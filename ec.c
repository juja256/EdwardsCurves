
#include "gf.h"
#include "ec.h"
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <stdio.h>

void EcWDoubleProj(Ec* ecc, const EcPointProj*, EcPointProj*);
void EcEdDoubleProj(Ec* ecc, const EcPointProj*, EcPointProj*);
void EcEdAddProj(Ec* ecc, const EcPointProj*, const EcPointProj*, EcPointProj*);
void EcWAddProj(Ec* ecc, const EcPointProj*, const EcPointProj*, EcPointProj*);

void PRNGInit(PRNG* generator, unsigned char* seed, int seed_len) {
    memcpy( generator->state, seed, seed_len );
    srand(time(NULL));
    generator->state[0] = 1;
    PRNGRun(generator);
}

void PRNGRun(PRNG* generator) {
    generator->state[0] ^= rand();
}

void PRNGGenerateSequence(PRNG* generator, int bit_len, unsigned char* dest) {
    int w = (bit_len % 64 == 0) ? bit_len/64 : bit_len/64+1; 
    memset(dest, 0, w*8);
    for (int i=0; i<bit_len; i++) {
        PRNGRun(generator);
        dest[i/8] ^= ( generator->state[0] & 1 ) << (i%8);
    }
}

int EcPointCmp(Ec* ecc, const EcPoint* A, const EcPoint* B) {
    return GFCmp(ecc, A->x, B->x) || GFCmp(ecc, A->y, B->y);
}

void EcCopy(Ec* ecc, EcPoint* dest, const EcPoint* src) {
    copy(dest->x, src->x, ecc->wordLen);
    copy(dest->y, src->y, ecc->wordLen);
}

void EcCopyProj(Ec* ecc, EcPointProj* dest, const EcPointProj* src) {
    copy(dest->X, src->X, ecc->wordLen);
    copy(dest->Y, src->Y, ecc->wordLen);
    copy(dest->Z, src->Z, ecc->wordLen);
}

void EcIdentityPoint(Ec* ecc, EcPoint* P) {
    if (ecc->isEdwards) {
        EcCopy(ecc, P, &uPEd);
    }
    else {
        copy(P->x, zero, ecc->wordLen);
        copy(P->y, zero, ecc->wordLen);
    }
}

void EcIdentityPointProj(Ec* ecc, EcPointProj* P) {
    if (ecc->isEdwards) {
        EcCopyProj(ecc, P, &uPPEd);
    }
    else {
        EcCopyProj(ecc, P, &uPPW);
    }
}

void EcDestroy(Ec* ecc) {
    if (ecc->T != NULL)
        free(ecc->T);
}

void BaseEcInit(Ec* ecc, u64 bitLen, const BigInt p, const EcPoint* bp, const BigInt n)
{
    srand(time(NULL));
    ecc->bitLen = bitLen;
    ecc->wordLen = (bitLen%64 == 0) ? (bitLen / 64) : (bitLen / 64) + 1;
    memset(ecc->p, 0, sizeof(BigInt));
    memset(ecc->n, 0, sizeof(BigInt));
    copy(ecc->n, n, ecc->wordLen);

    ecc->GFMul = GFMul_Cmn;
    ecc->GFSqr = GFSqr_Cmn;

    copy(ecc->p, p, ecc->wordLen);

    
    if (ecc->isEdwards) {
        ecc->curve_id = ED_NOT_STANDARD;
    }
    else {
        ecc->curve_id = FIPS_NOT_STANDARD;
    }
    if (GFCmp(ecc, ecc->p, p192) == 0) {
        ecc->GFMul = GFMul_FIPS192;
        ecc->GFSqr = GFSqr_FIPS192;
        ecc->curve_id |= 0x192;
    }
    else if (GFCmp(ecc, ecc->p, p224) == 0) {
        ecc->GFMul = GFMul_FIPS224;
        ecc->GFSqr = GFSqr_FIPS224;
        ecc->curve_id |= 0x224;
    }
    else if (GFCmp(ecc, ecc->p, p256) == 0) {
        ecc->GFMul = GFMul_FIPS256;
        ecc->GFSqr = GFSqr_FIPS256;
        ecc->curve_id |= 0x256;
    }
    else if (GFCmp(ecc, ecc->p, p384) == 0) {
        ecc->GFMul = GFMul_FIPS384;
        ecc->GFSqr = GFSqr_FIPS384;
        ecc->curve_id |= 0x384;
    }
    else if (GFCmp(ecc, ecc->p, p521) == 0) {
        ecc->GFMul = GFMul_FIPS521;
        ecc->GFSqr = GFSqr_FIPS521;
        ecc->curve_id |= 0x521;
    }

    unsigned char seed[4];
    srand(time(NULL));
    seed[0] = rand();
    PRNGInit(&(ecc->prng), seed, 4);

    if (bp != NULL) {
        copy(ecc->BasePoint.x, bp->x, ecc->wordLen);
        copy(ecc->BasePoint.y, bp->y, ecc->wordLen);
    }
}

int EcEdInit(EcEd* ecc, u64 bitLen, const BigInt p, const EcPoint* bp, const BigInt n, const GFElement d) {
    ecc->isEdwards = 1;
    BaseEcInit(ecc, bitLen, p, bp, n);
    ecc->EcAdd = EcEdAddProj;
    ecc->EcDouble = EcEdDoubleProj;

    copy(ecc->d, d, ecc->wordLen);
    copy(ecc->a, unity, ecc->wordLen);
    EcScalarMulWindowedPrecomputation(ecc, &(ecc->BasePoint), &(ecc->T), WINDOW_SIZE);
    return 0;
}

int EcEdTwistedUAInit(EcEd* ecc, u64 bitLen, const BigInt p, const EcPoint* bp, const BigInt n, const GFElement d, u32 curve_id) {
    ecc->isEdwards = 1;
    BaseEcInit(ecc, bitLen, p, bp, n);
    ecc->curve_id = curve_id;
    ecc->hash_id = KUPYNA_HASH;
    ecc->EcAdd = EcEdAddProj;
    ecc->EcDouble = EcEdDoubleProj;

    copy(ecc->d, d, ecc->wordLen);
    copy(ecc->a, zero, ecc->wordLen);
    ecc->a[0] = 2; // UA curves
    EcScalarMulWindowedPrecomputation(ecc, &(ecc->BasePoint), &(ecc->T), WINDOW_SIZE);

    /* check if d, a both are not quadratic residues for the lack of convenience */
    GFElement pp, tmp;
    copy(pp, ecc->p, ecc->wordLen);
    div2(ecc->wordLen, pp);
    GFPow(ecc, ecc->d, pp, tmp);
    if (GFCmp(ecc, tmp, unity) == 0) return INVALID_D; // invalid d

    GFPow(ecc, ecc->a, pp, tmp);
    if (GFCmp(ecc, tmp, unity) == 0) return INVALID_A; // invalid a
    switch (ecc->bitLen)
    {
    case 256:
        ecc->max_msg_size = P256_MAX_MSG_SIZE;
        ecc->hash_out_size = P256_HASH_SIZE;
        break;
    case 384:
        ecc->max_msg_size = P384_MAX_MSG_SIZE;
        ecc->hash_out_size = P384_HASH_SIZE;
        break;
    case 512:
        ecc->max_msg_size = P512_MAX_MSG_SIZE;
        ecc->hash_out_size = P512_HASH_SIZE;
        break;
    case 768:
        ecc->max_msg_size = P768_MAX_MSG_SIZE;
        ecc->hash_out_size = P768_HASH_SIZE;
        break;
    default:
        ecc->max_msg_size = 0;
        ecc->hash_out_size = 0;
    }
    return 0;
}

int EcWInit(EcW* ecc, u64 bitLen, const BigInt p, const EcPoint* bp, const BigInt n, const GFElement a, const GFElement b) {  
    ecc->isEdwards = 0;
    BaseEcInit(ecc, bitLen, p, bp, n); 

    ecc->EcAdd = EcWAddProj;
    ecc->EcDouble = EcWDoubleProj;

    copy(ecc->a, a, ecc->wordLen);
    copy(ecc->b, b, ecc->wordLen);
    EcScalarMulWindowedPrecomputation(ecc, &(ecc->BasePoint), &(ecc->T), WINDOW_SIZE);
    return 0;
}

int EcIsIdentityPoint(Ec* ecc, const EcPoint* P) {
    if (ecc->isEdwards) {
        return EcPointCmp(ecc, P, &uPEd) == 0;
    }
    else {
        return EcPointCmp(ecc, P, &uPW) == 0;
    }
}

int EcEdCheckPointOnCurve(EcEd* ecc, const EcPoint* P) {
    GFElement x,y,z,xy;
    GFSqr(ecc, P->x, x);
    GFSqr(ecc, P->y, y);
    GFMul(ecc, x, y, xy);
    GFMulByD(ecc, xy);
    GFAdd(ecc, xy, unity, xy);

    GFMulByA(ecc, y, y);
    GFAdd(ecc, x, y, z);
    return !GFCmp(ecc, xy, z); //1 - ok
}

int EcWCheckPointOnCurve(EcW* ecc, const EcPoint* P) {
    GFElement l,r; //left,right
    GFSqr(ecc, P->y, l);
    GFSqr(ecc, P->x, r);
    GFAdd(ecc, r, ecc->a, r);
    GFMul(ecc, r, P->x, r); // now r = x(x^2 + a)
    GFAdd(ecc, r, ecc->b, r);
    return !GFCmp(ecc, l, r); //1 - ok
}

int EcCheckPointOnCurve(Ec* ecc, const EcPoint* P) {
    if (ecc->isEdwards) return EcEdCheckPointOnCurve(ecc, P);
    else return EcWCheckPointOnCurve(ecc, P);
}

int EcCheckPointInMainSubGroup(Ec* ecc, const EcPoint* P) {
    int r = EcCheckPointOnCurve(ecc, P);
    if (!r) return 0;
    EcPoint Q;
    EcScalarMul(ecc, P, ecc->n, &Q);
    return EcIsIdentityPoint(ecc, &Q);
}

void EcEdAddAf(EcEd* ecc, const EcPoint* A, const EcPoint* B, EcPoint* C) {
    GFElement z1, z2, z3, z4, z5, z6, z7;
    GFMul(ecc, A->x, B->x, z1); // z1 = x1 * x2
    GFMul(ecc, A->y, B->y, z2); // z2 = y1 * y2
    GFMul(ecc, z1, z2, z3); 
    GFMulByD(ecc, z3); // z3 = d * x1 * x2 * y1 * y2

    GFNeg(ecc, z3, z4); // z4 = - z3
    GFAdd(ecc, z3, unity, z3);

    GFInv(ecc, z3, z3);

    GFAdd(ecc, z4, unity, z4);

    GFInv(ecc, z4, z4);

    GFMul(ecc, A->x, B->y, z5); // z5 = x1 * y2
    GFMul(ecc, A->y, B->x, z6); // z6 = x2 * y1
    GFAdd(ecc, z5, z6, z5);
    GFMulByA(ecc, z2, z2);
    GFSub(ecc, z1, z2, z2);

    GFMul(ecc, z5, z3, C->y);
    GFMul(ecc, z2, z4, C->x);
}

int EcWAddAf(EcW* ecc, const EcPoint* A, const EcPoint* B, EcPoint* C)
{
    GFElement s,difX,difY;
    GFSub(ecc,A->x,B->x,difX); // difX = x1-x2
    GFSub(ecc,A->y,B->y,difY); // difY = y1-y2
    if (GFCmp(ecc, difX, zero) == 0) return INFINITY_POINT;
    //s = difY/difX
    GFInv(ecc, difX, difX);
    GFMul(ecc, difY, difX, s);
    //x3 = s^2 - x1 - x2
    GFSqr(ecc, s, difX);
    GFSub(ecc, difX, A->x, difX);
    GFSub(ecc, difX, B->x, difX); // x3 = difX

    //y3 = -y1 + s(x1 - x3)
    GFSub(ecc, A->x, difX, difY); //y3 = x1 - x3
    GFMul(ecc, difY, s ,difY); //y3 = y3 * s
    GFSub(ecc, difY, A->y, C->y); //y3 = y3 - y1
    copy(C->x, difX, ecc->wordLen);
    return NORMAL_POINT;
}

int EcAddAf(Ec* ecc, const EcPoint* A, const EcPoint* B, EcPoint* C) {
    if (ecc->isEdwards) { EcEdAddAf(ecc, A, B, C); return NORMAL_POINT; }
    else return EcWAddAf(ecc, A, B, C);
}

void EcEdDoubleAf(EcEd* ecc, const EcPoint* A, EcPoint* B) {
    GFElement z1, z2, z3, z4, z5;
    GFSqr(ecc, A->x, z1);
    GFSqr(ecc, A->y, z2);
    GFMul(ecc, A->x, A->y, z3);
    GFMulBy2(ecc, z3, z3);
    GFMul(ecc, z1, z2, z4);
    GFMulByD(ecc, z4);
    GFNeg(ecc, z4, z5);
    GFAdd(ecc, z4, unity, z4);
    GFInv(ecc, z4, z4);
    GFAdd(ecc, z5, unity, z5);
    GFInv(ecc, z5, z5);
    GFMulByA(ecc, z2, z2);
    GFSub(ecc, z1, z2, z2);
    GFMul(ecc, z4, z3, B->y);
    GFMul(ecc, z2, z5, B->x);
}

int EcWDoubleAf(EcW* ecc, const EcPoint* A, EcPoint* B)
{
    if (GFCmp(ecc, A->y, zero) == 0) return INFINITY_POINT;
    GFElement s,tmp,tmp2;
    copy(tmp2, A->x, ecc->wordLen);

    //s = (3*x1^2 + a)/2/y1
    GFSqr(ecc, A->x, s);
    GFMulBy2(ecc, s, tmp);  
    GFAdd(ecc, tmp, s, s); 

    GFAdd(ecc, s, ecc->a, s);

    GFMulBy2(ecc, A->y, tmp); 
    GFInv(ecc, tmp, tmp);

    GFMul(ecc, s, tmp, s);

    //x2 = s^2 - 2*x1
    GFSqr(ecc, s, tmp);
    GFMulBy2(ecc, A->x, B->x);
    GFSub(ecc, tmp, B->x, B->x);

    //y2 = -y1 + s(x1 - x2)
    GFSub(ecc, tmp2, B->x, tmp);
    GFMul(ecc, tmp, s, tmp);
    GFSub(ecc, tmp, A->y, B->y);
    return NORMAL_POINT;
}

int EcDoubleAf(Ec* ecc, const EcPoint* A, EcPoint* B) {
    if (ecc->isEdwards) { EcEdDoubleAf(ecc, A, B); return NORMAL_POINT; }
    else return EcWDoubleAf(ecc,A,B);
}

/* x^2 + y^2 = 1 + dx^2y^2 */
/* y^2 = (1 - x^2)/(1 - dx^2) */
void EcEdGenerateBasePoint(EcEd* ecc, EcPoint* bp) {
    PRNGGenerateSequence( &ecc->prng, ecc->bitLen, (u8*)(bp->x) );
    int isEqual = 0;
    GFElement t1, t2, t3, pp;
    copy(pp, ecc->p, ecc->wordLen);
    div2(ecc->wordLen, pp);

    do {
        GFSqr(ecc, bp->x, t1);
        GFNeg(ecc, t1, t2);
        GFAdd(ecc, t2, unity, t2);
        GFMulByD(ecc, t1);
        GFNeg(ecc, t1, t1);
        GFAdd(ecc, t1, ecc->a, t1);
        GFInv(ecc, t1, t1);
        GFMul(ecc, t1, t2, t1);
        GFPow(ecc, t1, pp, t2);
        isEqual = !GFCmp(ecc, t2, unity);
        if (!isEqual)
            GFAdd(ecc, bp->x, unity, bp->x);
        else break;
    } while(1);

    GFSqrt(ecc, t1, bp->y);
    EcDoubleAf(ecc, bp, bp);
    EcDoubleAf(ecc, bp, bp); /* Keeping Point in Main SubGroup by multiplication by factor (x4) */
}

//y^2 = x^3 + ax + b
void EcWGenerateBasePoint(EcW* ecc, EcPoint* bp) {
    PRNGGenerateSequence( &ecc->prng, ecc->bitLen, (u8*)(bp->x) );
    int y = 0;
    GFElement sum, pp,tmp;
    copy(pp, ecc->p, ecc->wordLen);
    div2(ecc->wordLen, pp);

    do {
        GFSqr(ecc,bp->x,sum);
        GFAdd(ecc,sum,ecc->a,sum);
        GFMul(ecc,sum,bp->x,sum);
        GFAdd(ecc,sum,ecc->b,sum);
        GFPow(ecc, sum, pp, tmp);
        y = !GFCmp(ecc, tmp, unity);
        if (!y)
            GFAdd(ecc, bp->x, unity, bp->x);
        else break;
    } while(1);

    GFSqrt(ecc, sum, bp->y);
}

void EcGenerateBasePoint(Ec* ecc, EcPoint* bp) {
    if (ecc->isEdwards) EcEdGenerateBasePoint(ecc, bp);
    else EcWGenerateBasePoint(ecc,bp);
}

void EcConvertAffineToProjective(Ec* ecc, const EcPoint* P, EcPointProj* Q) {
    PRNGGenerateSequence(&ecc->prng, ecc->bitLen, (u8*)(Q->Z) );
    GFMul(ecc, P->x, Q->Z, Q->X);
    GFMul(ecc, P->y, Q->Z, Q->Y);
}

void EcConvertProjectiveToAffine(Ec* ecc, const EcPointProj* P, EcPoint* Q) {
    GFElement Z_inv;
    if (GFCmp(ecc, P->Z, zero) == 0) {
        EcCopy(ecc, Q, &uPW);
    }
    GFInv(ecc, P->Z, Z_inv);
    GFMul(ecc, P->X, Z_inv, Q->x);
    GFMul(ecc, P->Y, Z_inv, Q->y);
}

/* According to http://elibrary.kubg.edu.ua/21879/1/A_Bessalov_Polytechnika_2017_FIT.pdf */
void EcEdAddProj(EcEd* ecc, const EcPointProj* P1, const EcPointProj* P2, EcPointProj* P3) {
    /* 10M + 1S */
    GFElement A, B, C, D, E, F, G, T;
    GFMul(ecc, P1->Z, P2->Z, A); // A = Z1Z2
    GFSqr(ecc, A, B); // B = A^2
    GFMul(ecc, P1->X, P2->X, C); // C = X1X2
    GFMul(ecc, P1->Y, P2->Y, D); // D = Y1Y2
    GFMul(ecc, C, D, E); 
    GFMulByD(ecc, E); // E = dCD
    GFSub(ecc, B, E, F); // F = B-E
    GFAdd(ecc, B, E, G); // G = B+E

    
    GFAdd(ecc, P1->X, P1->Y, T);
    GFAdd(ecc, P2->X, P2->Y, P3->Y);
    GFMul(ecc, P3->Y, T, P3->Y);
    GFSub(ecc, P3->Y, C, P3->Y);
    GFSub(ecc, P3->Y, D, P3->Y);
    GFMul(ecc, P3->Y, A, P3->Y);
    GFMul(ecc, P3->Y, F, P3->Y); // Y3 = AF((X1+Y1)(X2+Y2)-C-D)

    GFMulByA(ecc, D, D);
    GFSub(ecc, C, D, P3->X);
    GFMul(ecc, P3->X, A, P3->X);
    GFMul(ecc, P3->X, G, P3->X); // X3 = AG(C-D) 

    GFMul(ecc, F, G, P3->Z); // Z3 = FG
}

void EcWAddProj(EcW* ecc, const EcPointProj* A, const EcPointProj* B, EcPointProj* C)
{
    if ( GFCmp(ecc, A->Z, zero) == 0 ) { // Handling infinities
        EcCopyProj(ecc, C, B);
        return;
    }
    if ( GFCmp(ecc, B->Z, zero) == 0 ) { // Handling infinities
        EcCopyProj(ecc, C, A);
        return;
    }

    GFElement u,v,a,tmp,x1z2,sqrv,y1z2;
    //v = x2z1 - x1z2
    GFMul(ecc,A->X,B->Z,x1z2); 
    GFMul(ecc,A->Z,B->X,v);

    GFSub(ecc,v,x1z2,v);
    //u = y2z1 - y1z2
    GFMul(ecc,B->Y,A->Z,u);
    GFMul(ecc,B->Z,A->Y,y1z2);
    GFSub(ecc,u,y1z2,u);
    
    GFMul(ecc,A->Z,B->Z,tmp); // tmp = z1*z2
    GFSqr(ecc,v,sqrv); 
    //z3 = v^3*z1*z2
    GFMul(ecc,sqrv,v,C->Z);
    GFMul(ecc,C->Z,tmp,C->Z);
    //a = u^2*z1*z2 - v^2(v+2*x1*x2)
    GFSqr(ecc,u,a);
    GFMul(ecc,a,tmp,a);
    GFMulBy2(ecc,x1z2,tmp);
    GFAdd(ecc,tmp,v,tmp);
    GFMul(ecc,tmp,sqrv,tmp);
    GFSub(ecc,a,tmp,a);
    //x3 = va
    GFMul(ecc,v,a,C->X);
    //y3 = v^2(u*x1*z2 - v*y1*z2) - u*a 
    GFMul(ecc,u,x1z2,C->Y);
    //GFMul(ecc,A->Y,B->Z,tmp);
    GFMul(ecc,y1z2,v,tmp);
    GFSub(ecc,C->Y,tmp,C->Y);
    GFMul(ecc,C->Y,sqrv,C->Y);
    GFMul(ecc,u,a,tmp);
    GFSub(ecc,C->Y,tmp,C->Y);
}

void EcAddProj(Ec* ecc, const EcPointProj* A, const EcPointProj* B, EcPointProj* C) {
    ecc->EcAdd(ecc, A, B, C);
}


void EcEdDoubleProj(EcEd* ecc, const EcPointProj* P, EcPointProj* P2) {
    /* 4M + 3S */
    GFElement A,B,C,D,E,F,G;

    GFSqr(ecc, P->X, A); // A = X^2
    GFSqr(ecc, P->Y, B); // B = aY^2
    GFSqr(ecc, P->Z, C); // C = Z^2
    GFMulByA(ecc, B, B); // B = aY^2
    GFAdd(ecc, A, B, D); // D = A+B

    // G = 2XY
    GFMul(ecc, P->X, P->Y, G);
    GFMulBy2(ecc, G, G);
    
    GFSub(ecc, A, B, E); // E = A-B
    GFMulBy2(ecc, C, F); 
    GFSub(ecc, F, D, F); // F = 2C - A - B
    
    GFMul(ecc, F, G, P2->Y); // Y2 = FG 
    GFMul(ecc, D, E, P2->X); // X2 = DE
    GFMul(ecc, D, F, P2->Z); // Z3 = DF 
}

void EcWDoubleProj(EcW* ecc, const EcPointProj* A, EcPointProj* B)
{
    GFElement h,s,w,b,tmp;
    //w = a*z1^2 + 3*x1^2
    GFSqr(ecc,A->X,tmp); 
    GFMulBy2(ecc,tmp,w);
    GFAdd(ecc,w,tmp,w); //now w = 3*x1^2
    GFSqr(ecc,A->Z,tmp);
    GFMul(ecc,tmp,ecc->a,tmp);
    GFAdd(ecc,w,tmp,w);
    //s = y1*z1
    GFMul(ecc,A->Y,A->Z,s);
    //b = x1*y1*s
    GFMul(ecc,A->X,A->Y,b);
    GFMul(ecc,b,s,b);
    //h = w^2-8b
    GFSqr(ecc,w,h);

    GFMulBy2Power(ecc,b,3,tmp);

    GFSub(ecc,h,tmp,h);
    //x2 = 2hs
    GFMul(ecc,h,s,B->X);
    GFMulBy2(ecc,B->X,B->X);
    //y2 = w(4b - h) - 8(y1*s)^2
    GFMul(ecc,A->Y,s,tmp);
    GFMulBy2(ecc,tmp,tmp);
    GFSqr(ecc,tmp,tmp);
    GFMulBy2(ecc,tmp,tmp); // now tmp = 8*(y1*s)^2
    GFMulBy2Power(ecc,b,2,B->Y);

    GFSub(ecc,B->Y,h,B->Y);
    GFMul(ecc,w,B->Y,B->Y);
    GFSub(ecc,B->Y,tmp,B->Y);
    //z2 = 8s^3
    GFMulBy2(ecc,s,tmp);
    GFSqr(ecc,tmp,B->Z);
    GFMul(ecc,B->Z,tmp,B->Z);
}

void EcDoubleProj(Ec* ecc, const EcPointProj* A, EcPointProj* B) {
    ecc->EcDouble(ecc, A, B);
}

/* Easiest AddAndDouble ScalarMul */
void EcScalarMulNaive(Ec* ecc, const EcPointProj* A, const BigInt k, EcPointProj* B) {
    EcPointProj H;
    
    EcCopyProj(ecc, &H, A); // H := A
    EcIdentityPointProj(ecc, B);

    for (u32 i=0; i<ecc->bitLen; i++) {
        if (get_bit(k, i)) {
            ecc->EcAdd(ecc, B, &H, B);
        }
        ecc->EcDouble(ecc, &H, &H); 
    }
}

/* Constant-time AlwaysAddAndDouble Montgomery ScalarMul */
void EcScalarMulMontgomery(Ec* ecc, const EcPointProj* A, const BigInt k, EcPointProj* B) {
    EcPointProj H;
    
    EcCopyProj(ecc, &H, A); // H := A, H = P1
    EcIdentityPointProj(ecc, B);

    for (int i=ecc->bitLen-1; i>=0; i--) {
        if (get_bit(k, i) == 0) {
            ecc->EcAdd(ecc, B, &H, &H);
            ecc->EcDouble(ecc, B, B); 
        }
        else {
            ecc->EcAdd(ecc, B, &H, B);
            ecc->EcDouble(ecc, &H, &H); 
        }
    }
}


/* Fixed-window mostly constant-time ScalarMul */
void EcScalarMulWindowedPrecomputation(Ec* ecc, const EcPoint* A, EcPointProj** T, int windowSize) {
    *T = (EcPointProj*)malloc(sizeof(EcPointProj) * (1<<windowSize) );
    BigInt k;
    EcPointProj bpp;
    EcConvertAffineToProjective(ecc, A, &bpp);
    copy(k, zero, ecc->wordLen);
    for (int i=0; i<(1<<windowSize); i++) {
        EcScalarMulNaive(ecc, &bpp, k, &((*T)[i]));
        k[0]++;
    }
}

void EcScalarMulWindowed(Ec* ecc, const EcPointProj* T, int windowSize, const BigInt k, EcPointProj* B) {
    int m = (ecc->bitLen % windowSize == 0) ? (ecc->bitLen / windowSize) : (ecc->bitLen / windowSize) + 1;
    EcIdentityPointProj(ecc, B);

    for (int i=m-1; i>=0; i--) {
        for (int j=0;j<windowSize;j++){
            ecc->EcDouble(ecc, B, B);
        }
        int n = (k[ (i*windowSize) / 64 ] & ( (u64)( (1<<windowSize) - 1) << ((i*windowSize) % 64) )) >> ((i*windowSize) % 64);
        ecc->EcAdd(ecc, B, &(T[n]), B);
    }
}

typedef signed char BigIntNAF[sizeof(BigInt) * 8 + 1];

static inline char mods(BigInt d, int w) {
    if ( d[0] & (1<<w-1) ) return (char)(d[0] & ((1<<w)-1)) - (char)(1 << w);
    else return (char)(d[0] & ((1<<w)-1));
}

/* via Reitwiesner Algorithm */
static inline void ConvertBigIntTowNAF(Ec* ecc, const BigInt k, int w, BigIntNAF naf) {
    int i = 0;
    BigInt d;
    copy(d, k, ecc->wordLen);
    memset(naf, 0, sizeof(BigIntNAF));
    while (GFCmp(ecc, d, zero) != 0) {
        if ( (d[0] & 1) == 1) {
            naf[i] = mods(d, w);
            if (naf[i] >= 0) {
                sub_word(ecc->wordLen, d, naf[i], d);
            }
            else {
                add_word(ecc->wordLen, d, -naf[i], d);
            }
        }
        div2(ecc->wordLen, d);
        i++;
    }
}

void EcScalarMulwNAFPrecomputation(Ec* ecc, const EcPoint* A, EcPointProj** T, int windowSize) {
    *T = (EcPointProj*)malloc(sizeof(EcPointProj) * (1<<windowSize)/2 );
    BigInt k;
    EcPointProj bpp;
    EcConvertAffineToProjective(ecc, A, &bpp);
    copy(k, unity, ecc->wordLen);
    for (int i=0; i<(1<<windowSize)/2; i+=2) {
        EcScalarMulNaive(ecc, &bpp, k, &((*T)[i]));
        EcCopyProj(ecc, &((*T)[i+1]), &((*T)[i]));
        GFNeg(ecc, (*T)[i+1].Y, (*T)[i+1].Y);
        k[0]+=2;
    }
}

void EcScalarMulwNAF(Ec* ecc, const EcPointProj* T, int windowSize, const BigInt k, EcPointProj* B) {
    BigIntNAF naf;
    ConvertBigIntTowNAF(ecc, k, windowSize, naf);
    EcIdentityPointProj(ecc, B);

    for (int i=ecc->bitLen; i>=0; i--) {
        ecc->EcDouble(ecc, B, B);
        if (naf[i] != 0) {
            ecc->EcAdd(ecc, B, &(T[ abs(naf[i])/2 + (naf[i] < 0) ]), B);
        }
    }
}

void EcScalarMulByBasePoint(Ec* ecc, const BigInt k, EcPoint* B) {
    EcPointProj B_p;
    EcScalarMulWindowed(ecc, ecc->T, WINDOW_SIZE, k, &B_p);
    EcConvertProjectiveToAffine(ecc, &B_p, B);
}

void EcScalarMul(Ec* ecc, const EcPoint* A, const BigInt k, EcPoint* B) {
    EcPointProj A_p, B_p;
    EcConvertAffineToProjective(ecc, A, &A_p);
    EcScalarMulMontgomery(ecc, &A_p, k, &B_p);
    EcConvertProjectiveToAffine(ecc, &B_p, B);
}

/*
    Standard Edwards full curves:
        Ed-192:
    p = FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFF
    d = 6DBA6A
    n = 3FFFFFFFFFFFFFFFFFFFFFFFEA75D4027230DD4DFFDB0455
    x = 44F083BB00E51AD91A2743284D31F57EE5C84826FCC91F4B
    y = 15FC16E5870524E0DBBE9EC8BB9F066C02A02B1978D4E029

        Ed-224
    p = FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000000000000000000001
    d = 3608425
    n = 400000000000000000000000000020BBEC47CEDB34DD05BCB6B7E619
    x = C448CA02660F57204FF1BDE2B5CC3E25606A7460399FEA3DA9A06383
    y = 319117770D6FC7FE35F6A02905FE1F363156BD2E5B75BB89A64CAFAB

        Ed-256
    p = FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF
    d = 72A38
    n = 3FFFFFFFC00000003FFFFFFFFFFFFFFFBA76FA29C30CC3AA4954B53EDBE54D75
    x = 894F8283626AEE6848515DDDC3B8DBB3D5302DEE0EE75080D6753E4D39BA5AB2
    y = EA612346223F6480CBBAFA39DB95D54D21469DD3074A957EFDA4FD79FEB630B5

        Ed-384
    p = FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFF0000000000000000FFFFFFFF
    d = 12593A
    n = 4000000000000000000000000000000000000000000000005063576B5A9A0C3A23E9510EA680650B4884E63A2968DD71
    x = 1FC0E8E61F599813E376D11F7510D77F177C2F1CDE19FD14D63A2EC5EAD4D0DED1BD06676CCF365243BF3C0675A31B62
    y = F52B4FA352B257D7A102FA45C56A50CCBDB3DEC053D5610EDBD0188C11F321F28A43E2FC50395E4A8BD0029AE71D51AA

        Ed-521
    p = 1FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
    d = 16A
    n =  7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF46087BC4A294FCC80B3F45D8CEBFB21479BA651BA07DE913AD1D8392DE3FF8AF
    x =  749CC477889C4315E0097137C4704FA312C80BB1ED27658CB593229888C6091AB09D3AACFA68C24D60516EEEB707BB9D81DEBE5521819424F0E82014F35176A6DE
    y =  C3860A238AE5DB5A01BEFC4C26F6B041A6D6D2A16FA67D6336DC5D8084F702EBC5D741906190FCD3E42161594F7B58BD4B4502F331931EED042E9B77FC088A28C6

========================================================================================================

    Weierstrass NIST curves

    everywhere a = -3

        P-192
    p = FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFF
    n = FFFFFFFFFFFFFFFFFFFFFFFF99DEF836146BC9B1B4D22831
    a = FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFC
    b = 64210519E59C80E70FA7E9AB72243049FEB8DEECC146B9B1
    x = 188DA80EB03090F67CBF20EB43A18800F4FF0AFD82FF1012
    y = 07192B95FFC8DA78631011ED6B24CDD573F977A11E794811

        P-224
    p = FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000000000000000000001
    n = FFFFFFFFFFFFFFFFFFFFFFFFFFFF16A2E0B8F03E13DD29455C5C2A3D
    a = FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFE
    b = B4050A850C04B3ABF54132565044B0B7D7BFD8BA270B39432355FFB4
    x = B70E0CBD6BB4BF7F321390B94A03C1D356C21122343280D6115C1D21
    y = BD376388B5F723FB4C22DFE6CD4375A05A07476444D5819985007E34

        P-256
    p = FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF
    n = FFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551
    a = FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFC
    b = 5AC635D8AA3A93E7B3EBBD55769886BC651D06B0CC53B0F63BCE3C3E27D2604B
    x = 6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296
    y = 4FE342E2FE1A7F9B8EE7EB4A7C0F9E162BCE33576B315ECECBB6406837BF51F5

        P-384
    p = FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFF0000000000000000FFFFFFFF
    n = FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFC7634D81F4372DDF581A0DB248B0A77AECEC196ACCC52973
    a = FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFF0000000000000000FFFFFFFC
    b = B3312FA7E23EE7E4988E056BE3F82D19181D9C6EFE8141120314088F5013875AC656398D8A2ED19D2A85C8EDD3EC2AEF
    x = AA87CA22BE8B05378EB1C71EF320AD746E1D3B628BA79B9859F741E082542A385502F25DBF55296C3A545E3872760AB7
    y = 3617DE4A96262C6F5D9E98BF9292DC29F8F41DBD289A147CE9DA3113B5F0B8C00A60B1CE1D7E819D7A431D7C90EA0E5F

        P-521
    p = 1FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
    n = 1FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFA51868783BF2F966B7FCC0148F709A5D03BB5C9B8899C47AEBB6FB71E91386409
    a = 1FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFC
    b = 051953EB9618E1C9A1F929A21A0B68540EEA2DA725B99B315F3B8B489918EF109E156193951EC7E937B1652C0BD3BB1BF073573DF883D2C34F1EF451FD46B503F00
    x =  C6858E06B70404E9CD9E3ECB662395B4429C648139053FB521F828AF606B4D3DBAA14B5E77EFE75928FE1DC127A2FFA8DE3348B3C1856A429BF97E7E31C2E5BD66
    y = 11839296A789A3BC0045C8A5FB42C7D1BD998F54449579B446817AFBD17273E662C97EE72995EF42640C550B9013FAD0761353C7086A272C24088BE94769FD16650
*/
int EcInitStandardCurve(Ec* ecc, u64 bitLen, BOOL isEdwards) {
    GFElement p,n;
    EcPoint G;
    if (isEdwards) {
        GFElement d;
        switch (bitLen) {
            case 192:
            GFInitFromString(p,  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFF");
            GFInitFromString(n,  "3FFFFFFFFFFFFFFFFFFFFFFFEA75D4027230DD4DFFDB0455");
            GFInitFromString(d,  "6DBA6A");
            GFInitFromString(G.y,"44F083BB00E51AD91A2743284D31F57EE5C84826FCC91F4B");
            GFInitFromString(G.x,"15FC16E5870524E0DBBE9EC8BB9F066C02A02B1978D4E029");
            break;
            case 224:
            GFInitFromString(p,  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000000000000000000001");
            GFInitFromString(n,  "400000000000000000000000000020BBEC47CEDB34DD05BCB6B7E619");
            GFInitFromString(d,  "3608425");
            GFInitFromString(G.y,"C448CA02660F57204FF1BDE2B5CC3E25606A7460399FEA3DA9A06383");
            GFInitFromString(G.x,"319117770D6FC7FE35F6A02905FE1F363156BD2E5B75BB89A64CAFAB");
            break;
            case 256:
            GFInitFromString(p,  "FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF");
            GFInitFromString(n,  "3FFFFFFFC00000003FFFFFFFFFFFFFFFBA76FA29C30CC3AA4954B53EDBE54D75");
            GFInitFromString(d,  "72A38");
            GFInitFromString(G.y,"894F8283626AEE6848515DDDC3B8DBB3D5302DEE0EE75080D6753E4D39BA5AB2");
            GFInitFromString(G.x,"EA612346223F6480CBBAFA39DB95D54D21469DD3074A957EFDA4FD79FEB630B5");
            break;
            case 384:
            GFInitFromString(p,  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFF0000000000000000FFFFFFFF");
            GFInitFromString(n,  "4000000000000000000000000000000000000000000000005063576B5A9A0C3A23E9510EA680650B4884E63A2968DD71");
            GFInitFromString(d,  "12593A");
            GFInitFromString(G.y,"1FC0E8E61F599813E376D11F7510D77F177C2F1CDE19FD14D63A2EC5EAD4D0DED1BD06676CCF365243BF3C0675A31B62");
            GFInitFromString(G.x,"F52B4FA352B257D7A102FA45C56A50CCBDB3DEC053D5610EDBD0188C11F321F28A43E2FC50395E4A8BD0029AE71D51AA");
            break;
            case 521:
            GFInitFromString(p,  "1FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF");
            GFInitFromString(n,  "7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF46087BC4A294FCC80B3F45D8CEBFB21479BA651BA07DE913AD1D8392DE3FF8AF");
            GFInitFromString(d,  "16A");
            GFInitFromString(G.y,"1C4A8AA03837A73FC86D95F308BD9B738E2BFAB5B8ECCEB0EF90F02B5E90D1DF28D470C4F531212CD0E68F4E925E5ED82A74AB63335A9DB3E31650C6767EBA00681");
            GFInitFromString(G.x,"11365BAC31463E811ECB7FEF04F65DE5B8B15ADAFD72CC74EC804B9408EB31F0E44CE67EA36144A06B07D105B862E2493B556740F343DB87866E15C7BE9F3813666");
            break;
            default:
            return -1;
        }
        return EcEdInit(ecc, bitLen, p, &G, n, d);
    }
    else
    {
        GFElement a,b;
        switch (bitLen) {
            case 192:
            GFInitFromString(p, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFF");
            GFInitFromString(n, "FFFFFFFFFFFFFFFFFFFFFFFF99DEF836146BC9B1B4D22831");
            GFInitFromString(a, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFC");
            GFInitFromString(b, "64210519E59C80E70FA7E9AB72243049FEB8DEECC146B9B1");
            GFInitFromString(G.x, "188DA80EB03090F67CBF20EB43A18800F4FF0AFD82FF1012");
            GFInitFromString(G.y, "07192B95FFC8DA78631011ED6B24CDD573F977A11E794811");
            break;
            case 224:
            GFInitFromString(p, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000000000000000000001");
            GFInitFromString(n, "FFFFFFFFFFFFFFFFFFFFFFFFFFFF16A2E0B8F03E13DD29455C5C2A3D");
            GFInitFromString(a, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFE");
            GFInitFromString(b, "B4050A850C04B3ABF54132565044B0B7D7BFD8BA270B39432355FFB4");
            GFInitFromString(G.x, "B70E0CBD6BB4BF7F321390B94A03C1D356C21122343280D6115C1D21");
            GFInitFromString(G.y, "BD376388B5F723FB4C22DFE6CD4375A05A07476444D5819985007E34");
            break;
            case 256:
            GFInitFromString(p, "FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF");
            GFInitFromString(n, "FFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551");
            GFInitFromString(a, "FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFC");
            GFInitFromString(b, "5AC635D8AA3A93E7B3EBBD55769886BC651D06B0CC53B0F63BCE3C3E27D2604B");
            GFInitFromString(G.x, "6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296");
            GFInitFromString(G.y, "4FE342E2FE1A7F9B8EE7EB4A7C0F9E162BCE33576B315ECECBB6406837BF51F5");
            break;
            case 384:
            GFInitFromString(p, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFF0000000000000000FFFFFFFF");
            GFInitFromString(n, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFC7634D81F4372DDF581A0DB248B0A77AECEC196ACCC52973");
            GFInitFromString(a, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFF0000000000000000FFFFFFFC");
            GFInitFromString(b, "B3312FA7E23EE7E4988E056BE3F82D19181D9C6EFE8141120314088F5013875AC656398D8A2ED19D2A85C8EDD3EC2AEF");
            GFInitFromString(G.x, "AA87CA22BE8B05378EB1C71EF320AD746E1D3B628BA79B9859F741E082542A385502F25DBF55296C3A545E3872760AB7");
            GFInitFromString(G.y, "3617DE4A96262C6F5D9E98BF9292DC29F8F41DBD289A147CE9DA3113B5F0B8C00A60B1CE1D7E819D7A431D7C90EA0E5F");
            break;
            case 521:
            GFInitFromString(p,   "1FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF");
            GFInitFromString(n,   "1FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFA51868783BF2F966B7FCC0148F709A5D03BB5C9B8899C47AEBB6FB71E91386409");
            GFInitFromString(a,   "1FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFC");
            GFInitFromString(b,   "051953EB9618E1C9A1F929A21A0B68540EEA2DA725B99B315F3B8B489918EF109E156193951EC7E937B1652C0BD3BB1BF073573DF883D2C34F1EF451FD46B503F00");
            GFInitFromString(G.x,  "C6858E06B70404E9CD9E3ECB662395B4429C648139053FB521F828AF606B4D3DBAA14B5E77EFE75928FE1DC127A2FFA8DE3348B3C1856A429BF97E7E31C2E5BD66");
            GFInitFromString(G.y, "11839296A789A3BC0045C8A5FB42C7D1BD998F54449579B446817AFBD17273E662C97EE72995EF42640C550B9013FAD0761353C7086A272C24088BE94769FD16650");
            break;
            default:
            return -1;
        }
        return EcWInit(ecc, bitLen, p, &G, n, a, b);
    }
}

int  EcInitStandardCurveById(Ec* ecc, u32 id) {
    u32 type = (id & 0x0000F000) >> 12;
    u32 bitLen = (0xF & id) + ((0xF0 & id) >> 4)*10 + ((0xF00 & id) >> 8)*100;
    if (type == 0xF) return EcInitStandardCurve(ecc, bitLen, 0);
    else if (type == 0xE) return EcInitStandardCurve(ecc, bitLen, 1);

    /* Ukrainian twisted curves */
    GFElement p,n;
    EcPoint G;
    GFElement d;
    switch (id) {
        /* 256bit curves */
        case UA_256_1:
        GFInitFromString(p,  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFE4D");
        GFInitFromString(n,  "4000000000000000000000000000000029E26087789BC2815BDFF97093543CCF");
        GFInitFromString(d,  "18");
        GFInitFromString(G.y,"742F27A268641C9D7DDF69892BE3DF3D8F9CC52260B89A4953C8379C7C0A212B");
        GFInitFromString(G.x,"91F5D0E7E2D417E3108B13B075CDC7756045F8424479FCFE8F23D27250A0883F");
        break;
        case UA_256_2:
        GFInitFromString(p,  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFDB5");
        GFInitFromString(n,  "3FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFB99844FF43B7A9435320A8D19D9C043F");
        GFInitFromString(d,  "2A0");
        GFInitFromString(G.y,"324B41DE7F36BD8AE59F911EBA978E5954D25ED40C1A2F944BD2D7F25D25792A");
        GFInitFromString(G.x,"2845CBB60A38C03D484BA09942E63D8F3D7407F17409F7F5E07E90A934D1EA17");
        break;
        case UA_256_3:
        GFInitFromString(p,  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFC65");
        GFInitFromString(n,  "400000000000000000000000000000001E1D14B6ED53005C107434A962152593");
        GFInitFromString(d,  "6B");
        GFInitFromString(G.y,"9CE446260CB1DAB3A777A255036FF042F616AA43D60150CB818989491EC015FA");
        GFInitFromString(G.x,"62E70EAB243C5CB882E35E8BB99268736430B4E7E479D97A6D978874EBE5B1B7");
        break;
        case UA_256_4:
        GFInitFromString(p,  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFAED");
        GFInitFromString(n,  "4000000000000000000000000000000026AD0CE5A034E9970A0131EDD7919957");
        GFInitFromString(d,  "44");
        GFInitFromString(G.y,"1261ED0B35AA220C974EC82BD3748E1E97A5F7CA982E6CF6F67B0167B8BFA1C6");
        GFInitFromString(G.x,"DABB64DC63418A9366DB2CE2A3BF1A9F032A2AA3C714B98F26423A106B6B7729");
        break;
        case UA_256_5:
        GFInitFromString(p,  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF9FD");
        GFInitFromString(n,  "40000000000000000000000000000000412E994309D80E83AFA4642085233F05");
        GFInitFromString(d,  "170");
        GFInitFromString(G.y,"CC4A085C63F3F5A59306AF8B234F1DFA9A177A62922967D7250B99F797F4FAF9");
        GFInitFromString(G.x,"D9D411C794BC01F7A07DEC0ECD314F8DE16720623DF9BA1AE4E08339E067E746");
        break;

        /* 384bit curves */
        case UA_384_1:
        GFInitFromString(p,  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF79D");
        GFInitFromString(n,  "4000000000000000000000000000000000000000000000000AF905B73674AC7D4AF38C53331DC208A517DCB3F340EECF");
        GFInitFromString(d,  "214");
        GFInitFromString(G.y,"A9027207E88074F4AFA3D44D4590DD04BAFBF6AE3D321091F500C783F4707940B7F5EBDD93325C5391843F9A78526BB1");
        GFInitFromString(G.x,"F5FC151B6264CB53A4B879AA9A1F4A5156BBF063B56AAA912617C0E4CEFF15D2DF497D9AEF12374A22D22C3B402B2C71");
        break;
        case UA_384_2:
        GFInitFromString(p,  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF115");
        GFInitFromString(n,  "400000000000000000000000000000000000000000000000231C4F8D91A6B595EAA0789F9CCFF3C7FA9F05F6C028EACF");
        GFInitFromString(d,  "A8");
        GFInitFromString(G.y,"89FD5124BBCFC2FBFA908A0A2F8D46E9A443EA0D34A8101CC28EA068C32EEA2A49466ECEDD25E2DABECBDE0B016C8ACD");
        GFInitFromString(G.x,"45DF45F8020CA03A0DE417410BDA8BFB795A653450104321344EFB9A1C5086B25A970EC79ED51DCA9D362AFF9A86F528");
        break;
        case UA_384_3:
        GFInitFromString(p,  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEF45");
        GFInitFromString(n,  "400000000000000000000000000000000000000000000000169BD1B903298994262820E2A3E73B4AA252FB924D3974A1");
        GFInitFromString(d,  "43");
        GFInitFromString(G.y,"20520FB9FCABB124725D8FF3AD5B2004DC296B245A1D57C689962B2A9D0EABCFD94CC4157E7CD0809055711EB35DEB5A");
        GFInitFromString(G.x,"25DA67A2BCA7B317923BD7537E3F781E433EB28F0224070701BBDF3DF0D63D26BA5409747D7AB030E2BB6FE78593C0C8");
        break;
        case UA_384_4:
        GFInitFromString(p,  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFE395");
        GFInitFromString(n,  "40000000000000000000000000000000000000000000000050BFE5AEE2F9642045CBA3B5673EC7020317A99B7A2A7F97");
        GFInitFromString(d,  "304");
        GFInitFromString(G.y,"8242F49FA453B08E85B033A3399C0EECFB57522496A9DFE1FC8423C009951649DD1FF90F9B49AA9C6E28AD374A67B632");
        GFInitFromString(G.x,"F108BBC769E3358598060F9C1C5EB79C6492860F58967AC7D7E50B23B21767B3D8A360C346B2AA22A801E0DD652CAD2A");
        break;
        case UA_384_5:
        GFInitFromString(p,  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFE2D5");
        GFInitFromString(n,  "4000000000000000000000000000000000000000000000007269FD2518A6AC30A3FA9F3931A347FA8BAC4843ECF19B5F");
        GFInitFromString(d,  "114");
        GFInitFromString(G.y,"535E607ABF5739C3E5A946F9AD44C899C33BF00B04CAE695080D622B28DA64CA0B55AB041BF971382D956739A3AAE176");
        GFInitFromString(G.x,"E44A4A622DAA298B94A167989E99ADFF06180D11EAA618386D7B173C82C4C1E9822C47B7C7F5B26663C8BF3074C4762C");
        break;

        /* 512bit curves */
        case UA_512_1:
        GFInitFromString(p,  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFC95");
        GFInitFromString(n,  "400000000000000000000000000000000000000000000000000000000000000028A3CE52209E2BD4952882D5574165192C46C0D0311FEA6BF9FECE70EE63B59F");
        GFInitFromString(d,  "10D");
        GFInitFromString(G.y,"53A0D50CC63C9219762F451978AEF214DBCFCC3A5CB5EF27124991A86B42B3A1A832724A0E6B930FDD1DA2E27A540D6B675E4422C444F529C508F0BAE7D0A85");
        GFInitFromString(G.x,"5230A1EE747050A072BD7319741586EA520388B6B53094571C821A2FC9A9E83D56665346B5DB04C43E75261DBDA512728FAAFAC48AE9260A5A184E2933E3A400");
        break;
        case UA_512_2:
        GFInitFromString(p,  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF8DD");
        GFInitFromString(n,  "3FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFDC46E59CE4CBBB6248CFE7BB582E1F490FD5F9F74433AE745EFE6351254BCDA5");
        GFInitFromString(d,  "179");
        GFInitFromString(G.y,"A51291E6175265F4C5755445CF271E51E9248C728D01CC388C3938E477F3A94C69F21873F7CB97C64AE7D33A1056B41BB0E95436F88A03555C72911F313331EF");
        GFInitFromString(G.x,"A9297AC52D51EF4745E5E558F6217B302052A39544DFE6DC04FE427850225E232BC379E11F6C51FDB0C9D0C2511BC3E945E25B9829B00C4DD92843328606ADFB");
        break;
        case UA_512_3:
        GFInitFromString(p,  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF27D");
        GFInitFromString(n,  "40000000000000000000000000000000000000000000000000000000000000002DD8135044F2810D554443F2E2BD624D526C8798402AB63951BFA991E6087F7B");
        GFInitFromString(d,  "AF");
        GFInitFromString(G.y,"81177507F05F730D44D26B9D14B0FFAC6558F6834314F0414A9D6A9E48608C5F08507D8273AF1326ABD50A49B901462F3BE8005B63FC8FF851B853FF48B31B7B");
        GFInitFromString(G.x,"A8B673C87778EE5D31B9A404334AB9B5A572C8EB536C2C91443F112F8961B1E1BBE4E7AC851E8AA509F0B9AD2D9C203A816B67E0796AD7FF29726D262EA957B8");
        break;
        case UA_512_4:
        GFInitFromString(p,  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEE85");
        GFInitFromString(n,  "4000000000000000000000000000000000000000000000000000000000000000153ED5CAB52ECFB1E3832493DD7A1C06617A6790626B88DFAC45D806A147B99B");
        GFInitFromString(d,  "3C9");
        GFInitFromString(G.y,"AB0076944360C99551AE3AB73E24473AD55987DFD6378DD3A6D89A4B2C88576F837341D4141E791B512483B6AEED65E09F68A64CE5B3B78AF16924634C2BCF13");
        GFInitFromString(G.x,"54123C034D18682F739A011463EE02332B35825602C853A8EDE6DD1657894A47C2D6E1D2EF39DC434504051DB60A690FAFFA10B9053E37268D3AF352F18B423A");
        break;
        case UA_512_5:
        GFInitFromString(p,  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFE64D");
        GFInitFromString(n,  "400000000000000000000000000000000000000000000000000000000000000017F1EF80AF949ECB3AC390D869E51ACD0E05576A5B2F6A96416B8365C69B687B");
        GFInitFromString(d,  "55");
        GFInitFromString(G.y,"B950D70E217E95B5AD76310129D34F735AE5CDB08DC80BEFD144C83B9E1159726BCFFDE8E806C46852123673FD81541FB0CEDF946E055CFAA4FBDA6AAF2BC51");
        GFInitFromString(G.x,"FF87729EF99568229871AC58E64867874334F44B071650877091558128573484A502DA64BDBCC347084B79A4F4A895D73D47FE32D0512E0796F1788C0A3084CE");
        break;

        /* 768bit curves */
        case UA_768_1:
        GFInitFromString(p,  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFDA45");
        GFInitFromString(n,  "400000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000443A2EC4A0EE2AA34439FF2239B2E0A486E21A23CCBF0266B649F83B9C234F6A0FF7785DACE4C8F581D4D96C68345E43");
        GFInitFromString(d,  "25C");
        GFInitFromString(G.y,"8E507CB7494E04A4624223110A0FDBD63CDC40376E4786FF960B8D39B25E949C96EA40BBC4A59CB9DD902CCFFEE583BAA7D0FCE6E0D9F75742B08A13E29068E39DCE12816E52FA1411D89C76476012C979D93A863268D0D096CC98EE9A894D23");
        GFInitFromString(G.x,"7A4301399505745D99DC7F4FC4E931300E4DC193DF1E7BA1B579ABDE4D2E32BEBFB593B148406544D16D5585DF369A7AF65C2D6258F48FC8C5D0E42D51EBE8B250B3B1DAB47FC3698FCAAA86C62CCFE83E887850D0995A58D41471AB9285BED8");
        break;
        case UA_768_2:
        GFInitFromString(p,  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFD555");
        GFInitFromString(n,  "4000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000002D53DC8D2F1FC7CB83E33D94852D3EE414186D796A86D35A580E5DBACCC74E280BB2C0E44353D97CE1BFD24110542DCD");
        GFInitFromString(d,  "636");
        GFInitFromString(G.y,"E9B840CAAD47C510B54EFFF4785577F18308F04D91CA92838072C162E18824BA395E2902993040DB7C2584DB1ADDE382A2DFAEC5FF23DD7814AB592452FC387BC032182A0DCA6BC7C6C089016173824516DAF4D6AD3CED3794140A11E4B2D106");
        GFInitFromString(G.x,"FDFFC92CE4BE1B293E4A7B2704709B63B285F97FF072A28118A7A927B695F91871C3EA91B0503FCDD0C566CC4F05DBB38589569C711A2E17DB7E70FE34AA92B8E03DA1E977DFC73D6743D8A487CEE06B9F960116980B66630CA3AF257DF435BC");
        break;
        case UA_768_3:
        GFInitFromString(p,  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFCD5D");
        GFInitFromString(n,  "40000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000019307FDF78EEA1D089460A6BD798E03A0E6ACEFB1335695599829531A05EA69449290BE0671BE46140E0F45945DCB667");
        GFInitFromString(d,  "2F7");
        GFInitFromString(G.y,"C5A9126B7300C113878313AE14F28FA3F0EDE0976353F8EAF1951ABC09C708DFB73586A6B04CAC4A243E6B3A213F1E131ADF53335A05B4EF879EA145F24028C52896B383114723E0268864CD18A2FF72B6819F1E355E4BAF46C4E961FC205B41");
        GFInitFromString(G.x,"CAB4B564E23DA773CBAB0253493D7A3B86B3BB0E37AC4016C01FBDE375DC3385047E5151F4F4F033C9B32600F1F7E6C3C1B6AA6E00F4E0919063CEE8658D9BFD8C199A63AAFDE6E71C53782600E72EE4F6C92C75FD6AF4FEDF385CCB8F5AC0AF");
        break;
        case UA_768_4:
        GFInitFromString(p,  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFC1A5");
        GFInitFromString(n,  "3FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF9E89FC12B0FA8CF0C1B7F9863B4F6BE4401CE25414828C3ED6B47FDE5AD47D69ADB14B80F47BEF0640DEB581217C9C0F");
        GFInitFromString(d,  "463");
        GFInitFromString(G.y,"41812102C2AED66F1A5E26FA0EFC9CBEB5B1B506D7EEDD42509975CCE5309741B38321A96BD9E83C2922335AB72E25DD085FA6FC3769F583D228BF9763083826EE8D3B6C222EEF6A26EDD625935C85DB76C111105F9673A0260D2EB28779D7A5");
        GFInitFromString(G.x,"A9798A3C1EF33A7EEE16ACBF1B82F6F44240DA21E7CEAC3BD619C666FA51EB25B5B4260E095EC74204FEC48FE786F7AC074776D596F56FA049905B51741E9C6688951EF053661C4A4418BF874F04CB0E55CFB958731EC1019C72829776A33DE");
        break;
        case UA_768_5:
        GFInitFromString(p,  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFA83D");
        GFInitFromString(n,  "4000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000005616FD85445E620B96CE7B9971B100DE2945CF2FC4359CF343B1CD6C8C58C87B26689C5528D69695A2B88620AFBF775F");
        GFInitFromString(d,  "75");
        GFInitFromString(G.y,"48001FA5FCC83A5B1012DED30675CE30F3E1D4037F15ED8849D948BADCB52412DFF3A4E436B67CD9CD0A3F58C515770C9C618BB42DA89EA9D3FF885D4802495F28B4B412B0C286314A8D22EF1974526A3DBE4F09E920A64BB195B4430788A246");
        GFInitFromString(G.x,"EE0B664BACF7FD2F1A6403C3D496BFAA51DCD180BB32862FA937C8DFB449997DE6D2EDCEE3490421A159A3B7031D85B2A4958541EAE286A61F0CF5963D98FCE7E2AFF75576E803D9D1222695E233FC467F30D95B3028B4DFB33A41856CE330B");
        break;

        default:
        return -4;
    }
    return EcEdTwistedUAInit(ecc, bitLen, p, &G, n, d, id);
    
}

void EcDump(Ec* ecc, char* buf) {
    char e[60];
    switch (ecc->curve_id)
    {
    case FIPS_192:
        strcpy(buf, "FIPS P-192"); break;
    case FIPS_224:
        strcpy(buf, "FIPS P-224"); break;
    case FIPS_384:
        strcpy(buf, "FIPS P-384"); break;
    case FIPS_521:
        strcpy(buf, "FIPS P-521"); break;
    case ED_192:
        strcpy(buf, "ED P-192"); break;
    case ED_224:
        strcpy(buf, "ED P-224"); break;
    case ED_384:
        strcpy(buf, "ED P-384"); break;
    case ED_521:
        strcpy(buf, "ED P-521"); break;
    case UA_256_1:
        strcpy(buf, "UA P-256/1"); break;
    case UA_256_2:
        strcpy(buf, "UA P-256/2"); break;
    case UA_256_3:
        strcpy(buf, "UA P-256/3"); break;
    case UA_256_4:
        strcpy(buf, "UA P-256/4"); break;
    case UA_256_5:
        strcpy(buf, "UA P-256/5"); break;
    case UA_384_1:
        strcpy(buf, "UA P-384/1"); break;
    case UA_384_2:
        strcpy(buf, "UA P-384/2"); break;
    case UA_384_3:
        strcpy(buf, "UA P-384/3"); break;
    case UA_384_4:
        strcpy(buf, "UA P-384/4"); break;
    case UA_384_5:
        strcpy(buf, "UA P-384/5"); break;
    case UA_512_1:
        strcpy(buf, "UA P-512/1"); break;
    case UA_512_2:
        strcpy(buf, "UA P-512/2"); break;
    case UA_512_3:
        strcpy(buf, "UA P-512/3"); break;
    case UA_512_4:
        strcpy(buf, "UA P-512/4"); break;
    case UA_512_5:
        strcpy(buf, "UA P-512/5"); break;
    case UA_768_1:
        strcpy(buf, "UA P-768/1"); break;
    case UA_768_2:
        strcpy(buf, "UA P-768/2"); break;
    case UA_768_3:
        strcpy(buf, "UA P-768/3"); break;
    case UA_768_4:
        strcpy(buf, "UA P-768/4"); break;
    case UA_768_5:
        strcpy(buf, "UA P-768/5"); break;
    default:
        sprintf(buf, "Unknown curve: ID: %d, IsEdwards: %d, BitLength: %d", ecc->curve_id, ecc->isEdwards, ecc->bitLen);
        break;
    }

}