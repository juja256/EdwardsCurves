
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
    generator->state[0] ^= 1;
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
    return 0;
}

int EcWInit(EcW* ecc, u64 bitLen, const BigInt p, const EcPoint* bp, const BigInt n, const GFElement a, const GFElement b) {  
    ecc->isEdwards = 0;
    BaseEcInit(ecc, bitLen, p, bp, n); 

    ecc->EcAdd = EcWAddProj;
    ecc->EcDouble = EcWDoubleProj;

    copy(ecc->a, a, ecc->wordLen);
    copy(ecc->b, b, ecc->wordLen);
    
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
    GFElement x,y,z;
    GFSqr(ecc, P->x, x);
    GFSqr(ecc, P->y, y);
    GFAdd(ecc, x, y, z);
    GFMul(ecc, x, y, x);
    //GFMul(ecc, x, ecc->d, x);
    GFMulByD(ecc, x);
    GFAdd(ecc, x, unity, x);
    return !GFCmp(ecc, x, z); //1 - ok
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
    //GFMul(ecc, z3, ecc->d, z3); 
    GFMulByD(ecc, z3); // z3 = d * x1 * x2 * y1 * y2

    GFNeg(ecc, z3, z4); // z4 = - z3
    GFAdd(ecc, z3, unity, z3);

    GFInv(ecc, z3, z3);

    GFAdd(ecc, z4, unity, z4);

    GFInv(ecc, z4, z4);

    GFMul(ecc, A->x, B->y, z5); // z5 = x1 * y2
    GFMul(ecc, A->y, B->x, z6); // z6 = x2 * y1
    GFAdd(ecc, z5, z6, z5);
    GFSub(ecc, z2, z1, z2);

    GFMul(ecc, z5, z3, C->x);
    GFMul(ecc, z2, z4, C->y);
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
    //GFMul(ecc, z4, ecc->d, z4);
    GFMulByD(ecc, z4);
    GFNeg(ecc, z4, z5);
    GFAdd(ecc, z4, unity, z4);
    GFInv(ecc, z4, z4);
    GFAdd(ecc, z5, unity, z5);
    GFInv(ecc, z5, z5);
    GFSub(ecc, z2, z1, z2);
    GFMul(ecc, z4, z3, B->x);
    GFMul(ecc, z2, z5, B->y);
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
    int y = 0;
    GFElement t1, t2, t3, pp;
    copy(pp, ecc->p, ecc->wordLen);
    div2(ecc->wordLen, pp);

    do {
        GFSqr(ecc, bp->x, t1);
        GFNeg(ecc, t1, t2);
        GFAdd(ecc, t2, unity, t2);
        //GFMul(ecc, t1, ecc->d, t1);
        GFMulByD(ecc, t1);
        GFNeg(ecc, t1, t1);
        GFAdd(ecc, t1, unity, t1);
        GFInv(ecc, t1, t1);
        GFMul(ecc, t1, t2, t1);
        GFPow(ecc, t1, pp, t2);
        y = !GFCmp(ecc, t2, unity);
        if (!y)
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

/*void EcEdScalarMul(EcEd* ecc, const EcPoint* A, const BigInt k, EcPoint* B) {
    EcPoint P, H;
    EcCopy(ecc, &P, &uPEd);
    EcCopy(ecc, &H, A);

    for (u32 i=0; i<ecc->bitLen; i++) {
        if (get_bit(k, i)) {
            EcEdAdd(ecc, &P, &H, &P);
        }
        EcEdDouble(ecc, &H, &H); 
    }
    EcCopy(ecc, B, &P);
}

int EcWScalarMul(EcW* ecc, const EcPoint* A, const BigInt k, EcPoint* B)
{
    int s = NORMAL_POINT;
    
    int hb = bigint_bit_len(ecc->wordLen, k) - 1;

    EcPoint P;
    EcCopy(ecc, &P, A);
    for (int i = hb - 1;i >= 0;i--)
    {
        s &= EcWDouble(ecc, &P, &P); 
        if (get_bit(k, i))
            s &= EcWAdd(ecc, &P, A, &P);
    }
    EcCopy(ecc, B, &P);
    return s;
}

int EcScalarMul(Ec* ecc, const EcPoint* A, const BigInt k, EcPoint* B) {
    if (ecc->isEdwards) { EcEdScalarMul(ecc, A, k, B); return NORMAL_POINT; }
    else return EcWScalarMul(ecc, A, k, B);
}*/

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

void EcEdAddProj(EcEd* ecc, const EcPointProj* A, const EcPointProj* B, EcPointProj* C) {
    GFElement a,b,c,d,e,e1,e2,e3,f;
    GFMul(ecc, A->Z, B->Z, a);
    GFSqr(ecc, a, b);
    GFMul(ecc, A->X, B->X, c);
    GFMul(ecc, A->Y, B->Y, d);
    GFAdd(ecc, A->X, A->Y, e);
    GFAdd(ecc, B->X, B->Y, e1);
    GFMul(ecc, e, e1, e);
    GFSub(ecc, e, c, e);
    GFSub(ecc, e, d, e);
    GFMul(ecc, c, d, f);

    GFMulByD(ecc, f);
    GFSub(ecc, b, f, e2);
    GFMul(ecc, e, a, C->X);
    GFMul(ecc, C->X, e2, C->X);

    GFAdd(ecc, b, f, e3);
    GFSub(ecc, d, c, C->Y);
    GFMul(ecc, C->Y, e3, C->Y);
    GFMul(ecc, C->Y, a, C->Y);

    GFMul(ecc, e2, e3, C->Z);
}

void EcWAddProj(EcW* ecc, const EcPointProj* A, const EcPointProj* B, EcPointProj* C)
{
    if ( GFCmp(ecc, A->Z, &zero) == 0 ) { // Handling infinities
        EcCopyProj(ecc, C, B);
        return;
    }
    if ( GFCmp(ecc, B->Z, &zero) == 0 ) { // Handling infinities
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

void EcEdDoubleProj(EcEd* ecc, const EcPointProj* A, EcPointProj* B) {
    GFElement a,b,c,d,e,e1,e2,e3,f;
    GFSqr(ecc, A->Z, a);
    GFSqr(ecc, a, b);
    GFSqr(ecc, A->X, c);
    GFSqr(ecc, A->Y, d);

    GFMul(ecc, A->X, A->Y, e);
    GFMulBy2(ecc, e, e);

    GFMul(ecc, c, d, f);
    //GFMul(ecc, f, ecc->d, f);
    GFMulByD(ecc, f);
    GFSub(ecc, b, f, e2);
    GFMul(ecc, e, a, B->X);
    GFMul(ecc, B->X, e2, B->X);

    GFAdd(ecc, b, f, e3);
    GFSub(ecc, d, c, B->Y);
    GFMul(ecc, B->Y, e3, B->Y);
    GFMul(ecc, B->Y, a, B->Y);

    GFMul(ecc, e2, e3, B->Z);
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
    GFMulBy2(ecc,b,tmp);
    GFMulBy2(ecc,tmp,tmp);
    GFMulBy2(ecc,tmp,tmp);//now tmp = 8b
    GFSub(ecc,h,tmp,h);
    //x2 = 2hs
    GFMul(ecc,h,s,B->X);
    GFMulBy2(ecc,B->X,B->X);
    //y2 = w(4b - h) - 8(y1*s)^2
    GFMul(ecc,A->Y,s,tmp);
    GFMulBy2(ecc,tmp,tmp);
    GFSqr(ecc,tmp,tmp);
    GFMulBy2(ecc,tmp,tmp); // now tmp = 8*(y1*s)^2
    GFMulBy2(ecc,b,B->Y);
    GFMulBy2(ecc,B->Y,B->Y);
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

static inline void EcScalarMulNaive(Ec* ecc, const EcPointProj* A, const BigInt k, EcPointProj* B) {
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

static inline void EcScalarMulMontgomery(Ec* ecc, const EcPointProj* A, const BigInt k, EcPointProj* B) {
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

void EcScalarMulProj(Ec* ecc, const EcPointProj* A, const BigInt k, EcPointProj* B) {
    EcScalarMulMontgomery(ecc, A, k, B);
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
            GFInitFromString(G.x,"44F083BB00E51AD91A2743284D31F57EE5C84826FCC91F4B");
            GFInitFromString(G.y,"15FC16E5870524E0DBBE9EC8BB9F066C02A02B1978D4E029");
            break;
            case 224:
            GFInitFromString(p,  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000000000000000000001");
            GFInitFromString(n,  "400000000000000000000000000020BBEC47CEDB34DD05BCB6B7E619");
            GFInitFromString(d,  "3608425");
            GFInitFromString(G.x,"C448CA02660F57204FF1BDE2B5CC3E25606A7460399FEA3DA9A06383");
            GFInitFromString(G.y,"319117770D6FC7FE35F6A02905FE1F363156BD2E5B75BB89A64CAFAB");
            break;
            case 256:
            GFInitFromString(p,  "FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF");
            GFInitFromString(n,  "3FFFFFFFC00000003FFFFFFFFFFFFFFFBA76FA29C30CC3AA4954B53EDBE54D75");
            GFInitFromString(d,  "72A38");
            GFInitFromString(G.x,"894F8283626AEE6848515DDDC3B8DBB3D5302DEE0EE75080D6753E4D39BA5AB2");
            GFInitFromString(G.y,"EA612346223F6480CBBAFA39DB95D54D21469DD3074A957EFDA4FD79FEB630B5");
            break;
            case 384:
            GFInitFromString(p,  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFF0000000000000000FFFFFFFF");
            GFInitFromString(n,  "4000000000000000000000000000000000000000000000005063576B5A9A0C3A23E9510EA680650B4884E63A2968DD71");
            GFInitFromString(d,  "12593A");
            GFInitFromString(G.x,"1FC0E8E61F599813E376D11F7510D77F177C2F1CDE19FD14D63A2EC5EAD4D0DED1BD06676CCF365243BF3C0675A31B62");
            GFInitFromString(G.y,"F52B4FA352B257D7A102FA45C56A50CCBDB3DEC053D5610EDBD0188C11F321F28A43E2FC50395E4A8BD0029AE71D51AA");
            break;
            case 521:
            GFInitFromString(p,  "1FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF");
            GFInitFromString(n,  "7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF46087BC4A294FCC80B3F45D8CEBFB21479BA651BA07DE913AD1D8392DE3FF8AF");
            GFInitFromString(d,  "16A");
            GFInitFromString(G.x,"1C4A8AA03837A73FC86D95F308BD9B738E2BFAB5B8ECCEB0EF90F02B5E90D1DF28D470C4F531212CD0E68F4E925E5ED82A74AB63335A9DB3E31650C6767EBA00681");
            GFInitFromString(G.y,"11365BAC31463E811ECB7FEF04F65DE5B8B15ADAFD72CC74EC804B9408EB31F0E44CE67EA36144A06B07D105B862E2493B556740F343DB87866E15C7BE9F3813666");
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

void EcDump(Ec* ecc, char* buf) {
    char e[60];
    if (ecc->isEdwards) strcpy(e, "Edwards"); else strcpy(e, "Weierstrass");

    int l = sprintf(buf, "P-%d(%s)", ecc->bitLen, e);
}