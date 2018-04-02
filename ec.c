
#include "gf.h"
#include "ec.h"
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <stdio.h>

int EcPointCmp(const Ec* ecc, const EcPoint* A, const EcPoint* B) {
    return GFCmp(ecc, A->x, B->x) || GFCmp(ecc, A->y, B->y);
}

void EcCopy(const Ec* ecc, EcPoint* dest, const EcPoint* src) {
    copy(dest->x, src->x, ecc->wordLen);
    copy(dest->y, src->y, ecc->wordLen);
}

void EcCopyProj(const Ec* ecc, EcPointProj* dest, const EcPointProj* src) {
    copy(dest->X, src->X, ecc->wordLen);
    copy(dest->Y, src->Y, ecc->wordLen);
    copy(dest->Z, src->Z, ecc->wordLen);
}

void BaseEcInit(Ec* ecc,const EcPoint* bp,u64 bitLen, const BigInt n)
{
    srand(time(NULL));
    ecc->bitLen = bitLen;
    ecc->wordLen = (bitLen%64 == 0) ? (bitLen / 64) : (bitLen / 64) + 1;
    copy(ecc->n, n, ecc->wordLen);
    ecc->n[ecc->wordLen] = 0;
    copy(ecc->BasePoint.x, bp->x, ecc->wordLen);
    copy(ecc->BasePoint.y, bp->y, ecc->wordLen);

    switch (bitLen) {
        case 192:
            copy(ecc->p, p192, ecc->wordLen);
            ecc->GFMul = GFMul_FIPS192;
            ecc->GFSqr = GFSqr_FIPS192;
            break;
        case 224:
            copy(ecc->p, p224, ecc->wordLen);
            ecc->GFMul = GFMul_FIPS224;
            ecc->GFSqr = GFSqr_FIPS224;
            break;

        case 256:
            copy(ecc->p, p256, ecc->wordLen);
            ecc->GFMul = GFMul_FIPS256;
            ecc->GFSqr = GFSqr_FIPS256;
            break;

        case 384: 
            copy(ecc->p, p384, ecc->wordLen);
            ecc->GFMul = GFMul_FIPS384;
            ecc->GFSqr = GFSqr_FIPS384;
            break;
    }
    ecc->p[ecc->wordLen] = 0;

    BigInt two;
    memset(two, 0, ecc->wordLen * 8); two[0] = 2;
    sub(ecc->wordLen, ecc->p, two, ecc->p_min_two);

    unsigned char seed[4];
    srand(time(NULL));
    seed[0] = rand();
    PRNGInit(&(ecc->prng), seed, 4);
}

int EcEdInit(EcEd* ecc, const EcPoint* bp, u64 bitLen, const BigInt n, const GFElement d) {
    if ( (bitLen != 192) && (bitLen != 224) && (bitLen != 256) && (bitLen != 384)) {
        return -1;
    }
    BaseEcInit(ecc,bp,bitLen,n);
    ecc->isEdwards = 1;
    ecc->cofactor = 4; // Cofactor = 4 for all full Edwards curves
    copy(ecc->d, d, ecc->wordLen);
    return 0;
}

int EcWInit(EcW* ecc, const EcPoint* bp, u64 bitLen, const BigInt n, const GFElement a,const GFElement b)
{   
    if ( (bitLen != 192) && (bitLen != 224) && (bitLen != 256) && (bitLen != 384)) {
        return -1;
    }
    BaseEcInit(ecc,bp,bitLen,n);
    ecc->isEdwards = 0;
    ecc->cofactor = 1; // Cofactor = 1 for all NIST Recommended curves with a = -3
    copy(ecc->a, a, ecc->wordLen);
    copy(ecc->b, b, ecc->wordLen);
    return 0;
}

int EcEdCheckPointOnCurve(const EcEd* ecc, const EcPoint* P) {
    GFElement x,y,z;
    GFSqr(ecc, P->x, x);
    GFSqr(ecc, P->y, y);
    GFAdd(ecc, x, y, z);
    GFMul(ecc, x, y, x);
    GFMul(ecc, x, ecc->d, x);
    GFAdd(ecc, x, unity, x);
    return !GFCmp(ecc, x, z); //1 - ok
}

int EcWCheckPointOnCurve(const EcW* ecc, const EcPoint* P) {
    GFElement l,r; //left,right
    GFSqr(ecc, P->y, l);
    GFSqr(ecc, P->x, r);
    GFAdd(ecc, r, ecc->a, r);
    GFMul(ecc, r, P->x, r); // now r = x(x^2 + a)
    GFAdd(ecc, r, ecc->b, r);
    return !GFCmp(ecc, l, r); //1 - ok
}

int EcCheckPointOnCurve(const Ec* ecc, const EcPoint* P) {
    if (ecc->isEdwards) return EcEdCheckPointOnCurve(ecc, P);
    else return EcWCheckPointOnCurve(ecc, P);
}

int EcCheckPointInMainSubGroup(const Ec* ecc, const EcPoint* P) {
    int r = EcCheckPointOnCurve(ecc, P);
    if (!r) return 0;
    EcPoint Q;
    int s = EcScalarMul(ecc, P, ecc->n, &Q);
    if ( ecc->isEdwards && (EcPointCmp(ecc, &Q, &uP) == 0) ) return 1;
    else if ( !(ecc->isEdwards) && (s == 0) ) return 1;
    return 0;
}


/* x^2 + y^2 = 1 + dx^2y^2 */
/* y^2 = (1 - x^2)/(1 - dx^2) */
void EcEdGenerateBasePoint(const EcEd* ecc, EcPoint* bp) {
    randomize(ecc->wordLen, bp->x);
    if (ecc->bitLen % 64 == 32) { // P-224
        bp->x[ecc->wordLen-1] &= 0xFFFFFFFF;
    }
    int y = 0;
    GFElement t1, t2, t3, pp;
    copy(pp, ecc->p, ecc->wordLen);
    div2(ecc->wordLen, pp);

    do {
        GFSqr(ecc, bp->x, t1);
        GFNeg(ecc, t1, t2);
        GFAdd(ecc, t2, unity, t2);
        GFMul(ecc, t1, ecc->d, t1);
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

    tonelli_shanks_sqrt(ecc, t1, bp->y);
    EcDouble(ecc, bp, bp);
    EcDouble(ecc, bp, bp); /* Keeping Point in Main SubGroup by multiplication by factor (x4) */
}

//y^2 = x^3 + ax + b
void EcWGenerateBasePoint(const EcW* ecc, EcPoint* bp) {
    randomize(ecc->wordLen, bp->x);
    if (ecc->bitLen % 64 == 32) { // P-224
        bp->x[ecc->wordLen-1] &= 0xFFFFFFFF;
    }
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

    tonelli_shanks_sqrt(ecc, sum, bp->y);
}

void EcGenerateBasePoint(const Ec* ecc, EcPoint* bp) {
    if (ecc->isEdwards) EcEdGenerateBasePoint(ecc, bp);
    else EcWGenerateBasePoint(ecc,bp);
}

void EcEdAdd(const EcEd* ecc, const EcPoint* A, const EcPoint* B, EcPoint* C) {
    GFElement z1, z2, z3, z4, z5, z6, z7;
    GFMul(ecc, A->x, B->x, z1); // z1 = x1 * x2
    GFMul(ecc, A->y, B->y, z2); // z2 = y1 * y2
    GFMul(ecc, z1, z2, z3); 
    GFMul(ecc, z3, ecc->d, z3); // z3 = d * x1 * x2 * y1 * y2
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

int EcWAdd(const EcW* ecc, const EcPoint* A, const EcPoint* B, EcPoint* C)
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

int EcAdd(const Ec* ecc, const EcPoint* A, const EcPoint* B, EcPoint* C) {
    if (ecc->isEdwards) { EcEdAdd(ecc, A, B, C); return NORMAL_POINT; }
    else return EcWAdd(ecc, A, B, C);
}

void EcEdDouble(const EcEd* ecc, const EcPoint* A, EcPoint* B) {
    GFElement z1, z2, z3, z4, z5;
    GFSqr(ecc, A->x, z1);
    GFSqr(ecc, A->y, z2);
    GFMul(ecc, A->x, A->y, z3);
    GFMulBy2(ecc, z3, z3);
    GFMul(ecc, z1, z2, z4);
    GFMul(ecc, z4, ecc->d, z4);
    GFNeg(ecc, z4, z5);
    GFAdd(ecc, z4, unity, z4);
    GFInv(ecc, z4, z4);
    GFAdd(ecc, z5, unity, z5);
    GFInv(ecc, z5, z5);
    GFSub(ecc, z2, z1, z2);
    GFMul(ecc, z4, z3, B->x);
    GFMul(ecc, z2, z5, B->y);
}

int EcWDouble(const EcW* ecc, const EcPoint* A, EcPoint* B)
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

int EcDouble(const Ec* ecc, const EcPoint* A, EcPoint* B) {
    if (ecc->isEdwards) { EcEdDouble(ecc, A, B); return NORMAL_POINT; }
    else return EcWDouble(ecc,A,B);
}

void EcEdScalarMul(const EcEd* ecc, const EcPoint* A, const BigInt k, EcPoint* B) {
    EcPoint P, H;
    copy_point(&P, &uP, ecc->wordLen);
    copy_point(&H, A,  ecc->wordLen);

    for (u32 i=0; i<ecc->bitLen; i++) {
        if (get_bit(k, i)) {
            EcEdAdd(ecc, &P, &H, &P);
        }
        EcEdDouble(ecc, &H, &H); 
    }
    copy_point(B, &P, ecc->wordLen);
}

int EcWScalarMul(const EcW* ecc, const EcPoint* A, const BigInt k, EcPoint* B)
{
    int s = NORMAL_POINT;
    
    int hb = bigint_bit_len(ecc->wordLen, k) - 1;

    EcPoint P;
    copy_point(&P, A,  ecc->wordLen);
    for (int i = hb - 1;i >= 0;i--)
    {
        s &= EcWDouble(ecc, &P, &P); 
        if (get_bit(k, i))
            s &= EcWAdd(ecc, &P, A, &P);
    }
    copy_point(B, &P, ecc->wordLen);
    return s;
}

int EcScalarMul(const Ec* ecc, const EcPoint* A, const BigInt k, EcPoint* B) {
    if (ecc->isEdwards) { EcEdScalarMul(ecc, A, k, B); return NORMAL_POINT; }
    else return EcWScalarMul(ecc, A, k, B);
}

void EcConvertAffineToProjective(const Ec* ecc, const EcPoint* P, EcPointProj* Q) {
    randomize(ecc->wordLen, Q->Z);
    GFMul(ecc, P->x, Q->Z, Q->X);
    GFMul(ecc, P->y, Q->Z, Q->Y);
}

void EcConvertProjectiveToAffine(const Ec* ecc, const EcPointProj* P, EcPoint* Q) {
    GFElement Z_inv;
    GFInv(ecc, P->Z, Z_inv);
    GFMul(ecc, P->X, Z_inv, Q->x);
    GFMul(ecc, P->Y, Z_inv, Q->y);
}

void EcEdAddProj(const EcEd* ecc, const EcPointProj* A, const EcPointProj* B, EcPointProj* C) {
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
    GFMul(ecc, f, ecc->d, f);

    GFSub(ecc, b, f, e2);
    GFMul(ecc, e, a, C->X);
    GFMul(ecc, C->X, e2, C->X);

    GFAdd(ecc, b, f, e3);
    GFSub(ecc, d, c, C->Y);
    GFMul(ecc, C->Y, e3, C->Y);
    GFMul(ecc, C->Y, a, C->Y);

    GFMul(ecc, e2, e3, C->Z);
}

int EcWAddProj(const EcW* ecc, const EcPointProj* A, const EcPointProj* B, EcPointProj* C)
{
    GFElement u,v,a,tmp,x1z2,sqrv,y1z2;
    //v = x2z1 - x1z2
    GFMul(ecc,A->X,B->Z,x1z2); 
    GFMul(ecc,A->Z,B->X,v);
    if(!GFCmp(ecc,v,x1z2)) //=> z3 = 0, to infinity and beyond
        return INFINITY_POINT;
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
    return NORMAL_POINT;
}

int EcAddProj(const Ec* ecc, const EcPointProj* A, const EcPointProj* B, EcPointProj* C) {
    if (ecc->isEdwards) {EcEdAddProj(ecc, A, B, C); return NORMAL_POINT;}
    else return EcWAddProj(ecc, A, B, C);
}

void EcEdDoubleProj(const EcEd* ecc, const EcPointProj* A, EcPointProj* B) {
    GFElement a,b,c,d,e,e1,e2,e3,f;
    GFSqr(ecc, A->Z, a);
    GFSqr(ecc, a, b);
    GFSqr(ecc, A->X, c);
    GFSqr(ecc, A->Y, d);

    GFMul(ecc, A->X, A->Y, e);
    GFMulBy2(ecc, e, e);

    GFMul(ecc, c, d, f);
    GFMul(ecc, f, ecc->d, f);

    GFSub(ecc, b, f, e2);
    GFMul(ecc, e, a, B->X);
    GFMul(ecc, B->X, e2, B->X);

    GFAdd(ecc, b, f, e3);
    GFSub(ecc, d, c, B->Y);
    GFMul(ecc, B->Y, e3, B->Y);
    GFMul(ecc, B->Y, a, B->Y);

    GFMul(ecc, e2, e3, B->Z);
}

int EcWDoubleProj(const EcW* ecc, const EcPointProj* A, EcPointProj* B)
{
    if(!GFCmp(ecc,A->Y,zero))
        return INFINITY_POINT;
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
    return NORMAL_POINT;
}

int EcDoubleProj(const Ec* ecc, const EcPointProj* A, EcPointProj* B) {
    if (ecc->isEdwards) {EcEdDoubleProj(ecc, A, B); return NORMAL_POINT;}
    else return EcWDoubleProj(ecc, A, B);
}

void EcEdScalarMulProj(const EcEd* ecc, const EcPointProj* A, const BigInt k, EcPointProj* B) {
    EcPointProj P, H;
    EcConvertAffineToProjective(ecc, &uP, &P);
    
    EcCopyProj(ecc, &H, A); // H := A
    EcCopyProj(ecc, B, &P); // B := O

    for (u32 i=0; i<ecc->bitLen; i++) {
        if (get_bit(k, i)) {
            EcEdAddProj(ecc, B, &H, B);
        }
        EcEdDoubleProj(ecc, &H, &H); 
    }
}

int EcWScalarMulProj(const EcW* ecc, const EcPointProj* A, const BigInt k, EcPointProj* B) {
    EcPointProj H;
    
    EcCopyProj(ecc, B, A); // B := A
    EcCopyProj(ecc, &H, A); // H := A

    int s = NORMAL_POINT;
    
    int hb = bigint_bit_len(ecc->wordLen, k) - 1;

    for (int i = hb - 1;i >= 0; i--)
    {
        s &= EcWDoubleProj(ecc, B, B); 
        if (get_bit(k, i))
            s &= EcWAddProj(ecc, B, &H, B);
    }
    return s;
}

int EcScalarMulProj(const Ec* ecc, const EcPointProj* A, const BigInt k, EcPointProj* B) {
    if (ecc->isEdwards) {EcEdScalarMulProj(ecc, A, k, B); return NORMAL_POINT;}
    else return EcWScalarMulProj(ecc, A, k, B);
}



/*
    Standard Edwards full curves:
        P-192:
    p = FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFF
    d = 6DBA6A
    n = 3FFFFFFFFFFFFFFFFFFFFFFFEA75D4027230DD4DFFDB0455
    x = 44F083BB00E51AD91A2743284D31F57EE5C84826FCC91F4B
    y = 15FC16E5870524E0DBBE9EC8BB9F066C02A02B1978D4E029

        P-224
    p = FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000000000000000000001
    d = 3608425
    n = 400000000000000000000000000020BBEC47CEDB34DD05BCB6B7E619
    x = C448CA02660F57204FF1BDE2B5CC3E25606A7460399FEA3DA9A06383
    y = 319117770D6FC7FE35F6A02905FE1F363156BD2E5B75BB89A64CAFAB

        P-256
    p = FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF
    d = 72A38
    n = 3FFFFFFFC00000003FFFFFFFFFFFFFFFBA76FA29C30CC3AA4954B53EDBE54D75
    x = 894F8283626AEE6848515DDDC3B8DBB3D5302DEE0EE75080D6753E4D39BA5AB2
    y = EA612346223F6480CBBAFA39DB95D54D21469DD3074A957EFDA4FD79FEB630B5

        P-384
    p = FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFF0000000000000000FFFFFFFF
    d = 12593A
    n = 4000000000000000000000000000000000000000000000005063576B5A9A0C3A23E9510EA680650B4884E63A2968DD71
    x = 1FC0E8E61F599813E376D11F7510D77F177C2F1CDE19FD14D63A2EC5EAD4D0DED1BD06676CCF365243BF3C0675A31B62
    y = F52B4FA352B257D7A102FA45C56A50CCBDB3DEC053D5610EDBD0188C11F321F28A43E2FC50395E4A8BD0029AE71D51AA

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

*/
int EcInitStandardCurve(Ec* ecc, u64 bitLen, BOOL isEdwards) {
    GFElement p,n;
    EcPoint G;
    if (isEdwards) {
        GFElement d;
        switch (bitLen) {
            case 192:
            GFInitFromString(p, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFF");
            GFInitFromString(n, "3FFFFFFFFFFFFFFFFFFFFFFFEA75D4027230DD4DFFDB0455");
            GFInitFromString(d, "6DBA6A");
            GFInitFromString(G.x, "44F083BB00E51AD91A2743284D31F57EE5C84826FCC91F4B");
            GFInitFromString(G.y, "15FC16E5870524E0DBBE9EC8BB9F066C02A02B1978D4E029");
            break;
            case 224:
            GFInitFromString(p, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000000000000000000001");
            GFInitFromString(n, "400000000000000000000000000020BBEC47CEDB34DD05BCB6B7E619");
            GFInitFromString(d, "3608425");
            GFInitFromString(G.x, "C448CA02660F57204FF1BDE2B5CC3E25606A7460399FEA3DA9A06383");
            GFInitFromString(G.y, "319117770D6FC7FE35F6A02905FE1F363156BD2E5B75BB89A64CAFAB");
            break;
            case 256:
            GFInitFromString(p, "FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF");
            GFInitFromString(n, "3FFFFFFFC00000003FFFFFFFFFFFFFFFBA76FA29C30CC3AA4954B53EDBE54D75");
            GFInitFromString(d, "72A38");
            GFInitFromString(G.x, "894F8283626AEE6848515DDDC3B8DBB3D5302DEE0EE75080D6753E4D39BA5AB2");
            GFInitFromString(G.y, "EA612346223F6480CBBAFA39DB95D54D21469DD3074A957EFDA4FD79FEB630B5");
            break;
            case 384:
            GFInitFromString(p, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFF0000000000000000FFFFFFFF");
            GFInitFromString(n, "4000000000000000000000000000000000000000000000005063576B5A9A0C3A23E9510EA680650B4884E63A2968DD71");
            GFInitFromString(d, "12593A");
            GFInitFromString(G.x, "1FC0E8E61F599813E376D11F7510D77F177C2F1CDE19FD14D63A2EC5EAD4D0DED1BD06676CCF365243BF3C0675A31B62");
            GFInitFromString(G.y, "F52B4FA352B257D7A102FA45C56A50CCBDB3DEC053D5610EDBD0188C11F321F28A43E2FC50395E4A8BD0029AE71D51AA");
            break;
            default:
            return -1;
        }
        return EcEdInit(ecc, &G, bitLen, n, d);
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
            default:
            return -1;
        }
        return EcWInit(ecc, &G, bitLen, n, a, b);
    }
}

void EcDump(const Ec* ecc, char* buf) {
    char e[60];
    if (ecc->isEdwards) strcpy(e, "Edwards"); else strcpy(e, "Weierstrass");

    int l = sprintf(buf, "P-%d(%s)", ecc->bitLen, e);
}