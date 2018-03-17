
#include "gf.h"
#include "ec.h"
#include <stdlib.h>
#include <time.h>
#include <string.h>


int EcEdInit(EcEd* ecc, const EcPoint* bp, u64 bitLen, const BigInt n, const GFElement d) {
    srand(time(NULL));

    if ( (bitLen != 192) && (bitLen != 224) && (bitLen != 256) && (bitLen != 384)) {
        return -1;
    }
    ecc->isEdwards = 1;
    ecc->bitLen = bitLen;
    ecc->wordLen = (bitLen%64 == 0) ? (bitLen / 64) : (bitLen / 64) + 1;
    copy(ecc->d, d, ecc->wordLen);
    copy(ecc->n, n, ecc->wordLen);
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

    BigInt two;
    memset(two, 0, ecc->wordLen * 8); two[0] = 2;
    sub(ecc->wordLen, ecc->p, two, ecc->p_min_two);
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

int EcCheckPointOnCurve(const Ec* ecc, const EcPoint* P) {
    if (ecc->isEdwards) return EcEdCheckPointOnCurve(ecc, P);
    return -1;
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
}

void EcGenerateBasePoint(const Ec* ecc, EcPoint* bp) {
    if (ecc->isEdwards) EcEdGenerateBasePoint(ecc, bp);
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

void EcAdd(const Ec* ecc, const EcPoint* A, const EcPoint* B, EcPoint* C) {
    if (ecc->isEdwards) EcEdAdd(ecc, A, B, C);
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

void EcDouble(const Ec* ecc, const EcPoint* A, EcPoint* B) {
    if (ecc->isEdwards) EcEdDouble(ecc, A, B);
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

void EcScalarMul(const Ec* ecc, const EcPoint* A, const BigInt k, EcPoint* B) {
    if (ecc->isEdwards) EcEdScalarMul(ecc, A, k, B);
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

void EcAddProj(const Ec* ecc, const EcPointProj* A, const EcPointProj* B, EcPointProj* C) {
    if (ecc->isEdwards) EcEdAddProj(ecc, A, B, C);
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

void EcDoubleProj(const Ec* ecc, const EcPointProj* A, EcPointProj* B) {
    if (ecc->isEdwards) EcEdDoubleProj(ecc, A, B);
}

void EcEdScalarMulProj(const EcEd* ecc, const EcPoint* A, const BigInt k, EcPoint* B) {
    EcPointProj P, H;
    EcConvertAffineToProjective(ecc, &uP, &P);
    EcConvertAffineToProjective(ecc, A, &H);

    for (u32 i=0; i<ecc->bitLen; i++) {
        if (get_bit(k, i)) {
            EcEdAddProj(ecc, &P, &H, &P);
        }
        EcEdDoubleProj(ecc, &H, &H); 
    }
    EcConvertProjectiveToAffine(ecc, &P, B);
}

void EcScalarMulProj(const Ec* ecc, const EcPoint* A, const BigInt k, EcPoint* B) {
    if (ecc->isEdwards) EcEdScalarMulProj(ecc, A, k, B);
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
*/
int EcInitStandardCurve(Ec* ecc, u64 bitLen, BOOL isEdwards) {
    if (isEdwards) {
        GFElement p,d,n;
        EcPoint G;
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
    return -2;
}