#include "eced.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

#define MSB_M 0x8000000000000000
#define HEX_FORMAT "%.16llX"
static u64 unity[] = {0x1, 0,0,0,0,0};

static EcPoint uP = { {0,0,0,0,0,0}, {1,0,0,0,0,0} };

static u64 p192[] = { 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFE, 0xFFFFFFFFFFFFFFFF };

static u64 p224[] = { 1, 0xFFFFFFFF00000000, 0xFFFFFFFFFFFFFFFF, 0x00000000FFFFFFFF };

static u64 p256[] = { 0xFFFFFFFFFFFFFFFF, 0x00000000FFFFFFFF, 0x0000000000000000, 0xFFFFFFFF00000001 };

static u64 p386[] = { 0x00000000FFFFFFFF, 0xFFFFFFFF00000000, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF };

static inline u64 get_bit(GFElement a, u64 num) {
    return a[num/64] & (1UL << (num % 64));
}

static inline void copy(GFElement a, GFElement b, int len) {
    for (u64 i=0;i<len;i++)
        a[i] = b[i];
}

static inline void copy_point(EcPoint* P, EcPoint* Q, int len) {
    copy(P->x, Q->x, len);
    copy(P->y, Q->y, len);
}

static inline u64 add(u64 n, u64* a, u64* b, u64* c) {
    u64 msb_a, msb_b, carry = 0;

    for (u64 i=0;i < n; i++) {
        msb_a = a[i] & MSB_M;
        msb_b = b[i] & MSB_M;
        c[i] = a[i] + b[i] + carry;
        carry = ( (msb_a && msb_b) || ((msb_a ^ msb_b) && !(MSB_M & c[i])) );
    }
    return carry;
}

static inline int inc(u64 len, GFElement n) {
    for (u32 i=0;i<len;i++) {
        if (n[i] == 0xFFFFFFFFFFFFFFFF) {
            n[i]++;
            continue;
        }
        else {
            n[i]++;
            return 0;
        } 
    }
    return 1;
}

static inline u64 sub(u64 n, u64* a, u64* b, u64* c) {
    u64 borrow = 0;

    for (u32 i=0; i<n; i++) {
        if ((a[i] >= (b[i] + borrow)) && !(b[i] == 0xFFFFFFFFFFFFFFFF)) {
            c[i] = a[i] - b[i] - borrow;
            borrow = 0;
        }
        else {
            c[i] = (0xFFFFFFFFFFFFFFFF + (a[i] - b[i] - borrow)) + 1;
            borrow = 1;
        }
    }
    return borrow;
}

static inline void _mul_raw(unsigned long long a, unsigned long long b, unsigned long long* low, unsigned long long* high) {
    __int128 r = (__int128)a * (__int128)b;
    *low = (unsigned long long)r;
    *high = r >> 64;
}

static inline u64 _add_raw(unsigned long long a, unsigned long long b, unsigned long long* c) {
    __int128 r = (__int128)a + (__int128)b;
    *c = (u64)r;
    return r >> 64;
}

static inline void mul_by_word(EcEd* ecc, u64* a, u64 d, u64* c) {
    u64 carry = 0, carry_tmp;
    for (u64 i=0; i < ecc->wordLen; i++) {
        carry_tmp = carry;
        _mul_raw(d, a[i], &(c[i]), &carry);
        carry += _add_raw(c[i], carry_tmp, &(c[i]));
    }
    c[ecc->wordLen] = carry;
}

static inline void mul(EcEd* ecc, u64* a, u64* b, u64* c) {
    GFElement tmp;
    memset(c, 0, 2*8*ecc->wordLen);
    for (u64 i=0; i < ecc->wordLen; i++) {
        mul_by_word(ecc, a, b[i], (u64*)tmp);
        add(ecc->wordLen + 1, c+i, tmp, c+i);
    }
}

static inline u64 word_bit_len(u64 n) {
    u32 c = 64;
    while (c) {
        if ((((u64)1 << (64-1)) & n) >> (64-1))
            return c;
        n <<= 1;
        --c;
    }
    return 0;
}

static inline void mul2(u64 n, GFElement a) {
    u64 buf = 0;
    for (u32 i = 0; i < n; i++) {
        u64 cur = a[i];
        a[i] <<= 1;
        a[i] ^= buf;
        buf = (cur & (0xFFFFFFFFFFFFFFFF << ( 63 ))) >> ( 63 );
    }
    a[n] = buf;
}


static inline void div2(u64 n, GFElement a) {
    u64 buf = 0;
    for (int i = n-1; i >= 0; i--) {
        u64 cur = a[i];
        a[i] >>= 1;
        a[i] ^= buf;
        buf = (cur & (u64)1) << 63;
    }
}

void GFInitFromString(GFElement a, const char* str) {
    u64 s_len = strlen(str);
    u64 tmp;

    memset(a, 0, 8*6);

    for (int i = s_len-1; i >= 0; i--) {
        if ((str[i] >= '0') && (str[i] <= '9')) { // 0,1,2,...
            tmp = str[i] - 48;
        }
        else if ((str[i] >= 'A') && (str[i] <= 'F')) { // A,B,...
            tmp = str[i] - 55;
        }
        else {
            return;
        }
        a[(s_len - 1 - i) / 16] ^= ( tmp << (((s_len - 1 - i) % 16)*4) );
    }
}

void GFDump(EcEd* ecc, GFElement a) {
    for (int i=ecc->wordLen-1; i>=0; i--)
        printf(HEX_FORMAT, a[i]);
    printf("\n");
}

/* FIPS-192 Fp: p = 2^192 - 2^64 - 1*/

void GFNeg(EcEd* ecc, GFElement a, _out_ GFElement c) {
    sub(ecc->wordLen, ecc->p, (u64*)a, (u64*)c);
}

void GFAdd(EcEd* ecc, GFElement a, GFElement b, _out_ GFElement c) {
    u64 carry = 0;
    if (ecc->bitLen % 64 == 0)
        carry = add(ecc->wordLen, a,b,c);
    else {
        add(ecc->wordLen, a,b,c);
        carry = (c[ ecc->wordLen -1 ] & (1UL<<32) );
    }
    if (carry) {
        sub(ecc->wordLen, c, ecc->p, c);
    }
    if ((c[2] == 0xFFFFFFFFFFFFFFFF) && (c[1] == 0xFFFFFFFFFFFFFFFF)) {
        sub(ecc->wordLen, c, ecc->p, c);
    }
}

void GFSub(EcEd* ecc, GFElement a, GFElement b, _out_ GFElement c) {
    GFElement tmp;
    copy(tmp, a, ecc->wordLen);
    sub(ecc->wordLen, ecc->p, (u64*)b, (u64*)c);
    GFAdd(ecc, tmp, c, c);
}

void GFMul_FIPS192(EcEd* ecc, GFElement a, GFElement b, _out_ GFElement c) {
    GFElement tmp;
    mul(ecc, a, b, tmp);
    /*for (int i=2*ecc->wordLen-1;i>=0; i--) {
        printf(HEX_FORMAT, tmp[i]);
    }
    printf("\n");*/
    copy(c, tmp, 2*ecc->wordLen);
    //GFDump(ecc, c);
    u64 carry = 0;
    tmp[0] = c[3];
    tmp[1] = c[3];
    tmp[2] = 0;
    carry += add(3, c, tmp, c);
    //GFDump(ecc, c);
    tmp[0] = 0;
    tmp[1] = c[4];
    tmp[2] = c[4];
    carry = add(3, c, tmp, c);
    if (carry) sub(ecc->wordLen, c, ecc->p, c);
    //GFDump(ecc, c);
    tmp[0] = c[5];
    tmp[1] = c[5];
    tmp[2] = c[5];
    carry = add(3, c, tmp, c);
    //GFDump(ecc, c);
    if (carry) sub(ecc->wordLen, c, ecc->p, c);
    //GFDump(ecc, c);
    if ((c[2] == 0xFFFFFFFFFFFFFFFF) && (c[1] == 0xFFFFFFFFFFFFFFFF)) {
        sub(ecc->wordLen, c, ecc->p, c);
    }
    c[3] = 0;
}

void GFSqr_FIPS192(EcEd* ecc, GFElement a, _out_ GFElement b) {
    GFMul_FIPS192(ecc, a, a, b);
}

/* FIPS-224 Fp: p = 2^224 - 2^96 - 1 */

void GFMul_FIPS224(EcEd* ecc, GFElement a, GFElement b, _out_ GFElement c) {
    GFElement res, tmp;
    mul(ecc, a, b, res);
    copy(c, res, ecc->wordLen);

    u64 carry = 0, borrow = 0;
    c[3] = c[3] & 0xFFFFFFFF;

    tmp[0] = 0;
    ((u32*)tmp)[2] = 0;
    ((u32*)tmp)[3] = ((u32*)res)[7];
    ((u32*)tmp)[4] = ((u32*)res)[8];
    ((u32*)tmp)[5] = ((u32*)res)[9];
    ((u32*)tmp)[6] = ((u32*)res)[10];
    ((u32*)tmp)[7] = 0;

    add(4, tmp, c, c);

    ((u32*)tmp)[3] = ((u32*)res)[11];
    ((u32*)tmp)[4] = ((u32*)res)[12];
    ((u32*)tmp)[5] = ((u32*)res)[13];
    ((u32*)tmp)[6] = 0;

    add(4, tmp, c, c);

    ((u32*)tmp)[0] = ((u32*)res)[11];
    ((u32*)tmp)[1] = ((u32*)res)[12];
    ((u32*)tmp)[2] = ((u32*)res)[13];
    ((u32*)tmp)[3] = 0;
    tmp[2] = 0;
    tmp[3] = 0;
    sub(4, c, tmp, c);

    ((u32*)tmp)[0] = ((u32*)res)[7];
    ((u32*)tmp)[1] = ((u32*)res)[8];
    ((u32*)tmp)[2] = ((u32*)res)[9];
    ((u32*)tmp)[3] = ((u32*)res)[10];
    ((u32*)tmp)[4] = ((u32*)res)[11];
    ((u32*)tmp)[5] = ((u32*)res)[12];
    ((u32*)tmp)[6] = ((u32*)res)[13];

    sub(4, c, tmp, c);
    if (MSB_M & c[3]) {
        add(4, c, p224, c);
    }
    else if (((u32*)c)[7] == 1) {
        sub(4, c, p224, c);
    }
}

void GFSqr_FIPS224(EcEd* ecc, GFElement a, _out_ GFElement b) {
    GFMul_FIPS224(ecc, a, a, b);
}

void GFMulBy2(EcEd* ecc, GFElement a, GFElement b) {
    u64 carry = 0;
    copy(b, a, ecc->wordLen);
    b[ecc->wordLen] = 0;
    mul2(ecc->wordLen, b);
    if (ecc->bitLen % 64 == 0) {
        carry = b[ecc->wordLen] & 1;
        //printf("CARRY %d!", carry);
    }
    else {
        carry = (b[ ecc->wordLen -1 ] & (1UL<<32) );
    }
    if (carry) {
        sub(ecc->wordLen, b, ecc->p, b);
    }
}

void GFPow(EcEd* ecc, GFElement a, BigInt n, GFElement b) {
    GFElement tmp, tmp2;
    copy(tmp, a, ecc->wordLen);
    memset(tmp2, 0, 8*ecc->wordLen);
    tmp2[0] = 1;
    for (u64 i=0; i<ecc->bitLen; i++) {
        if (get_bit(n, i)) {
            ecc->GFMul(ecc, tmp2, tmp, tmp2);
            //printf("%d : ", i);
            //GFDump(ecc, tmp2);
        }
        ecc->GFSqr(ecc, tmp, tmp);
    }
    copy(b, tmp2, ecc->wordLen);
}

void GFInv(EcEd* ecc, GFElement a, GFElement b) {
    GFPow(ecc, a, ecc->p_min_two, b);
}

int GFCmp(EcEd* ecc, GFElement a, GFElement b) {
    for (int i=ecc->wordLen - 1; i>=0; i--) {
        if (a[i] != b[i]) return 1;
    }
    return 0;
}

int EcEdInit(EcEd* ecc, EcPoint* bp, u64 bitLen, BigInt n, GFElement d) {
    srand(time(NULL));
    if ( (bitLen != 192) && (bitLen != 224) ) {
        return -1;
    }
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
    }

    BigInt two;
    memset(two, 0, ecc->wordLen * 8); two[0] = 2;
    sub(ecc->wordLen, ecc->p, two, ecc->p_min_two);
    return 0;
}

/* x^2 + y^2 = 1 + dx^2y^2 */
/* y^2 = (1 - x^2)/(1 - dx^2) */
void EcEdGenerateBasePoint(EcEd* ecc, EcPoint* bp) {
    randomize(ecc->wordLen, bp->x);
    int y = 0;
    GFElement t1, t2, t3, pp;
    copy(pp, ecc->p, ecc->wordLen);
    div2(ecc->wordLen, pp);

    do {
        ecc->GFSqr(ecc, bp->x, t1);
        GFNeg(ecc, t1, t2);
        GFAdd(ecc, t2, unity, t2);
        ecc->GFMul(ecc, t1, ecc->d, t1);
        GFNeg(ecc, t1, t1);
        GFAdd(ecc, t1, unity, t1);
        GFInv(ecc, t1, t1);
        ecc->GFMul(ecc, t1, t2, t1);
        GFPow(ecc, t1, pp, t2);
        y = !GFCmp(ecc, t2, unity);
    } while(y == 0);
    
}

void EcEdAdd(EcEd* ecc, EcPoint* A, EcPoint* B, _out_ EcPoint* C) {
    GFElement z1, z2, z3, z4, z5, z6, z7;
    ecc->GFMul(ecc, A->x, B->x, z1); // z1 = x1 * x2
    //printf("z1\n");
    //GFDump(ecc, z1);
    ecc->GFMul(ecc, A->y, B->y, z2); // z2 = y1 * y2
    //printf("z2\n");
    //GFDump(ecc, z2);
    ecc->GFMul(ecc, z1, z2, z3); 
    ecc->GFMul(ecc, z3, ecc->d, z3); // z3 = d * x1 * x2 * y1 * y2
    //printf("z3\n");
    //GFDump(ecc, z3);

    GFNeg(ecc, z3, z4); // z4 = - z3
    //printf("z4\n");
    //GFDump(ecc, z4);
    //inc(ecc->wordLen, z3);
    GFAdd(ecc, z3, unity, z3);
    //printf("z3\n");
    //GFDump(ecc, z3);

    GFInv(ecc, z3, z3);
    //printf("z3\n");
    //GFDump(ecc, z3);

    GFAdd(ecc, z4, unity, z4);
    //printf("z4\n");
    //GFDump(ecc, z4);

    GFInv(ecc, z4, z4);
    //printf("z4\n");
    //GFDump(ecc, z4);

    ecc->GFMul(ecc, A->x, B->y, z5); // z5 = x1 * y2
    //printf("z5\n");
    //GFDump(ecc, z5);
    ecc->GFMul(ecc, A->y, B->x, z6); // z6 = x2 * y1
    //printf("z6\n");
    //GFDump(ecc, z6);
    GFAdd(ecc, z5, z6, z5);
    GFSub(ecc, z2, z1, z2);
    //printf("z2 z4\n");

    //GFDump(ecc, z2);
    //GFDump(ecc, z4);
    ecc->GFMul(ecc, z5, z3, C->x);
    ecc->GFMul(ecc, z2, z4, C->y);
}

void EcEdDouble(EcEd* ecc, EcPoint* A, EcPoint* B) {
    GFElement z1, z2, z3, z4, z5;
    ecc->GFSqr(ecc, A->x, z1);
    ecc->GFSqr(ecc, A->y, z2);
    ecc->GFMul(ecc, A->x, A->y, z3);
    GFMulBy2(ecc, z3, z3);
    ecc->GFMul(ecc, z1, z2, z4);
    ecc->GFMul(ecc, z4, ecc->d, z4);
    GFNeg(ecc, z4, z5);
    GFAdd(ecc, z4, unity, z4);
    GFInv(ecc, z4, z4);
    GFAdd(ecc, z5, unity, z5);
    GFInv(ecc, z5, z5);
    GFSub(ecc, z2, z1, z2);
    ecc->GFMul(ecc, z4, z3, B->x);
    ecc->GFMul(ecc, z2, z5, B->y);
}

void EcEdScalarMulOrdinary(EcEd* ecc, EcPoint* A, BigInt k, _out_ EcPoint* B) {
    EcPoint P, H;
    copy_point(&P, &uP, ecc->wordLen);
    copy_point(&H, A,  ecc->wordLen);

    for (u32 i=0; i<ecc->bitLen; i++) {
        if (get_bit(k, i)) {
            EcEdAdd(ecc, &P, &H, &P);
            /*printf("%d : ", i);
            GFDump(ecc, P.x);
            GFDump(ecc, P.y);
            printf("\n");*/
        }
        EcEdDouble(ecc, &H, &H); 
        /*printf("%d_ : ", i);
        GFDump(ecc, H.x);
        GFDump(ecc, H.y);
        printf("\n");*/
    }

    copy_point(B, &P, ecc->wordLen);
}

static inline void randomize(u64 len, GFElement a) {
    for (u64 i=0; i<len; i++)
        a[i] = ((u64)rand() << 32) ^ (u64)rand();    
}

void EcEdConvertAffineToProjective(EcEd* ecc, EcPoint* P, EcPointProj* Q) {
    randomize(ecc->wordLen, Q->Z);
    ecc->GFMul(ecc, P->x, Q->Z, Q->X);
    ecc->GFMul(ecc, P->y, Q->Z, Q->Y);
}

void EcEdConvertProjectiveToAffine(EcEd* ecc, EcPointProj* P, EcPoint* Q) {
    GFElement Z_inv;
    GFInv(ecc, P->Z, Z_inv);
    ecc->GFMul(ecc, P->X, Z_inv, Q->x);
    ecc->GFMul(ecc, P->Y, Z_inv, Q->y);
}

void EcEdAddProj(EcEd* ecc, EcPointProj* A, EcPointProj* B, _out_ EcPointProj* C) {
    GFElement a,b,c,d,e,e1,e2,e3,f;
    ecc->GFMul(ecc, A->Z, B->Z, a);
    ecc->GFSqr(ecc, a, b);
    ecc->GFMul(ecc, A->X, B->X, c);
    ecc->GFMul(ecc, A->Y, B->Y, d);
    GFAdd(ecc, A->X, A->Y, e);
    GFAdd(ecc, B->X, B->Y, e1);
    ecc->GFMul(ecc, e, e1, e);
    GFSub(ecc, e, c, e);
    GFSub(ecc, e, d, e);
    ecc->GFMul(ecc, c, d, f);
    ecc->GFMul(ecc, f, ecc->d, f);

    GFSub(ecc, b, f, e2);
    ecc->GFMul(ecc, e, a, C->X);
    ecc->GFMul(ecc, C->X, e2, C->X);

    GFAdd(ecc, b, f, e3);
    GFSub(ecc, d, c, C->Y);
    ecc->GFMul(ecc, C->Y, e3, C->Y);
    ecc->GFMul(ecc, C->Y, a, C->Y);

    ecc->GFMul(ecc, e2, e3, C->Z);
}

void EcEdDoubleProj(EcEd* ecc, EcPointProj* A, _out_ EcPointProj* B) {
    GFElement a,b,c,d,e,e1,e2,e3,f;
    ecc->GFSqr(ecc, A->Z, a);
    ecc->GFSqr(ecc, a, b);
    ecc->GFSqr(ecc, A->X, c);
    ecc->GFSqr(ecc, A->Y, d);

    ecc->GFMul(ecc, A->X, A->Y, e);
    GFMulBy2(ecc, e, e);

    ecc->GFMul(ecc, c, d, f);
    ecc->GFMul(ecc, f, ecc->d, f);

    GFSub(ecc, b, f, e2);
    ecc->GFMul(ecc, e, a, B->X);
    ecc->GFMul(ecc, B->X, e2, B->X);

    GFAdd(ecc, b, f, e3);
    GFSub(ecc, d, c, B->Y);
    ecc->GFMul(ecc, B->Y, e3, B->Y);
    ecc->GFMul(ecc, B->Y, a, B->Y);

    ecc->GFMul(ecc, e2, e3, B->Z);
}

void EcEdScalarMul(EcEd* ecc, EcPoint* A, BigInt k, _out_ EcPoint* B) {
    EcPointProj P, H;
    EcEdConvertAffineToProjective(ecc, &uP, &P);
    EcEdConvertAffineToProjective(ecc, A, &H);

    for (u32 i=0; i<ecc->bitLen; i++) {
        if (get_bit(k, i)) {
            EcEdAddProj(ecc, &P, &H, &P);
            /*printf("%d : ", i);
            GFDump(ecc, P.x);
            GFDump(ecc, P.y);
            printf("\n");*/
        }
        EcEdDoubleProj(ecc, &H, &H); 
        /*printf("%d_ : ", i);
        GFDump(ecc, H.x);
        GFDump(ecc, H.y);
        printf("\n");*/
    }

    EcEdConvertProjectiveToAffine(ecc, &P, B);
}