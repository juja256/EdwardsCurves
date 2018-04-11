#include "gf.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

#ifdef _WIN64
#include <intrin.h>
#pragma intrinsic(_umul128) 
#else
#include <x86intrin.h>
#endif // _WIN64


#define MAX_U64    0xFFFFFFFFFFFFFFFF
#define MSB_M      0x8000000000000000
#define HEX_FORMAT "%.16llX"

const u64 unity[] = { 0x1, 0, 0, 0, 0, 0, 0, 0, 0 };

const u64 zero[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0  };

const EcPoint uP  = { { 0, 0, 0, 0, 0, 0, 0, 0, 0 }, { 1, 0, 0, 0, 0, 0, 0, 0, 0} };

const u64 p192[]  = { 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFE, 0xFFFFFFFFFFFFFFFF };

const u64 p256[] = { 0xFFFFFFFFFFFFFFFF, 0x00000000FFFFFFFF, 0x0000000000000000, 0xFFFFFFFF00000001};

const u64 p224[]  = { 0x0000000000000001, 0xFFFFFFFF00000000, 0xFFFFFFFFFFFFFFFF, 0x00000000FFFFFFFF };

const u64 p384[]  = { 0x00000000FFFFFFFF, 0xFFFFFFFF00000000, 0xFFFFFFFFFFFFFFFE, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF };

const u64 p521[] = { 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 
                     0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0x1FF };             

inline u64 get_bit(const u64* a, u64 num) {
    return a[num/64] & ((u64)1 << (num % 64));
}

inline void copy(u64* a, const u64* b, int len) {
    for (int i=0;i<len;i++)
        a[i] = b[i];
}

inline u64 add(u64 n, const u64* a, const u64* b, u64* c) {
    u64 msb_a, msb_b, carry = 0;

    for (u32 i=0;i < n; i++) {
        msb_a = a[i] & MSB_M;
        msb_b = b[i] & MSB_M;
        c[i] = a[i] + b[i] + carry;
        carry = ( (msb_a && msb_b) || ((msb_a ^ msb_b) && !(MSB_M & c[i])) );
    }
    return carry;
}

inline u64 sub(u64 n, const u64* a, const u64* b, u64* c) {
    u64 borrow = 0;
    for (int i=0; i<n; i++) {
        u64 t_a = a[i];
        u64 t_b = b[i];
        c[i] = t_a - t_b - borrow;
        borrow = ( (~t_a) & (c[i] | t_b) | (c[i] & t_b) ) >> (63);
    }
    return borrow;
}

inline void _mul_raw(u64 a, u64 b, u64* low, u64* high) {
#ifdef _WIN64
    *low = _umul128(a, b, high);
#else
    __int128 r = (__int128)a * (__int128)b;
    *low = (unsigned long long)r;
    *high = r >> 64;
#endif // _WIN64
}

inline u64 _add_raw(u64 a, u64 b, u64* c) {
#ifdef _WIN64
    u64 msb_a, msb_b;
    msb_a = a & MSB_M;
    msb_b = b & MSB_M;
    *c = a + b;
    return ((msb_a && msb_b) || ((msb_a ^ msb_b) && !(MSB_M & *c)));
#else 
    __int128 r = (__int128)a + (__int128)b;
    *c = (u64)r;
    return r >> 64;
#endif 
}

inline void mul_by_word(u64 n, const u64* a, u64 d, u64* c) {
    u64 carry = 0, carry_tmp;
    for (int i=0; i < n; i++) {
        carry_tmp = carry;
        _mul_raw(d, a[i], &(c[i]), &carry);
        carry += _add_raw(c[i], carry_tmp, &(c[i]));
    }
    c[n] = carry;
}

inline void mul(u64 n, const u64* a, const u64* b, u64* c) {
    GFElement tmp;
    memset(c, 0, 2*8*n);
    for (u64 i=0; i < n; i++) {
        mul_by_word(n, a, b[i], (u64*)tmp);
        add(n + 1, c+i, tmp, c+i);
    }
}

inline int word_bit_len(u64 n) {
    int c = 64;
    while (c) {
        if ((((u64)1 << (64-1)) & n) >> (64-1))
            return c;
        n <<= 1;
        --c;
    }
    return 0;
}

inline int bigint_bit_len(u64 nWords, const u64* a) {
    int bit_len = nWords * 64;
    int i=nWords-1;
    do {
        bit_len-=64;
    } while ((i>=0) && (a[i--] == 0));

    bit_len += word_bit_len(a[i+1]);
    return bit_len;
}

inline void shl(u64 n, const u64* a, u64* res, u64 bits) {
    u64 buf = 0;
    int chk = bits / 64;
    bits = bits % 64;
    u64 cur;

    if (bits) {
        for (int i = 0; i < n; i++) {
            cur = a[i];
            res[i+chk] = (cur << bits) ^ buf;
            buf = (cur & (MAX_U64 << ( 64-bits ))) >> ( 64-bits );
        }
    }
    else {
        for (int i = 0; i < n; i++) {
            cur = a[i];
            res[i+chk] = cur;
        }
    }

    for (int i = 0; i < chk; i++) {
        res[i] = 0;
    }
    res[n+chk] = buf;
}

inline void shr(u64 n, const u64* a, u64* res, u64 bits) {
    u64 buf = 0;
    int chk = bits / 64;
    bits = bits % 64;
    u64 cur;

    if (bits) {
        for (int i = n-1+chk; i >= chk; i--) {
            cur = a[i];
            res[i-chk] = (cur >> bits) ^ buf;
            buf = (cur & (MAX_U64 >> ( 64-bits ))) << (64 - bits);
        }
    }
    else {
        for (int i = n-1+chk; i >= chk; i--) {
            cur = a[i];
            res[i-chk] = cur;
        }
    }
    for (int i = n; i<n+chk; i++) {
        res[i] = 0;
    }
}

static inline void dump(u64 n, const u64* a) {
    for (int i=n-1; i>=0; i--)
        printf(HEX_FORMAT, a[i]);
    printf("\n");
}

inline int cmp(u64 n, const u64* a, const u64* b) {
    for (int i = n - 1; i >= 0; i--)
    {
        if (a[i] > b[i]) return 1;
        else if (a[i] < b[i]) return -1;
    }
    return 0;
}

void zero_int(u64 n, u64* a) {
    for (int i=0; i<n; i++) a[i] = 0;
}

void add_mod(u64 n, const BigInt a, const BigInt b, const BigInt m, BigInt res) {
    res[n] = add(n, a, b, res);
    if (cmp(n+1, res, m) == 1) sub(n+1, res, m, res);
}

void mul_mod(u64 n, const BigInt a, const BigInt b, const BigInt m, BigInt res) {
    /* new multiplication with reduction by division */
    VeryBigInt d;
    mul(n, a, b, d);
    divide(n, d, m, NULL, res);
    
    /* old multiplication using only additions */
    /*
    u64 b_len = bigint_bit_len(n, b);
    BigInt mm, r;
    zero_int(n, r);
    
    copy(mm, a, n);
    for (int i=0; i<b_len; i++) {
        if (get_bit(b, i)) add_mod(n, r, mm, m, r);
        add_mod(n, mm, mm, m, mm);
    }
    copy(res, r, n);
    */    
}

void exp_mod(u64 n, const BigInt a, const BigInt p, const BigInt m, BigInt res) {
    u64 b_len = bigint_bit_len(n, p);
    VeryBigInt mm, r;
    zero_int(n+1, r);
    r[0] = 1;
    copy(mm, a, n);
    for (int i=0; i<b_len; i++) {
        if (get_bit(p, i)) mul_mod(n, r, mm, m, r);
        mul_mod(n, mm, mm, m, mm);
    }
    copy(res, r, n);
}

void imul(u64 n, const u64* a, const u64* b, u64* c) {
    /* c := a * b, where b could be signed value */
    int b_isneg = b[n] & MSB_M;
    BigInt bb;
    if (b_isneg) {
        sub(n+1, zero, b, bb);
    }
    else {
        copy(bb, b, n+1);
    }
    
    mul(n+1, a, bb, c);
    
    if (b_isneg) {
        sub(2*n, zero, c, c);
    }
}

void inv_mod(u64 n, const BigInt a, const BigInt m, BigInt res) {
    /* Old realization of inversion using powering to p-2 */
    /*
    BigInt mm;
    copy(mm, m, n);
    mm[0] -= 2; 
    exp_mod(n, a, mm, m, res);
    */

    /* Compute a^{-1} mod m with Extended Euclidean Algorithm */
    VeryBigInt q, tmp, r;
    BigInt t, newt, newr;
    zero_int(2*n, r);
    zero_int(n+1, t); // t := 0
    zero_int(n+1, newt); newt[0] = 1; // newt := 1
    copy(r, m, n); // r := m
    copy(newr, a, n); // newr := a

    while(cmp(n, newr, zero) != 0) {
        divide(n, r, newr, q, tmp); // q := r div newr, tmp := r mod newr
        copy(r, newr, n); // r := newr
        copy(newr, tmp, n); // newr := tmp 

        imul(n, q, newt, tmp); // tmp := q*newt

        sub(n+1, t, tmp, tmp); // tmp := t - tmp
        copy(t, newt, n+1);
        copy(newt, tmp, n+1);
    }
    
    if (t[n] & MSB_M) {
        add(n, t, m, res);
    }
    else {
        copy(res, t, n);
    }
}

void divide(u64 n, const u64* a, const u64* b, u64* quotient, u64* reminder) {
    VeryBigInt q;
    VeryBigInt tmp;
    VeryBigInt r;

    copy(r, a, 2*n);
    zero_int(2*n, q);
    zero_int(2*n, tmp);

    int k = bigint_bit_len(2*n, a);
    int t = bigint_bit_len(n, b);

    if (k < t) {
        if (quotient) copy(quotient, q, 2*n);
        copy(reminder, r, n);
        return;
    }

    k = k-t;
    shl(n, b, tmp, k);
    while (k >= 0) {
        if (sub(2*n, r, tmp, r) == 0) {
            q[ k/64 ] |= (u64)1 << (k % 64);
        }
        else {
            add(2*n, r, tmp, r);
        }
        
        div2(2*n, tmp);
        k--;
    }

    if (quotient) copy(quotient, q, 2*n);
    copy(reminder, r, n);
}

void GFSqrt(const Ec* ecc, const GFElement a, GFElement r) {
    // p = Q*2^s
    GFElement Q, pp, z, tmp, c, t, R, b;
    copy(z, unity, ecc->wordLen);
    z[0]++;
    copy(Q, ecc->p, ecc->wordLen);
    u64 M, S = 1;
    div2(ecc->wordLen, Q);
    copy(pp, Q, ecc->wordLen);
    //p-1 = 2^s*Q
    while ((Q[0] & 1) == 0) {
        div2(ecc->wordLen, Q);
        S++;
    }
    //finding z - non quadratic residue
    while(1) {
        GFPow(ecc, z, pp, tmp); 
        if (GFCmp(ecc, tmp, unity))
            break;
        else 
            GFAdd(ecc, z, unity, z);
    }
    
    M = S;
    GFPow(ecc, z, Q, c);
    GFPow(ecc, a, Q, t);
    GFAdd(ecc, Q, unity, tmp);
    div2(ecc->wordLen, tmp);
    u64 i=1;
    GFPow(ecc, a, tmp, R); // tmp = (Q+1)/2

    while (1) {
        if (!GFCmp(ecc, t, unity)) {
            copy(r, R, ecc->wordLen);
            break;
        } 
        copy(tmp,t,ecc->wordLen); //copying value of t, we'll need it later
        for (i=1; i<M; i++) {
            GFSqr(ecc, t, t);
            if (!GFCmp(ecc, t, unity))
                break;
        }
        copy(b, c, ecc->wordLen);
        for (u64 j=0; j<M-i-1; j++) {
            GFSqr(ecc, b, b);
        }
        M = i;
        GFSqr(ecc, b, c);
        GFMul(ecc, tmp, c, t); //original value of t multiplied by c^2 
        GFMul(ecc, R, b, R);
    }
    
}

void GFInitFromString(GFElement a, const char* str) {
    u64 s_len = strlen(str);
    u64 tmp;

    memset(a, 0, sizeof(GFElement));

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

void GFDump(const Ec* ecc, const GFElement a) {
    for (int i=ecc->wordLen-1; i>=0; i--)
        printf(HEX_FORMAT, a[i]);
    printf("\n");
}


void GFNeg(const Ec* ecc, const GFElement a, GFElement c) {
    sub(ecc->wordLen, ecc->p, (u64*)a, (u64*)c);
}

void GFAdd(const Ec* ecc, const GFElement a, const GFElement b, GFElement c) {
    u64 carry = 0;
    switch (ecc->bitLen) {
        case 192:
        case 256:
        case 384: {
            carry = add(ecc->wordLen, a,b,c);
            if (carry) {
                sub(ecc->wordLen, c, ecc->p, c);
            }
            if(GFCmp(ecc,c,ecc->p)!=-1) {
                sub(ecc->wordLen, c, ecc->p, c);
            }
            break;
        }
        case 224: {
            add(ecc->wordLen, a, b, c);
            carry = (c[ecc->wordLen - 1] & ((u64)1 << 32));
            if (carry) {
                sub(ecc->wordLen, c, ecc->p, c);
            }
            if(GFCmp(ecc,c,ecc->p)!=-1) {
                sub(ecc->wordLen, c, ecc->p, c);
            }
            break;
        }
        case 521: {
            add(ecc->wordLen, a, b, c);
            carry = (c[ecc->wordLen - 1] & ((u64)1 << 9));
            /*if (carry) {
                sub(ecc->wordLen, c, ecc->p, c);
            }*/
            if(GFCmp(ecc, c, ecc->p) != -1)
                sub(ecc->wordLen, c, ecc->p, c);
            break;
        }
    }
}

void GFSub(const Ec* ecc, const GFElement a, const GFElement b, GFElement c) {
    GFElement tmp;
    copy(tmp, a, ecc->wordLen);                           
    sub(ecc->wordLen, ecc->p, (u64*)b, (u64*)c);
    GFAdd(ecc, tmp, c, c);
}

/* FIPS-192 Fp: p = 2^192 - 2^64 - 1*/

void GFMul_FIPS192(const Ec* ecc, const GFElement a, const GFElement b, GFElement res) {
    BigInt tmp, c;
    mul(3, a, b, c);

    u64 carry = 0;

    // S_1
    tmp[0] = c[3];
    tmp[1] = c[3];
    tmp[2] = 0;
    
    carry += add(3, c, tmp, res);
    
    // S_2
    tmp[0] = 0;
    tmp[1] = c[4];
    tmp[2] = c[4];

    carry = add(3, res, tmp, res);
    if (carry) sub(3, res, ecc->p, res);
    
    // S_3
    tmp[0] = c[5];
    tmp[1] = c[5];
    tmp[2] = c[5];

    carry = add(3, res, tmp, res);
    if (carry) sub(3, res, ecc->p, res);

    res[3] = 0;

    if(GFCmp(ecc, res, ecc->p) != -1)
        sub(3, res, ecc->p, res);
}

void GFSqr_FIPS192(const Ec* ecc, const GFElement a, GFElement b) {
    GFMul_FIPS192(ecc, a, a, b);
}

/* FIPS-224 Fp: p = 2^224 - 2^96 + 1 */

void GFMul_FIPS224(const Ec* ecc, const GFElement a, const GFElement b, GFElement res) {
    BigInt tmp, c;
    mul(ecc->wordLen, a, b, c);

    copy(res, c, 4);

    u64 carry = 0;
    res[3] &= 0xFFFFFFFF;

    // S_1
    tmp[0] = 0;
    ((u32*)tmp)[2] = 0;
    ((u32*)tmp)[3] = ((u32*)c)[7];
    ((u32*)tmp)[4] = ((u32*)c)[8];
    ((u32*)tmp)[5] = ((u32*)c)[9];
    ((u32*)tmp)[6] = ((u32*)c)[10];
    ((u32*)tmp)[7] = 0;

    add(4, tmp, res, res);
    
    // S_2
    ((u32*)tmp)[3] = ((u32*)c)[11];
    ((u32*)tmp)[4] = ((u32*)c)[12];
    ((u32*)tmp)[5] = ((u32*)c)[13];
    ((u32*)tmp)[6] = 0;

    add(4, tmp, res, res);

    // D_1
    ((u32*)tmp)[0] = ((u32*)c)[11];
    ((u32*)tmp)[1] = ((u32*)c)[12];
    ((u32*)tmp)[2] = ((u32*)c)[13];
    ((u32*)tmp)[3] = 0;
    tmp[2] = 0;
    tmp[3] = 0;
    
    sub(4, res, tmp, res);

    // D_2
    ((u32*)tmp)[0] = ((u32*)c)[7];
    ((u32*)tmp)[1] = ((u32*)c)[8];
    ((u32*)tmp)[2] = ((u32*)c)[9];
    ((u32*)tmp)[3] = ((u32*)c)[10];
    ((u32*)tmp)[4] = ((u32*)c)[11];
    ((u32*)tmp)[5] = ((u32*)c)[12];
    ((u32*)tmp)[6] = ((u32*)c)[13];

    sub(4, res, tmp, res);
    
	while(((int*)res)[7] < 0) { 
        add(4, res, ecc->p, res);
    }
    while(((int*)res)[7] > 0) { 
        sub(4, res, ecc->p, res);
    }

    if(GFCmp( ecc, res, ecc->p) != -1) {
        sub(4, res, ecc->p, res);
    }
}

void GFSqr_FIPS224(const Ec* ecc, const GFElement a, GFElement b) {
    GFMul_FIPS224(ecc, a, a, b);
}

/*p = 2^256 – 2^224 + 2^192 + 2^96 – 1*/
void GFMul_FIPS256(const Ec* ecc,const GFElement a,const GFElement b, GFElement res) 
{
    GFElement tmp, c;

    mul(4, a, b, c);


    int carry; //signed, its important
    tmp[4] = 0;
    tmp[3] = c[7];
    tmp[2] = c[6];
    ((u32*)tmp)[3] = ((u32*)c)[11];
    ((u32*)tmp)[2] = 0;
    tmp[0] = 0;
    mul2(4, tmp);
    carry = tmp[4]; //in case 2*tmp > 2^256
    carry += add(4, c, tmp, res);

    tmp[4] = 0;
    ((u32*)tmp)[7] = 0;
    ((u32*)tmp)[6] = ((u32*)c)[15];
    ((u32*)tmp)[5] = ((u32*)c)[14];
    ((u32*)tmp)[4] = ((u32*)c)[13];
    ((u32*)tmp)[3] = ((u32*)c)[12];
    mul2(4, tmp);
    carry += tmp[4];
    carry += add(4, res, tmp, res);

    tmp[3] = c[7];
    tmp[2] = 0;
    ((u32*)tmp)[3] = 0;
    ((u32*)tmp)[2] = ((u32*)c)[10];
    tmp[0] = c[4];
    carry += add(4, res, tmp, res);

    ((u32*)tmp)[7] = ((u32*)c)[8];
    ((u32*)tmp)[6] = ((u32*)c)[13];
    tmp[2] = c[7];
    ((u32*)tmp)[3] = ((u32*)c)[13];
    ((u32*)tmp)[2] = ((u32*)c)[11];
    ((u32*)tmp)[1] = ((u32*)c)[10];
    ((u32*)tmp)[0] = ((u32*)c)[9];
    carry += add(4, res, tmp, res);

    ((u32*)tmp)[7] = ((u32*)c)[10];
    ((u32*)tmp)[6] = ((u32*)c)[8];
    tmp[2] = 0;
    ((u32*)tmp)[3] = 0;
    ((u32*)tmp)[2] = ((u32*)c)[13];
    ((u32*)tmp)[1] = ((u32*)c)[12];
    ((u32*)tmp)[0] = ((u32*)c)[11];
    carry -= sub(4, res, tmp, res);

    ((u32*)tmp)[7] = ((u32*)c)[11];
    ((u32*)tmp)[6] = ((u32*)c)[9];
    tmp[2] = 0;
    tmp[1] = c[7];
    tmp[0] = c[6];
    carry -= sub(4, res, tmp, res);

    ((u32*)tmp)[7] = ((u32*)c)[12];
    ((u32*)tmp)[6] = 0;
    ((u32*)tmp)[5] = ((u32*)c)[10];
    ((u32*)tmp)[4] = ((u32*)c)[9];
    ((u32*)tmp)[3] = ((u32*)c)[8];
    ((u32*)tmp)[2] = ((u32*)c)[15];
    ((u32*)tmp)[1] = ((u32*)c)[14];
    ((u32*)tmp)[0] = ((u32*)c)[13];

    carry-=sub(4, res, tmp, res);

    ((u32*)tmp)[7] = ((u32*)c)[13];
    ((u32*)tmp)[6] = 0;
    tmp[2] = c[5];
    ((u32*)tmp)[3] = ((u32*)c)[9];
    ((u32*)tmp)[2] = 0;
    tmp[0] = c[7];

    carry-=sub(4, res, tmp, res);

    while(carry > 0) 
    {
        sub(4, res, ecc->p, res);
        carry--;
    }
    while(carry < 0) 
    {
        add(4, res, ecc->p, res);
        carry++;
    }
    
    if(GFCmp(ecc, res, ecc->p) != -1)
        sub(4, res, ecc->p, res);
}

void GFSqr_FIPS256(const Ec* ecc,const GFElement a,  GFElement  b) {
    GFMul_FIPS256(ecc, a, a, b);
}

/* FIPS-384 Fp: p = p^384 - 2^128 - 2^96 +  2^32 - 1 */
void GFMul_FIPS384(const Ec* ecc, const GFElement a, const GFElement b, GFElement res) {
    BigInt tmp;
    VeryBigInt c;
    mul(6, a, b, c);

    int carry = 0;

    // 2*S_1
    tmp[0] = 0;
    tmp[1] = 0;
    ((u32*)tmp)[4] = ((u32*)c)[21];
    ((u32*)tmp)[5] = ((u32*)c)[22];
    ((u32*)tmp)[6] = ((u32*)c)[23];
    ((u32*)tmp)[7] = 0;
    tmp[4] = 0;
    tmp[5] = 0;

    mul2(6, tmp);
    carry = tmp[6];
    carry += add(6, tmp, c, res);

    // S_2
    tmp[0] = c[6];
    tmp[1] = c[7];
    tmp[2] = c[8];
    tmp[3] = c[9];
    tmp[4] = c[10];
    tmp[5] = c[11];

    carry += add(6, tmp, res, res);

    // S_3
    ((u32*)tmp)[0] = ((u32*)c)[21];
    ((u32*)tmp)[1] = ((u32*)c)[22];
    ((u32*)tmp)[2] = ((u32*)c)[23];
    ((u32*)tmp)[3] = ((u32*)c)[12];
    ((u32*)tmp)[4] = ((u32*)c)[13];
    ((u32*)tmp)[5] = ((u32*)c)[14];
    ((u32*)tmp)[6] = ((u32*)c)[15];
    ((u32*)tmp)[7] = ((u32*)c)[16];
    ((u32*)tmp)[8] = ((u32*)c)[17];
    ((u32*)tmp)[9] = ((u32*)c)[18];
    ((u32*)tmp)[10] = ((u32*)c)[19];
    ((u32*)tmp)[11] = ((u32*)c)[20];

    carry += add(6, res, tmp, res);

    // S_4
    ((u32*)tmp)[0] = 0;
    ((u32*)tmp)[1] = ((u32*)c)[23];
    ((u32*)tmp)[2] = 0;
    ((u32*)tmp)[3] = ((u32*)c)[20];
    ((u32*)tmp)[4] = ((u32*)c)[12];
    ((u32*)tmp)[5] = ((u32*)c)[13];
    ((u32*)tmp)[6] = ((u32*)c)[14];
    ((u32*)tmp)[7] = ((u32*)c)[15];
    ((u32*)tmp)[8] = ((u32*)c)[16];
    ((u32*)tmp)[9] = ((u32*)c)[17];
    ((u32*)tmp)[10] = ((u32*)c)[18];
    ((u32*)tmp)[11] = ((u32*)c)[19];

    carry += add(6, res, tmp, res);

    // S_5
    tmp[0] = 0;
    tmp[1] = 0;
    tmp[2] = c[10];
    tmp[3] = c[11];
    tmp[4] = 0;
    tmp[5] = 0;

    carry += add(6, res, tmp, res);

    // S_6
    ((u32*)tmp)[0] = ((u32*)c)[20];
    ((u32*)tmp)[1] = 0;
    ((u32*)tmp)[2] = 0;
    ((u32*)tmp)[3] = ((u32*)c)[21];
    ((u32*)tmp)[4] = ((u32*)c)[22];
    ((u32*)tmp)[5] = ((u32*)c)[23];	
    tmp[3] = 0;
    tmp[4] = 0;
    tmp[5] = 0;

    carry += add(6, res, tmp, res);

    // D_1
    ((u32*)tmp)[0] = ((u32*)c)[23];
    ((u32*)tmp)[1] = ((u32*)c)[12];
    ((u32*)tmp)[2] = ((u32*)c)[13];
    ((u32*)tmp)[3] = ((u32*)c)[14];
    ((u32*)tmp)[4] = ((u32*)c)[15];
    ((u32*)tmp)[5] = ((u32*)c)[16];
    ((u32*)tmp)[6] = ((u32*)c)[17];
    ((u32*)tmp)[7] = ((u32*)c)[18];
    ((u32*)tmp)[8] = ((u32*)c)[19];
    ((u32*)tmp)[9] = ((u32*)c)[20];
    ((u32*)tmp)[10] = ((u32*)c)[21];
    ((u32*)tmp)[11] = ((u32*)c)[22];

    carry -= sub(6, res, tmp, res);

    // D_2
    ((u32*)tmp)[0] = 0;
    ((u32*)tmp)[1] = ((u32*)c)[20];
    ((u32*)tmp)[2] = ((u32*)c)[21];
    ((u32*)tmp)[3] = ((u32*)c)[22];
    ((u32*)tmp)[4] = ((u32*)c)[23];
    ((u32*)tmp)[5] = 0;
    tmp[3] = 0;
    tmp[4] = 0;
    tmp[5] = 0;
    
    carry -= sub(6, res, tmp, res);

    // D_3
    tmp[0] = 0;
    ((u32*)tmp)[2] = 0;
    ((u32*)tmp)[3] = ((u32*)c)[23];
    ((u32*)tmp)[4] = ((u32*)c)[23];
    ((u32*)tmp)[5] = 0;
    tmp[3] = 0;
    tmp[4] = 0;
    tmp[5] = 0;
    
    carry -= sub(6, res, tmp, res);

    while(carry > 0) {
        sub(6, res, ecc->p, res);
        carry--;
    }
    while(carry++ < 0) {
        add(6, res, ecc->p, res);
        carry++;
    }
    //in case of 0<c<p256
    if(GFCmp(ecc, res, ecc->p) != -1)
        sub(6, res, ecc->p, res);
}

void GFSqr_FIPS384(const Ec* ecc, const GFElement a, GFElement b) {
    GFMul_FIPS384(ecc, a, a, b);
}

void GFMul_FIPS521(const Ec* ecc, const GFElement a, const GFElement b, GFElement res) {
    VeryBigInt c, tmp;

    mul(ecc->wordLen, a, b, c);
    shr(ecc->wordLen, c, tmp, 521);

    c[8] &= MAX_U64 >> (64-9);
    tmp[8] &= MAX_U64 >> (64-9);

    add(9, c, tmp, res);

    //if ( res[8] & (1 << 9) ) sub(9, res, ecc->p, res);
    if(GFCmp(ecc, res, ecc->p) != -1)
        sub(9, res, ecc->p, res);
}

void GFSqr_FIPS521(const Ec* ecc, const GFElement a, GFElement c) {
    GFMul_FIPS521(ecc, a, a, c);
}

void GFMul_Cmn(const Ec* ecc, const GFElement a, const GFElement b, GFElement c) {
    mul_mod(ecc->wordLen, a, b, ecc->p, c);
}

void GFSqr_Cmn(const Ec* ecc, const GFElement a, GFElement c) {
    GFMul_Cmn(ecc, a, a, c);
}

void GFMulBy2Power(const Ec* ecc, const GFElement a, int pp, GFElement b) {
    u64 carry = 0;
    copy(b, a, ecc->wordLen);
    b[ecc->wordLen] = 0;
    shl(ecc->wordLen, b, b, pp);

    if (ecc->bitLen % 64 == 0) {
        carry = b[ecc->wordLen] & ((1<<pp)-1);
        while(carry) {
            sub(ecc->wordLen+1, b, ecc->p, b);
            carry = b[ecc->wordLen] & ((1<<pp)-1);
        }
    }
    else {
        carry = b[ ecc->wordLen -1 ] & ((u64)(-1)<< ( ecc->bitLen % 64 ) );
        while(carry) {
            sub(ecc->wordLen, b, ecc->p, b);
            carry = b[ ecc->wordLen -1 ] & ((u64)(-1)<< (ecc->bitLen % 64) );
        }
    }
}

void GFMulByD(const EcEd* ecc, GFElement a) {
    if (ecc->d_len < 64) {
        VeryBigInt t;
        mul_by_word(ecc->wordLen, a, ecc->d[0], t);
        divide(ecc->wordLen, t, ecc->p, NULL, a);
    }
    else {
        GFMul(ecc, a, ecc->d, a);
    }
}

void GFPow(const Ec* ecc, const GFElement a, const BigInt n, GFElement b) {
    GFElement tmp, tmp2;
    copy(tmp, a, ecc->wordLen);
    memset(tmp2, 0, 8*ecc->wordLen);
    tmp2[0] = 1;
    for (u64 i=0; i<ecc->bitLen; i++) {
        if (get_bit(n, i)) {
            GFMul(ecc, tmp2, tmp, tmp2);
        }
        GFSqr(ecc, tmp, tmp);
    }
    copy(b, tmp2, ecc->wordLen);
}

void GFInv(const Ec* ecc, const GFElement a, GFElement b) {
    //GFPow(ecc, a, ecc->p_min_two, b); 
    
    /* inv_mod(Euclidean alg.) works slightly better than powering */
    inv_mod(ecc->wordLen, a, ecc->p, b);
}

void GFMul(const Ec* ecc, const GFElement a, const GFElement b, GFElement c) {
    ecc->GFMul(ecc, a, b, c);
}

void GFSqr(const Ec* ecc, const GFElement a, GFElement c) {
    ecc->GFSqr(ecc, a, c);
}

int GFCmp(const Ec* ecc, const GFElement a, const GFElement b) {
    for (int i=ecc->wordLen - 1; i>=0; i--)
    {
        if (a[i] > b[i]) return 1;
        else if (a[i] < b[i]) return -1;
    }
    return 0;
}