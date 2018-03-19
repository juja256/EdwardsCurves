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

const u64 unity[] = { 0x1, 0, 0, 0, 0, 0 };

const EcPoint uP  = { { 0, 0, 0, 0, 0, 0 }, { 1, 0, 0, 0, 0, 0} };

const u64 p192[]  = { 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFE, 0xFFFFFFFFFFFFFFFF };

const u64 p256[] = { 0xFFFFFFFFFFFFFFFF, 0x00000000FFFFFFFF, 0x0000000000000000, 0xFFFFFFFF00000001};

const u64 p224[]  = { 0x0000000000000001, 0xFFFFFFFF00000000, 0xFFFFFFFFFFFFFFFF, 0x00000000FFFFFFFF };

const u64 p384[]  = { 0x00000000FFFFFFFF, 0xFFFFFFFF00000000, 0xFFFFFFFFFFFFFFFE, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF };
                                          

inline u64 get_bit(const GFElement a, u64 num) {
    return a[num/64] & ((u64)1 << (num % 64));
}

inline void copy(GFElement a, const GFElement b, int len) {
    for (int i=0;i<len;i++)
        a[i] = b[i];
}

inline void copy_point(EcPoint* P, const EcPoint* Q, int len) {
    copy(P->x, Q->x, len);
    copy(P->y, Q->y, len);
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

inline void mul_by_word(const Ec* ecc, const u64* a, u64 d, u64* c) {
    u64 carry = 0, carry_tmp;
    for (int i=0; i < ecc->wordLen; i++) {
        carry_tmp = carry;
        _mul_raw(d, a[i], &(c[i]), &carry);
        carry += _add_raw(c[i], carry_tmp, &(c[i]));
    }
    c[ecc->wordLen] = carry;
}

inline void mul(const Ec* ecc, const u64* a, const u64* b, u64* c) {
    GFElement tmp;
    memset(c, 0, 2*8*ecc->wordLen);
    for (u64 i=0; i < ecc->wordLen; i++) {
        mul_by_word(ecc, a, b[i], (u64*)tmp);
        add(ecc->wordLen + 1, c+i, tmp, c+i);
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

inline int bigint_bit_len(u64 nWords, const BigInt a) {
    int bit_len = nWords * 64;
    int i=nWords-1;
    do {
        bit_len-=64;
    } while (a[i--] == 0);

    bit_len += word_bit_len(a[i+1]);
    return bit_len;
}

inline void shl(u64 n, const GFElement a, GFElement res, u64 bits) {
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

inline void shr(u64 n, const GFElement a, GFElement res, u64 bits) {
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



inline void dump(u64 n, const BigInt a) {
    for (int i=n-1; i>=0; i--)
        printf(HEX_FORMAT, a[i]);
    printf("\n");
}

inline int cmp(u64 n, const BigInt a, const BigInt b) {
    for (int i = n - 1; i >= 0; i--)
    {
        if (a[i] > b[i]) return 1;
        else if (a[i] < b[i]) return -1;
    }
    return 0;
}

void basic_reduction(u64 n, const BigInt a, const BigInt p, BigInt res) {
    BigInt tmp;
    BigInt m;

    copy(res, a, 2*n);

    int k = bigint_bit_len(n, p);
    int t;
    

    while (cmp(2*n, res, p) != -1) {
        t = bigint_bit_len(2*n, res);
        shl(n, p, m, t-k);
        if (cmp(2*n, m, res) == 1) {
            t--;
            shl(n, p, m, t-k);
        }
        sub(2*n, res, m, res);
        printf("%d %d\n", k,t);
    }
}

inline void randomize(u64 len, GFElement a) {
    for (u64 i=0; i<len; i++)
        a[i] = ((u64)rand() << 32) ^ (u64)rand();    
}

void tonelli_shanks_sqrt(const Ec* ecc, const GFElement a, GFElement r) {
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
    }
}

void GFSub(const Ec* ecc, const GFElement a, const GFElement b, GFElement c) {
    GFElement tmp;
    copy(tmp, a, ecc->wordLen);                           
    sub(ecc->wordLen, ecc->p, (u64*)b, (u64*)c);
    GFAdd(ecc, tmp, c, c);
}

/* FIPS-192 Fp: p = 2^192 - 2^64 - 1*/

void GFMul_FIPS192(const Ec* ecc, const GFElement a, const GFElement b, GFElement c) {
    GFElement tmp;
    mul(ecc, a, b, tmp);
    copy(c, tmp, 2*ecc->wordLen);
    u64 carry = 0;

    // S_1
    tmp[0] = c[3];
    tmp[1] = c[3];
    tmp[2] = 0;
    
    carry += add(3, c, tmp, c);
    
    // S_2
    tmp[0] = 0;
    tmp[1] = c[4];
    tmp[2] = c[4];

    carry = add(3, c, tmp, c);
    if (carry) sub(ecc->wordLen, c, ecc->p, c);
    
    // S_3
    tmp[0] = c[5];
    tmp[1] = c[5];
    tmp[2] = c[5];

    carry = add(3, c, tmp, c);
    if (carry) sub(ecc->wordLen, c, ecc->p, c);

    c[3] = 0;

    if(GFCmp(ecc,c,p192)!=-1)
        sub(ecc->wordLen, c, p192, c);
}

void GFSqr_FIPS192(const Ec* ecc, const GFElement a, GFElement b) {
    GFMul_FIPS192(ecc, a, a, b);
}

/* FIPS-224 Fp: p = 2^224 - 2^96 + 1 */

void GFMul_FIPS224(const Ec* ecc, const GFElement a, const GFElement b, GFElement c) {
    GFElement res, tmp,ta,tb;
    copy(ta, a, ecc->wordLen);
    copy(tb, b, ecc->wordLen);
    ta[3] = ta[3] & 0xFFFFFFFF; //truncated a
    tb[3] = tb[3] & 0xFFFFFFFF; //truncated b
    mul(ecc, ta, tb, res);
    copy(c, res, ecc->wordLen);

    u64 carry = 0;
    c[3] = c[3] & 0xFFFFFFFF;

    // S_1
    tmp[0] = 0;
    ((u32*)tmp)[2] = 0;
    ((u32*)tmp)[3] = ((u32*)res)[7];
    ((u32*)tmp)[4] = ((u32*)res)[8];
    ((u32*)tmp)[5] = ((u32*)res)[9];
    ((u32*)tmp)[6] = ((u32*)res)[10];
    ((u32*)tmp)[7] = 0;

    add(ecc->wordLen, tmp, c, c);
    
    // S_2
    ((u32*)tmp)[3] = ((u32*)res)[11];
    ((u32*)tmp)[4] = ((u32*)res)[12];
    ((u32*)tmp)[5] = ((u32*)res)[13];
    ((u32*)tmp)[6] = 0;

    add(ecc->wordLen, tmp, c, c);

    // D_1
    ((u32*)tmp)[0] = ((u32*)res)[11];
    ((u32*)tmp)[1] = ((u32*)res)[12];
    ((u32*)tmp)[2] = ((u32*)res)[13];
    ((u32*)tmp)[3] = 0;
    tmp[2] = 0;
    tmp[3] = 0;
    
    sub(ecc->wordLen, c, tmp, c);

    // D_2
    ((u32*)tmp)[0] = ((u32*)res)[7];
    ((u32*)tmp)[1] = ((u32*)res)[8];
    ((u32*)tmp)[2] = ((u32*)res)[9];
    ((u32*)tmp)[3] = ((u32*)res)[10];
    ((u32*)tmp)[4] = ((u32*)res)[11];
    ((u32*)tmp)[5] = ((u32*)res)[12];
    ((u32*)tmp)[6] = ((u32*)res)[13];

    sub(ecc->wordLen, c, tmp, c);
    
	while(MSB_M & c[3]) { // maybe just 'if' is enough
        add(ecc->wordLen, c, p224, c);
    }
    while(((u32*)c)[7] >= 1) { //and here
        sub(ecc->wordLen, c, p224, c);
    }

    if(GFCmp(ecc,c,p224)!=-1) {
        sub(ecc->wordLen, c, p224, c);
    }
}

void GFSqr_FIPS224(const Ec* ecc, const GFElement a, GFElement b) {
    GFMul_FIPS224(ecc, a, a, b);
}

/*p = 2^256 – 2^224 + 2^192 + 2^96 – 1*/
void GFMul_FIPS256(const Ec* ecc,const GFElement a,const GFElement b,  GFElement c) 
{
    GFElement tmp,res;
    mul(ecc, a, b, res);
    int len = ecc->wordLen; 
    copy(c, res, ecc->wordLen);

    int carry; //signed, its important
    tmp[len] = 0;

    tmp[3] = res[7];
    tmp[2] = res[6];
    ((u32*)tmp)[3] = ((u32*)res)[11];
    ((u32*)tmp)[2] = 0;
    tmp[0] = 0;
    mul2(ecc->wordLen,tmp);
    carry = tmp[len]; //in case 2*tmp > 2^256
    carry+=add(len,c,tmp,c);

    tmp[len] = 0;
    ((u32*)tmp)[7] = 0;
    ((u32*)tmp)[6] = ((u32*)res)[15];
    ((u32*)tmp)[5] = ((u32*)res)[14];
    ((u32*)tmp)[4] = ((u32*)res)[13];
    ((u32*)tmp)[3] = ((u32*)res)[12];
    mul2(ecc->wordLen,tmp);
    carry+=tmp[len];
    carry+=add(len,c,tmp,c);

    tmp[3] = res[7];
    tmp[2] = 0;
    ((u32*)tmp)[3] = 0;
    ((u32*)tmp)[2] = ((u32*)res)[10];
    tmp[0] = res[4];
    carry+=add(len,c,tmp,c);

    ((u32*)tmp)[7] = ((u32*)res)[8];
    ((u32*)tmp)[6] = ((u32*)res)[13];
    tmp[2] = res[7];
    ((u32*)tmp)[3] = ((u32*)res)[13];
    ((u32*)tmp)[2] = ((u32*)res)[11];
    ((u32*)tmp)[1] = ((u32*)res)[10];
    ((u32*)tmp)[0] = ((u32*)res)[9];
    carry+=add(len,c,tmp,c);

    ((u32*)tmp)[7] = ((u32*)res)[10];
    ((u32*)tmp)[6] = ((u32*)res)[8];
    tmp[2] = 0;
    ((u32*)tmp)[3] = 0;
    ((u32*)tmp)[2] = ((u32*)res)[13];
    ((u32*)tmp)[1] = ((u32*)res)[12];
    ((u32*)tmp)[0] = ((u32*)res)[11];
    carry-=sub(len, c, tmp, c);

    ((u32*)tmp)[7] = ((u32*)res)[11];
    ((u32*)tmp)[6] = ((u32*)res)[9];
    tmp[2] = 0;
    tmp[1] = res[7];
    tmp[0] = res[6];
    carry-=sub(len, c, tmp, c);

    ((u32*)tmp)[7] = ((u32*)res)[12];
    ((u32*)tmp)[6] = 0;
    ((u32*)tmp)[5] = ((u32*)res)[10];
    ((u32*)tmp)[4] = ((u32*)res)[9];
    ((u32*)tmp)[3] = ((u32*)res)[8];
    ((u32*)tmp)[2] = ((u32*)res)[15];
    ((u32*)tmp)[1] = ((u32*)res)[14];
    ((u32*)tmp)[0] = ((u32*)res)[13];
    carry-=sub(len, c, tmp, c);

    ((u32*)tmp)[7] = ((u32*)res)[13];
    ((u32*)tmp)[6] = 0;
    tmp[2] = res[5];
    ((u32*)tmp)[3] = ((u32*)res)[9];
    ((u32*)tmp)[2] = 0;
    tmp[0] = res[7];
    carry-=sub(len, c, tmp, c);

    
    while(carry>0) //in case of c > 2^256
    {
        sub(len,c,p256,c);
        carry--;
    }
    while(carry<0) //in case of c<0
    {
        add(len,c,p256,c);
        carry++;
    }
    
    if(GFCmp(ecc,c,p256)!=-1)
        sub(len, c, p256, c);
}

void GFSqr_FIPS256(const Ec* ecc,const GFElement a,  GFElement  b) {
    GFMul_FIPS256(ecc, a, a, b);
}

/* FIPS-384 Fp: p = p^384 - 2^128 - 2^96 +  2^32 - 1 */
void GFMul_FIPS384(const Ec* ecc, const GFElement a, const GFElement b, GFElement c) {
    GFElement res, tmp;
    mul(ecc, a, b, res);
    copy(c, res, 2*ecc->wordLen);

    int carry = 0;

    // 2*S_1
    tmp[0] = 0;
    tmp[1] = 0;
    ((u32*)tmp)[4] = ((u32*)res)[21];
    ((u32*)tmp)[5] = ((u32*)res)[22];
    ((u32*)tmp)[6] = ((u32*)res)[23];
    ((u32*)tmp)[7] = 0;
    tmp[4] = 0;
    tmp[5] = 0;

    mul2(ecc->wordLen, tmp);
    carry = tmp[ecc->wordLen];
    carry+=add(ecc->wordLen, tmp, c, c);

    // S_2
    ((u32*)tmp)[0] = ((u32*)res)[12];
    ((u32*)tmp)[1] = ((u32*)res)[13];
    ((u32*)tmp)[2] = ((u32*)res)[14];
    ((u32*)tmp)[3] = ((u32*)res)[15];
    ((u32*)tmp)[4] = ((u32*)res)[16];
    ((u32*)tmp)[5] = ((u32*)res)[17];
    ((u32*)tmp)[6] = ((u32*)res)[18];
    ((u32*)tmp)[7] = ((u32*)res)[19];
    ((u32*)tmp)[8] = ((u32*)res)[20];
    ((u32*)tmp)[9] = ((u32*)res)[21];
    ((u32*)tmp)[10] = ((u32*)res)[22];
    ((u32*)tmp)[11] = ((u32*)res)[23];

    carry+= add(ecc->wordLen, tmp, c, c);

    // S_3
    ((u32*)tmp)[0] = ((u32*)res)[21];
    ((u32*)tmp)[1] = ((u32*)res)[22];
    ((u32*)tmp)[2] = ((u32*)res)[23];
    ((u32*)tmp)[3] = ((u32*)res)[12];
    ((u32*)tmp)[4] = ((u32*)res)[13];
    ((u32*)tmp)[5] = ((u32*)res)[14];
    ((u32*)tmp)[6] = ((u32*)res)[15];
    ((u32*)tmp)[7] = ((u32*)res)[16];
    ((u32*)tmp)[8] = ((u32*)res)[17];
    ((u32*)tmp)[9] = ((u32*)res)[18];
    ((u32*)tmp)[10] = ((u32*)res)[19];
    ((u32*)tmp)[11] = ((u32*)res)[20];

    carry += add(ecc->wordLen, c, tmp, c);

    // S_4
    ((u32*)tmp)[0] = 0;
    ((u32*)tmp)[1] = ((u32*)res)[23];
    ((u32*)tmp)[2] = 0;
    ((u32*)tmp)[3] = ((u32*)res)[20];
    ((u32*)tmp)[4] = ((u32*)res)[12];
    ((u32*)tmp)[5] = ((u32*)res)[13];
    ((u32*)tmp)[6] = ((u32*)res)[14];
    ((u32*)tmp)[7] = ((u32*)res)[15];
    ((u32*)tmp)[8] = ((u32*)res)[16];
    ((u32*)tmp)[9] = ((u32*)res)[17];
    ((u32*)tmp)[10] = ((u32*)res)[18];
    ((u32*)tmp)[11] = ((u32*)res)[19];

    carry += add(ecc->wordLen, c, tmp, c);

    // S_5
    tmp[0] = 0;
    tmp[1] = 0;
    ((u32*)tmp)[4] = ((u32*)res)[20];
    ((u32*)tmp)[5] = ((u32*)res)[21];
    ((u32*)tmp)[6] = ((u32*)res)[22];
    ((u32*)tmp)[7] = ((u32*)res)[23];
    tmp[4] = 0;
    tmp[5] = 0;

    carry+=add(ecc->wordLen, c, tmp, c);

    // S_6
    ((u32*)tmp)[0] = ((u32*)res)[20];
    ((u32*)tmp)[1] = 0;
    ((u32*)tmp)[2] = 0;
    ((u32*)tmp)[3] = ((u32*)res)[21];
    ((u32*)tmp)[4] = ((u32*)res)[22];
    ((u32*)tmp)[5] = ((u32*)res)[23];	
    tmp[3] = 0;
    tmp[4] = 0;
    tmp[5] = 0;

    carry+=add(ecc->wordLen, c, tmp, c);

    // D_1
    ((u32*)tmp)[0] = ((u32*)res)[23];
    ((u32*)tmp)[1] = ((u32*)res)[12];
    ((u32*)tmp)[2] = ((u32*)res)[13];
    ((u32*)tmp)[3] = ((u32*)res)[14];
    ((u32*)tmp)[4] = ((u32*)res)[15];
    ((u32*)tmp)[5] = ((u32*)res)[16];
    ((u32*)tmp)[6] = ((u32*)res)[17];
    ((u32*)tmp)[7] = ((u32*)res)[18];
    ((u32*)tmp)[8] = ((u32*)res)[19];
    ((u32*)tmp)[9] = ((u32*)res)[20];
    ((u32*)tmp)[10] = ((u32*)res)[21];
    ((u32*)tmp)[11] = ((u32*)res)[22];

    carry-=sub(ecc->wordLen, c, tmp, c);

    // D_2
    ((u32*)tmp)[0] = 0;
    ((u32*)tmp)[1] = ((u32*)res)[20];
    ((u32*)tmp)[2] = ((u32*)res)[21];
    ((u32*)tmp)[3] = ((u32*)res)[22];
    ((u32*)tmp)[4] = ((u32*)res)[23];
    ((u32*)tmp)[5] = 0;
    tmp[3] = 0;
    tmp[4] = 0;
    tmp[5] = 0;
    
    carry-=sub(ecc->wordLen, c, tmp, c);

    // D_3
    tmp[0] = 0;
    ((u32*)tmp)[2] = 0;
    ((u32*)tmp)[3] = ((u32*)res)[23];
    ((u32*)tmp)[4] = ((u32*)res)[23];
    ((u32*)tmp)[5] = 0;
    tmp[3] = 0;
    tmp[4] = 0;
    tmp[5] = 0;
    
    carry-=sub(ecc->wordLen, c, tmp, c);

    while(carry>0) //in case of c > 2^256
    {
        sub(ecc->wordLen,c,ecc->p,c);
        carry--;
    }
    while(carry<0) //in case of c<0
    {
        add(ecc->wordLen,c,ecc->p,c);
        carry++;
    }
    //in case of 0<c<p256
    if(GFCmp(ecc,c,ecc->p)!=-1)
        sub(ecc->wordLen, c, ecc->p, c);
}

void GFSqr_FIPS384(const Ec* ecc, const GFElement a, GFElement b) {
    GFMul_FIPS384(ecc, a, a, b);
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
        carry = b[ ecc->wordLen -1 ] & ((u64)(-1)<<32);
        while(carry) {
            sub(ecc->wordLen, b, ecc->p, b);
            carry = b[ ecc->wordLen -1 ] & ((u64)(-1)<<32);
        }
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
    GFPow(ecc, a, ecc->p_min_two, b);
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