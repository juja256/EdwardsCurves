#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

#include "const.h"


static inline u64 get_bit(const GFElement a, u64 num) {
    return a[num/64] & ((u64)1 << (num % 64));
}

static inline void copy(GFElement a, const GFElement b, int len) {
    for (int i=0;i<len;i++)
        a[i] = b[i];
}

static inline void copy_point(EcPoint* P, const EcPoint* Q, int len) {
    copy(P->x, Q->x, len);
    copy(P->y, Q->y, len);
}

static inline u64 add(u64 n, const u64* a, const u64* b, u64* c) {
    u64 msb_a, msb_b, carry = 0;

    for (u32 i=0;i < n; i++) {
        msb_a = a[i] & MSB_M;
        msb_b = b[i] & MSB_M;
        c[i] = a[i] + b[i] + carry;
        carry = ( (msb_a && msb_b) || ((msb_a ^ msb_b) && !(MSB_M & c[i])) );
    }
    return carry;
}

static inline u64 sub(u64 n, const u64* a, const u64* b, u64* c) {
    u64 borrow = 0;
    for (int i=0; i<n; i++) {
        u64 t_a = a[i];
        u64 t_b = b[i];
        c[i] = t_a - t_b - borrow;
        borrow = ( (~t_a) & (c[i] | t_b) | (c[i] & t_b) ) >> (63);
    }
    return borrow;
}

static inline void _mul_raw(u64 a, u64 b, u64* low, u64* high) {
#ifdef _WIN64
	*low = _umul128(a, b, high);
#else
	__int128 r = (__int128)a * (__int128)b;
	*low = (unsigned long long)r;
	*high = r >> 64;
#endif // _WIN64
}

static inline u64 _add_raw(u64 a, u64 b, u64* c) {
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

static inline void mul_by_word(const EcEd* ecc, const u64* a, u64 d, u64* c) {
    u64 carry = 0, carry_tmp;
    for (int i=0; i < ecc->wordLen; i++) {
        carry_tmp = carry;
        _mul_raw(d, a[i], &(c[i]), &carry);
        carry += _add_raw(c[i], carry_tmp, &(c[i]));
    }
    c[ecc->wordLen] = carry;
}

static inline void mul(const EcEd* ecc, const u64* a, const u64* b, u64* c) {
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
    for (int i = 0; i < n; i++) {
        u64 cur = a[i];
        a[i] <<= 1;
        a[i] ^= buf;
        buf = (cur & (MAX_U64 << ( 63 ))) >> ( 63 );
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

static inline void randomize(u64 len, GFElement a) {
    for (u64 i=0; i<len; i++)
        a[i] = ((u64)rand() << 32) ^ (u64)rand();    
}

static inline void tonelli_shanks_sqrt(const EcEd* ecc, const GFElement a, GFElement r) {
    // p = Q*2^s
    GFElement Q, pp, z, tmp, c, t, R, b;
    copy(z, unity, ecc->wordLen);
    z[0]++;
    copy(Q, ecc->p, ecc->wordLen);
    u64 M, S = 1;
    div2(ecc->wordLen, Q);
    copy(pp, Q, ecc->wordLen);

    while ((Q[0] & 1) == 0) {
        div2(ecc->wordLen, Q);
        S++;
    }

    while(1) {
        GFPow(ecc, z, pp, tmp); 
        if (!GFCmp(ecc, tmp, unity))
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
    GFPow(ecc, a, tmp, R);

    while (1) {
        if (!GFCmp(ecc, t, unity)) {
            copy(r, R, ecc->wordLen);
            break;
        } 
        for (i=1; i<S; i++) {
            GFSqr(ecc, t, t);
            if (!GFCmp(ecc, t, unity));
                break;
        }
        copy(b, c, ecc->wordLen);
        for (u64 j=0; j<M-i; j++) {
            GFSqr(ecc, b, b);
        }
        M = i;
        GFSqr(ecc, b, c);
        GFMul(ecc, t, c, t);
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

void GFDump(const EcEd* ecc, const GFElement a) {
    for (int i=ecc->wordLen-1; i>=0; i--)
        printf(HEX_FORMAT, a[i]);
    printf("\n");
}


void GFNeg(const EcEd* ecc, const GFElement a, GFElement c) {
    sub(ecc->wordLen, ecc->p, (u64*)a, (u64*)c);
}

void GFAdd(const EcEd* ecc, const GFElement a, const GFElement b, GFElement c) {
    u64 carry = 0;
    
	if (ecc->bitLen == 192 || ecc->bitLen == 256 || ecc->bitLen == 384) {
		carry = add(ecc->wordLen, a, b, c);
		if (carry) {
			sub(ecc->wordLen, c, ecc->p, c);
		}
		if (GFCmp(ecc, c, ecc->p) != -1) {
			sub(ecc->wordLen, c, ecc->p, c);
		}
	}
	else if (ecc->bitLen == 224) {
		add(ecc->wordLen, a, b, c);
		carry = (c[ecc->wordLen - 1] & ((u64)1 << 32));
		if (carry) {
			sub(ecc->wordLen, c, ecc->p, c);
		}
		if (GFCmp(ecc, c, ecc->p) == 1) {
			sub(ecc->wordLen, c, ecc->p, c);
		}
	}
}

void GFSub(const EcEd* ecc, const GFElement a, const GFElement b, GFElement c) {
    GFElement tmp;
    copy(tmp, a, ecc->wordLen);                           
    sub(ecc->wordLen, ecc->p, (u64*)b, (u64*)c);
    GFAdd(ecc, tmp, c, c);
}

/* FIPS-192 Fp: p = 2^192 - 2^64 - 1*/

void GFMul_FIPS192(const EcEd* ecc, const GFElement a, const GFElement b, GFElement c) {
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
    if ((c[2] == MAX_U64) && (c[1] == MAX_U64)) {
        sub(ecc->wordLen, c, ecc->p, c);
    }
    c[3] = 0;
}

void GFSqr_FIPS192(const EcEd* ecc, const GFElement a, GFElement b) {
    GFMul_FIPS192(ecc, a, a, b);
}

/* FIPS-224 Fp: p = 2^224 - 2^96 - 1 */

void GFMul_FIPS224(const EcEd* ecc, const GFElement a, const GFElement b, GFElement c) {
    GFElement res, tmp;
    mul(ecc, a, b, res);
    copy(c, res, ecc->wordLen);

    u64 carry = 0;
    c[ecc->wordLen-1] = c[ecc->wordLen-1] & 0xFFFFFFFF;

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
	((u32*)tmp)[0] = ((u32*)res)[7];
	((u32*)tmp)[1] = ((u32*)res)[8];
	((u32*)tmp)[2] = ((u32*)res)[9];
	((u32*)tmp)[3] = ((u32*)res)[10];
	((u32*)tmp)[4] = ((u32*)res)[11];
	((u32*)tmp)[5] = ((u32*)res)[12];
	((u32*)tmp)[6] = ((u32*)res)[13];

	sub(ecc->wordLen, c, tmp, c);

	// D_2
    ((u32*)tmp)[0] = ((u32*)res)[11];
    ((u32*)tmp)[1] = ((u32*)res)[12];
    ((u32*)tmp)[2] = ((u32*)res)[13];
    ((u32*)tmp)[3] = 0;
    tmp[2] = 0;
    tmp[3] = 0;
    
	sub(ecc->wordLen, c, tmp, c);

    
	if (MSB_M & c[3]) {
        add(ecc->wordLen, c, ecc->p, c);
    }
    else if (((u32*)c)[7] == 1) {
        sub(ecc->wordLen, c, ecc->p, c);
    }
	if (GFCmp(ecc, c, p224) != -1) {
		sub(ecc->wordLen, c, ecc->p, c);
	}
}

void GFSqr_FIPS224(const EcEd* ecc, const GFElement a, GFElement b) {
    GFMul_FIPS224(ecc, a, a, b);
}

/* FIPS-256 Fp: p = 2^256 - 2^224 + 2^192 + 2^96 - 1 */

void GFMul_FIPS256(const EcEd* ecc, const GFElement a, const GFElement b, GFElement c) {

}

void GFSqr_FIPS256(const EcEd* ecc, const GFElement a, GFElement b) {
	GFMul_FIPS256(ecc, a, a, b);
}

/* FIPS-384 Fp: p = p^384 - 2^128 - 2^96 +  2^32 - 1 */

void GFMul_FIPS384(const EcEd* ecc, const GFElement a, const GFElement b, GFElement c) {
	GFElement res, tmp;
	mul(ecc, a, b, res);
	copy(c, res, 2*ecc->wordLen);

	u64 carry = 0;

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
	add(ecc->wordLen, tmp, c, c);
	if (carry) sub(ecc->wordLen, c, ecc->p, c);

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

	carry = add(ecc->wordLen, tmp, c, c);
	if (carry) sub(ecc->wordLen, c, ecc->p, c);

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

	carry = add(ecc->wordLen, c, tmp, c);
	if (carry) sub(ecc->wordLen, c, ecc->p, c);


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

	carry = add(ecc->wordLen, c, tmp, c);
	if (carry) sub(ecc->wordLen, c, ecc->p, c);

	// S_5
	tmp[0] = 0;
	tmp[1] = 0;
	((u32*)tmp)[4] = ((u32*)res)[20];
	((u32*)tmp)[5] = ((u32*)res)[21];
	((u32*)tmp)[6] = ((u32*)res)[22];
	((u32*)tmp)[7] = ((u32*)res)[23];
	tmp[4] = 0;
	tmp[5] = 0;

	add(ecc->wordLen, c, tmp, c);
	if (carry) sub(ecc->wordLen, c, ecc->p, c);

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

	add(ecc->wordLen, c, tmp, c);
	if (carry) sub(ecc->wordLen, c, ecc->p, c);

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

	sub(ecc->wordLen, c, tmp, c);

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
	
	sub(ecc->wordLen, c, tmp, c);

	// D_3
	tmp[0] = 0;
	((u32*)tmp)[2] = 0;
	((u32*)tmp)[3] = ((u32*)res)[23];
	((u32*)tmp)[4] = ((u32*)res)[23];
	((u32*)tmp)[5] = 0;
	tmp[3] = 0;
	tmp[4] = 0;
	tmp[5] = 0;
	
	sub(ecc->wordLen, c, tmp, c);

	if (MSB_M & c[5]) {
		add(6, c, p384, c);
	}
	/*else if (((u32*)c)[7] == 1) {
		sub(ecc->wordLen, c, p384, c);
	}*/
}

void GFSqr_FIPS384(const EcEd* ecc, const GFElement a, GFElement b) {
	GFMul_FIPS384(ecc, a, a, b);
}

void GFMulBy2(const EcEd* ecc, const GFElement a, GFElement b) {
    u64 carry = 0;
    copy(b, a, ecc->wordLen);
    b[ecc->wordLen] = 0;
    mul2(ecc->wordLen, b);
    if (ecc->bitLen % 64 == 0) {
        carry = b[ecc->wordLen] & 1;
    }
    else {
        carry = (b[ ecc->wordLen -1 ] & ((u64)1<<32) );
    }
    if (carry) {
        sub(ecc->wordLen, b, ecc->p, b);
    }
}

void GFPow(const EcEd* ecc, const GFElement a, const BigInt n, GFElement b) {
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

void GFInv(const EcEd* ecc, const GFElement a, GFElement b) {
    GFPow(ecc, a, ecc->p_min_two, b);
}

void GFMul(const EcEd* ecc, const GFElement a, const GFElement b, GFElement c) {
    ecc->GFMul(ecc, a, b, c);
}

void GFSqr(const EcEd* ecc, const GFElement a, GFElement c) {
    ecc->GFSqr(ecc, a, c);
}

int GFCmp(const EcEd* ecc, const GFElement a, const GFElement b) {
    for (int i=ecc->wordLen - 1; i>=0; i--) {
		if (a[i] > b[i]) return 1;
		if (a[i] < b[i]) return -1;
    }
    return 0;
}

int EcEdInit(EcEd* ecc, const EcPoint* bp, u64 bitLen, const BigInt n, const GFElement d) {
    srand(time(NULL));
    if ( (bitLen != 192) && (bitLen != 224) && (bitLen != 256) && (bitLen != 384)) {
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
		case 256:

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
    return !GFCmp(ecc, x, z);
}

/* x^2 + y^2 = 1 + dx^2y^2 */
/* y^2 = (1 - x^2)/(1 - dx^2) */
void EcEdGenerateBasePoint(const EcEd* ecc, EcPoint* bp) {
    randomize(ecc->wordLen, bp->x);
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

void EcEdScalarMulOrdinary(const EcEd* ecc, const EcPoint* A, const BigInt k, EcPoint* B) {
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

void EcEdConvertAffineToProjective(const EcEd* ecc, const EcPoint* P, EcPointProj* Q) {
    randomize(ecc->wordLen, Q->Z);
    GFMul(ecc, P->x, Q->Z, Q->X);
    GFMul(ecc, P->y, Q->Z, Q->Y);
}

void EcEdConvertProjectiveToAffine(const EcEd* ecc, const EcPointProj* P, EcPoint* Q) {
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

void EcEdScalarMul(const EcEd* ecc, const EcPoint* A, const BigInt k, EcPoint* B) {
    EcPointProj P, H;
    EcEdConvertAffineToProjective(ecc, &uP, &P);
    EcEdConvertAffineToProjective(ecc, A, &H);

    for (u32 i=0; i<ecc->bitLen; i++) {
        if (get_bit(k, i)) {
            EcEdAddProj(ecc, &P, &H, &P);
        }
        EcEdDoubleProj(ecc, &H, &H); 
    }
    EcEdConvertProjectiveToAffine(ecc, &P, B);
}
