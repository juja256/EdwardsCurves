#ifndef GF_H
#define GF_H

#include "ec.h"

extern const u64 unity[6];
extern const EcPoint uP;
extern const u64 p192[3];
extern const u64 p224[4];
extern const u64 p256[4];
extern const u64 p384[6];

/* Common Arithmetics */
void tonelli_shanks_sqrt(const Ec* ecc, const GFElement a, GFElement r);
void randomize(u64 len, GFElement a);
void shr(u64 n, const GFElement a, GFElement res, u64 bits);
#define div2(n, a) shr((n), (a), (a), 1)
void shl(u64 n, const GFElement a, GFElement res, u64 bits);
#define mul2(n, a) shl((n), (a), (a), 1)
u64 get_bit(const GFElement a, u64 num);
void copy(GFElement a, const GFElement b, int len);
void copy_point(EcPoint* P, const EcPoint* Q, int len);
u64 add(u64 n, const u64* a, const u64* b, u64* c);
u64 sub(u64 n, const u64* a, const u64* b, u64* c);
void _mul_raw(u64 a, u64 b, u64* low, u64* high);
u64 _add_raw(u64 a, u64 b, u64* c);
void mul_by_word(const Ec* ecc, const u64* a, u64 d, u64* c);
void mul(const Ec* ecc, const u64* a, const u64* b, u64* c);
int word_bit_len(u64 n);
int bigint_bit_len(u64 nWords, const BigInt a);
void basic_reduction(u64 n, const BigInt a, const BigInt p, BigInt res);

/* Galois' Fields Arithmetics  */
void GFInitFromString(GFElement a, const char* str);
void GFDump(const Ec* ecc, const GFElement a);
void GFAdd(const Ec* ecc, const GFElement a, const GFElement b, GFElement c);
void GFSub(const Ec* ecc, const GFElement a, const GFElement b, GFElement c);
void GFNeg(const Ec* ecc, const GFElement a, GFElement c);
void GFPow(const Ec* ecc, const GFElement a, const BigInt n, GFElement b);
void GFInv(const Ec* ecc, const GFElement a, GFElement b);
int  GFCmp(const Ec* ecc, const GFElement a, const GFElement b);
void GFMul(const Ec* ecc, const GFElement a, const GFElement b, GFElement c);
void GFSqr(const Ec* ecc, const GFElement a, GFElement c);
void GFMulBy2Power(const Ec* ecc, const GFElement a, int pp, GFElement b);
#define GFMulBy2(e, a, b) GFMulBy2Power(e, a, 1, b);

void GFMul_FIPS192(const Ec* ecc, const GFElement a, const GFElement b, GFElement c);
void GFMul_FIPS224(const Ec* ecc, const GFElement a, const GFElement b, GFElement c);
void GFMul_FIPS256(const Ec* ecc, const GFElement a, const GFElement b, GFElement c);
void GFMul_FIPS384(const Ec* ecc, const GFElement a, const GFElement b, GFElement c);

void GFSqr_FIPS192(const Ec* ecc, const GFElement a, GFElement c);
void GFSqr_FIPS224(const Ec* ecc, const GFElement a, GFElement c);
void GFSqr_FIPS256(const Ec* ecc, const GFElement a, GFElement c);
void GFSqr_FIPS384(const Ec* ecc, const GFElement a, GFElement c);


#endif /* GF_H */
