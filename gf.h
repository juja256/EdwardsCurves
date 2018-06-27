#ifndef GF_H
#define GF_H

#include "ec.h"

#ifdef __cplusplus
extern "C" {
#endif


extern const u64 unity[9];
extern const EcPoint uPEd;
extern const EcPoint uPW;
extern const EcPointProj uPPEd;
extern const EcPointProj uPPW;
extern const u64 p192[3];
extern const u64 p224[4];
extern const u64 p256[4];
extern const u64 p384[6];
extern const u64 p521[9];
extern const u64 zero[20];

/* Common Arithmetics */

void shr(u64 n, const u64* a, u64* res, u64 bits);
#define div2(n, a) shr((n), (a), (a), 1)
void shl(u64 n, const u64* a, u64* res, u64 bits);
#define mul2(n, a) shl((n), (a), (a), 1)
u64 get_bit(const u64* a, u64 num);
void copy(u64* a, const u64* b, int len);
u64 add(u64 n, const u64* a, const u64* b, u64* c);
u64 sub(u64 n, const u64* a, const u64* b, u64* c);
void _mul_raw(u64 a, u64 b, u64* low, u64* high);
u64 _add_raw(u64 a, u64 b, u64* c);
void mul_by_word(u64 n, const u64* a, u64 d, u64* c);
void mul(u64 n, const u64* a, const u64* b, u64* c);
void imul(u64 n, const u64* a, const u64* b, u64* c);

int word_bit_len(u64 w);
int bigint_bit_len(u64 n, const u64* a);

void divide(u64 n, const u64* a, const u64* b, u64* quotient, u64* reminder);
int cmp(u64 n, const u64* a, const u64* b);


void add_mod(u64 n, const BigInt a, const BigInt b, const BigInt m, BigInt res);
void mul_mod(u64 n, const BigInt a, const BigInt b, const BigInt m, BigInt res);
void exp_mod(u64 n, const BigInt a, const BigInt b, const BigInt m, BigInt res);
void inv_mod(u64 n, const BigInt a, const BigInt m, BigInt res);


/* Galois' Fields Arithmetics  */
void GFSqrt(const Ec* ecc, const GFElement a, GFElement r); // via Tonelli-Shanks
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
#define GFMulBy2(e, a, b) GFAdd(e, a, a, b);
void GFMulByD(const EcEd* ecc, GFElement a);

void GFMul_Cmn(const Ec* ecc, const GFElement a, const GFElement b, GFElement c);
void GFMul_FIPS192(const Ec* ecc, const GFElement a, const GFElement b, GFElement c);
void GFMul_FIPS224(const Ec* ecc, const GFElement a, const GFElement b, GFElement c);
void GFMul_FIPS256(const Ec* ecc, const GFElement a, const GFElement b, GFElement c);
void GFMul_FIPS384(const Ec* ecc, const GFElement a, const GFElement b, GFElement c);
void GFMul_FIPS521(const Ec* ecc, const GFElement a, const GFElement b, GFElement c);

void GFSqr_Cmn(const Ec* ecc, const GFElement a, GFElement c);
void GFSqr_FIPS192(const Ec* ecc, const GFElement a, GFElement c);
void GFSqr_FIPS224(const Ec* ecc, const GFElement a, GFElement c);
void GFSqr_FIPS256(const Ec* ecc, const GFElement a, GFElement c);
void GFSqr_FIPS384(const Ec* ecc, const GFElement a, GFElement c);
void GFSqr_FIPS521(const Ec* ecc, const GFElement a, GFElement c);

#ifdef __cplusplus
}
#endif

#endif /* GF_H */
