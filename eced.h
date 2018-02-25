#ifndef ECED_H
#define ECED_H

#ifdef _WIN64
#include <intrin.h>
#pragma intrinsic(_umul128) 
#else
#include <x86intrin.h>
#endif // _WIN64

typedef unsigned long long u64;
typedef unsigned u32;
typedef unsigned char u8;

typedef u64 GFElement[12];
typedef u64 BigInt[6];

typedef struct _EcPoint {
    GFElement x;
    GFElement y;
} EcPoint;

typedef struct _EcPointProj {
    GFElement X;
    GFElement Y;
    GFElement Z;
} EcPointProj;

typedef struct _EcEd EcEd;

typedef void TGFAddFunc(const EcEd*, const GFElement, const GFElement, GFElement);
typedef void TGFMulFunc(const EcEd*, const GFElement, const GFElement, GFElement);
typedef void TGFSubFunc(const EcEd*, const GFElement, const GFElement, GFElement);
typedef void TGFSqrFunc(const EcEd*, const GFElement, GFElement);

typedef struct _EcEd {
    GFElement d;
    BigInt n;
    EcPoint BasePoint; 
    BigInt p;
    BigInt p_min_two;
    u64 bitLen;
    u64 wordLen;
    TGFMulFunc* GFMul;
    TGFSqrFunc* GFSqr;
} EcEd;

void GFInitFromString(GFElement a, const char* str);
void GFDump(const EcEd* ecc, const GFElement a);
void GFAdd(const EcEd* ecc, const GFElement a, const GFElement b, GFElement c);
void GFSub(const EcEd* ecc, const GFElement a, const GFElement b, GFElement c);
void GFPow(const EcEd* ecc, const GFElement a, const BigInt n, GFElement b);
void GFInv(const EcEd* ecc, const GFElement a, GFElement b);
int  GFCmp(const EcEd* ecc, const GFElement a, const GFElement b);
void GFMul(const EcEd* ecc, const GFElement a, const GFElement b, GFElement c);
void GFSqr(const EcEd* ecc, const GFElement a, GFElement c);

int  EcEdInit(EcEd* ecc, const EcPoint* bp, u64 bitLen, const BigInt n, const GFElement d);
void EcEdGenerateBasePoint(const EcEd* ecc, EcPoint* bp);
int  EcEdCheckPointOnCurve(const EcEd* ecc, const EcPoint* P);

void EcEdAdd(const EcEd* ecc, const EcPoint* A, const EcPoint* B, EcPoint* C);
void EcEdDouble(const EcEd* ecc, const EcPoint* A, EcPoint* B);

void EcEdAddProj(const EcEd* ecc, const EcPointProj* A, const EcPointProj* B, EcPointProj* C);
void EcEdDoubleProj(const EcEd* ecc, const EcPointProj* A, EcPointProj* B);

void EcEdScalarMul(const EcEd* ecc, const EcPoint* A, const BigInt k, EcPoint* B);
void EcEdScalarMulOrdinary(const EcEd* ecc, EcPoint* A, const BigInt k, EcPoint* B);


#endif /* ECED_H */
