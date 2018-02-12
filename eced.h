#ifndef ECED_H
#define ECED_H

#define _out_
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

typedef void TGFAddFunc(EcEd*, GFElement, GFElement, GFElement);
typedef void TGFMulFunc(EcEd*, GFElement, GFElement, GFElement);
typedef void TGFSubFunc(EcEd*, GFElement, GFElement, GFElement);
typedef void TGFSqrFunc(EcEd*, GFElement, GFElement);

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
void GFDump(EcEd* ecc, GFElement a);
void GFAdd(EcEd* ecc, GFElement a, GFElement b, _out_ GFElement c);
void GFSub(EcEd* ecc, GFElement a, GFElement b, _out_ GFElement c);
void GFPow(EcEd* ecc, GFElement a, BigInt n, GFElement b);
void GFInv(EcEd* ecc, GFElement a, GFElement b);
int  GFCmp(EcEd* ecc, GFElement a, GFElement b);
void GFMul(EcEd* ecc, GFElement a, GFElement b, _out_ GFElement c);
void GFSqr(EcEd* ecc, GFElement a, _out_ GFElement c);

int  EcEdInit(EcEd* ecc, EcPoint* bp, u64 bitLen, BigInt n, GFElement d);
void EcEdGenerateBasePoint(EcEd* ecc, EcPoint* bp);
int  EcEdCheckPointOnCurve(EcEd* ecc,EcPoint*);

void EcEdAdd(EcEd* ecc, EcPoint* A, EcPoint* B, _out_ EcPoint* C);
void EcEdDouble(EcEd* ecc, EcPoint* A, _out_ EcPoint* B);

void EcEdAddProj(EcEd* ecc, EcPointProj* A, EcPointProj* B, _out_ EcPointProj* C);
void EcEdDoubleProj(EcEd* ecc, EcPointProj* A, _out_ EcPointProj* B);

void EcEdScalarMul(EcEd* ecc, EcPoint* A, BigInt k, _out_ EcPoint* B);
void EcEdScalarMulOrdinary(EcEd* ecc, EcPoint* A, BigInt k, _out_ EcPoint* B);


#endif