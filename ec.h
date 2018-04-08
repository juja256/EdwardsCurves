#ifndef EC_H
#define EC_H

#define NORMAL_POINT 1
#define INFINITY_POINT 0

typedef unsigned long long u64;
typedef unsigned u32;
typedef unsigned char u8;

typedef u64 GFElement[10];
typedef u64 BigInt[10]; // up to 640 bits
typedef u64 VeryBigInt[18]; // up to 1152 bits

typedef struct _EcPoint {
    GFElement x;
    GFElement y;
} EcPoint;

typedef struct _EcPointProj {
    GFElement X;
    GFElement Y;
    GFElement Z;
} EcPointProj;

typedef int BOOL;
typedef struct _Ec Ec;

typedef void TGFAddFunc(const Ec*, const GFElement, const GFElement, GFElement);
typedef void TGFMulFunc(const Ec*, const GFElement, const GFElement, GFElement);
typedef void TGFSubFunc(const Ec*, const GFElement, const GFElement, GFElement);
typedef void TGFSqrFunc(const Ec*, const GFElement, GFElement);

#define PRNG_STATE_LEN 256 
typedef struct {
    unsigned char state[PRNG_STATE_LEN];
} PRNG;

void PRNGInit(PRNG* generator, unsigned char* seed, int len);
void PRNGRun(PRNG* generator); 
void PRNGGenerateSequence(PRNG* generator, int bit_len, unsigned char* dest);

/* Elliptic curve in Edwards:
    x^2 + y^2 = 1 + d*x^2*y^2
   or in Weierstrass form:
    y^2 = x^3 + a*x + b 
*/
typedef struct _Ec {
    GFElement d, a, b;
    BOOL isEdwards;

    u64 bitLen;
    u64 wordLen;

    BigInt n;
    u64 cofactor;
    EcPoint BasePoint;

    BigInt p;
    BigInt p_min_two;
    
    TGFMulFunc* GFMul;
    TGFSqrFunc* GFSqr;
    PRNG prng;
} Ec;

typedef Ec EcEd;
typedef Ec EcW;


int  EcEdInit(EcEd* ecc, u64 bitLen, const BigInt p, const EcPoint* bp, const BigInt n, const GFElement d);
int  EcWInit(EcW* ecc, u64 bitLen, const BigInt p, const EcPoint* bp, const BigInt n, const GFElement a, const GFElement b);

int  EcInitStandardCurve(Ec* ecc, u64 bitLen, BOOL isEdwards);

void EcGenerateBasePoint(Ec* ecc, EcPoint* bp);

int  EcCheckPointInMainSubGroup(Ec* ecc, const EcPoint* P);
int  EcCheckPointOnCurve(Ec* ecc, const EcPoint* P);

void EcDump(Ec* ecc, char* buf);

/* Arithmetic on Elliptic Curves */
int  EcPointCmp(Ec* ecc, const EcPoint* A, const EcPoint* B);
void EcCopy(Ec* ecc, EcPoint* dest, const EcPoint* src);
void EcCopyProj(Ec* ecc, EcPointProj* dest, const EcPointProj* src);

void EcConvertAffineToProjective(Ec* ecc, const EcPoint* P, EcPointProj* Q);
void EcConvertProjectiveToAffine(Ec* ecc, const EcPointProj* P, EcPoint* Q);

int EcAdd(Ec* ecc, const EcPoint* A, const EcPoint* B, EcPoint* C);
int EcDouble(Ec* ecc, const EcPoint* A, EcPoint* B);

int EcAddProj(Ec* ecc, const EcPointProj* A, const EcPointProj* B, EcPointProj* C);
int EcDoubleProj(Ec* ecc, const EcPointProj* A, EcPointProj* B);

int EcScalarMulProj(Ec* ecc, const EcPointProj* A, const BigInt k, EcPointProj* B);
int EcScalarMul(Ec* ecc, const EcPoint* A, const BigInt k, EcPoint* B);


#endif /* EC_H */