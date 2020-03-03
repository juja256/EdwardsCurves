#ifndef EC_H
#define EC_H

#ifdef __cplusplus
extern "C" {
#endif

#define NORMAL_POINT 1
#define INFINITY_POINT 0

typedef unsigned long long u64;
typedef unsigned u32;
typedef unsigned char u8;

/* actually 521bit number needs only 9 words, but there is one additional word for signed numbers support */
typedef u64 GFElement[13]; 
typedef u64 BigInt[13]; // up to 832 bits
typedef u64 VeryBigInt[26]; // up to 2*832 bits

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
#define TRUE 1
#define FALSE 0

typedef struct _Ec Ec;

typedef void TGFAddFunc(const Ec*, const GFElement, const GFElement, GFElement);
typedef void TGFMulFunc(const Ec*, const GFElement, const GFElement, GFElement);
typedef void TGFSubFunc(const Ec*, const GFElement, const GFElement, GFElement);
typedef void TGFSqrFunc(const Ec*, const GFElement, GFElement);

typedef void TEcAdd(Ec*, const EcPointProj*, const EcPointProj*, EcPointProj*);
typedef void TEcDouble(Ec*, const EcPointProj*, EcPointProj*);

#define PRNG_STATE_LEN 256 
typedef struct {
    unsigned char state[PRNG_STATE_LEN];
} PRNG;

void PRNGInit(PRNG* generator, unsigned char* seed, int len);
void PRNGRun(PRNG* generator); 
void PRNGGenerateSequence(PRNG* generator, int bit_len, unsigned char* dest);

#define INVALID_D -2
#define INVALID_A -3

/* Elliptic curve in Edwards:
    x^2 + a*y^2 = 1 + d*x^2*y^2
   or in Weierstrass form:
    y^2 = x^3 + a*x + b 
*/

#define ED_192 0xE192
#define ED_224 0xE224
#define ED_256 0xE256
#define ED_384 0xE384
#define ED_521 0xE521
#define ED_NOT_STANDARD 0xE000

#define FIPS_192 0xF192
#define FIPS_224 0xF224
#define FIPS_256 0xF256
#define FIPS_384 0xF384
#define FIPS_521 0xF521
#define FIPS_NOT_STANDARD 0xF000

#define UA_256_1 0x17561256 // UA 256/1
#define UA_256_2 0x27561256 // UA 256/2
#define UA_256_3 0x37561256 // UA 256/3
#define UA_256_4 0x47561256 // UA 256/4
#define UA_256_5 0x57561256 // UA 256/5

#define UA_384_1 0x17561384 // UA 384/1
#define UA_384_2 0x27561384 // UA 384/2
#define UA_384_3 0x37561384 // UA 384/3
#define UA_384_4 0x47561384 // UA 384/4
#define UA_384_5 0x57561384 // UA 384/5

#define UA_512_1 0x17561512 // UA 512/1
#define UA_512_2 0x27561512 // UA 512/2
#define UA_512_3 0x37561512 // UA 512/3
#define UA_512_4 0x47561512 // UA 512/4
#define UA_512_5 0x57561512 // UA 512/5

#define KUPYNA_HASH 0x01

#define P256_HASH_SIZE 32
#define P384_HASH_SIZE 64
#define P512_HASH_SIZE 64

#define HASH_ID_SIZE 8

#define P256_MAX_MSG_SIZE (256-P256_HASH_SIZE-HASH_ID_SIZE)
#define P384_MAX_MSG_SIZE (384-P384_HASH_SIZE-HASH_ID_SIZE)
#define P512_MAX_MSG_SIZE (512-P512_HASH_SIZE-HASH_ID_SIZE)

#define WINDOW_SIZE 4

typedef struct _Ec {
    GFElement d, a, b;
    BOOL isEdwards;

    u64 bitLen;
    u64 wordLen;

    BigInt n;
    u64 cofactor;
    EcPoint BasePoint;

    BigInt p;
    u64 curve_id;
    
    TGFMulFunc* GFMul;
    TGFSqrFunc* GFSqr;

    TEcAdd* EcAdd;
    TEcDouble* EcDouble;

    PRNG prng;

    EcPointProj* T; // for precomputations

    u32 hash_id;

    unsigned max_msg_size;
    unsigned hash_out_size;
} Ec;

typedef Ec EcEd;
typedef Ec EcW;


int  EcEdInit(EcEd* ecc, u64 bitLen, const BigInt p, const EcPoint* bp, const BigInt n, const GFElement d);
int  EcEdTwistedUAInit(EcEd* ecc, u64 bitLen, const BigInt p, const EcPoint* bp, const BigInt n, const GFElement d, u32 curve_id);
int  EcWInit(EcW* ecc, u64 bitLen, const BigInt p, const EcPoint* bp, const BigInt n, const GFElement a, const GFElement b);

int  EcInitStandardCurve(Ec* ecc, u64 bitLen, BOOL isEdwards);
int  EcInitStandardCurveById(Ec* ecc, u32 id);

void EcGenerateBasePoint(Ec* ecc, EcPoint* bp);

int  EcCheckPointInMainSubGroup(Ec* ecc, const EcPoint* P);
int  EcCheckPointOnCurve(Ec* ecc, const EcPoint* P);

void EcDump(Ec* ecc, char* buf);

void EcDestroy(Ec* ecc);

/* Arithmetic on Elliptic Curves */
int  EcPointCmp(Ec* ecc, const EcPoint* A, const EcPoint* B);
void EcCopy(Ec* ecc, EcPoint* dest, const EcPoint* src);
void EcCopyProj(Ec* ecc, EcPointProj* dest, const EcPointProj* src);

void EcConvertAffineToProjective(Ec* ecc, const EcPoint* P, EcPointProj* Q);
void EcConvertProjectiveToAffine(Ec* ecc, const EcPointProj* P, EcPoint* Q);

void EcAdd(Ec* ecc, const EcPoint* A, const EcPoint* B, EcPoint* C);
void EcDouble(Ec* ecc, const EcPoint* A, EcPoint* B);

void EcAddProj(Ec* ecc, const EcPointProj* A, const EcPointProj* B, EcPointProj* C);
void EcDoubleProj(Ec* ecc, const EcPointProj* A, EcPointProj* B);


/* 
Scalar Multiplications(all in the projective coordinates):
    - AddAndDouble naive:
    Result: Q = kP
    Q <- 0P
    H <- 1P
    for i from 0 to m do
        if k[i] = 1
            Q <- Q + H
        H <- 2H

    - Montgomery constant-time:
    Result: Q = kP
    Q <- 0P
    H <- 1P
    for i from 0 to m do
        if k[i] = 0
            H <- H + Q
            Q <- 2Q
        else 
            Q <- H + Q
            H <- 2H 
    - Windowed constant time for Edwards:
    Precomputation: 0P, 1P, ... (2^w - 1)P
    Result: Q = kP
    Q <- 0P
    for i from m/w downto 0 do
        Q <- (2^w)Q
        Q <- Q + k[i]P

    - wNAF(not tested) 
*/

void EcScalarMulWindowedPrecomputation(Ec* ecc, const EcPoint* A, EcPointProj** T, int windowSize);
void EcScalarMulWindowed(Ec* ecc, const EcPointProj* T, int windowSize, const BigInt k, EcPointProj* B);


void EcScalarMulwNAFPrecomputation(Ec* ecc, const EcPoint* A, EcPointProj** T, int windowSize);
void EcScalarMulwNAF(Ec* ecc, const EcPointProj* T, int windowSize, const BigInt k, EcPointProj* B);

void EcScalarMulNaive(Ec* ecc, const EcPointProj* A, const BigInt k, EcPointProj* B);

void EcScalarMulMontgomery(Ec* ecc, const EcPointProj* A, const BigInt k, EcPointProj* B);

/* Fast windowed ScalarMul */
void EcScalarMulByBasePoint(Ec* ecc, const BigInt k, EcPoint* B);

/* Not so fast, but constant time Montgomery ladder */
void EcScalarMul(Ec* ecc, const EcPoint * A, const BigInt k, EcPoint * B);

#ifdef __cplusplus
}
#endif

#endif /* EC_H */