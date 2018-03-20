#include "ec.h"
#include "gf.h"

#include <stdio.h>

#ifdef _WIN64
#include <intrin.h>
#pragma intrinsic(_umul128) 
#else
#include <x86intrin.h>
#endif // _WIN64

void test_fips(u64 bit_len, int isEdwards) {
    Ec cur;
    EcPoint G, H, A, B, Z;
    BigInt n, d, p;
    GFElement e1, e2 ,e3, e4, X;

    int r = EcInitStandardCurve(&cur, bit_len, isEdwards);

    printf("Init status: %d %d %d\n", r, cur.bitLen, cur.wordLen);

    int isOnCurve = EcCheckPointOnCurve(&cur, &(cur.BasePoint));
    printf("Base point, on curve: %d\n", isOnCurve);

    u64 s, e;
    EcScalarMul(&cur, &(cur.BasePoint), cur.n, &B);
    s = __rdtsc();
    r = EcScalarMul(&cur, &(cur.BasePoint), cur.n, &B);
    e = __rdtsc();
    GFDump(&cur, B.x);
    GFDump(&cur, B.y);
    printf("Scalar Mul(Aff.), status: %d, time: %d\n", r, e-s);
    
    EcScalarMulProj(&cur, &(cur.BasePoint), cur.n, &B);
    s = __rdtsc();
    EcScalarMulProj(&cur, &(cur.BasePoint), cur.n, &B);
    e = __rdtsc();
    GFDump(&cur, B.x);
    GFDump(&cur, B.y);
    printf("Scalar Mul(Proj.) time: %d\n", e-s);

    EcGenerateBasePoint(&cur, &B);
    isOnCurve = EcCheckPointOnCurve(&cur, &B);
    printf("Generated point, on curve: %d\n", isOnCurve);
    GFDump(&cur, B.x);
    GFDump(&cur, B.y);
}

void test_eddsa() {
    Ec cur;
    EcPoint G, H, A, B, Z;
    BigInt n, d, p;
    GFElement e1, e2 ,e3, e4, X, a;
    GFInitFromString(G.x, "44F083BB00E51AD91A2743284D31F57EE5C84826FCC91F4B");
    GFInitFromString(G.y, "15FC16E5870524E0DBBE9EC8BB9F066C02A02B1978D4E029");
    GFInitFromString(X,   "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF0000000000000000");
    GFInitFromString(d,   "6DBA6A");
    GFInitFromString(p,   "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFF");
    GFInitFromString(n,   "3FFFFFFFFFFFFFFFFFFFFFFFEA75D4027230DD4DFFDB0455");
    GFInitFromString(Z.x, "00");
    GFInitFromString(Z.y, "01");
    int r = EcEdInit(&cur, &G, 192, n, d);

    GFInitFromString(p, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFF");
    GFInitFromString(a, "FFFFFFFFFFFFFFFFF0FFFFFFFFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFF");
    //GFDump(&cur, a);
    printf( "BitLen: %d\n",  bigint_bit_len(3, d) );
    printf("Reduction result: \n");
    basic_reduction(3, a, p, a);
    GFDump(&cur, a);
}

int main() {

	printf("------- Testing P-192 (Edwards) -------\n");
    test_fips(192, 1);
    printf("------- Testing P-192 (Weierstrass) -------\n");
    test_fips(192, 0);

    printf("------- Testing P-224 (Edwards) -------\n");
    test_fips(224, 1);
    printf("------- Testing P-224 (Weierstrass) -------\n");
    test_fips(224, 0);

    printf("------- Testing P-256 (Edwards) -------\n");
    test_fips(256, 1);
    printf("------- Testing P-256 (Weierstrass) -------\n");
    test_fips(256, 0);

	printf("------- Testing P-384 (Edwards) -------\n");
	test_fips(384, 1);
    printf("------- Testing P-384 (Weierstrass) -------\n");
	test_fips(384, 0);

    //printf("------- Testing ECDSA -------\n");
    //test_eddsa();
    #ifdef _WIN64
    system("pause");
    #endif
    return 0;
}
