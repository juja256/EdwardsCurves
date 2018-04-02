#include "ec.h"
#include "gf.h"
#include "dsa.h"

#include <stdio.h>
#include <string.h>
#ifdef _WIN64
#include <intrin.h>
#pragma intrinsic(_umul128) 
#include <windows.h>
#else
#include <x86intrin.h>
#include <time.h>
double GetTickCount(void) 
{
  struct timespec now;
  if (clock_gettime(CLOCK_MONOTONIC, &now))
    return 0;
  return now.tv_sec * 1000.0 + now.tv_nsec / 1000000.0;
}
#endif // _WIN64



void test_fips(u64 bit_len, int isEdwards) {
    Ec cur;
    EcPoint G, H, A, Z;
    EcPointProj B;
    BigInt n, d, p;
    GFElement e3, e4, X;

    int r = EcInitStandardCurve(&cur, bit_len, isEdwards);

    printf("Init status: %d %d %d\n", r, cur.bitLen, cur.wordLen);

    int isOnCurve = EcCheckPointInMainSubGroup(&cur, &(cur.BasePoint));
    printf("Base point, in main subgroup: %d\n", isOnCurve);

    double s1, e1, s2, e2;

    EcScalarMul(&cur, &(cur.BasePoint), cur.n, &Z);
    s1 = GetTickCount();
    r = EcScalarMul(&cur, &(cur.BasePoint), cur.n, &Z);
    e1 = GetTickCount();
    GFDump(&cur, Z.x);
    GFDump(&cur, Z.y);
    printf("Scalar Mul(Aff.), status: %d, time: %lf\n", r, e1-s1);
    
    EcConvertAffineToProjective(&cur, &(cur.BasePoint), &B);
    s2 = GetTickCount();
    r = EcScalarMulProj(&cur, &B, cur.n, &B);
    e2 = GetTickCount();
    EcConvertProjectiveToAffine(&cur, &B, &G);
    GFDump(&cur, G.x);
    GFDump(&cur, G.y);

    printf("Scalar Mul(Proj.), status: %d, time: %lf\n",r, e2-s2);

    EcGenerateBasePoint(&cur, &H);
    isOnCurve = EcCheckPointInMainSubGroup(&cur, &H);
    printf("Generated point, in main subgroup: %d\n", isOnCurve);
    GFDump(&cur, H.x);
    GFDump(&cur, H.y);
}

void test_eddsa(u64 bit_len, int isEdwards) {
    Ec cur;
    char name[60];
    
    EcInitStandardCurve(&cur, bit_len, isEdwards);
    EcDump(&cur, name);
    printf("%s\n", name);

    EcSignature sig;
    BigInt hash, key;
    EcPoint Q;
    memset(hash, 0, 7*8);
    strcpy((char*)hash, "Gopher" );
    EcDsaGenerateKey(&cur, key, &Q);
    EcDsaSign(&cur, key, hash, &sig);
    GFDump(&cur, sig.r);
    GFDump(&cur, sig.s);
    int s = EcDsaVerify(&cur, &Q, hash, &sig);
    printf("Signature Verification status: %d\n", s);
}

int main() {

	/*printf("------- Testing P-192 (Edwards) -------\n");
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
	test_fips(384, 0);*/

    printf("------- Testing ECDSA -------\n");
    test_eddsa(192, 1);
    test_eddsa(192, 0);
    //test_eddsa(224, 1);
    //test_eddsa(224, 0);
    test_eddsa(256, 1);
    test_eddsa(256, 0);
    test_eddsa(384, 1);
    test_eddsa(384, 0);
    #ifdef _WIN64
    system("pause");
    #endif
    return 0;
}
