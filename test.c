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



void test_ariphmetic(u64 bit_len, int isEdwards) {
    Ec cur;
    EcPoint G, H, Z;
    EcPointProj B, A;
    BigInt n, d, p;
    GFElement e3, e4, X;
    char name[60];
    int r = EcInitStandardCurve(&cur, bit_len, isEdwards);
    EcDump(&cur, name);
    printf("--- %s ---\n", name);

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
    for (int i=0;i<10;i++)
        EcScalarMulProj(&cur, &B, cur.n, &A);
    e2 = GetTickCount();

    r = EcScalarMulProj(&cur, &B, cur.n, &B);
    
    EcConvertProjectiveToAffine(&cur, &B, &G);
    GFDump(&cur, G.x);
    GFDump(&cur, G.y);

    printf("Scalar Mul(Proj.), status: %d, time: %lf\n",r, (e2-s2)/10);

    EcGenerateBasePoint(&cur, &H);
    isOnCurve = EcCheckPointInMainSubGroup(&cur, &H);
    printf("Generated point, in main subgroup: %d\n", isOnCurve);
    GFDump(&cur, H.x);
    GFDump(&cur, H.y);
}

void test_eddsa(u64 bit_len, int isEdwards) {
    double s1, e1, s2, e2;
    Ec cur;
    char name[60];
    
    EcInitStandardCurve(&cur, bit_len, isEdwards);
    EcDump(&cur, name);
    printf("--- %s ---\n", name);

    EcSignature sig;
    BigInt hash, key;
    EcPoint Q;
    memset(hash, 0, 7*8);
    strcpy((char*)hash, "Gopher" );
    s1 = GetTickCount();
    EcDsaGenerateKey(&cur, key, &Q);
    e1 = GetTickCount();
    EcDsaSign(&cur, key, hash, &sig);
    s2 = GetTickCount();
    int s = EcDsaVerify(&cur, &Q, hash, &sig);
    e2 = GetTickCount();
    printf("Generated ECDSA keypair (key, Q), time: %lf\n", e1-s1);
    GFDump(&cur, key);
    GFDump(&cur, Q.x);
    GFDump(&cur, Q.y);
    printf("Digital signature (r, s) for hash: %s, time: %lf\n", (char*)hash, s2-e1);
    GFDump(&cur, sig.r);
    GFDump(&cur, sig.s);
    printf("Signature Verification status: %d, time: %lf\n", s, e2-s2);
}

int main() {
    
	printf("------- Testing EC Arithmetic -------\n");
    test_ariphmetic(192, 0);
    test_ariphmetic(192, 1);
    test_ariphmetic(224, 0);
    test_ariphmetic(224, 1);
    test_ariphmetic(256, 0);
    test_ariphmetic(256, 1);
    test_ariphmetic(384, 0);
    test_ariphmetic(384, 1);
    test_ariphmetic(521, 0);
    test_ariphmetic(521, 1);

    printf("------- Testing ECDSA -------\n");
    test_eddsa(192, 0);
    test_eddsa(192, 1);
    test_eddsa(224, 0);
    test_eddsa(224, 1);
    test_eddsa(256, 0);
    test_eddsa(256, 1);
    test_eddsa(384, 0);
    test_eddsa(384, 1);
    test_eddsa(521, 0);
    test_eddsa(521, 1);
    #ifdef _WIN64
    system("pause");
    #endif
    return 0;
}
