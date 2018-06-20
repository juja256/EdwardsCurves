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

#include "sha3.h"
#include <assert.h>


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
    
    EcConvertAffineToProjective(&cur, &(cur.BasePoint), &B);
    s2 = GetTickCount();
    for (int i=0;i<10;i++)
        EcScalarMulProj(&cur, &B, cur.n, &A);
    e2 = GetTickCount();

    EcScalarMulProj(&cur, &B, cur.n, &B);
    
    EcConvertProjectiveToAffine(&cur, &B, &G);
    GFDump(&cur, G.x);
    GFDump(&cur, G.y);

    printf("Scalar Mul, time: %lf\n", (e2-s2)/10);

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

static const char* hexConvTab = "0123456789abcdef";

void convertToHex(const BYTE* inBuf, int size, char* outBuf) {
    for (int i=0; i<size; i++) {
        outBuf[2*i+1] = hexConvTab[inBuf[i] & 0x0F];
        outBuf[2*i] = hexConvTab[(inBuf[i] & 0xF0) >> 4];
    }
    outBuf[2*size] = '\0';
}

BYTE ord(char a) {
    for (int i=0; i<16; i++) 
        if (a == hexConvTab[i]) return i;
    return 0;
}

void convertFromHex(const char* inBuf, BYTE* outBuf) {
    int arrLen = strlen(inBuf) / 2;
    for (int i=0; i<arrLen; i++) 
        outBuf[i] = (ord(inBuf[i]) << 4) | ord(inBuf[i+1]);
}

void test_sha3() {
    const char* empty = "";
    const char* empty224 = "6b4e03423667dbb73b6e15454f0eb1abd4597f9a1b078e3f5b5a6bc7";
    const char* empty256 = "a7ffc6f8bf1ed76651c14756a061d662f580ff4de43b49fa82d80a4b80f8434a";
    const char* empty384 = "0c63a75b845e4f7d01107d852e4c2485c51a50aaaa94fc61995e71bbee983a2ac3713831264adb47fb6bd1e058d5f004";
    const char* empty512 = "a69f73cca23a9ac5c8b567dc185a756e97c982164fe25859e0d1dcc1475c80a615b2123af1f5f94c11e3e9402c3ac558f500199d95b6d3e301758586281dcd26";
    const char* tqbfjotld = "The quick brown fox jumps over the lazy dog";
    const char* tqbfjotld224 = "d15dadceaa4d5d7bb3b48f446421d542e08ad8887305e28d58335795";
    const char* tqbfjotld256 = "69070dda01975c8c120c3aada1b282394e7f032fa9cf32f4cb2259a0897dfc04";
    const char* tqbfjotld384 = "7063465e08a93bce31cd89d2e3ca8f602498696e253592ed26f07bf7e703cf328581e1471a7ba7ab119b1a9ebdf8be41";
    const char* tqbfjotld512 = "01dedd5de4ef14642445ba5f5b97c15e47b9ad931326e4b0727cd94cefc44fff23f07bf543139939b49128caf436dc1bdee54fcb24023a08d9403f9b4bf0d450";

    Sha3Engine sha3;
    BYTE digest[512/8];
    char hexDigest[1024/8 + 1];

    /* empty tests */
    SHA3Sum(224, (const BYTE*)empty, 0, digest);
    convertToHex(digest, 224/8, hexDigest);
    printf("sha3-224(''): %s\n", hexDigest);
    assert( strcmp(empty224, hexDigest) == 0 );

    SHA3Sum(256, (const BYTE*)empty, 0, digest);
    convertToHex(digest, 256/8, hexDigest);
    printf("sha3-256(''): %s\n", hexDigest);
    assert( strcmp(empty256, hexDigest) == 0 );
    
    SHA3Sum(384, (const BYTE*)empty, 0, digest);
    convertToHex(digest, 384/8, hexDigest);
    printf("sha3-384(''): %s\n", hexDigest);
    assert( strcmp(empty384, hexDigest) == 0 );

    SHA3Sum(512, (const BYTE*)empty, 0, digest);
    convertToHex(digest, 512/8, hexDigest);
    printf("sha3-512(''): %s\n", hexDigest);
    assert( strcmp(empty512, hexDigest) == 0 );

    /* tqbfjotld tests */
    SHA3Sum(224, (const BYTE*)tqbfjotld, strlen(tqbfjotld), digest);
    convertToHex(digest, 224/8, hexDigest);
    printf("sha3-224('%s'): %s\n", tqbfjotld, hexDigest);
    assert( strcmp(tqbfjotld224, hexDigest) == 0 );

    SHA3Sum(256, (const BYTE*)tqbfjotld, strlen(tqbfjotld), digest);
    convertToHex(digest, 256/8, hexDigest);
    printf("sha3-256('%s'): %s\n", tqbfjotld, hexDigest);
    assert( strcmp(tqbfjotld256, hexDigest) == 0 );
    
    SHA3Sum(384, (const BYTE*)tqbfjotld, strlen(tqbfjotld), digest);
    convertToHex(digest, 384/8, hexDigest);
    printf("sha3-384('%s'): %s\n",tqbfjotld, hexDigest);
    assert( strcmp(tqbfjotld384, hexDigest) == 0 );

    SHA3Sum(512, (const BYTE*)tqbfjotld, strlen(tqbfjotld), digest);
    convertToHex(digest, 512/8, hexDigest);
    printf("sha3-512('%s'): %s\n",tqbfjotld, hexDigest);
    assert( strcmp(tqbfjotld512, hexDigest) == 0 );

    #define BUF_SIZE (1024*1024)
    #define CYCLES 100
    BYTE largeBuf[BUF_SIZE];
    double s1, e1;
    SHA3Init(&sha3, 224);
    s1 = GetTickCount();
    for (int i=0; i<CYCLES; i++)
        SHA3Update(&sha3, largeBuf, BUF_SIZE);
    e1 = GetTickCount();
    printf("sha3-%d rate: %lf Mb/sec\n", 224, 1000* (CYCLES) / (e1-s1));


    SHA3Init(&sha3, 256);
    s1 = GetTickCount();
    for (int i=0; i<CYCLES; i++)
        SHA3Update(&sha3, largeBuf, BUF_SIZE);
    e1 = GetTickCount();
    printf("sha3-%d rate: %lf Mb/sec\n", 256, 1000* (CYCLES) / (e1-s1));

    SHA3Init(&sha3, 384);
    s1 = GetTickCount();
    for (int i=0; i<CYCLES; i++)
        SHA3Update(&sha3, largeBuf, BUF_SIZE);
    e1 = GetTickCount();
    printf("sha3-%d rate: %lf Mb/sec\n", 384, 1000* (CYCLES) / (e1-s1));

    SHA3Init(&sha3, 512);
    s1 = GetTickCount();
    for (int i=0; i<CYCLES; i++)
        SHA3Update(&sha3, largeBuf, BUF_SIZE);
    e1 = GetTickCount();
    printf("sha3-%d rate: %lf Mb/sec\n", 512, 1000* (CYCLES) / (e1-s1));

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

    test_sha3();
    #ifdef _WIN64
    system("pause");
    #endif
    return 0;
}
