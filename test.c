#include "eced.h"
#include <stdio.h>


void test_fips192() {
    EcEd cur;
    EcPoint G, H, A, B, Z;
    BigInt n, d, p;
    GFElement e1, e2 ,e3, e4, X;
    GFInitFromString(G.x, "AE709D07B2D112CECD4A7AE103757F2C101D054ACB1A0F17");
    GFInitFromString(H.y, "8BDF8D5A994AAB1E1455889D28E29CA3EC68548D3F7CA32F");
    GFInitFromString(G.y, "8BDF8D5A994AAB1E1455889D28E29CA3EC68548D3F7CA32F");
    GFInitFromString(X,   "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF0000000000000000");
    GFInitFromString(d,   "28453E");
    GFInitFromString(p,   "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFF");
    GFInitFromString(n,   "4000000000000000000000001AEAD1229D137F564D7FF6D5");
    GFInitFromString(Z.x, "00");
    GFInitFromString(Z.y, "01");
    int r = EcEdInit(&cur, &G, 192, n, d);
    printf("Init status: %d %d %d\n", r, cur.bitLen, cur.wordLen);
    
    GFDump(&cur, G.x);
    GFDump(&cur, G.y);

    GFAdd(&cur, G.x, G.y, e1);
    GFDump(&cur, e1);
    GFSub(&cur, G.x, G.y, e2);
    GFDump(&cur, e2);
    GFSub(&cur, G.y, G.x, e3);
    GFDump(&cur, G.x);
    cur.GFMul(&cur, G.x, G.y, e4);
    GFDump(&cur, e4);

    EcEdAdd(&cur, &G, &Z, &B);

    GFDump(&cur, B.x);
    GFDump(&cur, B.y);
    GFInv(&cur, G.x, e4);
    GFDump(&cur, e4);
    cur.GFSqr(&cur, X, e4);
    GFDump(&cur, e4);

    e1[0]=5; e1[1] = 0; e1[2] = 0; e1[3] = 0;
    GFPow(&cur, G.x, e1, e4);
    GFDump(&cur, e4);

    GFSub(&cur, cur.p, G.x, H.x);
    GFDump(&cur, H.x);
    GFDump(&cur, H.y);
    EcEdAdd(&cur, &G, &H, &B);
    GFDump(&cur, B.x);
    GFDump(&cur, B.y);

    u64 s, e;
    EcEdScalarMul(&cur, &G, cur.n, &B);
    s = __rdtsc();
    EcEdScalarMul(&cur, &G, cur.n, &B);
    e = __rdtsc();
    GFDump(&cur, B.x);
    GFDump(&cur, B.y);
    printf("Scalar Mul(Hom.) time: %d\n", e-s);
    
    EcEdScalarMulOrdinary(&cur, &G, cur.n, &B);
    s = __rdtsc();
    EcEdScalarMulOrdinary(&cur, &G, cur.n, &B);
    e = __rdtsc();
    GFDump(&cur, B.x);
    GFDump(&cur, B.y);
    printf("Scalar Mul(Ord.) time: %d\n", e-s);

    EcEdGenerateBasePoint(&cur, &B);
    int isOnCurve = EcEdCheckPointOnCurve(&cur, &B);
    printf("Generated point, on curve: %d\n", isOnCurve);
    GFDump(&cur, B.x);
    GFDump(&cur, B.y);
}

int main() {
    test_fips192();
    #ifdef _WIN64
    system("pause");
    #endif
    return 0;
}
