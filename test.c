#include "eced.h"
#include <stdio.h>

/*
    Curves for testing:

        P-192:
    p = FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFF
    d = 6DBA6A
    n = 3FFFFFFFFFFFFFFFFFFFFFFFEA75D4027230DD4DFFDB0455
    x = 44F083BB00E51AD91A2743284D31F57EE5C84826FCC91F4B
    y = 15FC16E5870524E0DBBE9EC8BB9F066C02A02B1978D4E029

        P-224
    p = FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000000000000000000001
    d = 3608425
    n = 400000000000000000000000000020BBEC47CEDB34DD05BCB6B7E619
    x = C448CA02660F57204FF1BDE2B5CC3E25606A7460399FEA3DA9A06383
    y = 319117770D6FC7FE35F6A02905FE1F363156BD2E5B75BB89A64CAFAB

        P-256
    p = FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF
    d = 72A38
    n = 3FFFFFFFC00000003FFFFFFFFFFFFFFFBA76FA29C30CC3AA4954B53EDBE54D75
    x = 894F8283626AEE6848515DDDC3B8DBB3D5302DEE0EE75080D6753E4D39BA5AB2
    y = EA612346223F6480CBBAFA39DB95D54D21469DD3074A957EFDA4FD79FEB630B5

        P-384
    p = FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFF0000000000000000FFFFFFFF
    d = 12593A
    n = 4000000000000000000000000000000000000000000000005063576B5A9A0C3A23E9510EA680650B4884E63A2968DD71
    x = 1FC0E8E61F599813E376D11F7510D77F177C2F1CDE19FD14D63A2EC5EAD4D0DED1BD06676CCF365243BF3C0675A31B62
    y = F52B4FA352B257D7A102FA45C56A50CCBDB3DEC053D5610EDBD0188C11F321F28A43E2FC50395E4A8BD0029AE71D51AA
*/

void test_fips192() {
    EcEd cur;
    EcPoint G, H, A, B, Z;
    BigInt n, d, p;
    GFElement e1, e2 ,e3, e4, X;
    GFInitFromString(G.x, "44F083BB00E51AD91A2743284D31F57EE5C84826FCC91F4B");
    GFInitFromString(G.y, "15FC16E5870524E0DBBE9EC8BB9F066C02A02B1978D4E029");
    GFInitFromString(X,   "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF0000000000000000");
    GFInitFromString(d,   "6DBA6A");
    GFInitFromString(p,   "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFF");
    GFInitFromString(n,   "3FFFFFFFFFFFFFFFFFFFFFFFEA75D4027230DD4DFFDB0455");
    GFInitFromString(Z.x, "00");
    GFInitFromString(Z.y, "01");
    int r = EcEdInit(&cur, &G, 192, n, d);
    printf("Init status: %d %d %d\n", r, cur.bitLen, cur.wordLen);
    
    GFDump(&cur, G.x);
    GFDump(&cur, G.y);

    GFAdd(&cur, p, p, e1);
    int isOnCurve = EcEdCheckPointOnCurve(&cur, &G);
    printf("Base point, on curve: %d\n", isOnCurve);

    GFAdd(&cur, G.x, G.y, e1);
    GFDump(&cur, e1);
    GFSub(&cur, p, p, e2);
    GFDump(&cur, e2);
    GFSub(&cur, G.y, G.x, e3);
    GFDump(&cur, e3);
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
    isOnCurve = EcEdCheckPointOnCurve(&cur, &B);
    printf("Generated point, on curve: %d\n", isOnCurve);
    GFDump(&cur, B.x);
    GFDump(&cur, B.y);
}

void test_fips224() {
    EcEd cur;
    EcPoint G, H, A, B, Z;
    BigInt n, d, p;
    GFElement e1, e2 ,e3, e4, X;
	
	
	GFInitFromString(G.x, "C448CA02660F57204FF1BDE2B5CC3E25606A7460399FEA3DA9A06383");
    GFInitFromString(G.y, "319117770D6FC7FE35F6A02905FE1F363156BD2E5B75BB89A64CAFAB");
    GFInitFromString(X,   "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF0000000000000000");
    GFInitFromString(d,   "3608425");
    GFInitFromString(p,   "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000000000000000000001");
    GFInitFromString(n,   "400000000000000000000000000020BBEC47CEDB34DD05BCB6B7E619");
    GFInitFromString(Z.x, "00");
    GFInitFromString(Z.y, "01");
    int r = EcEdInit(&cur, &G, 224, n, d);
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
    printf("on curve %d\n",EcEdCheckPointOnCurve(&cur, &B));

    GFDump(&cur, B.x);
    GFDump(&cur, B.y);
    GFInv(&cur, G.x, e4);
    GFDump(&cur, e4);
    cur.GFSqr(&cur, X, e4);
    GFDump(&cur, e4);

    e1[0]=5; e1[1] = 0; e1[2] = 0; e1[3] = 0;
    GFPow(&cur, G.x, e1, e4);
    GFDump(&cur, e4);

    u64 s, e;
    EcEdScalarMul(&cur, &G, cur.n, &B);
    s = __rdtsc();
    EcEdScalarMul(&cur, &G, cur.n, &B);
    e = __rdtsc();
    GFDump(&cur, B.x);
    GFDump(&cur, B.y);
    printf("on curve %d\n",EcEdCheckPointOnCurve(&cur, &B));
    printf("Scalar Mul(Hom.) time: %d\n", e-s);
    
    EcEdScalarMulOrdinary(&cur, &G, cur.n, &B);
    s = __rdtsc();
    EcEdScalarMulOrdinary(&cur, &G, cur.n, &B);
    e = __rdtsc();
    GFDump(&cur, B.x);
    GFDump(&cur, B.y);
    printf("on curve %d\n",EcEdCheckPointOnCurve(&cur, &B));
    printf("Scalar Mul(Ord.) time: %d\n", e-s);

    EcEdGenerateBasePoint(&cur, &B);
    int isOnCurve = EcEdCheckPointOnCurve(&cur, &B);
    printf("Generated point, on curve: %d\n", isOnCurve);
    GFDump(&cur, B.x);
    GFDump(&cur, B.y);
}

void test_fips256() {

    EcEd cur;
    EcPoint G, H, A, B, Z;
    BigInt n, d, p;
    GFElement e1, e2 ,e3, e4, X;
    GFInitFromString(G.x, "894F8283626AEE6848515DDDC3B8DBB3D5302DEE0EE75080D6753E4D39BA5AB2");
    GFInitFromString(H.y, "EA612346223F6480CBBAFA39DB95D54D21469DD3074A957EFDA4FD79FEB630B5");
    GFInitFromString(G.y, "EA612346223F6480CBBAFA39DB95D54D21469DD3074A957EFDA4FD79FEB630B5");
    GFInitFromString(X,   "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF0000000000000000");
    GFInitFromString(d,   "72A38");
    GFInitFromString(p,   "FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF");
    GFInitFromString(n,   "3FFFFFFFC00000003FFFFFFFFFFFFFFFBA76FA29C30CC3AA4954B53EDBE54D75");
    GFInitFromString(Z.x, "00");
    GFInitFromString(Z.y, "01");
    int r = EcEdInit(&cur, &G, 256, n, d);
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
    printf("on curve %d\n",EcEdCheckPointOnCurve(&cur, &B));

    GFDump(&cur, B.x);
    GFDump(&cur, B.y);
    GFInv(&cur, G.x, e4);
    GFDump(&cur, e4);
    cur.GFSqr(&cur, X, e4);
    GFDump(&cur, e4);

    e1[0]=5; e1[1] = 0; e1[2] = 0; e1[3] = 0; e1[4] = 0;
    GFPow(&cur, G.x, e1, e4);
    GFDump(&cur, e4);

    GFSub(&cur, cur.p, G.x, H.x);
    GFDump(&cur, H.x);
    GFDump(&cur, H.y);
    EcEdAdd(&cur, &G, &H, &B);
    printf("on curve %d\n",EcEdCheckPointOnCurve(&cur, &B));
    GFDump(&cur, B.x);
    GFDump(&cur, B.y);

    u64 s, e;
    EcEdScalarMul(&cur, &G, cur.n, &B);
    s = __rdtsc();
    EcEdScalarMul(&cur, &G, cur.n, &B);
    e = __rdtsc();
    printf("on curve %d\n",EcEdCheckPointOnCurve(&cur, &B));
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

void test_fips384() {
	EcEd cur;
	EcPoint G, H, A, B, Z;
	BigInt n, d, p;
	GFElement e1, e2, e3, e4, X;
	GFInitFromString(G.x, "1FC0E8E61F599813E376D11F7510D77F177C2F1CDE19FD14D63A2EC5EAD4D0DED1BD06676CCF365243BF3C0675A31B62");
	GFInitFromString(G.y, "F52B4FA352B257D7A102FA45C56A50CCBDB3DEC053D5610EDBD0188C11F321F28A43E2FC50395E4A8BD0029AE71D51AA");
	GFInitFromString(X, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF0000000000000000");
	GFInitFromString(d, "12593A");
	GFInitFromString(p, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFF0000000000000000FFFFFFFF");
	GFInitFromString(n, "4000000000000000000000000000000000000000000000005063576B5A9A0C3A23E9510EA680650B4884E63A2968DD71");
	GFInitFromString(Z.x, "00");
	GFInitFromString(Z.y, "01");
	int r = EcEdInit(&cur, &G, 384, n, d);
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
    printf("on curve %d\n",EcEdCheckPointOnCurve(&cur, &B));

	GFDump(&cur, B.x);
	GFDump(&cur, B.y);
	GFInv(&cur, G.x, e4);
	GFDump(&cur, e4);
	cur.GFSqr(&cur, X, e4);
	GFDump(&cur, e4);

	e1[0] = 5; e1[1] = 0; e1[2] = 0; e1[3] = 0;
	GFPow(&cur, G.x, e1, e4);
	GFDump(&cur, e4);

	u64 s, e;
	EcEdScalarMul(&cur, &G, cur.n, &B);
	s = __rdtsc();
	EcEdScalarMul(&cur, &G, cur.n, &B);
	e = __rdtsc();
	GFDump(&cur, B.x);
	GFDump(&cur, B.y);
	printf("Scalar Mul(Hom.) time: %d\n", e - s);
    printf("on curve %d\n",EcEdCheckPointOnCurve(&cur, &B));

	EcEdScalarMulOrdinary(&cur, &G, cur.n, &B);
	s = __rdtsc();
	EcEdScalarMulOrdinary(&cur, &G, cur.n, &B);
	e = __rdtsc();
	GFDump(&cur, B.x);
	GFDump(&cur, B.y);
	printf("Scalar Mul(Ord.) time: %d\n", e - s);

	EcEdGenerateBasePoint(&cur, &B);
	int isOnCurve = EcEdCheckPointOnCurve(&cur, &B);
	printf("Generated point, on curve: %d\n", isOnCurve);
	GFDump(&cur, B.x);
	GFDump(&cur, B.y);
}

void test()
{
	EcEd cur;
	EcPoint G, H, A, B, Z;
	BigInt n, d, p;
	GFElement e1, e2, e3, e4, X;
	GFInitFromString(G.x, "1FC0E8E61F599813E376D11F7510D77F177C2F1CDE19FD14D63A2EC5EAD4D0DED1BD06676CCF365243BF3C0675A31B62");
	GFInitFromString(G.y, "F52B4FA352B257D7A102FA45C56A50CCBDB3DEC053D5610EDBD0188C11F321F28A43E2FC50395E4A8BD0029AE71D51AA");
	GFInitFromString(d,   "12593A");
	GFInitFromString(p,   "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFF0000000000000000FFFFFFFF");
	GFInitFromString(n,   "4000000000000000000000000000000000000000000000005063576B5A9A0C3A23E9510EA680650B4884E63A2968DD71");
	GFInitFromString(Z.x, "00");
	GFInitFromString(Z.y, "01");
	int r = EcEdInit(&cur, &G, 384, n, d);

	printf("a+b: ");
	GFAdd(&cur, G.x, G.y , e1);
	GFDump(&cur, e1);
	printf("\n");

	printf("a-b: ");
	GFSub(&cur, G.x, G.y, e2);
	GFDump(&cur, e2);
	printf("\n");


	printf("a*b: ");
	GFMul(&cur, G.x, G.y, e3);
	GFDump(&cur, e3);
	printf("\n");

	int isOnCurve = EcEdCheckPointOnCurve(&cur, &G);
	printf("Generated point, on curve: %d\n", isOnCurve);

}

void test2()
{
	EcEd cur;
	EcPoint G, H, A, B, Z;
	BigInt n, d, p;
	GFElement e1, e2, e3, e4, X;

	GFInitFromString(G.x, "C448CA02660F57204FF1BDE2B5CC3E25606A7460399FEA3DA9A06383");
	GFInitFromString(G.y, "419117770D6FC7FE35F6A02905FE1F363156BD2E5B75BB89A64CAFAB");
	GFInitFromString(X, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF0000000000000000");
	GFInitFromString(d, "3608425");
	GFInitFromString(p, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000000000000000000001");
	GFInitFromString(n, "400000000000000000000000000020BBEC47CEDB34DD05BCB6B7E619");
	GFInitFromString(Z.x, "00");
	GFInitFromString(Z.y, "01");
	int r = EcEdInit(&cur, &G, 224, n, d);

	printf("a+b: ");
	GFAdd(&cur, G.x, G.y, e1);
	GFDump(&cur, e1);
	printf("\n");

	printf("a-b: ");
	GFSub(&cur, G.x, G.y, e2);
	GFDump(&cur, e2);
	printf("\n");


	printf("a*b: ");
	GFMul(&cur, G.x, G.y, e3);
	GFDump(&cur, e3);
	printf("\n");

	int isOnCurve = EcEdCheckPointOnCurve(&cur, &G);
	printf("Generated point, on curve: %d\n", isOnCurve);
}

extern void basic_reduction(u64 n, const BigInt a, const BigInt p, BigInt res);
extern int bigint_bit_len(u64 nWords, const BigInt a);

void test_eddsa() {
    EcEd cur;
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

	printf("------- Testing P-192 -------\n");
    test_fips192();

    printf("------- Testing P-224 -------\n");
    test_fips224();

    printf("------- Testing P-256 -------\n");
    test_fips256();
	
	printf("------- Testing P-384 -------\n");
	test_fips384();

    printf("------- Testing ECDSA -------\n");
    test_eddsa();
    #ifdef _WIN64
    system("pause");
    #endif
    return 0;
}
