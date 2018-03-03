#include "eced.h"
#include <stdio.h>


EcEd curveInit(int modLength)
{
	EcEd cur;
	EcPoint G;
	BigInt n, d, p;

	switch(modLength)
	{
	case 192:
		GFInitFromString(G.x, "44F083BB00E51AD91A2743284D31F57EE5C84826FCC91F4B");
		GFInitFromString(G.y, "15FC16E5870524E0DBBE9EC8BB9F066C02A02B1978D4E029");
		GFInitFromString(d,   "6DBA6A");
		GFInitFromString(p,   "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFF");
		GFInitFromString(n,   "3FFFFFFFFFFFFFFFFFFFFFFFEA75D4027230DD4DFFDB0455");
		break;
	case 224:
		GFInitFromString(G.x, "C448CA02660F57204FF1BDE2B5CC3E25606A7460399FEA3DA9A06383");
		GFInitFromString(G.y, "319117770D6FC7FE35F6A02905FE1F363156BD2E5B75BB89A64CAFAB");
		GFInitFromString(d,   "3608425");
		GFInitFromString(p,   "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000000000000000000001");
		GFInitFromString(n,   "400000000000000000000000000020BBEC47CEDB34DD05BCB6B7E619");
		break;
	case 256:
		GFInitFromString(G.x, "894F8283626AEE6848515DDDC3B8DBB3D5302DEE0EE75080D6753E4D39BA5AB2");
		GFInitFromString(G.y, "EA612346223F6480CBBAFA39DB95D54D21469DD3074A957EFDA4FD79FEB630B5");
		GFInitFromString(d,   "72A38");
		GFInitFromString(p,   "FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF");
		GFInitFromString(n,   "3FFFFFFFC00000003FFFFFFFFFFFFFFFBA76FA29C30CC3AA4954B53EDBE54D75");
		break;
	case 384:
		GFInitFromString(G.x, 
			"1FC0E8E61F599813E376D11F7510D77F177C2F1CDE19FD14D63A2EC5EAD4D0DED1BD06676CCF365243BF3C0675A31B62");
		GFInitFromString(G.y, 
			"F52B4FA352B257D7A102FA45C56A50CCBDB3DEC053D5610EDBD0188C11F321F28A43E2FC50395E4A8BD0029AE71D51AA");
		GFInitFromString(d,           
			"12593A");
		GFInitFromString(p,        
			"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFF0000000000000000FFFFFFFF");
		GFInitFromString(n,           
			"4000000000000000000000000000000000000000000000005063576B5A9A0C3A23E9510EA680650B4884E63A2968DD71");
		break;
	default:
		break;
	}

	int r = EcEdInit(&cur, &G, modLength, n, d);
	printf("Init status: %d %d %d\n", r, cur.bitLen, cur.wordLen);

	return cur;
}

// modLength = [192,224,256,382]
void testCurve(int modLength)
{
	printf("Testing P-%d\n\n",modLength);

	EcEd cur = curveInit(modLength);
	
	EcPoint H, A, B, Z;
	EcPointProj PP;
	GFElement e1, e2, e3, e4, X;

	GFInitFromString(X, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF0000000000000000");
	GFInitFromString(Z.x, "00");
	GFInitFromString(Z.y, "01");
	
	printf("BasePoint.x, BasePoint.y: \n");
	GFDump(&cur, cur.BasePoint.x);
	GFDump(&cur, cur.BasePoint.y);

	int isOnCurve = EcEdCheckPointOnCurve(&cur, &cur.BasePoint);
	printf("Base point, on curve: %d\n", isOnCurve);

	printf("BasePoint.x + BasePoint.y: \n");
	GFAdd(&cur, cur.BasePoint.x, cur.BasePoint.y, e1);
	GFDump(&cur, e1);
	
	printf("BasePoint.y - BasePoint.x: \n");
	GFSub(&cur, cur.BasePoint.y, cur.BasePoint.x, e3);
	GFDump(&cur, e3);
	
	printf("BasePoint.x * BasePoint.y: \n");
	cur.GFMul(&cur, cur.BasePoint.x, cur.BasePoint.y, e4);
	GFDump(&cur, e4);

	printf("p(mod) + p(mod): \n");
	GFAdd(&cur, cur.p, cur.p, e2);
	GFDump(&cur, e2);

	printf("p(mod) - p(mod): \n");
	GFSub(&cur, cur.p, cur.p, e2);
	GFDump(&cur, e2);

	printf("BasePoint + (0,1) -> B.x,B.y: \n");
	EcEdAdd(&cur, &cur.BasePoint, &Z, &B);
	GFDump(&cur, B.x);
	GFDump(&cur, B.y);

	printf("Inv(BasePoint.x): \n");
	GFInv(&cur, cur.BasePoint.x, e4);
	GFDump(&cur, e4);
	
	printf("Sqr(BasePoint): \n");
	cur.GFSqr(&cur, X, e4);
	GFDump(&cur, e4);

	e1[0] = 5; e1[1] = 0; e1[2] = 0; e1[3] = 0;
	printf("Pow(BasePoint.x,e1): \n");
	GFPow(&cur, cur.BasePoint.x, e1, e4);
	GFDump(&cur, e4);

	u64 s, e;

	printf("ScalarMul(BasePoint,n) -> B.X,B.Y,B.Z: \n");
	EcEdScalarMul(&cur, &cur.BasePoint, cur.n, &B);
	s = __rdtsc();
	EcEdScalarMul(&cur, &cur.BasePoint, cur.n, &B);
	e = __rdtsc();
	GFDump(&cur, B.x);
	GFDump(&cur, B.y);
	printf("Scalar Mul(Hom.) time: %d\n", e - s);

	printf("ScalarMulOrdinary(BasePoint,n) -> B.x, B.y: \n");
	EcEdScalarMulOrdinary(&cur, &cur.BasePoint, cur.n, &B);
	s = __rdtsc();
	EcEdScalarMulOrdinary(&cur, &cur.BasePoint, cur.n, &B);
	e = __rdtsc();
	GFDump(&cur, B.x);
	GFDump(&cur, B.y);
	printf("Scalar Mul(Ord.) time: %d\n", e - s);

	printf("Gegerating point.......\n");
	EcEdGenerateBasePoint(&cur, &B);
	isOnCurve = EcEdCheckPointOnCurve(&cur, &B);
	printf("Generated point, on curve: %d\n", isOnCurve);
	GFDump(&cur, B.x);
	GFDump(&cur, B.y);
}


int main() {
	
	testCurve(224);
	
#ifdef _WIN64
    system("pause");
#endif
    return 0;
}