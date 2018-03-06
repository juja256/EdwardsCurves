#include "eddsa.h"
#include <stdlib.h>
#include <time.h>


int EdDsaGenerateKey(const EcEd* ecc, BigInt key, EcPoint* Q) {
    srand(time(NULL));
    for (u64 i=0; i<ecc->wordLen; i++) {
        key[i] = rand() | ((u64)(rand()) << 32);
    }
    if (ecc->bitLen % 64 == 32) key[ecc->wordLen-1] &= 0xFFFFFFFF;
    if (GFCmp(ecc, key, ecc->n) == 1) {
        GFSub(ecc, key, ecc->n, key);
    }
    EcEdScalarMul(ecc, &(ecc->BasePoint), key, Q);
    return 0;
}

