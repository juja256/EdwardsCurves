#ifndef EDDSA_H
#define EDDSA_H

#include "eced.h"

typedef struct {
    BigInt r;
    BigInt s;
} EcSignature;

int EdDsaGenerateKey(EcEd* ecc, BigInt* key, EcPoint* Q);
int EdDsaSign(EcEd* ecc, BigInt key, BigInt hash, EcSignature* signature);
int EdDsaVerify(EcEd* ecc, EcPoint* Q, BigInt hash, EcSignature* signature);

#endif /* EDDSA_H */
