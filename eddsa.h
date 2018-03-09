#ifndef EDDSA_H
#define EDDSA_H

#include "eced.h"

typedef struct {
    BigInt r;
    BigInt s;
} EcSignature;

int EdDsaGenerateKey(const EcEd* ecc, BigInt key, EcPoint* Q);
int EdDsaSign(const EcEd* ecc, const BigInt key, const BigInt hash, EcSignature* signature);
int EdDsaVerify(const EcEd* ecc, const EcPoint* Q, const BigInt hash, const EcSignature* signature);

#endif /* EDDSA_H */
