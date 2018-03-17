#ifndef DSA_H
#define DSA_H

#include "ec.h"

typedef struct {
    BigInt r;
    BigInt s;
} EcSignature;

int DsaGenerateKey(const Ec* ecc, BigInt key, EcPoint* Q);
int DsaSign(const Ec* ecc, const BigInt key, const BigInt hash, EcSignature* signature);
int DsaVerify(const Ec* ecc, const EcPoint* Q, const BigInt hash, const EcSignature* signature);

#endif /* DSA_H */
