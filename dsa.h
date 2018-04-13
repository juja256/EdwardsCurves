#ifndef DSA_H
#define DSA_H

#include "ec.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    BigInt r;
    BigInt s;
} EcSignature;

#define VER_OK 0
#define VER_BROKEN_KEY 1
#define VER_BROKEN_SIGNATURE 2

int EcDsaGenerateKey(Ec* ecc, BigInt key, EcPoint* Q);
int EcDsaSign(Ec* ecc, const BigInt key, const BigInt hash, EcSignature* signature);
int EcDsaVerify(Ec* ecc, const EcPoint* Q, const BigInt hash, const EcSignature* signature);


#ifdef __cplusplus
}
#endif

#endif /* DSA_H */
