#ifndef EDDH_H
#define EDDH_H

#include "ec.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    GFElement r;
    GFElement t;
} Ciphertext;


int EcDhStartKeyNegotiation(Ec* ecc, EcPoint* Q_B, BigInt d_A, EcPoint* P_enc);
int EcDhEndKeyNegotiation(Ec* ecc, EcPoint* P_enc, BigInt d_A, EcPoint* P_secret);

int UaKemGeneratePrivateKey(Ec* ecc, BigInt d, EcPoint* Q);
int UaKemEncrypt(Ec* ecc, EcPoint* Q, unsigned char* msg, unsigned size, Ciphertext* C);
int UaKemDecrypt(Ec* ecc, BigInt d, EcPoint* Q, Ciphertext* C, unsigned char* msg);

#ifdef __cplusplus
}
#endif
#endif /* EDDH_H */
