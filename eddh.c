#include "eddh.h"

#define SUCCESS 0

int EdDhStartKeyNegotioation(EcEd* ecc, EcPoint* Q_B, BigInt d_A, EcPoint* P_enc) {
    EcEdScalarMul(ecc, Q_B, d_A, P_enc);
    return SUCCESS;
}

int EdDhObtainSecretPoint(EcEd* ecc, BigInt d_A, EcPoint* P_enc, EcPoint* P_secret) {
    EcEdScalarMul(ecc, P_enc, d_A, P_secret);
    return SUCCESS;
}
