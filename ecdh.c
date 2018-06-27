#include "ecdh.h"
#include "ec.h"

#define DH_SUCCESS 0
#define DH_FAIL -1


int EcDhStartKeyNegotiation(Ec* ecc, EcPoint* Q_B, BigInt d_A, EcPoint* P_enc) {
    if (!EcCheckPointInMainSubGroup(ecc, Q_B)) return DH_FAIL;
    EcScalarMul(ecc, Q_B, d_A, P_enc);
    return DH_SUCCESS;

}

int EcDhEndKeyNegotiation(Ec* ecc, EcPoint* P_enc, BigInt d_A, EcPoint* P_secret) {
    if (!EcCheckPointInMainSubGroup(ecc, P_enc)) return DH_FAIL;
    EcScalarMul(ecc, P_enc, d_A, P_secret);
    return DH_SUCCESS;
}