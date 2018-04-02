#include "ecdh.h"
#include "ec.h"

#define DH_SUCCESS 0
#define DH_FAIL -1


int EcDhStartKeyNegotiation(Ec* ecc, EcPoint* Q_B, BigInt d_A, EcPoint* P_enc) {
    if (!EcCheckPointInMainSubGroup(ecc, Q_B)) return DH_FAIL;
    EcPointProj P_p;
    EcConvertAffineToProjective(ecc, Q_B, &P_p);
    EcScalarMulProj(ecc, &P_p, d_A, &P_p);
    EcConvertProjectiveToAffine(ecc, &P_p, P_enc);
    return DH_SUCCESS;

}

int EcDhEndKeyNegotiation(Ec* ecc, BigInt d_A, EcPoint* P_enc, EcPoint* P_secret) {
    if (!EcCheckPointInMainSubGroup(ecc, P_enc)) return DH_FAIL;
    EcPointProj P_p;
    EcConvertAffineToProjective(ecc, P_enc, &P_p);
    EcScalarMulProj(ecc, &P_p, d_A, &P_p);
    EcConvertProjectiveToAffine(ecc, &P_p, P_secret);
    return DH_SUCCESS;
}