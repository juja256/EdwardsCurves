#ifndef EDDH_H
#define EDDH_H

#include "ec.h"

int EcDhStartKeyNegotiation(Ec* ecc, EcPoint* Q_B, BigInt d_A, EcPoint* P_enc);
int EcDhEndKeyNegotiation(Ec* ecc, BigInt d_A, EcPoint* P_enc, EcPoint* P_secret);

#endif /* EDDH_H */
