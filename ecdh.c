#include "ecdh.h"
#include "ec.h"
#include "gf.h"
#include "kupyna.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#define DH_SUCCESS 0
#define DH_FAIL -1

#define ENCRYPT_SUCCESS 0
#define INVALID_CURVE -1
#define MESSAGE_TOO_LONG -2
#define HASH_ERROR -3;

#define DECRYPT_SUCCESS 0
#define INVALID_R -4
#define INVALID_C -5
#define HASH_INTEGRITY_ERROR -6

static void InvertBytes(unsigned char* a, unsigned sz) {
    unsigned char tmp;
    for (int i=0; i<sz/2; i++) {
        tmp = a[i];
        a[i] = a[sz-1-i];
        a[sz-1-i] = tmp;
    }
}

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

int UaKemGeneratePrivateKey(Ec* ecc, BigInt d, EcPoint* Q) {
    #ifdef UAKEM_DEBUG
    copy(d, zero, ecc->wordLen);
    d[0] = 0x25;
    #else
    PRNGGenerateSequence(&ecc->prng, ecc->bitLen, (unsigned char*)d);
    while (GFCmp(ecc, d, ecc->n) == 1) {
        GFSub(ecc, d, ecc->n, d);
    }
    #endif
    EcScalarMulByBasePoint(ecc, d, Q);
    #ifdef UAKEM_DEBUG
    printf("Q: \n");
    GFDump(ecc, Q->x);
    GFDump(ecc, Q->y);
    #endif
    return 0;
}

int UaKemEncrypt(Ec* ecc, EcPoint* Q, unsigned char* msg, unsigned size, Ciphertext* C) {
    if (ecc->curve_id & 0x0FFFF000 != 0x07561000) return INVALID_CURVE;
    if (size*8 > ecc->max_msg_size) return MESSAGE_TOO_LONG;

    BigInt M;
    kupyna_t kupyna_engine;
    if ((ecc->hash_id != KUPYNA_HASH) || (KupynaInit(ecc->hash_out_size, &kupyna_engine) != 0)) return HASH_ERROR;
    
    copy(M, zero, ecc->wordLen); // prepare zero padded
    memcpy(M, msg, size); // copy message to lower part
    unsigned short ML = size*8l; // calculate length
    memcpy((uint8_t*)M + ecc->max_msg_size/8, &ML, 2); // set length
    //InvertBytes((uint8_t*)M, ecc->max_msg_size/8 + 2);
    KupynaHash(&kupyna_engine, (uint8_t*)M, ecc->max_msg_size+16, (uint8_t*)M + 2 + ecc->max_msg_size/8); // set hash
    *((uint8_t*)M + ecc->bitLen/8 - 1) = (uint8_t)ecc->hash_id; // set hash id
    #ifdef UAKEM_DEBUG
    printf("M(m_size=%d): \n", ecc->max_msg_size);
    GFDump(ecc, M);
    #endif
    BigInt epsilon;
    GFElement t;
    EcPoint R, T;
    #ifndef UAKEM_DEBUG
    PRNGGenerateSequence(&ecc->prng, ecc->bitLen-2, (unsigned char*)epsilon); // generate one time epsilon
    #else
    copy(epsilon, zero, ecc->wordLen);
    epsilon[0] = 7;
    #endif
    EcScalarMulByBasePoint(ecc, epsilon, &R); // R = εP
    #ifdef UAKEM_DEBUG
    printf("R: \n");
    GFDump(ecc, R.x);
    GFDump(ecc, R.y);
    #endif
    EcScalarMul(ecc, Q, epsilon, &T);
    #ifdef UAKEM_DEBUG
    printf("T: \n");
    GFDump(ecc, T.x);
    GFDump(ecc, T.y);
    #endif
    GFMul(ecc, T.x, M, t); // t = xT * M mod p
    copy(C->r, R.x, ecc->wordLen);
    copy(C->t, t, ecc->wordLen);
    #ifdef UAKEM_DEBUG
    printf("C: \n");
    GFDump(ecc, C->r);
    GFDump(ecc, C->t);
    #endif
    return ENCRYPT_SUCCESS;
}

int UaKemDecrypt(Ec* ecc, BigInt e, EcPoint* Q, Ciphertext* C, unsigned char* msg, unsigned short* mlen) {
    if (!GFCmp(ecc, zero, C->r) || !GFCmp(ecc, unity, C->r)) return INVALID_R;
    GFElement x,y;
    GFSqr(ecc, C->r, x);
    GFNeg(ecc, x, y);
    GFAdd(ecc, y, unity, y); // y = 1-x^2
    GFMulByD(ecc, x);
    GFNeg(ecc, x, x);
    GFAdd(ecc, x, ecc->a, x); // x = a - dx^2
    GFInv(ecc, x, x);
    GFMul(ecc, x, y, y);
    if (GFLegendreSymbol(ecc, y) == -1) return INVALID_R;

    EcPoint R, T;
    copy(R.x, C->r, ecc->wordLen);
    GFSqrt(ecc, y, R.y);
    EcScalarMul(ecc, &R, e, &T); // T = eR 
    GFInv(ecc, T.x, x);
    GFMul(ecc, C->t, x, x); // x = t*x_inv

    kupyna_t kupyna_engine;
    uint8_t hash[16];
    if ((((uint8_t*)x)[ecc->bitLen/8 - 1] != KUPYNA_HASH) || (ecc->hash_id != KUPYNA_HASH) || (KupynaInit(ecc->hash_out_size, &kupyna_engine) != 0)) return HASH_ERROR;
    KupynaHash(&kupyna_engine, (uint8_t*)x, ecc->max_msg_size+16, hash);
    *mlen = *(unsigned short*)((uint8_t*)x + ecc->max_msg_size/8);
    if (memcmp(hash, (uint8_t*)x + 2 + ecc->max_msg_size/8, ecc->hash_out_size/8) != 0) return HASH_INTEGRITY_ERROR;
    memcpy(msg, x, *mlen/8);
    return DECRYPT_SUCCESS;
}