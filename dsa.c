#include "dsa.h"
#include "ec.h"
#include "gf.h"
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <stdio.h>

int EcDsaGenerateKey(Ec* ecc, BigInt key, EcPoint* Q) {
    PRNGGenerateSequence( &(ecc->prng), ecc->bitLen, (unsigned char*)key );
    while (GFCmp(ecc, key, ecc->n) == 1) {
        GFSub(ecc, key, ecc->n, key);
    }
    EcScalarMul(ecc, &(ecc->BasePoint), key, Q);
    return 0;
}

int EcDsaSign(Ec* ecc, const BigInt key, const BigInt hash, EcSignature* signature) {
    BigInt k;
    EcPoint P;
    
    gen_k:
    PRNGGenerateSequence( &(ecc->prng), ecc->bitLen, (unsigned char*)k ); // generate k

    while (GFCmp(ecc, k, ecc->n) != -1) {
        GFSub(ecc, k, ecc->n, k);
    }

    EcScalarMul(ecc, &(ecc->BasePoint), k, &P); // (x1, y1) = k*P

    if (GFCmp(ecc, P.x, zero) == 0) {
        goto gen_k; 
    }

    copy(signature->r, P.x, ecc->wordLen);
    
    while (GFCmp(ecc, signature->r, ecc->n) != -1) {
        GFSub(ecc, signature->r, ecc->n, signature->r);
    } // r = x1 mod n
    
    /* s = k_inv*(hash + key*r) mod n */
    mul_mod( ecc->wordLen, key, signature->r, ecc->n, signature->s );
    add_mod( ecc->wordLen, hash, signature->s, ecc->n, signature->s);
    inv_mod( ecc->wordLen, k, ecc->n, k); 
    mul_mod( ecc->wordLen, signature->s, k, ecc->n, signature->s);
    if (GFCmp(ecc, signature->s, zero) == 0) goto gen_k;
    return 0;
}

int EcDsaVerify(Ec* ecc, const EcPoint* Q, const BigInt hash, const EcSignature* signature) {
    //if ( !EcCheckPointInMainSubGroup(ecc, Q) ) return VER_BROKEN_KEY;
    if ( ( GFCmp(ecc, signature->r, unity) != 1 ) || ( GFCmp(ecc, signature->r, ecc->n) != -1 ) 
        || ( GFCmp(ecc, signature->s, unity) != 1 ) || ( GFCmp(ecc, signature->s, ecc->n) != -1 ) ) return VER_BROKEN_SIGNATURE;

    BigInt u1, u2, s, v;
    inv_mod(ecc->wordLen, signature->s, ecc->n, s);
    mul_mod(ecc->wordLen, s, hash, ecc->n, u1); // u1 = s_inv * hash mod n
    mul_mod(ecc->wordLen, s, signature->r, ecc->n, u2); // u2 = s_inv * r mod n

    EcPointProj P_p, Q_p;
    EcPoint P;
    EcConvertAffineToProjective(ecc, &(ecc->BasePoint), &P_p);
    EcConvertAffineToProjective(ecc, Q, &Q_p);
    EcScalarMulNaive(ecc, &Q_p, u2, &Q_p);
    EcScalarMulNaive(ecc, &P_p, u1, &P_p);
    EcAddProj(ecc, &P_p, &Q_p, &P_p); // P = u1*G + u2*Q
    EcConvertProjectiveToAffine(ecc, &P_p, &P);

    copy(v, P.x, ecc->wordLen);
    while (GFCmp(ecc, v, ecc->n) != -1) {
        GFSub(ecc, v, ecc->n, v);
    } // v = P.x mod n

    /* v =? r */
    if ( GFCmp( ecc, v, signature->r ) != 0 ) return VER_BROKEN_SIGNATURE;
    return VER_OK;
}
