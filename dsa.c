#include "dsa.h"
#include "ec.h"
#include "gf.h"
#include <stdlib.h>
#include <time.h>
#include <string.h>

/* Magic 2048bit Number M=pq:
0x82780CC55BAE0479B0F478FE1F648D21180C71EF09655F1103F21B98765A4A926034967321666C59056461D85616746BAF309D393692C7EBCE285DCC5F865868E0C9F2048DFD87240229B039BAA07BE52FCBBFF95DF7E8BD2189A00B3D4832510D15BE80354A6560CED52426F7404EBCC3285392E894037FB02B01B19E47CDAA2DD029EC78A56963EDF8A6EA9E670F964E96C3ED6912DB8EBDCDC8E959F36ECF37CD464B9153D2FDD2EAA12E5982B0A33B448F290FA31868FE48BFF5339C9AC8697A2C040AF823922AD9EB9807E45912EBDDD15BB5AD3AB04C21C8774E41D42486ABFDC7EECC392FE0629372B63A4334AA23A2F721F7AF131201E43B9CE2C787 
*/
static const u64 M[] = { 
0x1201e43b9ce2c787,
0xaa23a2f721f7af13,
0xe0629372b63a4334,
0x86abfdc7eecc392f,
0x4c21c8774e41d424,
0xebddd15bb5ad3ab0,
0x2ad9eb9807e45912,
0x697a2c040af82392,
0xfe48bff5339c9ac8,
0x3b448f290fa31868,
0xd2eaa12e5982b0a3,
0x37cd464b9153d2fd,
0xbdcdc8e959f36ecf,
0x4e96c3ed6912db8e,
0xedf8a6ea9e670f96,
0x2dd029ec78a56963,
0xb02b01b19e47cdaa,
0xc3285392e894037f,
0xced52426f7404ebc,
0xd15be80354a6560,
0x2189a00b3d483251,
0x2fcbbff95df7e8bd,
0x229b039baa07be5,
0xe0c9f2048dfd8724,
0xce285dcc5f865868,
0xaf309d393692c7eb,
0x56461d85616746b,
0x6034967321666c59,
0x3f21b98765a4a92,
0x180c71ef09655f11,
0xb0f478fe1f648d21,
0x82780cc55bae0479
 };

int BBSInit(BBS2048_Generator* generator, BigInt seed) {
    copy( generator->state, seed, 2048/64 );
    mul_mod(32, generator->state, generator->state, M, generator->state);
}

int BBSGenerateSequence(BBS2048_Generator* generator, u64 bit_len, BigInt dest) {
    int w = (bit_len % 64 == 0) ? bit_len/64 : bit_len/64+1; 

    memset(dest, 0, w);
    for (int i=0; i<bit_len; i++) {
        mul_mod(32, generator->state, generator->state, M, generator->state);
        dest[i/64] ^= ( generator->state[0] & 1 ) << (i%64);
    }
}

int EcDsaGenerateKey(Ec* ecc, BigInt key, EcPoint* Q) {
    BigInt seed;
    EcPointProj Q_p;

    BBSGenerateSequence( &(ecc->prng), ecc->bitLen, key );

    if (GFCmp(ecc, key, ecc->n) == 1) {
        GFSub(ecc, key, ecc->n, key);
    }

    EcConvertAffineToProjective(ecc, &(ecc->BasePoint), &Q_p);
    EcScalarMulProj(ecc, &Q_p, key, &Q_p);
    EcConvertProjectiveToAffine(ecc, &Q_p, Q);
    return 0;
}

int EcDsaSign(Ec* ecc, const BigInt key, const BigInt hash, EcSignature* signature) {
    BigInt k;
    EcPoint P;
    EcPointProj Q_p;
    
    gen_k:
    BBSGenerateSequence( &ecc->prng, ecc->bitLen, k ); // generate k

    if (GFCmp(ecc, k, ecc->n) == 1) {
        GFSub(ecc, k, ecc->n, k);
    }

    EcConvertAffineToProjective(ecc, &(ecc->BasePoint), &Q_p);
    EcScalarMulProj(ecc, &Q_p, k, &Q_p); // (x1, y1) = k*P
    EcConvertProjectiveToAffine(ecc, &Q_p, &P);

    if (GFCmp(ecc, P.x, ecc->n) == 1) {
        GFSub(ecc, k, ecc->n, signature->r);
    }
    else if (GFCmp(ecc, P.x, zero) == 0) {
        goto gen_k; 
    }
    else {
        copy(signature->r, P.x, ecc->wordLen);
    } // r = x1 mod n

    /* s = k_inv*(hash + key*r) mod n */
    mul_mod( ecc->wordLen, key, signature->r, ecc->n, signature->s );
    add_mod( ecc->wordLen, hash, signature->s, ecc->n, signature->s);
    inv_mod( ecc->wordLen, k, ecc->n, k); 
    mul_mod( ecc->wordLen, signature->s, k, ecc->n, signature->s);
    if (GFCmp(ecc, signature->s, zero) == 0) goto gen_k;
    return 0;
}

int EcDsaVerify(const Ec* ecc, const EcPoint* Q, const BigInt hash, const EcSignature* signature) {
    if ( !EcCheckPointInMainSubGroup(ecc, Q) ) return VER_BROKEN_KEY;
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
    EcScalarMulProj(ecc, &Q_p, u2, &Q_p);
    EcScalarMulProj(ecc, &P_p, u1, &P_p);
    EcAddProj(ecc, &P_p, &Q_p, &P_p); // P = u1*G + u2*Q
    EcConvertProjectiveToAffine(ecc, &P_p, &P);

    if (GFCmp(ecc, P.x, ecc->n) == 1) {
        GFSub(ecc, P.x, ecc->n, v);
    }
    else {
        copy(v, P.x, ecc->wordLen);
    } // v = P.x mod n

    /* v =? r */
    if ( GFCmp( ecc, v, signature->r ) != 0 ) return VER_BROKEN_SIGNATURE;
    return VER_OK;
}
