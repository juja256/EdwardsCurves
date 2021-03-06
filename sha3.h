#ifndef SHA3_H
#define SHA3_H


#ifdef __cplusplus
extern "C" {
#endif

typedef unsigned char BYTE;

typedef unsigned long long KeccakState[25];

typedef void RoundFunc(KeccakState*);

typedef struct {
    int c;
    int r;
    int delayed;
    KeccakState state;
    RoundFunc* keccak_p;
} KeccakSpoonge;

typedef KeccakSpoonge Sha3Engine;

void Keccak_p( KeccakState* state );
void SpoongeInit( KeccakSpoonge* spoonge, int capacity, int rate, RoundFunc* rnd);
void SpoongeAbsorb( KeccakSpoonge* spoonge, const BYTE* inBuf, unsigned size );
void SpoongeSqueeze( KeccakSpoonge* spoonge, BYTE* outBuf, unsigned size );

#define INVALID_HASH_LEN -1
#define INIT_SUCCESS 0;

int SHA3Init( Sha3Engine* state, unsigned digestSize);
void SHA3Update( Sha3Engine* state, const BYTE* inBuf, unsigned size );
void SHA3Final( Sha3Engine* state, const BYTE* inBuf, unsigned size );
void SHA3GetDigest( Sha3Engine* state, BYTE* digest );
int SHA3Sum( unsigned digestSize, const BYTE* inBuf, unsigned size, BYTE* digest );

#ifdef __cplusplus
}
#endif

#endif