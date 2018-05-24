#ifndef SHA3_H
#define SHA3_H

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

typedef KeccakSpoonge SHA3State;

void Keccak_theta( KeccakState* state);
void Keccak_rho(KeccakState* state);
void Keccak_pi(KeccakState* state);
void Keccak_chi(KeccakState* state);
void Keccak_iota(KeccakState* state, int rnd);

void Keccak_p( KeccakState* state );
void SpoongeInit( KeccakSpoonge* spoonge, int capacity, int rate, RoundFunc* rnd);
void SpoongeAbsorb( KeccakSpoonge* spoonge, BYTE* inBuf, unsigned size );
void SpoongeSqueeze( KeccakSpoonge* spoonge, BYTE* outBuf, unsigned size );

void SHA3Init( SHA3State* state, unsigned digestSize);
void SHA3Update( SHA3State* state, BYTE inBuf, unsigned size );
void SHA3Final( SHA3State* state, BYTE inBuf, unsigned size );
void SHA3GetDigest( SHA3State* state, BYTE* digest );
void SHA3Sum( unsigned digestSize, BYTE* inBuf, unsigned size, BYTE* digest );

#endif