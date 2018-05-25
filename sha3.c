#include "sha3.h"
#include <string.h>
#include <stdio.h>

#define w 64
#define WORD unsigned long long

static inline void xor(BYTE* acc, BYTE* other, unsigned size) {
    for (unsigned i=0; i<size; i++) {
        acc[i] ^= other[i];
    }
}

#define ROT_L( X, Y ) (( X << Y ) | ( X >> (64 - Y) ))

static void dump(KeccakState st) {
    BYTE* bb = (BYTE*)st;
    for (int i=0; i<200; i++)
        printf("%02X ", bb[i]);
    printf("\n");
}

#define A(x,y,z) (( (*state)[5*(y)+(x)] & (1 << (z))) >> (z)) 
#define Parity(x,z) (A(x, 0, z) ^ A(x, 1, z) ^ A(x, 2, z) ^ A(x, 3, z) ^ A(x, 4, z))
#define Axor(x,y,z,b) (*state)[5*(y)+(x)] ^= ((b) << (z))
#define Mod(x,n) ( (x < 0) ? (x+n) : (x%n) )
#define Lane(x, y) (*state)[5*(y)+(x)]


static BYTE rho_offsets[5][5] = 
{ 
    {0, 1, 62, 28, 27},
    {36, 44, 6, 55, 20},
    {3, 10, 43, 25, 39},
    {41, 45, 15, 21, 8},
    {18, 2, 61, 56, 14} 
};

static WORD RC[24] = {
    0x0000000000000001, 0x0000000000008082,
    0x800000000000808A, 0x8000000080008000,
    0x000000000000808B, 0x0000000080000001,
    0x8000000080008081, 0x8000000000008009,
    0x000000000000008A, 0x0000000000000088,
    0x0000000080008009, 0x000000008000000A,
    0x000000008000808B, 0x800000000000008B,
    0x8000000000008089, 0x8000000000008003,
    0x8000000000008002, 0x8000000000000080,
    0x000000000000800A, 0x800000008000000A,
    0x8000000080008081, 0x8000000000008080,
    0x0000000080000001, 0x8000000080008008,
};
// A[x,y,z] = S[w(5y + x) + z]

void Keccak_theta(KeccakState* state) {
    WORD C[5];
    WORD D[5];

    C[0] = Lane(0, 0) ^ Lane(0, 1) ^ Lane(0, 2) ^ Lane(0, 3) ^ Lane(0, 4);
    C[1] = Lane(1, 0) ^ Lane(1, 1) ^ Lane(1, 2) ^ Lane(1, 3) ^ Lane(1, 4);
    C[2] = Lane(2, 0) ^ Lane(2, 1) ^ Lane(2, 2) ^ Lane(2, 3) ^ Lane(2, 4);
    C[3] = Lane(3, 0) ^ Lane(3, 1) ^ Lane(3, 2) ^ Lane(3, 3) ^ Lane(3, 4);
    C[4] = Lane(4, 0) ^ Lane(4, 1) ^ Lane(4, 2) ^ Lane(4, 3) ^ Lane(4, 4);

    D[0] = C[4] ^ ROT_L( C[1], 1 );
    D[1] = C[0] ^ ROT_L( C[2], 1 );
    D[2] = C[1] ^ ROT_L( C[3], 1 );
    D[3] = C[2] ^ ROT_L( C[4], 1 );
    D[4] = C[3] ^ ROT_L( C[0], 1 );
    
    for (int i=0; i<5; i++) {
        for (int j=0; j<5; j++) {
            Lane(i, j) ^= D[i];
        }
    }
}

void Keccak_rho(KeccakState* state) {
    for (int x=0; x<5; x++) {
        for (int y=0; y<5; y++) {
            Lane(x,y) = ROT_L( Lane(x,y), rho_offsets[x][y]); 
        }
    }
}

void Keccak_pi(KeccakState* state) {
    static KeccakState state_tmp;
    for (int x=0; x<5; x++) {
        for (int y=0; y<5; y++) {
            state_tmp[x+5*y] = Lane( Mod( (x+3*y), 5), x );
        }
    }
    memcpy(*state, state_tmp, sizeof(state_tmp));
}

void Keccak_chi(KeccakState* state) {
    WORD a1, a2, a3, a4, a0;

    for (int y=0; y<5; y++) {
        a0 = Lane(0, y);
        a1 = Lane(1, y);
        a2 = Lane(2, y);
        a3 = Lane(3, y);
        a4 = Lane(4, y);

        Lane(0, y) ^= (~a1) & a2;
        Lane(1, y) ^= (~a2) & a3;
        Lane(2, y) ^= (~a3) & a4;
        Lane(3, y) ^= (~a4) & a0;
        Lane(4, y) ^= (~a0) & a1;
    }
}

void Keccak_iota(KeccakState* state, int rnd) {
    (*state)[0] ^= RC[rnd];
}

void Keccak_p( KeccakState* state ) {
    for (int i=0; i<24; i++) {
        
        Keccak_theta(state);
        //printf("\ntheta ");
        //dump(*state);
        Keccak_rho(state);
        Keccak_pi(state);
        Keccak_chi(state);
        Keccak_iota(state, i);
    }
}

void SpoongeInit( KeccakSpoonge* spoonge, int capacity, int rate, RoundFunc* rnd) {
    spoonge->c = capacity/8;
    spoonge->r = rate/8,
    spoonge->keccak_p = rnd;
    spoonge->delayed = 0;
    for (int i=0; i<25; i++) spoonge->state[i] = 0;
}

void SpoongeAbsorb( KeccakSpoonge* spoonge, BYTE* inBuf, unsigned size ) {
    
    int delayedSize = (size < spoonge->r - spoonge->delayed) ? size : spoonge->r - spoonge->delayed;
    BYTE* statePtr = (BYTE*)spoonge->state;
    
    xor(statePtr + spoonge->delayed, inBuf, delayedSize);
    size -= delayedSize;
    inBuf += delayedSize;

    spoonge->delayed += delayedSize;

    if (spoonge->delayed == spoonge->r) {
        spoonge->keccak_p(&spoonge->state);
        spoonge->delayed = 0;
    }
    if (size == 0) {
        return;
    }

    for (int i=0; i<size / spoonge->r; i++) {
        xor(statePtr, inBuf + i*spoonge->r, spoonge->r);
        spoonge->keccak_p(&spoonge->state);
    }

    inBuf += size;
    xor(statePtr, inBuf + size, size % spoonge->r);
    spoonge->delayed = size % spoonge->r;
}

void SpoongeSqueeze( KeccakSpoonge* spoonge, BYTE* outBuf, unsigned size ) {
    int firstSize = (size <= spoonge->r) ? size : spoonge->r;
    
    memcpy(outBuf, spoonge->state, firstSize);

    int restSize = ((size - firstSize) % spoonge->r == 0) ? (size - firstSize)/spoonge->r : (size - firstSize)/spoonge->r + 1;
    for (int i=0; i<restSize; i++) {
        spoonge->keccak_p(&spoonge->state);
        int toBeCopied = (i == restSize-1) ? size % spoonge->r : spoonge->r; 
        memcpy(outBuf + firstSize + i*spoonge->r, spoonge->state, toBeCopied );
    }    
}

int SHA3Init( Sha3Engine* state, unsigned digestSize) {
    if (! ((digestSize == 224) || (digestSize == 256) || (digestSize == 384) || (digestSize == 512)) ) return INVALID_HASH_LEN;
    SpoongeInit(state, 2*digestSize, 1600 - 2*digestSize, Keccak_p);
    return INIT_SUCCESS;
}

void SHA3Update( Sha3Engine* state, BYTE* inBuf, unsigned size ) {
    SpoongeAbsorb(state, inBuf, size);
}

void SHA3Final( Sha3Engine* state, BYTE* inBuf, unsigned size ) {
    SpoongeAbsorb(state, inBuf, size);
    
    BYTE* statePtr = (BYTE*)state->state;

    /* SHA3 padding */
    statePtr[state->delayed] ^= 0x06;
    statePtr[state->r-1] ^= 0x80;
    
    state->keccak_p(&state->state);
}

void SHA3GetDigest( Sha3Engine* state, BYTE* digest ) {
    memcpy(digest, state->state, state->c/2);
}

int SHA3Sum( unsigned digestSize, BYTE* inBuf, unsigned size, BYTE* digest ) {
    Sha3Engine sha3;
    int r = SHA3Init(&sha3, digestSize);
    SHA3Final(&sha3, inBuf, size);
    SHA3GetDigest(&sha3, digest);
    return r;
}