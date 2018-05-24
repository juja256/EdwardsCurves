#include "sha3.h"
#include <string.h>

static inline void xor(BYTE* acc, BYTE* other, unsigned size) {
    for (unsigned i=0; i<size; i++) {
        acc[i] ^= other[i];
    }
}

void Keccak_p( KeccakState* state ) {
    for (int i=0; i<24; i++) {
        Keccak_theta(state);
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

    if (size == 0) {
        spoonge->delayed += delayedSize;
        return;
    }

    for (int i=0; i<size / spoonge->r; i++) {
        xor(statePtr, inBuf + i*spoonge->r, spoonge->r);
        spoonge->keccak_p(spoonge->state);
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
        spoonge->keccak_p(spoonge->state);
        int toBeCopied = (i == restSize-1) ? size % spoonge->r : spoonge->r; 
        memcpy(outBuf + firstSize + i*spoonge->r, spoonge->state, toBeCopied );
    }    
}

void SHA3Init( SHA3State* state, unsigned digestSize);
void SHA3Update( SHA3State* state, BYTE inBuf, unsigned size );
void SHA3Final( SHA3State* state, BYTE inBuf, unsigned size );
void SHA3GetDigest( SHA3State* state, BYTE* digest );
void SHA3Sum( unsigned digestSize, BYTE* inBuf, unsigned size, BYTE* digest );