#include "eced.h"

#define MAX_U64    0xFFFFFFFFFFFFFFFF
#define MSB_M      0x8000000000000000
#define HEX_FORMAT "%.16llX"

static const u64 unity[] = { 0x1, 0, 0, 0, 0, 0 };

static const EcPoint uP = { { 0, 0, 0, 0, 0, 0 }, { 1, 0, 0, 0, 0, 0 } };

static const u64 p192[] = { 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFE, 0xFFFFFFFFFFFFFFFF };

static const u64 p224[] = { 0x0000000000000001, 0xFFFFFFFF00000000, 0xFFFFFFFFFFFFFFFF, 0x00000000FFFFFFFF };

static const u64 p256[] = { 0xFFFFFFFFFFFFFFFF, 0x00000000FFFFFFFF, 0x0000000000000000, 0xFFFFFFFF00000001 };

static const u64 p384[] = { 0x00000000FFFFFFFF, 0xFFFFFFFF00000000, 0xFFFFFFFFFFFFFFFE, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF };