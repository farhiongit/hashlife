#include "bitl.h"
#define TYPE intbig_t
#define FORMAT "%'+V"
#define EQ(a, b) (intbig_cmp ((a), (b)) == 0)
#define LT(a, b) (intbig_cmp ((a), (b)) < 0)
#define ZERO LL_TO_LLL (0)
#define ONE LL_TO_LLL (1)
#define ADD(a, b) intbig_add ((a), (b))
#define MINUS(a) intbig_opposite (a)

