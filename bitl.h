// Big interger toy library
#ifndef __BITL_H__
#  define __BITL_H__

#  include <stdlib.h>

// Big integers are UBI_NB_BITS bits long.
#  define UBI_NB_BITS 256
#  define ULL_NB_BITS (sizeof (unsigned long long int) * 8)
#  define UBI_LENGTH (UBI_NB_BITS / ULL_NB_BITS)

//---------------------------------------------------
// unsigned long long int
//---------------------------------------------------
typedef struct
{
  unsigned long long int array[UBI_LENGTH];
} uintbig_t;

// Converts an unsigned long long integer into an unsigned big integer.
#  define ULL_TO_ULLL(ull) (uintbig_t){ { ull } }
// Converts an unsigned big integer into an unsigned long long integer (with possible loss).
#  define ULLL_TO_ULL(ulll) ((ulll).array[0])
// Zero
extern uintbig_t UINTBIG_ZERO;
// The maximal unsigned big integer.
extern uintbig_t UINTBIG_MAX;

// Compares two unsigned big integers.
// Returns a negative, null or positive integer whether a is lower, equal or higher than b.
int uintbig_cmp (uintbig_t a, uintbig_t b);

// Returns 1 if unsigned big integer is equal to 0.
int uintbig_is_zero (uintbig_t a);

// Returns ~a.
uintbig_t uintbig_swapbits (uintbig_t a);

// Returns a & b.
uintbig_t uintbig_and (uintbig_t a, uintbig_t b);

// Returns a | b.
uintbig_t uintbig_or (uintbig_t a, uintbig_t b);

// Returns a ^ b.
uintbig_t uintbig_xor (uintbig_t a, uintbig_t b);

// Returns a + b.
uintbig_t uintbig_add (uintbig_t a, uintbig_t b);

// Returns a - b.
uintbig_t uintbig_sub (uintbig_t a, uintbig_t b);

// Returns a << shift.
uintbig_t uintbig_shiftleft (uintbig_t a, size_t shift);

// Returns a >> shift.
uintbig_t uintbig_shiftright (uintbig_t a, size_t shift);

//---------------------------------------------------
// signed long long int
//---------------------------------------------------
typedef uintbig_t intbig_t;

// 0.
extern intbig_t INTBIG_ZERO;
// The minimal signed big integer.
extern intbig_t INTBIG_MIN;
// The maximal signed big integer.
extern intbig_t INTBIG_MAX;
// Converts a signed big integer into a long long integer (with possible loss).
#  define LLL_TO_LL(lll) ((lll).array[0])

// Returns 1 if a is negative.
int intbig_is_negative (intbig_t a);

// Returns 1 if a is positive.
int intbig_is_positive (intbig_t a);

// Returns ~a.
intbig_t intbig_swapbits (intbig_t a);

// Returns -a
intbig_t intbig_opposite (intbig_t a);

// Returns |a|.
intbig_t intbig_abs (intbig_t a);

// Converts a long long integer into a signed big integer.
intbig_t LL_TO_LLL (signed long long int ll);

// Converts an unsigned big integer into a signed big integer (with possible loss).
intbig_t ULLL_TO_LLL (uintbig_t ua);

// Compares two signed big integers.
// Returns a negative, null or positive integer whether a is lower, equal or higher than b.
int intbig_cmp (intbig_t a, intbig_t b);

// Returns 1 if a is zero.
int intbig_is_zero (intbig_t a);

// Returns a + b.
intbig_t intbig_add (intbig_t a, intbig_t b);

// Returns a - b.
intbig_t intbig_sub (intbig_t a, intbig_t b);

//---------------------------------------------------
// printf extension (%V for intbig_t, %U for uintbig_t)
//---------------------------------------------------
#  define PRIINTBIG "V"
#  define PRIUINTBIG "U"

// A call to this function prepares fprintf family of functions for formatted output to stream.
// Use conversion specifier U for an unsigned big integer, V for a signed big integer. Flag characters + and ' are supported.
// Big integers are formatted in base 2^64.
// E.g.                     : printf ("%+'V ; %+'V\n", INTBIG_MIN, INTBIG_MAX);
// yields (in french locale): -(9 223 372 036 854 775 808 x 2^192) ; +(9 223 372 036 854 775 807 x 2^192 + 18 446 744 073 709 551 615 x 2^128 + 18 446 744 073 709 551 615 x 2^64 + 18 446 744 073 709 551 615)
void xintbig_printf_init (void);

#endif
