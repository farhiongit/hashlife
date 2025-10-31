// Big interger toy library
#include <stdlib.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <locale.h>
#include <printf.h>

#include "bitl.h"

//---------------------------------------------------
// unsigned long long int
//---------------------------------------------------
uintbig_t UINTBIG_ZERO = ULL_TO_ULLL (0);
uintbig_t UINTBIG_MAX = { {ULLONG_MAX, ULLONG_MAX, ULLONG_MAX, ULLONG_MAX}
};

int
uintbig_cmp (uintbig_t a, uintbig_t b)
{
  for (size_t i = UBI_LENGTH; i > 0; i--)
    if (a.array[i - 1] > b.array[i - 1])
      return 1;
    else if (a.array[i - 1] < b.array[i - 1])
      return -1;
  return 0;
}

int
uintbig_is_zero (uintbig_t a)
{
  for (size_t i = 0; i < UBI_LENGTH; i++)
    if (a.array[i])
      return 0;
  return 1;
}

uintbig_t
uintbig_swapbits (uintbig_t a)
{
  for (size_t j = 0; j < UBI_LENGTH; j++)
    a.array[j] = ~a.array[j];
  return a;
}

uintbig_t
uintbig_and (uintbig_t a, uintbig_t b)
{
  for (size_t j = 0; j < UBI_LENGTH; j++)
    a.array[j] &= b.array[j];
  return a;
}

uintbig_t
uintbig_or (uintbig_t a, uintbig_t b)
{
  for (size_t j = 0; j < UBI_LENGTH; j++)
    a.array[j] |= b.array[j];
  return a;
}

uintbig_t
uintbig_xor (uintbig_t a, uintbig_t b)
{
  for (size_t j = 0; j < UBI_LENGTH; j++)
    a.array[j] ^= b.array[j];
  return a;
}

uintbig_t
uintbig_add (uintbig_t a, uintbig_t b)
{
  uintbig_t carry = UINTBIG_ZERO;

  for (size_t i = 0; i < UBI_LENGTH; i++)
  {
    // if a.array[i] + b.array[i] <= ULLONG_MAX
    if (a.array[i] <= ULLONG_MAX - b.array[i])
      a.array[i] += b.array[i];
    else                        // a.array[i] > ULLONG_MAX - b.array[i]
    {
      a.array[i] -= (ULLONG_MAX - b.array[i]) + 1ULL;   // --> 0 <= a.array[i] <= ULLONG_MAX - 1
      if (i < UBI_LENGTH - 1)
        carry.array[i + 1] = 1ULL;
      /* else wrap around */
    }
    if (a.array[i] <= ULLONG_MAX - carry.array[i])
      a.array[i] += carry.array[i];
    else
    {
      a.array[i] -= (ULLONG_MAX - carry.array[i]) + 1ULL;
      if (i < UBI_LENGTH - 1)
        carry.array[i + 1] = 1ULL;
      /* else wrap around */
    }
  }
  return a;
}

uintbig_t
uintbig_sub (uintbig_t a, uintbig_t b)
{
  uintbig_t carry = UINTBIG_ZERO;

  for (size_t i = 0; i < UBI_LENGTH; i++)
  {
    if (a.array[i] >= b.array[i])
      a.array[i] -= b.array[i];
    else                        // if (a.array[i] < b.array[i])
    {
      a.array[i] += (ULLONG_MAX - b.array[i]) + 1ULL;
      if (i < UBI_LENGTH - 1)
        carry.array[i + 1] = 1ULL;
      /* else wrap around */
    }
    if (a.array[i] >= carry.array[i])
      a.array[i] -= carry.array[i];
    else
    {
      a.array[i] += (ULLONG_MAX - carry.array[i]) + 1ULL;
      if (i < UBI_LENGTH - 1)
        carry.array[i + 1] = 1ULL;
      /* else wrap around */
    }
  }
  return a;
}

uintbig_t
uintbig_shiftleft (uintbig_t a, size_t shift)
{
  if (uintbig_is_zero (a) || shift >= UBI_NB_BITS)
    return UINTBIG_ZERO;

  for (; shift >= ULL_NB_BITS; shift -= ULL_NB_BITS)
  {
    for (size_t j = UBI_LENGTH; j > 1; j--)
      a.array[j - 1] = a.array[j - 2];
    a.array[0] = 0;
  }

  static unsigned long long int ULLONG_LAST = ~((~0ULL) >> 1);
  // a << shift
  for (size_t i = 0; i < shift; i++)
  {
    if (a.array[UBI_LENGTH - 1] & ULLONG_LAST)
    {
      /* a is truncated */
    }
    for (size_t j = UBI_LENGTH; j > 1; j--)
    {
      a.array[j - 1] <<= 1;
      if (a.array[j - 2] & ULLONG_LAST)
        a.array[j - 1] |= 1ULL;
    }
    a.array[0] <<= 1;
  }
  return a;
}

uintbig_t
uintbig_shiftright (uintbig_t a, size_t shift)
{
  if (uintbig_is_zero (a) || shift >= UBI_NB_BITS)
    return UINTBIG_ZERO;

  for (; shift >= ULL_NB_BITS; shift -= ULL_NB_BITS)
  {
    for (size_t j = 0; j < UBI_LENGTH - 1; j++)
      a.array[j] = a.array[j + 1];
    a.array[UBI_LENGTH - 1] = 0;
  }

  static unsigned long long int ULLONG_LAST = ~((~0ULL) >> 1);
  // a >> shift
  for (size_t i = 0; i < shift; i++)
  {
    for (size_t j = 0; j < UBI_LENGTH - 1; j++)
    {
      a.array[j] >>= 1;
      if (a.array[j + 1] & 1ULL)
        a.array[j] |= ULLONG_LAST;
    }
    a.array[UBI_LENGTH - 1] >>= 1;
  }
  return a;
}

//---------------------------------------------------
// signed long long int
//---------------------------------------------------
intbig_t INTBIG_ZERO = (intbig_t) { {0}
};
intbig_t INTBIG_MAX = { {ULLONG_MAX, ULLONG_MAX, ULLONG_MAX, LLONG_MAX}
};
intbig_t INTBIG_MIN = { {~ULLONG_MAX, ~ULLONG_MAX, ~ULLONG_MAX, ~(typeof (INTBIG_MIN.array[0])) LLONG_MAX}
};

int
intbig_is_negative (intbig_t a)
{
  static unsigned long long int ULLONG_LAST = ~((~0ULL) >> 1);
  if (a.array[UBI_LENGTH - 1] & ULLONG_LAST)
    return 1;
  else
    return 0;
}

int
intbig_is_positive (intbig_t a)
{
  return !intbig_is_negative (a) && !intbig_is_zero (a);
}

intbig_t
intbig_swapbits (intbig_t a)
{
  for (size_t j = 0; j < UBI_LENGTH; j++)
    a.array[j] = ~a.array[j];
  return a;
}

intbig_t
intbig_opposite (intbig_t a)
{
  // By convention, the opposite of a is coded as ~(a-1).
  // and therefore, the opposite of ~a+1 is coded as a.
  // As the opposite of the opposite is neutral, and opposition is bijective, ~a+1 is coded as the opposite of a.
  if (intbig_is_negative (a))
    return uintbig_add (intbig_swapbits (a), ULL_TO_ULLL (1));
  else                          // also works if a is equal to 0.
    return intbig_swapbits (uintbig_sub (a, ULL_TO_ULLL (1)));
}

intbig_t
intbig_abs (intbig_t a)
{
  if (intbig_is_negative (a))
    return uintbig_add (intbig_swapbits (a), ULL_TO_ULLL (1));
  else
    return a;
}

intbig_t
LL_TO_LLL (signed long long int ll)
{
  if (ll >= 0)
    return ULL_TO_ULLL ((long long unsigned int) ll);
  else
    return intbig_opposite (ULL_TO_ULLL ((long long unsigned int) (-ll)));
}

intbig_t
ULLL_TO_LLL (uintbig_t ua)
{
  // Might truncate.
  ua.array[UBI_LENGTH - 1] &= ~(0ULL) >> 1;
  return ua;
}

int
intbig_cmp (intbig_t a, intbig_t b)
{
  int s;
  if ((s = intbig_is_negative (b)) != intbig_is_negative (a))
    return s ? 1 : -1;
  else
    return uintbig_cmp (a, b);  // Also works if a and b are negative.
}

int
intbig_is_zero (intbig_t a)
{
  return uintbig_is_zero (a);
}

intbig_t
intbig_add (intbig_t a, intbig_t b)
{
  return uintbig_add (a, b);
}

intbig_t
intbig_sub (intbig_t a, intbig_t b)
{
  return uintbig_sub (a, b);
}

//---------------------------------------------------
// printf extension (%V for intbig_t, %U for uintbig_t)
//---------------------------------------------------

// see register_printf_specifier
// https://sourceware.org/bugzilla/attachment.cgi?id=3874&action=view
// http://www.gnu.org/software/libc/manual/html_node/Customizing-Printf.html#Customizing-Printf

static int PA_INTBIG;
static wchar_t SPEC_INTBIG = L'V';

static int PA_UINTBIG;
static wchar_t SPEC_UINTBIG = L'U';

#define DEFINE_PRINTER(name, printer)                         \
static int                                                    \
uintbig_to_##name (FILE * f, uintbig_t a, char flag)          \
{                                                             \
  (void) f;                                                   \
  char fmt0[] = "%0llu";                                      \
  char fmt1[] = "2^%u";                                       \
  char fmtn[] = "%0llu x 2^%u";                               \
  if (flag == '\'')                                           \
    fmtn[1] = fmt0[1] = flag;                                 \
                                                              \
  int ret = 0;                                                \
  ret += printer (f, "(");                                    \
  int more = 0;                                               \
  for (size_t j = UBI_LENGTH; j > 1; j--)                     \
    if (a.array[j - 1])                                       \
    {                                                         \
      if (more)                                               \
        ret += printer (f, " + ");                            \
      if (a.array[j - 1] == 1)                                \
        ret += printer (f, fmt1, (j - 1) * ULL_NB_BITS);      \
      else                                                    \
        ret += printer (f, fmtn, a.array[j - 1], (j - 1) * ULL_NB_BITS);  \
      more = 1;                                               \
    }                                                         \
                                                              \
  if (!more || a.array[0])                                    \
  {                                                           \
    if (more)                                                 \
      ret += printer (f, " + ");                              \
    ret += printer (f, fmt0, a.array[0]);                     \
  }                                                           \
  ret += printer (f, ")");                                    \
                                                              \
  return ret;                                                 \
}                                                             \
                                                              \
static int                                                    \
intbig_to_##name (FILE * f, intbig_t a, char flag_sign, char flag_group)  \
{                                                             \
  (void) f;                                                   \
  int ret = 0;                                                \
                                                              \
  if (intbig_is_negative (a))                                 \
  {                                                           \
    ret += printer (f, "%c", '-');                            \
    a = intbig_opposite (a);                                  \
  }                                                           \
  else if (intbig_is_zero (a) || flag_sign == ' ')            \
    ret += printer (f, " ");                                  \
  else if (flag_sign == '+')                                  \
    ret += printer (f, "+");                                  \
                                                              \
  ret += uintbig_to_##name (f, a, flag_group);                \
                                                              \
  return ret;                                                 \
}                                                             \
struct __useless_struct_##name##_IMPL

// The macro DEFINE_PRINTER avoids the use of a pointer to a variadic function fsinkf (file, ...) = snprinf (0, 0, ...)
DEFINE_PRINTER (stream, fprintf);
#define FSINKF(file, ...) snprintf (0, 0, __VA_ARGS__)
DEFINE_PRINTER (sink, FSINKF);

static int
intbig_printf (FILE *stream, const struct printf_info *info, const void *const args[])
{
  // info->spec == SPEC_INTBIG
  intbig_t a;
  memcpy (&a, *((void *const *const *) args)[0], sizeof (a));   // uh !
  int ret = intbig_to_sink (0, a, info->showsign ? '+' : info->space ? ' ' : 0,
                            info->group ? '\'' : 0);
  char *str = calloc ((size_t) (ret + 1), sizeof (*str));
  FILE *buffer = fmemopen (str, (size_t) (ret + 1), "w");
  intbig_to_stream (buffer, a, info->showsign ? '+' : info->space ? ' ' : 0, info->group ? '\'' : 0);
  fclose (buffer);
  ret = fprintf (stream, "%*s", (info->left ? -info->width : info->width), str);
  free (str);
  return ret;
}

static int
uintbig_printf (FILE *stream, const struct printf_info *info, const void *const args[])
{
  // info->spec == SPEC_UINTBIG
  uintbig_t a;
  memcpy (&a, *((void *const *const *) args)[0], sizeof (a));   // uh ! again
  int ret = uintbig_to_sink (0, a, info->group ? '\'' : 0);
  char *str = calloc ((size_t) (ret + 1), sizeof (*str));
  FILE *buffer = fmemopen (str, (size_t) (ret + 1), "w");
  uintbig_to_stream (buffer, a, info->group ? '\'' : 0);
  fclose (buffer);
  ret = fprintf (stream, "%*s", (info->left ? -info->width : info->width), str);
  free (str);
  return ret;
}

static int
intbig_print_arginfo_size (const struct printf_info *info, size_t n, int *argtypes, int *size)
{
  (void) info;
  /* We always take exactly one argument and this is big integer */
  if (n > 0)
    argtypes[0] = PA_INTBIG;
  *size = (int) sizeof (intbig_t);
  return 1;
}

static int
uintbig_print_arginfo_size (const struct printf_info *info, size_t n, int *argtypes, int *size)
{
  (void) info;
  /* We always take exactly one argument and this is big integer */
  if (n > 0)
    argtypes[0] = PA_UINTBIG;
  *size = (int) sizeof (uintbig_t);
  return 1;
}

static void
uintbig_va (void *mem, va_list *ap)
{
  uintbig_t v = va_arg (*ap, uintbig_t);
  memcpy (mem, &v, sizeof (v));
}

static void
intbig_va (void *mem, va_list *ap)
{
  intbig_t v = va_arg (*ap, intbig_t);
  memcpy (mem, &v, sizeof (v));
}

void
xintbig_printf_init (void)
{
  static int XINTBIG_PRINTF_INIT_DONE;

  if (XINTBIG_PRINTF_INIT_DONE)
    return;

  XINTBIG_PRINTF_INIT_DONE = 1;

  PA_INTBIG = register_printf_type (intbig_va);
  PA_UINTBIG = register_printf_type (uintbig_va);
  register_printf_specifier (SPEC_UINTBIG, uintbig_printf, uintbig_print_arginfo_size);
  register_printf_specifier (SPEC_INTBIG, intbig_printf, intbig_print_arginfo_size);
}

//---------------------------------------------------
// Unit testing
//---------------------------------------------------
#ifdef TUEXE
int
main (void)
{
  setlocale (LC_ALL, "");

  xintbig_printf_init ();

  printf ("%'llu;%'llu\n", 0ULL, ULLONG_MAX);
  printf ("%+'lli;%+'lli\n", LLONG_MIN, ~LLONG_MIN);
  printf ("%+'lli;%+'lli\n", LLONG_MAX, ~LLONG_MAX);
  printf ("%+'lli;%+'lli\n", 0LL, ~0LL);
  printf ("%+'lli;%+'lli\n", 10LL, ~10LL);
  printf ("%+'lli;%+'lli\n", -10LL, ~-10LL);

#  define PU(a) printf ("[%'100" PRIUINTBIG "]\n", a)
  uintbig_t a = UINTBIG_ZERO;
  a = ULL_TO_ULLL (1);

  PU (a);
  PU (a = uintbig_shiftleft (a, 70 + 64));
  PU (a = uintbig_shiftright (a, 69 + 64));

  uintbig_t b = ULL_TO_ULLL (~0ULL);
  PU (b);
  PU (a = uintbig_add (a, b));
  PU (a = uintbig_sub (a, b));

  PU ((a = (uintbig_t)
       {
       {
       ~0ULL, ~0ULL, 1, 0}}));
  PU (uintbig_add (a, ULL_TO_ULLL (1)));

#  define PS(a) printf ("[%'+100" PRIINTBIG "]\n", a)
  intbig_t sa;

  printf ("%+'V -- %+'V\n", INTBIG_MIN, INTBIG_MAX);

  PS (uintbig_sub (ULL_TO_ULLL (17), ULL_TO_ULLL (117)));
  PS (a = uintbig_sub (ULL_TO_ULLL (ULLONG_MAX), ULL_TO_ULLL (10058117ULL + LLONG_MAX)));
  PS (a = uintbig_sub (a, ULL_TO_ULLL (10058117ULL + (ULLONG_MAX - LLONG_MAX))));

  PS (sa = INTBIG_ZERO);
  PS (sa = intbig_opposite (sa));

  PS (sa = LL_TO_LLL (-1LL));
  PS (sa = intbig_opposite (sa));
  PS (sa = intbig_opposite (sa));

  PS (sa = LL_TO_LLL (1LL << 62));
  PS (sa = intbig_opposite (sa));
  PS (sa = intbig_opposite (sa));

  printf ("%i\n", intbig_cmp (INTBIG_ZERO, INTBIG_ZERO));
  printf ("%i\n", intbig_cmp (INTBIG_ZERO, LL_TO_LLL (1LL)));
  printf ("%i\n", intbig_cmp (LL_TO_LLL (1LL), INTBIG_ZERO));
  printf ("%i\n", intbig_cmp (INTBIG_ZERO, intbig_opposite (LL_TO_LLL (1LL))));
  printf ("%i\n", intbig_cmp (intbig_opposite (LL_TO_LLL (1LL)), INTBIG_ZERO));

  intbig_t sb = LL_TO_LLL (1LL << 62);
  PS (sb);
  PS (sa = intbig_add (sa, sb));
  PS (sa = intbig_sub (sa, sb));
  PS (sa = intbig_sub (sa, sb));
  PS (sa = intbig_sub (sa, sb));
  PS (sa = intbig_add (sa, sb));
  PS (sa = intbig_add (sa, sb));
  PS (sa = intbig_add (sa, LL_TO_LLL (-187496325)));
  PS (intbig_sub (ULLL_TO_LLL (uintbig_shiftleft (intbig_abs (sa), 151)), LL_TO_LLL (14789)));
  PS (intbig_opposite (intbig_sub (ULLL_TO_LLL (uintbig_shiftleft (intbig_abs (sa), 151)), LL_TO_LLL (14789))));
}
#endif
