// Hashlife Conway's game of life implementation
// This code implements the HashLife algorithms described by R. Gosper in 1984 in his article
// [GOSPER] Exploiting regularities in large cellular spaces, R.Wm. GOSPER, Physica D (Nonlinear Phenomena, volume 10) (1984) 75--80, North-Holland, Amsterdam.

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <locale.h>
#include <string.h>
#include <limits.h>
#include <regex.h>
#include <ctype.h>
#include <unistd.h>
#include "bitl.h"
#include "hgolbi.h"

#ifdef DEBUG
#  define ASSERT(...) do { if (!(__VA_ARGS__)) { fprintf (stderr, "Asserton is false at line %i.\n", __LINE__); abort(); } } while (0)
#  define PRINT(...) do { fprintf (stderr, __VA_ARGS__); } while (0)
#else
#  define ASSERT(...)
#  define PRINT(...)
#endif

static unsigned char NB_BITS[UINT16_MAX + 1] = { 0 };

// The rules of the game of life are implemented in this function (and nowhere else).
// [GOSPER] "Life is a two state, nine-neighborhood rule applied on an ordinary, two-dimensional grid"
//          "(...) a 4 by 4 cell [which] doesn't know its RESULT, (...) computes it by brute force,
//           , i.e. by applying the Life rule [4] to the nine-neighborhoods of each of its four central cells (bits)."
static unsigned char
next2x2 (uint16_t field4x4, unsigned char S, unsigned char B)
{
  /* *INDENT-OFF* */
  static struct { uint16_t N/*eighbours*/, C/*entral*/; } F[4] =
    {
     {.N = 1 << 0 | 1 << 1 | 1 << 2 | 1 << 4 | 1 <<   6 | 1 <<   8 | 1 <<   9 | 1 << 0xa, .C = 1 <<   5},  // -> 0
     {.N = 1 << 1 | 1 << 2 | 1 << 3 | 1 << 5 | 1 <<   7 | 1 <<   9 | 1 << 0xa | 1 << 0xb, .C = 1 <<   6},  // -> 1
     {.N = 1 << 4 | 1 << 5 | 1 << 6 | 1 << 8 | 1 << 0xa | 1 << 0xc | 1 << 0xd | 1 << 0xe, .C = 1 <<   9},  // -> 2
     {.N = 1 << 5 | 1 << 6 | 1 << 7 | 1 << 9 | 1 << 0xb | 1 << 0xd | 1 << 0xe | 1 << 0xf, .C = 1 << 0xa},  // -> 3
    };
  /* *INDENT-ON* */

  /* *INDENT-OFF* */
  /*
                 (xxxx)    (0123)
          (01)   (x01x)    (4567)
     hr = (23) = (x23x) = f(89ab)
                 (xxxx)    (cdef)
  */
  /* *INDENT-ON* */
  unsigned char hr = 0;
  for (size_t i = 0; i < sizeof (F) / sizeof (*F); i++)
    if (((field4x4 & F[i].C) ?  // Is the central cell of the all nine field (5, 6, 9 or a) alive ?
         S :                    // If alive, then, might survive ;
         B                      // otherwise, might generate.
        ) & (1 << NB_BITS[field4x4 & F[i].N]))
      hr |= (unsigned char) (1 << i);

  return hr;
}

static void
global_init (void)
{
  if (NB_BITS[1])               // NB_BITS[1] is non zero after global initialization, 0 before.
    return;

  NB_BITS[0] = 0;
  for (size_t i = 1; i < sizeof (NB_BITS) / sizeof (*NB_BITS); i++)
    NB_BITS[i] = (unsigned char) (NB_BITS[i >> 1] + (i & 1));
}

/*
    0--> x  (east)
    |
    v
    
    y
 (south)
*/
#define NB_QUADRANTS 4
typedef enum
{ NW = 0, NE = 1, SW = 2, SE = 3 } Quadrant;

// [GOSPER] "There are two key components - (1) a hash mechanism and (2) macro-cells"
// (2) Macro-cells
// [GOSPER] "A macro-cell represents a 2^n by 2^n block of (...) cells, where n is any non-negative integer."
//          "a macro-cell of size 2^n (n > 0) requires just five units.
//           These hold (pointers to) the four macro-cells of size 2^(n-1) which comprise the four quadrants,
//           and, if we are lucky (and n > 1), the RESULT, also of quadrant size.
//           The entire structure and evolution of an initial configuration [universe] will be encoded in the interlinkings
//           of macro-cells, which are computed as we probe its future."
typedef struct sMacrocell *MacrocellId;
struct sMacrocell
{
  uintbig_t population;         // Number of active cells in the macrocell.
  uintbig_t nb_instances;       // Macrocells are shared resources. This counter keeps track of the number of instances of macrocell in the universe (past, present and future).
  MacrocellId quadrant[NB_QUADRANTS];   // The four quadrants (NW, NE, SW, SE) that compose a macrocell.
  MacrocellId result;           // RESULT for the macrocell.
};

// (1) A hash mecanism. The Macrocell comparator is a perfect hash function for the hash-mecanism.
// Two macrocells compare equal if and only if they point to the same 4 quadrants
// (a == b) === !macrocell_lt (a, b) && !macrocell_lt (b, a)
static int
macrocell_lt (MacrocellId a, MacrocellId b)
{
  for (size_t i = 0; i < NB_QUADRANTS; i++)
    if (a->quadrant[i] == b->quadrant[i])
      continue;
    else if (a->quadrant[i] < b->quadrant[i])
      return 1;
    else                        //if (a->quadrant[i] > b->quadrant[i])
      return 0;

  return 0;
}

// [GOSPER] "At the bottom (...) are the 2^0 by 2^0 (i.e. 1 by 1) cells, of which there are at most two, since Life is a two state automaton"
// ON is a static leaf, shared between all universes.
#define ON (&_on)
static struct sMacrocell _on = { {{1}}, {{0}}, {0}, 0 };        // A leaf (cell is on)

#define QUERY (&_query)
static struct sMacrocell _query = { 0 };        // A fake leaf

static uintbig_t
macrocell_get_population (MacrocellId m, unsigned int depth)
{
  if (m == 0)
    return UINTBIG_ZERO;
  if (depth == 0)
    return m->population;

  uintbig_t nb = UINTBIG_ZERO;
  for (Quadrant q = 0; q < NB_QUADRANTS; q++)
    nb = uintbig_add (nb, macrocell_get_population (m->quadrant[q], depth - 1));

  return nb;
}

static int
xypos_lt (XYPos a, XYPos b)
{
  int dx = intbig_cmp (a.x, b.x);
  int dy = intbig_cmp (a.y, b.y);
  return dy < 0 || (dy == 0 && dx < 0);
}

#include "set_impl.h"
// A position in a 2D space
DECLARE_SET (XYPos);
DEFINE_OPERATORS (XYPos);
DEFINE_SET (XYPos);

typedef struct
{
  unsigned int height;
  uintbig_t xmin, ymin, tbase;
} SpaceTimeRegion;

DECLARE_SET (SpaceTimeRegion);
DEFINE_OPERATORS (SpaceTimeRegion);
DEFINE_SET (SpaceTimeRegion);

#ifdef DEBUG
static void
show_overlapping (unsigned int height, uintbig_t tbase, uintbig_t instant, uintbig_t xmin, uintbig_t ymin, Window window)
{
  uintbig_t size = uintbig_sl (ULL_TO_ULLL (1), height);
  uintbig_t quarter_size = uintbig_sr (size, 2);
  intbig_t sxmin, sxmax, symin, symax, diff;
  diff = uintbig_sub (UINTBIG_MAX, INTBIG_MAX);
  sxmin = uintbig_sub (xmin, diff);
  sxmax = intbig_sub (intbig_add (sxmin, size), LL_TO_LLL (1));
  symin = uintbig_sub (ymin, diff);
  symax = intbig_sub (intbig_add (symin, size), LL_TO_LLL (1));
  PRINT ("height %u:\n", height);
  PRINT (" (x) ([%+'V ; %+'V] ^ [%+'V ; %+'V]) ^\n", window.NWvertex.x, window.SEvertex.x, sxmin, sxmax);
  PRINT (" (y) ([%+'V ; %+'V] ^ [%+'V ; %+'V]) ^\n", window.NWvertex.y, window.SEvertex.y, symin, symax);
  PRINT (" (t) (%'U ^ [%'U ; %'U])\n", instant, tbase, uintbig_add (tbase, quarter_size));
}
#else
#  define show_overlapping(...)  do {} while (0)
#endif

static int
time_overlap (unsigned int height, uintbig_t tbase, uintbig_t instant)
{
  ASSERT (height);
  if (uintbig_cmp (instant, tbase) < 0)
    return 0;
  uintbig_t deltat = uintbig_sub (instant, tbase);
  uintbig_t quarter_size = height >= 2 ? uintbig_sl (ULL_TO_ULLL (1), height - 2) : UINTBIG_ZERO;
  if (uintbig_cmp (deltat, quarter_size) > 0)
    return 0;
  return 1;
}

static int
space_overlap (unsigned int height, uintbig_t xmin, uintbig_t ymin, Window window)
{
  uintbig_t wxmin = uintbig_add (uintbig_sub (UINTBIG_MAX, INTBIG_MAX), window.NWvertex.x);
  uintbig_t wxmax = uintbig_add (uintbig_sub (UINTBIG_MAX, INTBIG_MAX), window.SEvertex.x);
  uintbig_t wymin = uintbig_add (uintbig_sub (UINTBIG_MAX, INTBIG_MAX), window.NWvertex.y);
  uintbig_t wymax = uintbig_add (uintbig_sub (UINTBIG_MAX, INTBIG_MAX), window.SEvertex.y);
  uintbig_t xmax = uintbig_sub (uintbig_add (xmin, uintbig_sl (ULL_TO_ULLL (1), height)), ULL_TO_ULLL (1));
  uintbig_t ymax = uintbig_sub (uintbig_add (ymin, uintbig_sl (ULL_TO_ULLL (1), height)), ULL_TO_ULLL (1));

  if (uintbig_cmp (wxmax, xmin) < 0 || uintbig_cmp (wxmin, xmax) > 0 || uintbig_cmp (wymax, ymin) < 0 || uintbig_cmp (wymin, ymax) > 0)
    return 0;
  else
    return 1;
}

static uintbig_t
macrocell_get_cells_in_window (MacrocellId m, unsigned int height, uintbig_t xmin, uintbig_t ymin, Window window, SET (XYPos) *cells)
{
  if (!m)
    return UINTBIG_ZERO;

  if (!space_overlap (height, xmin, ymin, window))
    return UINTBIG_ZERO;

  uintbig_t ret = UINTBIG_ZERO;
  if (height)
  {
    uintbig_t half_size = uintbig_sl (ULL_TO_ULLL (1), height - 1);
    for (Quadrant q = 0; q < NB_QUADRANTS; q++)
      ret = uintbig_add (ret,
                         macrocell_get_cells_in_window (m->quadrant[q], height - 1,
                                                        uintbig_add (xmin, (q == NE || q == SE) ? half_size : UINTBIG_ZERO),
                                                        uintbig_add (ymin, (q == SW || q == SE) ? half_size : UINTBIG_ZERO), window, cells));
  }
  else
  {
    SET_ADD (cells, ((XYPos)
                     {
                     .x = (intbig_t) uintbig_sub (xmin, uintbig_sub (UINTBIG_MAX, INTBIG_MAX)),.y = (intbig_t) uintbig_sub (ymin, uintbig_sub (UINTBIG_MAX, INTBIG_MAX))}
             ));
    ret = ULL_TO_ULLL (1);
  }
  PRINT ("Collected %'U cells in macrocell.\n", ret);
  return ret;
}

  // Declare the use of template for sets (ordered lists) of MacrocellId.
#include "set_impl.h"
DECLARE_SET (MacrocellId);
DEFINE_OPERATORS (MacrocellId);
DEFINE_SET (MacrocellId);

static MacrocellId
macrocell_fetch_pattern (MacrocellId m, SET (MacrocellId) *set)
{
  if (!m)
    return 0;
  int empty = 1;
  for (Quadrant q = 0; q < NB_QUADRANTS && empty; q++)
    if (m->quadrant[q])
      empty = 0;

  if (empty)
  {
    free (m);
    return 0;
  }
  else
  {
    SNODE (MacrocellId) * node;
    node = SET_FIND (set, m);
    ASSERT (node);
    if (*SNODE_KEY (node) != m)
      free (m);
    return *SNODE_KEY (node);
  }
}

// [GOSPER] "The hash mechanism prevents the recomputation of indistinguishable scenarios."
//          "a macrocell is never created if one having the same quadrants already exists."
static MacrocellId
macrocell_patternify (MacrocellId m, SET (MacrocellId) *set)
{
  if (!m)
    return 0;
  int empty = 1;
  for (Quadrant q = 0; q < NB_QUADRANTS && empty; q++)
    if (m->quadrant[q])
      empty = 0;

  SNODE (MacrocellId) * node;
  if (empty)
  {
    // If the macrocell m is an empty region, it can be forgotten.
    free (m);
    return 0;
  }
  else if ((node = SET_FIND (set, m)))
  {
    // [GOSPER] "When the algorithm tries to group four quadrants to form a pre-existing macro-cell, the hash
    //           mechanism notices the coincidence and returns the old cell instead of a new one. Most importantly,
    //           this old cell may already know its RESULT"
    // Aggregation:
    // Is there another identical (twin) macrocell at the same height ?
    // If the macrocell is not an empty region but is already registered in the set, it can be forgotten: forget m.
    if (*SNODE_KEY (node) != m)
    {
      free (m);                 // The previously created macrocell m has a twin, it is not needed anymore and is destroyed.
      (*SNODE_KEY (node))->nb_instances = uintbig_add ((*SNODE_KEY (node))->nb_instances, ULL_TO_ULLL (1));     // The already existing pattern is used once more.
    }
    return *SNODE_KEY (node);   // The previously created macrocell m is aggregated with a preexisting one.
  }
  else                          // newmc is not anywhere else in the universe
  {
    // If the macrocell is not an empty region and is not registered yet, it must be registered in the set.
    // Add it from the path into the set of patterns.
    m->result = QUERY;          // Makes sure that an extraneous result is not registered.
    m->nb_instances = ULL_TO_ULLL (1);
    m->population = macrocell_get_population (m, 1);
    SET_ADD (set, m);           // register the macrocell in the ordered set.
    return m;
  }
}

// Declare the use of template for lists of Level.
// Register of macrocells for a given level.
typedef struct sLevel Level;
struct sLevel
{
  SET (MacrocellId) * macrocells;       // An ordered set of unique Macrocell references (defined as the address of the macrocells)
  // [GOSPER] "This will usually require fewer macro-cells than you might think, due to two restrictions on when a macro-cell
  //           can be created. First, a macrocell is never created if one having the same quadrants already exists."
  //           This applies recursively to the quadrants. At the bottom of the recursion are the 2^0 by 2^0 (i.e. 1 by 1) cells,
  //           of which there are at most two, since Life is a two state automaton."
  //          "The second restriction on macro-cell creation is implicit in the algorithm: a macro-cell of size 2^n,
  //           (n/> 2) can only be created when its x, y, and time coordinates (relative to all its parent macro-cells)
  //           are multiples of 2^(n-2). Thus, proliferations of cells is limited by indistinguishability when they are
  //           small, and by infrequency of creation when they are large."
};

#include "list_impl.h"
DECLARE_LIST (Level);
DEFINE_OPERATORS (Level);
DEFINE_LIST (Level);

// ---------------------- Universes
typedef struct sUniverse
{
  unsigned int height;
  uintbig_t x0, y0;
  MacrocellId root;             // The macrocell at the top of the universe
  /* *INDENT-OFF* */
  LIST (Level) * listOfLevels;   // Register of macrocells at every level of the Universe.
  /* *INDENT-ON* */
  unsigned char RESULT4x4[UINT16_MAX + 1];      // RESULT for 4x4 macro-cells for rule B/S.
  unsigned char S;              // Pattern for survival.
  unsigned char B;              // Pattern for birth.
} Universe;

static void
universe_init (Universe *pUniverse)
{
  global_init ();

  static Universe UNIVERSE_INITIALIZER = { 0 };
  *pUniverse = UNIVERSE_INITIALIZER;
  pUniverse->S = 1 << 2 | 1 << 3;       // Pattern for survival (if 2 ou 3 neighbours).
  pUniverse->B = 1 << 3;        // Pattern for birth (if 3 neighbours).
  for (size_t i = 0; i < sizeof (pUniverse->RESULT4x4) / sizeof (*(pUniverse->RESULT4x4)); i++)
    pUniverse->RESULT4x4[i] = next2x2 ((uint16_t) i, pUniverse->S, pUniverse->B);
}

Universe *
universe_create (void)
{
  xintbig_printf_init ();
  Universe *pUniverse = calloc (1, sizeof (*pUniverse));
  if (pUniverse)
    universe_init (pUniverse);
  return pUniverse;
}

void
universe_reinitialize (Universe *pUniverse)
{
  if (pUniverse->listOfLevels)
  {
    while (LIST_SIZE (pUniverse->listOfLevels))
    {
      SET (MacrocellId) * s = LNODE_VALUE (LIST_LAST (pUniverse->listOfLevels))->macrocells;
      while (SET_SIZE (s))
      {
        MacrocellId m = *SNODE_KEY (SET_LAST (s));
        if (LIST_SIZE (pUniverse->listOfLevels) > 1)    // Leafs are static or allocated by the caller of set_cell and should not be free'd here.
          free (m);
        SET_REMOVE (s, SET_LAST (s));
      }
      SET_DESTROY (s);
      LIST_REMOVE (pUniverse->listOfLevels, LIST_LAST (pUniverse->listOfLevels));
    }
    ASSERT (pUniverse->listOfLevels);
    LIST_DESTROY (pUniverse->listOfLevels);
    pUniverse->listOfLevels = 0;
  }
  universe_init (pUniverse);
}

void
universe_destroy (Universe *pUniverse)
{
  universe_reinitialize (pUniverse);
  free (pUniverse);
}

// [GOSPER] "Macro-cells mechanize the information-compression of the spacetime behavior of configurations."
//          "Each macro-cell seeks to determine its RESULT, namely the concentric 2^(n-1) by 2^(n-1) macro-cell
//           which the parent macro-cell exclusively determines after 2^(n-2) time steps."
// [GOSPER] "the algorithm is indifferent to the x, y, and time coordinates and even the sizes of the macro-cells.
//           Their spatial coordinates are implicit in the quadrant structure of their owners, and their
//           time coordinates are implicit in the RESULT structure."
static MacrocellId
universe_get_RESULT (Universe *pUniverse, MacrocellId m, unsigned int height)
{
  if (m == 0 || height < 2)
    return 0;

  // [GOSPER] "A macro-cell at time 0 is the top stratum of a patch of earth. Successively deeper strata hold future
  //           time-slices (...) down to the RESULT stratum, which is the flat bottom of a hole with sides of slope 1,
  //           with the initial macro-cell as its base. (the light cone of the initial macro-cell.)
  //           The depth of the hole is 1/4 of its width, and half the width of the bottom."
  // Compute RESULT, S/4 generations ahead of m, where S is the size of m.
  /* Macrocell m (.), of size S  ->  RESULT (H), of size S/2, S/4 generations ahead of m.

     ........                         ........
     ........                         ........
     ........                         ..HHHH..
     ........                    ->   ..HHHH..
     ........                         ..HHHH..
     ........                         ..HHHH..
     ........                         ........
     ........                         ........
   */

  int empty;
  empty = 1;
  for (Quadrant q = 0; q < NB_QUADRANTS && empty; q++)
    if (m->quadrant[q])
      empty = 0;
  if (empty)
    return m->result = 0;

  ASSERT (LIST_SIZE (pUniverse->listOfLevels) > height);
  Level *l0 = 0;
  Level *l1 = 0;
  Level *l2 = 0;

  l0 = LNODE_VALUE (LIST_INDEX (pUniverse->listOfLevels, height));
  ASSERT (l0 && l0->macrocells);
  ASSERT (m == macrocell_patternify (m, l0->macrocells));       // m should have been registered already.
  m = macrocell_patternify (m, l0->macrocells);
  // [GOSPER] "If the queried macro-cell already knows its RESULT (from having computed it previously), it just returns it."
  if (m->result != QUERY)       // If the RESULT was already computed, return it.
    return m->result;

  l1 = LNODE_VALUE (LIST_INDEX (pUniverse->listOfLevels, height - 1));
  ASSERT (l1 && l1->macrocells);
  MacrocellId result = 0;
  if (height == 2)
  {
    // [GOSPER] "The smallest macro-cell which can have a RESULT is 4 by 4."
    /*
       mmmm
       mHHm
       mHHm
       mmmm
     */

    uint16_t field4x4 = 0;
    for (Quadrant i = 0; i < NB_QUADRANTS; i++)
      if (m->quadrant[i])
        for (Quadrant j = 0; j < NB_QUADRANTS; j++)
          if (m->quadrant[i]->quadrant[j] == ON)
            field4x4 |= (uint16_t) (1 << ((i == NW ? 0 : i == NE ? 2 : i == SW ? 8 : 10) + (j == NW ? 0 : j == NE ? 1 : j == SW ? 4 : 5)));

    // [GOSPER] "If a 4 by 4 cell doesn't know its RESULT, it computes it by brute force."
    unsigned char hr2 = pUniverse->RESULT4x4[field4x4];

    result = calloc (1, sizeof (*result));      // Set all quadrants to the state 0.
    for (Quadrant i = 0; i < NB_QUADRANTS; i++)
      if ((hr2 >> i) & 1)
      {
        // [GOSPER] "At the bottom of.the recursion are the 2^0 by 2^0 (i.e. 1 by 1) cells, of which there are at most two,
        //           since Life [universe] is a two state automaton."
        result->quadrant[i] = ON;
        result->population = uintbig_add (result->population, ULL_TO_ULLL (1));
      }
  }                             // if (height == 2)
  else                          // if (height != 2)
  {
    // [GOSPER] "reuse of many RESULTs in the construction of larger RESULTs (which represent larger time-steps)."
    // [GOSPER] "Larger cells determine their RESULTs by a (...) recursion which involves combining a total of thirteen
    //           separate RESULTs of quadrants, and other quadrant-sized macro-cells formed by grouping
    //           RESULTS and regrouping quarants of quadrants."
    /* Macrocell m, of size S:
       ........                            ........
       ........                            .cccccc.
       ........                 cell[36]   .cccccc.
       ........                 36 cells   .cccccc.
       ........                            .cccccc.
       ........                            .cccccc.
       ........                            .cccccc.
       ........                            ........
     */
    MacrocellId cell[36] = { 0 };       // height - 3

    /* 1 to 4 */
    // [GOSPER] "By taking the RESULTs of the four quadrants, we reach half of the desired depth,
    //           but there remain dikes covering 5/9 of this halfway bottom."
    /* Get RESULT (h) of quadrants (q) of macrocell m.
       q is half the size (S/2) of m, therefore h is S/2/4 = S/8 generations ahead of m.

       qqqq....      ....qqqq      ........      ........
       qhhq....      ....qhhq      ........      ........
       qhhq....      ....qhhq      ........      ........
       qqqq....      ....qqqq      ........      ........
       ........      ........      qqqq....      ....qqqq
       ........      ........      qhhq....      ....qhhq
       ........      ........      qhhq....      ....qhhq
       ........      ........      qqqq....      ....qqqq
     */
    /* *INDENT-OFF* */
    struct
    {
      MacrocellId *cell[NB_QUADRANTS];
    } unit_1_4[NB_QUADRANTS] =
    {
      {{&cell[ 0], &cell[ 1], &cell[ 6], &cell[ 7]}},  // NW
      {{&cell[ 4], &cell[ 5], &cell[10], &cell[11]}},  // NE
      {{&cell[24], &cell[25], &cell[30], &cell[31]}},  // SW
      {{&cell[28], &cell[29], &cell[34], &cell[35]}}   // SE
    };
    /* *INDENT-ON* */

    for (Quadrant u = 0; u < NB_QUADRANTS; u++)
      if ((result = universe_get_RESULT (pUniverse, m->quadrant[u], height - 1)))       // height - 2
        for (Quadrant i = 0; i < NB_QUADRANTS; i++)
          *unit_1_4[u].cell[i] = result->quadrant[i];   // quadrants of result are already paternified.

    /* 5 to 9 */
    // [GOSPER] "To excavate these dikes, five artificial, shifted "quadrants" must be constructed from quadrants' quadrants,
    //           and then RESULTed. This will involve the reexcavation of some thin air, but at little cost, since
    //           the reexcavated cells will remember their RESULTs."
    /* Get RESULT (h) of quadrants (q) of macrocell m
       q is half the size (S/2) of m, therefore h is S/2/4 = S/8 generations ahead of m.

       ..qqqq..      ........      ........      ........      ........
       ..qhhq..      ........      ........      ........      ........
       ..qhhq..      ........      qqqq....      ....qqqq      ..qqqq..
       ..qqqq..      ........      qhhq....      ....qhhq      ..qhhq..
       ........      ..qqqq..      qhhq....      ....qhhq      ..qhhq..
       ........      ..qhhq..      qqqq....      ....qqqq      ..qqqq..
       ........      ..qhhq..      ........      ........      ........
       ........      ..qqqq..      ........      ........      ........
     */
    MacrocellId mtemp_5_9[5];   // height - 1
    for (unsigned int u = 0; u < 5; u++)
      mtemp_5_9[u] = calloc (1, sizeof (*mtemp_5_9[u]));
    /* *INDENT-OFF* */
    struct
    {
      MacrocellId *cell[NB_QUADRANTS];
    } unit_5_9[5] =
    {
      {{&cell[ 2], &cell[ 3], &cell[ 8], &cell[ 9]}},  // 5
      {{&cell[26], &cell[27], &cell[32], &cell[33]}},  // 6
      {{&cell[12], &cell[13], &cell[18], &cell[19]}},  // 7
      {{&cell[16], &cell[17], &cell[22], &cell[23]}},  // 8
      {{&cell[14], &cell[15], &cell[20], &cell[21]}}   // 9
    };
    /* *INDENT-ON* */
    /* 5 */
    if (m->quadrant[NW])
    {
      mtemp_5_9[0]->quadrant[NW] = m->quadrant[NW]->quadrant[NE];
      mtemp_5_9[0]->quadrant[SW] = m->quadrant[NW]->quadrant[SE];
    }
    if (m->quadrant[NE])
    {
      mtemp_5_9[0]->quadrant[NE] = m->quadrant[NE]->quadrant[NW];
      mtemp_5_9[0]->quadrant[SE] = m->quadrant[NE]->quadrant[SW];
    }
    /* 6 */
    if (m->quadrant[SW])
    {
      mtemp_5_9[1]->quadrant[NW] = m->quadrant[SW]->quadrant[NE];
      mtemp_5_9[1]->quadrant[SW] = m->quadrant[SW]->quadrant[SE];
    }
    if (m->quadrant[SE])
    {
      mtemp_5_9[1]->quadrant[NE] = m->quadrant[SE]->quadrant[NW];
      mtemp_5_9[1]->quadrant[SE] = m->quadrant[SE]->quadrant[SW];
    }
    /* 7 */
    if (m->quadrant[NW])
    {
      mtemp_5_9[2]->quadrant[NW] = m->quadrant[NW]->quadrant[SW];
      mtemp_5_9[2]->quadrant[NE] = m->quadrant[NW]->quadrant[SE];
    }
    if (m->quadrant[SW])
    {
      mtemp_5_9[2]->quadrant[SW] = m->quadrant[SW]->quadrant[NW];
      mtemp_5_9[2]->quadrant[SE] = m->quadrant[SW]->quadrant[NE];
    }
    /* 8 */
    if (m->quadrant[NE])
    {
      mtemp_5_9[3]->quadrant[NW] = m->quadrant[NE]->quadrant[SW];
      mtemp_5_9[3]->quadrant[NE] = m->quadrant[NE]->quadrant[SE];
    }
    if (m->quadrant[SE])
    {
      mtemp_5_9[3]->quadrant[SW] = m->quadrant[SE]->quadrant[NW];
      mtemp_5_9[3]->quadrant[SE] = m->quadrant[SE]->quadrant[NE];
    }
    /* 9 */
    if (m->quadrant[NW])
      mtemp_5_9[4]->quadrant[NW] = m->quadrant[NW]->quadrant[SE];
    if (m->quadrant[NE])
      mtemp_5_9[4]->quadrant[NE] = m->quadrant[NE]->quadrant[SW];
    if (m->quadrant[SW])
      mtemp_5_9[4]->quadrant[SW] = m->quadrant[SW]->quadrant[NE];
    if (m->quadrant[SE])
      mtemp_5_9[4]->quadrant[SE] = m->quadrant[SE]->quadrant[NW];

    /* 5 to 9 */
    for (unsigned int u = 0; u < 5; u++)
    {
      // The quadrants of mtemp_5_9[u] have been paternified already.
      mtemp_5_9[u] = macrocell_patternify (mtemp_5_9[u], l1->macrocells);

      if ((result = universe_get_RESULT (pUniverse, mtemp_5_9[u], height - 1))) // height - 2
        for (Quadrant i = 0; i < NB_QUADRANTS; i++)
          *unit_5_9[u].cell[i] = result->quadrant[i];   // quadrants of result are already paternified.
    }

    /* Therefore, we have 36 already paternified cells (h), S/8 generations ahead of m.
       ........
       .hhhhhh.
       .hhhhhh.
       .hhhhhh.
       .hhhhhh.
       .hhhhhh.
       .hhhhhh.
       ........
     */

    /* 10 to 13 */
    // [GOSPER] "We are now on the halfway bottom, composed of nine subresults, which are then grouped in fours
    //           to form four overlapping squares. The grouping of the RESULTs of these four squares is the grand RESULT."
    /* Get RESULT (H) of quadrants (x) of cells h
       x is half the size (S/2) of m, therefore H is S/2/4 = S/8 generations ahead of x, part of h, itself S/8 generations ahead of m.
       Therefore, H is S/8 + S/8 = S/4 generations ahead of m.

       ........      ........      ........      ........
       .xxxxhh.      .hhxxxx.      .hhhhhh.      .hhhhhh.
       .xHHxhh.      .hhxHHx.      .hhhhhh.      .hhhhhh.
       .xHHxhh.      .hhxHHx.      .xxxxhh.      .hhxxxx.
       .xxxxhh.      .hhxxxx.      .xHHxhh.      .hhxHHx.
       .hhhhhh.      .hhhhhh.      .xHHxhh.      .hhxHHx.
       .hhhhhh.      .hhhhhh.      .xxxxhh.      .hhxxxx.
       ........      ........      ........      ........
     */

    /* *INDENT-OFF* */
    struct
    {
      MacrocellId cell[NB_QUADRANTS /*j*/][NB_QUADRANTS /*k*/];
    } unit_10_13[NB_QUADRANTS] =
    {
      { // NW (j)
        {{cell[ 0], cell[ 1], cell[ 6], cell[ 7]},  // NW (k)
         {cell[ 2], cell[ 3], cell[ 8], cell[ 9]},  // NE
         {cell[12], cell[13], cell[18], cell[19]},  // SW
         {cell[14], cell[15], cell[20], cell[21]}}  // SE
      },
      { // NE
        {{cell[ 2], cell[ 3], cell[ 8], cell[ 9]},  // NW
         {cell[ 4], cell[ 5], cell[10], cell[11]},  // NE
         {cell[14], cell[15], cell[20], cell[21]},  // SW
         {cell[16], cell[17], cell[22], cell[23]}}  // SE
      },
      { // SW
        {{cell[12], cell[13], cell[18], cell[19]},  // NW
         {cell[14], cell[15], cell[20], cell[21]},  // NE
         {cell[24], cell[25], cell[30], cell[31]},  // SW
         {cell[26], cell[27], cell[32], cell[33]}}  // SE
      },
      { // SE
        {{cell[14], cell[15], cell[20], cell[21]},  // NW
         {cell[16], cell[17], cell[22], cell[23]},  // NE
         {cell[26], cell[27], cell[32], cell[33]},  // SW
         {cell[28], cell[29], cell[34], cell[35]}}  // SE
      }
    };
    /* *INDENT-ON* */

    l2 = LNODE_VALUE (LIST_INDEX (pUniverse->listOfLevels, height - 2));
    ASSERT (l2 && l2->macrocells);

    result = calloc (1, sizeof (*result));

    for (Quadrant u = 0; u < NB_QUADRANTS; u++)
    {
      MacrocellId mtemp = calloc (1, sizeof (*mtemp));
      for (Quadrant j = 0; j < NB_QUADRANTS; j++)
      {
        MacrocellId qtemp = calloc (1, sizeof (*qtemp));
        for (Quadrant k = 0; k < NB_QUADRANTS; k++)
          qtemp->quadrant[k] = unit_10_13[u].cell[j][k];        // Quadrants of qtemp are already paternified.
        mtemp->quadrant[j] = macrocell_patternify (qtemp, l2->macrocells);      // Quadrants of mtemp are paternified.
      }
      mtemp = macrocell_patternify (mtemp, l1->macrocells);     // mtemp is paternified.
      result->quadrant[u] = universe_get_RESULT (pUniverse, mtemp, height - 1); // Quadrants of result are paternified, at height - 2
    }

    /* Done: hresult H, half the size of m, is S/4 generations ahead of m:
       ........
       ........
       ..HHHH..
       ..HHHH..
       ..HHHH..
       ..HHHH..
       ........
       ........
     */
  }                             // if (height != 2)

  // Register result (suppose here that m has already been registered as a pattern)
  return m->result = macrocell_patternify (result, l1->macrocells);
}

// Test whether or not a universe is closed.
static int
universe_is_closed (Universe *pUniverse)
{
  /* A universe is NOT closed (therefore open) if surrouded by empty space (.):
     ........
     ........
     ..xxxx..
     ..xxxx..
     ..xxxx..
     ..xxxx..
     ........
     ........
   */
  if (pUniverse->root == 0)
    return 0;
  ASSERT (pUniverse->height != 0);
  MacrocellId r = pUniverse->root;
  if (pUniverse->height == 1)
  {
    for (Quadrant i = 0; i < NB_QUADRANTS; i++)
      if (r->quadrant[i])
        // Universe is not surrounded by empty space.
        return 1;
  }
  else if (pUniverse->height >= 2)
  {
    for (Quadrant i = 0; i < NB_QUADRANTS; i++)
      if (r->quadrant[i])
        for (Quadrant j = 1; j < NB_QUADRANTS; j++)
          if (r->quadrant[i]->quadrant[(NB_QUADRANTS - 1 - i + j) % NB_QUADRANTS])
            // Universe is not surrounded by empty space.
            return 1;
  }

  // Universe is surrounded by empty space.
  return 0;
}

static void
universe_expand (Universe *pUniverse)
{
  if (!pUniverse->root)
    return;
  ASSERT (pUniverse->height != 0);

  PRINT ("*** Expand universe from height %i to ", pUniverse->height);

  int do_not_expand = 1;
  for (Quadrant i = 0; i < NB_QUADRANTS && do_not_expand; i++)
    if (pUniverse->root->quadrant[i])
      do_not_expand = 0;
  if (do_not_expand)
    return;

  MacrocellId newroot = calloc (1, sizeof (*newroot));
  newroot->result = QUERY;
  newroot->nb_instances = ULL_TO_ULLL (1);      // The root is unique.
  newroot->population = pUniverse->root->population;
  for (Quadrant i = 0; i < NB_QUADRANTS; i++)
    if (pUniverse->root->quadrant[i])
    {
      newroot->quadrant[i] = calloc (1, sizeof (*newroot->quadrant[i]));
      newroot->quadrant[i]->quadrant[NB_QUADRANTS - 1 - i] = pUniverse->root->quadrant[i];
      newroot->quadrant[i]->nb_instances = ULL_TO_ULLL (1);     // Each of the 4 quadrants newroot->quadrant[i] is unique
      newroot->quadrant[i]->result = QUERY;
      newroot->quadrant[i]->population = pUniverse->root->quadrant[i]->population;
    }

  ASSERT (LIST_SIZE (pUniverse->listOfLevels) >= pUniverse->height + 1);
  Level *pl = LNODE_VALUE (LIST_INDEX (pUniverse->listOfLevels, pUniverse->height));
  ASSERT (pl && pl->macrocells);
  SNODE (MacrocellId) * node = SET_FIND (pl->macrocells, pUniverse->root);
  ASSERT (node && *SNODE_KEY (node) == pUniverse->root);
  SET_REMOVE (pl->macrocells, node);
  for (Quadrant i = 0; i < NB_QUADRANTS; i++)
    if (newroot->quadrant[i])
      newroot->quadrant[i] = macrocell_patternify (newroot->quadrant[i], pl->macrocells);
  free (pUniverse->root);
  pUniverse->x0 = uintbig_sub (pUniverse->x0, uintbig_sl (ULL_TO_ULLL (1), pUniverse->height - 1));
  pUniverse->y0 = uintbig_sub (pUniverse->y0, uintbig_sl (ULL_TO_ULLL (1), pUniverse->height - 1));
  pUniverse->height++;
  pUniverse->root = newroot;
  if (LIST_SIZE (pUniverse->listOfLevels) == pUniverse->height)
  {
    Level l = {
      0
    };
    l.macrocells = SET_CREATE (MacrocellId, macrocell_lt);
    LIST_APPEND (pUniverse->listOfLevels, l);
  }
  pl = LNODE_VALUE (LIST_INDEX (pUniverse->listOfLevels, pUniverse->height));
  ASSERT (pl && pl->macrocells);
  pUniverse->root = macrocell_patternify (pUniverse->root, pl->macrocells);

  PRINT ("%i.\n", pUniverse->height);
}

static int
universe_contains (Universe *pUniverse, intbig_t sx, intbig_t sy)
{
  if (!pUniverse->root)
    return 0;
  uintbig_t x = uintbig_add (uintbig_sub (UINTBIG_MAX, INTBIG_MAX), sx);
  uintbig_t y = uintbig_add (uintbig_sub (UINTBIG_MAX, INTBIG_MAX), sy);
  if (uintbig_cmp (x, pUniverse->x0) < 0 || uintbig_cmp (y, pUniverse->y0) < 0)
    return 0;
  if (!uintbig_is_zero (uintbig_sr (uintbig_sub (x, pUniverse->x0), pUniverse->height)) || !uintbig_is_zero (uintbig_sr (uintbig_sub (y, pUniverse->y0), pUniverse->height)))
    return 0;
  return 1;
}

// [GOSPER] "The entire structure and evolution of an initial configuration will be encoded in the interlinkings of macro-cells"
static MacrocellId
universe_cell_accessor (Universe *pUniverse, intbig_t sx, intbig_t sy, MacrocellId leaf)
{
  uintbig_t x = uintbig_add (uintbig_sub (UINTBIG_MAX, INTBIG_MAX), sx);
  uintbig_t y = uintbig_add (uintbig_sub (UINTBIG_MAX, INTBIG_MAX), sy);
  MacrocellId oldleaf = 0;
  ASSERT (!leaf || (!leaf->quadrant[NW] && !leaf->quadrant[NE] && !leaf->quadrant[SW] && !leaf->quadrant[SE])); // leaf is a leaf
  if ((leaf == 0 || leaf == QUERY) && !universe_contains (pUniverse, sx, sy))   // Cannot remove a cell from or query an unexisting position
    return 0;
  // Insert the cell in the tree
  else if (pUniverse->root == 0)
  {
    // If the universe is empty, create it ex nihilo.
    ASSERT (leaf);
    pUniverse->height = 1;
    pUniverse->x0 = uintbig_sl (uintbig_sr (x, pUniverse->height), pUniverse->height);
    pUniverse->y0 = uintbig_sl (uintbig_sr (y, pUniverse->height), pUniverse->height);
    Quadrant q = ULLL_TO_ULL (uintbig_add (uintbig_and (uintbig_sr (x, pUniverse->height - 1), ULL_TO_ULLL (1)),
                                           uintbig_sl (uintbig_and (uintbig_sr (y, pUniverse->height - 1),
                                                                    ULL_TO_ULLL (1)), 1)));
    pUniverse->root = calloc (1, sizeof (*pUniverse->root));
    pUniverse->root->nb_instances = ULL_TO_ULLL (1);    // The root is unique.
    pUniverse->root->quadrant[q] = leaf;
    pUniverse->root->result = QUERY;    // Unneeded as height < 2
    pUniverse->root->population = ULL_TO_ULLL (1);
    ASSERT (!pUniverse->listOfLevels);
    pUniverse->listOfLevels = LIST_CREATE (Level);
    Level l0 = {
      0
    };
    l0.macrocells = SET_CREATE (MacrocellId, macrocell_lt);
    SET_ADD (l0.macrocells, leaf);
    Level l1 = {
      0
    };
    l1.macrocells = SET_CREATE (MacrocellId, macrocell_lt);
    pUniverse->root = macrocell_patternify (pUniverse->root, l1.macrocells);
    ASSERT (LIST_SIZE (pUniverse->listOfLevels) == 0);
    LIST_APPEND (pUniverse->listOfLevels, l0);
    LIST_APPEND (pUniverse->listOfLevels, l1);
  }                             // if (!pUniverse->root)
  else                          // if (pUniverse->root)
  {
    // Phase 1: Universe expansion. Expand the universe to include the cell at position (x, y)
    // Build tree upward, to widen the universe to include (x, y)
    while (!universe_contains (pUniverse, sx, sy))
      universe_expand (pUniverse);

    ASSERT (pUniverse->height > 0);
    ASSERT (LIST_SIZE (pUniverse->listOfLevels) >= pUniverse->height + 1);
    // Phase 2: Find the target leaf in the tree.
    struct
    {
      MacrocellId oldmc, newmc;
      Quadrant q;               // quadrant q of oldmc which is modified
    } *path = calloc (pUniverse->height, sizeof (*path));
    // Top to bottom : create quadrants down to the position (x, y), starting from the root at the top.
    MacrocellId m = pUniverse->root;    // (*) At step 0, the macrocell m is the root. It is unique by construction.
    for (size_t h = pUniverse->height; h >= 1; h--)
    {
      ASSERT (h >= 1);
      Quadrant q = ULLL_TO_ULL (uintbig_add (uintbig_and (uintbig_sr (uintbig_sub (x, pUniverse->x0), h - 1), ULL_TO_ULLL (1)),
                                             uintbig_sl (uintbig_and (uintbig_sr (uintbig_sub (y, pUniverse->y0), h - 1), ULL_TO_ULLL (1)), 1)));
      path[h - 1].oldmc = m;
      // Creates a copy of m
      path[h - 1].newmc = calloc (1, sizeof (*path[h - 1].newmc));
      path[h - 1].newmc->result = QUERY;
      if (m)
      {
        *path[h - 1].newmc = *m;        // Copy of m
        m = m->quadrant[q];
      }
      path[h - 1].q = q;
      if (h < pUniverse->height)
        path[h].newmc->quadrant[path[h].q] = path[h - 1].newmc; // Point to the copy instead of the original.
    }                           // for (size_t h = pUniverse->height; h >= 1; h--)

    // Here, m is the bottom leaf
    oldleaf = m;
    if (m == leaf || leaf == QUERY)
    {
      // No need for change
      for (size_t h = pUniverse->height; h >= 1; h--)
        free (path[h - 1].newmc);
      free (path);
      return oldleaf;
    }

    path[0].newmc->quadrant[path[0].q] = leaf;
    path[0].newmc->population = macrocell_get_population (path[0].newmc, 1);
    // Phase 3: space contraction.
    // Bottom to top : aggregate identical macrocells of the same height, from level 1 to level pUniverse->height - 1
    // Level 0 contains only the static state (ON) and does not need to be aggregated.
    // Level pUniverse->height contains only the root and does not need to be aggregated either.
    for (size_t h = 1; h <= pUniverse->height; h++)
    {
      // Try to aggregate m at height h.
      // path has pUniverse->height elements (from level 1 (in path[0] to level pUniverse->height (in path[pUniverse->height - 1]) included)
      ASSERT (LIST_SIZE (pUniverse->listOfLevels) > h);
      Level *l = LNODE_VALUE (LIST_INDEX (pUniverse->listOfLevels, h));
      ASSERT (l && l->macrocells);
      SNODE (MacrocellId) * node;
      if (path[h - 1].oldmc)
      {
        ASSERT (!uintbig_is_zero (path[h - 1].oldmc->nb_instances));    // The olc macrocell should have been previously registered.
        path[h - 1].oldmc->nb_instances = uintbig_sub (path[h - 1].oldmc->nb_instances, ULL_TO_ULLL (1));
      }

      // [GOSPER] "a macrocell is never created if one having the same quadrants already exists.
      //           This applies recursively to the quadrants."
      // Register the new macrocell newmc (if necessary).
      if (h < pUniverse->height)
        path[h].newmc->quadrant[path[h].q] = macrocell_patternify (path[h - 1].newmc, l->macrocells);
      else
        pUniverse->root = macrocell_patternify (path[h - 1].newmc, l->macrocells);
      // Update the population of the new cell
      if (h < pUniverse->height)
        path[h].newmc->population = macrocell_get_population (path[h].newmc, 1);

      // if oldmc is not in the universe anymore, it can be forgotten.
      if (path[h - 1].oldmc && uintbig_is_zero (path[h - 1].oldmc->nb_instances))
      {
        node = SET_FIND (l->macrocells, path[h - 1].oldmc);
        ASSERT (node);
        ASSERT (SNODE_KEY (node));
        free (*SNODE_KEY (node));
        SET_REMOVE (l->macrocells, node);
      }
    }                           // for (size_t h = 1; h <= pUniverse->height; h++)
    free (path);
  }                             // if (pUniverse->height || pUniverse->x0 != x || pUniverse->y0 != y)

  if (!pUniverse->root)
    universe_reinitialize (pUniverse);
  return oldleaf;
}

void
universe_cell_set (Universe *pUniverse, intbig_t x, intbig_t y)
{
  universe_cell_accessor (pUniverse, x, y, ON);
}

void
universe_cell_unset (Universe *pUniverse, intbig_t x, intbig_t y)
{
  universe_cell_accessor (pUniverse, x, y, 0);
}

int
universe_cell_is_set (Universe *pUniverse, intbig_t x, intbig_t y)
{
  return universe_cell_accessor (pUniverse, x, y, QUERY) ? 1 : 0;
}

#define IFNOTEXIT(cond, ...) \
do { \
  if (!(cond)) \
  { \
    fprintf (stderr, "" __VA_ARGS__); \
    fprintf (stderr, "\n"); \
    exit (EXIT_FAILURE); \
  } \
} while(0)

uintbig_t
universe_RLE_readfile (Universe *pUniverse, FILE *f, intbig_t x, intbig_t y, int header)
{
  universe_reinitialize (pUniverse);

  if (header)
  {
    char *line = 0;
    size_t length;
    do
    {
      if (getline (&line, &length, f) < 0)
      {
        free (line);
        line = 0;
      }
    }
    while (line && *line == '#');

    IFNOTEXIT (line, "Missing header line");

    regex_t regvar = { 0 };
    IFNOTEXIT (regcomp (&regvar, " *([[:alnum:]]+) *= *([^ ,]+) *,?", REG_EXTENDED | REG_ICASE) == 0, "Invalid ERE. Comma separated parameters of the form 'var=value' expected.");
    regmatch_t match[3];
    for (size_t offset = 0; regexec (&regvar, line + offset, sizeof (match) / sizeof (*match), match, REG_NOTBOL | REG_NOTEOL) == 0; offset += (size_t) match[0].rm_eo)
    {
      if (!strncmp ("rule", line + offset + match[1].rm_so, (size_t) (match[1].rm_eo - match[1].rm_so)))
      {
        regex_t regrule = { 0 };
        IFNOTEXIT (regcomp (&regrule, "B([[:digit:]]+)/S([[:digit:]]+)", REG_EXTENDED | REG_ICASE) == 0, "Invalid ERE");
        regmatch_t matchBS[3];
        IFNOTEXIT (regexec (&regrule, line + offset + match[2].rm_so,
                            sizeof (matchBS) / sizeof (*matchBS), matchBS, REG_NOTBOL | REG_NOTEOL) == 0, "Invalid format for 'rule'. Format 'rule=Bnnn/Snnn' expected.");

        pUniverse->B = pUniverse->S = 0;
        for (const char *c = line + offset + match[2].rm_so + matchBS[1].rm_so; c < line + offset + match[2].rm_so + matchBS[1].rm_eo; c++)
        {
          IFNOTEXIT (isdigit (*c), "Invalid number '%c' for rule B", *c);
          pUniverse->B |= (unsigned char) (1 << (*c - '0'));
        }
        for (const char *c = line + offset + match[2].rm_so + matchBS[2].rm_so; c < line + offset + match[2].rm_so + matchBS[2].rm_eo; c++)
        {
          IFNOTEXIT (isdigit (*c), "Invalid number '%c' for rule S", *c);
          pUniverse->S |= (unsigned char) (1 << (*c - '0'));
        }
        regfree (&regrule);
      }
    }
    free (line);
    regfree (&regvar);
  }

  // Take B and S into account.
  for (size_t i = 0; i < sizeof (pUniverse->RESULT4x4) / sizeof (*(pUniverse->RESULT4x4)); i++)
    pUniverse->RESULT4x4[i] = next2x2 ((uint16_t) i, pUniverse->S, pUniverse->B);

  long unsigned int counter = 1;
  for (int c = 0; (c = fgetc (f)) && c != '!' && c != EOF;)
    switch (c)
    {
      case 'O':                // alive cell
      case 'X':                // alive cell
      case 'o':                // alive cell
      case 'x':                // alive cell
        for (; counter >= 1; counter--)
        {
          universe_cell_set (pUniverse, x, y);
          ASSERT (intbig_cmp (x, INTBIG_MAX) < 0);
          x = intbig_add (x, LL_TO_LLL (1));
        }
        counter = 1;
        break;
      case '.':                // dead cell
      case 'b':                // dead cell
      case 'B':                // dead cell
        for (; counter >= 1; counter--)
        {
          ASSERT (intbig_cmp (x, INTBIG_MAX) < 0);
          x = intbig_add (x, LL_TO_LLL (1));
        }
        counter = 1;
        break;
      case '$':                // end of line
        for (; counter >= 1; counter--)
        {
          ASSERT (intbig_cmp (y, INTBIG_MIN) > 0);
          y = intbig_sub (y, LL_TO_LLL (1));
        }
        x = INTBIG_ZERO;
        counter = 1;
        break;
      default:
        if (isdigit (c))
        {
          ungetc (c, f);
          errno = 0;
          IFNOTEXIT (fscanf (f, "%lu", &counter) == 1, "Invalid character '%c'", c);
        }
        else
          IFNOTEXIT (isblank (c) || iscntrl (c), "Invalid character '%c'", c);
        break;
    }

  return macrocell_get_population (pUniverse->root, 1);
}

static int
print_cell (SNODE (XYPos) *n, void *arg)
{
  XYPos pos = *SNODE_KEY (n);
  Explorer *explorer = arg;
  if (explorer->extractor.foreach)
    explorer->extractor.foreach (explorer->universe, explorer->spacetime, pos.x, pos.y, explorer->extractor.context);
  return EXIT_SUCCESS;
}

// [GOSPER] "To SHOW the intersection of such a slab with the spacetime, teach the macro-cell classes to check
//           whether their future cone intersects the slab, and if so, propagate the SHOW message, along with
//           appropriate x, y, and time offsets, to the quadrants and RESULTs, which are computed as necessary.
//           Then teach the 1 by 1s to signal the querying window if they get a SHOW message with time and
//           space coordinates all 0."
// [GOSPER] "SHOW is only logarithmic in the time coordinate."
static MacrocellId
universe_show_RESULT (Universe *pUniverse, MacrocellId m, SpaceTimeRegion offset, Explorer *pE, SET (XYPos) *found_cells, SET (SpaceTimeRegion) *already_explored)
{
  if (m == 0 || offset.height < 2)
    return 0;

  int empty;
  empty = 1;
  for (Quadrant q = 0; q < NB_QUADRANTS && empty; q++)
    if (m->quadrant[q])
      empty = 0;
  if (empty)
    return 0;

  // [GOSPER] "(...) most future cones (and, by the geometry, the cones of all their components, recursively,)
  //           will not intersect non-gigantic windows."
  // No need to explore deeper in this macrocell.
  if (!time_overlap (offset.height, offset.tbase, pE->spacetime.time.instant) || !space_overlap (offset.height, offset.xmin, offset.ymin, pE->spacetime.space.window))
    return m->result;

  // Exploring m comes to exploring 13 overlapping quadrants (recursively).
  // Therefore, we must avoid to explore twice the same spacetime region of the universe (this spares a LOT of CPU time.)
  if (SET_FIND (already_explored, offset))
    return m->result;
  else
    SET_ADD (already_explored, offset);

  // m should be registered already.
  ASSERT (LIST_SIZE (pUniverse->listOfLevels) > offset.height &&
          LNODE_VALUE (LIST_INDEX (pUniverse->listOfLevels, offset.height)) &&
          LNODE_VALUE (LIST_INDEX (pUniverse->listOfLevels, offset.height))->macrocells
          && m == macrocell_fetch_pattern (m, LNODE_VALUE (LIST_INDEX (pUniverse->listOfLevels, offset.height))->macrocells));
  ASSERT (m->result != QUERY);  // The RESULT should have been computed already.

  uintbig_t deltat = uintbig_sub (pE->spacetime.time.instant, offset.tbase);
  // The instant of m corresponds to the request instant of pE
  if (uintbig_is_zero (deltat))
  {
    // For debugging purpose only.
    show_overlapping (offset.height, offset.tbase, pE->spacetime.time.instant, offset.xmin, offset.ymin, pE->spacetime.space.window);
    macrocell_get_cells_in_window (m, offset.height, offset.xmin, offset.ymin, pE->spacetime.space.window, found_cells);
    return m->result;
  }
  // The instant of m->hresult corresponds to the request instant of pE
  uintbig_t quarter_size = uintbig_sl (ULL_TO_ULLL (1), offset.height - 2);
  if (uintbig_cmp (deltat, quarter_size) == 0)
  {
    if (m->result)
    {
      // For debugging purpose only.
      show_overlapping (offset.height - 1, pE->spacetime.time.instant, pE->spacetime.time.instant,
                        uintbig_add (offset.xmin, quarter_size), uintbig_add (offset.ymin, quarter_size), pE->spacetime.space.window);
      macrocell_get_cells_in_window (m->result, offset.height - 1,
                                     uintbig_add (offset.xmin, quarter_size), uintbig_add (offset.ymin, quarter_size), pE->spacetime.space.window, found_cells);
    }
    return m->result;
  }

  // For debugging purpose only.
  show_overlapping (offset.height, offset.tbase, pE->spacetime.time.instant, offset.xmin, offset.ymin, pE->spacetime.space.window);

  // Explore RESULT, S/4 generations ahead of m, where S is the size of m.
  /* Macrocell m (.), of size S  ->  RESULT (H), of size S/2, S/4 generations ahead of m.

     ........                         ........
     ........                         ........
     ........                         ..HHHH..
     ........                    ->   ..HHHH..
     ........                         ..HHHH..
     ........                         ..HHHH..
     ........                         ........
     ........                         ........
   */

  // HERE, 0 < deltat < quarter_size

  Level *l1 = LNODE_VALUE (LIST_INDEX (pUniverse->listOfLevels, offset.height - 1));
  ASSERT (l1 && l1->macrocells);

  /* Macrocell m, of size S:
     ........                            ........
     ........                            .cccccc.
     ........                 cell[36]   .cccccc.
     ........                 36 cells   .cccccc.
     ........                            .cccccc.
     ........                            .cccccc.
     ........                            .cccccc.
     ........                            ........
   */
  MacrocellId cell[36] = { 0 }; // height - 3

  /* 1 to 4 */
  /* Get RESULT (h) of quadrants (q) of macrocell m.
     q is half the size (S/2) of m, therefore h is S/2/4 = S/8 generations ahead of m.
     (NW)          (NE)          (SW)          (SE)
     qqqq....      ....qqqq      ........      ........
     qhhq....      ....qhhq      ........      ........
     qhhq....      ....qhhq      ........      ........
     qqqq....      ....qqqq      ........      ........
     ........      ........      qqqq....      ....qqqq
     ........      ........      qhhq....      ....qhhq
     ........      ........      qhhq....      ....qhhq
     ........      ........      qqqq....      ....qqqq
   */
    /* *INDENT-OFF* */
    struct
    {
      MacrocellId *cell[NB_QUADRANTS];    // cells h
    } unit_1_4[NB_QUADRANTS] =
    {
      {{&cell[ 0], &cell[ 1], &cell[ 6], &cell[ 7]}},  // NW
      {{&cell[ 4], &cell[ 5], &cell[10], &cell[11]}},  // NE
      {{&cell[24], &cell[25], &cell[30], &cell[31]}},  // SW
      {{&cell[28], &cell[29], &cell[34], &cell[35]}}   // SE
    };
    /* *INDENT-ON* */

  for (Quadrant u = 0; u < NB_QUADRANTS; u++)
  {
    SpaceTimeRegion r2 = offset;
    r2.height = offset.height - 1;
    if (u == NE || u == SE)
      r2.xmin = uintbig_add (uintbig_add (offset.xmin, quarter_size), quarter_size);
    if (u == SE || u == SW)
      r2.ymin = uintbig_add (uintbig_add (offset.ymin, quarter_size), quarter_size);
    if (universe_show_RESULT (pUniverse, m->quadrant[u], r2, pE, found_cells, already_explored))        // height - 2
      for (Quadrant i = 0; i < NB_QUADRANTS; i++)
        *unit_1_4[u].cell[i] = m->quadrant[u]->result->quadrant[i];
  }

  /* 5 to 9 */
  /* Get RESULT (h) of quadrants (q) of macrocell m
     q is half the size (S/2) of m, therefore h is S/2/4 = S/8 generations ahead of m.
     (5)           (6)           (7)           (8)           (9)
     ..qqqq..      ........      ........      ........      ........
     ..qhhq..      ........      ........      ........      ........
     ..qhhq..      ........      qqqq....      ....qqqq      ..qqqq..
     ..qqqq..      ........      qhhq....      ....qhhq      ..qhhq..
     ........      ..qqqq..      qhhq....      ....qhhq      ..qhhq..
     ........      ..qhhq..      qqqq....      ....qqqq      ..qqqq..
     ........      ..qhhq..      ........      ........      ........
     ........      ..qqqq..      ........      ........      ........
   */
  MacrocellId mtemp_5_9[5];     // height - 1
  for (unsigned int u = 0; u < 5; u++)
    mtemp_5_9[u] = calloc (1, sizeof (*mtemp_5_9[u]));

  /* *INDENT-OFF* */
  struct
  {
    MacrocellId *cell[NB_QUADRANTS];
  } unit_5_9[5] =
  {
    {{&cell[ 2], &cell[ 3], &cell[ 8], &cell[ 9]}},  // 5
    {{&cell[26], &cell[27], &cell[32], &cell[33]}},  // 6
    {{&cell[12], &cell[13], &cell[18], &cell[19]}},  // 7
    {{&cell[16], &cell[17], &cell[22], &cell[23]}},  // 8
    {{&cell[14], &cell[15], &cell[20], &cell[21]}}   // 9
  };
  /* *INDENT-ON* */
  /* 5 */
  if (m->quadrant[NW])
  {
    mtemp_5_9[0]->quadrant[NW] = m->quadrant[NW]->quadrant[NE];
    mtemp_5_9[0]->quadrant[SW] = m->quadrant[NW]->quadrant[SE];
  }
  if (m->quadrant[NE])
  {
    mtemp_5_9[0]->quadrant[NE] = m->quadrant[NE]->quadrant[NW];
    mtemp_5_9[0]->quadrant[SE] = m->quadrant[NE]->quadrant[SW];
  }
  /* 6 */
  if (m->quadrant[SW])
  {
    mtemp_5_9[1]->quadrant[NW] = m->quadrant[SW]->quadrant[NE];
    mtemp_5_9[1]->quadrant[SW] = m->quadrant[SW]->quadrant[SE];
  }
  if (m->quadrant[SE])
  {
    mtemp_5_9[1]->quadrant[NE] = m->quadrant[SE]->quadrant[NW];
    mtemp_5_9[1]->quadrant[SE] = m->quadrant[SE]->quadrant[SW];
  }
  /* 7 */
  if (m->quadrant[NW])
  {
    mtemp_5_9[2]->quadrant[NW] = m->quadrant[NW]->quadrant[SW];
    mtemp_5_9[2]->quadrant[NE] = m->quadrant[NW]->quadrant[SE];
  }
  if (m->quadrant[SW])
  {
    mtemp_5_9[2]->quadrant[SW] = m->quadrant[SW]->quadrant[NW];
    mtemp_5_9[2]->quadrant[SE] = m->quadrant[SW]->quadrant[NE];
  }
  /* 8 */
  if (m->quadrant[NE])
  {
    mtemp_5_9[3]->quadrant[NW] = m->quadrant[NE]->quadrant[SW];
    mtemp_5_9[3]->quadrant[NE] = m->quadrant[NE]->quadrant[SE];
  }
  if (m->quadrant[SE])
  {
    mtemp_5_9[3]->quadrant[SW] = m->quadrant[SE]->quadrant[NW];
    mtemp_5_9[3]->quadrant[SE] = m->quadrant[SE]->quadrant[NE];
  }
  /* 9 */
  if (m->quadrant[NW])
    mtemp_5_9[4]->quadrant[NW] = m->quadrant[NW]->quadrant[SE];
  if (m->quadrant[NE])
    mtemp_5_9[4]->quadrant[NE] = m->quadrant[NE]->quadrant[SW];
  if (m->quadrant[SW])
    mtemp_5_9[4]->quadrant[SW] = m->quadrant[SW]->quadrant[NE];
  if (m->quadrant[SE])
    mtemp_5_9[4]->quadrant[SE] = m->quadrant[SE]->quadrant[NW];

  /* 5 to 9 */
  for (unsigned int u = 0; u < 5; u++)
  {
    mtemp_5_9[u] = macrocell_fetch_pattern (mtemp_5_9[u], l1->macrocells);

    SpaceTimeRegion r2 = offset;
    r2.height = offset.height - 1;
    switch (u)
    {
      case 0:
        r2.xmin = uintbig_add (offset.xmin, quarter_size);
        break;
      case 1:
        r2.xmin = uintbig_add (offset.xmin, quarter_size);
        r2.ymin = uintbig_add (uintbig_add (offset.ymin, quarter_size), quarter_size);
        break;
      case 2:
        r2.ymin = uintbig_add (offset.ymin, quarter_size);
        break;
      case 3:
        r2.xmin = uintbig_add (uintbig_add (offset.xmin, quarter_size), quarter_size);
        r2.ymin = uintbig_add (offset.ymin, quarter_size);
        break;
      case 4:
        r2.xmin = uintbig_add (offset.xmin, quarter_size);
        r2.ymin = uintbig_add (offset.ymin, quarter_size);
        break;
      default:
    }

    if (universe_show_RESULT (pUniverse, mtemp_5_9[u], r2, pE, found_cells, already_explored))  // height - 2
      for (Quadrant i = 0; i < NB_QUADRANTS; i++)
        *unit_5_9[u].cell[i] = mtemp_5_9[u]->result->quadrant[i];
  }

  /* Therefore, we have 36 cells h, S/8 generations ahead of m.
     ........
     .hhhhhh.
     .hhhhhh.
     .hhhhhh.
     .hhhhhh.
     .hhhhhh.
     .hhhhhh.
     ........
   */

  /* 10 to 13 */
  /* Get RESULT (H) of quadrants (x) of cells h of m (.)
     x is half the size (S/2) of m, therefore H is S/2/4 = S/8 generations ahead of x, part of h, itself S/8 generations ahead of m.
     Therefore, H is S/8 + S/8 = S/4 generations ahead of m.
     (NW)          (NE)          (SW)          (SE)
     ........      ........      ........      ........
     .xxxxhh.      .hhxxxx.      .hhhhhh.      .hhhhhh.
     .xHHxhh.      .hhxHHx.      .hhhhhh.      .hhhhhh.
     .xHHxhh.      .hhxHHx.      .xxxxhh.      .hhxxxx.
     .xxxxhh.      .hhxxxx.      .xHHxhh.      .hhxHHx.
     .hhhhhh.      .hhhhhh.      .xHHxhh.      .hhxHHx.
     .hhhhhh.      .hhhhhh.      .xxxxhh.      .hhxxxx.
     ........      ........      ........      ........
   */

    /* *INDENT-OFF* */
    struct
    {
      MacrocellId cell[NB_QUADRANTS /*j*/][NB_QUADRANTS /*k*/];
    } unit_10_13[NB_QUADRANTS] =
    {
      { // NW (j)
        {{cell[ 0], cell[ 1], cell[ 6], cell[ 7]},  // NW (k)
         {cell[ 2], cell[ 3], cell[ 8], cell[ 9]},  // NE
         {cell[12], cell[13], cell[18], cell[19]},  // SW
         {cell[14], cell[15], cell[20], cell[21]}}  // SE
      },
      { // NE
        {{cell[ 2], cell[ 3], cell[ 8], cell[ 9]},  // NW
         {cell[ 4], cell[ 5], cell[10], cell[11]},  // NE
         {cell[14], cell[15], cell[20], cell[21]},  // SW
         {cell[16], cell[17], cell[22], cell[23]}}  // SE
      },
      { // SW
        {{cell[12], cell[13], cell[18], cell[19]},  // NW
         {cell[14], cell[15], cell[20], cell[21]},  // NE
         {cell[24], cell[25], cell[30], cell[31]},  // SW
         {cell[26], cell[27], cell[32], cell[33]}}  // SE
      },
      { // SE
        {{cell[14], cell[15], cell[20], cell[21]},  // NW
         {cell[16], cell[17], cell[22], cell[23]},  // NE
         {cell[26], cell[27], cell[32], cell[33]},  // SW
         {cell[28], cell[29], cell[34], cell[35]}}  // SE
      }
    };
    /* *INDENT-ON* */

  Level *l2 = LNODE_VALUE (LIST_INDEX (pUniverse->listOfLevels, offset.height - 2));
  ASSERT (l2 && l2->macrocells);

  uintbig_t eighth_size = uintbig_sr (quarter_size, 1);
  for (Quadrant u = 0; u < NB_QUADRANTS; u++)
  {
    MacrocellId mtemp = calloc (1, sizeof (*mtemp));
    MacrocellId qtemp[NB_QUADRANTS];
    for (Quadrant j = 0; j < NB_QUADRANTS; j++)
    {
      qtemp[j] = calloc (1, sizeof (*qtemp[j]));
      for (Quadrant k = 0; k < NB_QUADRANTS; k++)
        qtemp[j]->quadrant[k] = unit_10_13[u].cell[j][k];
      mtemp->quadrant[j] = macrocell_fetch_pattern (qtemp[j], l2->macrocells);
    }
    mtemp = macrocell_fetch_pattern (mtemp, l1->macrocells);
    SpaceTimeRegion r2 = offset;
    r2.tbase = uintbig_add (offset.tbase, eighth_size);
    r2.xmin = uintbig_add (offset.xmin, eighth_size);
    r2.ymin = uintbig_add (offset.ymin, eighth_size);
    r2.height = offset.height - 1;
    switch (u)
    {
      case NW:
        break;
      case NE:
        r2.xmin = uintbig_add (r2.xmin, quarter_size);
        break;
      case SW:
        r2.ymin = uintbig_add (r2.ymin, quarter_size);
        break;
      case SE:
        r2.xmin = uintbig_add (r2.xmin, quarter_size);
        r2.ymin = uintbig_add (r2.ymin, quarter_size);
        break;
      default:
    }
    universe_show_RESULT (pUniverse, mtemp, r2, pE, found_cells, already_explored);     // height - 2
  }                             // for (Quadrant u = 0; u < NB_QUADRANTS; u++)

  /* Done: hresult H, half the size of m, is S/4 generations ahead of m:
     ........
     ........
     ..HHHH..
     ..HHHH..
     ..HHHH..
     ..HHHH..
     ........
     ........
   */

  return m->result;
}

// [GOSPER] "The following algorithm will let us freely explore the future space-times of large initial configurations,
//           provided that they are sufficiently repetitious, both structurally and behaviorally"
// [GOSPER] "let the observation drive the computing"
// [GOSPER] "one wishes to have one or more "windows" into the spacetime, rectangular slabs of cells, one time unit thick."
// Returns the number of cells in the explorer.space.window at instant explorer.time.instant.
uintbig_t
universe_explore (Universe *pUniverse, Explorer explorer)
{
  if (intbig_cmp (explorer.spacetime.space.window.NWvertex.x, explorer.spacetime.space.window.SEvertex.x) >= 0)
  {
    explorer.spacetime.space.window.NWvertex.x = INTBIG_MIN;
    explorer.spacetime.space.window.SEvertex.x = INTBIG_MAX;
  }
  if (intbig_cmp (explorer.spacetime.space.window.NWvertex.y, explorer.spacetime.space.window.SEvertex.y) >= 0)
  {
    explorer.spacetime.space.window.NWvertex.y = INTBIG_MIN;
    explorer.spacetime.space.window.SEvertex.y = INTBIG_MAX;
  }

  explorer.universe = pUniverse;

  if (explorer.extractor.preaction)
    explorer.extractor.preaction (pUniverse, explorer.spacetime, explorer.extractor.context);

  SET (XYPos) * found_cells = SET_CREATE (XYPos, xypos_lt);
  SET (SpaceTimeRegion) * already_explored = SET_CREATE (SpaceTimeRegion);

  uintbig_t population;
  if (pUniverse->root == 0)
  {                             /* do nothing */
  }
  else if (uintbig_is_zero (explorer.spacetime.time.instant))
    macrocell_get_cells_in_window (pUniverse->root, pUniverse->height, pUniverse->x0, pUniverse->y0, explorer.spacetime.space.window, found_cells);
  else
  {
    unsigned int min_height = 2;
    for (uintbig_t t = uintbig_sub (explorer.spacetime.time.instant, ULL_TO_ULLL (1)); !uintbig_is_zero (t); t = uintbig_sr (t, 1))
      min_height++;

    // [GOSPER] "outermost SHOW method ensures that the configuration being probed is surroundedby enough vacuum so that its
    //           future cone entirely contains the probe window, no matter how large or remote in space or time. Thus,
    //           there are never any edge effects."
    // Make sure the universe m is high enough and surrounded by sufficient empty space so that its horizon is not reachable at light speed.
    while (pUniverse->height < min_height || universe_is_closed (pUniverse))
      universe_expand (pUniverse);

    uintbig_t quarter_size = uintbig_sl (ULL_TO_ULLL (1), pUniverse->height - 2);
    ASSERT (uintbig_cmp (explorer.spacetime.time.instant, quarter_size) <= 0);

    // Compute the forecast of the whole universe and its reachable neighborhood (ie within reachable distance at speed of light).
    // Compute 4 RESULTs, S/4 generations ahead of m, where S is the size of m.
    // Macrocell m (. and x), of size S  ->  RESULTs, of size S/2, S/4 generations ahead of m.

    // Step 1 : Create 4 shifted universes, uSE, uSW, uNW, uNE, of size S, respectively centered around each quadrant of m.
    // 'x' = might be populated, '.' = empty space of m, '-' = empty space of the shifted universe.
    //   Universe                          uSE               uSW               uNW               uNE
    //      m                           --------          --------
    //                                  --------          --------
    //   ........                       --......          ......--
    //   ........                       --......          ......--
    //   ..xxxx..                       --..xxxx          xxxx..--          xxxx..--          --..xxxx
    //   ..xxxx..                    -> --..xxxx          xxxx..--          xxxx..--          --..xxxx
    //   ..xxxx..                       --..xxxx          xxxx..--          xxxx..--          --..xxxx
    //   ..xxxx..                       --..xxxx          xxxx..--          xxxx..--          --..xxxx
    //   ........                                                           ......--          --......
    //   ........                                                           ......--          --......
    //                                                                      --------          --------
    //                                                                      --------          --------
    // uSE.se.nw = m.nw.se    uSW.sw.... =        uNW.nw.... =       uNE.ne.... = ...
    // uSE.se.ne = m.ne.sw
    // uSE.se.sw = m.sw.ne
    // uSE.se.se = m.se.nw
    // Step 2 : Compute the four RESULTs for each of theses universes.

    Level *pl0 = LNODE_VALUE (LIST_INDEX (pUniverse->listOfLevels, pUniverse->height));
    ASSERT (pl0 && pl0->macrocells);
    Level *pl1 = LNODE_VALUE (LIST_INDEX (pUniverse->listOfLevels, pUniverse->height - 1));
    ASSERT (pl1 && pl1->macrocells);
    MacrocellId shifted[NB_QUADRANTS];
    for (Quadrant u = 0; u < NB_QUADRANTS; u++)
    {
      shifted[u] = calloc (1, sizeof (*shifted[u]));
      shifted[u]->quadrant[u] = calloc (1, sizeof (*shifted[u]->quadrant[u]));
      for (Quadrant q = 0; q < NB_QUADRANTS; q++)
        if (pUniverse->root->quadrant[q])
          shifted[u]->quadrant[u]->quadrant[q] = pUniverse->root->quadrant[q]->quadrant[NB_QUADRANTS - q - 1];
      shifted[u]->quadrant[u] = macrocell_patternify (shifted[u]->quadrant[u], pl1->macrocells);
      shifted[u] = macrocell_patternify (shifted[u], pl0->macrocells);
      PRINT ("Compute quadrant %i...\n", u + 1);
      universe_get_RESULT (pUniverse, shifted[u], pUniverse->height);

      uintbig_t xmin, ymin;
      switch (u)
      {
        case NW:
          xmin = uintbig_add (pUniverse->x0, quarter_size);
          ymin = uintbig_add (pUniverse->y0, quarter_size);
          break;
        case SW:
          xmin = uintbig_add (pUniverse->x0, quarter_size);
          ymin = uintbig_sub (pUniverse->y0, quarter_size);
          break;
        case NE:
          xmin = uintbig_sub (pUniverse->x0, quarter_size);
          ymin = uintbig_add (pUniverse->y0, quarter_size);
          break;
        case SE:
          xmin = uintbig_sub (pUniverse->x0, quarter_size);
          ymin = uintbig_sub (pUniverse->y0, quarter_size);
          break;
        default:
      }
      PRINT ("Explore quadrant %i...\n", u + 1);
      SpaceTimeRegion r = {.xmin = xmin,.ymin = ymin,.tbase = UINTBIG_ZERO,.height = pUniverse->height
      };
      universe_show_RESULT (pUniverse, shifted[u], r, &explorer, found_cells, already_explored);
    }
  }

  population = ULL_TO_ULLL (SET_SIZE (found_cells));
  SET_TRAVERSE (found_cells, print_cell, &explorer);
  SET_DESTROY (found_cells);
  SET_DESTROY (already_explored);

  if (explorer.extractor.postaction)
    explorer.extractor.postaction (pUniverse, explorer.spacetime, population, explorer.extractor.context);

  return population;
}
