#ifndef __HGOLBI_H__
#  define __HGOLBI_H__

#  include "bitl.h"

// Abstract data type Universe.
struct sUniverse;
typedef struct sUniverse Universe;

// An instant in time.
typedef struct
{
  uintbig_t instant;
} Time;
// A window in a 2D space rectangle, defined by its opposite vertices at North-West and South-East.
typedef struct
{
  intbig_t y, x;
} XYPos;
typedef struct
{
  XYPos NWvertex, SEvertex;
} Window;
// A window in space.
typedef struct
{
  Window window;
} Space;
typedef struct
{
  Space space;
  Time time;
} SpaceTime;
// An extractor.
typedef struct
{
  int (*preaction) (Universe *, SpaceTime st, void * context);
  int (*foreach) (Universe *, SpaceTime st, intbig_t x, intbig_t y, void * context);
  int (*postaction) (Universe *, SpaceTime st, uintbig_t num_cells, void * context);
  void *context;         // For user purpose.
} Extractor;
// An explorer of universe in space and time.
typedef struct
{
  SpaceTime spacetime;
  Extractor extractor;
  Universe *universe;    // Reserved, do not use.
} Explorer;

// Create a universe.
// Calls xintbig_printf_init () for convenience.
Universe *universe_create (void);

// Reinitialize a universe before reuse.
void universe_reinitialize (Universe * pUniverse);

// Destroy a universe.
void universe_destroy (Universe * pUniverse);

// Add a cell in universe at position (x, y).
// x and y can be initialized with LL_TO_LLL (l) where l is an long long integer.
void universe_cell_set (Universe * pUniverse, intbig_t x, intbig_t y);

// Remove a cell in universe from position (x, y).
// x and y can be initialized with LL_TO_LLL (l) where l is an long long integer.
void universe_cell_unset (Universe * pUniverse, intbig_t x, intbig_t y);

// Check if a cell is set at position (x, y) in universe.
// x and y can be initialized with LL_TO_LLL (l) where l is an long long integer.
int universe_cell_is_set (Universe * pUniverse, intbig_t x, intbig_t y);

// Initialize the universe from a RLE file f.
uintbig_t universe_RLE_readfile (Universe * pUniverse, FILE * f, intbig_t x, intbig_t y, int header);

// Explore all cells at generation explorer.time.instant in window explorer.space.window.
// The callback functions of explorer.extractor are :
// - preaction is called before exploration with context as argument.
// - foreach is called for each found cell, with its position (x and y) and context as arguments.
// - postaction is called after exploration with the number of cells and context as arguments.
// Each of these functions have pUniverse and explorer.spacetime as first arguments.
uintbig_t universe_explore (Universe * pUniverse, Explorer explorer);

#endif
