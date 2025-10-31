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
  void (*preaction) (Universe *, SpaceTime st, void *context);
  void (*foreach) (Universe *, SpaceTime st, intbig_t x, intbig_t y, void *context);
  void (*postaction) (Universe *, SpaceTime st, uintbig_t num_cells, void *context);
  void *context;                // For user purpose.
} Extractor;
// An explorer of universe in space and time.
typedef struct
{
  SpaceTime spacetime;
  Extractor extractor;
  Universe *universe;           // Reserved, do not use.
} Explorer;

// Create a universe and return the pointer to it.
// The rule of game of life B3/S23 is used by default.
// Calls xintbig_printf_init () for convenience.
Universe *universe_create (void);

// Set the rule which format must be "Bnnn/Snnn"
int universe_set_BLE_rules (Universe * pUniverse, const char *rules);

// Destroy a universe.
void universe_destroy (Universe * pUniverse);

// Add a cell in universe at position (x, y).
// x and y can be initialised with LL_TO_LLL (l) where l is an long long integer.
void universe_cell_set (Universe * pUniverse, intbig_t x, intbig_t y);

// Remove a cell in universe from position (x, y).
// x and y can be initialised with LL_TO_LLL (l) where l is an long long integer.
void universe_cell_unset (Universe * pUniverse, intbig_t x, intbig_t y);

// Check if a cell is set (returns 1) or unset (returns 0) at position (x, y) in universe.
// x and y can be initialised with LL_TO_LLL (l) where l is an long long integer.
int universe_cell_is_set (Universe * pUniverse, intbig_t x, intbig_t y);

// Initialise the universe from a RLE file f.
// Returns the number of cells initialised in the universe.
// x and y are the coordinates of the North-West corner of the universe.
// x and y can be initialised with LL_TO_LLL (l) where l is an long long integer.
// header should be 1 if the RLE contains a header line (possibly preceded by commented lines starting with #), 0 otherwise.
size_t universe_RLE_readfile (Universe * pUniverse, FILE * f, intbig_t x, intbig_t y, int header);

// Explore all cells in window explorer.space.window at generation explorer.time.instant.
// Returns the number of cells found in the window.
// The callback functions of explorer.extractor are used as following :
// - preaction is called before exploration.
// - foreach is called for each found cell, with its position (x and y) as 3rd and 4th argument.
// - postaction is called after exploration with the number of found cells as 3rd argument.
// Each of these functions have pUniverse and explorer.spacetime as first two arguments and explorer.extractor.context as last argument.
uintbig_t universe_explore (Universe * pUniverse, Explorer explorer);

#endif
