#include <stdio.h>
#include <locale.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include "hgolbi.h"

#define IFNOTEXIT(cond, ...) \
do { \
  if (!(cond)) \
  { \
    fprintf (stderr, "" __VA_ARGS__); \
    fprintf (stderr, "\n"); \
    exit (EXIT_FAILURE); \
  } \
} while(0)

static void
preaction (Universe * universe, SpaceTime st, void * arg)
{
  (void) arg;
  (void) universe;
  printf ("Cells in universe within window [%1$+'V ; %2$+'V] x [%3$+'V ; %4$+'V] at generation %5$'U:\n",
          st.space.window.NWvertex.x, st.space.window.SEvertex.x,
          st.space.window.NWvertex.y, st.space.window.SEvertex.y,
          st.time.instant);
}

static void
extractor (Universe * universe, SpaceTime st, intbig_t x, intbig_t y, void *arg)
{
  (void) arg;
  (void) universe;
  printf ("- cell at position (%+'12V, %+'12V) at time %'12U\n", x, y, st.time.instant);
}

static void
postaction (Universe * universe, SpaceTime st, uintbig_t numcells, void *arg)
{
  (void) arg;
  (void) universe;
  printf ("There are %1$'U cells in universe within window [%2$+'V ; %3$+'V] x [%4$+'V ; %5$+'V] at generation %6$'U.\n",
          numcells, st.space.window.NWvertex.x, st.space.window.SEvertex.x,
          st.space.window.NWvertex.y, st.space.window.SEvertex.y,
          st.time.instant);
}

#define EXPLORE(mytime)\
  do {\
    e.spacetime.time.instant = (mytime);\
    universe_explore (pUniverse, e);\
  } while(0)

static void
TU (Universe *pUniverse, Explorer e)
{
#ifndef DEBUG
  // Glider
  universe_cell_set (pUniverse, LL_TO_LLL (0), LL_TO_LLL (0));
  universe_cell_set (pUniverse, LL_TO_LLL (1), LL_TO_LLL (0));
  universe_cell_set (pUniverse, LL_TO_LLL (2), LL_TO_LLL (0));
  universe_cell_set (pUniverse, LL_TO_LLL (2), LL_TO_LLL (1));
  universe_cell_set (pUniverse, LL_TO_LLL (1), LL_TO_LLL (2));
  // Block
  universe_cell_set (pUniverse, LL_TO_LLL (10), LL_TO_LLL (10));
  universe_cell_set (pUniverse, LL_TO_LLL (10), LL_TO_LLL (11));
  universe_cell_set (pUniverse, LL_TO_LLL (11), LL_TO_LLL (10));
  universe_cell_set (pUniverse, LL_TO_LLL (11), LL_TO_LLL (11));

  for (unsigned long int i = 0; i < 16; i += 1)
    EXPLORE (ULL_TO_ULLL (i));
  for (unsigned long int i = 16; i; i -= 1)
    EXPLORE (ULL_TO_ULLL (i));
  EXPLORE (ULL_TO_ULLL (10000));
  EXPLORE (ULL_TO_ULLL (20000));
  EXPLORE (ULL_TO_ULLL (20000000));
  EXPLORE (ULL_TO_ULLL (ULONG_MAX));
  universe_reinitialize (pUniverse);
#endif
}

int
main (int argc, char *const argv[])
{
  setlocale (LC_ALL, "");

  Explorer e = { 0 };
  e.extractor.preaction = preaction;
  e.extractor.foreach = extractor;
  e.extractor.postaction = postaction;
  Universe *pUniverse = universe_create ();
  IFNOTEXIT (pUniverse, "Memory allocation error.");

  // Parsing command line, e.g.: ./hgolbi_example -t1/2 -x-9/10,3/4 -y-5/6,7/8 </dev/null
  uintbig_t t = UINTBIG_ZERO;
  intbig_t xmin, ymin, xmax, ymax;
  xmin = ymin = xmax = ymax = INTBIG_ZERO;
  char *ptr = 0;
  char *endptr = 0;
  char *ntokptr = 0;
  char *ctokptr = 0;
  intbig_t *v = 0;
  int opt, sign;
  const char *num_sep = "_";
  const char *coord_sep = ",";
  const char *optstring = ":Ut:x:y:";
  while ((opt = getopt (argc, argv, optstring)) != -1)
    switch (opt)
    {
      case 'U':
        TU (pUniverse, e);
        break;
      case 't':
        for (t = UINTBIG_ZERO, ptr = optarg; (ptr = strtok_r (ptr, num_sep, &ntokptr)); ptr = 0)
        {
          errno = 0;
          t = uintbig_add (uintbig_sl (t, ULL_NB_BITS), ULL_TO_ULLL (strtoull (ptr, &endptr, 10)));
          IFNOTEXIT (optarg && *optarg && *endptr == 0, "Invalid number %s for option '-%c'", ptr, opt);
        }
        break;
      case 'x':
      case 'y':
        for (v = (opt == 'x' ? &xmin : &ymin), ptr = optarg; (ptr = strtok_r (ptr, coord_sep, &ctokptr));
             v = (opt == 'x' ? &xmax : &ymax), ptr = 0)
        {
          for (sign = 0, *v = INTBIG_ZERO; (ptr = strtok_r (ptr, num_sep, &ntokptr)); ptr = 0)
          {
            errno = 0;
            if (sign == 0)
            {
              if (intbig_is_negative (*v = LL_TO_LLL (strtoll (ptr, &endptr, 10))))
              {
                sign = -1;
                *v = intbig_opposite (*v);
              }
              else
                sign = 1;
            }
            else
              *v = uintbig_add (uintbig_sl (*v, ULL_NB_BITS), ULL_TO_ULLL (strtoull (ptr, &endptr, 10)));
            IFNOTEXIT (ptr && *ptr && *endptr == 0, "Invalid number %s for option '-%c'", ptr, opt);
          }
          if (sign == -1)
            *v = intbig_opposite (*v);
        }
        break;
    }
  e.spacetime.space.window.NWvertex.x = xmin;
  e.spacetime.space.window.NWvertex.y = ymin;
  e.spacetime.space.window.SEvertex.x = xmax;
  e.spacetime.space.window.SEvertex.y = ymax;

  FILE *f = 0;
  if (argc > optind)
    IFNOTEXIT (f = fopen (argv[optind], "r"), "Can not read file '%s'", argv[optind]);
  else
  {
    f = stdin;

    int c;
    if ((c = fgetc (f)) != EOF)
      ungetc (c, f);
    else
    {
#define RULE "rule=B3/S23\n"
      //const char *pattern = RULE "9bo12b$7bobo12b$6bobo13b$2o3bo2bo11b2o$2o4bobo11b2o$7bobo12b$9bo!";  // Queen bee shuttle, period 30:
      //const char *pattern = RULE "24bo11b$22bobo11b$12b2o6b2o12b2o$11bo3bo4b2o12b2o$2o8bo5bo3b2o14b$2o8bo3bob2o4bobo11b$10bo5bo7bo11b$11bo3bo20b$12b2o!";      // Gosper glider gun
      //const char *pattern = RULE ".XX$XX$.X";  // R-pentomino, stabilizes at generation 1103 with 116 cells, including one escaped glider at generation 69:
      //const char *pattern = RULE "10X";        // Pentadecathlon (period 15):
      //const char *pattern = RULE "3X"; // Blinker:
      //const char *pattern = RULE "xx$xx20$.o.$..o$ooo"; // Block + Glider:
      const char *pattern = RULE "bo5b$3bo3b$2o2b3o!";  // Acorn, takes 5206 generations to stabilize to 633 cells, including 13 escaped gliders:
      //const char *pattern = RULE "ooo$.o.";    // Tee or Tetromino, stabilizes to 12 cells in a 9x9 square at 10th generation.
      //const char *pattern = RULE "......o$oo$.o...ooo"; // Die-hard, eventually disappears after 130 generations
      //const char *pattern = RULE "......X$....X.XX$....X.X$....X$..X$X.X";     // Infinite growth, block-laying switch engine that leaves behind two-by-two still life blocks as its translates itself across the game's universe.
      //const char *pattern = RULE "77bo$77bo$77bo21$3o20$3bo$3bo$3bo5$20b3o$9b3o10bo$22bo$21bo!";        // 18-cell 40514-generation methuselah. The stable pattern that results from 40514M (excluding 70 escaping gliders) has 3731 cells and consists of 248 blinkers (including 21 traffic lights), 218 blocks, 163 beehives (including nine honey farms), 56 loaves, 39 boats, 10 ships, nine tubs, five ponds, four beacons, two toads, one barge, one eater 1 and one long boat.
      //const char *pattern = RULE "10001o!";
      //const char *pattern = RULE "15366bo$15366bo$15364boo$15363bo$15363bo$15363bo$15363bo6$15393bo$" "15392boo$15390bobbo$$15390bobo$15391bo133$15568boo$15569boo$15569bo29$" "15554bo$15553bobo$15555bo$15556bo507$59722boo$59721boo$59722bo29$" "59737bo$59736bobo$59736bo$59735bo13907$bo3bo$bbobo$obbo$o$o21$33bo$32b" "o$31bo$32bo$33bo$29b3o!";       // Metacatacryst, exhibits quadratic growth.

      IFNOTEXIT (f = fmemopen ((void *) pattern, strlen (pattern), "r"), "Can not read pattern");
    }
  }

  printf ("%'U cells have been read from the RLE pattern.\n", universe_RLE_readfile (pUniverse, f, INTBIG_ZERO, INTBIG_ZERO, 1));
  fclose (f);

  EXPLORE (UINTBIG_ZERO);
  if (!uintbig_is_zero (t))
    EXPLORE (t);

  universe_destroy (pUniverse);
  printf ("Done.\n");

  return EXIT_SUCCESS;
}
