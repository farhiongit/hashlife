An implementation of the Hashlife algorithm
====================

This code, in `hgolbi.c`, implements in C language the
[Hash-life algorithm proposed by R. Gosper in 1984](https://doi.org/10.1016%2F0167-2789%2884%2990251-3)
to explore the evolution of the Conway's Game Of Life.

It offers :

- a clean public user interface defined in `hgolbi.h`
- a reader of universe confifurations in RLE format

It can handle almost infinite universe
(up to 10^77^ x 10^77^, using a lightweight implementation of big integers coded on 256 bits, in `bitl.h` and `bitl.c`.)

However, it makes use of template lists and sets.
This makes the code an order of magniture slower than [the highly optimized golly](http://golly.sourceforge.net/).

The code `hgolbi.c` is commented, and refers to excerpts from the Gosper's paper.

This one is a tough one, so, hold on !

# General scheme

The paper "Exploiting regularities in large cellular spaces" written by R. Gosper in 1984 describes (with some details hidden between lines) an algorithm based on three main ideas:

- space contraction
- time contraction
- spacetime exploration

## Structure of the universe, and space contraction

The universe is structured by a tree of encapsulated sub-structures.

- First, the universe, a square in which cells will evolve following the rules of the Game of Life, is structured as a tree of quadrants:
  the universe is divided in 4 equal squares (named macro-cells by Gosper), each of theses parts divided in 4 squares as well,
  and so recursively down to squares of size 1.
- Two (or more) identical squares (in which cells are equally located inside the squares) share the same reference.
  This acts as a space contraction, where a universe with a repeated structure are modelized by very few macro-cells.
  Two macro-cells are identical if their four macro-cells are theirselves identical (and therefore share the same reference).
  A hash-mecanism, implemented by the function `macrocell_lt`, let identify identical macro-cells.
- Macrocells that don't contain any cell do all reference the same NULL macrocell.
  Therefore, the size of the tree only depend logarithmically on the number of cells in the universe,
  and not on the surface of the universe.

The function `universe_cell_set` let add cells in the universe, one by one, building the tree structure and
joining equivalent squares (by calls to `macrocell_patternify`).

## Universe evolution and time contraction

Once the universe has been populated with cells, it defines its initial configuration at time 0.
Those cells will evolve following the rules of the Game Of Life, for each time step.

In order to reduce (exponentially) the complexity (number of operations) to calculate the evolution of the universe:

- Two identical parts of the universe (square or macro-cells), even occuring at different times,
  and surrounded by the same environnement, will evolve the same way. The evolution will therefore be computed only once.
  For this, each macrocell records the result of its evolution after _t_ steps (where _t_ is one fourth the size of the macrocell),
  of its central cells (a square half the size of the macrocell).
- The evolution of a macro-cell at instant 2 x _t_ can be deduced, without extra computation,
  from the evolution at instant _t_ of 13 half-sized sub-macrocells (by a call to the function `universe_get_RESULT`).
- This is applied recursively down to the smallest macro-cells of size 2 by 2.
- Macro-cells of size 2 by 2 are therefore the only ones computed by direct application of the rules of the Game Of Life (by `next2x2`).

## Universe exploration

The function `universe_explore` will forecast the evolution of the universe at any instant of time.

- First, it makes sure the universe is large enough (by calling `universe_is_closed`) so that no cells will reach
  the border of the universe and escape between evolution in time (by sufficient calls to `universe_expand`).
- It decomposes the universe in four overlapping universes of same size _S_, where _S_ is the size of the initial universe. 
- It computes the evolution of these four universes at instant _S_/4 by calling `universe_get_RESULT`.
- It explores the configuration of the universe at any time between 0 and _S_/4 by calling `universe_show_RESULT`,
  by requesting information of their evolution to 13 half-sized sub-macro-cells (applied recursively.)

# Objects
## Universe
`Universe` is a public abstract data type (which internal definition is hidden from the user public interface).
Internally, it contains the reference to the top macro-cell, in which all cells are included.

## Explorer
`Explorer` is a public object type used by the function `universe_explore`.
It defines the region of the universe and the instant in time to explore, as well as the actions to perform before,
during and after the exploration.

```
Explorer
|\ spacetime
||\ space
|| \ window, region in space to explore
|| |\ NWvertex
|| | \ x, type: intbig_t
|| | \ y, type: intbig_t
||  \ SEvertex
||   \ x, type: intbig_t
||   \ y, type: intbig_t
| \ time, instant in time to explore
|  \ instant, type uintbig_t
 \ extractor, actions to perform
  \ preaction, type: void (*) (Universe *, SpaceTime, void *)
  \ foreach, type: void (*) (Universe *, SpaceTime, intbig_t, intbig_t, void *)
  \ postaction, type: void (*) (Universe *, SpaceTime, uintbig_t, void *)
  \ context, type void *
```

## Macro-cell
A `Macrocell` is a private object composed of:
- the references to its four sub-macro-cells (north-west, north-east, south-west and south-east squares.)
- the reference to the result of its evolution, itself a cenral sub-macro-cells.

```
Macrocell
\ quadrant, Array of 4 pointers to the sub-macrocells quadrants of type *Macrocell
\ result,  Pointer to the macrocell RESULT, of type *Macrocell
\ population, number of cells in the macro-cell, of type uintbig_t
\ nb_instances, number of instances of this macro-cell in space-time, of type uintbig_t
```

# Functions
## Global function
`next2x2` computes the evolution at next step of the 4 central cells of a 4 by 4 square, applying the rules of the Game Of Life.

## Main functions on Universe
### Public functions
- `universe_create` creates an empty new universe.
  It calls `xintbig_printf_init ()` for convenience.
- `universe_destroy` destroys a previously created universe (and releases all associated resources.)
- `universe_cell_set` add a cell into the universe at a given position at initial time 0.
- `universe_cell_unset` removes a cell from the universe at a given position at initial time 0.
- `universe_cell_is_set` controls the presence of a cell in the universe at a given position at initial time 0.
- `universe_set_BLE_rules` re-initialises and re-configures the initial time 0 of a universe for reuse and sets rule.
- `universe_RLE_readfile` re initialises and re-configures the initial time 0 of a universe from a given RLE file (see [Run Length Encoded](https://conwaylife.com/wiki/Run_Length_Encoded)).
- `universe_explore` find cells in a given region of space and at a given instant of time.

### Private functions
- `universe_get_RESULT` computes the RESULT of a macro-cell.
- `universe_is_closed` checks wether or not the universe is surrounded by empty space.
- `universe_expand` adds empty space around universe ; the size of the universe is doubled.
- `universe_cell_accessor` drills down into the tree of macro-cells of the universe to a given position.
- `universe_show_RESULT` finds cells from a macro-cell at a given time.

## Main functions on Macrocell
### Private functions
- `macrocell_get_cells_in_window`
- `macrocell_fetch_pattern`
- `macrocell_patternify`

# Usage

1. Include `hgolbi.h`.
1. Create a universe with `universe_create`.
1. Populate the universe:
    - with succesive calls to `universe_cell_set`.
    - from a file containing a RLE pattern with `universe_RLE_readfile`.
1. Declare actionners for exploration (at least `Explorer.extractor.foreach`).
1. Explore the universe through space and time (at least defining `Explorer.spacetime.time.instant`) with calls to `universe_explore`.
1. Destroy the universe with `universe_destroy`.

# Example

`hgolbi_example.c` is an example of usage of the hash-life algorithm.

To build it, type: `make`

`hgolbi.c` makes use of template lists and sets, which files can be found [here](https://github.com/farhiongit/Ctemplates).

`hgolbi_example` accepts options in command line (see `Makefile` for an example):

- `-t value` to set an instant in time,where `value` is a 2^64^ based non negative integer number.
- `-x min,max` to set a x-slice in space,where `min` and `max` are 2^64^ based signed integer numbers.
- `-y min,max` to set a y-slice in space,where `min` and `max` are 2^64^ based signed integer numbers.

2^64^ based integer numbers are composed of four or less components separated by `_`.
E.g.,

- `11_7` is equal to 11 x 2^64^ + 7.
- `-111_0_0` is equal to -111 x 2^128^.


**Have fun !**
