/*Daala video codec
Copyright (c) 2006-2013 Daala project contributors.  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.*/

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>

/*Encoder-only motion compensation routines.*/

/*TODO:
 - Develop a real encoding and measure real bits.
 - Thresholds for DP.
   + How do we calculate them?
   + How do they propagate between frames (block sizes change)*/


typedef struct od_mv_err_node od_mv_err_node;

#include "logging.h"
#include "mcenc.h"

typedef int od_offset[2];
typedef int od_pattern[8];

/*The number of bits to reduce chroma SADs by, if used.*/
#define OD_MC_CHROMA_SCALE (2)

/*The subdivision level of a MV in the mesh, given its position
   (mod OD_MVB_DELTA0).*/
static const int OD_MC_LEVEL[OD_MVB_DELTA0][OD_MVB_DELTA0] = {
  { 0, 6, 4, 6, 2, 6, 4, 6 },
  { 6, 5, 6, 5, 6, 5, 6, 5 },
  { 4, 6, 3, 6, 4, 6, 3, 6 },
  { 6, 5, 6, 5, 6, 5, 6, 5 },
  { 2, 6, 4, 6, 1, 6, 4, 6 },
  { 6, 5, 6, 5, 6, 5, 6, 5 },
  { 4, 6, 3, 6, 4, 6, 3, 6 },
  { 6, 5, 6, 5, 6, 5, 6, 5 }
};

/*Ancestor lists for a vertex.
  These are stored as lists of offsets to the vertices in the domain.
  All the ancestors at one level must be listed before ancestors at any higher
   level, and there must be no duplicates.
  By convention, we order each level top-to-bottom, then left-to-right.
  Level 0 ancestors are not included, as they cannot be decimated.
  We do, however, draw a square of level 0 vertices in the ASCII-art diagrams
   to provide some context on where each vertex lies.*/

/*Lists for level 2 vertices.*/
static const od_offset OD_ANCESTORS2[2][2] = {
  /*         1
             |
             |
             |
         0---2---0
         |   |   |
         |   |   |
         |   |   |
         |   1   |
         |       |
         |       |
         |       |
         0-------0 */
  { { 0, -4}, { 0,  4} },
  /*     0-------0
         |       |
         |       |
         |       |
     1---2---1   |
         |       |
         |       |
         |       |
         0-------0 */
  { {-4,  0}, { 4,  0} }
};

/*Lists for level 3 vertices.*/
static const od_offset OD_ANCESTORS3[4][5] = {
  /*         1
             |
             |
             |
         0---2---0
         |  /|   |
         | 3 |   |
         |/  |   |
     1---2---1   |
         |       |
         |       |
         |       |
         0-------0 */
  { { 2, -2}, {-2, 2}, { 2, -6}, {-6,  2}, { 2,  2} },
  /*         1
             |
             |
             |
         0---2---0
         |   |\  |
         |   | 3 |
         |   |  \|
         |   1---2---1
         |       |
         |       |
         |       |
         0-------0 */
  { {-2, -2}, { 2, 2}, {-2, -6}, {-2,  2}, { 6,  2} },
  /*     0-------0
         |       |
         |       |
         |       |
     1---2---1   |
         |\  |   |
         | 3 |   |
         |  \|   |
         0---2---0
             |
             |
             |
             1 */
  { {-2, -2}, { 2, 2}, {-6, -2}, { 2, -2}, { 2,  6} },
  /*     0-------0
         |       |
         |       |
         |       |
         |   1---2---1
         |   |  /|
         |   | 3 |
         |   |/  |
         0---2---0
             |
             |
             |
             1 */
  { { 2, -2}, {-2, 2}, {-2, -2}, { 6, -2}, {-2,  6} }
};

/*Lists for level 4 vertices.*/
static const od_offset OD_ANCESTORS4[8][9] = {
  /* 1---2---1
          \  |
           3 |
           |\|
         0-4-2---0
         | |/|   |
         | 3 |   |
         |/  |   |
     1---2---1   |
         |       |
         |       |
         |       |
         0-------0 */
  {
    { 0, -2}, { 0,  2}, {-2, -4}, { 2,  0}, {-2,  4},
    {-6, -4}, { 2, -4}, {-6,  4}, { 2,  4}
  },
  /*         1---2---1
             |  /
             | 3
             |/|
         0---2-4-0
         |   |\| |
         |   | 3 |
         |   |  \|
         |   1---2---1
         |       |
         |       |
         |       |
         0-------0 */
  {
    { 0, -2}, { 0,  2}, { 2, -4}, {-2,  0}, { 2,  4},
    {-2, -4}, { 6, -4}, {-2,  4}, { 6,  4}
  },
  /* 1       1
     |       |
     |       |
     |       |
     2   0---2---0
     |\  |  /|   |
     | 3-4-3 |   |
     |  \|/  |   |
     1---2---1   |
         |       |
         |       |
         |       |
         0-------0 */
  {
    {-2,  0}, { 2,  0}, {-4, -2}, { 4, -2}, { 0,  2},
    {-4, -6}, { 4, -6}, {-4,  2}, { 4,  2}
  },
  /*         1
             |
             |
             |
         0---2---0
         |  /|\  |
         | 3-4-3 |
         |/  |  \|
     1---2---1---2---1
         |       |
         |       |
         |       |
         0-------0 */
  {
    {-2,  0}, { 2,  0}, { 0, -2}, {-4,  2}, { 4,  2},
    { 0, -6}, {-8,  2}, { 0,  2}, { 8,  2}
  },
  /*         1
             |
             |
             |
         0---2---0
         |  /|   |
         | 3 |   |
         |/| |   |
     1---2-4-1   |
         |\| |   |
         | 3 |   |
         |  \|   |
         0---2---0
             |
             |
             |
             1 */
  {
    { 0, -2}, { 0,  2}, { 2, -4}, {-2,  0}, { 2,  4},
    { 2, -8}, {-6,  0}, { 2,  0}, { 2,  8}
  },
  /*         1
             |
             |
             |
         0---2---0
         |   |\  |
         |   | 3 |
         |   | |\|
         |   1-4-2---1
         |   | |/|
         |   | 3 |
         |   |/  |
         0---2---0
             |
             |
             |
             1 */
  {
    { 0, -2}, { 0,  2}, {-2, -4}, { 2,  0}, {-2,  4},
    {-2, -8}, {-2,  0}, { 6,  0}, {-2,  8}
  },
  /*     0-------0
         |       |
         |       |
         |       |
     1---2---1   |
     |  /|\  |   |
     | 3-4-3 |   |
     |/  |  \|   |
     2   0---2---0
     |       |
     |       |
     |       |
     1       1 */
  {
    {-2,  0}, { 2,  0}, { 0, -2}, {-4,  2}, { 4,  2},
    {-4, -2}, { 4, -2}, {-4,  6}, { 4,  6}
  },
  /*     0-------0
         |       |
         |       |
         |       |
     1---2---1---2---1
         |\  |  /|
         | 3-4-3 |
         |  \|/  |
         0---2---0
             |
             |
             |
             1 */
  {
    {-2,  0}, { 2,  0}, {-4, -2}, { 4, -2}, { 0,  2},
    {-8, -2}, { 0, -2}, { 8, -2}, { 0,  6}
  }
};

/*Lists for level 5 vertices.*/
static const od_offset OD_ANCESTORS5[16][14] = {
  /* 1---2---1
     |    \  |
     |     3 |
     |     |\|
     2   0-4-2---0
     |\  |5|/|   |
     | 3-4-3 |   |
     |  \|/  |   |
     1---2---1   |
         |       |
         |       |
         |       |
         0-------0 */
  {
    { 1, -1}, {-1,  1}, { 1, -3}, {-3,  1}, { 1,  1}, {-1, -5}, {-5, -1},
    { 3, -1}, {-1,  3}, {-5, -5}, { 3, -5}, {-5,  3}, { 3,  3}
  },
  /* 1---2---1
          \  |
           3 |
           |\|
         0-4-2---0
         | |5|\  |
         | 3-4-3 |
         |/  |  \|
     1---2---1---2---1
         |       |
         |       |
         |       |
         0-------0 */
  {
    {-1, -1}, { 1,  1}, {-1, -3}, {-1,  1}, { 3,  1}, {-3, -5}, { 1, -1},
    {-3,  3}, { 5,  3}, {-7, -5}, { 1, -5}, {-7,  3}, { 1,  3}, { 9,  3}
  },
  /*         1---2---1
             |  /
             | 3
             |/|
         0---2-4-0
         |  /|5| |
         | 3-4-3 |
         |/  |  \|
     1---2---1---2---1
         |       |
         |       |
         |       |
         0-------0 */
  {
    { 1, -1}, {-1,  1}, { 1, -3}, {-3,  1}, { 1,  1}, { 3, -5}, {-1, -1},
    {-5,  3}, { 3,  3}, {-1, -5}, { 7, -5}, {-9,  3}, {-1,  3}, { 7,  3}
  },
  /*         1---2---1
             |  /    |
             | 3     |
             |/|     |
         0---2-4-0---2
         |   |\|5|  /|
         |   | 3-4-3 |
         |   |  \|/  |
         |   1---2---1
         |       |
         |       |
         |       |
         0-------0 */
  {
    {-1, -1}, { 1,  1}, {-1, -3}, {-1,  1}, { 3,  1}, { 1, -5}, {-3, -1},
    { 5, -1}, { 1,  3}, {-3, -5}, { 5, -5}, {-3,  3}, { 5,  3}
  },
  /* 1       1
     |       |
     |       |
     |       |
     2   0---2---0
     |\  |  /|   |
     | 3-4-3 |   |
     |  \|5| |   |
     1---2-4-1   |
         |\| |   |
         | 3 |   |
         |  \|   |
         0---2---0
             |
             |
             |
             1 */
  {
    {-1, -1}, { 1,  1}, {-3, -1}, { 1, -1}, { 1,  3}, {-5, -3}, { 3, -3},
    {-1,  1}, { 3,  5}, {-5, -7}, { 3, -7}, {-5,  1}, { 3,  1}, { 3,  9}
  },
  /*         1
             |
             |
             |
         0---2---0
         |  /|\  |
         | 3-4-3 |
         |/|5|  \|
     1---2-4-1---2---1
         |\| |   |
         | 3 |   |
         |  \|   |
         0---2---0
             |
             |
             |
             1 */
  {
    { 1, -1}, {-1,  1}, {-1, -1}, { 3, -1}, {-1,  3}, { 1, -3}, {-3,  1},
    { 5,  1}, { 1,  5}, { 1, -7}, {-7,  1}, { 1,  1}, { 9,  1}, { 1,  9}
  },
  /*         1
             |
             |
             |
         0---2---0
         |  /|\  |
         | 3-4-3 |
         |/  |5|\|
     1---2---1-4-2---1
         |   | |/|
         |   | 3 |
         |   |/  |
         0---2---0
             |
             |
             |
             1 */
  {
    {-1, -1}, { 1,  1}, {-3, -1}, { 1, -1}, { 1,  3}, {-1, -3}, {-5,  1},
    { 3,  1}, {-1,  5}, {-1, -7}, {-9,  1}, {-1,  1}, { 7,  1}, {-1,  9}
  },
  /*         1       1
             |       |
             |       |
             |       |
         0---2---0   2
         |   |\  |  /|
         |   | 3-4-3 |
         |   | |5|/  |
         |   1-4-2---1
         |   | |/|
         |   | 3 |
         |   |/  |
         0---2---0
             |
             |
             |
             1 */
  {
    { 1, -1}, {-1,  1}, {-1, -1}, { 3, -1}, {-1,  3}, {-3, -3}, { 5, -3},
    { 1,  1}, {-3,  5}, {-3, -7}, { 5, -7}, {-3,  1}, { 5,  1}, {-3,  9}
  },
  /*         1
             |
             |
             |
         0---2---0
         |  /|   |
         | 3 |   |
         |/| |   |
     1---2-4-1   |
     |  /|5| |   |
     | 3-4-3 |   |
     |/  |  \|   |
     2   0---2---0
     |       |
     |       |
     |       |
     1       1 */
  {
    { 1, -1}, {-1,  1}, { 1, -3}, {-3,  1}, { 1,  1}, { 3, -5}, {-1, -1},
    {-5,  3}, { 3,  3}, { 3, -9}, {-5, -1}, { 3, -1}, {-5,  7}, { 3,  7}
  },
  /*         1
             |
             |
             |
         0---2---0
         |  /|   |
         | 3 |   |
         |/| |   |
     1---2-4-1---2---1
         |\|5|  /|
         | 3-4-3 |
         |  \|/  |
         0---2---0
             |
             |
             |
             1 */
  {
    {-1, -1}, { 1,  1}, {-1, -3}, {-1,  1}, { 3,  1}, { 1, -5}, {-3, -1},
    { 5, -1}, { 1,  3}, { 1, -9}, {-7, -1}, { 1, -1}, { 9, -1}, { 1,  7}
  },
  /*         1
             |
             |
             |
         0---2---0
         |   |\  |
         |   | 3 |
         |   | |\|
     1---2---1-4-2---1
         |\  |5|/|
         | 3-4-3 |
         |  \|/  |
         0---2---0
             |
             |
             |
             1 */
  {
    { 1, -1}, {-1,  1}, { 1, -3}, {-3,  1}, { 1,  1}, {-1, -5}, {-5, -1},
    { 3, -1}, {-1,  3}, {-1, -9}, {-9, -1}, {-1, -1}, { 7, -1}, {-1,  7}
  },
  /*         1
             |
             |
             |
         0---2---0
         |   |\  |
         |   | 3 |
         |   | |\|
         |   1-4-2---1
         |   | |5|\  |
         |   | 3-4-3 |
         |   |/  |  \|
         0---2---0   2
             |       |
             |       |
             |       |
             1       1 */
  {
    {-1, -1}, { 1,  1}, {-1, -3}, {-1,  1}, { 3,  1}, {-3, -5}, { 1, -1},
    {-3,  3}, { 5,  3}, {-3, -9}, {-3, -1}, { 5, -1}, {-3,  7}, { 5,  7}
  },
  /*     0-------0
         |       |
         |       |
         |       |
     1---2---1   |
     |  /|\  |   |
     | 3-4-3 |   |
     |/  |5|\|   |
     2   0-4-2---0
     |     |/|
     |     3 |
     |    /  |
     1---2---1 */
  {
    {-1, -1}, { 1,  1}, {-3, -1}, { 1, -1}, { 1,  3}, {-1, -3}, {-5,  1},
    { 3,  1}, {-1,  5}, {-5, -3}, { 3, -3}, {-5,  5}, { 3,  5}
  },
  /*     0-------0
         |       |
         |       |
         |       |
     1---2---1---2---1
         |\  |  /|
         | 3-4-3 |
         | |5|/  |
         0-4-2---0
           |/|
           3 |
          /  |
     1---2---1 */
  {
    { 1, -1}, {-1,  1}, {-1, -1}, { 3, -1}, {-1,  3}, {-3, -3}, { 5, -3},
    { 1,  1}, {-3,  5}, {-7, -3}, { 1, -3}, { 9, -3}, {-7,  5}, { 1,  5}
  },
  /*     0-------0
         |       |
         |       |
         |       |
     1---2---1---2---1
         |\  |  /|
         | 3-4-3 |
         |  \|5| |
         0---2-4-0
             |\|
             | 3
             |  \
             1---2---1 */
  {
    {-1, -1}, { 1,  1}, {-3, -1}, { 1, -1}, { 1,  3}, {-5, -3}, { 3, -3},
    {-1,  1}, { 3,  5}, {-9, -3}, {-1, -3}, { 7, -3}, {-1,  5}, { 7,  5}
  },
  /*     0-------0
         |       |
         |       |
         |       |
         |   1---2---1
         |   |  /|\  |
         |   | 3-4-3 |
         |   |/|5|  \|
         0---2-4-0   2
             |\|     |
             | 3     |
             |  \    |
             1---2---1 */
  {
    { 1, -1}, {-1,  1}, {-1, -1}, { 3, -1}, {-1,  3}, { 1, -3}, {-3,  1},
    { 5,  1}, { 1,  5}, {-3, -3}, { 5, -3}, {-3,  5}, { 5,  5}
  }
};

/*Lists for level 6 vertices.*/
static const od_offset OD_ANCESTORS6[32][20] = {
  /* 1---2---1
     |  / \  |
     | 3-4-3 |
     |/   5|\|
     2   064-2---0
     |\  |5|/|   |
     | 3-4-3 |   |
     |  \|/  |   |
     1---2---1   |
         |       |
         |       |
         |       |
         0-------0 */
  {
    { 0, -1}, { 0,  1}, {-1, -2}, { 1,  0}, {-1,  2}, {-3, -2}, { 1, -2},
    {-3,  2}, { 1,  2}, {-1, -4}, {-5,  0}, { 3,  0}, {-1,  4},
    {-5, -4}, { 3, -4}, {-5,  4}, { 3,  4}
  },
  /* 1---2---1---2---1
          \  |  /
           3-4-3
           |5|/
         0-462---0
         | |5|\  |
         | 3-4-3 |
         |/  |  \|
     1---2---1---2---1
         |       |
         |       |
         |       |
         0-------0 */
  {
    { 0, -1}, { 0,  1}, { 1, -2}, {-1,  0}, { 1,  2}, {-1, -2}, { 3, -2},
    {-1,  2}, { 3,  2}, {-3, -4}, { 5, -4}, { 1,  0}, {-3,  4}, { 5,  4},
    {-7, -4}, { 1, -4}, { 9, -4}, {-7,  4}, { 1,  4}, { 9,  4}
  },
  /* 1---2---1---2---1
          \  |  /
           3-4-3
            \|5|
         0---264-0
         |  /|5| |
         | 3-4-3 |
         |/  |  \|
     1---2---1---2---1
         |       |
         |       |
         |       |
         0-------0 */
  {
    { 0, -1}, { 0,  1}, {-1, -2}, { 1,  0}, {-1,  2}, {-3, -2}, { 1, -2},
    {-3,  2}, { 1,  2}, {-5, -4}, { 3, -4}, {-1,  0}, {-5,  4}, { 3,  4},
    {-9, -4}, {-1, -4}, { 7, -4}, {-9,  4}, {-1,  4}, { 7,  4}
  },
  /*         1---2---1
             |  /|\  |
             | 3-4-3 |
             |/|5   \|
         0---2-460   2
         |   |\|5|  /|
         |   | 3-4-3 |
         |   |   |/  |
         |   1---2---1
         |       |
         |       |
         |       |
         0-------0 */
  {
    { 0, -1}, { 0,  1}, { 1, -2}, {-1,  0}, { 1,  2}, {-1, -2}, { 3, -2},
    {-1,  2}, { 3,  2}, { 1, -4}, {-3,  0}, { 5,  0}, { 1,  4},
    {-3, -4}, { 5, -4}, {-3,  4}, { 5,  4}
  },
  /* 1---2---1
     |  / \  |
     | 3   3 |
     |/|   |\|
     2 4 0-4-2---0
     |\|565|/|   |
     | 3-4-3 |   |
     |  \|/  |   |
     1---2---1   |
         |       |
         |       |
         |       |
         0-------0 */
  {
    {-1,  0}, { 1,  0}, {-2, -1}, { 2, -1}, { 0,  1}, {-2, -3}, { 2, -3},
    {-2,  1}, { 2,  1}, { 0, -5}, {-4, -1}, { 4, -1}, { 0,  3},
    {-4, -5}, { 4, -5}, {-4,  3}, { 4,  3}
  },
  /* 1---2---1
     |    \  |
     |     3 |
     |     |\|
     2   0-4-2---0
     |\  |565|\  |
     | 3-4-3-4-3 |
     |  \|/  |  \|
     1---2---1---2---1
         |       |
         |       |
         |       |
         0-------0 */
  {
    {-1,  0}, { 1,  0}, { 0, -1}, {-2,  1}, { 2,  1}, { 0, -3}, {-4,  1},
    { 0,  1}, { 4,  1}, {-2, -5}, {-6, -1}, { 2, -1}, {-2,  3}, { 6,  3},
    {-6, -5}, { 2, -5}, {-6,  3}, { 2,  3}, {10,  3}
  },
  /* 1---2---1---2---1
          \  |  /
           3 | 3
           |\|/|
         0-4-2-4-0
         | |565| |
         | 3-4-3 |
         |/  |  \|
     1---2---1---2---1
         |       |
         |       |
         |       |
         0-------0 */
  {
    {-1,  0}, { 1,  0}, {-2, -1}, { 2, -1}, { 0,  1}, {-2, -3}, { 2, -3},
    {-2,  1}, { 2,  1}, {-4, -5}, { 4, -5}, { 0, -1}, {-4,  3}, { 4,  3},
    {-8, -5}, { 0, -5}, { 8, -5}, {-8,  3}, { 0,  3}, { 8,  3}
  },
  /*         1---2---1
             |  /    |
             | 3     |
             |/|     |
         0---2-4-0   2
         |  /|565|  /|
         | 3-4-3-4-3 |
         |/  |  \|/  |
     1---2---1---2---1
         |       |
         |       |
         |       |
         0-------0 */
  {
    {-1,  0}, { 1,  0}, { 0, -1}, {-2,  1}, { 2,  1}, { 0, -3}, {-4,  1},
    { 0,  1}, { 4,  1}, { 2, -5}, {-2, -1}, { 6, -1}, {-6,  3}, { 2,  3},
    {-2, -5}, { 6, -5}, {-10, 3}, {-2,  3}, { 6,  3}
  },
  /* 1---2---1
     |    \  |
     |     3 |
     |     |\|
     2   0-4-2---0
     |\  |5|/|   |
     | 3-463 |   |
     |  \|5| |   |
     1---2-4-1   |
         |\| |   |
         | 3 |   |
         |  \|   |
         0---2---0
             |
             |
             |
             1 */
  {
    { 0, -1}, { 0,  1}, { 1, -2}, {-1,  0}, { 1,  2}, { 1, -4}, {-3,  0},
    { 1,  0}, { 1,  4}, {-1, -6}, {-5, -2}, { 3, -2}, {-1,  2}, { 3,  6},
    {-5, -6}, { 3, -6}, {-5,  2}, { 3,  2}, { 3, 10}
  },
  /* 1---2---1
          \  |
           3 |
           |\|
         0-4-2---0
         | |5|\  |
         | 364-3 |
         |/|5|  \|
     1---2-4-1---2---1
         |\| |   |
         | 3 |   |
         |  \|   |
         0---2---0
             |
             |
             |
             1 */
  {
    { 0, -1}, { 0,  1}, {-1, -2}, { 1,  0}, {-1,  2}, {-1, -4}, {-1,  0},
    { 3,  0}, {-1,  4}, {-3, -6}, { 1, -2}, {-3,  2}, { 5,  2}, { 1,  6},
    {-7, -6}, { 1, -6}, {-7,  2}, { 1,  2}, { 9,  2}, { 1, 10}
  },
  /*         1---2---1
             |  /
             | 3
             |/|
         0---2-4-0
         |  /|5| |
         | 3-463 |
         |/  |5|\|
     1---2---1-4-2---1
         |   | |/|
         |   | 3 |
         |   |/  |
         0---2---0
             |
             |
             |
             1 */
  {
    { 0, -1}, { 0,  1}, { 1, -2}, {-1,  0}, { 1,  2}, { 1, -4}, {-3,  0},
    { 1,  0}, { 1,  4}, { 3, -6}, {-1, -2}, {-5,  2}, { 3,  2}, {-1,  6},
    {-1, -6}, { 7, -6}, {-9,  2}, {-1,  2}, { 7,  2}, {-1, 10}
  },
  /*         1---2---1
             |  /    |
             | 3     |
             |/|     |
         0---2-4-0   2
         |   | |5|  /|
         |   | 364-3 |
         |   | |5|/  |
         |   1-4-2---1
         |   | |/|
         |   | 3 |
         |   |/  |
         0---2---0
             |
             |
             |
             1 */
  {
    { 0, -1}, { 0,  1}, {-1, -2}, { 1,  0}, {-1,  2}, {-1, -4}, {-1,  0},
    { 3,  0}, {-1,  4}, { 1, -6}, {-3, -2}, { 5, -2}, { 1,  2}, {-3,  6},
    {-3, -6}, { 5, -6}, {-3,  2}, { 5,  2}, {-3, 10}
  },
  /* 1       1
     |       |
     |       |
     |       |
     2   0---2---0
     |\  |  /|   |
     | 3-4-3 |   |
     | |565| |   |
     1-4-2-4-1   |
     | |/|\| |   |
     | 3 | 3 |   |
     |/  |  \|   |
     2   0---2---0
     |       |
     |       |
     |       |
     1       1 */
  {
    {-1,  0}, { 1,  0}, { 0, -1}, {-2,  1}, { 2,  1}, {-2, -1}, { 2, -1},
    {-2,  3}, { 2,  3}, {-4, -3}, { 4, -3}, { 0,  1}, {-4,  5}, { 4,  5},
    {-4, -7}, { 4, -7}, {-4,  1}, { 4,  1}, {-4,  9}, { 4,  9}
  },
  /* 1       1
     |       |
     |       |
     |       |
     2   0---2---0
     |\  |  /|\  |
     | 3-4-3-4-3 |
     |  \|565|  \|
     1---2-4-1---2---1
         |\| |   |
         | 3 |   |
         |  \|   |
         0---2---0
             |
             |
             |
             1 */
  {
    {-1,  0}, { 1,  0}, {-2, -1}, { 2, -1}, { 0,  1}, {-4, -1}, { 0, -1},
    { 4, -1}, { 0,  3}, {-6, -3}, { 2, -3}, {-2,  1}, { 6,  1}, { 2,  5},
    {-6, -7}, { 2, -7}, {-6,  1}, { 2,  1}, {10,  1}, { 2,  9}
  },
  /*         1
             |
             |
             |
         0---2---0
         |  /|\  |
         | 3-4-3 |
         |/|565|\|
     1---2-4-1-4-2---1
         |\| | |/|
         | 3 | 3 |
         |  \|/  |
         0---2---0
             |
             |
             |
             1 */
  {
    {-1,  0}, { 1,  0}, { 0, -1}, {-2,  1}, { 2,  1}, {-2, -1}, { 2, -1},
    {-2,  3}, { 2,  3}, { 0, -3}, {-4,  1}, { 4,  1}, { 0,  5},
    { 0, -7}, {-8,  1}, { 0,  1}, { 8,  1}, { 0,  9}
  },
  /*         1       1
             |       |
             |       |
             |       |
         0---2---0   2
         |  /|\  |  /|
         | 3-4-3-4-3 |
         |/  |565|/  |
     1---2---1-4-2---1
         |   | |/|
         |   | 3 |
         |   |/  |
         0---2---0
             |
             |
             |
             1 */
  {
    {-1,  0}, { 1,  0}, {-2, -1}, { 2, -1}, { 0,  1}, {-4, -1}, { 0, -1},
    { 4, -1}, { 0,  3}, {-2, -3}, { 6, -3}, {-6,  1}, { 2,  1}, {-2,  5},
    {-2, -7}, { 6, -7}, {-10, 1}, {-2,  1}, { 6,  1}, {-2,  9}
  },
  /* 1       1
     |       |
     |       |
     |       |
     2   0---2---0
     |\  |  /|   |
     | 3-4-3 |   |
     |  \|5| |   |
     1---264-1   |
     |  /|5| |   |
     | 3-4-3 |   |
     |/  |  \|   |
     2   0---2---0
     |       |
     |       |
     |       |
     1       1 */
  {
    { 0, -1}, { 0,  1}, {-1, -2}, { 1,  0}, {-1,  2}, {-3, -2}, { 1, -2},
    {-3,  2}, { 1,  2}, {-5, -4}, { 3, -4}, {-1,  0}, {-5,  4}, { 3,  4},
    {-5, -8}, { 3, -8}, {-5,  0}, { 3,  0}, {-5,  8}, { 3,  8}
  },
  /*         1
             |
             |
             |
         0---2---0
         |  /|\  |
         | 3-4-3 |
         |/|5|  \|
     1---2-461---2---1
         |\|5|  /|
         | 3-4-3 |
         |  \|/  |
         0---2---0
             |
             |
             |
             1 */
  {
    { 0, -1}, { 0,  1}, { 1, -2}, {-1,  0}, { 1,  2}, {-1, -2}, { 3, -2},
    {-1,  2}, { 3,  2}, { 1, -4}, {-3,  0}, { 5,  0}, { 1,  4},
    { 1, -8}, {-7,  0}, { 1,  0}, { 9,  0}, { 1,  8}
  },
  /*         1
             |
             |
             |
         0---2---0
         |  /|\  |
         | 3-4-3 |
         |/  |5|\|
     1---2---164-2---1
         |\  |5|/|
         | 3-4-3 |
         |  \|/  |
         0---2---0
             |
             |
             |
             1 */
  {
    { 0, -1}, { 0,  1}, {-1, -2}, { 1,  0}, {-1,  2}, {-3, -2}, { 1, -2},
    {-3,  2}, { 1,  2}, {-1, -4}, {-5,  0}, { 3,  0}, {-1,  4},
    {-1, -8}, {-9,  0}, {-1,  0}, { 7,  0}, {-1,  8}
  },
  /*         1       1
             |       |
             |       |
             |       |
         0---2---0   2
         |   |\  |  /|
         |   | 3-4-3 |
         |   | |5|/  |
         |   1-462---1
         |   | |5|\  |
         |   | 3-4-3 |
         |   |/  |  \|
         0---2---0   2
             |       |
             |       |
             |       |
             1       1 */
  {
    { 0, -1}, { 0,  1}, { 1, -2}, {-1,  0}, { 1,  2}, {-1, -2}, { 3, -2},
    {-1,  2}, { 3,  2}, {-3, -4}, { 5, -4}, { 1,  0}, {-3,  4}, { 5,  4},
    {-3, -8}, { 5, -8}, {-3,  0}, { 5,  0}, {-3,  8}, { 5,  8}
  },
  /* 1       1
     |       |
     |       |
     |       |
     2   0---2---0
     |\  |  /|   |
     | 3 | 3 |   |
     | |\|/| |   |
     1-4-2-4-1   |
     | |565| |   |
     | 3-4-3 |   |
     |/  |  \|   |
     2   0---2---0
     |       |
     |       |
     |       |
     1       1 */
  {
    {-1,  0}, { 1,  0}, {-2, -1}, { 2, -1}, { 0,  1}, {-2, -3}, { 2, -3},
    {-2,  1}, { 2,  1}, {-4, -5}, { 4, -5}, { 0, -1}, {-4,  3}, { 4,  3},
    {-4, -9}, { 4, -9}, {-4, -1}, { 4, -1}, {-4,  7}, { 4,  7}
  },
  /*         1
             |
             |
             |
         0---2---0
         |  /|   |
         | 3 |   |
         |/| |   |
     1---2-4-1---2---1
     |  /|565|  /|
     | 3-4-3-4-3 |
     |/  |  \|/  |
     2   0---2---0
     |       |
     |       |
     |       |
     1       1 */
  {
    {-1,  0}, { 1,  0}, { 0, -1}, {-2,  1}, { 2,  1}, { 0, -3}, {-4,  1},
    { 0,  1}, { 4,  1}, { 2, -5}, {-2, -1}, { 6, -1}, {-6,  3}, { 2,  3},
    { 2, -9}, {-6, -1}, { 2, -1}, {10, -1}, {-6,  7}, { 2,  7}
  },
  /*         1
             |
             |
             |
         0---2---0
         |  /|\  |
         | 3 | 3 |
         |/| | |\|
     1---2-4-1-4-2---1
         |\|565|/|
         | 3-4-3 |
         |  \|/  |
         0---2---0
             |
             |
             |
             1 */
  {
    {-1,  0}, { 1,  0}, {-2, -1}, { 2, -1}, { 0,  1}, {-2, -3}, { 2, -3},
    {-2,  1}, { 2,  1}, { 0, -5}, {-4, -1}, { 4, -1}, { 0,  3},
    { 0, -9}, {-8, -1}, { 0, -1}, { 8, -1}, { 0,  7}
  },
  /*         1
             |
             |
             |
         0---2---0
         |   |\  |
         |   | 3 |
         |   | |\|
     1---2---1-4-2---1
         |\  |565|\  |
         | 3-4-3-4-3 |
         |  \|/  |  \|
         0---2---0   2
             |       |
             |       |
             |       |
             1       1 */
  {
    {-1,  0}, { 1,  0}, { 0, -1}, {-2,  1}, { 2,  1}, { 0, -3}, {-4,  1},
    { 0,  1}, { 4,  1}, {-2, -5}, {-6, -1}, { 2, -1}, {-2,  3}, { 6,  3},
    {-2, -9}, {-10,-1}, {-2, -1}, { 6, -1}, {-2,  7}, { 6,  7}
  },
  /*         1
             |
             |
             |
         0---2---0
         |  /|   |
         | 3 |   |
         |/| |   |
     1---2-4-1   |
     |  /|5| |   |
     | 3-463 |   |
     |/  |5|\|   |
     2   0-4-2---0
     |     |/|
     |     3 |
     |    /  |
     1---2---1 */
  {
    { 0, -1}, { 0,  1}, { 1, -2}, {-1,  0}, { 1,  2}, { 1, -4}, {-3,  0},
    { 1,  0}, { 1,  4}, { 3, -6}, {-1, -2}, {-5,  2}, { 3,  2}, {-1,  6},
    {3, -10}, {-5, -2}, { 3, -2}, {-5,  6}, { 3,  6}
  },
  /*         1
             |
             |
             |
         0---2---0
         |  /|   |
         | 3 |   |
         |/| |   |
     1---2-4-1---2---1
         |\|5|  /|
         | 364-3 |
         | |5|/  |
         0-4-2---0
           |/|
           3 |
          /  |
     1---2---1 */
  {
    { 0, -1}, { 0,  1}, {-1, -2}, { 1,  0}, {-1,  2}, {-1, -4}, {-1,  0},
    { 3,  0}, {-1,  4}, { 1, -6}, {-3, -2}, { 5, -2}, { 1,  2}, {-3,  6},
    {1, -10}, {-7, -2}, { 1, -2}, { 9, -2}, {-7,  6}, { 1,  6}
  },
  /*         1
             |
             |
             |
         0---2---0
         |   |\  |
         |   | 3 |
         |   | |\|
     1---2---1-4-2---1
         |\  |5|/|
         | 3-463 |
         |  \|5| |
         0---2-4-0
             |\|
             | 3
             |  \
             1---2---1 */
  {
    { 0, -1}, { 0,  1}, { 1, -2}, {-1,  0}, { 1,  2}, { 1, -4}, {-3,  0},
    { 1,  0}, { 1,  4}, {-1, -6}, {-5, -2}, { 3, -2}, {-1,  2}, { 3,  6},
    {-1,-10}, {-9, -2}, {-1, -2}, { 7, -2}, {-1,  6}, { 7,  6}
  },
  /*         1
             |
             |
             |
         0---2---0
         |   |\  |
         |   | 3 |
         |   | |\|
         |   1-4-2---1
         |   | |5|\  |
         |   | 364-3 |
         |   |/|5|  \|
         0---2-4-0   2
             |\|     |
             | 3     |
             |  \    |
             1---2---1 */
  {
    { 0, -1}, { 0,  1}, {-1, -2}, { 1,  0}, {-1,  2}, {-1, -4}, {-1,  0},
    { 3,  0}, {-1,  4}, {-3, -6}, { 1, -2}, {-3,  2}, { 5,  2}, { 1,  6},
    {-3,-10}, {-3, -2}, { 5, -2}, {-3,  6}, { 5,  6}
  },
  /*     0-------0
         |       |
         |       |
         |       |
     1---2---1   |
     |  /|\  |   |
     | 3-4-3 |   |
     |/|565|\|   |
     2 4 0-4-2---0
     |\|   |/|
     | 3   3 |
     |  \ /  |
     1---2---1 */
  {
    {-1,  0}, { 1,  0}, { 0, -1}, {-2,  1}, { 2,  1}, {-2, -1}, { 2, -1},
    {-2,  3}, { 2,  3}, { 0, -3}, {-4,  1}, { 4,  1}, { 0,  5},
    {-4, -3}, { 4, -3}, {-4,  5}, { 4,  5}
  },
  /*     0-------0
         |       |
         |       |
         |       |
     1---2---1---2---1
     |  /|\  |  /|
     | 3-4-3-4-3 |
     |/  |565|/  |
     2   0-4-2---0
     |     |/|
     |     3 |
     |    /  |
     1---2---1 */
  {
    {-1,  0}, { 1,  0}, {-2, -1}, { 2, -1}, { 0,  1}, {-4, -1}, { 0, -1},
    { 4, -1}, { 0,  3}, {-2, -3}, { 6, -3}, {-6,  1}, { 2,  1}, {-2,  5},
    {-6, -3}, { 2, -3}, {10, -3}, {-6,  5}, { 2,  5}
  },
  /*     0-------0
         |       |
         |       |
         |       |
     1---2---1---2---1
         |\  |  /|
         | 3-4-3 |
         | |565| |
         0-4-2-4-0
           |/|\|
           3 | 3
          /  |  \
     1---2---1---2---1*/
  {
    {-1,  0}, { 1,  0}, { 0, -1}, {-2,  1}, { 2,  1}, {-2, -1}, { 2, -1},
    {-2,  3}, { 2,  3}, {-4, -3}, { 4, -3}, { 0,  1}, {-4,  5}, { 4,  5},
    {-8, -3}, { 0, -3}, { 8, -3}, {-8,  5}, { 0,  5}, { 8,  5}
  },
  /*     0-------0
         |       |
         |       |
         |       |
     1---2---1---2---1
         |\  |  /|\  |
         | 3-4-3-4-3 |
         |  \|565|  \|
         0---2-4-0   2
             |\|     |
             | 3     |
             |  \    |
             1---2---1 */
  {
    {-1,  0}, { 1,  0}, {-2, -1}, { 2, -1}, { 0,  1}, {-4, -1}, { 0, -1},
    { 4, -1}, { 0,  3}, {-6, -3}, { 2, -3}, {-2,  1}, { 6,  1}, { 2,  5},
    {-10,-3}, {-2, -3}, { 6, -3}, {-2,  5}, { 6,  5}
  }
};

/*The number of unique ancestors in each list in the grid pattern.
  Up through and including level 4, this number is constant for all vertices in
   the same level, but for levels beyond that, it depends on the location in
   the grid.*/
static const int OD_NANCESTORS[OD_MVB_DELTA0][OD_MVB_DELTA0] = {
  {  0, 17,  9, 20,  2, 20,  9, 17 },
  { 17, 13, 19, 14, 20, 14, 19, 13 },
  {  9, 19,  5, 20,  9, 20,  5, 19 },
  { 20, 14, 20, 14, 18, 14, 20, 14 },
  {  2, 20,  9, 18,  0, 18,  9, 20 },
  { 20, 14, 20, 14, 18, 14, 20, 14 },
  {  9, 19,  5, 20,  9, 20,  5, 19 },
  { 17, 13, 19, 14, 20, 14, 19, 13 }
};

/*The lists for each vertex in the grid pattern.*/
static const od_offset *OD_ANCESTORS[OD_MVB_DELTA0][OD_MVB_DELTA0] = {
  {
                 NULL,  OD_ANCESTORS6[0],  OD_ANCESTORS4[0],  OD_ANCESTORS6[1],
     OD_ANCESTORS2[0],  OD_ANCESTORS6[2],  OD_ANCESTORS4[1],  OD_ANCESTORS6[3]
  },
  {
     OD_ANCESTORS6[4],  OD_ANCESTORS5[0],  OD_ANCESTORS6[5],  OD_ANCESTORS5[1],
     OD_ANCESTORS6[6],  OD_ANCESTORS5[2],  OD_ANCESTORS6[7],  OD_ANCESTORS5[3]
  },
  {
     OD_ANCESTORS4[2],  OD_ANCESTORS6[8],  OD_ANCESTORS3[0],  OD_ANCESTORS6[9],
     OD_ANCESTORS4[3], OD_ANCESTORS6[10],  OD_ANCESTORS3[1], OD_ANCESTORS6[11]
  },
  {
    OD_ANCESTORS6[12],  OD_ANCESTORS5[4], OD_ANCESTORS6[13],  OD_ANCESTORS5[5],
    OD_ANCESTORS6[14],  OD_ANCESTORS5[6], OD_ANCESTORS6[15],  OD_ANCESTORS5[7]
  },
  {
     OD_ANCESTORS2[1], OD_ANCESTORS6[16],  OD_ANCESTORS4[4], OD_ANCESTORS6[17],
                 NULL, OD_ANCESTORS6[18],  OD_ANCESTORS4[5], OD_ANCESTORS6[19]
  },
  {
    OD_ANCESTORS6[20],  OD_ANCESTORS5[8], OD_ANCESTORS6[21],  OD_ANCESTORS5[9],
    OD_ANCESTORS6[22], OD_ANCESTORS5[10], OD_ANCESTORS6[23], OD_ANCESTORS5[11]
  },
  {
     OD_ANCESTORS4[6], OD_ANCESTORS6[24],  OD_ANCESTORS3[2], OD_ANCESTORS6[25],
     OD_ANCESTORS4[7], OD_ANCESTORS6[26],  OD_ANCESTORS3[3], OD_ANCESTORS6[27]
  },
  {
    OD_ANCESTORS6[28], OD_ANCESTORS5[12], OD_ANCESTORS6[29], OD_ANCESTORS5[13],
    OD_ANCESTORS6[30], OD_ANCESTORS5[14], OD_ANCESTORS6[31], OD_ANCESTORS5[15]
  }
};

/*Computes the Sum of Absolute Differences: slow path.*/
int32_t od_mc_compute_sad8_c(const unsigned char *src, int systride,
 const unsigned char *ref, int dystride, int w, int h) {
  int i;
  int j;
  int ret;
  ret = 0;
  for (j = 0; j < h; j++) {
    for (i = 0; i < w; i++) {
      ret += abs(ref[i] - src[i]);
    }
    src += systride;
    ref += dystride;
  }
  return ret;
}

int32_t od_mc_compute_sad8_4x4_c(const unsigned char *src, int systride,
 const unsigned char *ref, int dystride) {
  return od_mc_compute_sad8_c(src, systride, ref, dystride, 4, 4);
}

int32_t od_mc_compute_sad8_8x8_c(const unsigned char *src, int systride,
 const unsigned char *ref, int dystride) {
  return od_mc_compute_sad8_c(src, systride, ref, dystride, 8, 8);
}

int32_t od_mc_compute_sad8_16x16_c(const unsigned char *src, int systride,
 const unsigned char *ref, int dystride) {
  return od_mc_compute_sad8_c(src, systride, ref, dystride, 16, 16);
}

int32_t od_mc_compute_sad8_32x32_c(const unsigned char *src, int systride,
 const unsigned char *ref, int dystride) {
  return od_mc_compute_sad8_c(src, systride, ref, dystride, 32, 32);
}

int32_t od_mc_compute_sad16_c(const unsigned char *src, int systride,
 const unsigned char *ref, int dystride, int w, int h) {
  int i;
  int j;
  int32_t ret;
  ret = 0;
  for (j = 0; j < h; j++) {
    for (i = 0; i < w; i++) {
      ret += abs(((int16_t *)ref)[i] - ((int16_t *)src)[i]);
    }
    src += systride;
    ref += dystride;
  }
  return ret + (1 << OD_COEFF_SHIFT >> 1) >> OD_COEFF_SHIFT;
}

int32_t od_mc_compute_sad16_4x4_c(const unsigned char *src, int systride,
 const unsigned char *ref, int dystride) {
  return od_mc_compute_sad16_c(src, systride, ref, dystride, 4, 4);
}

int32_t od_mc_compute_sad16_8x8_c(const unsigned char *src, int systride,
 const unsigned char *ref, int dystride) {
  return od_mc_compute_sad16_c(src, systride, ref, dystride, 8, 8);
}

int32_t od_mc_compute_sad16_16x16_c(const unsigned char *src, int systride,
 const unsigned char *ref, int dystride) {
  return od_mc_compute_sad16_c(src, systride, ref, dystride, 16, 16);
}

int32_t od_mc_compute_sad16_32x32_c(const unsigned char *src, int systride,
 const unsigned char *ref, int dystride) {
  return od_mc_compute_sad16_c(src, systride, ref, dystride, 32, 32);
}

static void od_mc_hadamard_1d(int32_t *diff,
 int n, int stride0, int stride1){
  int i;
  if (n == 4) {
    /*4x4 is the base case as it's our smallest block size.
      Perform the low-level 2x2 butterflies.*/
    int32_t a;
    int32_t b;
    int32_t c;
    int32_t d;
    for (i = 0; i < 4; i++) {
      a = diff[0*stride1] + diff[1*stride1];
      b = diff[0*stride1] - diff[1*stride1];
      c = diff[2*stride1] + diff[3*stride1];
      d = diff[2*stride1] - diff[3*stride1];
      diff[0*stride1] = a + c;
      diff[2*stride1] = a - c;
      diff[1*stride1] = b + d;
      diff[3*stride1] = b - d;
      diff += stride0;
    }
  }
  else {
    /*Recursive case for 8x8, 16x16, etc.
      Subdivide then combine.*/
    int n2;
    int j;
    int32_t *ptr0;
    int32_t *ptr1;
    n2 = n >> 1;
    od_mc_hadamard_1d(diff, n2, stride0, stride1);
    od_mc_hadamard_1d(diff + n2*stride0, n2, stride0, stride1);
    od_mc_hadamard_1d(diff + n2*stride1, n2, stride0, stride1);
    od_mc_hadamard_1d(diff + n2*stride0 + n2*stride1, n2, stride0, stride1);
    ptr0 = diff;
    ptr1 = diff + stride1*n2;
    for (i = 0; i < n; i++){
      for (j = 0; j < n2*stride1; j+=stride1){
        int32_t temp;
        temp = ptr0[j] - ptr1[j];
        ptr0[j] += ptr1[j];
        ptr1[j] = temp;
      }
      ptr0 += stride0;
      ptr1 += stride0;
    }
  }
}

static int od_mc_compute_satd8(int32_t *work, int ln,
 const unsigned char *src, int systride,
 const unsigned char *ref, int rystride) {
  int n;
  int i;
  int j;
  int32_t satd;
  int32_t *p;
  n = 1 << ln;
  p = work;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      p[j] = src[j] - ref[j];
    }
    p += n;
    src += systride;
    ref += rystride;
  }
  /*Horizontal transform.*/
  od_mc_hadamard_1d(work, n, n, 1);
  /*Vertical transform.*/
  od_mc_hadamard_1d(work, n, 1, n);
  satd = 0;
  for (i = 0; i < n*n; i++) satd += abs(work[i]);
  return satd + (1 << ln >> 1) >> ln;
}

static int od_mc_compute_satd16(int32_t *work, int ln,
 const unsigned char *src, int systride,
 const unsigned char *ref, int rystride) {
  int n;
  int i;
  int j;
  int32_t satd;
  int32_t *p;
  n = 1 << ln;
  p = work;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      p[j] = ((int16_t *)src)[j] - ((int16_t *)ref)[j];
    }
    p += n;
    src += systride;
    ref += rystride;
  }
  /*Horizontal transform.*/
  od_mc_hadamard_1d(work, n, n, 1);
  /*Vertical transform.*/
  od_mc_hadamard_1d(work, n, 1, n);
  satd = 0;
  for (i = 0; i < n*n; i++) satd += abs(work[i]);
  return satd + (1 << ln + OD_COEFF_SHIFT >> 1) >> ln + OD_COEFF_SHIFT;
}

int32_t od_mc_compute_satd8_4x4_c(const unsigned char *src, int systride,
 const unsigned char *ref, int rystride) {
  int32_t work[4*4];
  return od_mc_compute_satd8(work, 2, src, systride, ref, rystride);
}

int32_t od_mc_compute_satd16_4x4_c(const unsigned char *src, int systride,
 const unsigned char *ref, int rystride) {
  int32_t work[4*4];
  return od_mc_compute_satd16(work, 2, src, systride, ref, rystride);
}

int32_t od_mc_compute_satd8_8x8_c(const unsigned char *src, int systride,
 const unsigned char *ref, int rystride) {
  int32_t work[8*8];
  return od_mc_compute_satd8(work, 3, src, systride, ref, rystride);
}

int32_t od_mc_compute_satd16_8x8_c(const unsigned char *src, int systride,
 const unsigned char *ref, int rystride) {
  int32_t work[8*8];
  return od_mc_compute_satd16(work, 3, src, systride, ref, rystride);
}

int32_t od_mc_compute_satd8_16x16_c(const unsigned char *src, int systride,
 const unsigned char *ref, int rystride) {
  int32_t work[16*16];
  return od_mc_compute_satd8(work, 4, src, systride, ref, rystride);
}

int32_t od_mc_compute_satd16_16x16_c(const unsigned char *src, int systride,
 const unsigned char *ref, int rystride) {
  int32_t work[16*16];
  return od_mc_compute_satd16(work, 4, src, systride, ref, rystride);
}

int32_t od_mc_compute_satd8_32x32_c(const unsigned char *src, int systride,
 const unsigned char *ref, int rystride) {
  int32_t work[32*32];
  return od_mc_compute_satd8(work, 5, src, systride, ref, rystride);
}

int32_t od_mc_compute_satd16_32x32_c(const unsigned char *src, int systride,
 const unsigned char *ref, int rystride) {
  int32_t work[32*32];
  return od_mc_compute_satd16(work, 5, src, systride, ref, rystride);
}

int32_t od_mc_compute_satd8_64x64_c(const unsigned char *src, int systride,
 const unsigned char *ref, int rystride) {
  int32_t work[64*64];
  return od_mc_compute_satd8(work, 6, src, systride, ref, rystride);
}

int32_t od_mc_compute_satd16_64x64_c(const unsigned char *src, int systride,
 const unsigned char *ref, int rystride) {
  int32_t work[64*64];
  return od_mc_compute_satd16(work, 6, src, systride, ref, rystride);
}

/*Computes the SAD of the input image against the given predictor.*/
static int32_t od_enc_sad(od_enc_ctx *enc, const unsigned char *p,
 int pystride, int pxstride, int pli, int x, int y, int log_blk_sz) {
  od_state *state;
  od_img_plane *iplane;
  unsigned char *src;
  int clipx;
  int clipy;
  int clipw;
  int cliph;
  int w;
  int h;
  int xdec = od_plane_info_tab[enc->input_img.pixel_format].xdec[pli];
  int ydec = od_plane_info_tab[enc->input_img.pixel_format].ydec[pli];
  state = &enc->state;
  iplane = enc->input_img.planes + pli;
  /*Compute the block dimensions in the target image plane.*/
  x >>= xdec;
  y >>= ydec;
  w = 1 << (log_blk_sz - xdec);
  h = 1 << (log_blk_sz - ydec);
  /*Clip the block against the active picture region.*/
  clipx = -x;
  if (clipx > 0) {
    w -= clipx;
    p += clipx*pxstride;
    x += clipx;
  }
  clipy = -y;
  if (clipy > 0) {
    h -= clipy;
    p += clipy*pystride;
    y += clipy;
  }
  clipw = ((state->info.pic_width + (1 << xdec) - 1) >> xdec)
   - x;
  w = OD_MINI(w, clipw);
  cliph = ((state->info.pic_height + (1 << ydec) - 1) >> ydec)
   - y;
  h = OD_MINI(h, cliph);
  /*OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
   "[%i, %i]x[%i, %i]", x, y, w, h));*/
  /*Compute the SAD.*/
  src = iplane->data + y*iplane->ystride + x*iplane->xstride;
  if (w == 4 && h == 4) {
    return (*enc->opt_vtbl.mc_compute_sad_4x4)(src, iplane->ystride,
     p, pystride);
  }
  if (w == 8 && h == 8) {
    return (*enc->opt_vtbl.mc_compute_sad_8x8)(src, iplane->ystride,
     p, pystride);
  }
  if (w == 16 && h == 16) {
    return (*enc->opt_vtbl.mc_compute_sad_16x16)(src, iplane->ystride,
     p, pystride);
  }
  if (w == 32 && h == 32) {
    return (*enc->opt_vtbl.mc_compute_sad_32x32)(src, iplane->ystride,
     p, pystride);
  }
  /*Default C implementation.*/
  if (pxstride == 1) {
    return od_mc_compute_sad8_c(src, iplane->ystride, p, pystride, w, h);
  }
  return od_mc_compute_sad16_c(src, iplane->ystride, p, pystride, w, h);
}

/*Computes the SATD of the input image block against the given predictor.*/
static int32_t od_enc_satd(od_enc_ctx *enc, const unsigned char *p,
 int pystride, int pxstride, int pli, int x, int y, int log_blk_sz) {
  od_state *state;
  od_img_plane *iplane;
  unsigned char *src;
  int clipx;
  int clipy;
  int clipw;
  int cliph;
  int w;
  int h;
  int xdec = od_plane_info_tab[enc->input_img.pixel_format].xdec[pli];
  int ydec = od_plane_info_tab[enc->input_img.pixel_format].ydec[pli];
  state = &enc->state;
  iplane = enc->input_img.planes + pli;
  /*Compute the block dimensions in the target image plane.*/
  x >>= xdec;
  y >>= ydec;
  w = 1 << (log_blk_sz - xdec);
  h = 1 << (log_blk_sz - ydec);
  /*Clip the block against the active picture region.*/
  clipx = -x;
  if (clipx > 0) {
    w -= clipx;
    p += clipx*pxstride;
    x += clipx;
  }
  clipy = -y;
  if (clipy > 0) {
    h -= clipy;
    p += clipy*pystride;
    y += clipy;
  }
  clipw = ((state->info.pic_width + (1 << xdec) - 1) >> xdec)
   - x;
  w = OD_MINI(w, clipw);
  cliph = ((state->info.pic_height + (1 << ydec) - 1) >> ydec)
   - y;
  h = OD_MINI(h, cliph);
  /*OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
   "[%i, %i]x[%i, %i]", x, y, w, h));*/
  /*Compute the SATD.*/
  src = iplane->data + y*iplane->ystride + x*iplane->xstride;
  if (w == 4 && h == 4) {
    return (*enc->opt_vtbl.mc_compute_satd_4x4)(src, iplane->ystride,
     p, pystride);
  }
  if (w == 8 && h == 8) {
    return (*enc->opt_vtbl.mc_compute_satd_8x8)(src, iplane->ystride,
     p, pystride);
  }
  if (w == 16 && h == 16) {
    return (*enc->opt_vtbl.mc_compute_satd_16x16)(src, iplane->ystride,
     p, pystride);
  }
  if (w == 32 && h == 32) {
    return (*enc->opt_vtbl.mc_compute_satd_32x32)(src, iplane->ystride,
     p, pystride);
  }
  if (w == 64 && h == 64) {
    return (*enc->opt_vtbl.mc_compute_satd_64x64)(src, iplane->ystride,
     p, pystride);
  }
  /*If not square SATD (on boundary always), run sad for now.
    TODO: Try padding 0's for undefined area of difference image,
    then apply square SATD.*/
  if(pxstride == 1){
    return od_mc_compute_sad8_c(src, iplane->ystride, p, pystride, w, h);
  }
  return od_mc_compute_sad16_c(src, iplane->ystride, p, pystride, w, h);
}

static int od_mv_est_init_impl(od_mv_est_ctx *est, od_enc_ctx *enc) {
  int nhmvbs;
  int nvmvbs;
  int log_mvb_sz;
  int vx;
  int vy;
  if (OD_UNLIKELY(!est)) {
    return OD_EFAULT;
  }
  OD_CLEAR(est, 1);
  est->enc = enc;
  nhmvbs = enc->state.nhmvbs;
  nvmvbs = enc->state.nvmvbs;
  for (log_mvb_sz = 0; log_mvb_sz < OD_LOG_MVB_DELTA0 ; log_mvb_sz++) {
    est->sad_cache[log_mvb_sz] = (od_sad4 **)od_malloc_2d(nvmvbs >> log_mvb_sz,
     nhmvbs >> log_mvb_sz, sizeof(est->sad_cache[log_mvb_sz][0][0]));
    if (OD_UNLIKELY(!est->sad_cache[log_mvb_sz])) return OD_EFAULT;
  }
  est->mvs = (od_mv_node **)od_calloc_2d(nvmvbs + 1, nhmvbs + 1,
   sizeof(est->mvs[0][0]));
  if (OD_UNLIKELY(!est->mvs)) {
    return OD_EFAULT;
  }
  est->refine_grid = (od_mv_grid_pt **)od_malloc_2d(nvmvbs + 1, nhmvbs + 1,
   sizeof(est->refine_grid[0][0]));
  if (OD_UNLIKELY(!est->refine_grid)) {
    return OD_EFAULT;
  }
  est->dp_nodes = (od_mv_dp_node *)malloc(
   sizeof(od_mv_dp_node)*(OD_MAXI(nhmvbs, nvmvbs) + 1));
  if (OD_UNLIKELY(!est->dp_nodes)) {
    return OD_EFAULT;
  }
  est->row_counts =
   (unsigned *)malloc(sizeof(*est->row_counts)*(nvmvbs + 1));
  if (OD_UNLIKELY(!est->row_counts)) {
    return OD_EFAULT;
  }
  est->col_counts =
   (unsigned *)malloc(sizeof(*est->col_counts)*(nhmvbs + 1));
  if (OD_UNLIKELY(!est->col_counts)) {
    return OD_EFAULT;
  }
  for (vy = 0; vy <= nvmvbs; vy++) {
    for (vx = 0; vx <= nhmvbs; vx++) {
      est->mvs[vy][vx].vx = vx;
      est->mvs[vy][vx].vy = vy;
      est->mvs[vy][vx].heapi = -1;
      enc->state.mv_grid[vy][vx].valid = 1;
    }
  }
  est->dec_heap = (od_mv_node **)malloc(
   sizeof(*est->dec_heap)*(nvmvbs + 1)*(nhmvbs + 1));
  if (OD_UNLIKELY(!est->dec_heap)) {
    return OD_EFAULT;
  }
  /*Set to UCHAR_MAX so that od_mv_est_clear_hit_cache initializes hit_cache.*/
  est->hit_bit = UCHAR_MAX;
  est->mv_res_min = 0;
  est->flags = OD_MC_USE_CHROMA;
  return OD_SUCCESS;
}

static void od_mv_est_clear(od_mv_est_ctx *est) {
  int log_mvb_sz;
  free(est->dec_heap);
  free(est->col_counts);
  free(est->row_counts);
  free(est->dp_nodes);
  od_free_2d(est->refine_grid);
  od_free_2d(est->mvs);
  for (log_mvb_sz = OD_LOG_MVB_DELTA0; log_mvb_sz-- > 0; ) {
    od_free_2d(est->sad_cache[log_mvb_sz]);
  }
}

static int od_mv_est_init(od_mv_est_ctx *est, od_enc_ctx *enc) {
  int ret;
  ret = od_mv_est_init_impl(est, enc);
  if (OD_UNLIKELY(ret < 0)) {
    od_mv_est_clear(est);
  }
  return ret;
}

/*STAGE 1: INITIAL MV ESTIMATES (via EPZS^2).*/

/*The amount to right shift the minimum error by when inflating it for
   computing the second maximum SAD threshold.*/
#define OD_MC_THRESH2_SCALE_BITS (3)

/*The vector offsets in the X direction for each search site in the various
   patterns.*/
static const int OD_SITE_DX[13] = {
  -1, 0, 1, -1, 0, 1, -1, 0, 1, -2, 0, 2, 0
};
/*The vector offsets in the Y direction for each search site in the various
   patterns.*/
static const int OD_SITE_DY[13] = {
  -1, -1, -1, 0, 0, 0, 1, 1, 1, 0, -2, 0, 2
};

/*Set to do the initial search with a square search pattern instead of the 3-D
   predict hexagon state machine.
  We still need the square pattern tables below even if this is unset, since
   we may use them for MV refinement in stages 3 and 4*/
#undef OD_USE_SQUARE_SEARCH

/*The number of sites to search of each boundary condition in the square
   pattern.
  Bit flags for the boundary conditions are as follows:
  1: -32 == dx
  2: dx == 31
  4: -32 == dy
  8: dy == 31*/
static const int OD_SQUARE_NSITES[11] = { 8, 5, 5, 0, 5, 3, 3, 0, 5, 3, 3 };
/*The list of sites to search for each boundary condition in the square
   pattern.*/
static const od_pattern OD_SQUARE_SITES[11] = {
  /* -32 < dx < 31,   -32 < dy < 31*/
  { 0, 1, 2, 3, 5, 6, 7, 8 },
  /*-32 == dx,        -32 < dy < 31*/
  { 1, 2, 5, 7, 8 },
  /*       dx == 31,  -32 < dy < 31*/
  { 0, 1, 3, 6, 7 },
  /*-32 == dx == 31,  -32 < dy < 31*/
  { -1 },
  /* -32 < dx < 31,  -32 == dy*/
  { 3, 5, 6, 7, 8 },
  /*-32 == dx,       -32 == dy*/
  { 5, 7, 8 },
  /*       dx == 31, -32 == dy*/
  { 3, 6, 7 },
  /*-32 == dx == 31, -32 == dy*/
  { -1 },
  /* -32 < dx < 31,         dy == 31*/
  { 0, 1, 2, 3, 5 },
  /*-32 == dx,              dy == 31*/
  { 1, 2, 5 },
  /*       dx == 31,        dy == 31*/
  { 0, 1, 3 }
};

/*The number of sites to search of each boundary condition in the diamond
   pattern.
  Bit flags for the boundary conditions are as follows:
  1: -32 == dx
  2: dx == 31
  4: -32 == dy
  8: dy == 31*/
static const int OD_DIAMOND_NSITES[11] = { 4, 3, 3, 0, 3, 2, 2, 0, 3, 2, 2 };
/*The list of sites to search for each boundary condition in the square
   pattern.*/
static const od_pattern OD_DIAMOND_SITES[11] = {
  /* -32 < dx < 31,   -32 < dy < 31*/
  { 1, 3, 5, 7 },
  /*-32 == dx,        -32 < dy < 31*/
  { 1, 5, 7 },
  /*       dx == 31,  -32 < dy < 31*/
  { 1, 3, 7 },
  /*-32 == dx == 31,  -32 < dy < 31*/
  { -1 },
  /* -32 < dx < 31,  -32 == dy*/
  { 3, 5, 7 },
  /*-32 == dx,       -32 == dy*/
  { 5, 7 },
  /*       dx == 31, -32 == dy*/
  { 3, 7 },
  /*-32 == dx == 31, -32 == dy*/
  { -1 },
  /* -32 < dx < 31,         dy == 31*/
  { 1, 3, 5 },
  /*-32 == dx,              dy == 31*/
  { 1, 5 },
  /*       dx == 31,        dy == 31*/
  { 1, 3 }
};

#if !defined(OD_USE_SQUARE_SEARCH)

/*The number of sites to search of each boundary condition in the horizontal
   hex pattern.
  Bit flags for the boundary conditions are as follows:
  1: -32 == dx
  2: dx == 31
  4: -32 == dy
  8: dy == 31*/
static const int OD_HHEX_NSITES[11] = { 6, 3, 3, 0, 4, 2, 2, 0, 4, 2, 2 };
/*The list of sites to search for each boundary condition in the horizontal
   hex pattern.*/
static const od_pattern OD_HHEX_SITES[11] = {
  /* -32 < dx < 31,   -32 < dy < 31*/
  { 0, 2, 6, 8, 9, 11 },
  /*-32 == dx,        -32 < dy < 31*/
  { 2, 8, 11 },
  /*       dx == 31,  -32 < dy < 31*/
  { 0, 6, 9 },
  /*-32 == dx == 31,  -32 < dy < 31*/
  { -1 },
  /* -32 < dx < 31,  -32 == dy*/
  { 6, 8, 9, 11 },
  /*-32 == dx,       -32 == dy*/
  { 8, 11 },
  /*       dx == 31, -32 == dy*/
  { 6, 9 },
  /*-32 == dx == 31, -32 == dy*/
  { -1 },
  /* -32 < dx < 31,         dy == 31*/
  { 0, 2, 9, 11 },
  /*-32 == dx,              dy == 31*/
  { 2, 11 },
  /*     dx == 31,          dy == 31*/
  { 0, 9 }
};

/*The number of sites to search of each boundary condition in the vertical hex
   pattern.
  Bit flags for the boundary conditions are as follows:
  1: -32 == dx
  2: dx == 31
  4: -32 == dy
  8: dy == 31*/
static const int OD_VHEX_NSITES[11] = { 6, 4, 4, 0, 3, 2, 2, 0, 3, 2, 2 };
/*The list of sites to search for each boundary condition in the vertical hex
   pattern.*/
static const od_pattern OD_VHEX_SITES[11] = {
  /* -32 < dx < 31,   -32 < dy < 31*/
  { 0, 2, 6, 8, 10, 12 },
  /*-32 == dx,        -32 < dy < 31*/
  { 2, 8, 10, 12 },
  /*       dx == 31,  -32 < dy < 31*/
  { 0, 6, 10, 12 },
  /*-32 == dx == 31,  -32 < dy < 31*/
  { -1 },
  /* -32 < dx < 31,  -32 == dy*/
  { 6, 8, 12 },
  /*-32 == dx,       -32 == dy*/
  { 8, 12 },
  /*       dx == 31, -32 == dy*/
  { 6, 12 },
  /*-32 == dx == 31, -32 == dy*/
  { -1 },
  /* -32 < dx < 31,         dy == 31*/
  { 0, 2, 10 },
  /*-32 == dx,              dy == 31*/
  { 2, 10 },
  /*       dx == 31,        dy == 31*/
  { 0, 10 }
};

#endif

#if defined(OD_USE_SQUARE_SEARCH)
/*The search state indicating we found a local minimum.*/
# define OD_SEARCH_STATE_DONE (1)

/*The number of sites in the pattern to use for each state.*/
static const int *const OD_SEARCH_NSITES[1] = {
  OD_SQUARE_NSITES
};

/*The sites in the pattern to use for each state.*/
static const od_pattern *const OD_SEARCH_SITES[1] = {
  OD_SQUARE_SITES
};

/*The successor state given the current state and the terminating site.*/
static const int OD_SEARCH_STATES[1][13] = {
  /*Just use a square pattern for the whole search.*/
  { 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, -1, -1, -1 }
};
#else
/*The state machine used to select which patterns to use in the gradient
   descent BMA search.
  These are based on the procedure described in
  @ARTICLE{TP06,
    author="Tsung-Han Tsai and Yu-Nan Pan",
    title="A Novel 3-D Predict Hexagon Search Algorithm for Fast Block Motion
     Estimation on H.264 Video Coding",
    journal="{IEEE} Transactions on Circuits and Systems for Video Technology",
    volume=16,
    number=12,
    pages="1542--1549",
    month=Dec,
    year=2006
  }*/

/*The search state indicating we found a local minimum.*/
# define OD_SEARCH_STATE_DONE (6)

/*The number of sites in the pattern to use for each state.*/
static const int *const OD_SEARCH_NSITES[6] = {
  OD_DIAMOND_NSITES,
  OD_DIAMOND_NSITES,
  OD_DIAMOND_NSITES,
  OD_HHEX_NSITES,
  OD_VHEX_NSITES,
  OD_DIAMOND_NSITES
};

/*The sites in the pattern to use for each state.*/
static const od_pattern *const OD_SEARCH_SITES[6] = {
  OD_DIAMOND_SITES,
  OD_DIAMOND_SITES,
  OD_DIAMOND_SITES,
  OD_HHEX_SITES,
  OD_VHEX_SITES,
  OD_DIAMOND_SITES
};

/*The successor state given the current state and the terminating site.*/
static const int OD_SEARCH_STATES[6][13] = {
  /*Start with a small diamond in the first step.*/
  { -1, 2, -1, 1, 6, 1, -1, 2, -1, -1, -1, -1, -1 },
  /*Use a small diamond for the second step, too, but remember if we took a
     horizontal step in the first step...*/
  { -1, 3, -1, 3, 6, 3, -1, 3, -1, -1, -1, -1, -1 },
  /*...or a vertical one.*/
  { -1, 4, -1, 4, 6, 4, -1, 4, -1, -1, -1, -1, -1 },
  /*Then switch to a horizontal hex pattern for all remaining steps...*/
  { 3, -1, 3, -1, 5, -1, 3, -1, 3, 3, -1, 3, -1 },
  /*...or a vertical hex pattern, depending.*/
  { 4, -1, 4, -1, 5, -1, 4, -1, 4, -1, 4, -1, 4 },
  /*And revert back to a small diamond for the last step.*/
  { -1, 6, -1, 6, 6, 6, -1, 6, -1, -1, -1, -1, -1 }
};
#endif

/*Clear the cache of motion vectors we've examined.*/
static void od_mv_est_clear_hit_cache(od_mv_est_ctx *est) {
  if (++est->hit_bit == UCHAR_MAX + 1) {
    memset(est->hit_cache, 0, sizeof(est->hit_cache));
    est->hit_bit = 1;
  }
}

/*Test if a motion vector has been examined.*/
static int od_mv_est_is_hit(od_mv_est_ctx *est, int mvx, int mvy) {
  return est->hit_cache[mvy + OD_MC_SEARCH_RANGE][mvx + OD_MC_SEARCH_RANGE] 
   == est->hit_bit;
}

/*Mark a motion vector examined.*/
static void od_mv_est_set_hit(od_mv_est_ctx *est, int mvx, int mvy) {
  est->hit_cache[mvy + OD_MC_SEARCH_RANGE][mvx + OD_MC_SEARCH_RANGE] = 
   (unsigned char)est->hit_bit;
}

/*Estimated rate (in units of OD_BITRES) of the >=3 part of a MV component of a
   given length.
  There is no theoretical basis for the current values, they were
   supposed to estimate the rate of a complete MV component (not just the >=3
   part) encoded using a Huffman code.
  TODO: This should be replaced with more realistic estimates (so far, attempts
  to do so have resulted in worse quality).*/
static const int OD_MV_GE3_EST_RATE[256] = {
    8,  32,  32,  48,  48,  48,  48,  64,
   64,  64,  64,  64,  64,  64,  64,  80,
   80,  80,  80,  80,  80,  80,  80,  80,
   80,  80,  80,  80,  80,  80,  80,  96,
   96,  96,  96,  96,  96,  96,  96,  96,
   96,  96,  96,  96,  96,  96,  96,  96,
   96,  96,  96,  96,  96,  96,  96,  96,
   96,  96,  96,  96,  96,  96,  96, 112,
  112, 112, 112, 112, 112, 112, 112, 112,
  112, 112, 112, 112, 112, 112, 112, 112,
  112, 112, 112, 112, 112, 112, 112, 112,
  112, 112, 112, 112, 112, 112, 112, 112,
  112, 112, 112, 112, 112, 112, 112, 112,
  112, 112, 112, 112, 112, 112, 112, 112,
  112, 112, 112, 112, 112, 112, 112, 112,
  112, 112, 112, 112, 112, 112, 112, 128,
  128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 144
};

/*Estimate the number of bits that will be used to encode the given MV and its
   predictor.*/
static int od_mv_est_cand_bits(od_mv_est_ctx *est, int equal_mvs,
 int dx, int dy, int predx, int predy, int ref, int ref_pred) {
  int ox;
  int oy;
  int id;
  int cost;
  int sign_cost;
  cost = 0;
  sign_cost = 1 << OD_BITRES;
  ox = dx - predx;
  oy = dy - predy;
  id = OD_MINI(abs(oy), 3)*4 + OD_MINI(abs(ox), 3);
  cost += ((ox != 0) + (oy != 0))*sign_cost;
  cost += est->mv_small_rate_est[equal_mvs][id];
  if (abs(ox) >= 3) {
    cost += OD_MV_GE3_EST_RATE[OD_MINI(abs(ox) - 3, 255)];
  }
  if (abs(oy) >= 3) {
    cost += OD_MV_GE3_EST_RATE[OD_MINI(abs(oy) - 3, 255)];
  }
  if (ref_pred != ref) {
    cost += 2 << OD_BITRES;
  }
  return cost;
}

/*Estimate the number of bits that will be used to encode the given MV grid
  point, including reference index.*/
static int od_mv_est_bits(od_mv_est_ctx *est, int vx, int vy, int mv_res) {
  od_state *state;
  int pred[2];
  int equal_mvs;
  int mv_rate;
  int ref_pred;
  od_mv_grid_pt *mvg;
  state = &est->enc->state;
  mvg = state->mv_grid[vy] + vx;
  equal_mvs = od_state_get_predictor(state, pred, vx, vy,
   OD_MC_LEVEL[vy & OD_MVB_MASK][vx & OD_MVB_MASK], mv_res, mvg->ref);
  ref_pred = od_mc_get_ref_predictor(state,  vx, vy,
   OD_MC_LEVEL[vy & OD_MVB_MASK][vx & OD_MVB_MASK]);
  mv_rate = od_mv_est_cand_bits(est, equal_mvs,
   mvg->mv[0] >> mv_res, mvg->mv[1] >> mv_res, pred[0], pred[1],
   mvg->ref, ref_pred);
  return mv_rate;
}

#if defined(OD_LOGGING_ENABLED)
# define OD_MV_EST_LOG_PRED(est, vx, vy, mv_res) \
   od_mv_est_log_pred(est, vx, vy, mv_res)
static void od_mv_est_log_pred(od_mv_est_ctx *est, int vx, int vy,
 int mv_res) {
  od_state *state;
  int pred[2];
  int equal_mvs;
  int ref_pred;
  od_mv_grid_pt *mvg;
  state = &est->enc->state;
  mvg = state->mv_grid[vy] + vx;
  equal_mvs = od_state_get_predictor(state, pred, vx, vy,
   OD_MC_LEVEL[vy & OD_MVB_MASK][vx & OD_MVB_MASK], mv_res, mvg->ref);
  ref_pred = od_mc_get_ref_predictor(state,  vx, vy,
   OD_MC_LEVEL[vy & OD_MVB_MASK][vx & OD_MVB_MASK]);
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
   "Predictor was: (%i, %i)@%i, %i equal_mvs",
   pred[0], pred[1], ref_pred, equal_mvs));
}
#else
# define OD_MV_EST_LOG_PRED(est, vx, vy, mv_res)
#endif

/*Computes the SAD of a whole-pel BMA block with the given parameters.*/
static int32_t od_mv_est_bma_sad(od_mv_est_ctx *est,
 int ref, int bx, int by, int mvx, int mvy, int log_mvb_sz) {
  daala_enc_ctx *enc;
  od_state *state;
  od_img_plane *iplane;
  int32_t ret;
  int refi;
  int dx;
  int dy;
  int xdec;
  int ydec;

  enc = est->enc;
  state = &enc->state;
  refi = state->ref_imgi[ref];
  xdec = od_plane_info_tab[state->ref_imgs[refi].pixel_format].xdec[0];
  ydec  = od_plane_info_tab[state->ref_imgs[refi].pixel_format].ydec[0];
  iplane = state->ref_imgs[refi].planes + 0;
  OD_ASSERT(xdec == 0 && ydec == 0);
  dx = bx + mvx;
  dy = by + mvy;
  ret = od_enc_sad(est->enc,
   iplane->data + dy*iplane->ystride + dx*iplane->xstride,
   iplane->ystride, iplane->xstride,
   0, bx, by, log_mvb_sz + OD_LOG_MVBSIZE_MIN);
  if (est->flags & OD_MC_USE_CHROMA) {
    int pli;
    unsigned char *ref_img;
    for (pli = 1; pli < enc->input_img.nplanes; pli++) {
      xdec = od_plane_info_tab[state->ref_imgs[refi].pixel_format].xdec[pli];
      ydec = od_plane_info_tab[state->ref_imgs[refi].pixel_format].ydec[pli];
      iplane = state->ref_imgs[refi].planes + pli;
      OD_ASSERT(((bx + (1 << iplane->xdec) - 1) & ~((1 << xdec) - 1))
       == bx);
      OD_ASSERT(((by + (1 << iplane->ydec) - 1) & ~((1 << ydec) - 1))
       == by);
      /*If the input chroma plane is sub-sampled, then the candidate block with
         subpel position for BMA search is interpolated at block level.*/
      ref_img = iplane->data
       + (by >> ydec)*iplane->ystride
       + (bx >> xdec)*iplane->xstride;
      (*state->opt_vtbl.mc_predict1fmv)
       (state, state->mc_buf[4], ref_img, iplane->ystride,
       mvx << (3 - xdec), mvy << (3 - ydec),
       log_mvb_sz + OD_LOG_MVBSIZE_MIN - xdec,
       log_mvb_sz + OD_LOG_MVBSIZE_MIN - ydec);
      /*Then, calculate SAD between a target block and the subpel interpolated
         MC block.*/
      ret += od_enc_sad(est->enc, state->mc_buf[4],
       iplane->xstride << (log_mvb_sz + OD_LOG_MVBSIZE_MIN - xdec),
       iplane->xstride,
       pli, bx, by, log_mvb_sz + OD_LOG_MVBSIZE_MIN) >> OD_MC_CHROMA_SCALE;
    }
  }
  return ret;
}

/*Computes the SAD of a block with the given parameters.*/
static int32_t od_mv_est_sad(od_mv_est_ctx *est,
 int vx, int vy, int oc, int s, int log_mvb_sz) {
  od_state *state;
  int32_t ret;
  int xstride;
  state = &est->enc->state;
  xstride = est->enc->input_img.planes[0].xstride;
  od_state_pred_block_from_setup(state, state->mc_buf[4], OD_MVBSIZE_MAX,
   0, vx, vy, oc, s, log_mvb_sz);
  ret = est->compute_distortion(est->enc, state->mc_buf[4],
   OD_MVBSIZE_MAX*xstride, xstride,
   0, vx << OD_LOG_MVBSIZE_MIN, vy << OD_LOG_MVBSIZE_MIN,
   log_mvb_sz + OD_LOG_MVBSIZE_MIN);
  if (est->flags & OD_MC_USE_CHROMA) {
    int pli;
    for (pli = 1; pli < est->enc->input_img.nplanes; pli++) {
      od_state_pred_block_from_setup(state, state->mc_buf[4], OD_MVBSIZE_MAX,
       pli, vx, vy, oc, s, log_mvb_sz);
      ret += est->compute_distortion(est->enc, state->mc_buf[4],
       OD_MVBSIZE_MAX*xstride, xstride,
       pli, vx << OD_LOG_MVBSIZE_MIN, vy << OD_LOG_MVBSIZE_MIN,
       log_mvb_sz + OD_LOG_MVBSIZE_MIN) >> OD_MC_CHROMA_SCALE;
    }
  }
  return ret;
}

#if defined(OD_ENABLE_ASSERTIONS) || defined(OD_LOGGING_ENABLED)
/*Checks to make sure our current mv_rate and sad values are correct.
  This is used for debugging only.*/
void od_mv_est_check_rd_block_state(od_mv_est_ctx *est,
 int vx, int vy, int log_mvb_sz) {
  od_state *state;
  int half_mvb_sz;
  state = &est->enc->state;
  half_mvb_sz = 1 << log_mvb_sz >> 1;
  if (log_mvb_sz > 0
   && state->mv_grid[vy + half_mvb_sz][vx + half_mvb_sz].valid) {
    od_mv_est_check_rd_block_state(est, vx, vy, log_mvb_sz - 1);
    od_mv_est_check_rd_block_state(est,
     vx + half_mvb_sz, vy, log_mvb_sz - 1);
    od_mv_est_check_rd_block_state(est,
     vx, vy + half_mvb_sz, log_mvb_sz - 1);
    od_mv_est_check_rd_block_state(est,
     vx + half_mvb_sz, vy + half_mvb_sz, log_mvb_sz - 1);
  }
  else {
    od_mv_node *block;
    int32_t sad;
    int oc;
    int s;
    block = est->mvs[vy] + vx;
    if (block->log_mvb_sz != log_mvb_sz) {
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_WARN,
       "Failure at node (%i, %i): log_mvb_sz should be %i (is %i)",
       vx, vy, log_mvb_sz, block->log_mvb_sz));
      OD_ASSERT(block->log_mvb_sz == log_mvb_sz);
    }
    if (log_mvb_sz < OD_LOG_MVB_DELTA0) {
      int mask;
      int s1vx;
      int s1vy;
      int s3vx;
      int s3vy;
      mask = (1 << (log_mvb_sz + 1)) - 1;
      oc = !!(vx & mask);
      if (vy & mask) oc = 3 - oc;
      if (block->oc != oc) {
        OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_WARN,
         "Failure at node (%i, %i): oc should be %i (is %i)",
         vx, vy, oc, block->oc));
        OD_ASSERT(block->oc == oc);
      }
      s1vx = vx + (OD_VERT_DX[(oc + 1) & 3] << log_mvb_sz);
      s1vy = vy + (OD_VERT_DY[(oc + 1) & 3] << log_mvb_sz);
      s3vx = vx + (OD_VERT_DX[(oc + 3) & 3] << log_mvb_sz);
      s3vy = vy + (OD_VERT_DY[(oc + 3) & 3] << log_mvb_sz);
      s = state->mv_grid[s1vy][s1vx].valid |
       state->mv_grid[s3vy][s3vx].valid << 1;
    }
    else {
      oc = 0;
      s = 3;
    }
    if (block->s != s) {
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_WARN,
       "Failure at node (%i, %i): s should be %i (is %i)",
       vx, vy, s, block->s));
      OD_ASSERT(block->s == s);
    }
    sad = od_mv_est_sad(est, vx, vy, oc, s, log_mvb_sz);
    if (block->sad != sad) {
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_WARN,
       "Failure at node (%i, %i): sad should be %i (is %i)",
       vx, vy, sad, block->sad));
      OD_ASSERT(block->sad == sad);
    }
  }
}

/*Checks to make sure our current mv_rate and sad values are correct.
  This is used for debugging only.*/
void od_mv_est_check_rd_state(od_mv_est_ctx *est, int mv_res) {
  od_state *state;
  int nhmvbs;
  int nvmvbs;
  int vx;
  int vy;
# if !defined(OD_ENABLE_ASSERTIONS)
  if (!od_logging_active(OD_LOG_MOTION_ESTIMATION, OD_LOG_WARN)) return;
# endif
  state = &est->enc->state;
  nhmvbs = state->nhmvbs;
  nvmvbs = state->nvmvbs;
  for (vy = 0; vy < nvmvbs; vy += OD_MVB_DELTA0) {
    for (vx = 0; vx < nhmvbs; vx += OD_MVB_DELTA0) {
      od_mv_est_check_rd_block_state(est, vx, vy, OD_LOG_MVB_DELTA0);
    }
  }
  for (vy = 0; vy < nvmvbs; vy++) {
    for (vx = 0; vx < nhmvbs; vx++) {
      od_mv_grid_pt *mvg;
      od_mv_node *mv;
      int mv_rate;
      int level;
      int mvb_sz;
      mvg = state->mv_grid[vy] + vx;
      if (!mvg->valid) continue;
      level = OD_MC_LEVEL[vy & OD_MVB_MASK][vx & OD_MVB_MASK];
      mvb_sz = 1 << ((OD_MC_LEVEL_MAX - level) >> 1);
      if (level & 1) {
        OD_ASSERT(state->mv_grid[vy - mvb_sz][vx - mvb_sz].valid
         && state->mv_grid[vy - mvb_sz][vx + mvb_sz].valid
         && state->mv_grid[vy + mvb_sz][vx - mvb_sz].valid
         && state->mv_grid[vy + mvb_sz][vx + mvb_sz].valid);
      }
      else {
        OD_ASSERT((vy - mvb_sz < 0 || state->mv_grid[vy - mvb_sz][vx].valid)
         && (vx - mvb_sz < 0 || state->mv_grid[vy][vx - mvb_sz].valid)
         && (vy + mvb_sz > nvmvbs || state->mv_grid[vy + mvb_sz][vx].valid)
         && (vx + mvb_sz > nhmvbs || state->mv_grid[vy][vx + mvb_sz].valid));
      }
      mv = est->mvs[vy] + vx;
      mv_rate = od_mv_est_bits(est, vx, vy, mv_res);
      if (mv_rate != mv->mv_rate) {
        OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_WARN,
         "Failure at node (%i, %i): mv_rate should be %i (is %i)",
         vx, vy, mv_rate, mv->mv_rate));
        OD_MV_EST_LOG_PRED(est, vx, vy, mv_res);
        OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_WARN,
         "MV was: (%i, %i)@%i", mvg->mv[0] >> mv_res,
         mvg->mv[1] >> mv_res, mvg->ref));
        OD_ASSERT(mv_rate == mv->mv_rate);
      }
    }
  }
}
#endif

#if defined(OD_DUMP_IMAGES) && defined(OD_ANIMATE)
static const unsigned char OD_YCbCr_MVCAND[3] = { 210, 16, 214 };
#endif

static void od_mv_est_init_mv(od_mv_est_ctx *est, int ref, int vx, int vy,
 int must_update) {
  static const od_mv_node ZERO_NODE;
  od_state *state;
  od_mv_grid_pt *mvg;
  od_mv_node *mv;
  const od_mv_node *cneighbors[4];
  const od_mv_node *pneighbors[4];
  int32_t t2;
  int32_t best_sad;
  int32_t best_cost;
  int best_rate;
  int cands[6][2];
  int best_vec[2];
  int nhmvbs;
  int nvmvbs;
  int level;
  int log_mvb_sz;
  int mvb_sz;
  int bx;
  int by;
  int ncns;
  int equal_mvs;
  int bxmin;
  int bymin;
  int bxmax;
  int bymax;
  int mvxmin;
  int mvxmax;
  int mvymin;
  int mvymax;
  int candx;
  int candy;
  int ref_pred;
  int pred[2];
  int ci;
  int previous_cost;
#if defined(OD_DUMP_IMAGES) && defined(OD_ANIMATE)
  int animating;
  int x0;
  int y0;
#endif
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
   "Initial search for MV (%i, %i):", vx, vy));
  state = &est->enc->state;
  mv = est->mvs[vy] + vx;
  level = OD_MC_LEVEL[vy & OD_MVB_MASK][vx & OD_MVB_MASK];
  log_mvb_sz = (OD_MC_LEVEL_MAX - level) >> 1;
  mvb_sz = 1 << log_mvb_sz;
  mvg = state->mv_grid[vy] + vx;
#if defined(OD_DUMP_IMAGES) && defined(OD_ANIMATE)
  animating = daala_granule_basetime(state, state->cur_time) == ANI_FRAME;
  if (animating) {
    od_state_mc_predict(state, state->ref_imgs
     + state->ref_imgi[OD_FRAME_SELF]);
    od_encode_fill_vis(est->enc);
    x0 = (vx << (OD_LOG_MVBSIZE_MIN + 1)) + (OD_BUFFER_PADDING << 1);
    y0 = (vy << (OD_LOG_MVBSIZE_MIN + 1)) + (OD_BUFFER_PADDING << 1);
  }
#endif
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG, "Level %i (%ix%i block)",
   level, mvb_sz << OD_LOG_MVBSIZE_MIN, mvb_sz << OD_LOG_MVBSIZE_MIN));
  bx = vx << OD_LOG_MVBSIZE_MIN;
  by = vy << OD_LOG_MVBSIZE_MIN;
  if (mvg->valid) {
    /* Current rate may need to be updated because of a reference change
       during this BMA pass */
    mv->mv_rate = od_mv_est_bits(est, vx, vy, 2);
  }
  /*Each grid point can contribute to 4 blocks around it at the current block
     size, _except_ at the border of the frame, since we do not construct a
     prediction for blocks completely outside the frame.
    We do not try to take into account that further splitting might restrict
     the area this MV influences, since we do not know how far we will
     ultimately split.
    So this MV can affect a block of pixels bounded by
     [bxmin, bmax) x [bymin, bymax), and MVs within that area must point no
     farther than OD_UMV_CLAMP pixels outside of the frame.*/
  bxmin = OD_MAXI(bx - (mvb_sz << OD_LOG_MVBSIZE_MIN), 0);
  mvxmin = OD_MAXI(bxmin - OD_MC_SEARCH_RANGE, -OD_UMV_CLAMP) - bxmin;
  bxmax = OD_MINI(bx + (mvb_sz << OD_LOG_MVBSIZE_MIN), state->frame_width);
  mvxmax = OD_MINI(bxmax + OD_MC_SEARCH_RANGE - 1,
   state->frame_width + OD_UMV_CLAMP) - bxmax;
  bymin = OD_MAXI(by - (mvb_sz << OD_LOG_MVBSIZE_MIN), 0);
  mvymin = OD_MAXI(bymin - OD_MC_SEARCH_RANGE, -OD_UMV_CLAMP) - bymin;
  bymax = OD_MINI(by + (mvb_sz << OD_LOG_MVBSIZE_MIN), state->frame_height);
  mvymax = OD_MINI(bymax + OD_MC_SEARCH_RANGE - 1,
   state->frame_height + OD_UMV_CLAMP) - bymax;
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
   "(%i, %i): Search range: [%i, %i]x[%i, %i]",
   bx, by, mvxmin, mvymin, mvxmax, mvymax));
  /*For the purposes of the BMA search, we match against a block of size
     mvb_sz << LOG_MVBSIZE_MIN centered on the current grid point:

      (bx - (mvb_sz << (LOG_MVBSIZE_MIN - 1)),
       by - (mvb_sz << (LOG_MVBSIZE_MIN - 1)))
          |
          V
          +---------------+
          |               |
          |               |
          |               |
          |       +       |
          |    (bx,by)    |
          |               |
          |               |
          +---------------+

          <--------------->
      mvb_sz << LOG_MVBSIZE_MIN

    Adjust (bx, by) here to point to the upper-left corner of that block.*/
  bx -= mvb_sz << (OD_LOG_MVBSIZE_MIN - 1);
  by -= mvb_sz << (OD_LOG_MVBSIZE_MIN - 1);
  ncns = 4;
  nhmvbs = state->nhmvbs;
  nvmvbs = state->nvmvbs;
  equal_mvs = od_state_get_predictor(state, pred, vx, vy, level, 2, ref);
  candx = OD_CLAMPI(mvxmin, OD_DIV2_RE(pred[0]), mvxmax);
  candy = OD_CLAMPI(mvymin, OD_DIV2_RE(pred[1]), mvymax);
  ref_pred = od_mc_get_ref_predictor(state, vx, vy, level);
  /*Find additional candidates.*/
  if (level == 0) {
    if (vy >= mvb_sz) {
      cneighbors[0] = vx >= mvb_sz ?
       est->mvs[vy - mvb_sz] + vx - mvb_sz : &ZERO_NODE;
      cneighbors[1] = est->mvs[vy - mvb_sz] + vx;
      cneighbors[2] = vx + mvb_sz <= nhmvbs ?
       est->mvs[vy - mvb_sz] + vx + mvb_sz : &ZERO_NODE;
      pneighbors[0] = est->mvs[vy - mvb_sz] + vx;
    }
    else {
      cneighbors[2] = cneighbors[1] = cneighbors[0] = &ZERO_NODE;
      pneighbors[0] = &ZERO_NODE;
    }
    cneighbors[3] = vx >= mvb_sz ? est->mvs[vy] + vx - mvb_sz : &ZERO_NODE;
    pneighbors[1] = vx >= mvb_sz ? est->mvs[vy] + vx - mvb_sz : &ZERO_NODE;
    pneighbors[2] = vx + mvb_sz <= nhmvbs ?
     est->mvs[vy] + vx + mvb_sz : &ZERO_NODE;
    pneighbors[3] = vy + mvb_sz <= nvmvbs ?
     est->mvs[vy + OD_MVB_DELTA0] + vx : &ZERO_NODE;
  }
  else {
    if (level & 1) {
      pneighbors[0] = est->mvs[vy - mvb_sz] + vx - mvb_sz;
      pneighbors[1] = est->mvs[vy - mvb_sz] + vx + mvb_sz;
      pneighbors[2] = est->mvs[vy + mvb_sz] + vx - mvb_sz;
      pneighbors[3] = est->mvs[vy + mvb_sz] + vx + mvb_sz;
      memcpy(cneighbors, pneighbors, sizeof(cneighbors));
    }
    else {
      pneighbors[0] = vy >= mvb_sz ? est->mvs[vy - mvb_sz] + vx : &ZERO_NODE;
      pneighbors[1] = vx >= mvb_sz ? est->mvs[vy] + vx - mvb_sz : &ZERO_NODE;
      pneighbors[2] = vx + mvb_sz <= nhmvbs ?
       est->mvs[vy] + vx + mvb_sz : &ZERO_NODE;
      pneighbors[3] = vy + mvb_sz <= nvmvbs ?
       est->mvs[vy + mvb_sz] + vx : &ZERO_NODE;
      cneighbors[0] = pneighbors[0];
      cneighbors[1] = pneighbors[1];
      /*NOTE: Only one of these candidates can be excluded at a time, so
         there will always be at least 3.*/
      if (vx > 0 && vx + mvb_sz > ((vx + OD_MVB_MASK) & ~OD_MVB_MASK)) ncns--;
      else cneighbors[2] = pneighbors[2];
      if (vy > 0 && vy + mvb_sz > ((vy + OD_MVB_MASK) & ~OD_MVB_MASK)) ncns--;
      else cneighbors[ncns - 1] = pneighbors[3];
    }
  }
  /*Spatially correlated predictors (from the current frame):*/
  for (ci = 0; ci < ncns; ci++) {
    cands[ci][0] = OD_CLAMPI(mvxmin, cneighbors[ci]->bma_mvs[0][ref][0],
     mvxmax);
    cands[ci][1] = OD_CLAMPI(mvymin, cneighbors[ci]->bma_mvs[0][ref][1],
     mvymax);
  }
  od_mv_est_clear_hit_cache(est);
#if defined(OD_DUMP_IMAGES) && defined(OD_ANIMATE)
  if (animating) {
    od_img_draw_line(&est->enc->vis_img, x0, y0,
     x0 + (candx << 1), y0 + (candy << 1), OD_YCbCr_MVCAND);
  }
#endif
  best_sad = od_mv_est_bma_sad(est, ref, bx, by, candx, candy, log_mvb_sz);
  best_rate = od_mv_est_cand_bits(est, equal_mvs,
   candx << 1, candy << 1, pred[0], pred[1], ref, ref_pred);
  best_cost = (best_sad << OD_ERROR_SCALE) + best_rate*est->lambda;
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
   "Median predictor: (%i, %i)   Cost: %i", candx, candy, best_cost));
  od_mv_est_set_hit(est, candx, candy);
  best_vec[0] = candx;
  best_vec[1] = candy;
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
   "Threshold: %i", est->thresh1[log_mvb_sz]));
  if (best_sad > est->thresh1[log_mvb_sz]) {
    int32_t sad;
    int32_t cost;
    int rate;
    /*Compute the early termination threshold for set B.*/
    t2 = mv->bma_sad;
    for (ci = 0; ci < ncns; ci++) {
      int log_cnb_sz;
      int clevel;
      int cvx;
      int cvy;
      cvx = cneighbors[ci]->vx;
      cvy = cneighbors[ci]->vy;
      clevel = OD_MC_LEVEL[cvy & OD_MVB_MASK][cvx & OD_MVB_MASK];
      log_cnb_sz = (OD_MC_LEVEL_MAX - clevel) >> 1;
      t2 = OD_MINI(t2,
       cneighbors[ci]->bma_sad >> ((log_cnb_sz - log_mvb_sz) << 1));
    }
    t2 = t2 + (t2 >> OD_MC_THRESH2_SCALE_BITS) + est->thresh2_offs[log_mvb_sz];
    /*Constant velocity predictor:*/
    cands[ncns][0] =
     OD_CLAMPI(mvxmin, mv->bma_mvs[1][ref][0], mvxmax);
    cands[ncns][1] =
     OD_CLAMPI(mvymin, mv->bma_mvs[1][ref][1], mvymax);
    ncns++;
    /*Zero predictor.*/
    cands[ncns][0] = 0;
    cands[ncns][1] = 0;
    ncns++;
    /*Examine the candidates in Set B.*/
    for (ci = 0; ci < ncns; ci++) {
      candx = cands[ci][0];
      candy = cands[ci][1];
      if (od_mv_est_is_hit(est, candx, candy)) {
        OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
         "Set B predictor %i: (%i, %i) ...Skipping.", ci, candx, candy));
        continue;
      }
      od_mv_est_set_hit(est, candx, candy);
#if defined(OD_DUMP_IMAGES) && defined(OD_ANIMATE)
      if (animating) {
        od_img_draw_line(&est->enc->vis_img, x0, y0,
         x0 + (candx << 1), y0 + (candy << 1), OD_YCbCr_MVCAND);
      }
#endif
      sad = od_mv_est_bma_sad(est, ref, bx, by, candx, candy, log_mvb_sz);
      rate = od_mv_est_cand_bits(est, equal_mvs,
       candx << 1, candy << 1, pred[0], pred[1], ref, ref_pred);
      cost = (sad << OD_ERROR_SCALE) + rate*est->lambda;
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
       "Set B predictor %i: (%i, %i)    Cost: %i", ci, candx, candy, cost));
      if (cost < best_cost) {
        best_sad = sad;
        best_rate = rate;
        best_cost = cost;
        best_vec[0] = candx;
        best_vec[1] = candy;
      }
    }
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG, "Threshold: %i", t2));
    if (best_sad > t2) {
      /*Constant velocity predictors from the previous frame:*/
      for (ci = 0; ci < 4; ci++) {
        cands[ci][0] =
         OD_CLAMPI(mvxmin, pneighbors[ci]->bma_mvs[1][ref][0], mvxmax);
        cands[ci][1] =
         OD_CLAMPI(mvymin, pneighbors[ci]->bma_mvs[1][ref][1], mvymax);
      }
      /*The constant acceleration predictor:*/
      cands[4][0] = OD_CLAMPI(mvxmin,
       OD_DIV_ROUND_POW2(mv->bma_mvs[1][ref][0]*est->mvapw[ref][0]
       - mv->bma_mvs[2][ref][0]*est->mvapw[ref][1], 16, 0x8000), mvxmax);
      cands[4][1] = OD_CLAMPI(mvymin,
       OD_DIV_ROUND_POW2(mv->bma_mvs[1][ref][1]*est->mvapw[ref][0]
       - mv->bma_mvs[2][ref][1]*est->mvapw[ref][1], 16, 0x8000), mvymax);
      /*Examine the candidates in Set C.*/
      for (ci = 0; ci < 5; ci++) {
        candx = cands[ci][0];
        candy = cands[ci][1];
        if (od_mv_est_is_hit(est, candx, candy)) {
          OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
           "Set C predictor %i: (%i, %i) ...Skipping.", ci, candx, candy));
          continue;
        }
        /*if (od_mv_est_is_hit(est, candx, candy)) continue;*/
        od_mv_est_set_hit(est, candx, candy);
#if defined(OD_DUMP_IMAGES) && defined(OD_ANIMATE)
        if (animating) {
          od_img_draw_line(&est->enc->vis_img, x0, y0,
           x0 + (candx << 1), y0 + (candy << 1), OD_YCbCr_MVCAND);
        }
#endif
        sad = od_mv_est_bma_sad(est, ref, bx, by, candx, candy, log_mvb_sz);
        rate = od_mv_est_cand_bits(est, equal_mvs,
         candx << 1, candy << 1, pred[0], pred[1], ref, ref_pred);
        cost = (sad << OD_ERROR_SCALE) + rate*est->lambda;
        OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
         "Set C predictor %i: (%i, %i)    Cost: %i", ci, candx, candy, cost));
        if (cost < best_cost) {
          best_sad = sad;
          best_rate = rate;
          best_cost = cost;
          best_vec[0] = candx;
          best_vec[1] = candy;
        }
      }
      /*Use the same threshold for Set C as in Set B.*/
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG, "Threshold: %i", t2));
      if (best_sad > t2) {
        const int *pattern;
        int mvstate;
        int best_site;
        int nsites;
        int sitei;
        int site;
        int b;
        /*Gradient descent pattern search.*/
        mvstate = 0;
        for (;;) {
          best_site = 4;
          b = (best_vec[0] <= mvxmin) | (best_vec[0] >= mvxmax) << 1 |
           (best_vec[1] <= mvymin) << 2 | (best_vec[1] >= mvymax) << 3;
          pattern = OD_SEARCH_SITES[mvstate][b];
          nsites = OD_SEARCH_NSITES[mvstate][b];
          for (sitei = 0; sitei < nsites; sitei++) {
            site = pattern[sitei];
            candx = best_vec[0] + OD_SITE_DX[site];
            candy = best_vec[1] + OD_SITE_DY[site];
            /*For the large search patterns, our simple mechanism to move
               bounds checking out of the inner loop doesn't work (it would
               need 2 more bits, or 4 times as much table storage, and require
               4 extra compares, when there are often fewer than 4 sites).
              If the displacement is larger than +/-1 in any direction (which
               happens when site > 8), check the bounds explicitly.*/
            if (site > 8 && (candx < mvxmin || candx > mvxmax
             || candy < mvymin || candy > mvymax)) {
              continue;
            }
            if (od_mv_est_is_hit(est, candx, candy)) {
              OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
               "Pattern search %i: (%i, %i) ...Skipping.",
               site, candx, candy));
              continue;
            }
            od_mv_est_set_hit(est, candx, candy);
#if defined(OD_DUMP_IMAGES) && defined(OD_ANIMATE)
            if (animating) {
              od_img_draw_line(&est->enc->vis_img, x0, y0,
               x0 + (candx << 1), y0 + (candy << 1), OD_YCbCr_MVCAND);
            }
#endif
            sad = od_mv_est_bma_sad(est,
             ref, bx, by, candx, candy, log_mvb_sz);
            rate = od_mv_est_cand_bits(est, equal_mvs,
             candx << 1, candy << 1, pred[0], pred[1], ref, ref_pred);
            cost = (sad << OD_ERROR_SCALE) + rate*est->lambda;
            OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
             "Pattern search %i: (%i, %i)    Cost: %i",
             site, candx, candy, cost));
            if (cost < best_cost) {
              best_sad = sad;
              best_rate = rate;
              best_cost = cost;
              best_site = site;
            }
          }
          mvstate = OD_SEARCH_STATES[mvstate][best_site];
          best_vec[0] += OD_SITE_DX[best_site];
          best_vec[1] += OD_SITE_DY[best_site];
          if (mvstate == OD_SEARCH_STATE_DONE) break;
        }
      }
    }
  }
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
   "Finished. Best vector: (%i, %i)  Best cost %i",
   best_vec[0], best_vec[1], best_cost));
#if defined(OD_DUMP_IMAGES) && defined(OD_ANIMATE)
  if (animating) {
    char iter_label[16];
    const od_offset *anc;
    od_mv_grid_pt *amvg;
    int nanc;
    int ai;
    int ax;
    int ay;
    mvg->valid = 1;
    nanc = OD_NANCESTORS[vy & OD_MVB_MASK][vx & OD_MVB_MASK];
    anc = OD_ANCESTORS[vy & OD_MVB_MASK][vx & OD_MVB_MASK];
    for (ai = 0; ai < nanc; ai++) {
      ax = vx + anc[ai][0];
      if (ax < 0 || ax > ((vx + OD_MVB_MASK) & ~OD_MVB_MASK)) continue;
      ay = vy + anc[ai][1];
      if (ay < 0 || ay > ((vy + OD_MVB_MASK) & ~OD_MVB_MASK)) continue;
      amvg = state->mv_grid[ay] + ax;
      amvg->valid = 1;
    }
    sprintf(iter_label, "ani%08i", est->enc->ani_iter++);
    od_state_dump_img(state, &est->enc->vis_img, iter_label);
  }
#endif
  mv->bma_mvs[0][ref][0] = best_vec[0];
  mv->bma_mvs[0][ref][1] = best_vec[1];
  /*previous_cost is our previous best cost from a previous pass of phase 1.*/
  previous_cost = (mv->bma_sad << OD_ERROR_SCALE) + mv->mv_rate*est->lambda;
  if (must_update || (best_cost < previous_cost)) {
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
     "Found a better SAD then previous best."));
    mvg->mv[0] = best_vec[0] << 3;
    mvg->mv[1] = best_vec[1] << 3;
    mvg->ref = ref;
    mvg->valid = 1;
    mv->bma_sad = best_sad;
    mv->mv_rate = best_rate;
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
     "Initialized MV (%2i, %2i): (%3i, %3i), SAD: %i, ref: %i",
     vx, vy, best_vec[0], best_vec[1], best_sad, ref));
#if defined(OD_ENABLE_ASSERTIONS) || defined(OD_ENABLE_LOGGING)
    {
      OD_ASSERT(mvg->mv[0] <= mvxmax << 3);
      OD_ASSERT(mvg->mv[0] >= mvxmin << 3);
      OD_ASSERT(mvg->mv[1] <= mvymax << 3);
      OD_ASSERT(mvg->mv[1] >= mvymin << 3);
      mv->mv_rate = od_mv_est_bits(est, vx, vy, 2);
      if (mv->mv_rate != best_rate) {
        OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_ERR,
         "Failure in MV rate init: %i != %i",
         mv->mv_rate, best_rate));
        OD_MV_EST_LOG_PRED(est, vx, vy, 2);
        OD_ASSERT(mv->mv_rate == best_rate);
      }
    }
#endif
  }
}

static void od_mv_est_init_mvs(od_mv_est_ctx *est, int ref, int must_update) {
  od_state *state;
  int nhmvbs;
  int nvmvbs;
  int vx;
  int vy;
  state = &est->enc->state;
  nhmvbs = state->nhmvbs;
  nvmvbs = state->nvmvbs;
  /*Move the motion vector predictors back a frame.*/
  for (vy = 0; vy <= nvmvbs; vy++) {
    for (vx = 0; vx <= nhmvbs; vx++) {
      od_mv_node *mv;
      mv = est->mvs[vy] + vx;
      OD_MOVE(mv->bma_mvs + 1, mv->bma_mvs + 0, 2);
    }
  }
  /*We initialize MVs a MVB at a time for cache coherency.
    Proceeding level-by-level would involve less branching and less complex
     code, but the SADs dominate.
    Initialization MUST proceed in order by level to ensure the necessary
     predictors are available.
    Order within a level does not matter, except for level 0.
    Level 0 is the only level to use predictors outside the current MVB,
     and must proceed in raster order, one row/column _ahead_ of the
     rest of the MVB (so that we have all four corners available to predict
     the lower levels in the current MVB).*/
  for (vx = 0; vx <= nhmvbs; vx += OD_MVB_DELTA0) {
    od_mv_est_init_mv(est, ref, vx, 0, must_update);
  }
  for (vy = 0; vy < nvmvbs; vy += OD_MVB_DELTA0) {
    od_mv_est_init_mv(est, ref, 0, vy + OD_MVB_DELTA0, must_update);
    for (vx = 0; vx < nhmvbs; vx += OD_MVB_DELTA0) {
      int log_mvb_sz;
      int level;
      /*Level 0 vertex.*/
      od_mv_est_init_mv(est, ref, vx + OD_MVB_DELTA0, vy + OD_MVB_DELTA0,
       must_update);
      /*All other levels.*/
      for(log_mvb_sz = OD_LOG_MVB_DELTA0, level = 1;
       log_mvb_sz-- > 0 && est->level_max >= level; level++) {
        int cx;
        int cy;
        int mvb_sz;
        mvb_sz = 1 << log_mvb_sz;
        /*Odd level vertices.*/
        for (cy = vy + mvb_sz; cy < vy + OD_MVB_DELTA0; cy += 2*mvb_sz) {
          for( cx = vx + mvb_sz; cx < vx + OD_MVB_DELTA0; cx += 2*mvb_sz) {
            od_mv_est_init_mv(est, ref, cx, cy, must_update);
          }
        }
        level++;
        if (est->level_max < level) break;
        /*Even level vertices.*/
        /*Add even-level vertices on the top/left edges of the frame as extra
           vertices in the first row/column of MVBs.
          Unlike other vertices on the edges of an MVB, they can use parents to
           the right/below them as predictors (or otherwise they would have no
           predictors).*/
        /*Skip the cy == vy row unless we're at the top of the frame.*/
        for (cy = vy + mvb_sz*!!vy; cy <= vy + OD_MVB_DELTA0; cy += mvb_sz) {
          /*Even level vertices appear in a quincunx pattern.
            We want to start every other row at an mvb_sz offset, and also to
             skip the first column on the rows flush with the edge of the block
             unless we're on the left edge of the whole frame.*/
          for( cx = vx + (cy & mvb_sz ? 2*mvb_sz*!!vx : mvb_sz);
           cx <= vx + OD_MVB_DELTA0; cx += 2*mvb_sz) {
            od_mv_est_init_mv(est, ref, cx, cy, must_update);
          }
        }
      }
    }
  }
}

/*STAGE 2: DECIMATION.*/

/*Merging domains.
  These are stored as lists of offsets to the vertices in the domain.
  Note that vertices in the merging domain must appear in order from finest
   scale (largest level) to coarsest (smallest level).
  Each list ends with the vertex (0, 0), the actual vertex be decimated.*/

/*Level 6 vertex:
                        6*/
static const od_offset OD_MERGEDOM6[1] = {
  { 0, 0 },
};

/*Level 5 vertex:
                        6
                        |
                      6-5-6
                        |
                        6*/
static const od_offset OD_MERGEDOM5[5] = {
  {  0, -1 }, { -1,  0 }, {  1,  0 }, {  0,  1 }, {  0,  0 }
};

/*Level 4 vertex:
                      6   6
                      |   |
                    6-5-6-5-6
                      | | |
                      6-4-6
                      | | |
                    6-5-6-5-6
                      |   |
                      6   6*/
static const od_offset OD_MERGEDOM4[17] = {
  { -1, -2 }, {  1, -2 }, { -2, -1 }, {  0, -1 }, {  2, -1 }, { -1,  0 },
  {  1,  0 }, { -2,  1 }, {  0,  1 }, {  2,  1 }, { -1,  2 }, {  1,  2 },
  { -1, -1 }, {  1, -1 }, { -1,  1 }, {  1,  1 }, {  0,  0 }
};

/*Level 3 vertex:
                      6   6
                      |   |
                    6-5-6-5-6
                      | | |
                  6   6-4-6   6
                  |   | | |   |
                6-5-6-5-6-5-6-5-6
                  | | | | | | |
                  6-4-6-3-6-4-6
                  | | | | | | |
                6-5-6-5-6-5-6-5-6
                  |   | | |   |
                  6   6-4-6   6
                      | | |
                    6-5-6-5-6
                      |   |
                      6   6*/
static const od_offset OD_MERGEDOM3[49] = {
  { -1, -4 }, {  1, -4 }, { -2, -3 }, {  0, -3 }, {  2, -3 }, { -3, -2 },
  { -1, -2 }, {  1, -2 }, {  3, -2 }, { -4, -1 }, { -2, -1 }, {  0, -1 },
  {  2, -1 }, {  4, -1 }, { -3,  0 }, { -1,  0 }, {  1,  0 }, {  3,  0 },
  { -4,  1 }, { -2,  1 }, {  0,  1 }, {  2,  1 }, {  4,  1 }, { -3,  2 },
  { -1,  2 }, {  1,  2 }, {  3,  2 }, { -2,  3 }, {  0,  3 }, {  2,  3 },
  { -1,  4 }, {  1,  4 }, { -1, -3 }, {  1, -3 }, { -3, -1 }, { -1, -1 },
  {  1, -1 }, {  3, -1 }, { -3,  1 }, { -1,  1 }, {  1,  1 }, {  3,  1 },
  { -1,  3 }, {  1,  3 }, {  0, -2 }, { -2,  0 }, {  2,  0 }, {  0,  2 },
  {  0,  0 }
};

/*Level 2 vertex:
                  6   6   6   6
                  |   |   |   |
                6-5-6-5-6-5-6-5-6
                  | | |   | | |
              6   6-4-6   6-4-6   6
              |   | | |   | | |   |
            6-5-6-5-6-5-6-5-6-5-6-5-6
              | | | | | | | | | | |
              6-4-6-3-6-4-6-3-6-4-6
              | | | | | | | | | | |
            6-5-6-5-6-5-6-5-6-5-6-5-6
              |   | | | | | | |   |
              6   6-4-6-2-6-4-6   6
              |   | | | | | | |   |
            6-5-6-5-6-5-6-5-6-5-6-5-6
              |   | | | | | | | | |
              6-4-6-3-6-4-6-3-6-4-6
              | | | | | | | | | | |
            6-5-6-5-6-5-6-5-6-5-6-5-6
              |   | | |   | | |   |
              6   6-4-6   6-4-6   6
                  | | |   | | |
                6-5-6-5-6-5-6-5-6
                  |   |   |   |
                  6   6   6   6*/
static const od_offset OD_MERGEDOM2[125] = {
  { -3, -6 }, { -1, -6 }, {  1, -6 }, {  3, -6 }, { -4, -5 }, { -2, -5 },
  {  0, -5 }, {  2, -5 }, {  4, -5 }, { -5, -4 }, { -3, -4 }, { -1, -4 },
  {  1, -4 }, {  3, -4 }, {  5, -4 }, { -6, -3 }, { -4, -3 }, { -2, -3 },
  {  0, -3 }, {  2, -3 }, {  4, -3 }, {  6, -3 }, { -5, -2 }, { -3, -2 },
  { -1, -2 }, {  1, -2 }, {  3, -2 }, {  5, -2 }, { -6, -1 }, { -4, -1 },
  { -2, -1 }, {  0, -1 }, {  2, -1 }, {  4, -1 }, {  6, -1 }, { -5,  0 },
  { -3,  0 }, { -1,  0 }, {  1,  0 }, {  3,  0 }, {  5,  0 }, { -6,  1 },
  { -4,  1 }, { -2,  1 }, {  0,  1 }, {  2,  1 }, {  4,  1 }, {  6,  1 },
  { -5,  2 }, { -3,  2 }, { -1,  2 }, {  1,  2 }, {  3,  2 }, {  5,  2 },
  { -6,  3 }, { -4,  3 }, { -2,  3 }, {  0,  3 }, {  2,  3 }, {  4,  3 },
  {  6,  3 }, { -5,  4 }, { -3,  4 }, { -1,  4 }, {  1,  4 }, {  3,  4 },
  {  5,  4 }, { -4,  5 }, { -2,  5 }, {  0,  5 }, {  2,  5 }, {  4,  5 },
  { -3,  6 }, { -1,  6 }, {  1,  6 }, {  3,  6 }, { -3, -5 }, { -1, -5 },
  {  1, -5 }, {  3, -5 }, { -5, -3 }, { -3, -3 }, { -1, -3 }, {  1, -3 },
  {  3, -3 }, {  5, -3 }, { -5, -1 }, { -3, -1 }, { -1, -1 }, {  1, -1 },
  {  3, -1 }, {  5, -1 }, { -5,  1 }, { -3,  1 }, { -1,  1 }, {  1,  1 },
  {  3,  1 }, {  5,  1 }, { -5,  3 }, { -3,  3 }, { -1,  3 }, {  1,  3 },
  {  3,  3 }, {  5,  3 }, { -3,  5 }, { -1,  5 }, {  1,  5 }, {  3,  5 },
  { -2, -4 }, {  2, -4 }, { -4, -2 }, {  0, -2 }, {  4, -2 }, { -2,  0 },
  {  2,  0 }, { -4,  2 }, {  0,  2 }, {  4,  2 }, { -2,  4 }, {  2,  4 },
  { -2, -2 }, {  2, -2 }, { -2,  2 }, {  2,  2 }, {  0,  0 }
};

/*Level 1 vertex:
                  6   6   6   6
                  |   |   |   |
                6-5-6-5-6-5-6-5-6
                  | | |   | | |
              6   6-4-6   6-4-6   6
              |   | | |   | | |   |
            6-5-6-5-6-5-6-5-6-5-6-5-6
              | | | | | | | | | | |
          6   6-4-6-3-6-4-6-3-6-4-6   6
          |   | | | | | | | | | | |   |
        6-5-6-5-6-5-6-5-6-5-6-5-6-5-6-5-6
          | | |   | | | | | | |   | | |
      6   6-4-6   6-4-6-2-6 4 6   6 4 6   6
      |   | | |   | | | | | | |   | | |   |
    6-5-6-5-6-5-6-5-6-5-6-5-6-5-6-5-6-5-6-5-6
      | | | | | | | | | | | | | | | | | | |
      6-4-6-3-6-4-6-3-6-4-6-3-6-4-6-3-6-4-6
      | | | | | | | | | | | | | | | | | | |
    6-5-6-5-6-5-6-5-6-5-6-5-6-5-6-5-6-5-6-5-6
      |   | | | | | | | | | | | | | | |   |
      6   6-4-6-2-6-4-6-1-6-4-6-2-6-4-6   6
      |   | | | | | | | | | | | | | | |   |
    6-5-6-5-6-5-6-5-6-5-6-5-6-5-6-5-6-5-6-5-6
      | | | | | | | | | | | | | | | | | | |
      6-4-6-3-6-4-6-3-6-4-6-3-6-4-6-3-6-4-6
      | | | | | | | | | | | | | | | | | | |
    6-5-6-5-6-5-6-5-6-5-6-5-6-5-6-5-6-5-6-5-6
      |   | | |   | | | | | | |   | | |   |
      6   6-4-6   6-4-6-2-6-4-6   6-4-6   6
          | | |   | | | | | | |   | | |
        6-5-6-5-6-5-6-5-6-5-6-5-6-5-6-5-6
          |   | | | | | | | | | | |   |
          6   6-4-6-3-6-4-6-3-6-4-6   6
              | | | | | | | | | | |
            6-5-6-5-6-5-6-5-6-5-6-5-6
              |   | | |   | | |   |
              6   6-4-6   6-4-6   6
                  | | |   | | |
                6-5-6-5-6-5-6-5-6
                  |   |   |   |
                  6   6   6   6*/
static const od_offset OD_MERGEDOM1[297] = {
  { -3,-10 }, { -1,-10 }, {  1,-10 }, {  3,-10 }, { -4, -9 }, { -2, -9 },
  {  0, -9 }, {  2, -9 }, {  4, -9 }, { -5, -8 }, { -3, -8 }, { -1, -8 },
  {  1, -8 }, {  3, -8 }, {  5, -8 }, { -6, -7 }, { -4, -7 }, { -2, -7 },
  {  0, -7 }, {  2, -7 }, {  4, -7 }, {  6, -7 }, { -7, -6 }, { -5, -6 },
  { -3, -6 }, { -1, -6 }, {  1, -6 }, {  3, -6 }, {  5, -6 }, {  7, -6 },
  { -8, -5 }, { -6, -5 }, { -4, -5 }, { -2, -5 }, {  0, -5 }, {  2, -5 },
  {  4, -5 }, {  6, -5 }, {  8, -5 }, { -9, -4 }, { -7, -4 }, { -5, -4 },
  { -3, -4 }, { -1, -4 }, {  1, -4 }, {  3, -4 }, {  5, -4 }, {  7, -4 },
  {  9, -4 }, {-10, -3 }, { -8, -3 }, { -6, -3 }, { -4, -3 }, { -2, -3 },
  {  0, -3 }, {  2, -3 }, {  4, -3 }, {  6, -3 }, {  8, -3 }, { 10, -3 },
  { -9, -2 }, { -7, -2 }, { -5, -2 }, { -3, -2 }, { -1, -2 }, {  1, -2 },
  {  3, -2 }, {  5, -2 }, {  7, -2 }, {  9, -2 }, {-10, -1 }, { -8, -1 },
  { -6, -1 }, { -4, -1 }, { -2, -1 }, {  0, -1 }, {  2, -1 }, {  4, -1 },
  {  6, -1 }, {  8, -1 }, { 10, -1 }, { -9,  0 }, { -7,  0 }, { -5,  0 },
  { -3,  0 }, { -1,  0 }, {  1,  0 }, {  3,  0 }, {  5,  0 }, {  7,  0 },
  {  9,  0 }, {-10,  1 }, { -8,  1 }, { -6,  1 }, { -4,  1 }, { -2,  1 },
  {  0,  1 }, {  2,  1 }, {  4,  1 }, {  6,  1 }, {  8,  1 }, { 10,  1 },
  { -9,  2 }, { -7,  2 }, { -5,  2 }, { -3,  2 }, { -1,  2 }, {  1,  2 },
  {  3,  2 }, {  5,  2 }, {  7,  2 }, {  9,  2 }, {-10,  3 }, { -8,  3 },
  { -6,  3 }, { -4,  3 }, { -2,  3 }, {  0,  3 }, {  2,  3 }, {  4,  3 },
  {  6,  3 }, {  8,  3 }, { 10,  3 }, { -9,  4 }, { -7,  4 }, { -5,  4 },
  { -3,  4 }, { -1,  4 }, {  1,  4 }, {  3,  4 }, {  5,  4 }, {  7,  4 },
  {  9,  4 }, { -8,  5 }, { -6,  5 }, { -4,  5 }, { -2,  5 }, {  0,  5 },
  {  2,  5 }, {  4,  5 }, {  6,  5 }, {  8,  5 }, { -7,  6 }, { -5,  6 },
  { -3,  6 }, { -1,  6 }, {  1,  6 }, {  3,  6 }, {  5,  6 }, {  7,  6 },
  { -6,  7 }, { -4,  7 }, { -2,  7 }, {  0,  7 }, {  2,  7 }, {  4,  7 },
  {  6,  7 }, { -5,  8 }, { -3,  8 }, { -1,  8 }, {  1,  8 }, {  3,  8 },
  {  5,  8 }, { -4,  9 }, { -2,  9 }, {  0,  9 }, {  2,  9 }, {  4,  9 },
  { -3, 10 }, { -1, 10 }, {  1, 10 }, {  3, 10 }, { -3, -9 }, { -1, -9 },
  {  1, -9 }, {  3, -9 }, { -5, -7 }, { -3, -7 }, { -1, -7 }, {  1, -7 },
  {  3, -7 }, {  5, -7 }, { -7, -5 }, { -5, -5 }, { -3, -5 }, { -1, -5 },
  {  1, -5 }, {  3, -5 }, {  5, -5 }, {  7, -5 }, { -9, -3 }, { -7, -3 },
  { -5, -3 }, { -3, -3 }, { -1, -3 }, {  1, -3 }, {  3, -3 }, {  5, -3 },
  {  7, -3 }, {  9, -3 }, { -9, -1 }, { -7, -1 }, { -5, -1 }, { -3, -1 },
  { -1, -1 }, {  1, -1 }, {  3, -1 }, {  5, -1 }, {  7, -1 }, {  9, -1 },
  { -9,  1 }, { -7,  1 }, { -5,  1 }, { -3,  1 }, { -1,  1 }, {  1,  1 },
  {  3,  1 }, {  5,  1 }, {  7,  1 }, {  9,  1 }, { -9,  3 }, { -7,  3 },
  { -5,  3 }, { -3,  3 }, { -1,  3 }, {  1,  3 }, {  3,  3 }, {  5,  3 },
  {  7,  3 }, {  9,  3 }, { -7,  5 }, { -5,  5 }, { -3,  5 }, { -1,  5 },
  {  1,  5 }, {  3,  5 }, {  5,  5 }, {  7,  5 }, { -5,  7 }, { -3,  7 },
  { -1,  7 }, {  1,  7 }, {  3,  7 }, {  5,  7 }, { -3,  9 }, { -1,  9 },
  {  1,  9 }, {  3,  9 }, { -2, -8 }, {  2, -8 }, { -4, -6 }, {  0, -6 },
  {  4, -6 }, { -6, -4 }, { -2, -4 }, {  2, -4 }, {  6, -4 }, { -8, -2 },
  { -4, -2 }, {  0, -2 }, {  4, -2 }, {  8, -2 }, { -6,  0 }, { -2,  0 },
  {  2,  0 }, {  6,  0 }, { -8,  2 }, { -4,  2 }, {  0,  2 }, {  4,  2 },
  {  8,  2 }, { -6,  4 }, { -2,  4 }, {  2,  4 }, {  6,  4 }, { -4,  6 },
  {  0,  6 }, {  4,  6 }, { -2,  8 }, {  2,  8 }, { -2, -6 }, {  2, -6 },
  { -6, -2 }, { -2, -2 }, {  2, -2 }, {  6, -2 }, { -6,  2 }, { -2,  2 },
  {  2,  2 }, {  6,  2 }, { -2,  6 }, {  2,  6 }, {  0, -4 }, { -4,  0 },
  {  4,  0 }, {  0,  4 }, {  0,  0 }
};

/*The merging domain for a vertex, indexed by (level - 1).*/
static const od_offset *OD_MERGEDOM[OD_MC_LEVEL_MAX] = {
  OD_MERGEDOM1,
  OD_MERGEDOM2,
  OD_MERGEDOM3,
  OD_MERGEDOM4,
  OD_MERGEDOM5,
  OD_MERGEDOM6
};

/*Error support regions.
  These are the blocks whose SAD will change after decimating a vertex at a
   given level, assuming no other vertices in the mesh have been decimated.
  Vertices in the figures at a higher level than the one removed illustrate one
   possible configuration (the first in raster order in a level-0 MVB); there
   may be others.*/
struct od_mv_err_node {
  int dx;
  int dy;
  int log_mvb_sz;
};

/*Level 6 support:
                      6-5-6
                      |/|\|
                      0-.-4
                      |\|/|
                      6-5-6*/
static const od_mv_err_node OD_ERRDOM6[4] = {
  { -1, -1, 0 }, {  0, -1, 0 }, { -1,  0, 0 }, {  0,  0, 0 }
};

/*Level 5 support:
                      6-5-6
                      |/|\|
                    6-0-.-4-6
                    |/|   |\|
                    5-.   .-5
                    |\|   |/|
                    6-4-.-3-6
                      |\|/|
                      6-5-6*/
static const od_mv_err_node OD_ERRDOM5[9] = {
  { -1, -2, 0 }, {  0, -2, 0 }, { -2, -1, 0 }, {  1, -1, 0 },
  { -2,  0, 0 }, {  1,  0, 0 }, { -1,  1, 0 }, {  0,  1, 0 },
  { -1, -1, 1 }
};

/*Level 4 support:
                    6-5-6-5-6
                    |/|\|/|\|
                  6-4-.-3-.-4-6
                  |/|  /|\  |\|
                  5-. / | \ .-5
                  |\|/  |  \|/|
                  6-0---.---2-6
                  |/|\  |  /|\|
                  5-. \ | / .-5
                  |\|  \|/  |/|
                  6-4-.-3-.-4-6
                    |\|/|\|/|
                    6-5-6-5-6*/
static const od_mv_err_node OD_ERRDOM4[20] = {
  { -2, -3, 0 }, { -1, -3, 0 }, {  0, -3, 0 }, {  1, -3, 0 },
  { -3, -2, 0 }, {  2, -2, 0 }, { -3, -1, 0 }, {  2, -1, 0 },
  { -3,  0, 0 }, {  2,  0, 0 }, { -3,  1, 0 }, {  2,  1, 0 },
  { -2,  2, 0 }, { -1,  2, 0 }, {  0,  2, 0 }, {  1,  2, 0 },
  { -2, -2, 1 }, {  0, -2, 1 }, { -2,  0, 1 }, {  0,  0, 1 }
};

/*Level 3 support:
                    6-5-6-5-6
                    |/|\|/|\|
                  6-4-.-3-.-4-6
                  |/|  /|\  |\|
                6-5-. / | \ .-5-6
                |/| |/  |  \| |\|
              6-4-.-0---.---2-.-4-6
              |/|  /|       |\  |\|
              5-. / |       | \ .-5
              |\|/  |       |  \|/|
              6-3---.       .---3-6
              |/|\  |       |  /|\|
              5-. \ |       | / .-5
              |\|  \|       |/  |/|
              6-4-.-2---.---1-.-4-6
                |\| |\  |  /| |/|
                6-5-. \ | / .-5-6
                  |\|  \|/  |/|
                  6-4-.-3-.-4-6
                    |\|/|\|/|
                    6-5-6-5-6*/
static const od_mv_err_node OD_ERRDOM3[37] = {
  { -2, -5, 0 }, { -1, -5, 0 }, {  0, -5, 0 }, {  1, -5, 0 },
  { -3, -4, 0 }, {  2, -4, 0 }, { -4, -3, 0 }, { -3, -3, 0 },
  {  2, -3, 0 }, {  3, -3, 0 }, { -5, -2, 0 }, {  4, -2, 0 },
  { -5, -1, 0 }, {  4, -1, 0 }, { -5,  0, 0 }, {  4,  0, 0 },
  { -5,  1, 0 }, {  4,  1, 0 }, { -4,  2, 0 }, { -3,  2, 0 },
  {  2,  2, 0 }, {  3,  2, 0 }, { -3,  3, 0 }, {  2,  3, 0 },
  { -2,  4, 0 }, { -1,  4, 0 }, {  0,  4, 0 }, {  1,  4, 0 },
  { -2, -4, 1 }, {  0, -4, 1 }, { -4, -2, 1 }, {  2, -2, 1 },
  { -4,  0, 1 }, {  2,  0, 1 }, { -2,  2, 1 }, {  0,  2, 1 },
  { -2, -2, 2 }
};

/*Level 2 support:
                6-5-6-5-6-5-6-5-6
                |/|\|/|\|/|\|/|\|
              6-4-.-3-.-4-.-3-.-4-6
              |/|  /|\  |  /|\  |\|
            6-5-. / | \ | / | \ .-5-6
            |/|\|/  |  \|/  |  \|/|\|
          6-4-.-2---.---1---.---2---4-6
          |/|  /|      /|\      |\  |\|
          5-. / |     / | \     | \ .-5
          |\|/  |    /  |  \    |  \|/|
          6-3---.   /   |   \   .---3-6
          |/|\  |  /    |    \  |  /|\|
          5-. \ | /     |     \ | / .-5
          |\|  \|/      |      \|/  |/|
          6-4---0-------.-------0---4-6
          |/|  /|\      |      /|\  |\|
          5-. / | \     |     / | \ .-5
          |\|/  |  \    |    /  |  \|/|
          6-3---.   \   |   /   .---3-6
          |/|\  |    \  |  /    |  /|\|
          5-. \ |     \ | /     | / .-5
          |\|  \|      \|/      |/  |/|
          6-4-.-2---.---1---.---2-.-4-6
            |\|/|\  |  /|\  |  /|\|/|
            6-5-. \ | / | \ | / .-5-6
              |\|  \|/  |  \|/  |/|
              6-4-.-3-.-4-.-3-.-4-6
                |\|/|\|/|\|/|\|/|
                6-5-6-5-6-5-6-5-6*/
static const od_mv_err_node OD_ERRDOM2[64] = {
  { -4, -7, 0 }, { -3, -7, 0 }, { -2, -7, 0 }, { -1, -7, 0 },
  {  0, -7, 0 }, {  1, -7, 0 }, {  2, -7, 0 }, {  3, -7, 0 },
  { -5, -6, 0 }, {  4, -6, 0 }, { -6, -5, 0 }, { -5, -5, 0 },
  {  4, -5, 0 }, {  5, -5, 0 }, { -7, -4, 0 }, {  6, -4, 0 },
  { -7, -3, 0 }, {  6, -3, 0 }, { -7, -2, 0 }, {  6, -2, 0 },
  { -7, -1, 0 }, {  6, -1, 0 }, { -7,  0, 0 }, {  6,  0, 0 },
  { -7,  1, 0 }, {  6,  1, 0 }, { -7,  2, 0 }, {  6,  2, 0 },
  { -7,  3, 0 }, {  6,  3, 0 }, { -6,  4, 0 }, { -5,  4, 0 },
  {  4,  4, 0 }, {  5,  4, 0 }, { -5,  5, 0 }, {  4,  5, 0 },
  { -4,  6, 0 }, { -3,  6, 0 }, { -2,  6, 0 }, { -1,  6, 0 },
  {  0,  6, 0 }, {  1,  6, 0 }, {  2,  6, 0 }, {  3,  6, 0 },
  { -4, -6, 1 }, { -2, -6, 1 }, {  0, -6, 1 }, {  2, -6, 1 },
  { -6, -4, 1 }, {  4, -4, 1 }, { -6, -2, 1 }, {  4, -2, 1 },
  { -6,  0, 1 }, {  4,  0, 1 }, { -6,  2, 1 }, {  4,  2, 1 },
  { -4,  4, 1 }, { -2,  4, 1 }, {  0,  4, 1 }, {  2,  4, 1 },
  { -4, -4, 2 }, {  0, -4, 2 }, { -4,  0, 2 }, {  0,  0, 2 }
};

/*Level 1 support:
                6-5-6-5-6-5-6-5-6
                |/|\|/|\|/|\|/|\|
              6-4-.-3-.-4-.-3-.-4-6
              |/|  /|\  |  /|\  |\|
            6-5-. / | \ | / | \ .-5-6
            |/|\|/  |  \|/  |  \|/|\|
          6-4-.-2---.---1---.---2-.-4-6
          |/|  /|      /|\      |\  |\|
        6-5-. / |     / | \     | \ .-5-6
        |/|\|/  |    /  |  \    |  \|/|\|
      6-4-.-3---.   /   |   \   .---3-.-4-6
      |/|  /|\  |  /    |    \  |  /|\  |\|
    6-5-. / | \ | /     |     \ | / | \ .-5-6
    |/|\|/  |  \|/      |      \|/  |  \|/|\|
  6-4-.-2---.---0-------.-------0---.---2-.-4-6
  |/|  /|      /|               |\      |\  |\|
  5-. / |     / |               | \     | \ .-5
  |\|/  |    /  |               |  \    |  \|/|
  6-3---.   /   |               |   \   .---3-6
  |/|\  |  /    |               |    \  |  /|\|
  5-. \ | /     |               |     \ | / .-5
  |\|  \|/      |               |      \|/  |/|
  6-4---1-------.               .-------1---4-6
  |/|  /|\      |               |      /|\  |\|
  5-. / | \     |               |     / | \ .-5
  |\|/  |  \    |               |    /  |  \|/|
  6-3---.   \   |               |   /   .---3-6
  |/|\  |    \  |               |  /    |  /|\|
  5-. \ |     \ |               | /     | / .-5
  |\|  \|      \|               |/      |/  |/|
  6-4-.-2---.---0-------.-------0---.---2-.-4-6
    |/|\|\  |  /|\      |      /|\  |  /|\|/|
    6-5-. \ | / | \     |     / | \ | / .-5-6
      |\|  \|/  |  \    |    /  |  \|/  |/|
      6-4-.-3---.   \   |   /   .---3-.-4-6
        |\|/|\  |    \  |  /    |  /|\|/|
        6-5-. \ |     \ | /     | / .-5-6
          |\|  \|      \|/      |/  |/|
          6-4-.-2---.---1---.---2-.-4-6
            |\|/|\  |  /|\  |  /|\|/|
            6-5-. \ | / | \ | / .-5-6
              |\|  \|/  |  \|/  |/|
              6-4-.-3-.-4-.-3-.-4-6
                |\|/|\|/|\|/|\|/|
                6-5-6-5-6-5-6-5-6*/
static const od_mv_err_node OD_ERRDOM1[105] = {
  { -4,-11, 0 }, { -3,-11, 0 }, { -2,-11, 0 }, { -1,-11, 0 },
  {  0,-11, 0 }, {  1,-11, 0 }, {  2,-11, 0 }, {  3,-11, 0 },
  { -5,-10, 0 }, {  4,-10, 0 }, { -6, -9, 0 }, { -5, -9, 0 },
  {  4, -9, 0 }, {  5, -9, 0 }, { -7, -8, 0 }, {  6, -8, 0 },
  { -8, -7, 0 }, { -7, -7, 0 }, {  6, -7, 0 }, {  7, -7, 0 },
  { -9, -6, 0 }, {  8, -6, 0 }, {-10, -5, 0 }, { -9, -5, 0 },
  {  8, -5, 0 }, {  9, -5, 0 }, {-11, -4, 0 }, { 10, -4, 0 },
  {-11, -3, 0 }, { 10, -3, 0 }, {-11, -2, 0 }, { 10, -2, 0 },
  {-11, -1, 0 }, { 10, -1, 0 }, {-11,  0, 0 }, { 10,  0, 0 },
  {-11,  1, 0 }, { 10,  1, 0 }, {-11,  2, 0 }, { 10,  2, 0 },
  {-11,  3, 0 }, { 10,  3, 0 }, {-10,  4, 0 }, { -9,  4, 0 },
  {  8,  4, 0 }, {  9,  4, 0 }, { -9,  5, 0 }, {  8,  5, 0 },
  { -8,  6, 0 }, { -7,  6, 0 }, {  6,  6, 0 }, {  7,  6, 0 },
  { -7,  7, 0 }, {  6,  7, 0 }, { -6,  8, 0 }, { -5,  8, 0 },
  {  4,  8, 0 }, {  5,  8, 0 }, { -5,  9, 0 }, {  4,  9, 0 },
  { -4, 10, 0 }, { -3, 10, 0 }, { -2, 10, 0 }, { -1, 10, 0 },
  {  0, 10, 0 }, {  1, 10, 0 }, {  2, 10, 0 }, {  3, 10, 0 },
  { -4,-10, 1 }, { -2,-10, 1 }, {  0,-10, 1 }, {  2,-10, 1 },
  { -6, -8, 1 }, {  4, -8, 1 }, { -8, -6, 1 }, { -6, -6, 1 },
  {  4, -6, 1 }, {  6, -6, 1 }, {-10, -4, 1 }, {  8, -4, 1 },
  {-10, -2, 1 }, {  8, -2, 1 }, {-10,  0, 1 }, {  8,  0, 1 },
  {-10,  2, 1 }, {  8,  2, 1 }, { -8,  4, 1 }, { -6,  4, 1 },
  {  4,  4, 1 }, {  6,  4, 1 }, { -6,  6, 1 }, {  4,  6, 1 },
  { -4,  8, 1 }, { -2,  8, 1 }, {  0,  8, 1 }, {  2,  8, 1 },
  { -4, -8, 2 }, {  0, -8, 2 }, { -8, -4, 2 }, {  4, -4, 2 },
  { -8,  0, 2 }, {  4,  0, 2 }, { -4,  4, 2 }, {  0,  4, 2 },
  { -4, -4, 3 }
};

/*The number of blocks in each decimated error domain.*/
static const int OD_NERRDOM[OD_MC_LEVEL_MAX] = { 105, 64, 37, 20, 9, 4 };
/*The error domain for a vertex, indexed by (level - 1).*/
static const od_mv_err_node *OD_ERRDOM[OD_MC_LEVEL_MAX] = {
  OD_ERRDOM1,
  OD_ERRDOM2,
  OD_ERRDOM3,
  OD_ERRDOM4,
  OD_ERRDOM5,
  OD_ERRDOM6
};

/*Returns a negative value, 0, or a positive value, depending on whether
  -dd1/dr1 is less, equal or greater than -dd2/dr2.*/
static int od_mv_dddr_cmp(int32_t dd1, int dr1,
 int32_t dd2, int dr2) {
  int64_t diff;
  /*dr == 0 (infinite slope) is a special case.
    We handle it as the limit of the non-infinite case.
    We remove all vertices where it is free to do so and it helps distortion
     immediately, and delay looking at vertices that hurt distortion but don't
     affect rate until they are the only thing left.
    Vertices that change neither rate nor distortion we are free to remove at
     any time, but we choose to do so up front, before vertices that do affect
     rate.
    Otherwise they would get sandwiched between "things that hurt distortion
     and hurt rate" and "things that hurt distortion and don't affect rate".
    This gives the following total ordering (from first to last):
                      ||                  ||         ||
      dr == 0, dd < 0 || dr == 0, dd == 0 || dr != 0 || dr == 0, dd > 0
                      ||                  ||         ||
    When dr != 0, items are sorted by (-dd/dr), whereas in the
     classes with dr == 0, they are sorted by dd directly (not negated).*/
  if (dr1 == 0) {
    return dr2 == 0 ? OD_SIGNI(dd1 - dd2) : (OD_SIGNI(dd1) << 1) - 1;
  }
  else if (dr2 == 0) return (OD_SIGNI(-dd2) << 1) + 1;
  diff = dd2*(int64_t)dr1 - dd1*(int64_t)dr2;
  return OD_SIGNI(diff);
}

/*Compare two nodes on the decimation heap.*/
static int od_mv_dec_cmp(od_mv_node *n1, od_mv_node *n2) {
  return od_mv_dddr_cmp(n1->dd, n1->dr, n2->dd, n2->dr);
}

/*Swap the two nodes on the decimation heap at indices p and q.*/
static void od_mv_dec_heap_swap(od_mv_node **heap, int p, int q) {
  od_mv_node *t;
  heap[p]->heapi = q;
  heap[q]->heapi = p;
  t = heap[p];
  heap[p] = heap[q];
  heap[q] = t;
}

/*Convert the list of nodes to be decimated to a heap.*/
static void od_mv_dec_heapify(od_mv_est_ctx *est) {
  od_mv_node **heap;
  int l;
  int r;
  int i;
  heap = est->dec_heap;
  l = est->dec_nheap >> 1;
  r = est->dec_nheap - 1;
  for (i = l; i-- > 0;) {
    int p;
    p = i;
    do {
      int q;
      q = (p << 1) + 1;
      if (q < r && od_mv_dec_cmp(heap[q], heap[q + 1]) >= 0) q++;
      if (od_mv_dec_cmp(heap[p], heap[q]) <= 0) break;
      od_mv_dec_heap_swap(heap, p, q);
      p = q;
    }
    while (p < l);
  }
}

/*Restore the heap structure at the given index by moving it down the heap.*/
static void od_mv_dec_heap_down(od_mv_est_ctx *est, int heapi) {
  od_mv_node **heap;
  int l;
  int r;
  int p;
  heap = est->dec_heap;
  l = est->dec_nheap >> 1;
  r = est->dec_nheap - 1;
  p = heapi;
  while (p < l) {
    int q;
    q = (p << 1) + 1;
    if (q < r && od_mv_dec_cmp(heap[q], heap[q + 1]) >= 0) q++;
    if (od_mv_dec_cmp(heap[p], heap[q]) <= 0) break;
    od_mv_dec_heap_swap(heap, p, q);
    p = q;
  }
}

/*Restore the heap structure at the given index by moving it up the heap.*/
static void od_mv_dec_heap_up(od_mv_est_ctx *est, int heapi) {
  od_mv_node **heap;
  int p;
  heap = est->dec_heap;
  p = heapi;
  while (p > 0) {
    int q;
    q = p;
    p = ((q + 1) >> 1) - 1;
    if (od_mv_dec_cmp(heap[p], heap[q]) <= 0) break;
    od_mv_dec_heap_swap(heap, p, q);
  }
}

/*Retrieve the item at the top of the heap.
  Returns NULL if there are no more nodes to decimate.*/
static od_mv_node *od_mv_dec_heap_delhead(od_mv_est_ctx *est) {
  od_mv_node *ret;
  if (est->dec_nheap <= 0) return NULL;
  ret = est->dec_heap[0];
  ret->heapi = -1;
  if (--est->dec_nheap > 0) {
    est->dec_heap[0] = est->dec_heap[est->dec_nheap];
    est->dec_heap[0]->heapi = 0;
    od_mv_dec_heap_down(est, 0);
  }
  return ret;
}

static void od_mv_dec_heap_del(od_mv_est_ctx *est, od_mv_node *node) {
  int heapi;
  heapi = node->heapi;
  if (heapi >= 0) {
    node->heapi = -1;
    est->dec_nheap--;
    if (est->dec_nheap > heapi) {
      est->dec_heap[heapi] = est->dec_heap[est->dec_nheap];
      est->dec_heap[heapi]->heapi = heapi;
      if (od_mv_dec_cmp(node, est->dec_heap[heapi]) >= 0) {
        od_mv_dec_heap_up(est, heapi);
      }
      else od_mv_dec_heap_down(est, heapi);
    }
    else est->dec_heap[est->dec_nheap] = NULL;
  }
}

/*Sets the dd and dr values of the given node, restoring the heap structure
   afterwards.*/
static void od_mv_dec_update(od_mv_est_ctx *est, od_mv_node *node,
 int dd, int dr) {
  int diff;
  diff = od_mv_dddr_cmp(dd, dr, node->dd, node->dr);
  node->dd = dd;
  node->dr = dr;
  if (node->heapi >= 0) {
    if (diff <= 0) od_mv_dec_heap_up(est, node->heapi);
    else od_mv_dec_heap_down(est, node->heapi);
  }
}

static void od_mv_est_init_nodes(od_mv_est_ctx *est) {
  od_state *state;
  od_mv_node *mv_row;
  od_mv_grid_pt *grid;
  int nhmvbs;
  int nvmvbs;
  int vx;
  int vy;
  state = &est->enc->state;
  nhmvbs = state->nhmvbs;
  nvmvbs = state->nvmvbs;
  for (vy = 0; vy <= nvmvbs; vy++) {
    mv_row = est->mvs[vy];
    grid = state->mv_grid[vy];
    for (vx = 0; vx <= nhmvbs; vx++) {
      int level;
      level = OD_MC_LEVEL[vy & OD_MVB_MASK][vx & OD_MVB_MASK];
      if (level <= est->level_max) {
        int flag_rate;
        /*While we're here, reset the MV state.*/
        OD_ASSERT(grid[vx].valid == 1);
        est->row_counts[vy]++;
        est->col_counts[vx]++;
        /*Inbetween the level limits, vertices require on average 2 bits to
           indicate the presence of children.
          TODO: Fix a more exact representation.
          TODO: Fix-up child flags for blocks outside the frame border.*/
        flag_rate = (est->level_min <= level && level < est->level_max) <<
         (1 + OD_BITRES);
        mv_row[vx].dr = -mv_row[vx].mv_rate - flag_rate;
      }
      else grid[vx].valid = 0;
    }
  }
}

/*Computes the SAD of all blocks at all scales with all possible edge
   splittings, using OBMC.
  These are what will drive the error of the adaptive subdivision process.*/
static void od_mv_est_calc_sads(od_mv_est_ctx *est) {
  od_state *state;
  int nhmvbs;
  int nvmvbs;
  int level_max;
  int level_min;
  int log_mvb_sz;
  int vx;
  int vy;
  int oc;
  int s;
  state = &est->enc->state;
  /*TODO: Interleaved evaluation would probably provide better cache
     coherency.*/
  nhmvbs = state->nhmvbs;
  nvmvbs = state->nvmvbs;
  level_max = est->level_max;
  level_min = est->level_min;
  for (log_mvb_sz = 0; log_mvb_sz < OD_LOG_MVB_DELTA0; log_mvb_sz++) {
    if (level_max >= OD_MC_LEVEL_MAX - 1 - 2*log_mvb_sz
     && level_min <= OD_MC_LEVEL_MAX - 2*log_mvb_sz) {
      od_sad4 **sad_cache;
      int smax;
      sad_cache = est->sad_cache[log_mvb_sz];
      smax = level_max >= OD_MC_LEVEL_MAX - 2*log_mvb_sz ? 4 : 1;
      for (vy = 0; vy < nvmvbs; vy++) {
        od_sad4 *sad_cache_row;
        od_mv_node *mv_row;
        sad_cache_row = sad_cache[vy];
        mv_row = est->mvs[vy << log_mvb_sz];
        for (vx = 0; vx < nhmvbs; vx++) {
          oc = (vx & 1) ^ ((vy & 1) << 1 | (vy & 1));
          for (s = 0; s < smax; s++) {
            sad_cache_row[vx][s] = od_mv_est_sad(est,
             vx << log_mvb_sz, vy << log_mvb_sz, oc, s, log_mvb_sz);
          }
          /*While we're here, fill in the block's setup state.*/
          if (level_max <= OD_MC_LEVEL_MAX - 2*log_mvb_sz) {
            mv_row[vx << log_mvb_sz].oc = oc;
            mv_row[vx << log_mvb_sz].log_mvb_sz = log_mvb_sz;
            mv_row[vx << log_mvb_sz].s = smax - 1;
            mv_row[vx << log_mvb_sz].sad = sad_cache_row[vx][smax - 1];
          }
        }
      }
    }
    nhmvbs >>= 1;
    nvmvbs >>= 1;
  }
  if (level_max <= 0) {
    for (vy = 0; vy < nvmvbs; vy++) {
      od_mv_node *mv_row;
      mv_row = est->mvs[vy << log_mvb_sz];
      for (vx = 0; vx < nhmvbs; vx++) {
        mv_row[vx << log_mvb_sz].oc = 0;
        mv_row[vx << log_mvb_sz].s = 3;
        mv_row[vx << log_mvb_sz].log_mvb_sz = log_mvb_sz;
        mv_row[vx << log_mvb_sz].sad = od_mv_est_sad(est,
         vx << OD_LOG_MVBSIZE_MIN, vy << OD_LOG_MVBSIZE_MIN, 0, 3, log_mvb_sz);
      }
    }
  }
}

static void od_mv_est_init_du(od_mv_est_ctx *est, int vx, int vy) {
  od_state *state;
  od_mv_node *dec;
  od_mv_node *merge;
  const od_mv_err_node *errdom;
  int nerrdom;
  const od_offset *mergedom;
  int nhmvbs;
  int nvmvbs;
  int level;
  int dlev;
  int log_mvb_sz_min;
  int log_mvb_sz;
  int di;
  int dvx;
  int dvy;
  int dx;
  int dy;
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
   "Computing du's for (%i, %i)", vx, vy));
  state = &est->enc->state;
  nhmvbs = state->nhmvbs;
  nvmvbs = state->nvmvbs;
  dec = est->mvs[vy] + vx;
  level = OD_MC_LEVEL[vy & OD_MVB_MASK][vx & OD_MVB_MASK];
  dlev = (OD_MC_LEVEL_MAX - est->level_max) >> 1;
  log_mvb_sz_min = (OD_MC_LEVEL_MAX - est->level_max + 1) >> 1;
  errdom = OD_ERRDOM[level - 1 + (dlev << 1)];
  nerrdom = OD_NERRDOM[level - 1 + (dlev << 1)];
  mergedom = OD_MERGEDOM[level - 1 + (dlev << 1)];
  dec->dd = 0;
  /*Subtract off the error before decimation.*/
  for (di = 0; di < nerrdom; di++) {
    dvx = vx + (errdom[di].dx << dlev);
    dvy = vy + (errdom[di].dy << dlev);
    if (dvx >= 0 && dvy >= 0 && dvx < nhmvbs && dvy < nvmvbs) {
      int mvb_sz;
      log_mvb_sz = errdom[di].log_mvb_sz + dlev;
      if (log_mvb_sz < log_mvb_sz_min) continue;
      mvb_sz = 1 << (log_mvb_sz - dlev);
      for (dy = 0; dy < mvb_sz; dy++) {
        for (dx = 0; dx < mvb_sz; dx++) {
          dec->dd -= est->mvs[dvy + (dy << dlev)][dvx + (dx << dlev)].sad;
          /*dec->dd -= est->sad_cache[dlev][(dvy >> dlev) + dy]
           [(dvx >> dlev) + dx][undecs];*/
          OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
           "Added error (%i, %i) [%ix%i]: %i", dvx + (dx << dlev),
           dvy + (dy << dlev), OD_MVBSIZE_MIN << dlev, OD_MVBSIZE_MIN << dlev,
           dec->dd));
        }
      }
    }
    else {
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
       "(%i, %i) outside [%i, %i]x[%i, %i]",
       dvx, dvy, 0, 0, nhmvbs, nvmvbs));
    }
  }
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
   "Subtracted initial error: %i", dec->dd));
  /*Decimate the vertices in the merging domain.
    Also sum up the rate changes while we do it.*/
  for (di = 0;; di++) {
    dvx = vx + (mergedom[di][0] << dlev);
    if (dvx < 0 || dvx > nhmvbs) continue;
    dvy = vy + (mergedom[di][1] << dlev);
    if (dvy < 0 || dvy > nvmvbs) continue;
    if (OD_MC_LEVEL[dvy & OD_MVB_MASK][dvx & OD_MVB_MASK] > est->level_max) {
      continue;
    }
    state->mv_grid[dvy][dvx].valid = 0;
    merge = est->mvs[dvy] + dvx;
    if (merge == dec) break;
    dec->dr += merge->dr;
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
     "Merged vertex (%2i, %2i), dr: %i", dvx, dvy, dec->dr));
  }
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
   "Merged vertex (%2i, %2i)", dvx, dvy));
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
   "Decimated vertices in merging domain."));
  /*Add in the error after decimation.*/
  for (di = 0; di < nerrdom; di++) {
    dvx = vx + (errdom[di].dx << dlev);
    dvy = vy + (errdom[di].dy << dlev);
    if (dvx >= 0 && dvy >= 0 && dvx < nhmvbs && dvy < nvmvbs) {
      log_mvb_sz = errdom[di].log_mvb_sz + dlev;
      if (log_mvb_sz < log_mvb_sz_min) continue;
      else if (log_mvb_sz < OD_LOG_MVB_DELTA0) {
        int mask;
        int s1vx;
        int s1vy;
        int s3vx;
        int s3vy;
        int oc;
        int s;
        mask = (1 << (log_mvb_sz + 1)) - 1;
        oc = !!(dvx & mask);
        if (dvy & mask) oc = 3 - oc;
        s1vx = dvx + (OD_VERT_DX[(oc + 1) & 3] << log_mvb_sz);
        s1vy = dvy + (OD_VERT_DY[(oc + 1) & 3] << log_mvb_sz);
        s3vx = dvx + (OD_VERT_DX[(oc + 3) & 3] << log_mvb_sz);
        s3vy = dvy + (OD_VERT_DY[(oc + 3) & 3] << log_mvb_sz);
        s = state->mv_grid[s1vy][s1vx].valid |
         state->mv_grid[s3vy][s3vx].valid << 1;
        dec->dd +=
         est->sad_cache[log_mvb_sz][dvy >> log_mvb_sz][dvx >> log_mvb_sz][s];
        OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
         "Added error (%i, %i) [%ix%i] {%i, %i}: %i", dvx, dvy,
         OD_MVBSIZE_MIN << log_mvb_sz, OD_MVBSIZE_MIN << log_mvb_sz,
         oc, s, dec->dd));
      }
      else {
        /*Cache the SAD for top-level blocks in the dd field, which is
           otherwise unused (since they cannot be decimated).*/
        est->mvs[dvy][dvx].dd = od_mv_est_sad(est,
         dvx, dvy, 0, 3, OD_LOG_MVB_DELTA0);
        dec->dd += est->mvs[dvy][dvx].dd;
        OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
         "Added error (%i, %i) [%ix%i]: %i", dvx, dvy,
         OD_MVBSIZE_MIN << log_mvb_sz, OD_MVBSIZE_MIN << log_mvb_sz, dec->dd));
      }
    }
  }
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
   "Total merging error: %i", dec->dd));
  /*Restore the vertices in the merging domain.*/
  for (di = 0;; di++) {
    dvx = vx + (mergedom[di][0] << dlev);
    if (dvx < 0 || dvx > nhmvbs) continue;
    dvy = vy + (mergedom[di][1] << dlev);
    if (dvy < 0 || dvy > nvmvbs) continue;
    if (OD_MC_LEVEL[dvy & OD_MVB_MASK][dvx & OD_MVB_MASK] > est->level_max) {
      continue;
    }
    state->mv_grid[dvy][dvx].valid = 1;
    if (dvx == vx && dvy == vy) break;
  }
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
   "Restored vertices in merging domain."));
  /*Add this node to the heap.*/
  dec->heapi = est->dec_nheap;
  est->dec_heap[est->dec_nheap++] = dec;
}

static void od_mv_est_init_dus(od_mv_est_ctx *est) {
  od_state *state;
  int nhmvbs;
  int nvmvbs;
  int log_mvb_sz;
  int level;
  int vx; 
  int vy;
  state = &est->enc->state;
  nhmvbs = state->nhmvbs;
  nvmvbs = state->nvmvbs;
  OD_CLEAR(est->row_counts, nvmvbs + 1);
  OD_CLEAR(est->col_counts, nhmvbs + 1);
  od_mv_est_init_nodes(est);
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG, "Finished MV bits."));
  od_mv_est_calc_sads(est);
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG, "Finished SADs."));
  /*Clear the merge heap.*/
  est->dec_nheap = 0;
  est->dec_heap[0] = NULL;
  /*The initialization is destructive to dr, and so must proceed by level from
     top to bottom.*/
  for (log_mvb_sz = OD_LOG_MVB_DELTA0, level = 1; log_mvb_sz-- > 0; level++) {
    int mvb_sz;
    mvb_sz = 1 << log_mvb_sz;
    /*Odd-level vertices.*/
    if (est->level_max < level) break;
    if (est->level_min < level) {
      for (vy = mvb_sz; vy <= nvmvbs; vy += 2*mvb_sz) {
        for (vx = mvb_sz; vx <= nhmvbs; vx += 2*mvb_sz) {
          od_mv_est_init_du(est, vx, vy);
        }
      }
    }
    /*Even-level vertices.*/
    level++;
    if (est->level_max < level) break;
    if (est->level_min < level) {
      for (vy = 0;; vy += mvb_sz) {
        for (vx = mvb_sz; vx <= nhmvbs; vx += 2*mvb_sz) {
          od_mv_est_init_du(est, vx, vy);
        }
        vy += mvb_sz;
        if (vy > nvmvbs) break;
        for (vx = 0; vx <= nhmvbs; vx += 2*mvb_sz) {
          od_mv_est_init_du(est, vx, vy);
        }
      }
    }
  }
  /*Make the node list into a proper heap.*/
  od_mv_dec_heapify(est);
}

static void od_mv_est_decimate(od_mv_est_ctx *est) {
  od_mv_node *dec;
  od_state *state;
  int nhmvbs;
  int nvmvbs;
  int dlev;
  int vx;
  int vy;
  od_mv_est_init_dus(est);
#if defined(OD_ENABLE_ASSERTIONS) || defined(OD_LOGGING_ENABLED)
  od_mv_est_check_rd_state(est, 2);
#endif
  state = &est->enc->state;
  nhmvbs = state->nhmvbs;
  nvmvbs = state->nvmvbs;
  dlev = (OD_MC_LEVEL_MAX - est->level_max) >> 1;
  for (;;) {
    const od_offset *mergedom;
    int level;
    int di;
#if defined(OD_DUMP_IMAGES) && defined(OD_ANIMATE)
    if (daala_granule_basetime(state, state->cur_time) == ANI_FRAME) {
      char iter_label[16];
      od_state_mc_predict(state, state->ref_imgs
       + state->ref_imgi[OD_FRAME_SELF]);
      od_encode_fill_vis(est->enc);
      sprintf(iter_label, "ani%08i", est->enc->ani_iter++);
      od_state_dump_img(state, &est->enc->vis_img, iter_label);
    }
#endif
    dec = od_mv_dec_heap_delhead(est);
    /*Stop if we've fully decimated the mesh, or if this decimation would not
       improve R-D performance at the current lambda.*/
    if (dec == NULL
     || dec->dr*est->lambda + (dec->dd << OD_ERROR_SCALE) > 0) {
      break;
    }
    level = OD_MC_LEVEL[dec->vy & OD_MVB_MASK][dec->vx & OD_MVB_MASK];
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
     "Merging node (%2i, %2i), level %i, dd %5i, dr %5i, dopt %5i:",
     dec->vx, dec->vy, level, dec->dd, dec->dr,
     dec->dr*est->lambda + (dec->dd << OD_ERROR_SCALE)));
    mergedom = OD_MERGEDOM[level - 1 + (dlev << 1)];
    for (di = 0;; di++) {
      od_mv_node *merge;
      od_mv_node *ancestor;
      od_mv_node *block;
      const od_offset *anc;
      int nanc;
      int ai;
      int ax;
      int ay;
      int bx;
      int by;
      int log_mvb_sz;
      int mask;
      /*Don't decimate vertices outside of the mesh.*/
      vx = dec->vx + (mergedom[di][0] << dlev);
      if (vx < 0 || vx > nhmvbs) continue;
      vy = dec->vy + (mergedom[di][1] << dlev);
      if (vy < 0 || vy > nvmvbs) continue;
      merge = est->mvs[vy] + vx;
      /*Don't decimate vertices that have already been decimated.*/
      if (!state->mv_grid[vy][vx].valid) {
        OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
         "Skipping node (%i, %i) (already merged).", vx, vy));
        continue;
      }
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
       "Merging node (%2i, %2i), dd %5i, dr %5i:",
       vx, vy, merge->dd, merge->dr));
      /*Update the deltas for this vertex in the merging domain.
        The simple rule applied below handles overlapped domains with an
         inclusion-exclusion approach.
        See Balmelli 2001 for details.*/
      nanc = OD_NANCESTORS[vy & OD_MVB_MASK][vx & OD_MVB_MASK];
      anc = OD_ANCESTORS[vy & OD_MVB_MASK][vx & OD_MVB_MASK];
      for (ai = 0; ai < nanc; ai++) {
        ax = vx + anc[ai][0];
        if (ax < 0 || ax > nhmvbs) continue;
        ay = vy + anc[ai][1];
        if (ay < 0 || ay > nvmvbs) continue;
        ancestor = est->mvs[ay] + ax;
        od_mv_dec_update(est, ancestor,
         ancestor->dd - merge->dd, ancestor->dr - merge->dr);
        OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
         "Updated ancestor (%2i, %2i) of (%2i, %2i): dd %5i, dr %5i",
         ax, ay, vx, vy, ancestor->dd, ancestor->dr));
      }
      state->mv_grid[vy][vx].valid = 0;
      od_mv_dec_heap_del(est, merge);
      est->row_counts[vy]--;
      est->col_counts[vx]--;
      level = OD_MC_LEVEL[vy & OD_MVB_MASK][vx & OD_MVB_MASK];
      log_mvb_sz = (OD_MC_LEVEL_MAX - level) >> 1;
      /*Account for quadrilaterals which may have only partially belonged to
         the merging domain (e.g., that would not have belonged were we using
         triangles).*/
      if (!(level & 1)) {
        static const int OD_CDX[4] = { -1, 1, -1, 1 };
        static const int OD_CDY[4] = { -1, -1, 1, 1 };
        int k;
        mask = (1 << (log_mvb_sz + 1)) - 1;
        for (k = 0; k < 4; k++) {
          int32_t ddd;
          int cx;
          int cy;
          int s;
          cx = vx + (OD_CDX[k] << log_mvb_sz);
          if (cx < 0 || cx > nhmvbs) continue;
          cy = vy + (OD_CDY[k] << log_mvb_sz);
          if (cy < 0 || cy > nvmvbs) continue;
          bx = vx + (OD_ERRDOM6[k].dx << log_mvb_sz);
          by = vy + (OD_ERRDOM6[k].dy << log_mvb_sz);
          block = est->mvs[by] + bx;
          by >>= log_mvb_sz;
          bx >>= log_mvb_sz;
          if (!state->mv_grid[cy][cx].valid) {
            block->s = 0;
            block->sad = est->sad_cache[log_mvb_sz][by][bx][0];
            /*If the opposing corner has already been decimated, the remaining
               adjustments have already been made.*/
            continue;
          }
          /*s is the split state of the error block with (vx, vy) decimated,
             and (cx, cy) undecimated.*/
          s = 1 << ((((k + 3) & 3) >> 1) ^ !!(vx & mask));
          block->s = s;
          block->sad = est->sad_cache[log_mvb_sz][by][bx][s];
          /*Replace the old decimation error change with the new one.*/
          ddd = est->sad_cache[log_mvb_sz][by][bx][0]
           - est->sad_cache[log_mvb_sz][by][bx][s ^ 3]
           + est->sad_cache[log_mvb_sz][by][bx][3]
           - est->sad_cache[log_mvb_sz][by][bx][s];
          OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
           "Checking opposing corner (%2i, %2i): ddd %i", cx, cy, ddd));
          /*This happens in regions of constant motion.*/
          if (ddd == 0) continue;
          ancestor = est->mvs[cy] + cx;
          od_mv_dec_update(est, ancestor, ancestor->dd + ddd, ancestor->dr);
          OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
           "Updated corner (%2i, %2i): dd %5i, dr %5i",
           cx, cy, ancestor->dd, ancestor->dr));
          /*Update the opposing corner's ancestors, which also, of
             necessity, must contain the affected quadrilateral, and must
             not have been decimated yet.*/
          nanc = OD_NANCESTORS[cy & OD_MVB_MASK][cx & OD_MVB_MASK];
          anc = OD_ANCESTORS[cy & OD_MVB_MASK][cx & OD_MVB_MASK];
          for (ai = 0; ai < nanc; ai++) {
            ax = cx + anc[ai][0];
            if (ax < 0 || ax > nhmvbs) continue;
            ay = cy + anc[ai][1];
            if (ay < 0 || ay > nvmvbs) continue;
            ancestor = est->mvs[ay] + ax;
            od_mv_dec_update(est, ancestor, ancestor->dd + ddd, ancestor->dr);
            OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
             "Updated ancestor (%2i, %2i): dd %5i, dr %5i",
             ax, ay, ancestor->dd, ancestor->dr));
          }
          /*Add back in the components that do not apply to the interior
             corner.*/
          ddd = -ddd;
          if (vx & mask) cx = vx;
          else cy = vy;
          OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
           "Checking interior corner (%2i, %2i): ddd %i", cx, cy, ddd));
          ancestor = est->mvs[cy] + cx;
          od_mv_dec_update(est, ancestor, ancestor->dd + ddd, ancestor->dr);
          OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
           "Updated corner (%2i, %2i): dd %5i, dr %5i",
           cx, cy, ancestor->dd, ancestor->dr));
          /*And update all the interior corner's ancestors, which also, of
             necessity, must contain the affected quadrilateral, and must not
             have been decimated yet.*/
          nanc = OD_NANCESTORS[cy & OD_MVB_MASK][cx & OD_MVB_MASK];
          anc = OD_ANCESTORS[cy & OD_MVB_MASK][cx & OD_MVB_MASK];
          for (ai = 0; ai < nanc; ai++) {
            ax = cx + anc[ai][0];
            if (ax < 0 || ax > nhmvbs) continue;
            ay = cy + anc[ai][1];
            if (ay < 0 || ay > nvmvbs) continue;
            ancestor = est->mvs[ay] + ax;
            od_mv_dec_update(est, ancestor, ancestor->dd + ddd, ancestor->dr);
            OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
             "Updated ancestor (%2i, %2i): dd %5i, dr %5i",
             ax, ay, ancestor->dd, ancestor->dr));
          }
        }
      }
      /*Otherwise, we eliminated several smaller blocks.
        Update the SAD and block setup for the larger block that took their
         place.*/
      else {
        int oc;
        bx = vx - (1 << log_mvb_sz);
        by = vy - (1 << log_mvb_sz);
        log_mvb_sz++;
        mask = (1 << (log_mvb_sz + 1)) - 1;
        oc = !!(bx & mask);
        if (by & mask) oc = 3 - oc;
        block = est->mvs[by] + bx;
        block->log_mvb_sz = log_mvb_sz;
        block->oc = oc;
        block->s = 3;
        if (log_mvb_sz < OD_LOG_MVB_DELTA0) {
          block->sad =
           est->sad_cache[log_mvb_sz][by >> log_mvb_sz][bx >> log_mvb_sz][3];
        }
        /*At the top level, we cached the SAD in the dd field.*/
        else block->sad = block->dd;
      }
      /*If we just decimated our target vertex, stop.*/
      if (merge == dec) break;
    }
  }
#if defined(OD_ENABLE_ASSERTIONS) || defined(OD_LOGGING_ENABLED)
  od_mv_est_check_rd_state(est, 2);
#endif
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG, "Finished merging."));
  /*if (dec != NULL) {
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
     "Node (%i, %i) dd %i, dr %i, dopt %i: not enough.", dec->vx, dec->vy,
     dec->dd, dec->dr, dec->dr*est->lambda + (dec->dd << OD_ERROR_SCALE)));
  }*/
  /*if (state->mv_grid[31][1].valid) {
    dec = est->mvs[31] + 1;
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
     "(%i, %i) remains. dd: %5i, dr: %2i, dopt: %6i.", dec->vx, dec->vy,
     dec->dd, dec->dr, dec->dr*est->lambda + (dec->dd << OD_ERROR_SCALE)));
  }*/
}

/*STAGE 3: Iterated Dynamic Programming.*/

/*The list of MVs that can be predicted by a level 0 MV, excluding those not
   yet considered by DP across rows.
    +-------+-------+
    |\      |      /
    | \     |     /
    |  \    |    /
    |   1---2---1
    |   |\ /|\ /
    |   | 3-4-3
    |   |/|565
    +---2-46*.. .   .
    |   |\|565
    |   | 3-4-3
    |   |/ \|/ \
    |   1---2---1
    |  /    |    \
    | /     |     \
    |/      |      \
    0-------0-------0 */
static const od_offset OD_ROW_PREDICTED0[24] = {
  /*These predicted MVs are changeable by future MVs in the DP path.*/
  {  4, -4 }, {  2, -2 }, {  1, -1 }, {  1,  1 }, {  2,  2 }, {  4,  4 },
  {  0,  8 }, {  8,  8 },
  /*The remaining ones are not.*/
  { -4, -4 }, {  0, -4 }, { -2, -2 }, {  0, -2 }, { -1, -1 }, {  0, -1 },
  { -4,  0 }, { -2,  0 }, { -1,  0 }, { -1,  1 }, {  0,  1 }, { -2,  2 },
  {  0,  2 }, { -4,  4 }, {  0,  4 }, { -8,  8 }
};
/*The list of MVs that can be predicted by a level 1 MV, excluding those
   not yet considered by DP across rows.
        +---2---+
        |\ /|\ /
        | 3-4-3
        |/|565
        2-46*.. .
        |\|565
        | 3-4-3
        |/ \|/ \
        +---2---+ */
static const od_offset OD_ROW_PREDICTED1[17] = {
  /*These predicted MVs are changeable by future MVs in the DP path.*/
  {  2, -2 }, {  1, -1 }, {  2,  2 }, {  1,  1 },
  /*The remaining ones are not.*/
  {  0, -4 }, { -2, -2 }, {  0, -2 }, { -1, -1 }, {  0, -1 }, { -4,  0 },
  { -2,  0 }, { -1,  0 }, { -1,  1 }, {  0,  1 }, { -2,  2 }, {  0,  2 },
  {  0,  4 }
};
/*The list of MVs that can be predicted by a level 2 MV, excluding those
   not yet considered by DP across rows.
          3-4-3
          |565
          46*..
          |565
          3-4-3 */
static const od_offset OD_ROW_PREDICTED2[14] = {
  /*These predicted MVs are changeable by future MVs in the DP path.*/
  {  2, -2 }, {  1, -1 }, {  1,  1 }, {  2,  2 },
  /*The remaining ones are not.*/
  { -2, -2 }, {  0, -2 }, { -1, -1 }, {  0, -1 }, { -2,  0 }, { -1,  0 },
  { -1,  1 }, {  0,  1 }, { -2,  2 }, {  0,  2 }
};
/*The list of MVs that can be predicted by a level 3 MV, excluding those
   not yet considered by DP across rows.
          +-4-+
          |565
          46*..
          |565
          +-4-+ */
static const od_offset OD_ROW_PREDICTED3[10] = {
  /*These predicted MVs are changeable by future MVs in the DP path.*/
  {  1, -1 }, {  1,  1 },
  /*The remaining ones are not.*/
  {  0, -2 }, { -1, -1 }, {  0, -1 }, { -2,  0 }, { -1,  0 }, { -1,  1 },
  {  0,  1 }, {  0,  2 }
};
/*The list of MVs that can be predicted by a level 4 MV, excluding those
   not yet considered by DP across rows.
           565
           6*.
           565 */
static const od_offset OD_ROW_PREDICTED4[7] = {
  /*These predicted MVs are changeable by future MVs in the DP path.*/
  {  1, -1 }, {  1,  1 },
  /*The remaining ones are not.*/
  { -1, -1 }, {  0, -1 }, { -1,  0 }, { -1,  1 }, {  0,  1 }
};
/*The list of MVs that can be predicted by a level 5 MV, excluding those
   not yet considered by DP across rows.
           +6+
           6*.
           +6+ */
static const od_offset OD_ROW_PREDICTED5[3] = {
  /*These predicted MVs are NOT changeable by future MVs in the DP path.*/
  {  0, -1 }, { -1,  0 }, {  0,  1 }
};

/*The list of MVs that can be predicted by a level 0 MV, excluding those not
   yet considered by DP across columns.
    +-------+-------+
    |\      |      /|
    | \     |     / |
    |  \    |    /  |
    |   1---2---1   |
    |   |\ /|\ /|   |
    |   | 3-4-3 |   |
    |   |/|565|\|   |
    +---2-46*64-2---0
    |   |\|5.5|/|   |
    |   | 3 . 3 |   |
    |   |/     \|   |
    |   1   .   1   |
    |  /         \  |
    | /           \ |
    |/             \|
    0       .       0 */
static const od_offset OD_COL_PREDICTED0[24] = {
  /*These predicted MVs are changeable by future MVs in the DP path.*/
  { -1,  1 }, {  1,  1 }, { -2,  2 }, {  2,  2 }, { -4,  4 }, {  4,  4 },
  {  8,  8 },
  /*The remaining ones are not.*/
  { -4, -4 }, {  0, -4 }, {  4, -4 }, { -2, -2 }, {  0, -2 }, {  2, -2 },
  { -1, -1 }, {  0, -1 }, {  1, -1 }, { -4,  0 }, { -2,  0 }, { -1,  0 },
  {  1,  0 }, {  2,  0 }, {  4,  0 }, {  8,  0 }, { -8,  8 }
};
/*The list of MVs that can be predicted by a level 1 MV, excluding those
   not yet considered by DP across columns.
        +---2---+
        |\ /|\ /|
        | 3-4-3 |
        |/|565|\|
        2-46*64-2
        |\|5.5|/|
        | 3 . 3 |
        |/     \|
        +   .   + */
static const od_offset OD_COL_PREDICTED1[17] = {
  /*These predicted MVs are changeable by future MVs in the DP path.*/
  { -1,  1 }, {  1,  1 }, { -2,  2 }, {  2,  2 },
  /*The remaining ones are not.*/
  {  0, -4 }, { -2, -2 }, {  0, -2 }, {  2, -2 }, { -1, -1 }, {  0, -1 },
  {  1, -1 }, { -4,  0 }, { -2,  0 }, { -1,  0 }, {  1,  0 }, {  2,  0 },
  {  4,  0 }
};
/*The list of MVs that can be predicted by a level 2 MV, excluding those
   not yet considered by DP across columns.
          3-4-3
          |565|
          46*64
          |5.5|
          3 . 3 */
static const od_offset OD_COL_PREDICTED2[14] = {
  /*These predicted MVs are changeable by future MVs in the DP path.*/
  { -1,  1 }, {  1,  1 }, { -2,  2 }, {  2,  2 },
  /*The remaining ones are not.*/
  { -2, -2 }, {  0, -2 }, {  2, -2 }, { -1, -1 }, {  0, -1 }, {  1, -1 },
  { -2,  0 }, { -1,  0 }, {  1,  0 }, {  2,  0 }
};
/*The list of MVs that can be predicted by a level 3 MV, excluding those
   not yet considered by DP across columns.
          +-4-+
          |565|
          46*64
          |5.5|
          + . + */
static const od_offset OD_COL_PREDICTED3[10] = {
  /*These predicted MVs are changeable by future MVs in the DP path.*/
  { -1,  1 }, {  1,  1 },
  /*The remaining ones are not.*/
  {  0, -2 }, { -1, -1 }, {  0, -1 }, {  1, -1 }, { -2,  0 }, { -1,  0 },
  {  1,  0 }, {  2,  0 }
};
/*The list of MVs that can be predicted by a level 4 MV, excluding those
   not yet considered by DP across columns.
           565
           6*6
           5.5 */
static const od_offset OD_COL_PREDICTED4[7] = {
  /*These predicted MVs are changeable by future MVs in the DP path.*/
  { -1,  1 }, {  1,  1 },
  /*The remaining ones are not.*/
  { -1, -1 }, {  0, -1 }, {  1, -1 }, { -1,  0 }, {  1,  0 }
};
/*The list of MVs that can be predicted by a level 5 MV, excluding those
   not yet considered by DP across columns.
           +6+
           6*6
           +.+ */
static const od_offset OD_COL_PREDICTED5[3] = {
  /*These predicted MVs are NOT changeable by future MVs in the DP path.*/
  {  0, -1 }, { -1,  0 }, {  1,  0 }
};

/*The number of predicted MVs in each list.*/
static const int OD_NPREDICTED[OD_MC_NLEVELS] = { 24, 17, 14, 10, 7, 3, 0 };
/*The number of changeable predicted MVs in each list.*/
static const int OD_NROW_PRED_CHANGEABLE[OD_MC_LEVEL_MAX] = {
  8, 4, 4, 2, 2, 0
};
/*The number of changeable predicted MVs in each list.*/
static const int OD_NCOL_PRED_CHANGEABLE[OD_MC_LEVEL_MAX] = {
  7, 4, 4, 2, 2, 0
};
/*The lists of offsets to predicted MVs for each level.*/
static const od_offset *const OD_ROW_PREDICTED[OD_MC_LEVEL_MAX] = {
  OD_ROW_PREDICTED0,
  OD_ROW_PREDICTED1,
  OD_ROW_PREDICTED2,
  OD_ROW_PREDICTED3,
  OD_ROW_PREDICTED4,
  OD_ROW_PREDICTED5
};
/*The lists of offsets to predicted MVs for each level.*/
static const od_offset *const OD_COL_PREDICTED[OD_MC_LEVEL_MAX] = {
  OD_COL_PREDICTED0,
  OD_COL_PREDICTED1,
  OD_COL_PREDICTED2,
  OD_COL_PREDICTED3,
  OD_COL_PREDICTED4,
  OD_COL_PREDICTED5
};

/*This should be the maximum value found in either OD_ROW_PRED_HIST_SIZE or
   OD_COL_PRED_HIST_SIZE.*/
#define OD_PRED_HIST_SIZE_MAX (16)

/*The amount of history to restore in the trellis state to ensure predicted MVs
   are evaluated correctly in row refinement.*/
static const int OD_ROW_PRED_HIST_SIZE[OD_MC_NLEVELS] = {
  16, 8, 4, 4, 2, 2, 1
};
/*The amount of history to restore in the trellis state to ensure predicted MVs
   are evaluated correctly in column refinement.*/
static const int OD_COL_PRED_HIST_SIZE[OD_MC_NLEVELS] = {
  8, 8, 4, 4, 2, 2, 1
};

/*Returns the boundary case indicating which motion vector range edges the
   current motion vector is abutting.
  vx: The horizontal position of the node.
  vy: The vertical position of the node.
  dx: The horizontal component of the motion vector.
  dy: The vertical component of the motion vector.
  dsz: The amount the vector is being adjusted by.
  log_blk_sz: The log base 2 of the maximum size of a block the vector can
               belong to.
  Return: A set of flags indicating the boundary conditions, after the
   documentation at OD_SQUARE_SITES.*/
static int od_mv_est_get_boundary_case(od_state *state,
 int vx, int vy, int dx, int dy, int dsz, int log_blk_sz) {
  int bxmin;
  int bymin;
  int bxmax;
  int bymax;
  int mvxmin;
  int mvxmax;
  int mvymin;
  int mvymax;
  int blk_sz;
  int bx;
  int by;
  blk_sz = 1 << log_blk_sz;
  bx = vx << OD_LOG_MVBSIZE_MIN;
  by = vy << OD_LOG_MVBSIZE_MIN;
  /*Each grid point can contribute to 4 blocks around it at the current block
     size, _except_ at the border of the frame, since we do not construct a
     prediction for blocks completely outside the frame.
    For simplicity, we do not try to take into account that further splitting
     might restrict the area this MV influences, even though we now know what
     that splitting is.
    So this MV can affect a block of pixels bounded by
     [bxmin, bmax) x [bymin, bymax), and MVs within that area must point no
     farther than OD_UMV_PADDING pixels outside of the frame.*/
  bxmin = OD_MAXI(bx - blk_sz, 0);
  mvxmin = (OD_MAXI(bxmin - OD_MC_SEARCH_RANGE, -OD_UMV_CLAMP) - bxmin) << 3;
  bxmax = OD_MINI(bx + blk_sz, state->frame_width);
  mvxmax = (OD_MINI(bxmax + OD_MC_SEARCH_RANGE - 1,
   state->frame_width + OD_UMV_CLAMP) - bxmax) << 3;
  bymin = OD_MAXI(by - blk_sz, 0);
  mvymin = (OD_MAXI(bymin - OD_MC_SEARCH_RANGE, -OD_UMV_CLAMP) - bymin) << 3;
  bymax = OD_MINI(by + blk_sz, state->frame_height);
  mvymax = (OD_MINI(bymax + OD_MC_SEARCH_RANGE - 1,
   state->frame_height + OD_UMV_CLAMP) - bymax) << 3;
  return (dx - dsz < mvxmin) | (dx + dsz > mvxmax) << 1 |
   (dy - dsz < mvymin) << 2 | (dy  + dsz > mvymax) << 3;
}

/*Computes the SAD of the specified block.*/
static int32_t od_mv_est_block_sad(od_mv_est_ctx *est,
 od_mv_node *block) {
  int32_t ret;
  ret = od_mv_est_sad(est, block->vx, block->vy,
   block->oc, block->s, block->log_mvb_sz);
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
   "Adding SAD (%3i, %3i) [%2ix%2i]: %6i", block->vx << OD_LOG_MVBSIZE_MIN,
   block->vy << OD_LOG_MVBSIZE_MIN, OD_MVBSIZE_MIN << block->log_mvb_sz,
   OD_MVBSIZE_MIN << block->log_mvb_sz, ret));
  return ret;
}

/*Gets the change in SAD for the blocks affected by the given DP node, using
   the current state of the grid.*/
static int32_t od_mv_dp_get_sad_change(od_mv_est_ctx *est,
 od_mv_dp_node *dp, int32_t block_sads[OD_DP_NBLOCKS_MAX]) {
  int bi;
  int32_t dd;
  dd = 0;
  for (bi = 0; bi < dp->nblocks; bi++) {
    od_mv_node *block;
    block = dp->blocks[bi];
    block_sads[bi] = od_mv_est_block_sad(est, block);
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
     "SAD change for block (%i, %i) [%ix%i]: %i - %i = %i", block->vx,
     block->vy, 1 << (block->log_mvb_sz + OD_LOG_MVBSIZE_MIN),
     1 << (block->log_mvb_sz + OD_LOG_MVBSIZE_MIN),
     block_sads[bi], block->sad, block_sads[bi] - block->sad));
    dd += block_sads[bi] - block->sad;
  }
  return dd;
}

/*Computes a rate adjustment for the predictors changed by following the given
   trellis path.
  As a side effect, enough of the trellis needed to evaluate that change is
   loaded into the MV grid.
  dp: The current DP node.
  cur_mv_rate: Returns the MV rate for the motion vector set by the current DP
                node.
  pred_mv_rates: Returns the MV rate for the motion vectors predicted by the
                  MV set by the current DP node.
  prevsi: The state index to follow in the previous DP node.
  mv_res: The motion vector resolution (0 = 1/8th pel to 2 = 1/2 pel).
  Return: The change in rate for the preceding MVs.*/
static int od_mv_dp_get_rate_change(od_mv_est_ctx *est, od_mv_dp_node *dp,
 int *cur_mv_rate, int pred_mv_rates[OD_DP_NPREDICTED_MAX],
 int prevsi, int mv_res) {
  od_mv_node *mv;
  od_mv_grid_pt *mvg;
  int pi;
  int dr;
  /*Move the state from the current trellis path into the grid.*/
  if (dp->min_predictor_node != NULL) {
    int pred_sis[OD_PRED_HIST_SIZE_MAX];
    int pred_si;
    int npreds;
    od_mv_dp_node *pred_dp;
    npreds = dp - dp->min_predictor_node;
    OD_ASSERT2(npreds <= OD_PRED_HIST_SIZE_MAX, "Too far back!");
    OD_LOG_PARTIAL((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG, "Restoring "));
    /*First, follow the trellis path backwards to find the state used in each
       node.*/
    pred_si = pred_sis[npreds - 1] = prevsi;
    for (pi = 2; pi <= npreds; pi++) {
      pred_dp = dp - pi;
      pred_si = pred_sis[npreds - pi] = pred_dp[1].states[pred_si].prevsi;
    }
    /*Then restore that state going FORWARDS.*/
    for (pred_dp = dp->min_predictor_node; pred_dp < dp; pred_dp++) {
      pred_si = pred_sis[pred_dp - dp->min_predictor_node];
      /*Restore the state for this MV itself.*/
      pred_dp->mv->mv_rate = pred_dp->states[pred_si].mv_rate;
      mvg = pred_dp->mvg;
      mvg->mv[0] = pred_dp->states[pred_si].mv[0];
      mvg->mv[1] = pred_dp->states[pred_si].mv[1];
      OD_LOG_PARTIAL((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
       "(%i, %i: %i)->(%i , %i) ",
       pred_dp->mv->vx, pred_dp->mv->vy, pred_si, mvg->mv[0], mvg->mv[1]));
      /*Restore the state for the MVs this one predicted.*/
      for (pi = 0; pi < pred_dp->npred_changeable; pi++) {
        pred_dp->predicted_mvs[pi]->mv_rate =
         pred_dp->states[pred_si].pred_mv_rates[pi];
      }
    }
    OD_LOG_PARTIAL((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG, "\n"));
  }
  /*Compute the new rate for the current MV.*/
  mv = dp->mv;
  mvg = dp->mvg;
  *cur_mv_rate = od_mv_est_bits(est, mv->vx, mv->vy, mv_res);
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
   "Current MV rate: %i - %i = %i",
   *cur_mv_rate, mv->mv_rate, *cur_mv_rate - mv->mv_rate));
  dr = *cur_mv_rate - mv->mv_rate;
  /*Compute the new rates for the MVs this one predicts.*/
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
   "Calculating predicted pred_mv_rates for node (%i, %i):",
   dp->mv->vx, dp->mv->vy));
  for (pi = 0; pi < dp->npredicted; pi++) {
    mv = dp->predicted_mvs[pi];
    mvg = dp->predicted_mvgs[pi];
    pred_mv_rates[pi] = od_mv_est_bits(est, mv->vx, mv->vy, mv_res);
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
     "Calculated predicted mv_rate of %i for (%i, %i)",
     pred_mv_rates[pi], mv->vx, mv->vy));
    OD_MV_EST_LOG_PRED(est, mv->vx, mv->vy, mv_res);
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
     "MV was: (%i, %i)",
     mvg->mv[0] >> mv_res, mvg->mv[1] >> mv_res));
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
     "Predicted MV (%i, %i) rate: %i - %i = %i", mv->vx, mv->vy,
     pred_mv_rates[pi], mv->mv_rate, pred_mv_rates[pi] - mv->mv_rate));
    dr += pred_mv_rates[pi] - mv->mv_rate;
  }
  return dr;
}

#if defined(OD_DUMP_IMAGES) && defined(OD_ANIMATE)
static const unsigned char OD_YCbCr_EDGE[3] = { 41, 240, 110 };

static void od_mv_dp_animate_encode(od_enc_ctx *enc,
 od_mv_dp_node *dp, int has_gap) {
  od_state *state;
  od_mv_dp_node *dp0;
  char iter_label[16];
  int active_states[OD_DP_NSTATES_MAX];
  int prev_active_states[OD_DP_NSTATES_MAX];
  int nactive_states;
  int nprev_active_states;
  int dp_state;
  int si;
  int d0vx;
  int d0vy;
  int d1vx;
  int d1vy;
  int x0;
  int y0;
  state = &enc->state;
  od_state_mc_predict(state, state->ref_imgs
   + state->ref_imgi[OD_FRAME_SELF]);
  od_encode_fill_vis(enc);
  /*Now, draw the current state of the DP.*/
  dp0 = dp;
  /*First draw the candidate edge labels for the active trellis paths.*/
  for (si = 0; si < dp0->nstates; si++) prev_active_states[si] = si;
  nprev_active_states = dp0->nstates;
  nactive_states = 0;
  do {
    if (nactive_states > 0) {
      d0vx = dp[0].mv->vx;
      d0vy = dp[0].mv->vy;
      x0 = (d0vx << (OD_LOG_MVBSIZE_MIN + 1)) + (OD_BUFFER_PADDING << 1);
      y0 = (d0vy << (OD_LOG_MVBSIZE_MIN + 1)) + (OD_BUFFER_PADDING << 1);
      if (nprev_active_states > 0) {
        int mvb_sz;
        int x1;
        int y1;
        d1vx = dp[1].mv->vx;
        d1vy = dp[1].mv->vy;
        x1 = (d1vx << (OD_LOG_MVBSIZE_MIN + 1)) + (OD_BUFFER_PADDING << 1);
        y1 = (d1vy << (OD_LOG_MVBSIZE_MIN + 1)) + (OD_BUFFER_PADDING << 1);
        od_img_draw_line(&enc->vis_img, x0, y0, x1, y1, OD_YCbCr_EDGE);
        if (d1vx - d0vx > 1) {
          mvb_sz = d1vx - d0vx;
          if (!has_gap || dp + 1 != dp0) mvb_sz >>= 1;
          if (!state->mv_grid[d0vy][d0vx + mvb_sz].valid) {
            if (d0vy >= mvb_sz
             && state->mv_grid[d0vy - mvb_sz][d0vx + mvb_sz].valid) {
              od_img_draw_line(&enc->vis_img,
               x0 + (mvb_sz << (OD_LOG_MVBSIZE_MIN + 1)),
               y0 - (mvb_sz << (OD_LOG_MVBSIZE_MIN + 1)),
               x0 + (mvb_sz << (OD_LOG_MVBSIZE_MIN + 1)), y1, OD_YCbCr_EDGE);
            }
            if (dp[0].mv->vy <= state->nvmvbs - mvb_sz
             && state->mv_grid[d0vy + mvb_sz][d0vx + mvb_sz].valid) {
              od_img_draw_line(&enc->vis_img,
               x0 + (mvb_sz << (OD_LOG_MVBSIZE_MIN + 1)),
               y0 + (mvb_sz << (OD_LOG_MVBSIZE_MIN + 1)),
               x0 + (mvb_sz << (OD_LOG_MVBSIZE_MIN + 1)), y1, OD_YCbCr_EDGE);
            }
          }
        }
        else if (d1vy - d0vy > 1) {
          mvb_sz = d1vy - d0vy;
          if (!has_gap || dp + 1 != dp0) mvb_sz >>= 1;
          if (!state->mv_grid[d0vy + mvb_sz][d0vx].valid) {
            if (d0vx >= mvb_sz
             && state->mv_grid[d0vy + mvb_sz][d0vx - mvb_sz].valid) {
              od_img_draw_line(&enc->vis_img,
               x0 - (mvb_sz << (OD_LOG_MVBSIZE_MIN + 1)),
               y0 + (mvb_sz << (OD_LOG_MVBSIZE_MIN + 1)),
               x1, y0 + (mvb_sz << (OD_LOG_MVBSIZE_MIN + 1)), OD_YCbCr_EDGE);
            }
            if (d0vx <= state->nhmvbs - mvb_sz
             && state->mv_grid[d0vy + mvb_sz][d0vx + mvb_sz].valid) {
              od_img_draw_line(&enc->vis_img,
               x0 + (mvb_sz << (OD_LOG_MVBSIZE_MIN + 1)),
               y0 + (mvb_sz << (OD_LOG_MVBSIZE_MIN + 1)),
               x1, y0 + (mvb_sz << (OD_LOG_MVBSIZE_MIN + 1)), OD_YCbCr_EDGE);
            }
          }
        }
      }
    }
    OD_COPY(active_states, prev_active_states, nprev_active_states);
    nactive_states = nprev_active_states;
    /*Follow the chain backwards to find the new active states.*/
    nprev_active_states = 0;
    for (si = 0; si < nactive_states; si++) {
      int sj;
      dp_state = dp[0].states[active_states[si]].prevsi;
      for (sj = 0;
       sj < nprev_active_states && prev_active_states[sj] != dp_state; sj++);
      if (sj >= nprev_active_states) {
        prev_active_states[nprev_active_states++] = dp_state;
      }
    }
  }
  while ((dp--)->states[0].prevsi >= 0);
  /*Now, draw all the candidate MVs in the active trellis paths.
    These two steps used to be together; now they're apart.
    Sorry for the mess that caused.*/
  /*Redraw the MVs, so they appear over the edge labels above.*/
  od_state_draw_mvs(enc);
  for (si = 0; si < dp0->nstates; si++) prev_active_states[si] = si;
  nprev_active_states = dp0->nstates;
  nactive_states = 0;
  dp = dp0;
  do {
    if (nactive_states > 0) {
      d0vx = dp[0].mv->vx;
      d0vy = dp[0].mv->vy;
      x0 = (d0vx << (OD_LOG_MVBSIZE_MIN + 1)) + (OD_BUFFER_PADDING << 1);
      y0 = (d0vy << (OD_LOG_MVBSIZE_MIN + 1)) + (OD_BUFFER_PADDING << 1);
      if (!has_gap || dp + 1 != dp0) {
        d1vx = dp[1].mv->vx;
        d1vy = dp[1].mv->vy;
        x0 = (d1vx << (OD_LOG_MVBSIZE_MIN + 1)) + (OD_BUFFER_PADDING << 1);
        y0 = (d1vy << (OD_LOG_MVBSIZE_MIN + 1)) + (OD_BUFFER_PADDING << 1);
        for (si = 0; si < nactive_states; si++) {
          dp_state = active_states[si];
          od_img_draw_line(&enc->vis_img, x0, y0,
           x0 + OD_DIV_ROUND_POW2(dp[1].states[dp_state].mv[0], 2, 2),
           y0 + OD_DIV_ROUND_POW2(dp[1].states[dp_state].mv[1], 2, 2),
           OD_YCbCr_MVCAND);
        }
      }
    }
    OD_COPY(active_states, prev_active_states, nprev_active_states);
    nactive_states = nprev_active_states;
    /*Follow the chain backwards to find the new active states.*/
    nprev_active_states = 0;
    for (si = 0; si < nactive_states; si++) {
      int sj;
      dp_state = dp[0].states[active_states[si]].prevsi;
      for (sj = 0;
       sj < nprev_active_states && prev_active_states[sj] != dp_state; sj++);
      if (sj >= nprev_active_states) {
        prev_active_states[nprev_active_states++] = dp_state;
      }
    }
  }
  while ((dp--)->states[0].prevsi >= 0);
  /*Draw the first state's MV's.*/
  d1vx = dp[1].mv->vx;
  d1vy = dp[1].mv->vy;
  x0 = (d1vx << (OD_LOG_MVBSIZE_MIN + 1)) + (OD_BUFFER_PADDING << 1);
  y0 = (d1vy << (OD_LOG_MVBSIZE_MIN + 1)) + (OD_BUFFER_PADDING << 1);
  for (si = 0; si < nactive_states; si++) {
    dp_state = active_states[si];
    od_img_draw_line(&enc->vis_img, x0, y0,
     x0 + OD_DIV_ROUND_POW2(dp[1].states[dp_state].mv[0], 2, 2),
     y0 + OD_DIV_ROUND_POW2(dp[1].states[dp_state].mv[1], 2, 2),
     OD_YCbCr_MVCAND);
  }
  sprintf(iter_label, "ani%08i", enc->ani_iter++);
  od_state_dump_img(state, &enc->vis_img, iter_label);
}
#endif

/*Row refinement.*/

static void od_mv_dp_row_init(od_mv_est_ctx *est,
 od_mv_dp_node *dp, int vx, int vy, od_mv_dp_node *prev_dp) {
  od_state *state;
  int nhmvbs;
  int nvmvbs;
  int level;
  int pred_hist;
  int npred;
  int nchangeable;
  int pi;
  state = &est->enc->state;
  dp->mv = est->mvs[vy] + vx;
  dp->mvg = state->mv_grid[vy] + vx;
  dp->original_mv[0] = dp->mvg->mv[0];
  dp->original_mv[1] = dp->mvg->mv[1];
  dp->original_mv_rate = dp->mv->mv_rate;
  nhmvbs = state->nhmvbs;
  nvmvbs = state->nvmvbs;
  /*Get the list of MVs we help predict.*/
  level = OD_MC_LEVEL[vy & OD_MVB_MASK][vx & OD_MVB_MASK];
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
   "Initializing node (%i, %i) [%i, %i] at level %i:",
   vx, vy, vx << OD_LOG_MVBSIZE_MIN, vy << OD_LOG_MVBSIZE_MIN, level));
  npred = nchangeable = 0;
  for (pi = 0; pi < OD_NPREDICTED[level]; pi++) {
    int px;
    int py;
    px = vx + OD_ROW_PREDICTED[level][pi][0];
    if (px < 0 || px > nhmvbs) continue;
    py = vy + OD_ROW_PREDICTED[level][pi][1];
    if (py < 0 || py > nvmvbs) continue;
    if (state->mv_grid[py][px].valid) {
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
       "Adding (%i, %i) as a PREDICTED MV.", px, py));
      dp->predicted_mvgs[npred] = state->mv_grid[py] + px;
      dp->predicted_mvs[npred] = est->mvs[py] + px;
      if (pi < OD_NROW_PRED_CHANGEABLE[level]) {
        OD_LOG((OD_LOG_MOTION_ESTIMATION,
         OD_LOG_DEBUG, "It is CHANGEABLE."));
        dp->original_mv_rates[npred] = est->mvs[py][px].mv_rate;
        nchangeable++;
      }
      npred++;
    }
  }
  dp->npredicted = npred;
  dp->npred_changeable = nchangeable;
  /*Now, figure out the earliest DP node that influences our own prediction,
     or that of one of the other MVs we predict.*/
  pred_hist = OD_ROW_PRED_HIST_SIZE[level];
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
   "Marking history up to %i back: %i>=%i",
   pred_hist, prev_dp != NULL ? prev_dp->mv->vx : -1, vx - pred_hist));
  if (prev_dp != NULL && prev_dp->mv->vx >= vx - pred_hist) {
    od_mv_dp_node *dp_pred;
    dp_pred = prev_dp;
    while (dp_pred->mv->vx > vx - pred_hist
     && dp_pred->states[0].prevsi >= 0) {
      dp_pred--;
    }
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
     "Stopped at (%i, %i) (%i <= %i? %i) (%i < 0? %i)",
     dp_pred->mv->vx, dp_pred->mv->vy, dp_pred->mv->vx, vx - pred_hist,
     dp_pred->mv->vx <= vx - pred_hist,
     dp_pred->states[0].prevsi, dp_pred->states[0].prevsi < 0));
    if (dp_pred->mv->vx < vx - pred_hist) {
      dp_pred++;
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
       "Too far, incrementing to (%i, %i).",
       dp_pred->mv->vx, dp_pred->mv->vy));
    }
    dp->min_predictor_node = dp_pred;
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
     "State will be restored back to (%i, %i).",
     dp_pred->mv->vx, dp_pred->mv->vy));
  }
  else dp->min_predictor_node = NULL;
}

/*Each DP node holds a list of blocks between it and the previous node that
   either MV affects.
  There can be up to 8 of them, e.g. the following grid,

  A---B---C                                      A--BB--C
  |\  |   |                                      |  ||  |
  | \ |   |                                      |  ||  |
  |  \|   |                                      E--DD--F
  |   D   |                                      A--DD--C
  |    \  |                                      |  ||  |
  |     \ |                                      |  ||  |
  |      \|                                      E--FE--F
  E-------F produces blocks that use these MVs:  E--FE--F
  |      /|                                      |  ||  |
  |     / |                                      |  ||  |
  |    /  |                                      H--GG--J
  |   G   |                                      E--GG--F
  |  /|   |                                      |  ||  |
  | / |   |                                      |  ||  |
  |/  |   |                                      H--II--J
  H---I---J

  If E and F are the MVs on the DP path, then at least one of the two affects
   all 8 blocks.*/

/*Finds which blocks to the left of the first node in a DP row are affected by
   that MV.*/
static void od_mv_dp_first_row_block_setup(od_mv_est_ctx *est,
 od_mv_dp_node *dp, int vx, int vy) {
  int nblocks;
  nblocks = 0;
  if (vx > 0) {
    od_state *state;
    int nvmvbs;
    int level;
    int log_mvb_sz;
    int mvb_sz;
    state = &est->enc->state;
    nvmvbs = state->nvmvbs;
    level = OD_MC_LEVEL[vy & OD_MVB_MASK][vx & OD_MVB_MASK];
    log_mvb_sz = (OD_MC_LEVEL_MAX - level) >> 1;
    mvb_sz = 1 << log_mvb_sz;
    /*DP chains begin either at the left of the frame or after a gap.
      vx > 0, so this is the first vertex after a gap.
      It's impossible to start a DP chain on an even level vertex
       after a gap as every even vertex has an odd vertex as parent
       left of it.
      Therefore this vertex is an odd level.
      This odd level vertex has no left even child, else it would not be
       the first vertex after a gap.
      There can be no split without that missing even child acting as
       parent for the odd vertex the next level down.
      Thus any block UL or LL of this vertex is unsplit.
      We need only check that the blocks exist (not off top or bottom
       of the frame.)*/
    if (vy >= mvb_sz) {
      OD_ASSERT(mvb_sz == 1
       || !state->mv_grid[vy - (mvb_sz >> 1)][vx - (mvb_sz >> 1)].valid);
      dp->blocks[nblocks++] = est->mvs[vy - mvb_sz] + vx - mvb_sz;
      OD_ASSERT(1 << dp->blocks[nblocks - 1]->log_mvb_sz == mvb_sz);
    }
    if (vy <= nvmvbs - mvb_sz) {
      OD_ASSERT(mvb_sz == 1
       || !state->mv_grid[vy + (mvb_sz >> 1)][vx - (mvb_sz >> 1)].valid);
      dp->blocks[nblocks++] = est->mvs[vy] + vx - mvb_sz;
      OD_ASSERT(1 << dp->blocks[nblocks - 1]->log_mvb_sz == mvb_sz);
    }
  }
  dp->nblocks = nblocks;
}

/*Finds which blocks to the left of each subsequent node in a DP row are
   affected by that MV or the previous one.*/
static void od_mv_dp_prev_row_block_setup(od_mv_est_ctx *est,
 od_mv_dp_node *dp, int vx, int vy) {
  od_state *state;
  int nvmvbs;
  int level;
  int prev_level;
  int log_mvb_sz;
  int prev_log_mvb_sz;
  int mvb_sz;
  int nblocks;
  state = &est->enc->state;
  nvmvbs = state->nvmvbs;
  level = OD_MC_LEVEL[vy & OD_MVB_MASK][vx & OD_MVB_MASK];
  log_mvb_sz = (OD_MC_LEVEL_MAX - level) >> 1;
  mvb_sz = 1 << log_mvb_sz;
  prev_level = OD_MC_LEVEL[vy & OD_MVB_MASK][(vx - mvb_sz) & OD_MVB_MASK];
  prev_log_mvb_sz = (OD_MC_LEVEL_MAX - prev_level) >> 1;
  nblocks = 0;
  if (level >= OD_MC_LEVEL_MAX - 1) {
    /*If the current MV is at the bottom of the hierarchy, then we always want
       one block up-left and down-left from it.*/
    if (vy > 0) {
      dp->blocks[nblocks++] = est->mvs[vy - 1] + vx - 1;
      OD_ASSERT(dp->blocks[nblocks - 1]->log_mvb_sz == 0);
      /*But if the _previous_ node was not, then
         a) level must be OD_MC_LEVEL_MAX, and
         b) If the MV up-left from us is not present, then the previous MV also
          affects one more block above us.*/
      if (prev_log_mvb_sz > 0 && !state->mv_grid[vy - 1][vx - 1].valid) {
        OD_ASSERT(level == OD_MC_LEVEL_MAX);
        dp->blocks[nblocks++] = est->mvs[vy - 2] + vx - 1;
        OD_ASSERT(dp->blocks[nblocks - 1]->log_mvb_sz == 0);
      }
    }
    if (vy < nvmvbs) {
      dp->blocks[nblocks++] = est->mvs[vy] + vx - 1;
      OD_ASSERT(dp->blocks[nblocks - 1]->log_mvb_sz == 0);
      /*The same case for a missing down-left MV.*/
      if (prev_log_mvb_sz > 0 && !state->mv_grid[vy + 1][vx - 1].valid) {
        OD_ASSERT(level == OD_MC_LEVEL_MAX);
        dp->blocks[nblocks++] = est->mvs[vy + 1] + vx - 1;
        OD_ASSERT(dp->blocks[nblocks - 1]->log_mvb_sz == 0);
      }
    }
  }
  else {
    int half_mvb_sz;
    int mvb_off;
    half_mvb_sz = mvb_sz >> 1;
    if (vy >= mvb_sz) {
      /*Figure out if the up-left block has been split at all.*/
      if (state->mv_grid[vy - half_mvb_sz][vx - half_mvb_sz].valid) {
        /*It has, now figure out how far down.*/
        mvb_off = half_mvb_sz;
        while (mvb_off > 1
         && state->mv_grid[vy - (mvb_off >> 1)][vx - (mvb_off >> 1)].valid) {
          mvb_off >>= 1;
        }
        dp->blocks[nblocks++] = est->mvs[vy - mvb_off] + vx - mvb_off;
        OD_ASSERT(1 << dp->blocks[nblocks - 1]->log_mvb_sz == mvb_off);
        /*If we were only partially split at this level, then this MV might
           be used for an extra block above.*/
        if (!state->mv_grid[vy - mvb_off][vx].valid) {
          dp->blocks[nblocks++] = est->mvs[vy - (mvb_off << 1)] + vx - mvb_off;
          OD_ASSERT(1 << dp->blocks[nblocks - 1]->log_mvb_sz == mvb_off);
        }
        /*Or for an extra block to the left.*/
        if (!state->mv_grid[vy][vx - mvb_off].valid) {
          dp->blocks[nblocks++] = est->mvs[vy - mvb_off] + vx - (mvb_off << 1);
          OD_ASSERT(1 << dp->blocks[nblocks - 1]->log_mvb_sz == mvb_off);
          /*And if so, if the _previous_ MV might have affected an extra block
             above it.*/
          if (!state->mv_grid[vy - mvb_off][vx - (mvb_off << 1)].valid) {
            dp->blocks[nblocks++] =
             est->mvs[vy - (mvb_off << 1)] + vx - (mvb_off << 1);
            OD_ASSERT(1 << dp->blocks[nblocks - 1]->log_mvb_sz == mvb_off);
          }
        }
      }
      else {
        /*If it has not been split, then this MV affects exactly one block on
           this side.*/
        dp->blocks[nblocks++] = est->mvs[vy - mvb_sz] + vx - mvb_sz;
        OD_ASSERT(1 << dp->blocks[nblocks - 1]->log_mvb_sz == mvb_sz);
        /*But if the previous node was at a higher level of the hierarchy, then
           a) level must be even, and
           b) If the MV up-left from us is not present, then the previous MV
            also affects one more block above us.*/
        if (prev_log_mvb_sz > log_mvb_sz
         && !state->mv_grid[vy - mvb_sz][vx - mvb_sz].valid) {
          OD_ASSERT(!(level & 1));
          dp->blocks[nblocks++] = est->mvs[vy - (mvb_sz << 1)] + vx - mvb_sz;
          OD_ASSERT(1 << dp->blocks[nblocks - 1]->log_mvb_sz == mvb_sz);
        }
      }
    }
    if (vy <= nvmvbs - mvb_sz) {
      /*Figure out if the down-left block has been split at all.*/
      if (state->mv_grid[vy + half_mvb_sz][vx - half_mvb_sz].valid) {
        /*It has, now figure out how far down.*/
        mvb_off = half_mvb_sz;
        while (mvb_off > 1
         && state->mv_grid[vy + (mvb_off >> 1)][vx - (mvb_off >> 1)].valid) {
          mvb_off >>= 1;
        }
        dp->blocks[nblocks++] = est->mvs[vy] + vx - mvb_off;
        OD_ASSERT(1 << dp->blocks[nblocks - 1]->log_mvb_sz == mvb_off);
        /*If we were only partially split at this level, then this MV might
           be used for an extra block below.*/
        if (!state->mv_grid[vy + mvb_off][vx].valid) {
          dp->blocks[nblocks++] = est->mvs[vy + mvb_off] + vx - mvb_off;
          OD_ASSERT(1 << dp->blocks[nblocks - 1]->log_mvb_sz == mvb_off);
        }
        /*Or for an extra block to the left.*/
        if (!state->mv_grid[vy][vx - mvb_off].valid) {
          dp->blocks[nblocks++] = est->mvs[vy] + vx - (mvb_off << 1);
          OD_ASSERT(1 << dp->blocks[nblocks - 1]->log_mvb_sz == mvb_off);
          /*And if so, if the _previous_ MV might have affected an extra block
             below it.*/
          if (!state->mv_grid[vy + mvb_off][vx - (mvb_off << 1)].valid) {
            dp->blocks[nblocks++] =
             est->mvs[vy + mvb_off] + vx - (mvb_off << 1);
            OD_ASSERT(1 << dp->blocks[nblocks - 1]->log_mvb_sz == mvb_off);
          }
        }
      }
      else {
        /*If it has not been split, then this MV affects exactly one block on
           this side.*/
        dp->blocks[nblocks++] = est->mvs[vy] + vx - mvb_sz;
        OD_ASSERT(1 << dp->blocks[nblocks - 1]->log_mvb_sz == mvb_sz);
        /*But if the previous node was at a higher level of the hierarchy, then
           a) level must be even, and
           b) If the MV down-left from us is not present, then the previous MV
            also affects one more block below us.*/
        if (prev_log_mvb_sz > log_mvb_sz
         && !state->mv_grid[vy + mvb_sz][vx - mvb_sz].valid) {
          OD_ASSERT(!(level & 1));
          dp->blocks[nblocks++] = est->mvs[vy + mvb_sz] + vx - mvb_sz;
          OD_ASSERT(1 << dp->blocks[nblocks - 1]->log_mvb_sz == mvb_sz);
        }
      }
    }
  }
  dp->nblocks = nblocks;
}

/*Finds which blocks to the right of the last node in a DP row are affected by
   that MV.
  The caller special-cases the case where there are no following blocks, so we
   don't need to check if we're on the edge of the frame (it has already done
   so).*/
static void od_mv_dp_last_row_block_setup(od_mv_est_ctx *est,
 od_mv_dp_node *dp, int vx, int vy) {
  od_state *state;
  int nvmvbs;
  int level;
  int log_mvb_sz;
  int mvb_sz;
  int nblocks;
  state = &est->enc->state;
  nvmvbs = state->nvmvbs;
  level = OD_MC_LEVEL[vy & OD_MVB_MASK][vx & OD_MVB_MASK];
  log_mvb_sz = (OD_MC_LEVEL_MAX - level) >> 1;
  mvb_sz = 1 << log_mvb_sz;
  nblocks = 0;
  /*The logic here is the same as in first_row, but at the end of a
     chain.
    The conclusion is the same; this is an odd vertex and the blocks
     to the right cannot be split.*/
  if (vy >= mvb_sz) {
    OD_ASSERT(mvb_sz == 1
     || !state->mv_grid[vy - (mvb_sz >> 1)][vx + (mvb_sz >> 1)].valid);
    dp->blocks[nblocks++] = est->mvs[vy - mvb_sz] + vx;
    OD_ASSERT(1 << dp->blocks[nblocks - 1]->log_mvb_sz == mvb_sz);
  }

  if (vy <= nvmvbs - mvb_sz) {
    OD_ASSERT(mvb_sz == 1
     || !state->mv_grid[vy + (mvb_sz >> 1)][vx + (mvb_sz >> 1)].valid);
    dp->blocks[nblocks++] = est->mvs[vy] + vx;
    OD_ASSERT(1 << dp->blocks[nblocks - 1]->log_mvb_sz == mvb_sz);
  }
  dp->nblocks = nblocks;
}

static void od_mv_dp_restore_row_state(od_mv_dp_node *dp) {
  od_mv_grid_pt *mvg;
  int pi;
  do {
    /*Restore the state for this MV itself.*/
    dp->mv->mv_rate = dp->original_mv_rate;
    mvg = dp->mvg;
    mvg->mv[0] = dp->original_mv[0];
    mvg->mv[1] = dp->original_mv[1];
    for (pi = 0; pi < dp->npred_changeable; pi++) {
      /*Restore the state for the MVs this one predicted.*/
      dp->predicted_mvs[pi]->mv_rate = dp->original_mv_rates[pi];
    }
  }
  while ((dp--)->states[0].prevsi >= 0);
}

static void od_mv_dp_install_row_state(od_mv_dp_node *dp, int prevsi) {
  od_mv_dp_node *dp0;
  od_mv_grid_pt *mvg;
  int nextsi;
  int si;
  int pi;
  int bi;
  /*We must install the state going FORWARDS, since the pred_mv_rates may have
     changed several times over the course of the trellis.
    Therefore, first we reverse all of the prevsi pointers to make them act
     like nextsi pointers.*/
  nextsi = -1;
  dp0 = dp;
  for (si = prevsi; si >= 0; si = prevsi) {
    dp--;
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
     "Node %i, prevsi: %i nextsi: %i", (int)(dp0 - dp), prevsi, nextsi));
    prevsi = dp->states[si].prevsi;
    dp->states[si].prevsi = nextsi;
    nextsi = si;
  }
  /*Now we traverse forward installing the rest of the state.*/
  for (si = nextsi; dp < dp0; dp++) {
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
     "Installing state %i for (%i, %i):", si, dp->mv->vx, dp->mv->vy));
    /*Install the state for this MV itself.*/
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
     "Installing current mv_rate for (%i, %i): %i",
     dp->mv->vx, dp->mv->vy, dp->states[si].mv_rate));
    dp->mv->mv_rate = dp->states[si].mv_rate;
    mvg = dp->mvg;
    mvg->mv[0] = dp->states[si].mv[0];
    mvg->mv[1] = dp->states[si].mv[1];
    /*Install the new block SADs.*/
    for (bi = 0; bi < dp->nblocks; bi++) {
      dp->blocks[bi]->sad = dp->states[si].block_sads[bi];
    }
    /*Install the state for the MVs this one predicted.*/
    for (pi = 0; pi < dp->npredicted; pi++) {
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
       "Installing predicted mv_rate for (%i, %i): %i",
       dp->predicted_mvs[pi]->vx, dp->predicted_mvs[pi]->vy,
       dp->states[si].pred_mv_rates[pi]));
      dp->predicted_mvs[pi]->mv_rate = dp->states[si].pred_mv_rates[pi];
    }
    si = dp->states[si].prevsi;
  }
}

static int32_t od_mv_est_refine_row(od_mv_est_ctx *est,
 int vy, int log_dsz, int mv_res, const int *pattern_nsites,
 const od_pattern *pattern) {
  od_state *state;
  od_mv_grid_pt *grid;
  od_mv_grid_pt *pmvg;
  od_mv_grid_pt *mvg;
  od_mv_dp_node *dp_node;
  od_mv_dp_state *cstate;
  od_mv_dp_state *pstate;
  int32_t dcost;
  int nhmvbs;
  int level;
  int log_mvb_sz;
  int mvb_sz;
  int nsites;
  int sitei;
  int site;
  int curx;
  int cury;
  int vx;
  int b;
  state = &est->enc->state;
  nhmvbs = state->nhmvbs;
  grid = state->mv_grid[vy];
  dcost = 0;
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
   "Refining row %i (%i)...", vy, vy << OD_LOG_MVBSIZE_MIN));
  for (vx = 0;; vx++) {
    int32_t block_sads[OD_DP_NSTATES_MAX][OD_DP_NBLOCKS_MAX];
    int32_t best_cost;
    int32_t cost;
    int32_t best_dd;
    int32_t dd;
    int cur_mv_rates[OD_DP_NSTATES_MAX];
    int pred_mv_rates[OD_DP_NSTATES_MAX][OD_DP_NPREDICTED_MAX];
    int best_dr;
    int dr;
    int best_si;
    int si;
    for (; vx <= nhmvbs && !grid[vx].valid; vx++);
    if (vx > nhmvbs) break;
    level = OD_MC_LEVEL[vy & OD_MVB_MASK][vx & OD_MVB_MASK];
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
     "Starting DP at vertex %i (%i), level %i",
     vx, vx << OD_LOG_MVBSIZE_MIN, level));
    log_mvb_sz = (OD_MC_LEVEL_MAX - level) >> 1;
    mvb_sz = 1 << log_mvb_sz;
    mvg = grid + vx;
    curx = mvg->mv[0];
    cury = mvg->mv[1];
    dp_node = est->dp_nodes;
    od_mv_dp_row_init(est, dp_node, vx, vy, NULL);
    od_mv_dp_first_row_block_setup(est, dp_node, vx, vy);
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG, "TESTING block SADs:"));
    od_mv_dp_get_sad_change(est, dp_node, block_sads[0]);
    /*Compute the set of states for the first node.*/
    b = od_mv_est_get_boundary_case(state, vx, vy, curx, cury,
     1 << log_dsz, log_mvb_sz + OD_LOG_MVBSIZE_MIN);
    nsites = pattern_nsites[b];
    for (sitei = 0, site = 4;; sitei++) {
      cstate = dp_node[0].states + sitei;
      cstate->mv[0] = curx + (OD_SITE_DX[site] << log_dsz);
      cstate->mv[1] = cury + (OD_SITE_DY[site] << log_dsz);
      cstate->prevsi = -1;
      mvg->mv[0] = cstate->mv[0];
      mvg->mv[1] = cstate->mv[1];
      cstate->dr = od_mv_dp_get_rate_change(est, dp_node,
       &cstate->mv_rate, cstate->pred_mv_rates, -1, mv_res);
      cstate->dd = od_mv_dp_get_sad_change(est, dp_node,
       cstate->block_sads);
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
       "State: %i (%g, %g)  dr: %i  dd: %i  dopt: %i",
       sitei, 0.125*cstate->mv[0], 0.125*cstate->mv[1], cstate->dr, cstate->dd,
       cstate->dr*est->lambda + (cstate->dd << OD_ERROR_SCALE)));
      if (sitei >= nsites) break;
      site = pattern[b][sitei];
    }
    dp_node[0].nstates = nsites + 1;
    pmvg = mvg;
    while (vx < nhmvbs) {
      /*Find the next available MV to advance to.*/
      if ((level & 1) && !grid[vx + mvb_sz].valid) {
        OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
         "Gap found at %i (%i), stopping", vx, vx << OD_LOG_MVBSIZE_MIN));
        break;
      }
      while (mvb_sz > 1 && grid[vx + (mvb_sz >> 1)].valid) mvb_sz >>= 1;
      vx += mvb_sz;
#if defined(OD_DUMP_IMAGES) && defined(OD_ANIMATE)
      if (daala_granule_basetime(state, state->cur_time) == ANI_FRAME) {
        od_mv_dp_restore_row_state(dp_node);
        od_mv_dp_animate_encode(est->enc, dp_node, 0);
      }
#endif
      level = OD_MC_LEVEL[vy & OD_MVB_MASK][vx & OD_MVB_MASK];
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
       "Continuing DP at vertex %i (%i), level %i",
       vx, vx << OD_LOG_MVBSIZE_MIN, level));
      log_mvb_sz = (OD_MC_LEVEL_MAX - level) >> 1;
      mvb_sz = 1 << log_mvb_sz;
      mvg = grid + vx;
      curx = mvg->mv[0];
      cury = mvg->mv[1];
      od_mv_dp_row_init(est, dp_node + 1, vx, vy, dp_node);
      od_mv_dp_prev_row_block_setup(est, dp_node + 1, vx, vy);
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG, "TESTING block SADs:"));
      if (od_logging_active(OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG)) {
        pmvg->mv[0] = dp_node[0].original_mv[0];
        pmvg->mv[1] = dp_node[0].original_mv[1];
        od_mv_dp_get_sad_change(est, dp_node + 1, block_sads[0]);
      }
      /*Compute the set of states for this node.*/
      b = od_mv_est_get_boundary_case(state,
       vx, vy, curx, cury, 1 << log_dsz, log_mvb_sz + OD_LOG_MVBSIZE_MIN);
      nsites = pattern_nsites[b];
      for (sitei = 0, site = 4;; sitei++) {
        cstate = dp_node[1].states + sitei;
        cstate->mv[0] = curx + (OD_SITE_DX[site] << log_dsz);
        cstate->mv[1] = cury + (OD_SITE_DY[site] << log_dsz);
        best_si = 0;
        best_dr = dp_node[0].states[0].dr;
        best_dd = dp_node[0].states[0].dd;
        best_cost = INT_MAX;
        mvg->mv[0] = cstate->mv[0];
        mvg->mv[1] = cstate->mv[1];
        for (si = 0; si < dp_node[0].nstates; si++) {
          pstate = dp_node[0].states + si;
          /*Get the rate change for this state using previous state si.
            This automatically loads the required bits of the trellis path into
             the grid, like the previous MV.*/
          cstate->dr = od_mv_dp_get_rate_change(est, dp_node + 1,
           cur_mv_rates + si, pred_mv_rates[si], si, mv_res);
          /*Test against the previous state.*/
          dr = pstate->dr + cstate->dr;
          dd = pstate->dd + od_mv_dp_get_sad_change(est,
           dp_node + 1, block_sads[si]);
          cost = dr*est->lambda + (dd << OD_ERROR_SCALE);
          OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
           "State: %2i (%7g, %7g) P.State: %i  dr: %3i  dd: %6i  dopt: %7i",
           sitei, 0.125*cstate->mv[0], 0.125*cstate->mv[1],
           si, dr, dd, cost));
          if (cost < best_cost) {
            best_si = si;
            best_cost = cost;
            best_dd = dd;
            best_dr = dr;
          }
        }
        OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
         "State: %2i  Best P.State: %i", sitei, best_si));
        cstate->prevsi = best_si;
        cstate->dr = best_dr;
        cstate->dd = best_dd;
        OD_COPY(cstate->block_sads, block_sads[best_si], dp_node[1].nblocks);
        cstate->mv_rate = cur_mv_rates[best_si];
        OD_COPY(cstate->pred_mv_rates, pred_mv_rates[best_si],
         dp_node[1].npredicted);
        if (sitei >= nsites) break;
        site = pattern[b][sitei];
      }
      dp_node[1].nstates = nsites + 1;
      dp_node++;
      pmvg = mvg;
    }
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
     "Finished DP at vertex %i (%i)",
     dp_node[0].mv->vx, dp_node[0].mv->vx << OD_LOG_MVBSIZE_MIN));
    best_si = 0;
    best_cost = INT_MAX;
    /*TODO: Once we stop optimizing at arbitrary places, we'll need to
       compute the rate change of MVs we didn't get to.*/
    dp_node[1].npredicted = dp_node[1].npred_changeable = 0;
    if (dp_node[0].mv->vx < nhmvbs) {
      od_mv_dp_last_row_block_setup(est, dp_node + 1, dp_node[0].mv->vx, vy);
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG, "TESTING block SADs:"));
      if (od_logging_active(OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG)) {
        pmvg->mv[0] = dp_node[0].original_mv[0];
        pmvg->mv[1] = dp_node[0].original_mv[1];
        od_mv_dp_get_sad_change(est, dp_node + 1, block_sads[0]);
      }
      for (si = 0; si < dp_node[0].nstates; si++) {
        pstate = dp_node[0].states + si;
        pmvg->mv[0] = pstate->mv[0];
        pmvg->mv[1] = pstate->mv[1];
        /*Test against the state with a following edge.*/
        dr = pstate->dr;
        dd = pstate->dd + od_mv_dp_get_sad_change(est,
         dp_node + 1, block_sads[si]);
        cost = dr*est->lambda + (dd << OD_ERROR_SCALE);
        OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
         "State: --  P.State: %i  dr: %3i  dd: %6i  dopt: %7i",
         si, dr, dd, cost));
        if (cost < best_cost) {
          best_si = si;
          best_cost = cost;
        }
      }
    }
    /*There are no blocks to accumulate SAD after this one, so pick the best
       state so far.*/
    else {
      dp_node[1].nblocks = 0;
      for (si = 0; si < dp_node[0].nstates; si++) {
        pstate = dp_node[0].states + si;
        cost = pstate->dr*est->lambda + (pstate->dd << OD_ERROR_SCALE);
        if (cost < best_cost) {
          best_si = si;
          best_cost = cost;
        }
      }
    }
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
     "Best P.State: %i dopt: %7i", best_si, best_cost));
    if (best_cost > 0) {
      /*Our optimal path is worse than what we started with!
        Restore the original state and give up.*/
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
       "Best cost (%7i) > 0! Optimization failed.", best_cost));
      od_mv_dp_restore_row_state(dp_node);
    }
    else {
      int bi;
#if defined(OD_DUMP_IMAGES) && defined(OD_ANIMATE)
      if (daala_granule_basetime(state, state->cur_time) == ANI_FRAME) {
        char iter_label[16];
        od_mv_dp_restore_row_state(dp_node);
        od_mv_dp_animate_encode(est->enc, dp_node, 0);
        od_mv_dp_install_row_state(dp_node + 1, best_si);
        od_state_mc_predict(state, state->ref_imgs
         + state->ref_imgi[OD_FRAME_SELF]);
        od_encode_fill_vis(est->enc);
        sprintf(iter_label, "ani%08i", est->enc->ani_iter++);
        od_state_dump_img(state, &est->enc->vis_img, iter_label);
      }
#endif
      /*Update the state along the optimal path.*/
      od_mv_dp_install_row_state(dp_node + 1, best_si);
      /*Store the SADs from this last node, too.*/
      for (bi = 0; bi < dp_node[1].nblocks; bi++) {
        dp_node[1].blocks[bi]->sad = block_sads[best_si][bi];
      }
      dcost += best_cost;
    }
  }
#if defined(OD_ENABLE_ASSERTIONS) || defined(OD_LOGGING_ENABLED)
  od_mv_est_check_rd_state(est, mv_res);
#endif
  return dcost;
}

/*Column refinement.*/

static void od_mv_dp_col_init(od_mv_est_ctx *est,
 od_mv_dp_node *dp, int vx, int vy, od_mv_dp_node *prev_dp) {
  od_state *state;
  int nhmvbs;
  int nvmvbs;
  int level;
  int pred_hist;
  int npred;
  int nchangeable;
  int pi;
  state = &est->enc->state;
  dp->mv = est->mvs[vy] + vx;
  dp->mvg = state->mv_grid[vy] + vx;
  dp->original_mv[0] = dp->mvg->mv[0];
  dp->original_mv[1] = dp->mvg->mv[1];
  dp->original_mv_rate = dp->mv->mv_rate;
  nhmvbs = state->nhmvbs;
  nvmvbs = state->nvmvbs;
  /*Get the list of MVs we help predict.*/
  level = OD_MC_LEVEL[vy & OD_MVB_MASK][vx & OD_MVB_MASK];
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
   "Initializing node (%i, %i) [%i, %i] at level %i:",
   vx, vy, vx << OD_LOG_MVBSIZE_MIN, vy << OD_LOG_MVBSIZE_MIN, level));
  npred = nchangeable = 0;
  for (pi = 0; pi < OD_NPREDICTED[level]; pi++) {
    int px;
    int py;
    px = vx + OD_COL_PREDICTED[level][pi][0];
    if (px < 0 || px > nhmvbs) continue;
    py = vy + OD_COL_PREDICTED[level][pi][1];
    if (py < 0 || py > nvmvbs) continue;
    if (state->mv_grid[py][px].valid) {
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
       "Adding (%i, %i) as a PREDICTED MV.", px, py));
      dp->predicted_mvgs[npred] = state->mv_grid[py] + px;
      dp->predicted_mvs[npred] = est->mvs[py] + px;
      if (pi < OD_NCOL_PRED_CHANGEABLE[level]) {
        OD_LOG((OD_LOG_MOTION_ESTIMATION,
         OD_LOG_DEBUG, "It is CHANGEABLE."));
        dp->original_mv_rates[npred] = est->mvs[py][px].mv_rate;
        nchangeable++;
      }
      npred++;
    }
  }
  dp->npredicted = npred;
  dp->npred_changeable = nchangeable;
  /*Now, figure out the earliest DP node that influences our own prediction,
     or that of one of the other MVs we predict.*/
  pred_hist = OD_COL_PRED_HIST_SIZE[level];
  if (prev_dp != NULL && prev_dp->mv->vy >= vy - pred_hist) {
    od_mv_dp_node *dp_pred;
    for (dp_pred = prev_dp; dp_pred->mv->vy > vy - pred_hist
     && dp_pred->states[0].prevsi >= 0; dp_pred--);
    if (dp_pred->mv->vy < vy - pred_hist) dp_pred++;
    dp->min_predictor_node = dp_pred;
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
     "State will be restored back to (%i, %i).",
     dp_pred->mv->vx, dp_pred->mv->vy));
  }
  else dp->min_predictor_node = NULL;
}

/*Finds which blocks above the first node in a DP column are affected by that
   MV.*/
static void od_mv_dp_first_col_block_setup(od_mv_est_ctx *est,
 od_mv_dp_node *dp, int vx, int vy) {
  int nblocks;
  nblocks = 0;
  if (vy > 0) {
    od_state *state;
    int nhmvbs;
    int level;
    int log_mvb_sz;
    int mvb_sz;
    state = &est->enc->state;
    nhmvbs = state->nhmvbs;
    level = OD_MC_LEVEL[vy & OD_MVB_MASK][vx & OD_MVB_MASK];
    log_mvb_sz = (OD_MC_LEVEL_MAX - level) >> 1;
    mvb_sz = 1 << log_mvb_sz;
    /*DP chains begin either at the top of the frame or after a gap.
      vy > 0, so this is the first vertex after a gap.
      It's impossible to start a DP chain on an even level vertex
       after a gap as every even vertex has an odd vertex as parent
       above it.
      Therefore this vertex is an odd level.
      This odd level vertex has no even child above, else it would not
       be the first vertex after a gap.
      There can be no split without that missing even child acting as
       parent for the odd vertex the next level down.
      Thus any block UL or UR of this vertex is unsplit.
      We need only check that the blocks exist (not off left or right
       of the frame.)*/
    if (vx >= mvb_sz) {
      OD_ASSERT(mvb_sz == 1
       || !state->mv_grid[vy - (mvb_sz >> 1)][vx - (mvb_sz >> 1)].valid);
      dp->blocks[nblocks++] = est->mvs[vy - mvb_sz] + vx - mvb_sz;
      OD_ASSERT(1 << dp->blocks[nblocks - 1]->log_mvb_sz == mvb_sz);
    }
    if (vx <= nhmvbs - mvb_sz) {
      OD_ASSERT(mvb_sz == 1
       || !state->mv_grid[vy - (mvb_sz >> 1)][vx + (mvb_sz >> 1)].valid);
      dp->blocks[nblocks++] = est->mvs[vy - mvb_sz] + vx;
      OD_ASSERT(1 << dp->blocks[nblocks - 1]->log_mvb_sz == mvb_sz);
    }
  }
  dp->nblocks = nblocks;
}

/*Finds which blocks above each subsequent node in a DP row are affected by
   that MV or the previous one.*/
static void od_mv_dp_prev_col_block_setup(od_mv_est_ctx *est,
 od_mv_dp_node *dp, int vx, int vy) {
  od_state *state;
  int nhmvbs;
  int level;
  int prev_level;
  int log_mvb_sz;
  int prev_log_mvb_sz;
  int mvb_sz;
  int nblocks;
  state = &est->enc->state;
  nhmvbs = state->nhmvbs;
  level = OD_MC_LEVEL[vy & OD_MVB_MASK][vx & OD_MVB_MASK];
  log_mvb_sz = (OD_MC_LEVEL_MAX - level) >> 1;
  mvb_sz = 1 << log_mvb_sz;
  prev_level = OD_MC_LEVEL[(vy - mvb_sz) & OD_MVB_MASK][vx & OD_MVB_MASK];
  prev_log_mvb_sz = (OD_MC_LEVEL_MAX - prev_level) >> 1;
  nblocks = 0;
  if (level >= OD_MC_LEVEL_MAX - 1) {
    /*If the current MV is at the bottom of the hierarchy, then we always want
       one block up-left and up-right from it.*/
    if (vx > 0) {
      dp->blocks[nblocks++] = est->mvs[vy - 1] + vx - 1;
      OD_ASSERT(dp->blocks[nblocks - 1]->log_mvb_sz == 0);
      /*But if the _previous_ node was not, then
         a) level must be OD_MC_LEVEL_MAX, and
         b) If the MV up-left from us is not present, then the previous MV also
          affects one more block to the left of us.*/
      if (prev_log_mvb_sz > 0 && !state->mv_grid[vy - 1][vx - 1].valid) {
        dp->blocks[nblocks++] = est->mvs[vy - 1] + vx - 2;
        OD_ASSERT(dp->blocks[nblocks - 1]->log_mvb_sz == 0);
      }
    }
    if (vx < nhmvbs) {
      dp->blocks[nblocks++] = est->mvs[vy - 1] + vx;
      OD_ASSERT(dp->blocks[nblocks - 1]->log_mvb_sz == 0);
      /*The same case for a missing up-right MV.*/
      if (prev_log_mvb_sz > 0 && !state->mv_grid[vy - 1][vx + 1].valid) {
        dp->blocks[nblocks++] = est->mvs[vy - 1] + vx + 1;
        OD_ASSERT(dp->blocks[nblocks - 1]->log_mvb_sz == 0);
      }
    }
  }
  else {
    int half_mvb_sz;
    int mvb_off;
    half_mvb_sz = mvb_sz >> 1;
    if (vx >= mvb_sz) {
      /*Figure out if the up-left block has been split at all.*/
      if (state->mv_grid[vy - half_mvb_sz][vx - half_mvb_sz].valid) {
        /*It has, now figure out how far down.*/
        mvb_off = half_mvb_sz;
        while (mvb_off > 1
         && state->mv_grid[vy - (mvb_off >> 1)][vx - (mvb_off >> 1)].valid) {
          mvb_off >>= 1;
        }
        dp->blocks[nblocks++] = est->mvs[vy - mvb_off] + vx - mvb_off;
        OD_ASSERT(1 << dp->blocks[nblocks - 1]->log_mvb_sz == mvb_off);
        /*If we were only partially split at this level, then this MV might
           be used for an extra block to the left.*/
        if (!state->mv_grid[vy][vx - mvb_off].valid) {
          dp->blocks[nblocks++] = est->mvs[vy - mvb_off] + vx - (mvb_off << 1);
          OD_ASSERT(1 << dp->blocks[nblocks - 1]->log_mvb_sz == mvb_off);
        }
        /*Or for an extra block above.*/
        if (!state->mv_grid[vy - mvb_off][vx].valid) {
          dp->blocks[nblocks++] = est->mvs[vy - (mvb_off << 1)] + vx - mvb_off;
          OD_ASSERT(1 << dp->blocks[nblocks - 1]->log_mvb_sz == mvb_off);
          /*And if so, if the _previous_ MV might have affected an extra block
             to the left of it.*/
          if (!state->mv_grid[vy - (mvb_off << 1)][vx - mvb_off].valid) {
            dp->blocks[nblocks++] =
             est->mvs[vy - (mvb_off << 1)] + vx - (mvb_off << 1);
            OD_ASSERT(1 << dp->blocks[nblocks - 1]->log_mvb_sz == mvb_off);
          }
        }
      }
      else {
        /*If it has not been split, then this MV affects exactly one block on
           this side.*/
        dp->blocks[nblocks++] = est->mvs[vy - mvb_sz] + vx - mvb_sz;
        OD_ASSERT(1 << dp->blocks[nblocks - 1]->log_mvb_sz == mvb_sz);
        /*But if the previous node was at a higher level of the hierarchy, then
           a) level must be even, and
           b) If the MV up-left from us is not present, then the previous MV
            also affects one more block to the left of us.*/
        if (prev_log_mvb_sz > log_mvb_sz
         && !state->mv_grid[vy - mvb_sz][vx - mvb_sz].valid) {
          dp->blocks[nblocks++] = est->mvs[vy - mvb_sz] + vx - (mvb_sz << 1);
          OD_ASSERT(1 << dp->blocks[nblocks - 1]->log_mvb_sz == mvb_sz);
        }
      }
    }
    if (vx <= nhmvbs - mvb_sz) {
      /*Figure out if the up-right block has been split at all.*/
      if (state->mv_grid[vy - half_mvb_sz][vx + half_mvb_sz].valid) {
        /*It has, now figure out how far down.*/
        mvb_off = half_mvb_sz;
        while (mvb_off > 1
         && state->mv_grid[vy - (mvb_off >> 1)][vx + (mvb_off >> 1)].valid) {
          mvb_off >>= 1;
        }
        dp->blocks[nblocks++] = est->mvs[vy - mvb_off] + vx;
        OD_ASSERT(1 << dp->blocks[nblocks - 1]->log_mvb_sz == mvb_off);
        /*If we were only partially split at this level, then this MV might
           be used for an extra block to the right.*/
        if (!state->mv_grid[vy][vx + mvb_off].valid) {
          dp->blocks[nblocks++] = est->mvs[vy - mvb_off] + vx + mvb_off;
          OD_ASSERT(1 << dp->blocks[nblocks - 1]->log_mvb_sz == mvb_off);
        }
        /*Or for an extra block above.*/
        if (!state->mv_grid[vy - mvb_off][vx].valid) {
          dp->blocks[nblocks++] = est->mvs[vy - (mvb_off << 1)] + vx;
          OD_ASSERT(1 << dp->blocks[nblocks - 1]->log_mvb_sz == mvb_off);
          /*And if so, if the _previous_ MV might have affected an extra block
             to the right of it.*/
          if (!state->mv_grid[vy - (mvb_off << 1)][vx + mvb_off].valid) {
            dp->blocks[nblocks++] =
             est->mvs[vy - (mvb_off << 1)] + vx + mvb_off;
            OD_ASSERT(1 << dp->blocks[nblocks - 1]->log_mvb_sz == mvb_off);
          }
        }
      }
      else {
        /*If it has not been split, then this MV affects exactly one block on
           this side.*/
        dp->blocks[nblocks++] = est->mvs[vy - mvb_sz] + vx;
        OD_ASSERT(1 << dp->blocks[nblocks - 1]->log_mvb_sz == mvb_sz);
        /*But if the previous node was at a higher level of the hierarchy, then
           a) level must be even, and
           b) If the MV up-right from us is not present, then the previous MV
            also affects one more block to the right of us.*/
        if (prev_log_mvb_sz > log_mvb_sz
         && !state->mv_grid[vy - mvb_sz][vx + mvb_sz].valid) {
          dp->blocks[nblocks++] = est->mvs[vy - mvb_sz] + vx + mvb_sz;
          OD_ASSERT(1 << dp->blocks[nblocks - 1]->log_mvb_sz == mvb_sz);
        }
      }
    }
  }
  dp->nblocks = nblocks;
}

/*Finds which blocks below the last node in a DP column are affected by that
   MV.
  The caller special-cases the case where there are no following blocks, so we
   don't need to check if we're on the edge of the frame (it has already done
   so).*/
static void od_mv_dp_last_col_block_setup(od_mv_est_ctx *est,
 od_mv_dp_node *dp, int vx, int vy) {
  od_state *state;
  int nhmvbs;
  int level;
  int log_mvb_sz;
  int mvb_sz;
  int nblocks;
  state = &est->enc->state;
  nhmvbs = state->nhmvbs;
  level = OD_MC_LEVEL[vy & OD_MVB_MASK][vx & OD_MVB_MASK];
  log_mvb_sz = (OD_MC_LEVEL_MAX - level) >> 1;
  mvb_sz = 1 << log_mvb_sz;
  nblocks = 0;
  /*The logic here is the same as in first_col, but at the end of a
     chain.
    The conclusion is the same; this is an odd vertex and the blocks
     below cannot be split.*/
  if (vx >= mvb_sz) {
    OD_ASSERT(mvb_sz == 1
     || !state->mv_grid[vy + (mvb_sz >> 1)][vx - (mvb_sz >> 1)].valid);
    dp->blocks[nblocks++] = est->mvs[vy] + vx - mvb_sz;
    OD_ASSERT(1 << dp->blocks[nblocks - 1]->log_mvb_sz == mvb_sz);
  }
  if (vx <= nhmvbs - mvb_sz) {
    OD_ASSERT(mvb_sz == 1
     || !state->mv_grid[vy + (mvb_sz >> 1)][vx + (mvb_sz >> 1)].valid);
    dp->blocks[nblocks++] = est->mvs[vy] + vx;
    OD_ASSERT(1 << dp->blocks[nblocks - 1]->log_mvb_sz == mvb_sz);
  }
  dp->nblocks = nblocks;
}

static void od_mv_dp_restore_col_state(od_mv_dp_node *dp) {
  od_mv_grid_pt *mvg;
  int pi;
  do {
    /*Restore the state for this MV itself.*/
    dp->mv->mv_rate = dp->original_mv_rate;
    mvg = dp->mvg;
    mvg->mv[0] = dp->original_mv[0];
    mvg->mv[1] = dp->original_mv[1];
    for (pi = 0; pi < dp->npred_changeable; pi++) {
      /*Restore the state for the MVs this one predicted.*/
      dp->predicted_mvs[pi]->mv_rate = dp->original_mv_rates[pi];
    }
  }
  while ((dp--)->states[0].prevsi >= 0);
}

static void od_mv_dp_install_col_state(od_mv_dp_node *dp, int prevsi) {
  od_mv_dp_node *dp0;
  od_mv_grid_pt *mvg;
  int            nextsi;
  int            si;
  int            pi;
  int            bi;
  /*We must install the state going FORWARDS, since the pred_mv_rates may have
     changed several times over the course of the trellis.
    Therefore, first we reverse all of the prevsi pointers to make them act
     like nextsi pointers.*/
  nextsi = -1;
  dp0 = dp;
  for (si = prevsi; si >= 0; si = prevsi) {
    dp--;
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
     "Node %i, prevsi: %i nextsi: %i", (int)(dp0 - dp), prevsi, nextsi));
    prevsi = dp->states[si].prevsi;
    dp->states[si].prevsi = nextsi;
    nextsi = si;
  }
  /*Now we traverse forward installing the rest of the state.*/
  for (si = nextsi; dp < dp0; dp++) {
    /*Install the state for this MV itself.*/
    dp->mv->mv_rate = dp->states[si].mv_rate;
    mvg = dp->mvg;
    mvg->mv[0] = dp->states[si].mv[0];
    mvg->mv[1] = dp->states[si].mv[1];
    /*Install the new block SADs.*/
    for (bi = 0; bi < dp->nblocks; bi++) {
      dp->blocks[bi]->sad = dp->states[si].block_sads[bi];
    }
    /*Install the state for the MVs this one predicted.*/
    for (pi = 0; pi < dp->npredicted; pi++) {
      dp->predicted_mvs[pi]->mv_rate = dp->states[si].pred_mv_rates[pi];
    }
    si = dp->states[si].prevsi;
  }
}

static int32_t od_mv_est_refine_col(od_mv_est_ctx *est,
 int vx, int log_dsz, int mv_res, const int *pattern_nsites,
 const od_pattern *pattern) {
  od_state *state;
  od_mv_grid_pt **grid;
  od_mv_grid_pt *pmvg;
  od_mv_grid_pt *mvg;
  od_mv_dp_node *dp_node;
  od_mv_dp_state *cstate;
  od_mv_dp_state *pstate;
  int32_t dcost;
  int nvmvbs;
  int level;
  int log_mvb_sz;
  int mvb_sz;
  int nsites;
  int sitei;
  int site;
  int curx;
  int cury;
  int vy;
  int b;
  state = &est->enc->state;
  nvmvbs = state->nvmvbs;
  grid = state->mv_grid;
  dcost = 0;
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
   "Refining column %i (%i)...", vx, vx << OD_LOG_MVBSIZE_MIN));
  for (vy = 0;; vy++) {
    int32_t block_sads[OD_DP_NSTATES_MAX][OD_DP_NBLOCKS_MAX];
    int32_t best_cost;
    int32_t cost;
    int32_t best_dd;
    int32_t dd;
    int cur_mv_rates[OD_DP_NSTATES_MAX];
    int pred_mv_rates[OD_DP_NSTATES_MAX][OD_DP_NPREDICTED_MAX];
    int best_dr;
    int dr;
    int best_si;
    int si;
    for (; vy <= nvmvbs && !grid[vy][vx].valid; vy++);
    if (vy > nvmvbs) break;
    level = OD_MC_LEVEL[vy & OD_MVB_MASK][vx & OD_MVB_MASK];
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
     "Starting DP at vertex %i (%i), level %i",
     vy, vy << OD_LOG_MVBSIZE_MIN, level));
    log_mvb_sz = (OD_MC_LEVEL_MAX - level) >> 1;
    mvb_sz = 1 << log_mvb_sz;
    mvg = grid[vy] + vx;
    curx = mvg->mv[0];
    cury = mvg->mv[1];
    dp_node = est->dp_nodes;
    od_mv_dp_col_init(est, dp_node, vx, vy, NULL);
    od_mv_dp_first_col_block_setup(est, dp_node, vx, vy);
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG, "TESTING block SADs:"));
    if (od_logging_active(OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG)) {
      od_mv_dp_get_sad_change(est, dp_node, block_sads[0]);
    }
    /*Compute the set of states for the first node.*/
    b = od_mv_est_get_boundary_case(state,
     vx, vy, curx, cury, 1 << log_dsz, log_mvb_sz + OD_LOG_MVBSIZE_MIN);
    nsites = pattern_nsites[b];
    for (sitei = 0, site = 4;; sitei++) {
      cstate = dp_node[0].states + sitei;
      cstate->mv[0] = curx + (OD_SITE_DX[site] << log_dsz);
      cstate->mv[1] = cury + (OD_SITE_DY[site] << log_dsz);
      cstate->prevsi = -1;
      mvg->mv[0] = cstate->mv[0];
      mvg->mv[1] = cstate->mv[1];
      cstate->dr = od_mv_dp_get_rate_change(est, dp_node,
       &cstate->mv_rate, cstate->pred_mv_rates, -1, mv_res);
      cstate->dd = od_mv_dp_get_sad_change(est, dp_node,
       cstate->block_sads);
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
       "State: %i  dr: %i  dd: %i  dopt: %i", sitei, cstate->dr, cstate->dd,
       cstate->dr*est->lambda + (cstate->dd << OD_ERROR_SCALE)));
      if (sitei >= nsites) break;
      site = pattern[b][sitei];
    }
    dp_node[0].nstates = nsites + 1;
    pmvg = mvg;
    while (vy < nvmvbs) {
      /*Find the next available MV to advance to.*/
      if ((level & 1) && !grid[vy + mvb_sz][vx].valid) {
        OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
         "Gap found at %i (%i), stopping", vy, vy << OD_LOG_MVBSIZE_MIN));
        break;
      }
      while (mvb_sz > 1 && grid[vy + (mvb_sz >> 1)][vx].valid) mvb_sz >>= 1;
      vy += mvb_sz;
#if defined(OD_DUMP_IMAGES) && defined(OD_ANIMATE)
      if (daala_granule_basetime(state, state->cur_time) == ANI_FRAME) {
        od_mv_dp_restore_col_state(dp_node);
        od_mv_dp_animate_encode(est->enc, dp_node, 0);
      }
#endif
      level = OD_MC_LEVEL[vy & OD_MVB_MASK][vx & OD_MVB_MASK];
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
       "Continuing DP at vertex %i (%i), level %i",
       vy, vy << OD_LOG_MVBSIZE_MIN, level));
      log_mvb_sz = (OD_MC_LEVEL_MAX - level) >> 1;
      mvb_sz = 1 << log_mvb_sz;
      mvg = grid[vy] + vx;
      curx = mvg->mv[0];
      cury = mvg->mv[1];
      od_mv_dp_col_init(est, dp_node + 1, vx, vy, dp_node);
      od_mv_dp_prev_col_block_setup(est, dp_node + 1, vx, vy);
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG, "TESTING block SADs:"));
      if (od_logging_active(OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG)) {
        pmvg->mv[0] = dp_node[0].original_mv[0];
        pmvg->mv[1] = dp_node[0].original_mv[1];
        od_mv_dp_get_sad_change(est, dp_node + 1, block_sads[0]);
      }
      /*Compute the set of states for this node.*/
      b = od_mv_est_get_boundary_case(state,
       vx, vy, curx, cury, 1 << log_dsz, log_mvb_sz + OD_LOG_MVBSIZE_MIN);
      nsites = pattern_nsites[b];
      for (sitei = 0, site = 4;; sitei++) {
        cstate = dp_node[1].states + sitei;
        cstate->mv[0] = curx + (OD_SITE_DX[site] << log_dsz);
        cstate->mv[1] = cury + (OD_SITE_DY[site] << log_dsz);
        best_si = 0;
        best_dr = dp_node[0].states[0].dr;
        best_dd = dp_node[0].states[0].dd;
        best_cost = INT_MAX;
        mvg->mv[0] = cstate->mv[0];
        mvg->mv[1] = cstate->mv[1];
        for (si = 0; si < dp_node[0].nstates; si++) {
          pstate = dp_node[0].states + si;
          /*Get the rate change for this state using previous state si.
            This automatically loads the required bits of the trellis path into
             the grid, like the previous MV.*/
          cstate->dr = od_mv_dp_get_rate_change(est, dp_node + 1,
           cur_mv_rates + si, pred_mv_rates[si], si, mv_res);
          /*Test against the previous state.*/
          dr = pstate->dr + cstate->dr;
          dd = pstate->dd + od_mv_dp_get_sad_change(est,
           dp_node + 1, block_sads[si]);
          cost = dr*est->lambda + (dd << OD_ERROR_SCALE);
          OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
           "State: %2i  P.State: %i  dr: %3i  dd: %6i  dopt: %7i",
           sitei, si, dr, dd, cost));
          if (cost < best_cost) {
            best_si = si;
            best_cost = cost;
            best_dd = dd;
            best_dr = dr;
          }
        }
        OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
         "State: %2i  Best P.State: %i", sitei, best_si));
        cstate->prevsi = best_si;
        cstate->dr = best_dr;
        cstate->dd = best_dd;
        OD_COPY(cstate->block_sads, block_sads[cstate->prevsi],
         dp_node[1].nblocks);
        cstate->mv_rate = cur_mv_rates[best_si];
        OD_COPY(cstate->pred_mv_rates, pred_mv_rates[best_si],
         dp_node[1].npredicted);
        if (sitei >= nsites) break;
        site = pattern[b][sitei];
      }
      dp_node[1].nstates = nsites + 1;
      dp_node++;
      pmvg = mvg;
    }
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
     "Finished DP at vertex %i (%i)",
     dp_node[0].mv->vy, dp_node[0].mv->vy << OD_LOG_MVBSIZE_MIN));
    best_si = 0;
    best_cost = INT_MAX;
    /*TODO: Once we stop optimizing at arbitrary places, we'll need to
       compute the rate change of MVs we didn't get to.*/
    dp_node[1].npredicted = dp_node[1].npred_changeable = 0;
    if (dp_node[0].mv->vy < nvmvbs) {
      od_mv_dp_last_col_block_setup(est, dp_node + 1, vx, dp_node[0].mv->vy);
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG, "TESTING block SADs:"));
      if (od_logging_active(OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG)) {
        pmvg->mv[0] = dp_node[0].original_mv[0];
        pmvg->mv[1] = dp_node[0].original_mv[1];
        od_mv_dp_get_sad_change(est, dp_node + 1, block_sads[0]);
      }
      for (si = 0; si < dp_node[0].nstates; si++) {
        pstate = dp_node[0].states + si;
        pmvg->mv[0] = pstate->mv[0];
        pmvg->mv[1] = pstate->mv[1];
        /*Test against the state with a following edge.*/
        dr = pstate->dr;
        dd = pstate->dd + od_mv_dp_get_sad_change(est,
         dp_node + 1, block_sads[si]);
        cost = dr*est->lambda + (dd << OD_ERROR_SCALE);
        OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
         "State: --  P.State: %i  dr: %3i  dd: %6i  dopt: %7i",
         si, dr, dd, cost));
        if (best_si < 0 || cost < best_cost) {
          best_si = si;
          best_cost = cost;
        }
      }
    }
    /*There are no blocks to accumulate SAD after this one, so pick the best
       state so far.*/
    else {
      dp_node[1].nblocks = 0;
      for (si = 0; si < dp_node[0].nstates; si++) {
        dp_node[1].nblocks = 0;
        pstate = dp_node[0].states + si;
        cost = pstate->dr*est->lambda + (pstate->dd << OD_ERROR_SCALE);
        if (best_si < 0 || cost < best_cost) {
          best_si = si;
          best_cost = cost;
        }
      }
    }
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
     "Best P.State: %i dopt: %7i", best_si, best_cost));
    if (best_cost > 0) {
      /*Our optimal path is worse than what we started with!
        Restore the original state and give up.*/
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
       "Best cost (%7i) > 0! Optimization failed.", best_cost));
      od_mv_dp_restore_col_state(dp_node);
    }
    else {
      int bi;
#if defined(OD_DUMP_IMAGES) && defined(OD_ANIMATE)
      if (daala_granule_basetime(state, state->cur_time) == ANI_FRAME) {
        char iter_label[16];
        od_mv_dp_restore_col_state(dp_node);
        od_mv_dp_animate_encode(est->enc, dp_node, 0);
        od_mv_dp_install_col_state(dp_node + 1, best_si);
        od_state_mc_predict(state, state->ref_imgs
         + state->ref_imgi[OD_FRAME_SELF]);
        od_encode_fill_vis(est->enc);
        sprintf(iter_label, "ani%08i", est->enc->ani_iter++);
        od_state_dump_img(state, &est->enc->vis_img, iter_label);
      }
#endif
      /*Update the state along the optimal path.*/
      od_mv_dp_install_col_state(dp_node + 1, best_si);
      /*Store the SADs from this last node, too.*/
      for (bi = 0; bi < dp_node[1].nblocks; bi++) {
        dp_node[1].blocks[bi]->sad = block_sads[best_si][bi];
      }
      dcost += best_cost;
    }
  }
#if defined(OD_ENABLE_ASSERTIONS) || defined(OD_LOGGING_ENABLED)
  od_mv_est_check_rd_state(est, mv_res);
#endif
  return dcost;
}

static int32_t od_mv_est_refine(od_mv_est_ctx *est, int log_dsz,
 int mv_res, const int *pattern_nsites, const od_pattern *pattern) {
  od_state *state;
  int32_t dcost;
  int nhmvbs;
  int nvmvbs;
  int vx;
  int vy;
  state = &est->enc->state;
  nhmvbs = state->nhmvbs;
  nvmvbs = state->nvmvbs;
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
   "Refining with displacements of %0g and 1/%i pel MV resolution.",
   (1 << log_dsz)*0.125, 1 << (3 - mv_res)));
  dcost = 0;
  for (vy = 0; vy <= nvmvbs; vy++) {
    if (est->row_counts[vy]) {
      dcost += od_mv_est_refine_row(est, vy, log_dsz, mv_res,
       pattern_nsites, pattern);
    }
  }
  for (vx = 0; vx <= nhmvbs; vx++) {
    if (est->col_counts[vx]) {
      dcost += od_mv_est_refine_col(est, vx, log_dsz, mv_res,
       pattern_nsites, pattern);
    }
  }
  return dcost;
}

/*STAGE 4: Sub-pel Refinement.*/

/*Stores the full-pel MVs for use by EPZS^2 in the next frame before sub-pel
   refinement.*/
void od_mv_est_update_fullpel_mvs(od_mv_est_ctx *est) {
  od_state *state;
  int nhmvbs;
  int nvmvbs;
  int vx;
  int vy;
  state = &est->enc->state;
  nhmvbs = state->nhmvbs;
  nvmvbs = state->nvmvbs;
  for (vy = 0; vy <= nvmvbs; vy++) {
    for (vx = 0; vx <= nhmvbs; vx++) {
      od_mv_grid_pt *mvg;
      od_mv_node *mv;
      mvg = state->mv_grid[vy] + vx;
      if (!mvg->valid) continue;
      mv = est->mvs[vy] + vx;
      mv->bma_mvs[0][mvg->ref][0] = mvg->mv[0] >> 3;
      mv->bma_mvs[0][mvg->ref][1] = mvg->mv[1] >> 3;
    }
  }
}

/*Sets the mv_rate of each node in the mesh, using the given MV resolution.
  Returns the change in rate.*/
int od_mv_est_update_mv_rates(od_mv_est_ctx *est, int mv_res) {
  od_state *state;
  int nhmvbs;
  int nvmvbs;
  int vx;
  int vy;
  int dr;
  state = &est->enc->state;
  nhmvbs = state->nhmvbs;
  nvmvbs = state->nvmvbs;
  dr = 0;
  for (vy = 0; vy <= nvmvbs; vy++) {
    for (vx = 0; vx <= nhmvbs; vx++) {
      od_mv_grid_pt *mvg;
      od_mv_node *mv;
      mvg = state->mv_grid[vy] + vx;
      if (!mvg->valid) continue;
      mv = est->mvs[vy] + vx;
      dr -= mv->mv_rate;
      mv->mv_rate = od_mv_est_bits(est, vx, vy, mv_res);
      dr += mv->mv_rate;
    }
  }
  return dr;
}

od_mv_est_ctx *od_mv_est_alloc(od_enc_ctx *enc) {
  od_mv_est_ctx *ret;
  ret = (od_mv_est_ctx *)malloc(sizeof(*ret));
  if (od_mv_est_init(ret, enc) < 0) {
    free(ret);
    return NULL;
  }
  return ret;
}

void od_mv_est_free(od_mv_est_ctx *est) {
  if (est != NULL) {
    od_mv_est_clear(est);
    free(est);
  }
}

void od_mv_est_reset_rd_block_state(od_mv_est_ctx *est,
 int vx, int vy, int log_mvb_sz) {
  od_state *state;
  int half_mvb_sz;
  state = &est->enc->state;
  half_mvb_sz = 1 << log_mvb_sz >> 1;
  if (log_mvb_sz > 0
   && state->mv_grid[vy + half_mvb_sz][vx + half_mvb_sz].valid) {
    od_mv_est_reset_rd_block_state(est, vx, vy, log_mvb_sz - 1);
    od_mv_est_reset_rd_block_state(est,
     vx + half_mvb_sz, vy, log_mvb_sz - 1);
    od_mv_est_reset_rd_block_state(est,
     vx, vy + half_mvb_sz, log_mvb_sz - 1);
    od_mv_est_reset_rd_block_state(est,
     vx + half_mvb_sz, vy + half_mvb_sz, log_mvb_sz - 1);
  }
  else {
    od_mv_node *block;
    int32_t sad;
    int oc;
    int s;
    block = est->mvs[vy] + vx;
    if (log_mvb_sz < OD_LOG_MVB_DELTA0) {
      int mask;
      int s1vx;
      int s1vy;
      int s3vx;
      int s3vy;
      mask = (1 << (log_mvb_sz + 1)) - 1;
      oc = !!(vx & mask);
      if (vy & mask) oc = 3 - oc;
      s1vx = vx + (OD_VERT_DX[(oc + 1) & 3] << log_mvb_sz);
      s1vy = vy + (OD_VERT_DY[(oc + 1) & 3] << log_mvb_sz);
      s3vx = vx + (OD_VERT_DX[(oc + 3) & 3] << log_mvb_sz);
      s3vy = vy + (OD_VERT_DY[(oc + 3) & 3] << log_mvb_sz);
      s = state->mv_grid[s1vy][s1vx].valid |
       state->mv_grid[s3vy][s3vx].valid << 1;
    }
    else {
      oc = 0;
      s = 3;
    }
    sad = od_mv_est_sad(est, vx, vy, oc, s, log_mvb_sz);
    block->sad = sad;
  }
}

void od_mv_subpel_refine(od_mv_est_ctx *est, int cost_thresh) {
  od_state *state;
  od_mv_grid_pt **grid;
  int32_t dcost;
  int32_t subpel_cost;
  int nhmvbs;
  int nvmvbs;
  int complexity;
  const int *pattern_nsites;
  const od_pattern *pattern;
  int mv_res;
  int best_mv_res;
  state = &est->enc->state;
  nhmvbs = state->nhmvbs;
  nvmvbs = state->nvmvbs;
  /*Save the fullpell MVs now for use by EPZS^2 on the next frame.
    We could also try rounding the results after refinement, I guess.
    I'm not sure it makes much difference*/
  od_mv_est_update_fullpel_mvs(est);
  complexity = est->enc->complexity;
  if (complexity >= OD_MC_SQUARE_SUBPEL_REFINEMENT_COMPLEXITY) {
    pattern_nsites = OD_SQUARE_NSITES;
    pattern = OD_SQUARE_SITES;
  }
  else {
    /*This speeds up each iteration by over 3x compared to a square pattern.*/
    pattern_nsites = OD_DIAMOND_NSITES;
    pattern = OD_DIAMOND_SITES;
  }
  do {
    dcost = od_mv_est_refine(est, 2, 2, pattern_nsites, pattern);
  }
  while (dcost < cost_thresh);
  for (best_mv_res = mv_res = 2; mv_res-- > est->mv_res_min;) {
    subpel_cost = od_mv_est_update_mv_rates(est, mv_res)*est->lambda;
    /*If the rate penalty for refining is small, bump the termination threshold
       down to make sure we actually get a decent improvement.
      We make sure not to let it get too small, however, so we're not here all
       day (a motion field of all (0, 0)'s would have a rate penalty of 0!).*/
    cost_thresh = OD_MAXI(cost_thresh,
     -OD_MAXI(subpel_cost, 16 << OD_ERROR_SCALE));
    OD_COPY(est->refine_grid[0], state->mv_grid[0],
     (nhmvbs + 1)*(nvmvbs + 1));
    do {
      dcost = od_mv_est_refine(est, mv_res, mv_res,
       pattern_nsites, pattern);
      subpel_cost += dcost;
    }
    while (dcost < cost_thresh);
    if (subpel_cost >= 0) {
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_INFO,
       "1/%i refinement FAILED:    dopt %7i", 1 << (3 - mv_res), subpel_cost));
      grid = est->refine_grid;
      est->refine_grid = state->mv_grid;
      state->mv_grid = grid;
      break;
    }
    else {
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_INFO,
       "1/%i refinement SUCCEEDED: dopt %7i", 1 << (3 - mv_res), subpel_cost));
      best_mv_res = mv_res;
    }
  }
  od_state_set_mv_res(state, best_mv_res);
}

void od_mv_est(od_mv_est_ctx *est, int lambda) {
  od_state *state;
  od_img_plane *iplane;
  int32_t dcost;
  int cost_thresh;
  int nhmvbs;
  int nvmvbs;
  int complexity;
  int use_satd;
  const int *pattern_nsites;
  const od_pattern *pattern;
  int log_mvb_sz;
  int pli;
  int i;
  int j;
  int xdec;
  int ydec;
  use_satd = est->enc->use_satd;
  state = &est->enc->state;
  nhmvbs = state->nhmvbs;
  nvmvbs = state->nvmvbs;
  xdec = od_plane_info_tab[est->enc->input_img.pixel_format].xdec[0];
  ydec = od_plane_info_tab[est->enc->input_img.pixel_format].ydec[0];

  iplane = est->enc->input_img.planes + 0;
  /*Sanitize user parameters*/
  est->level_min = OD_MINI(est->enc->params.mv_level_min,
   est->enc->params.mv_level_max);
  est->level_max = est->enc->params.mv_level_max;
  /*Rate estimations. Note that this does not depend on the previous frame: at
     this point, the probabilities have been reset by od_adapt_ctx_reset.*/
  for (i = 0; i < 5; i++) {
    for (j = 0; j < 16; j++) {
      est->mv_small_rate_est[i][j] = (int)((1 << OD_BITRES)
       *(OD_LOG2(est->enc->state.adapt.mv_small_cdf[i][15])
       - (OD_LOG2(est->enc->state.adapt.mv_small_cdf[i][j]
       - (j > 0 ? est->enc->state.adapt.mv_small_cdf[i][j - 1] : 0)))) + 0.5);
    }
  }
  /*If the luma plane is decimated for some reason, then our distortions will
     be smaller, so scale lambda appropriately.*/
  est->lambda = lambda >> (xdec + ydec);
  /*Compute termination thresholds for EPZS^2.*/
  for (log_mvb_sz = 0; log_mvb_sz < OD_NMVBSIZES; log_mvb_sz++) {
    est->thresh1[log_mvb_sz] =
     1 << 2*(log_mvb_sz + OD_LOG_MVBSIZE_MIN) >> (xdec + ydec);
  }
  /*If we're using the chroma planes, then our distortions will be larger.
    Compensate by increasing lambda and the termination thresholds.*/
  if (est->flags & OD_MC_USE_CHROMA) {
    for (pli = 1; pli < est->enc->input_img.nplanes; pli++) {
     xdec = od_plane_info_tab[est->enc->input_img.pixel_format].xdec[pli];
     ydec = od_plane_info_tab[est->enc->input_img.pixel_format].ydec[pli];
      est->lambda +=
       lambda >> (xdec + ydec + OD_MC_CHROMA_SCALE);
      for (log_mvb_sz = 0; log_mvb_sz < OD_NMVBSIZES; log_mvb_sz++) {
        est->thresh1[log_mvb_sz] += 1 << 2*(log_mvb_sz + OD_LOG_MVBSIZE_MIN) >>
         (xdec + ydec + OD_MC_CHROMA_SCALE);
      }
    }
  }
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG, "lambda: %i", est->lambda));
  for (log_mvb_sz = 0; log_mvb_sz < OD_NMVBSIZES; log_mvb_sz++) {
    est->thresh2_offs[log_mvb_sz] = est->thresh1[log_mvb_sz] >> 1;
  }
  /*Accelerated predictor weights.*/
  est->mvapw[OD_FRAME_PREV][0] = 0x20000;
  est->mvapw[OD_FRAME_PREV][1] = 0x10000;
  est->mvapw[OD_FRAME_GOLD][0] = 0x20000;
  est->mvapw[OD_FRAME_GOLD][1] = 0x10000;
  /*TODO: Constant velocity predictor weight.*/
#if defined(OD_DUMP_IMAGES) && defined(OD_ANIMATE)
  /*Set some initial state.
    This would get reset eventually by the algorithm in a more convenient
     place, but is needed earlier by the visualization.*/
  if (daala_granule_basetime(&est->enc->state, est->enc->state.cur_time) ==
   ANI_FRAME) {
    od_state_mvs_clear(&est->enc->state);
  }
#endif
  /*Use SAD for stages here after.*/
  est->compute_distortion = od_enc_sad;
  od_mv_est_init_mvs(est, OD_FRAME_PREV, 1);
  /* At very high lambdas, the signaling overhead of multiref is too high. */
  if (lambda < 150) {
    if (state->ref_imgi[OD_FRAME_GOLD] >= 0) {
      od_mv_est_init_mvs(est, OD_FRAME_GOLD, 0);
    }
  }
  od_mv_est_decimate(est);
  /*This threshold is somewhat arbitrary.
    Chen and Willson use 6000 (with SSD as an error metric).
    We would like something more dependent on the frame size.
    For CIF, there are a maximum of 6992 vertices in the mesh, which is pretty
     close to 6000.
    With a SAD error metric like we use, the square root of 6000 would be a
     more appropriate value, however that gives a PSNR improvement of less than
     0.01 dB, and requires almost twice as many iterations to achieve.*/
  cost_thresh = -nhmvbs*nvmvbs << OD_ERROR_SCALE;
  complexity = est->enc->complexity;
  if (complexity >= OD_MC_SQUARE_REFINEMENT_COMPLEXITY) {
    pattern_nsites = OD_SQUARE_NSITES;
    pattern = OD_SQUARE_SITES;
  }
  else {
    /*This speeds up each iteration by over 3x compared to a square pattern.*/
    pattern_nsites = OD_DIAMOND_NSITES;
    pattern = OD_DIAMOND_SITES;
  }
  do {
    dcost = 0;
    /*Logarithmic (telescoping) search.
      This is 3x more expensive than basic refinement, but can help escape
       local minima.*/
    if (complexity >= OD_MC_LOGARITHMIC_REFINEMENT_COMPLEXITY) {
      dcost += od_mv_est_refine(est, 5, 2, pattern_nsites, pattern);
      dcost += od_mv_est_refine(est, 4, 2, pattern_nsites, pattern);
    }
    dcost += od_mv_est_refine(est, 3, 2, pattern_nsites, pattern);
  }
  while (dcost < cost_thresh);
  if (use_satd) {
    /* The two #defines below apply to sub-pel ME only. */
# define OD_ME_SATD_THRESH_SCALE (0.7) /* 1.0 means the same as SAD. */
# define OD_ME_SATD_LAMBDA_SCALE (0.6)
    est->compute_distortion = od_enc_satd;/* Use SATD from this stage. */
    est->lambda = (int)(est->lambda*OD_ME_SATD_LAMBDA_SCALE);
    cost_thresh = (int)(cost_thresh*OD_ME_SATD_THRESH_SCALE);
    /* Reset sad values for each ME block since those are used as base
        when measuring the change in ME error in the next stage. */
    {
      int vx;
      int vy;
      for (vy = 0; vy < nvmvbs; vy += OD_MVB_DELTA0) {
        for (vx = 0; vx < nhmvbs; vx += OD_MVB_DELTA0) {
          od_mv_est_reset_rd_block_state(est, vx, vy, OD_LOG_MVB_DELTA0);
        }
      }
    }
  }
  od_mv_subpel_refine(est, cost_thresh);
  od_restore_fpu(state);
}
