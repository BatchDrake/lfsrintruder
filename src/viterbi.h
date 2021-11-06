/*

  viterbi.h: Viterbi decoder implementation for 1/n codes
  Copyright (C) 2016 Gonzalo José Carracedo Carballal

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as
  published by the Free Software Foundation, either version 3 of the
  License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this program.  If not, see
  <http://www.gnu.org/licenses/>

*/


#ifndef _VITERBI_H
#define _VITERBI_H

#include "types.h"

#define VITERBI_MAX_K          10
#define VITERBI_TRELLIS_LENGTH 5
#define VITERBI_MAX_N          16

#define VITERBI_INFINITY 0xffffffff

/*
 * Polynomial notation:
 *   MSB: current input
 *   LSB: oldest bit (the next one to be discarded)
 */
struct viterbi_params {
  unsigned int n;
  unsigned int K;
  const uint32_t *poly;

  BOOL (*on_data) (
      const uint8_t *path,
      unsigned int len,
      unsigned int errors,
      void *private);

  void *private;
};

struct viterbi_node {
  uint32_t prev; /* Previous path */
  uint32_t path_dist; /* Path distance so far */
};

struct viterbi {
  struct viterbi_params params;
  unsigned int state_count;
  unsigned int trellis_length;
  unsigned int p;
  uint32_t   code_mask;
  uint32_t   constraint_mask;
  uint8_t   *sequence; /* Reconstructed sequence */
  uint32_t  *code_dict; /* Converts state to codeword */
  uint32_t   best; /* Best state so far */
  struct viterbi_node **trellis;
};

typedef struct viterbi viterbi_t;

viterbi_t *viterbi_new(const struct viterbi_params *params);
BOOL viterbi_feed(viterbi_t *self, uint32_t codeword);
void viterbi_destroy(viterbi_t *self);

#endif /* _VITERBI_H */
