/*

  viterbi.c: Viterbi decoder implementation
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

#include "viterbi.h"

#include <assert.h>

void
viterbi_destroy(viterbi_t *self)
{
  unsigned int i;

  if (self->trellis != NULL) {
    for (i = 0; i < self->state_count; ++i)
      if (self->trellis[i] != NULL)
        free(self->trellis[i]);

    free(self->trellis);
  }

  if (self->sequence != NULL)
    free(self->sequence);

  if (self->code_dict != NULL)
    free(self->code_dict);

  free(self);
}

static void
print_binary(uint32_t word, unsigned int width)
{
  unsigned int i;

  for (i = 0; i < width; ++i)
    putchar('0' + ((word >> (width - i - 1)) & 1));
}

BOOL
viterbi_feed(viterbi_t *self, uint32_t codeword)
{
  uint32_t i;
  uint32_t min_path_dist;
  uint32_t dist, dist1;
  uint32_t path_dist, path_dist1;
  uint32_t prev;

  codeword &= self->code_mask;

  /* Traceback test */
  if (++self->p == self->trellis_length) {
    prev = self->best;
    do {
      --self->p;
      prev = self->trellis[prev & self->constraint_mask][self->p].prev;
      self->sequence[self->p] = prev >> (self->params.K - 1);
    } while (self->p != 0);

    /* Reset first state */
    for (i = 0; i < self->state_count; ++i)
      self->trellis[i][0].path_dist = i == self->best ? 0 : VITERBI_INFINITY;

    self->trellis[self->best][0].prev =
        self->trellis[self->best][self->trellis_length - 1].prev;

    /* Call data callback */
    if (!(self->params.on_data) (
        self->sequence,
        self->trellis_length - 1,
        self->trellis[self->best][self->trellis_length - 1].path_dist,
        self->params.private))
      return FALSE;

    self->p = 1;
  }

  min_path_dist = VITERBI_INFINITY;

  for (i = 0; i < self->state_count; ++i) {
    /* In this step: previous possible states are:
     *   Sp0 = ((i << 1) | 0) & constraint_mask: if oldest bit was 0
     *   Sp1 = ((i << 1) | 1) & constraint_mask: if oldest bit was 1
     *
     * We inspect which Sp has the smallest Hamming distance, and choose
     * it as previous.
     */
    prev = i << 1;

    /* Distance to codeword when oldest bit was 0 */
    dist  = popcount64(self->code_dict[prev | 0] ^ codeword);
    if ((path_dist = self->trellis[prev & self->constraint_mask][self->p - 1].path_dist)
        != VITERBI_INFINITY)
      path_dist += dist;

    dist1 = popcount64(self->code_dict[prev | 1] ^ codeword);
    if ((path_dist1 = self->trellis[(prev | 1) & self->constraint_mask][self->p - 1].path_dist)
        != VITERBI_INFINITY)
      path_dist1 += dist1;

    if (path_dist > path_dist1) {
      prev |= 1; /* Yes */
      path_dist = path_dist1;
    }

    self->trellis[i][self->p].path_dist = path_dist;
    self->trellis[i][self->p].prev      = prev;

    /* Update minimum */
    if (path_dist < min_path_dist) {
      min_path_dist = path_dist;
      self->best = i;
    }
  }

  return TRUE;
}

viterbi_t *
viterbi_new(const struct viterbi_params *params)
{
  viterbi_t *self = NULL;
  unsigned int i, j;
  uint32_t codeword;

  assert(params->K <= VITERBI_MAX_K);
  assert(params->n <= VITERBI_MAX_N);

  ALLOCATE(self, viterbi_t);

  self->params = *params;
  self->state_count = 1 << (params->K - 1);
  self->code_mask = (1 << params->n) - 1;
  self->constraint_mask = self->state_count - 1;
  self->trellis_length = params->K * VITERBI_TRELLIS_LENGTH;

  ALLOCATE_MANY(self->sequence,  self->trellis_length, uint8_t);
  ALLOCATE_MANY(self->code_dict, self->state_count * 2, uint32_t);
  ALLOCATE_MANY(self->trellis,   self->state_count, struct viterbi_node *);

  for (i = 0; i < self->state_count; ++i)
    ALLOCATE_MANY(
        self->trellis[i],
        self->trellis_length,
        struct viterbi_node);

  for (i = 0; i < self->state_count; ++i) {
    /* Populate trellis initial state */
    self->trellis[i][0].path_dist = i == 0 ? 0 : VITERBI_INFINITY;

    /* Populate code dictionary */
    for (j = 0; j < params->n; ++j)
      self->code_dict[i] |=
          (popcount64(params->poly[j] & i) & 1) << (params->n - j - 1);

    codeword = i | self->state_count;
    for (j = 0; j < params->n; ++j)
      self->code_dict[codeword] |=
          (popcount64(params->poly[j] & codeword) & 1) << (params->n - j - 1);
  }

  return self;

fail:
  if (self != NULL)
    viterbi_destroy(self);

  return NULL;
}
