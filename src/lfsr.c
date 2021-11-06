/*

  lfsr.c: Multiplicative descrambler
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "types.h"
#include "lfsr.h"

void
lfsr_destroy(lfsr_t *self)
{
  free(self);
}

char *
lfsr_get_poly(const lfsr_t *self)
{
  unsigned int i;
  char *prev = NULL;
  char *poly = NULL;

  for (i = 63; i >= 1; --i)
    if ((self->mask & (1ull << i)) != 0) {
      TRY(poly = strbuild("%sx^%d + ", prev == NULL ? "" : prev, i));
      if (prev != NULL)
        free(prev);
      prev = poly;
    }

  TRY(poly = strbuild("%s1", prev == NULL ? "" : prev));

  if (prev != NULL)
    free(prev);

  return poly;

fail:
  if (prev != NULL)
    free(prev);

  return NULL;
}

void
lfsr_reset(lfsr_t *self)
{
  self->reg = 0xffffffffffffffffull;
}

lfsr_t *
lfsr_new(const unsigned int *taps, unsigned int tap_len)
{
  lfsr_t *self = NULL;
  unsigned int i;

  if ((self = calloc(1, sizeof(lfsr_t))) == NULL)
    goto fail;

  for (i = 0; i < tap_len; ++i) {
    if (taps[i] >= LFSR_MAX_TAPS) {
      ERROR("Invalid tap %d\n", taps[i]);
      goto fail;
    }

    self->mask |= 1ull << taps[i];
    if (self->len < taps[i])
      self->len = taps[i];
  }

  lfsr_reset(self);

  self->cycle_len = (1ull << self->len) - 1;

  --self->len;

  return self;

fail:
  if (self != NULL)
    lfsr_destroy(self);

  return NULL;
}

static inline uint8_t
lfsr_core(lfsr_t *self, uint8_t input, unsigned int direction)
{
  uint8_t x = input & 1;
  uint8_t y = (popcount64(self->reg & self->mask) & 1) ^ x;
  uint8_t newbit = direction ? x : y;
  uint8_t output = direction ? y : self->reg & 1;

  self->reg = (self->reg >> 1) | (newbit << self->len);

  return y;
}

uint8_t
lfsr_scramble(lfsr_t *self, uint8_t input)
{
  return lfsr_core(self, input, 0);
}

uint8_t
lfsr_descramble(lfsr_t *self, uint8_t input)
{
  return lfsr_core(self, input, 1);
}

