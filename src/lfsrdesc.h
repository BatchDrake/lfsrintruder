/*

  lfsrdesc.h: LFSR description, as read from file
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


#ifndef _LFSRDESC_H
#define _LFSRDESC_H

#include "types.h"
#include "lfsr.h"

struct lfsrdesc {
  unsigned int *poly;
  size_t poly_size;
  lfsr_t *lfsr;
};

typedef struct lfsrdesc lfsrdesc_t;

static inline uint64_t
lfsrdesc_get_cycle_len(const lfsrdesc_t *desc)
{
  return lfsr_get_cycle_len(desc->lfsr);
}

lfsrdesc_t *lfsrdesc_new(const unsigned int *poly, size_t poly_size);
uint8_t *lfsrdesc_generate(lfsrdesc_t *desc, size_t len);
char *lfsrdesc_get_poly(const lfsrdesc_t *self);
void lfsrdesc_destroy(lfsrdesc_t *);

BOOL lfsrdesc_load_from_file(const char *path);

#endif /* _LFSRDESC_H */

