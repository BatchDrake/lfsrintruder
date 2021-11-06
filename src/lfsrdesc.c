/*

  lfsrdesc.h: LFSR description, as read from polynomial database file
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

#include "lfsrdesc.h"

#include <string.h>
#include <ctype.h>

PTR_LIST(lfsrdesc_t, desc);

void
lfsrdesc_destroy(lfsrdesc_t *self)
{
  if (self->poly != NULL)
    free(self->poly);

  if (self->lfsr != NULL)
    lfsr_destroy(self->lfsr);

  free(self);
}

lfsrdesc_t *
lfsrdesc_new(const unsigned int *polynomial, size_t poly_len)
{
  lfsrdesc_t *self = NULL;

  ALLOCATE(self, lfsrdesc_t);

  ALLOCATE_MANY(self->poly, poly_len, unsigned int);

  memcpy(self->poly, polynomial, poly_len * sizeof(unsigned int));
  self->poly_size = poly_len;

  CONSTRUCT(self->lfsr, lfsr, polynomial, poly_len);

  return self;

fail:
  if (self != NULL)
    lfsrdesc_destroy(self);

  return NULL;
}

uint8_t *
lfsrdesc_generate(lfsrdesc_t *self, size_t len)
{
  uint8_t *alloc = NULL;
  unsigned int i;

  ALLOCATE_MANY(alloc, len, uint8_t);

  lfsr_reset(self->lfsr);

  /* Empty pipeline */
  for (i = 0; i < 64; ++i)
    lfsr_scramble(self->lfsr, 0);

  for (i = 0; i < len; ++i)
    alloc[i] = lfsr_scramble(self->lfsr, 0);

  return alloc;

fail:
  if (alloc != NULL)
    free(alloc);

  return NULL;
}

char *
lfsrdesc_get_poly(const lfsrdesc_t *self)
{
  return lfsr_get_poly(self->lfsr);
}

BOOL
lfsrdesc_load_from_file(const char *path)
{
  FILE *fp = NULL;
  lfsrdesc_t *desc = NULL;
  BOOL ok = FALSE;
  BOOL primitive = TRUE;
  char *line = NULL;
  char *p;
  arg_list_t *args = NULL;
  unsigned int *taps = NULL, *tmp;
  unsigned int max_tap_len = 0;
  unsigned int i;

  TRY(fp = fopen(path, "r"));

  while ((line = fread_line(fp)) != NULL) {
    p = line;

    while (isspace(*p))
      ++p;

    if (*p == '\0') {
      primitive = TRUE; /* Empty line: polynomials are primitive again */
    } else if (strstr(p, "non-primitive") != NULL) {
      primitive = FALSE;
    } else if (*line != '#' && primitive) {
      TRY(args = csv_split_line(p));
      if (args->al_argc > max_tap_len) {
        TRY(tmp = realloc(taps, sizeof(unsigned int) * args->al_argc));
        max_tap_len = args->al_argc;
        taps = tmp;
      }

      for (i = 0; i < args->al_argc; ++i)
        TRY(sscanf(args->al_argv[i], "%u", taps + i) == 1);

      CONSTRUCT(desc, lfsrdesc, taps, args->al_argc);

      TRY(PTR_LIST_APPEND_CHECK(desc, desc) != -1);

      desc = NULL;

      free_al(args);
      args = NULL;
    }

    free(line);
    line = NULL;
  }

  ok = TRUE;

fail:
  if (desc != NULL)
    lfsrdesc_destroy(desc);

  if (taps != NULL)
    free(taps);

  if (args != NULL)
    free_al(args);

  if (line != NULL)
    free(line);

  if (fp != NULL)
    fclose(fp);

  return ok;
}

