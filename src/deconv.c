/*

  deconv.c: Decode convolutional code
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

#include "viterbi.h"

void
usage(const char *a0)
{
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "  %s K poly1 [poly2 [...]]\n", a0);
}

struct error_info {
  int total;
  int failed;
};

static BOOL
on_data(
    const uint8_t *path,
    unsigned int len,
    unsigned int errors,
    void *private)
{
  unsigned int i;
  struct error_info *info = (struct error_info *) private;

  if (errors > len / (VITERBI_TRELLIS_LENGTH - 1)) {
    fprintf(stderr, "Stream is extremely corrupted here! (%d errors)\n", errors);
    ++info->failed;
  } else if (errors > 0) {
    fprintf(
        stderr,
        "Warning: %d unexplainable bits in %d codewords\n",
        errors,
        len);
    ++info->failed;
  }

  for (i = 0; i < len; ++i)
    putchar('0' + path[i]);

  ++info->total;

  return TRUE;
}

int
main(int argc, char *argv[])
{
  unsigned int n, K;
  unsigned int i;
  uint8_t *c;
  uint32_t codeword;
  uint32_t *polies;
  viterbi_t *viterbi;
  struct error_info info = {0, 0};
  struct viterbi_params params;

  if (argc < 3) {
    fprintf(stderr, "%s: wrong number of arguments\n", argv[0]);
    usage(argv[0]);
    exit(EXIT_FAILURE);
  }

  if (sscanf(argv[1], "%u", &K) != 1) {
    fprintf(
        stderr,
        "%s: invalid constraint length \"%s\"\n",
        argv[0],
        argv[1]);
    usage(argv[0]);
    exit(EXIT_FAILURE);
  }

  n = argc - 2;

  ALLOCATE_MANY(c,      n, uint8_t);
  ALLOCATE_MANY(polies, n, unsigned int);

  for (i = 0; i < n; ++i) {
    if (sscanf(argv[2 + i], "%u", polies + i) != 1) {
      fprintf(stderr, "%s: invalid polinomial\n", argv[0], argv[2 + i]);
      usage(argv[0]);
      exit(EXIT_FAILURE);
    }
  }

  params.K = K;
  params.n = n;
  params.poly = polies;
  params.on_data = on_data;
  params.private = &info;

  CONSTRUCT(viterbi, viterbi, &params);

  while (fread(c, sizeof(uint8_t), n, stdin) == n) {
    codeword = 0;
    for (i = 0; i < n; ++i) {
      if (c[i] == '1') {
        codeword |= 1 << (n - i - 1);
      } else if (c[i] != '0') {
        fprintf(
            stderr,
            "%s: invalid character \0%o in input\n",
            argv[0],
            c[i]);
        goto fail;
      }
    }

    if (!viterbi_feed(viterbi, codeword)) {
      fprintf(stderr, "%s: Viterbi decoder refused to continue\n", argv[0]);
      goto fail;
    }
  }

  if (info.total == info.failed) {
    fprintf(stderr, "%s: all tracebacks failed. Decoding failed\n");
    goto fail;
  }

  exit(EXIT_SUCCESS);

fail:
  exit(EXIT_FAILURE);
}
