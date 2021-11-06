/*

  correlator.h: Fast correlator
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

#ifndef _CORRELATOR_H
#define _CORRELATOR_H

#include "lfsrdesc.h"

#include <fftw3.h>

struct correlator_candidate {
  lfsrdesc_t *desc;
  unsigned int offset;
  unsigned int phase;
};

struct correlator {
  uint8_t *data;               /* Copy of input data */
  fftwf_complex *data_freq; /* Data in frequency domain */
  fftwf_complex *seq_freq;   /* Sequence in frequency domain */
  fftwf_complex *xcorr;      /* Computed on each run */

  size_t N;

  fftwf_plan fft_plan;     /* FFT(seq_freq) --> seq_freq */
  fftwf_plan fft_plan_inv; /* IFFT(seq_freq) --> xcorr */

  PTR_LIST(struct correlator_candidate, candidate);

  float best_score;
};

typedef struct correlator correlator_t;

void correlator_destroy(correlator_t *self);

BOOL correlator_walk_candidates(
    correlator_t *self,
    BOOL (*callback) (const struct correlator_candidate *, void *),
    void *private);

BOOL correlator_run(correlator_t *corr);

correlator_t *correlator_new(const uint8_t *data, size_t N);

#endif /* _CORRELATOR_H */

