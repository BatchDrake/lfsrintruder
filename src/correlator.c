/*

  correlator.c: Fast correlator
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

#include <complex.h>
#include <math.h>

#include "correlator.h"

#include <string.h>
#include <sys/stat.h>

PTR_LIST_EXTERN(lfsrdesc_t, desc);

void
correlator_destroy(correlator_t *self)
{
  unsigned int i;

  for (i = 0; i < self->candidate_count; ++i)
    if (self->candidate_list[i] != NULL)
      free(self->candidate_list[i]);

  if (self->candidate_list != NULL)
    free(self->candidate_list);

  if (self->data != NULL)
    free(self->data);

  if (self->data_freq != NULL)
    free(self->data_freq);

  if (self->seq_freq != NULL)
    free(self->seq_freq);

  if (self->xcorr != NULL)
    free(self->xcorr);

  if (self->fft_plan_inv != NULL)
    fftwf_destroy_plan(self->fft_plan_inv);

  if (self->fft_plan != NULL)
    fftwf_destroy_plan(self->fft_plan);

  free(self);
}

static void
correlator_attempt_save(const char *path, const uint8_t *data, size_t N)
{
  FILE *fp;

  if ((fp = fopen(path, "w")) == NULL)
    return;

  while (N-- != 0)
    fputc(*data++ + '0', fp);

  fclose(fp);
}

static BOOL
correlator_save_candidate(
    const correlator_t *self,
    const uint8_t *seq,
    unsigned int offset,
    const char *name)
{
  char *path = NULL;
  FILE *fp = NULL;
  unsigned int i = 0;
  unsigned int hw = 0;
  unsigned int max_seq = 0;
  unsigned int curr_seq = 0;
  unsigned int flip_count = 0;
  uint8_t prev = 0;

  uint8_t b;
  BOOL ok = FALSE;

  if (access("candidates/", F_OK) == -1)
    TRY(mkdir("candidates", 0755) != -1);

  TRY(
      path = strbuild(
          "candidates/unscrambled-off%d-%s.log", name));

  TRY(fp = fopen(path, "w"));

  for (i = 0; i < self->N; ++i) {
    b = (seq[(i + offset) % self->N] ^ self->data[i]);
    if (i > 0) {
      if (b == prev) {
        if (++curr_seq > max_seq)
          max_seq = curr_seq;
      } else {
        ++flip_count;
        curr_seq = 0;
      }
    }

    prev = b;
    hw += b;
    fputc('0' + b, fp);
  }

  _DEBUG("  Hamming weight:   %d\n", hw);
  _DEBUG("  Longest sequence: %d\n", max_seq);
  _DEBUG("  Bit flip count:   %d\n", flip_count);

  ok = TRUE;

fail:
  if (path != NULL)
    free(path);

  if (fp != NULL)
    fclose(fp);

  return ok;
}

BOOL
correlator_walk_candidates(
    correlator_t *self,
    BOOL (*callback) (const struct correlator_candidate *, void *),
    void *private)
{
  unsigned int i;

  for (i = 0; i < self->candidate_count; ++i)
    if (!(callback) (self->candidate_list[i], private))
      return FALSE;

  return TRUE;
}

BOOL
correlator_register_candidate(
    correlator_t *self,
    lfsrdesc_t *desc,
    unsigned int offset)
{
  struct correlator_candidate *candidate = NULL;

  ALLOCATE(candidate, struct correlator_candidate);

  candidate->desc = desc;
  candidate->offset = offset;
  candidate->phase = offset % lfsrdesc_get_cycle_len(desc);

  TRY(PTR_LIST_APPEND_CHECK(self->candidate, candidate) != -1);

  return TRUE;

fail:
  if (candidate != NULL)
    free(candidate);

  return FALSE;
}

BOOL
correlator_run(correlator_t *self)
{
  unsigned int i, j;
  unsigned int max_j;
  float amp, max;
  char *poly = NULL;
  uint8_t *seq = NULL;
  float K = 1.f / self->N;
  BOOL ok = FALSE;

  _DEBUG("Running against %d polynomials\n", desc_count);

  self->best_score = 0;

  /* Run correlator on each polynomial */
  for (i = 0; i < desc_count; ++i) {
    /* Get polynomial desc */
    TRY(poly = lfsrdesc_get_poly(desc_list[i]));

    /* Generate float sequence */
    TRY(seq = lfsrdesc_generate(desc_list[i], self->N));

    for (j = 0; j < self->N; ++j)
      self->seq_freq[j] = 2 * K * (seq[j] - .5);

    fftwf_execute(self->fft_plan); /* Change to frequency */

    /* Multiply by data in frequency domain  */
    for (j = 0; j < self->N; ++j)
      self->seq_freq[j] *= conj(self->data_freq[j]);

    /* Compute inverse FFT */
    fftwf_execute(self->fft_plan_inv);

    max = 0;
    max_j = 0;
    for (j = 0; j < self->N; ++j) {
      amp = creal(self->xcorr[j] * conj(self->xcorr[j]));
      if (amp > max) {
        max = amp;
        max_j = j;
      }
    }

    if (max > self->best_score) {
      TRY(correlator_register_candidate(self, desc_list[i], max_j));
      self->best_score = max;

      _DEBUG(
          "Best score: %6.2f%% in %-5d (polynomial %s)\n",
          100.f * max,
          max_j,
          poly);

      correlator_save_candidate(self, seq, max_j, poly);
    }

    free(poly);
    poly = NULL;

    free(seq);
    seq = NULL;
  }

  ok = TRUE;

fail:
  if (poly != NULL)
    free(poly);

  if (seq != NULL)
    free(seq);

  return ok;
}

correlator_t *
correlator_new(const uint8_t *data, size_t N)
{
  correlator_t *new = NULL;
  fftwf_plan plan = NULL;
  BOOL ok = FALSE;
  float K;
  unsigned int i;

  ALLOCATE(new, correlator_t);

  ALLOCATE_MANY(new->data, N, uint8_t);
  ALLOCATE_FFT(new->data_freq, N);
  ALLOCATE_FFT(new->seq_freq, N);
  ALLOCATE_FFT(new->xcorr, N);

  memcpy(new->data, data, N * sizeof(uint8_t));

  new->N = N;

  correlator_attempt_save("input.log", data, N);

  /*
   * Compute some FFTs
   */

  K = 1. / N;
  for (i = 0; i < N; ++i)
    new->data_freq[i] = 2 * K * (data[i] - .5);

  _DEBUG("Computing FFT of data (%d bins)\n", N);

  TRY(plan = fftwf_plan_dft_1d(
      new->N,
      new->data_freq,
      new->data_freq,
      FFTW_FORWARD,
      FFTW_ESTIMATE));

  _DEBUG("Computing...\n");

  fftwf_execute(plan); /* In xcorr_freq: FFT of data */
  fftwf_destroy_plan(plan);
  plan = NULL;

  _DEBUG("Done\n");

  TRY(new->fft_plan = fftwf_plan_dft_1d(
      new->N,
      new->seq_freq,
      new->seq_freq,
      FFTW_FORWARD,
      FFTW_ESTIMATE));

  TRY(new->fft_plan_inv = fftwf_plan_dft_1d(
      new->N,
      new->seq_freq,
      new->xcorr,
      FFTW_BACKWARD,
      FFTW_ESTIMATE));

  ok = TRUE;

fail:
  if (plan != NULL)
    fftwf_destroy_plan(plan);

  if (!ok && new != NULL) {
    correlator_destroy(new);
    new = NULL;
  }

  return new;
}
