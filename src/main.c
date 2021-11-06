/*
 * main.c: entry point for lfsrintruder
 * Creation date: Tue Jan 22 10:15:41 2019
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <getopt.h>
#include <sys/stat.h>

#include "correlator.h"

#define OUTPUT_DIRECTORY "descrambled"

struct lfsr_params_hit {
  unsigned int offset;
  unsigned int hits;
};

struct lfsr_hit {
  lfsrdesc_t *desc;
  unsigned int hits;
  unsigned int max_offset_hits;
  PTR_LIST(struct lfsr_params_hit, params_hit);
};

PTR_LIST(struct lfsr_hit, hit);

void
lfsr_hit_destroy(struct lfsr_hit *hit)
{
  unsigned int i;

  for (i = 0; i < hit->params_hit_count; ++i)
    if (hit->params_hit_list[i] != NULL)
      free(hit->params_hit_list[i]);

  if (hit->params_hit_list != NULL)
    free(hit->params_hit_list);

  free(hit);
}

BOOL
lfsr_hit_push(struct lfsr_hit *self, unsigned int offset)
{
  unsigned int i;
  struct lfsr_params_hit *hit = NULL;

  for (i = 0; i < self->params_hit_count; ++i)
    if (self->params_hit_list[i]->offset == offset) {
      hit = self->params_hit_list[i];
      break;
    }

  if (hit == NULL) {
    ALLOCATE(hit, struct lfsr_params_hit);
    hit->offset = offset;
    TRY(PTR_LIST_APPEND_CHECK(self->params_hit, hit) != -1);
  }

  ++self->hits;
  ++hit->hits;

  if (hit->hits > self->max_offset_hits)
    self->max_offset_hits = hit->hits;

  return TRUE;

fail:
  if (hit != NULL)
    free(hit);

  return FALSE;
}

struct lfsr_hit *
lfsr_hit_new(lfsrdesc_t *desc)
{
  struct lfsr_hit *self = NULL;

  ALLOCATE(self, struct lfsr_hit);

  self->desc = desc;

  return self;

fail:
  if (self != NULL)
    lfsr_hit_destroy(self);

  return NULL;
}

struct lfsr_hit *
lfsr_hit_lookup(const lfsrdesc_t *desc)
{
  unsigned int i;

  for (i = 0; i < hit_count; ++i)
    if (hit_list[i]->desc == desc)
      return hit_list[i];

  return NULL;
}

BOOL
lfsr_hit_assert(
    lfsrdesc_t *desc,
    unsigned int offset)
{
  struct lfsr_hit *hit, *new_hit = NULL;
  unsigned int i = 0;

  if ((hit = lfsr_hit_lookup(desc)) == NULL) {
    CONSTRUCT(new_hit, lfsr_hit, desc);
    TRY(PTR_LIST_APPEND_CHECK(hit, new_hit) != -1);
    hit = new_hit;
    new_hit = NULL;
  }

  TRY(lfsr_hit_push(hit, offset));

  return TRUE;

fail:
  if (new_hit != NULL)
    lfsr_hit_destroy(new_hit);

  return FALSE;
}

static BOOL
on_candidate(const struct correlator_candidate *candidate, void *private)
{
  /* Record only polynomials whose cycle length is at least 31 */
  if (lfsrdesc_get_cycle_len(candidate->desc) >= 16)
    TRY(lfsr_hit_assert(candidate->desc, candidate->phase));

  return TRUE;

fail:
  return FALSE;
}

BOOL
lfsr_hit_descramble_file(
    const struct lfsr_hit *candidate,
    const char *input,
    unsigned int index,
    unsigned int offset)
{
  char *path = NULL;
  FILE *ofp = NULL;
  FILE *fp = NULL;

  uint8_t *seq = NULL;
  uint8_t c;
  unsigned int len = lfsrdesc_get_cycle_len(candidate->desc);
  unsigned int p = offset % len;
  BOOL ok = FALSE;

  if (access(OUTPUT_DIRECTORY, F_OK) == -1)
    TRY_EXCEPT(
        mkdir(OUTPUT_DIRECTORY, 0755) == -1,
        fprintf(
            stderr,
            "Failed to create output directory %s: %s\n",
            OUTPUT_DIRECTORY,
            strerror(errno)));

  TRY(path = strbuild("%s/descrambled-%06d.log", OUTPUT_DIRECTORY, index));

  TRY_EXCEPT(
      fp = fopen(input, "r"),
      fprintf(
          stderr,
          "Failed to open %s for reading: %s\n",
          input,
          strerror(errno)));

  TRY_EXCEPT(
      ofp = fopen(path, "w"),
      fprintf(
          stderr,
          "Failed to open %s for writing: %s\n",
          path,
          strerror(errno)));

  /* Generate a cycle */
  TRY(seq = lfsrdesc_generate(candidate->desc, len));

  while (fread(&c, 1, 1, fp) == 1)
    if (c == '0' || c == '1') {
      c -= '0';
      c ^= seq[p++];
      if (p == len)
        p = 0;
      fputc(c + '0', ofp);
    }
  ok = TRUE;

fail:
  if (path != NULL)
    free(path);

  if (ofp != NULL)
    fclose(ofp);

  if (fp != NULL)
    fclose(fp);

  if (seq != NULL)
    free(seq);

  return ok;
}

int
main(int argc, char *argv[], char *envp[])
{
  FILE *fp = NULL;
  correlator_t *corr = NULL;
  char *buffer = NULL;
  uint8_t c;
  unsigned int p = 0;
  unsigned int i, j;
  unsigned int files = 0;
  unsigned int max_hits = 0;
  unsigned int best_offset = 0;
  struct lfsr_hit *best_hit = NULL;
  char *poly;

  struct stat sbuf;

  if (argc < 2) {
    fprintf(
        stderr,
        "Wrong number of arguments. Usage:\n\t%s file1.log [file2.log [...]]\n",
        argv[0]);
    exit(EXIT_FAILURE);
  }

  if (!lfsrdesc_load_from_file("all-irredpoly.txt")) {
    fprintf(stderr, "%s: cannot load polynomials\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  for (i = 1; i < argc; ++i) {
    if (stat(argv[i], &sbuf) == -1) {
      fprintf(
          stderr,
          "%s: cannot stat %s: %s\n",
          argv[0],
          argv[1],
          strerror(errno));

      goto cleanup;
    }

    if (sbuf.st_size == 0) {
      fprintf(stderr, "%s: file %s is empty, skipping...\n", argv[0], argv[i]);
      goto cleanup;
    }

    ALLOCATE_MANY(buffer, sbuf.st_size, uint8_t);

    if ((fp = fopen(argv[i], "rb")) == NULL) {
      fprintf(
          stderr,
          "%s: cannot open %s: %s\n",
          argv[0],
          argv[1],
          strerror(errno));
      goto cleanup;
    }

    p = 0;
    while (fread(&c, 1, 1, fp) == 1)
      if (c == '0' || c == '1')
        buffer[p++] = c - '0';

    if ((corr = correlator_new(buffer, p)) == NULL) {
      fprintf(stderr, "%s: cannot correlate %d bytes\n", argv[0], p);
      goto cleanup;
    }

    TRY(correlator_run(corr));

    /* Everything went alright */
    TRY(correlator_walk_candidates(corr, on_candidate, NULL));

    ++files;

cleanup:
    if (corr != NULL) {
      correlator_destroy(corr);
      corr = NULL;
    }

    if (fp != NULL) {
      fclose(fp);
      fp = NULL;
    }

    if (buffer != NULL) {
      free(buffer);
      buffer = NULL;
    }
  }

  for (i = 0; i < hit_count; ++i) {
    if (files == 1 || hit_list[i]->hits > 1) {
      TRY(poly = lfsrdesc_get_poly(hit_list[i]->desc));
      printf("%3d/%d hits: %s\n", hit_list[i]->hits, files, poly);
      free(poly);
      for (j = 0; j < hit_list[i]->params_hit_count; ++j)
        printf(
            "      Offset %4d with %3d hits\n",
            hit_list[i]->params_hit_list[j]->offset,
            hit_list[i]->params_hit_list[j]->hits);

      putchar(10);

    }


    if (max_hits < hit_list[i]->max_offset_hits) {
      best_hit = hit_list[i];
      max_hits = hit_list[i]->max_offset_hits;
    }
  }

  if (best_hit != NULL) {
    TRY(poly = lfsrdesc_get_poly(best_hit->desc));
    printf(
        "\033[1mBEST MATCH: [%s] WITH %d/%d HITS\033[0m\n",
        poly,
        best_hit->hits,
        files);
    free(poly);

    for (i = 0; i < best_hit->params_hit_count; ++i)
      if (best_hit->params_hit_list[i]->hits == max_hits)
        best_offset = best_hit->params_hit_list[i]->offset;

    printf(
        "\033[1mBEST OFFSET: %d WITH %d/%d HITS\033[0m\n",
        best_offset,
        max_hits,
        best_hit->hits);

    files = 0;

    for (i = 1; i < argc; ++i)
      if (lfsr_hit_descramble_file(
        best_hit,
        argv[i],
        i,
        best_offset))
        ++files;

    printf(
        "\033[1mDESCRAMBLED %d FILES UNDER %s\033[0m\n",
        files,
        OUTPUT_DIRECTORY);
  } else {
    printf("%s: no candidate polynomials found. Shame :(\n", argv[0]);
  }

  return 0;

fail:
  exit(EXIT_FAILURE);
}

