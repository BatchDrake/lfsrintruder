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

#ifndef _TYPES_H
#define _TYPES_H

#include <util.h>
#include <stdint.h>

#ifndef _RELEASE
#  include <stdio.h>
#  define _DEBUG(fmt, arg...)                  \
  fprintf(                                      \
      stderr,                                   \
      "[i] " fmt,                               \
      ##arg)
#else
#  define _DEBUG(fmt, arg...)
#endif /* !_RELEASE */

#define TRY_EXCEPT(expr, catch_stmts)        \
  if (!(expr)) {                                \
    catch_stmts;                                \
    goto fail;                                  \
  }

#define TRY(expr)                            \
    TRY_EXCEPT(                              \
        expr,                                   \
        _DEBUG(                                \
            "%s:%d: exception in \"%s\"\n",     \
            __FILE__,                           \
            __LINE__,                           \
            STRINGIFY(expr)                     \
            ))

#define ALLOCATE_FFT(dest, n)         \
    TRY_EXCEPT(                              \
        dest = fftwf_alloc_complex(n),         \
        _DEBUG(                                \
            "%s:%d: failed to allocate FFT array of %d elements\n", \
            __FILE__,                           \
            __LINE__,                           \
            n))


#define ALLOCATE_MANY(dest, n, type)         \
    TRY_EXCEPT(                              \
        dest = calloc(n, sizeof(type)),         \
        _DEBUG(                                \
            "%s:%d: failed to allocate %d objects of type %s\n", \
            __FILE__,                           \
            __LINE__,                           \
            n,                                  \
            STRINGIFY(type)                     \
            ))

#define ALLOCATE(dest, type) ALLOCATE_MANY(dest, 1, type)

#define INIT(self, type, params...)          \
  TRY_EXCEPT(                                \
      JOIN(type, _init) (self, ##params),       \
      _DEBUG(                                  \
          "%s:%d: failed to call constructor of class \"%s\"\n", \
          __FILE__,                             \
          __LINE__,                             \
          STRINGIFY(params)))                   \


#define CONSTRUCT(dest, type, params...)     \
  TRY_EXCEPT(                                \
      dest = JOIN(type, _new) (params),       \
      _DEBUG(                                  \
          "%s:%d: failed to construct object of class %s\n", \
          __FILE__,                             \
          __LINE__,                             \
          STRINGIFY(type)))                   \

enum boolean {
  FALSE,
  TRUE
};

typedef enum boolean BOOL;

static inline unsigned int
popcount64(uint64_t b)
{
  b = (b & 0x5555555555555555ull) + (b >> 1 & 0x5555555555555555ull);
  b = (b & 0x3333333333333333ull) + (b >> 2 & 0x3333333333333333ull);
  b = b + (b >> 4) & 0x0F0F0F0F0F0F0F0Full;
  b = b + (b >> 8);
  b = b + (b >> 16);
  b = b + (b >> 32) & 0x0000007F;

  return (unsigned int) b;
}


#endif /* _TYPES_H */
