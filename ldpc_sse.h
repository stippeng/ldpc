/*****************************************************************

    Copyright (C) 2013 Stefan Grönroos

    Authors: Stefan Grönroos <stefan.gronroos@abo.fi>

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

********************************************************************/

#ifndef SSE_LDPC_H
#define SSE_LDPC_H

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <pthread.h>
#include <time.h>
#include "ldpc.h"
#include "helpers.h"

//We require at least SSSE3
#ifdef __SSE4__
#include <smmintrin.h>
#else
#include <tmmintrin.h>
#endif

/* If we process 16 codewords in parallel using SIMD instructions,
    we need 8 blocks of codewords to decode 128 codewords. */
#define LDPC_CODEWORD_BLOCKS 8

//Number of threads to be spawned for each checknode/bitnode update in each iteration.
#define LDPC_MAX_NUM_THREADS 128

typedef char i8_vec __attribute__ ((__vector_size__ (16)));
typedef unsigned char u8_vec __attribute__ ((__vector_size__ (16)));

typedef union {
    i8_vec v;
    char b[16];
} ldpc_llr_t;

typedef union {
    u8_vec v;
    unsigned char b[16];
} ldpc_bit_t;

typedef struct {
        int i_next;
} ldpc_edge_t;

typedef struct {
    ldpc_llr_t message[LDPC_CODEWORD_BLOCKS];
} ldpc_msg_t;

struct bn_update_args {
    int first_n;
    int num_n;
    ldpc_edge_t *hc;
    ldpc_llr_t *llr;
    int N;
    ldpc_msg_t *emsg;
    int *col_idx;
    unsigned short max_iter;
    pthread_barrier_t *barr_0;
    pthread_barrier_t *barr_1;
};

struct bn_update_bitval_args {
    int first_n;
    int num_n;
    ldpc_edge_t *hc;
    ldpc_bit_t *bitval;
    ldpc_llr_t *llr;
    int N;
    ldpc_msg_t *emsg;
    int *col_idx;
};

struct cn_update_args {
    int first_n;
    int num_n;
    ldpc_edge_t *hb;
    int M;
    ldpc_msg_t *emsg;
    int *row_idx;
    unsigned short max_iter;
    pthread_barrier_t *barr_0;
    pthread_barrier_t *barr_1;
};

struct check_satisfied_args {
    int first_n;
    int num_n;
    int M;
    ldpc_bit_t *bitval;
    ldpc_edge_t *hb;
    int *llr_map;
    int *row_idx;
    char *unsat_eq;
};


//#define IMAX(X,Y) (max(X,Y))
//#define IMIN(X,Y) (min(X,Y))
#define IMAX(X,Y) ((X >= Y) ? X : Y)
#define IMIN(X,Y) ((X < Y) ? X : Y)

#define CLAMP(VAL,MINV,MAXV) (IMIN(IMAX(VAL,MINV),MAXV))

#define ABS(X) (X >= 0 ? X : -X)

#define CHARTOFLOAT_H(X) (((float)(X)) / F2C_FACTOR)
#define CHARTOFLOAT(X) (char2float[X+127])

#ifdef __SSE4__
#define VEC_BLEND(a, b, mask) (_mm_blendv_epi8(a, b, mask))
#else
#define VEC_BLEND(a, b, mask) (_mm_or_si128(_mm_and_si128(mask, b), _mm_andnot_si128(mask, a)))
#endif


/* Initialize memory structures */
//void sse_ldpc_ms_load(ldpc_edge_t *hb_host, ldpc_edge_t *hc_host, int numEdges, int M, int N, int* gpu_row_idx, int *gpu_col_idx, int *gpu_llr_map);
ldpc_t *ldpc_init_sse(ldpc_param_t *param);

/* Decode a block of LLRs */
int ldpc_decode_sse(ldpc_t *h, char *llr, unsigned char *bitval);

/* Clean up memory */
void ldpc_destroy_sse(ldpc_t *h);

#endif // SSE_LDPC_H
