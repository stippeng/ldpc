/*****************************************************************

    Copyright (C) 2014 Stefan Grönroos

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

#ifndef LDPC_H
#define LDPC_H

#include <stdio.h>

typedef struct ldpc_t ldpc_t; //Decoder handle

/* Intermediate format for representing the LDPC H matrix, using linked lists */
typedef struct ldpc_ll_edge_t {
    int idx, row, col;
    struct ldpc_ll_edge_t *right;
    struct ldpc_ll_edge_t *down;
} ldpc_ll_edge_t;

typedef struct {
    int M,N,K, num_edges;
    ldpc_ll_edge_t **rows;
    ldpc_ll_edge_t **cols;
} ldpc_ll_matrix_t;

/* Decoder initialization parameters */
typedef struct ldpc_param_t {
    /* To initialize the decoder with a certain code,
     * supply an ldpc_ll_matrix_t containing the H matrix
     */
    ldpc_ll_matrix_t *h_matrix;

    unsigned short max_iter; /* Maximum number of LDPC decoder iterations */
    unsigned short num_threads; /* Number of simulataneous active worker threads */

} ldpc_param_t;

/********************
 * Decoder functions
*********************/

/*
 * Initialize the decoder resources using a populated ldpc_param_t structure.
 * This needs to be done before any call to ldpc_decode.
 */
ldpc_t* (*ldpc_init)(ldpc_param_t *param);

/*
 * Decode one batch of codewords.
 * As of now, this takes a batch of 128 encoded codewords (soft-bits, 8-bit value per bit), and produces
 * 128 decoded codewords in bitval (memory for bitval must be allocated).
 */
int (*ldpc_decode)(ldpc_t *h, char *llr_in, unsigned char *bitval);

/*
 * Free decoder resources
 */
void (*ldpc_destroy)(ldpc_t *h);

/* Give the required size of the input LLR array required by the decoder */
size_t (*ldpc_decoder_input_size)(ldpc_t *h);
/* Give the required size of the output buffer */
size_t (*ldpc_decoder_output_size)(ldpc_t *h);

/*******************
 * Encoder functions
 *******************/
int (*ldpc_encode)(ldpc_param_t *h, int len, char *input, char *output);


/*******************
 * Misc functions
 *******************/
void ldpc_ll_matrix_destroy(ldpc_ll_matrix_t *H);
void ldpc_param_init(ldpc_param_t *param);
void ldpc_param_destroy(ldpc_param_t *param);

#endif //LDPC_H
