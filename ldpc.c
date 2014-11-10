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

#include "ldpc.h"
#include <stdlib.h>
#include <string.h>

void ldpc_ll_matrix_destroy(ldpc_ll_matrix_t *H) {
    int i;
    ldpc_ll_edge_t *cur;
    ldpc_ll_edge_t *next;
    if (!H)
        return;

    for(i=0;i<H->M;i++)
    {
        cur = H->rows[i];
        while(cur)
        {
            next = cur->right;
            free(cur);
            cur = next;
        }
    }

    free(H->rows);
    free(H->cols);

    free(H);
    return;
}

void ldpc_param_init(ldpc_param_t *param)
{
    param->h_matrix = NULL;
    param->max_iter = 30; /* Default to 30 iterations */

    return;
}

void ldpc_param_destroy(ldpc_param_t *param)
{
    ldpc_ll_matrix_destroy(param->h_matrix);
}


/* For now, encodes a K*128 byte array (one bit per byte)
 * Warning! This function is not created for speed, but only for testing.
 *
 */
int ldpc_encode_c(ldpc_param_t *p, int len, char *input, char *output)
{
    int bit;
    ldpc_ll_edge_t *node;

    if (p->h_matrix && len == p->h_matrix->K)
    {
        /* The first K bits in the encoded codeword are identical to the input */
        memcpy(output, input, p->h_matrix->K*sizeof(char));

        /* Zero out the rest of the output data */
        memset( (void *)(output+p->h_matrix->K), 0, p->h_matrix->M*sizeof(char) );

        /* Calculate the parity bits */
        for (bit=0;bit<p->h_matrix->M;bit++)
        {
            output[p->h_matrix->K + bit] = 0;
            node = p->h_matrix->rows[bit];
            do {
                output[p->h_matrix->K + bit] ^= output[node->col];
                node = node->right;
            } while(node->right);
        }

    }
    else
        return -1;

    return 1;
}

int (*ldpc_encode)(ldpc_param_t *p, int len, char *input, char *output) = ldpc_encode_c;


