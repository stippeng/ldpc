/*****************************************************************
    Functions for parsing alist files describing an LDPC code.

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
#include <stdio.h>
#include <stdlib.h>

void ldpc_ll_add_edge(ldpc_ll_matrix_t *H, int row, int col, int idx)
{
    ldpc_ll_edge_t *curr;
    ldpc_ll_edge_t *node;

    node = (ldpc_ll_edge_t *)calloc(1, sizeof(ldpc_ll_edge_t));
    node->down = NULL;
    node->right = NULL;
    node->idx = idx;
    node->row = row;
    node->col = col;

    curr = H->rows[row];
    if (!curr)
        H->rows[row] = node;
    else
    {
        while (curr->right)
            curr = curr->right;
        curr->right = node;
    }

    curr = H->cols[col];
    if(!curr)
        H->cols[col] = node;
    else
    {
        while (curr->down)
            curr = curr->down;
        curr->down = node;
    }


}

/* Parse file "fname" and produce an ldpc_ll_matrix_t containing the LDPC code */
ldpc_ll_matrix_t *ldpc_alist_parse(char *fname) {
    FILE *fp = NULL;
    ldpc_ll_matrix_t *H = NULL;
    int N,M, max_cd, max_bd;
    int *cdegs = NULL;
    int *bdegs = NULL;
    int *row = NULL;
    int i,j,idx;
    int succ = 0;

    if (fp = fopen(fname, "r")) {
        if (fscanf(fp, "%d %d %d %d ", &M, &N, &max_cd, &max_bd ) == 4) {
            fprintf(stderr, "N: %d, M: %d, max_cd: %d, max_bd: %d\n", N, M, max_cd, max_bd);
            cdegs = (int *)calloc(M, sizeof(int));
            bdegs = (int *)calloc(N, sizeof(int));
            H = (ldpc_ll_matrix_t *)calloc(1, sizeof(ldpc_ll_matrix_t));
            H->rows = (ldpc_ll_edge_t **)calloc(M, sizeof(ldpc_ll_edge_t*));
            H->cols = (ldpc_ll_edge_t **)calloc(N, sizeof(ldpc_ll_edge_t*));
            H->N = N;
            H->M = M;
            H->K = N-M;


            row = (int *)calloc(max_cd, sizeof(int));

            for (i=0;i<M;i++)
                if (!fscanf(fp, "%d ", &cdegs[i]))
                    goto cleanup;
            for (i=0;i<N;i++)
                if (!fscanf(fp, "%d ", &bdegs[i]))
                    goto cleanup;

            idx = 0;
            for (i=0;i<M;i++)
            {
                for (j=0;j<max_cd;j++)
                    if (!fscanf(fp, "%d ", &row[j]))
                        goto cleanup;

                for (j=0;j<cdegs[i];j++)
                {
                    ldpc_ll_add_edge(H, i, row[j]-1, idx++);
                }
            }

            H->num_edges = idx;
            fprintf(stderr, "edges: %d\n", idx);
            succ = 1;
        } else {
            fprintf(stderr, "Alist file %s: bad format\n", fname);
            goto cleanup;
        }
    cleanup:
        if (fp) fclose(fp);
        if (cdegs) free(cdegs);
        if (bdegs) free(bdegs);
        if (row) free(row);
        if (!succ && H) ldpc_ll_matrix_destroy(H);

        return succ ? H : NULL;

    } else {
        fprintf(stderr, "Error opening alist file %s\n", fname);
        return NULL;
    }
}
