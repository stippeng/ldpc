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

#include "ldpc_sse.h"

struct ldpc_t {
    int M;
    int N;
    int K;

    ldpc_edge_t *hb;
    ldpc_edge_t *hc;
    ldpc_msg_t *edge_msg;

    int *row_idx;
    int *col_idx;

    int *llr_map;

    int num_edges;

    struct bn_update_args *bn_args;
    struct bn_update_bitval_args *bn_bv_args;
    struct cn_update_args *cn_args;
    struct check_satisfied_args *cs_args;

    char *unsatisfied_eqs;

    unsigned short max_iter;
    unsigned short num_threads;

    pthread_barrier_t *barr_0;
    pthread_barrier_t *barr_1;

};


void ldpc_init_messages_sse(ldpc_t *h)
{
    /* This seems to be faster than memset */
    __m128i zero = _mm_setzero_si128();
    ldpc_llr_t *init_p = (ldpc_llr_t *)h->edge_msg;

    for (ldpc_llr_t *p = init_p; p < init_p+(h->num_edges*LDPC_CODEWORD_BLOCKS); p++) {
        _mm_store_si128((__m128i *)p, zero);
    }

}


/* Bit node update without hard decision */
void *sse_ldpc_ms_bn_update(void *threadarg) {

    __m128i msg;
    __m128i m;
    __m128i xmmtmp;
    const __m128i mesfloor = _mm_set1_epi8(-127);
    int col_start;

    int index;
    int temp;
    struct bn_update_args *arg;
    arg = (struct bn_update_args *) threadarg;
    unsigned short iter = 0;

    xmmtmp = _mm_set1_epi8(1);

    while (iter++ < arg->max_iter)
    {
        for (int i=arg->first_n; i < arg->first_n + arg->num_n; i++) {
            for (int cw_block = 0; cw_block < LDPC_CODEWORD_BLOCKS;cw_block++) {

                col_start = arg->col_idx[i];
                index = col_start;
                temp = (i*LDPC_CODEWORD_BLOCKS + cw_block);

                m = _mm_load_si128((__m128i *)&arg->llr[temp].v);

                do {
                    msg = _mm_load_si128((__m128i *)&arg->emsg[index].message[cw_block].v);

                    m = _mm_adds_epi8(m, msg);

                    index = arg->hc[index].i_next;
                } while (index != col_start);

                do {
                    msg = _mm_load_si128((__m128i *)&arg->emsg[index].message[cw_block].v);

                    msg = _mm_subs_epi8(m, msg);
                    //Hack: We do not want -128, as that ruins correction performance
                    msg = _mm_adds_epi8(msg, xmmtmp);
                    //msg = _mm_min_epu8(msg,mesfloor);

                    _mm_store_si128((__m128i *)&arg->emsg[index].message[cw_block].v, msg);
                    index = arg->hc[index].i_next;
                } while (index != col_start);

            }
        }
        pthread_barrier_wait(arg->barr_0);
        pthread_barrier_wait(arg->barr_1);
    }

    pthread_exit(NULL);
}

/* Bit node update with hard decision */
void *sse_ldpc_ms_bn_update_bitval(void *threadarg) {

    __m128i msg;
    __m128i m;
    const __m128i xmmtmp = _mm_set1_epi8(1);
    const __m128i mesfloor = _mm_set1_epi8(-127);
    int col_start;

    int index;
    int temp;
    struct bn_update_bitval_args *arg;
    arg = (struct bn_update_bitval_args *) threadarg;

    for (int i=arg->first_n; i < arg->first_n + arg->num_n; i++) {
        for (int cw_block = 0; cw_block < LDPC_CODEWORD_BLOCKS;cw_block++) {

            col_start = arg->col_idx[i];
            index = col_start;
            temp = (i*LDPC_CODEWORD_BLOCKS + cw_block);

            m = _mm_load_si128((__m128i *)&arg->llr[temp].v);

            do {
                msg = _mm_load_si128((__m128i *)&arg->emsg[index].message[cw_block].v);

                m = _mm_adds_epi8(m, msg);

                index = arg->hc[index].i_next;
            } while (index != col_start);

            do {
                msg = _mm_load_si128((__m128i *)&arg->emsg[index].message[cw_block].v);

                msg = _mm_subs_epi8(m, msg);
                //Hack: We do not want -128, as that ruins correction performance
                msg = _mm_adds_epi8(msg, xmmtmp);
                //msg = _mm_min_epu8(msg,mesfloor);

                _mm_store_si128((__m128i *)&arg->emsg[index].message[cw_block].v, msg);
                index = arg->hc[index].i_next;
            } while (index != col_start);

            // Hard decision
            // bitval should be sign bit of m
            msg = _mm_srli_epi16(m, 7);
            msg = _mm_and_si128(msg, xmmtmp);
            _mm_store_si128((__m128i *)&arg->bitval[temp].v, msg);

        }
    }

    pthread_exit(NULL);
}

void *sse_ldpc_ms_cn_update(void *threadarg) {
    __m128i minLLR, nMinLLR, absol, minMsg,tmp1, tmp2, mask1, mask2, mask3, msg, sign, zero;

    int row_start;
    char counter;
    unsigned short iter = 0;

    struct cn_update_args *arg;
    arg = (struct cn_update_args *)threadarg;

    while (iter++ < arg->max_iter)
    {
        pthread_barrier_wait(arg->barr_0);
        for (int i=arg->first_n; i < arg->first_n + arg->num_n; i++) {
            for (int cw_block = 0; cw_block < LDPC_CODEWORD_BLOCKS;cw_block++) {
                minLLR = _mm_set1_epi8(127);
                nMinLLR = _mm_set1_epi8(127);
                sign =  _mm_set1_epi8(1);
                counter = 0;

                row_start = arg->row_idx[i];
                int index = row_start;

                do {
                    //msg = emsg[index].message[cw_block].b[cw];
                    msg = _mm_load_si128((__m128i *)&arg->emsg[index].message[cw_block].v);

                    //sign *= (msg >= 0 ? 1:-1);
                    sign = _mm_xor_si128(sign, msg);

                    //abs(msg)
                    absol = _mm_abs_epi8(msg);


                    mask1 = _mm_cmplt_epi8(absol, minLLR); //0xFF if less than minLLR
                    mask2 = _mm_cmplt_epi8(absol, nMinLLR); //0xFF if less than nMinLLR

                    tmp1 = VEC_BLEND(absol, minLLR, mask1);
                    nMinLLR = VEC_BLEND(nMinLLR, tmp1, mask2);

                    minLLR = VEC_BLEND(minLLR, absol, mask1);

                    tmp1 = _mm_set1_epi8(counter);
                    minMsg = VEC_BLEND(minMsg, tmp1, mask1);

                    index = arg->hb[index].i_next;
                    counter++;
                } while (index != row_start);

                counter = 0;
                zero = _mm_setzero_si128();
                sign = _mm_cmplt_epi8(sign,zero); //0xFF = -1 if < 0

                do {
                    msg = _mm_load_si128((__m128i *)&arg->emsg[index].message[cw_block].v);

                    mask1 = _mm_cmpeq_epi8(minMsg, _mm_set1_epi8(counter));
                    mask2 = _mm_cmplt_epi8(msg, zero);
                    mask3 = _mm_xor_si128(sign, mask2); //if sign*msg < 0 =>0xFF, else 0x00
                    mask3 = _mm_or_si128(mask3, _mm_set1_epi8(1));
                    tmp1 = _mm_sign_epi8(minLLR, mask3);
                    tmp2 = _mm_sign_epi8(nMinLLR, mask3);

                    msg = VEC_BLEND(tmp1, tmp2, mask1);

                    _mm_store_si128((__m128i *)&arg->emsg[index].message[cw_block].v, msg);
                    index = arg->hb[index].i_next;
                    counter++;
                } while (index != row_start);
            }
        }
        pthread_barrier_wait(arg->barr_1);
    }

    pthread_exit(NULL);

}

void *sse_ldpc_ms_check_unsatisfied(void *threadarg) {
    __m128i hard_val;
    __m128i sum;
    const __m128i zeros = _mm_setzero_si128();

    int row_start;
    char unsatisfied = 0;

    struct check_satisfied_args *arg;
    arg = (struct check_satisfied_args *)threadarg;

    for (int i=arg->first_n; i < arg->first_n + arg->num_n; i++) {
        row_start = arg->row_idx[i];
        for (int cw_block = 0; cw_block < LDPC_CODEWORD_BLOCKS;cw_block++) {
            sum = _mm_setzero_si128();

            int index = row_start;

            do {
                hard_val = _mm_load_si128((__m128i *) &arg->bitval[arg->llr_map[index]*LDPC_CODEWORD_BLOCKS+cw_block].v);
                sum = _mm_xor_si128(sum, hard_val);

                index = arg->hb[index].i_next;
            } while (index != row_start);
            if (! _mm_movemask_epi8(_mm_cmpeq_epi32(sum, zeros))) {
                unsatisfied = 1;
                goto unsat;
            }
        }
    }
    unsat:
    arg->unsat_eq[0] = unsatisfied;
    pthread_exit(NULL);

}

ldpc_t *ldpc_init_sse(ldpc_param_t *param)
{
    ldpc_t *h;
    ldpc_ll_edge_t *node;
    int i, first_idx;

    /* Set up LDPC code structures from supplied matrix */
    if (param->h_matrix)
    {
        h = (ldpc_t *)calloc(1, sizeof(ldpc_t));
        h->M = param->h_matrix->M;
        h->N = param->h_matrix->N;
        h->K = param->h_matrix->K;
        h->num_edges = param->h_matrix->num_edges;

        h->hb = (ldpc_edge_t *)_mm_malloc(h->num_edges * sizeof(ldpc_edge_t), 16);
        h->hc = (ldpc_edge_t *)_mm_malloc(h->num_edges * sizeof(ldpc_edge_t), 16);
        h->row_idx = (int *)_mm_malloc(h->M * sizeof(int), 16);
        h->col_idx = (int *)_mm_malloc(h->N * sizeof(int), 16);
        h->llr_map = (int *)_mm_malloc(h->num_edges * sizeof(int), 16);

        for(i=0;i<h->M;i++)
        {
            node = param->h_matrix->rows[i];
            first_idx = node->idx;
            h->row_idx[i] = node->idx;
            do {
                if (node->right)
                    h->hb[node->idx].i_next = node->right->idx;
                else /* No further edges in row. Link back to first edge */
                    h->hb[node->idx].i_next = first_idx;

                node = node->right;
            } while (node);
        }

        for(i=0;i<h->N;i++)
        {
            node = param->h_matrix->cols[i];
            first_idx = node->idx;
            h->col_idx[i] = node->idx;
            do {
                if (node->down)
                    h->hc[node->idx].i_next = node->down->idx;
                else
                    h->hc[node->idx].i_next = first_idx;

                h->llr_map[node->idx] = i;

                node = node->down;
            } while (node);
        }

    } else {
        fprintf(stderr, "No LDPC code supplied!\n");
        return NULL;
    }

    /* num_threads should divide both N and M */
    h->num_threads = CLAMP(param->num_threads, 1, h->M); /* num_threads must be between 1 and M */
    while (h->M % h->num_threads != 0 || h->N % h->num_threads != 0)
        --h->num_threads;

    fprintf(stderr, "Using %d simultaneous threads\n", h->num_threads);

    /* Allocate decoder memory */
    h->edge_msg = (ldpc_msg_t *)_mm_malloc(h->num_edges * sizeof(ldpc_msg_t), 64);

    h->bn_args = (struct bn_update_args *)malloc(h->num_threads*sizeof(struct bn_update_args));
    h->bn_bv_args = (struct bn_update_bitval_args *)malloc(h->num_threads*sizeof(struct bn_update_bitval_args));
    h->cn_args = (struct cn_update_args *)malloc(h->num_threads*sizeof(struct cn_update_args));
    h->cs_args = (struct check_satisfied_args *)malloc(h->num_threads*sizeof(struct check_satisfied_args));

    h->unsatisfied_eqs = (char *)_mm_malloc(h->num_threads * sizeof(char), 16);
    h->max_iter = CLAMP(param->max_iter, 0, 100);

    h->barr_0 = malloc(sizeof(pthread_barrier_t));
    h->barr_1 = malloc(sizeof(pthread_barrier_t));
    pthread_barrier_init(h->barr_0, NULL, h->num_threads*2);
    pthread_barrier_init(h->barr_1, NULL, h->num_threads*2);


    for(int i=0; i<h->num_threads; i++) {
        int first, num;
        num = h->N/h->num_threads;
        first = i*num;
        if (i == h->num_threads - 1) num = h->N-first;

        h->bn_args[i].first_n = first;
        h->bn_args[i].num_n = num;
        h->bn_args[i].hc = h->hc;
        h->bn_args[i].N = h->N;
        h->bn_args[i].emsg = h->edge_msg;
        h->bn_args[i].col_idx = h->col_idx;
        h->bn_args[i].max_iter = h->max_iter;
        h->bn_args[i].barr_0 = h->barr_0;
        h->bn_args[i].barr_1 = h->barr_1;

        h->bn_bv_args[i].first_n = first;
        h->bn_bv_args[i].num_n = num;
        h->bn_bv_args[i].hc = h->hc;
        h->bn_bv_args[i].N = h->N;
        h->bn_bv_args[i].emsg = h->edge_msg;
        h->bn_bv_args[i].col_idx = h->col_idx;

        num = h->M/h->num_threads;
        first = i*num;
        if (i == h->num_threads - 1) num = h->M-first;

        h->cn_args[i].first_n = first;
        h->cn_args[i].num_n = num;
        h->cn_args[i].hb = h->hb;
        h->cn_args[i].M = h->M;
        h->cn_args[i].emsg = h->edge_msg;
        h->cn_args[i].row_idx = h->row_idx;
        h->cn_args[i].max_iter = h->max_iter;
        h->cn_args[i].barr_0 = h->barr_0;
        h->cn_args[i].barr_1 = h->barr_1;

        h->cs_args[i].first_n = first;
        h->cs_args[i].num_n = num;
        h->cs_args[i].M = h->M;
        h->cs_args[i].llr_map = h->llr_map;
        h->cs_args[i].row_idx = h->row_idx;
        h->cs_args[i].unsat_eq = &h->unsatisfied_eqs[i];
        h->cs_args[i].hb = h->hb;
    }

    return h;
}

void ldpc_destroy_sse(ldpc_t *h) {
   _mm_free(h->hb);
   _mm_free(h->hc);
   _mm_free(h->edge_msg);
   _mm_free(h->row_idx);
   _mm_free(h->col_idx);
   _mm_free(h->llr_map);

   /* Free thread data */
   free(h->bn_args);
   free(h->bn_bv_args);
   free(h->cn_args);
   free(h->cs_args);

   free(h->barr_0);
   free(h->barr_1);

   free(h);
}

int ldpc_decode_sse(ldpc_t *h, char *llr, unsigned char *bitval) {
    ldpc_llr_t *llr_interl;
    ldpc_bit_t *bitval_interl;

    START_CLOCK(1);


    llr_interl = (ldpc_llr_t *) _mm_malloc(LDPC_CODEWORD_BLOCKS*h->N*sizeof(ldpc_llr_t), 64); //64-byte to align with cache lines

    for(int i = 0; i < h->N;i++) {
        for (int n=0;n<LDPC_CODEWORD_BLOCKS;n++)  {
            for (int k=0;k<16;k++) {
                llr_interl[i*LDPC_CODEWORD_BLOCKS+n].b[k] = llr[(16*n + k)*h->N + i];
            }
        }
    }


    bitval_interl = (ldpc_bit_t *)_mm_malloc(LDPC_CODEWORD_BLOCKS*h->N*sizeof(ldpc_bit_t), 64);

    /* Set up threads */
    for(int i=0; i<h->num_threads; i++) {
        h->bn_args[i].llr = llr_interl;

        h->bn_bv_args[i].llr = llr_interl;
        h->bn_bv_args[i].bitval = bitval_interl;

        h->cs_args[i].bitval = bitval_interl;



    }
    pthread_t thread_handles[h->num_threads*2];
    int rc;
    void *thread_status;

    BENCHMARK_NOW(1, "SSE LDPC intitialization time: ");

    START_CLOCK(1);
    ldpc_init_messages_sse(h);

    int unsat = 1;

    /* Bit node update
           Note: No early termination
        */
    for (int t=0; t<h->num_threads; t++) {
        rc = pthread_create(&thread_handles[t], NULL, sse_ldpc_ms_bn_update, (void *) &h->bn_args[t]);
        if (rc) {
            fprintf(stderr, "ERROR; return code from pthread_create() is %d\n", rc);
            exit(-1);
        }
    }

    /* Check node update */
    for (int t=0; t<h->num_threads; t++) {
        rc = pthread_create(&thread_handles[h->num_threads+t], NULL, sse_ldpc_ms_cn_update, (void *) &h->cn_args[t]);
        if (rc) {
            fprintf(stderr, "ERROR; return code from pthread_create() is %d\n", rc);
            exit(-1);
        }
    }

    // Join all the threads in order to ensure they're all finished
    for (int t=0; t<h->num_threads*2; t++) {
        rc = pthread_join(thread_handles[t], &thread_status);
        if (rc) {
            printf("ERROR; return code from pthread_join() is %d\n", rc);
            exit(-1);
        }
    }

    if(unsat != 0) {
        /* Bit node update with hard decision */
        for (int t=0; t<h->num_threads; t++) {
            rc = pthread_create(&thread_handles[t], NULL, sse_ldpc_ms_bn_update_bitval, (void *) &h->bn_bv_args[t]);
            if (rc) {
                fprintf(stderr, "ERROR; return code from pthread_create() is %d\n", rc);
                exit(-1);
            }
        }
        for (int t=0; t<h->num_threads; t++) {
            rc = pthread_join(thread_handles[t], &thread_status);
            if (rc) {
                printf("ERROR; return code from pthread_join() is %d\n", rc);
                exit(-1);
            }
        }
    }


    BENCHMARK_NOW(1, "SSE LDPC decode time: ");

    START_CLOCK(1)

    for(int n=0;n<h->K;n++) {
        for(int i=0;i<LDPC_CODEWORD_BLOCKS;i++) {
            for (int k=0;k<16;k++) {
                bitval[(i*16+k)*h->K + n] = (unsigned char) bitval_interl[n*LDPC_CODEWORD_BLOCKS + i].b[k];
            }
        }
    }

    BENCHMARK_NOW(1, "SSE LDPC copy back time: ");

    //Free memory
    _mm_free(llr_interl);
    _mm_free(bitval_interl);

    return 1;

}

size_t ldpc_decoder_input_size_sse(ldpc_t *h) {
    return h->N*LDPC_CODEWORD_BLOCKS*16*sizeof(char);
}

size_t ldpc_decoder_output_size_sse(ldpc_t *h) {
    return h->K*LDPC_CODEWORD_BLOCKS*16*sizeof(unsigned char);
}

/* Link architecture-dependent functions to function pointers defined in interface */
int (*ldpc_decode)(ldpc_t *h, char *llr_in, unsigned char *bitval) = ldpc_decode_sse;
ldpc_t * (*ldpc_init)(ldpc_param_t *param) = ldpc_init_sse;
void (*ldpc_destroy)(ldpc_t *h) = ldpc_destroy_sse;
size_t (*ldpc_decoder_input_size)(ldpc_t *h) = ldpc_decoder_input_size_sse;
size_t (*ldpc_decoder_output_size)(ldpc_t *h) = ldpc_decoder_output_size_sse;

