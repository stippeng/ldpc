/*****************************************************************
    Simple test program for generating a random bit sequence, encoding it,
    decoding it, and checking the resulting BER.

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

#include "ldpc.h" /* LDPC decoder interface */
#include "alist.h"
#include "ldpc.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#define ROUNDS 5

int main() {
    ldpc_ll_matrix_t *H;
    ldpc_param_t param;
    ldpc_t *decoder;
    char *input;
    char *enc;
    char *chan;
    unsigned char *dec;
    int r,cw;
    float ber;
    int sum, csum = 0;
    int rounds = 0;

    srand(time(NULL));

    /***
     * First, we load an LDPC code, i.e. its associated parity-check matrix
     * into a linked-list-based format from a .alist file.
     * The parity check matrix determines the number of data and parity bits in a
     * coded bit sequence (a codeword).
     * N (the number of columns in the parity check matrix) determines the total
     * number of bits in a codeword. M (number of rows in the matrix) determines
     * the number of parity (redundancy) bits. N-M = K is the number of data (input)
     * bits.
     ***/
    H = ldpc_alist_parse("matrices/dvbs2_LDPC_matrix_12.alist");
    if (!H)
        return 1;

    /***
     * An ldpc_param_t structure contains various parameters used to set up the
     * decoder instance. The most important parameter is a parity check matrix
     * in the ldpc_ll_matrix_t format, as produced by for example the ldpc_alist_parse
     * function.
     ***/
    ldpc_param_init(&param);
    param.h_matrix = H;
    param.max_iter = 30; //Number of decoder iterations to run
    param.num_threads = 16;

    /* Now we can initialize the decoder using the populated param structure. */
    decoder = ldpc_init(&param);

    /***
     * As of now, the decoder decodes 128 codewords at a time (H->K*128 bits),
     * so all data needs to be allocated as multiples of 128.
     *
     * Here we allocate memory for test input, encoded data, and decoded data.
     * The input data has length 128 x H->K bits in a one bit per byte format.
     * This is encoded into 128 x H->N (H->N > H->K!) including parity bits.
     *
     * In chan, we convert the encoded data to "soft" bits, where a negative value means
     * the bit is probably a one, while a positive probably means a zero.
     *
     * The decoder gives us the decoded 128 x H->K bits back.
     ***/

    input = (char *)malloc(128*H->K*sizeof(char));
    enc = (char *)malloc(128*H->N*sizeof(char));
    chan = (char *)malloc(128*H->N*sizeof(char));
    dec = (unsigned char *)malloc(128*H->K*sizeof(unsigned char));

    while (rounds++ < ROUNDS)
    {
        fprintf(stderr, "Round: %d\n", rounds);
        /* Generate a random bit sequence as input */
        for(r=0;r<128*H->K;r++)
            input[r] = rand()%2;

        /* The encoder takes one block of H->K bits and outputs H->N encoded bits */
        for (r=0;r<128;r++)
            ldpc_encode(&param, H->K, input+(r*H->K), enc+(r*H->N));

        /* Very simple "channel model" with soft bit generation :-) */
        for (r=0;r<128*H->N;r++)
        {
            chan[r] = enc[r] ? (-64+(rand()%76)) : 64-(rand()%76);
        }

        /***
         * Decode a 128 x H->N input sequence, into a 128 x H->K output bit sequence.
         * The decoder currently does NOT check wether the output bits are valid!
         ***/
        for (r=0;r<1;r++)
            ldpc_decode(decoder, chan, dec);

        /***
         * Calculate the bit-error-ratio (BER) by comparing to the
         * original input data.
         ***/
        sum = 0;
        for (cw=0;cw<128;cw++)
            for (r=0;r<H->K;r++)
                sum += input[cw*H->K + r] == dec[cw*H->K + r] ? 1 : 0;

        ber = 1-(((float)sum)/(float)(128*H->K));
        printf("BER: %f (%d)\n", ber, 128*H->K-sum);
        csum += sum;
    }

    ber = 1-(((float)csum)/(float)(128*H->K*ROUNDS));
    fprintf(stderr, "Total BER: %f (%d)\n", ber, 128*H->K*ROUNDS-csum);

    /* Free decoder resources, param (including H matrix), and other
     * allocated memory */

    ldpc_destroy(decoder);
    ldpc_param_destroy(&param);
    free(input);
    free(enc);
    free(dec);
    free(chan);


    return 0;
}
