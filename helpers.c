/*****************************************************************
    Some helper macros for benchmarking and CUDA functions

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

#include <stdio.h>
#include <stdlib.h>

// For benchmarking
#ifdef BENCHMARKING
struct timespec time1[5], time2;
struct timespec diff(struct timespec start, struct timespec end)
{
        struct timespec temp;
        if ((end.tv_nsec-start.tv_nsec)<0) {
                temp.tv_sec = end.tv_sec-start.tv_sec-1;
                temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
        } else {
                temp.tv_sec = end.tv_sec-start.tv_sec;
                temp.tv_nsec = end.tv_nsec-start.tv_nsec;
        }
        return temp;
}
#endif

/* Convert a character array to a one byte per bit format */
char *bits_to_bytes(int len, char *in)
{
    int i;
    unsigned char c;
    char *out = (char *)malloc(len*8*sizeof(char));

    for (i=0;i<len;i++)
    {
        for (c=0;c<8;c++)
        {
            out[i*8+c] = (in[i] >> c) & 1;
        }
    }
}
