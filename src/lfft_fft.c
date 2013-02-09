/* LambdaFFT is covered by the new BSD license
 *
 * Copyright (c) 2011-2013, Clemens Korner
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *  1. Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *  2. Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

/*! \file
 * This file contains functions to calculate the fft (fast fourier
 * transformation) and ifft (inverse fast fourier transfomration)
 * It uses the radix-2 algorithm for calculation.
 * The calculation is optimized for processors without a FPU.
 * float is only used during the initialization process.
 *
 * \author Clemens Korner
 * \version 0.1.0
 */

#include "lfft_fft.h"

#include "lfft_config.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/*!
 * Initializes the fft and allocates the necessary memory.
 * \param fft pointer to struct to be initialized
 * \param samples number of input samples; must be to the power of 2
 * \return 0: successful 1: samples is not to the power of 2
 */
static lfft_errno _lfft_fft_init(lfft_Fft * fft, uint16_t samples);

/*!
 * Prepares the input data for the calculation.
 * It uses only real input data, the imaginary part is set to 0.
 * \param fft initialized lfft_Fft struct
 * \param real real input data for fft
 * \param calculate_ifft true: calculate ifft false: calculate fft
 */
static void _lfft_fft(lfft_Fft * fft, const int32_t real[], bool calculate_ifft);

/*!
 * Prepares the input data for the calculation.
 * It uses real input and imaginary input data.
 * \param fft initialized lfft_Fft struct
 * \param real real input data for fft
 * \param imag imaginary input data for fft.
 * \param calculate_ifft true: calculate ifft false: calculate fft
 */
static void _lfft_fft_complex(lfft_Fft * fft, const int32_t real[], const int32_t imag[], bool calculate_ifft);

/*!
 * Prepares the input data for the calculation.
 * It uses only real input data, the imaginary part is set to 0.
 * \param fft initialized lfft_Fft struct
 * \param real real input data for fft
 * \param calculate_ifft true: calculate ifft false: calculate fft
 */
static void _lfft_fft_float(lfft_Fft * fft, const float real[], bool calculate_ifft);

/*!
 * Prepares the input data for the calculation.
 * It uses real input and imaginary input data.
 * \param fft initialized lfft_Fft struct
 * \param real real input data for fft
 * \param imag imaginary input data for fft.
 * \param calculate_ifft true: calculate ifft false: calculate fft
 */
static void _lfft_fft_complex_float(lfft_Fft * fft, const float real[], const float imag[], bool calculate_ifft);

/*!
 * Checks if x is to the power of 2.
 * \param x number to be checkeds
 * \return true if x is to the power of 2, false if not
 */
static bool _lfft_is_power_2(uint16_t x);

lfft_errno lfft_fft_new(lfft_Fft * fft, uint16_t samples)
{
    return _lfft_fft_init(fft, samples);
}

void lfft_fft_delete(lfft_Fft * fft)
{
    // deallocate memory
    free(fft->switching_table);
    free(fft->wk_real);
    free(fft->wk_imag);
    free(fft->result_real);
    free(fft->result_imag);
}

void lfft_fft(lfft_Fft * fft, const int32_t real[])
{
    _lfft_fft(fft, real, false);
}

void lfft_fft_complex(lfft_Fft * fft, const int32_t real[], const int32_t imag[])
{
    _lfft_fft_complex(fft, real, imag, false);
}

void lfft_fft_float(lfft_Fft * fft, const float real[])
{
    _lfft_fft_float(fft, real, false);
}

void lfft_fft_complex_float(lfft_Fft * fft, const float real[], const float imag[])
{
    _lfft_fft_complex_float(fft, real, imag, false);

}

void lfft_ifft(lfft_Fft * fft, const int32_t real[])
{
    _lfft_fft(fft, real, true);
}

void lfft_ifft_complex(lfft_Fft * fft, const int32_t real[], const int32_t imag[])
{
    _lfft_fft_complex(fft, real, imag, true);
}

void lfft_ifft_float(lfft_Fft * fft, const float real[])
{
    _lfft_fft_float(fft, real, true);
}

void lfft_ifft_complex_float(lfft_Fft * fft, const float real[], const float imag[])
{
    _lfft_fft_complex_float(fft, real, imag, true);

}

int32_t lfft_fft_result_real_at(lfft_Fft * fft, uint16_t n)
{
    return fft->result_real[n]>>LFFT_RHS_BITS;
}

int32_t lfft_fft_result_imag_at(lfft_Fft * fft, uint16_t n)
{
    return fft->result_imag[n]>>LFFT_RHS_BITS;
}

float lfft_fft_result_real_float_at(lfft_Fft * fft, uint16_t n)
{
    return ((float) fft->result_real[n])/(1<<LFFT_RHS_BITS);
}

float lfft_fft_result_imag_float_at(lfft_Fft * fft, uint16_t n)
{
    return ((float) fft->result_imag[n])/(1<<LFFT_RHS_BITS);
}

uint16_t lfft_fft_abs_at(lfft_Fft * fft, uint16_t n)
{
    uint32_t op1 = fft->result_real[n]>>(LFFT_RHS_BITS);
    uint32_t op2 = fft->result_imag[n]>>(LFFT_RHS_BITS);
    op1 *= op1;
    op2 *= op2;
    return lfft_isqrt(op1+op2);
}

uint16_t lfft_fft_abs_and_norm_at(lfft_Fft * fft, uint16_t n)
{
    uint32_t op1 = fft->result_real[n]>>(fft->steps+LFFT_RHS_BITS);
    uint32_t op2 = fft->result_imag[n]>>(fft->steps+LFFT_RHS_BITS);
    op1 *= op1;
    op2 *= op2;
    return lfft_isqrt(op1+op2);
}

uint16_t lfft_isqrt(uint32_t x)
{
    // this function is an adapted version of an isqrt algorithm
    // written by Wilco Dijkstra
    uint32_t a;
    uint32_t root = 0;

#define ISQRT_ITER(N) \
    a = root+(1<<(N)); \
    if(x >= a<<(N)) \
    { \
        x -= a<<(N); \
        root |= 2<<(N); \
    }

    ISQRT_ITER(15);
    ISQRT_ITER(14);
    ISQRT_ITER(13);
    ISQRT_ITER(12);
    ISQRT_ITER(11);
    ISQRT_ITER(10);
    ISQRT_ITER(9);
    ISQRT_ITER(8);
    ISQRT_ITER(7);
    ISQRT_ITER(6);
    ISQRT_ITER(5);
    ISQRT_ITER(4);
    ISQRT_ITER(3);
    ISQRT_ITER(2);
    ISQRT_ITER(1);
    ISQRT_ITER(0);

#undef ISQRT_ITER

    return root>>1;
}

void _lfft_fft_calculation(lfft_Fft * fft, bool calculate_ifft)
{
    uint16_t i;
    uint16_t j;
    uint16_t n_wk;
    uint16_t butterfly_counter;

    // temporary variables
    int32_t real_1;
    int32_t real_2;
    int32_t imag_1;
    int32_t imag_2;
    int32_t wk_real;
    int32_t wk_imag;

    uint16_t space_butterfly_operant = 1; // space between operants of butterfly graph
    uint16_t n_wk_counter            = fft->samples/2; // counter to calculate next n_wk

    for(i = 0; i < fft->steps; i++)
    {
        j = 0;
        n_wk = 0;
        butterfly_counter = 0;
        while(j < fft->samples)
        {
            // bitand with ones_mask to delete unwanted carry
            n_wk = fft->ones_mask&n_wk;

            real_1 = fft->result_real[j];
            imag_1 = fft->result_imag[j];
            imag_2 = fft->result_imag[j+space_butterfly_operant];
            real_2 = fft->result_real[j+space_butterfly_operant];

            // in both cases, fft and ifft, wk_real = fft->wk_real[n_wk]
            // cos(x) == cos(-x)
            wk_real = fft->wk_real[n_wk];
            // -sin(x) == sin(-x)
            if(calculate_ifft)
            {
                wk_imag = -fft->wk_imag[n_wk];

            }
            else
            {
                wk_imag = fft->wk_imag[n_wk];
            }

            // first butterfly operant
            fft->result_real[j] = real_1+
                    ((real_2*wk_real)>>LFFT_RHS_BITS)-
                    ((imag_2*wk_imag)>>LFFT_RHS_BITS);

            fft->result_imag[j] = imag_1+
                    ((imag_2*wk_real)>>LFFT_RHS_BITS)+
                    ((real_2*wk_imag)>>LFFT_RHS_BITS);

            // second butterfly operant
            fft->result_real[j+space_butterfly_operant] = real_1-
                    ((real_2*wk_real)>>LFFT_RHS_BITS)+
                    ((imag_2*wk_imag)>>LFFT_RHS_BITS);

            fft->result_imag[j+space_butterfly_operant] = imag_1-
                    ((imag_2*wk_real)>>LFFT_RHS_BITS)-
                    ((real_2*wk_imag)>>LFFT_RHS_BITS);

            n_wk += n_wk_counter;

            butterfly_counter++;
            j++;
            if(butterfly_counter == space_butterfly_operant)
            {
                butterfly_counter = 0;
                j += space_butterfly_operant;
            }
        }
        // n_wk_counter = n_wk_counter/2
        n_wk_counter >>= 1;
        // space_butterfly_operant = space_butterfly_operant*2
        space_butterfly_operant <<= 1;
    }

    if(calculate_ifft)
    {
        // divide the results by the fft size (fft->samples)
        // x/fft->samples == x>>fft->steps
        for(i = 0; i < fft->samples; i++)
        {
            fft->result_real[i] = fft->result_real[i]>>fft->steps;
            fft->result_imag[i] = fft->result_imag[i]>>fft->steps;
        }
    }
}

static lfft_errno _lfft_fft_init(lfft_Fft * fft, uint16_t samples)
{
    uint16_t i;
    uint16_t j;
    uint16_t new_place;

    // samples is not to the power of 2
    if(!_lfft_is_power_2(samples))
    {
        return 1;
    }

    fft->samples = samples;

    // log2(x) = log10(x)/log10(2)
    // make a round before casting to uint8_t, otherwise there are sometimes
    // rounding errors
    // round(x) = floor(x+0.5)
    fft->steps =  ((uint8_t) floor((log((float) fft->samples)/log(2.0f))+0.5f));
    // create binary mask
    // example: fft->samples = 8 --> fft->ones_mask = 0b11 */
    fft->ones_mask = (fft->samples>>1)-1;

    // allocate memory
    fft->switching_table = (uint16_t *) calloc(fft->samples, sizeof(uint16_t));
    fft->wk_real         = (int32_t *)  calloc(fft->samples/2, sizeof(int32_t));
    fft->wk_imag         = (int32_t *)  calloc(fft->samples/2, sizeof(int32_t));
    fft->result_real     = (int32_t *)  calloc(fft->samples, sizeof(int32_t));
    fft->result_imag     = (int32_t *)  calloc(fft->samples, sizeof(int32_t));

    // create array with ascending numbers from 0 to fft->samples-1
    for(i = 0; i < fft->samples; i++)
    {
        fft->switching_table[i] = i;
    }

    // calculate the switching table; making the content of the table bit-reverse
    for(i = 0; i < fft->samples; i++)
    {
        new_place = 0;
        for(j = 0; j < fft->steps; j++)
        {
            new_place <<= 1;
            new_place += fft->switching_table[i]&1;
            fft->switching_table[i] >>= 1;
        }
        fft->switching_table[i] = new_place;
    }

    // calculate different wk_n
    // the result is multiplied with 2^LFFT_RHS_BITS, to have the size some
    // decimal place bits in the least sigificant bits of the integer values
    // wk = e^(-j*2*pi*i/fft->samples)*(1<<LFFT_RHS_BITS)
    for(i = 0; i < fft->samples/2; i++)
    {
        fft->wk_real[i] = (int32_t) (cos((2.0f*M_PI*i)/fft->samples)*(1<<LFFT_RHS_BITS));
        fft->wk_imag[i] = (int32_t) ((-sin((2.0f*M_PI*i)/fft->samples))*(1<<LFFT_RHS_BITS));
    }

    return 0;
}

static void _lfft_fft(lfft_Fft * fft, const int32_t real[], bool calculate_ifft)
{
    uint16_t i;

    // reorder real input data and multiply with 2^LFFT_RHS_BITS
    for(i = 0; i < fft->samples; i++)
    {
        fft->result_real[i] = real[fft->switching_table[i]]<<LFFT_RHS_BITS;
    }

    // imaginary part is set to 0
    memset(fft->result_imag, 0, sizeof(fft->result_imag));

    _lfft_fft_calculation(fft, calculate_ifft);
}

static void _lfft_fft_complex(lfft_Fft * fft, const int32_t real[], const int32_t imag[], bool calculate_ifft)
{
    uint16_t i;

    // reorder real and imaginary input data and multiply with 2^LFFT_RHS_BITS
    for(i = 0; i < fft->samples; i++)
    {
        fft->result_real[i] = real[fft->switching_table[i]]<<LFFT_RHS_BITS;
        fft->result_imag[i] = imag[fft->switching_table[i]]<<LFFT_RHS_BITS;
    }

    _lfft_fft_calculation(fft, calculate_ifft);
}

static void _lfft_fft_float(lfft_Fft * fft, const float real[], bool calculate_ifft)
{
    uint16_t i;

    // reorder real input data and multiply with 2^LFFT_RHS_BITS
    for(i = 0; i < fft->samples; i++)
    {
        fft->result_real[i] = (int32_t) (real[fft->switching_table[i]]*(1<<LFFT_RHS_BITS));
    }

    // imaginary part is set to 0
    memset(fft->result_imag, 0, sizeof(fft->result_imag));

    _lfft_fft_calculation(fft, calculate_ifft);
}

static void _lfft_fft_complex_float(lfft_Fft * fft, const float real[], const float imag[], bool calculate_ifft)
{
    uint16_t i;

    // reorder real and imaginary input data and multiply with 2^LFFT_RHS_BITS
    for(i = 0; i < fft->samples; i++)
    {
        fft->result_real[i] = (int32_t) (real[fft->switching_table[i]]*(1<<LFFT_RHS_BITS));
        fft->result_imag[i] = (int32_t) (imag[fft->switching_table[i]]*(1<<LFFT_RHS_BITS));
    }

    _lfft_fft_calculation(fft, calculate_ifft);
}

static bool _lfft_is_power_2(uint16_t x)
{
    return x&&(!(x&(x-1)));
}
