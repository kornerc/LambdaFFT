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
 * This file contains functions to calculate the 2D-fft (2 dimensional fast fourier
 * transformation) and 2D-ifft (2 dimwnsional inverse fast fourier transfomration)
 * It uses the fft algorithm of core/lfft_fft.h to calcolate the result.
 *
 * \author Clemens Korner
 * \version 0.1.0
 */

#include "lfft_fft2.h"

#include <stdlib.h>
#include <string.h>

/*!
 * Prepares the input data for the calculation.
 * It uses only real input data, the imaginary part is set to 0.
 * \param fft2 initialized lfft_Fft2 struct
 * \param real 2D-real input data for fft
 * \param calculate_ifft2 true: calculate 2D-ifft false: calculate 2D-fft
 */
static void _lfft_fft2(lfft_Fft2 * fft2, int32_t ** real, bool calculate_ifft2);

/*!
 * Prepares the input data for the calculation.
 * It uses real input and imaginary input data.
 * \param fft2 initialized lfft_Fft2 struct
 * \param 2D-real real input data for fft
 * \param 2D-imag imaginary input data for fft.
 * \param calculate_ifft2 true: calculate 2D-ifft false: calculate 2D-fft
 */
static void _lfft_fft2_complex(lfft_Fft2 * fft2, int32_t ** real, int32_t ** imag, bool calculate_ifft2);

/*!
 * Prepares the input data for the calculation.
 * It uses only real input data, the imaginary part is set to 0.
 * \param fft2 initialized lfft_Fft2 struct
 * \param real 2D-real input data for fft
 * \param calculate_ifft2 true: calculate 2D-ifft false: calculate 2D-fft
 */
static void _lfft_fft2_float(lfft_Fft2 * fft2, float ** real, bool calculate_ifft2);

/*!
 * Prepares the input data for the calculation.
 * It uses real input and imaginary input data.
 * \param fft2 initialized lfft_Fft2 struct
 * \param real 2D-real input data for fft
 * \param imag 2D-imaginary input data for fft.
 * \param calculate_ifft2 true: calculate 2D-ifft false: calculate 2D-fft
 */
static void _lfft_fft2_complex_float(lfft_Fft2 * fft2, float ** real, float ** imag, bool calculate_ifft2);

lfft_errno lfft_fft2_new(lfft_Fft2 * fft2, uint16_t rows, uint16_t columns)
{
    uint16_t i;

    fft2->rows = rows;
    fft2->columns = columns;

    fft2->fft_rows    = (lfft_Fft *) malloc(sizeof(lfft_Fft));
    fft2->fft_columns = (lfft_Fft *) malloc(sizeof(lfft_Fft));

    // initialize fft's and check if they size is to the power of 2
    if(lfft_fft_new(fft2->fft_rows, fft2->columns)||lfft_fft_new(fft2->fft_columns, fft2->rows))
    {
        return 1;
    }

    fft2->result_real_temp = (int32_t **) malloc(fft2->columns*sizeof(int32_t *));
    fft2->result_imag_temp = (int32_t **) malloc(fft2->columns*sizeof(int32_t *));

    for(i = 0; i < fft2->columns; i++)
    {
        fft2->result_real_temp[i] = (int32_t *) malloc(fft2->rows*sizeof(int32_t));
        fft2->result_imag_temp[i] = (int32_t *) malloc(fft2->rows*sizeof(int32_t));
    }

    fft2->result_real = (int32_t **) malloc(fft2->rows*sizeof(int32_t *));
    fft2->result_imag = (int32_t **) malloc(fft2->rows*sizeof(int32_t *));

    for(i = 0; i < fft2->rows; i++)
    {
        fft2->result_real[i] = (int32_t *) malloc(fft2->columns*sizeof(int32_t));
        fft2->result_imag[i] = (int32_t *) malloc(fft2->columns*sizeof(int32_t));
    }

    return 0;
}

void lfft_fft2_delete(lfft_Fft2 * fft2)
{
    lfft_fft_delete(fft2->fft_rows);
    lfft_fft_delete(fft2->fft_columns);

    free(fft2->result_real);
    free(fft2->result_imag);
    free(fft2->result_real_temp);
    free(fft2->result_imag_temp);
}

void lfft_fft2(lfft_Fft2 * fft2, int32_t ** real)
{
    _lfft_fft2(fft2, real, false);
}

void lfft_fft2_complex(lfft_Fft2 * fft2, int32_t ** real, int32_t ** imag)
{
    _lfft_fft2_complex(fft2, real, imag, false);
}

void lfft_fft2_float(lfft_Fft2 * fft2, float ** real)
{
    _lfft_fft2_float(fft2, real, false);
}

void lfft_fft2_complex_float(lfft_Fft2 * fft2, float ** real, float ** imag)
{
    _lfft_fft2_complex_float(fft2, real, imag, false);
}

void lfft_ifft2(lfft_Fft2 * fft2, int32_t ** real)
{
    _lfft_fft2(fft2, real, true);
}

void lfft_ifft2_complex(lfft_Fft2 * fft2, int32_t ** real, int32_t ** imag)
{
    _lfft_fft2_complex(fft2, real, imag, true);
}

void lfft_ifft2_float(lfft_Fft2 * fft2, float ** real)
{
    _lfft_fft2_float(fft2, real, true);
}

void lfft_ifft2_complex_float(lfft_Fft2 * fft2, float ** real, float ** imag)
{
    _lfft_fft2_complex_float(fft2, real, imag, true);
}

int32_t lfft_fft2_result_real_at(lfft_Fft2 * fft2, uint16_t row, uint16_t column)
{
    return fft2->result_real[row][column]>>LFFT_RHS_BITS;
}

int32_t lfft_fft2_result_imag_at(lfft_Fft2 * fft2, uint16_t row, uint16_t column)
{
    return fft2->result_imag[row][column]>>LFFT_RHS_BITS;
}

float lfft_fft2_result_real_float_at(lfft_Fft2 * fft2, uint16_t row, uint16_t column)
{
    return ((float) fft2->result_real[row][column])/(1<<LFFT_RHS_BITS);
}

float lfft_fft2_result_imag_float_at(lfft_Fft2 * fft2, uint16_t row, uint16_t column)
{
    return ((float) fft2->result_imag[row][column])/(1<<LFFT_RHS_BITS);
}

void _lfft_fft2_calculation(lfft_Fft2 * fft2, bool calculate_ifft2)
{
    uint16_t i;
    uint16_t j;

    for(i = 0; i < fft2->rows; i++)
    {
        fft2->fft_rows->result_real = fft2->result_real[i];
        fft2->fft_rows->result_imag = fft2->result_imag[i];
        _lfft_fft_calculation(fft2->fft_rows, calculate_ifft2);
        // rotate the result and switch the values
        //   input data     rotated data   rotated and switched data
        //   |00 01 02 03|  |00 10 20 30|  |00 10 20 30|
        //   |10 11 12 13|  |01 11 21 31|  |02 12 22 32|
        //   |20 21 22 23|  |02 12 22 32|  |01 11 21 31|
        //   |30 31 32 33|  |03 13 23 33|  |03 13 23 33|
        for(j = 0; j < fft2->columns; j++)
        {
            fft2->result_real_temp[j][fft2->fft_columns->switching_table[i]] = fft2->fft_rows->result_real[j];
            fft2->result_imag_temp[j][fft2->fft_columns->switching_table[i]] = fft2->fft_rows->result_imag[j];
        }
    }

    for(i = 0; i < fft2->columns; i++)
    {
        fft2->fft_columns->result_real = fft2->result_real_temp[i];
        fft2->fft_columns->result_imag = fft2->result_imag_temp[i];
        _lfft_fft_calculation(fft2->fft_columns, calculate_ifft2);
        // rotate the result
        //   data temp      rotated data
        //   |00 01 02 03|  |00 10 20 30|
        //   |10 11 12 13|  |01 11 21 31|
        //   |20 21 22 23|  |02 12 22 32|
        //   |30 31 32 33|  |03 13 23 33|
        for(j = 0; j < fft2->rows; j++)
        {
            fft2->result_real[j][i] = fft2->fft_columns->result_real[j];
            fft2->result_imag[j][i] = fft2->fft_columns->result_imag[j];
        }
    }
}

static void _lfft_fft2(lfft_Fft2 * fft2, int32_t ** real, bool calculate_ifft2)
{
    uint16_t i;
    uint16_t j;

    for(i = 0; i < fft2->rows; i++)
    {
        // imaginary part is set to 0
        memset(fft2->result_imag[i], 0, sizeof(fft2->result_imag[0][0])*fft2->columns);
        for(j = 0; j < fft2->columns; j++)
        {
            // copy and reorder real part
            fft2->result_real[i][j] = real[i][fft2->fft_rows->switching_table[j]]<<LFFT_RHS_BITS;
        }
    }

    _lfft_fft2_calculation(fft2, calculate_ifft2);
}

static void _lfft_fft2_complex(lfft_Fft2 * fft2, int32_t ** real, int32_t ** imag, bool calculate_ifft2)
{
    uint16_t i;
    uint16_t j;

    for(i = 0; i < fft2->rows; i++)
    {
        for(j = 0; j < fft2->columns; j++)
        {
            // copy and reorder real part
            fft2->result_real[i][j] = real[i][fft2->fft_rows->switching_table[j]]<<LFFT_RHS_BITS;
            // copy and reorder imag part
            fft2->result_imag[i][j] = imag[i][fft2->fft_rows->switching_table[j]]<<LFFT_RHS_BITS;
        }
    }

    _lfft_fft2_calculation(fft2, calculate_ifft2);
}

static void _lfft_fft2_float(lfft_Fft2 * fft2, float ** real, bool calculate_ifft2)
{
    uint16_t i;
    uint16_t j;

    for(i = 0; i < fft2->rows; i++)
    {
        // imaginary part is set to 0
        memset(fft2->result_imag[i], 0, sizeof(fft2->result_imag[0][0])*fft2->columns);
        for(j = 0; j < fft2->columns; j++)
        {
            // copy and reorder real part
            fft2->result_real[i][j] = (int32_t) (real[i][fft2->fft_rows->switching_table[j]]*(1<<LFFT_RHS_BITS));
        }
    }

    _lfft_fft2_calculation(fft2, calculate_ifft2);
}

static void _lfft_fft2_complex_float(lfft_Fft2 * fft2, float ** real, float ** imag, bool calculate_ifft2)
{
    uint16_t i;
    uint16_t j;

    for(i = 0; i < fft2->rows; i++)
    {
        for(j = 0; j < fft2->columns; j++)
        {
            // copy and reorder real part
            fft2->result_real[i][j] = (int32_t) (real[i][fft2->fft_rows->switching_table[j]]*(1<<LFFT_RHS_BITS));
            // copy and reorder imag part
            fft2->result_imag[i][j] = (int32_t) (imag[i][fft2->fft_rows->switching_table[j]]*(1<<LFFT_RHS_BITS));
        }
    }

    _lfft_fft2_calculation(fft2, calculate_ifft2);
}

