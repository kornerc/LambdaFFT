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

#ifndef _LFFT_FFT_H
#define _LFFT_FFT_H

#ifdef __cplusplus
extern "C" {
#endif

#include "lfft_config.h"

typedef int8_t lfft_errno;

typedef struct _lfft_Fft
{
    uint16_t   samples; //!< number of samples
    uint8_t    steps; //!< number of steps
    uint16_t   ones_mask; //!< binary mask used for fft calculation

    uint16_t * switching_table; //!< table containing switching values for the input data
    int32_t  * wk_real; //!< real wk values
    int32_t  * wk_imag; //!< imaginary wk values

    int32_t  * result_real; //!< array where real part of the fft calcilation is saved
    int32_t  * result_imag; //!< array where imaginary part of the fft calcilation is saved
} lfft_Fft;

/*!
 * Initilizes the fft.
 * \param fft pointer to struct to be initialized
 * \param samples number of input samples; must be to the power of 2
 * \return 0: successful 1: samples is not to the power of 2
 */
lfft_errno lfft_fft_new(lfft_Fft * fft, uint16_t samples);

/*!
 * Deallocate the used memory.
 * \param fft initialized lfft_Fft struct
 */
void lfft_fft_delete(lfft_Fft * fft);

/*!
 * Calculates the fft with the radix-2 fft algorithm.
 * It uses only real input data, the imaginary part is set to 0.
 * \param fft initialized lfft_Fft struct
 * \param real real input data for fft.
 */
void lfft_fft(lfft_Fft * fft, const int32_t real[]);

/*!
 * Calculates the fft with the radix-2 fft algorithm.
 * It uses real input and imaginary input data.
 * \param fft initialized lfft_Fft struct
 * \param real real input data for fft.
 * \param imag imaginary input data for fft.
 */
void lfft_fft_complex(lfft_Fft * fft, const int32_t real[], const int32_t imag[]);

/*!
 * Calculates the fft with the radix-2 fft algorithm.
 * It uses only real input data, the imaginary part is set to 0.
 * \param fft initialized lfft_Fft struct
 * \param real real input data for fft.
 */
void lfft_fft_float(lfft_Fft * fft, const float real[]);

/*!
 * Calculates the fft with the radix-2 fft algorithm.
 * It uses real input and imaginary input data.
 * \param fft initialized lfft_Fft struct
 * \param real real input data for fft.
 * \param imag imaginary input data for fft.
 */
void lfft_fft_complex_float(lfft_Fft * fft, const float real[], const float imag[]);

/*!
 * Calculates the inverse-fft with the radix-2 fft algorithm.
 * It uses only real input data, the imaginary part is set to 0.
 * \param fft initialized lfft_Fft struct
 * \param real real input data for fft.
 */
void lfft_ifft(lfft_Fft * fft, const int32_t real[]);

/*!
 * Calculates the inverse-fft with the radix-2 fft algorithm.
 * It uses real input and imaginary input data.
 * \param fft initialized lfft_Fft struct
 * \param real real input data for fft.
 * \param imag imaginary input data for fft.
 */
void lfft_ifft_complex(lfft_Fft * fft, const int32_t real[], const int32_t imag[]);

/*!
 * Calculates the inverse-fft with the radix-2 fft algorithm.
 * It uses only real input data, the imaginary part is set to 0.
 * \param fft initialized lfft_Fft struct
 * \param real real input data for fft.
 */
void lfft_ifft_float(lfft_Fft * fft, const float real[]);

/*!
 * Calculates the inverse-fft with the radix-2 fft algorithm.
 * It uses real input and imaginary input data.
 * \param fft initialized lfft_Fft struct
 * \param real real input data for fft.
 * \param imag imaginary input data for fft.
 */
void lfft_ifft_complex_float(lfft_Fft * fft, const float real[], const float imag[]);

/*!
 * Returns an element of the real part of the fft
 * \param fft initialized lfft_Fft struct
 * \param n number of the requested element
 * \return real result of the fft
 */
int32_t lfft_fft_result_real_at(lfft_Fft * fft, uint16_t n);

/*!
 * Returns an element of the imag part of the fft
 * \param fft initialized lfft_Fft struct
 * \param  n number of the requested element
 * \return imaginary result of the fft
 */
int32_t lfft_fft_result_imag_at(lfft_Fft * fft, uint16_t n);

/*!
 * Returns an element of the real part of the fft
 * \param fft initialized lfft_Fft struct
 * \param n number of the requested element
 * \return real result of the fft
 */
float lfft_fft_result_real_float_at(lfft_Fft * fft, uint16_t n);

/*!
 * Returns an element of the imag part of the fft
 * \param fft initialized lfft_Fft struct
 * \param  n number of the requested element
 * \return imaginary result of the fft
 */
float lfft_fft_result_imag_float_at(lfft_Fft * fft, uint16_t n);

/*!
 * Calculates the absoulte value at n.
 * result = sqrt(real[n]^2+imag[n]^2)
 * \param fft initialized lfft_Fft struct
 * \param n element at n
 * \return absolute value at n
 */
uint16_t lfft_fft_abs_at(lfft_Fft * fft, uint16_t n);

/*!
 * Calculates the absoulte value and normalizes it at n.
 * result = sqrt((real[n]/samples)^2+(imag[n]/samples)^2)
 * \param fft initialized lfft_Fft struct
 * \param n element at n
 * \return absolute and normalized value at n
 */
uint16_t lfft_fft_abs_and_norm_at(lfft_Fft * fft, uint16_t n);

/*!
 * Calculates the integer square root.
 * \param x positive number
 * \return the integer square root of x
 */
uint16_t lfft_isqrt(uint32_t x);

/*!
 * Calculates the base algorithm of the fft.
 * You have to reorder and adjust the input data manually.
 * Don't use this function if it is possible to use the above functions.
 * \param fft initialized lfft_Fft struct
 * \param calculate_ifft false to calculate fft, true to calculate ifft
 */
void _lfft_fft_calculation(lfft_Fft * fft, bool calculate_ifft);


#ifdef __cplusplus
}
#endif

#endif /* _LFFT_FFT_H */
