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

//---------------------------------------------------------------------
// !!! Be aware, the 2D-FFT algorithm still needs some improvements !!!
//---------------------------------------------------------------------

#ifndef _LFFT_FFT2_H
#define	_LFFT_FFT2_H

#ifdef	__cplusplus
extern "C" {
#endif

#include "lfft_config.h"
#include "lfft_fft.h"

typedef struct _lfft_Fft2
{
    uint16_t rows; //!< number of rows
    uint16_t columns; //!< number of columns

    int32_t ** result_real_temp; //!< real result
    int32_t ** result_imag_temp; //!< imaginary result
    int32_t ** result_real; //!< temporary real data
    int32_t ** result_imag; //!< temporary imaginary data

    lfft_Fft * fft_rows; //!< fft to calculate the fft in a row
    lfft_Fft * fft_columns; //!< fft to calculate the fft in a column

} lfft_Fft2;

/*!
 * Initilizes the 2D-fft.
 * \param fft2 pointer to struct to be initialized
 * \param rows number of rows of the data
 * \param columns number of columns of the data
 * \return 0: successful 1: rows or columns is not to the power of 2
 */
lfft_errno lfft_fft2_new(lfft_Fft2 * fft2, uint16_t rows, uint16_t columns);

/*!
 * Deallocate the used memory.
 * \param fft2 initialized lfft_Fft2 struct
 */
void lfft_fft2_delete(lfft_Fft2 * fft2);

/*!
 * Calculates the 2D-fft with the fft algorithm of core/lfft_fft.h.
 * It uses only real input data, the imaginary part is set to 0.
 * \param fft2 initialized lfft_Fft2 struct
 * \param real 2D-real input data for fft.
 */
void lfft_fft2(lfft_Fft2 * fft2, int32_t ** real);

/*!
 * Calculates the 2D-fft with the fft algorithm of core/lfft_fft.h.
 * It uses real input and imaginary input data.
 * \param fft2 initialized lfft_Fft2 struct
 * \param real 2D-real input data for fft.
 * \param imag 2D-imaginary input data for fft.
 */
void lfft_fft2_complex(lfft_Fft2 * fft2, int32_t ** real, int32_t ** imag);

/*!
 * Calculates the 2D-fft with the fft algorithm of core/lfft_fft.h.
 * It uses only real input data, the imaginary part is set to 0.
 * \param fft2 initialized lfft_Fft struct
 * \param real 2D-real input data for fft.
 */
void lfft_fft2_float(lfft_Fft2 * fft2, float ** real);

/*!
 * Calculates the 2D-fft with the fft algorithm of core/lfft_fft.h.
 * It uses real input and imaginary input data.
 * \param fft2 initialized lfft_Fft2 struct
 * \param real 2D-real input data for fft.
 * \param imag 2D-imaginary input data for fft.
 */
void lfft_fft2_complex_float(lfft_Fft2 * fft2, float ** real, float ** imag);

/*!
 * Calculates the 2D-inverse-fft with the fft algorithm of core/lfft_fft.h.
 * It uses only real input data, the imaginary part is set to 0.
 * \param fft2 initialized lfft_Fft2 struct
 * \param real 2D-real input data for fft.
 */
void lfft_ifft2(lfft_Fft2 * fft2, int32_t ** real);

/*!
 * Calculates the 2D-inverse-fft with the fft algorithm of core/lfft_fft.h.
 * It uses real input and imaginary input data.
 * \param fft2 initialized lfft_Fft2 struct
 * \param real 2D-real input data for fft.
 * \param imag 2D-imaginary input data for fft.
 */
void lfft_ifft2_complex(lfft_Fft2 * fft2, int32_t ** real, int32_t ** imag);

/*!
 * Calculates the 2D-inverse-fft with the fft algorithm of core/lfft_fft.h.
 * It uses only real input data, the imaginary part is set to 0.
 * \param fft2 initialized lfft_Fft struct
 * \param real 2D-real input data for fft.
 */
void lfft_ifft2_float(lfft_Fft2 * fft2, float ** real);

/*!
 * Calculates the 2D-inverse-fft with the fft algorithm of core/lfft_fft.h.
 * It uses real input and imaginary input data.
 * \param fft2 initialized lfft_Fft2 struct
 * \param real 2D-real input data for fft.
 * \param imag 2D-imaginary input data for fft.
 */
void lfft_ifft2_complex_float(lfft_Fft2 * fft2, float ** real, float ** imag);

/*!
 * Returns an element of the real part of the 2D-fft
 * \param fft2 initialized lfft_Fft2 struct
 * \param row row of the requested element
 * \param column column of the requested element
 * \return real result of the 2D-fft
 */
int32_t lfft_fft2_result_real_at(lfft_Fft2 * fft2, uint16_t row, uint16_t column);

/*!
 * Returns an element of the imaginary part of the 2D-fft
 * \param fft2 initialized lfft_Fft2 struct
 * \param row row of the requested element
 * \param column column of the requested element
 * \return imaginary result of the 2D-fft
 */
int32_t lfft_fft2_result_imag_at(lfft_Fft2 * fft2, uint16_t row, uint16_t column);

/*!
 * Returns an element of the real part of the 2D-fft
 * \param fft2 initialized lfft_Fft2 struct
 * \param row row of the requested element
 * \param column column of the requested element
 * \return real result of the 2D-fft
 */
float lfft_fft2_result_real_float_at(lfft_Fft2 * fft2, uint16_t row, uint16_t column);

/*!
 * Returns an element of the imaginary part of the 2D-fft
 * \param fft2 initialized lfft_Fft2 struct
 * \param row row of the requested element
 * \param column column of the requested element
 * \return imaginary result of the 2D-fft
 */
float lfft_fft2_result_imag_float_at(lfft_Fft2 * fft2, uint16_t row, uint16_t column);

/*!
 * Calculates the base algorithm of the 2D-fft.
 * You have to reorder and adjust the input data manually.
 * Don't use this function if it is possible to use the above functions.
 * \param fft2 initialized lfft_Fft2 struct
 * \param calculate_ifft2 false to calculate 2D-fft, true to calculate 2D-ifft
 */
void _lfft_fft2_calculation(lfft_Fft2 * fft2, bool calculate_ifft2);

#ifdef	__cplusplus
}
#endif

#endif	/* _LFFT_FFT2_H */

