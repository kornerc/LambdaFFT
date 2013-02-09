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
 * This file contains settings of the fft.
 *
 * \author Clemens Korner
 * \version 0.1.0
 */

#ifndef _LFFT_CONFIG_H
#define _LFFT_CONFIG_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
// if your system doesn't support stdint.h,
// try to comment #include <stdint.h> and uncomment the following lines
// (check if the sizes of the datatypes match)
    //typedef signed char int8_t; // 8 Bit signed
    //typedef unsigned char uint8_t; // 8 Bit unsigned
    //typedef short int16_t; // 16 Bit signed
    //typedef unsigned short uint16_t; // 16 Bit unsigned
    //typedef int int32_t; // 32 Bit signed
    //typedef unsigned int uint32_t; // 32 Bit unsigned

#include <stdbool.h>
// if your system doesn't support stdbool.h,
// try to comment #include <stdbool.h> and uncomment the following lines
    //typedef int8_t _Bool;
    //#define bool _Bool
    //#define true 1
    //#define false 0

// use constants like M_PI under Windows
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif /* _USE_MATH_DEFINES */


/*!
 * bits for decimal places
 */
#define LFFT_RHS_BITS 8

#ifdef __cplusplus
}
#endif

#endif /* _LFFT_CONFIG_H */
