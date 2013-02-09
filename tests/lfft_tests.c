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

#include <inttypes.h>
#include <math.h>
#include <stdlib.h>

#include <check.h>

#include "lfft.h"
#include "core/lfft_fft.c"

START_TEST(test_f__lfft_is_power_2)
{
    uint16_t i;

    fail_unless(_lfft_is_power_2(0x0000) == false, "False assumption: _lfft_is_power_2(0x0000) == false\n", i);
    fail_unless(_lfft_is_power_2(0x0003) == false, "False assumption: _lfft_is_power_2(0x0003) == false\n", i);
    fail_unless(_lfft_is_power_2(0x0005) == false, "False assumption: _lfft_is_power_2(0x0005) == false\n", i);
    fail_unless(_lfft_is_power_2(0x0006) == false, "False assumption: _lfft_is_power_2(0x0006) == false\n", i);
    fail_unless(_lfft_is_power_2(0x0007) == false, "False assumption: _lfft_is_power_2(0x0007) == false\n", i);
    fail_unless(_lfft_is_power_2(0x0009) == false, "False assumption: _lfft_is_power_2(0x0009) == false\n", i);
    fail_unless(_lfft_is_power_2(0x000A) == false, "False assumption: _lfft_is_power_2(0x000A) == false\n", i);
    fail_unless(_lfft_is_power_2(0x000B) == false, "False assumption: _lfft_is_power_2(0x000B) == false\n", i);
    fail_unless(_lfft_is_power_2(0x000C) == false, "False assumption: _lfft_is_power_2(0x000C) == false\n", i);
    fail_unless(_lfft_is_power_2(0x000D) == false, "False assumption: _lfft_is_power_2(0x000D) == false\n", i);
    fail_unless(_lfft_is_power_2(0x000E) == false, "False assumption: _lfft_is_power_2(0x000E) == false\n", i);
    fail_unless(_lfft_is_power_2(0x000F) == false, "False assumption: _lfft_is_power_2(0x000F) == false\n", i);

    // test all numbers 2^x from 0x0000 to 0xFFFF
    for(i = 0; i < 8; ++i)
    {
        fail_unless(_lfft_is_power_2(1<<i) == true, "False assumption: _lfft_is_power_2(1<<%d) == true\n", i);
    }
}
END_TEST

START_TEST(test_f_lfft_isqrt)
{
#define ROUND(x) ((uint32_t)(floor(x)+0.5f))

// difference between lfft_isqrt and sqrt is maximal 1
#define TEST_LOOP(start, stop) \
for(i = start; i < stop; ++i) \
{ \
    lfft_sqrt_result = lfft_isqrt(i); \
    sqrt_result = ROUND(sqrt((float)i)); \
    fail_unless(abs(lfft_sqrt_result-sqrt_result) <= 1, "False assumption: lfft_isqrt(%"PRIu16")=%d | sqrt(%"PRIu16")=%d\n", i, lfft_sqrt_result, i, sqrt_result); \
}

    uint32_t i;
    uint32_t lfft_sqrt_result;
    uint32_t sqrt_result;

    TEST_LOOP(0, 0xFFFF);
    TEST_LOOP(0xFFFF0000, 0xFFFFFFFF);

#undef TEST_LOOP
#undef ROUND
}
END_TEST

Suite* a_suite()
{
    Suite *suite = suite_create ("lfft");
    TCase *tcase = tcase_create ("case");
    tcase_add_test(tcase, test_f__lfft_is_power_2);
    tcase_add_test(tcase, test_f_lfft_isqrt);
    //tcase_add_loop_test (tcase, test_f__lfft_is_power_2, 0x0000, 0x0002);
    suite_add_tcase(suite, tcase);

    return suite;
}

int main()
{
    int number_failed;
    Suite *suite = a_suite();
    SRunner *runner = srunner_create(suite);
    srunner_run_all(runner, CK_NORMAL);
    number_failed = srunner_ntests_failed(runner);
    srunner_free(runner);
    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}

