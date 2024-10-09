/******************************************************************************
 * BIKE -- Bit Flipping Key Encapsulation
 *
 * Copyright (c) 2021 Nir Drucker, Shay Gueron, Rafael Misoczki, Tobias Oder,
 * Tim Gueneysu, Jan Richter-Brockmann.
 * Contact: drucker.nir@gmail.com, shay.gueron@gmail.com,
 * rafaelmisoczki@google.com, tobias.oder@rub.de, tim.gueneysu@rub.de,
 * jan.richter-brockmann@rub.de.
 *
 * Permission to use this code for BIKE is granted.
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice,
 *   this list of conditions and the following disclaimer.
 *
 * * Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in the
 *   documentation and/or other materials provided with the distribution.
 *
 * * The names of the contributors may not be used to endorse or promote
 *   products derived from this software without specific prior written
 *   permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHORS ""AS IS"" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS CORPORATION OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 ******************************************************************************/

#ifndef __DEFS_H_INCLUDED__
#define __DEFS_H_INCLUDED__

#include <assert.h>
#include <math.h>

////////////////////////////////////////////
//         BIKE main parameters
///////////////////////////////////////////

// UNCOMMENT TO ENABLE BANDWIDTH OPTIMISATION FOR BIKE-3:
//#define BANDWIDTH_OPTIMIZED

// BIKE shared-secret size:
#define ELL_BITS  256ULL
#define ELL_SIZE (ELL_BITS/8)

////////////////////////////////////////////
// Implicit Parameters (do NOT edit below)
///////////////////////////////////////////

// select the max between a and b:
#define MAX(a,b) ((a)>(b))?(a):(b)

// Parameters for BGF Decoder:
#define tau 3

// Divide by the divider and round up to next integer:
#define DIVIDE_AND_CEIL(x, divider)  ((x/divider) + (x % divider == 0 ? 0 : 1ULL))

// Round the size to the nearest byte.
// SIZE suffix, is the number of bytes (uint8_t).
#define N_BITS   (R_BITS*2)
#define R_SIZE   DIVIDE_AND_CEIL(R_BITS, 8ULL)
#define N_SIZE   DIVIDE_AND_CEIL(N_BITS, 8ULL)
#define R_DQWORDS DIVIDE_AND_CEIL(R_SIZE, 16ULL)

inline double VAR_TH_FCT( int x ){
    double pi1 = (x + THR_X) / T1 / DV;
    double pi0 = (DV * 2 * x - THR_X) / (N_BITS - T1) / DV;
    assert(pi1 >= pi0);
    double T = ceil((log((N_BITS - T1) / T1) + DV * log((1 - pi0) / (1 - pi1))) / (log(pi1 / pi0) + log((1 - pi0) / (1 - pi1))));
    return fmax(T, (DV + 1) / 2.0);
}

////////////////////////////////////////////
//             Debug
///////////////////////////////////////////

#ifndef VERBOSE
#define VERBOSE 0
#endif

#if (VERBOSE == 3)
#define MSG(...)     { printf(__VA_ARGS__); }
#define DMSG(...)    MSG(__VA_ARGS__)
#define EDMSG(...)   MSG(__VA_ARGS__)
#define SEDMSG(...)  MSG(__VA_ARGS__)
#elif (VERBOSE == 2)
#define MSG(...)     { printf(__VA_ARGS__); }
#define DMSG(...)    MSG(__VA_ARGS__)
#define EDMSG(...)   MSG(__VA_ARGS__)
#define SEDMSG(...)
#elif (VERBOSE == 1)
#define MSG(...)     { printf(__VA_ARGS__); }
#define DMSG(...)    MSG(__VA_ARGS__)
#define EDMSG(...)
#define SEDMSG(...)
#else
#define MSG(...)     //{ printf(__VA_ARGS__); }
#define DMSG(...)
#define EDMSG(...)
#define SEDMSG(...)
#endif

////////////////////////////////////////////
//              Printing
///////////////////////////////////////////

// Show timer results in cycles.
#define RDTSC

//#define PRINT_IN_BE
//#define NO_SPACE
//#define NO_NEWLINE

////////////////////////////////////////////
//              Testing
///////////////////////////////////////////
#define NUM_OF_CODE_TESTS       100ULL
#define NUM_OF_ENCRYPTION_TESTS 100ULL

#endif //__TYPES_H_INCLUDED__

