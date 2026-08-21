/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    CRandomSFMT

\*---------------------------------------------------------------------------*/
#include "CRandomSFMT.H"

namespace Foam
{

void CRandomSFMT::RandomInit(int seed) {
   // Re-seed
   uint32_t i;                         // Loop counter
   uint32_t y = seed;                  // Temporary
   uint32_t statesize = SFMT_N*4;      // Size of state vector
   if (UseMother) statesize += 5;      // Add states for Mother-Of-All generator

   // Fill state vector with random numbers from seed
   ((uint32_t*)state)[0] = y;
   const uint32_t factor = 1812433253U;// Multiplication factor

   for (i = 1; i < statesize; i++) {
      y = factor * (y ^ (y >> 30)) + i;
      ((uint32_t*)state)[i] = y;
   }

   // Further initialization and period certification
   Init2();
}


// Functions used by CRandomSFMT::RandomInitByArray
static uint32_t func1(uint32_t x) {
    return (x ^ (x >> 27)) * 1664525U;
}

static uint32_t func2(uint32_t x) {
    return (x ^ (x >> 27)) * 1566083941U;
}

void CRandomSFMT::RandomInitByArray(int const seeds[], int NumSeeds) {
   // Seed by more than 32 bits
   uint32_t i, j, count, r, lag;

   if (NumSeeds < 0) NumSeeds = 0;

   const uint32_t size = SFMT_N*4; // number of 32-bit integers in state

   // Typecast state to uint32_t *
   uint32_t * sta = (uint32_t*)state;

   if (size >= 623) {
      lag = 11;} 
   else if (size >= 68) {
      lag = 7;}
   else if (size >= 39) {
      lag = 5;}
   else {
      lag = 3;
   }
   const uint32_t mid = (size - lag) / 2;

   if (static_cast<uint32_t>(NumSeeds) + 1 > size) {
      count = static_cast<uint32_t>(NumSeeds);
   }
   else {
      count = size - 1;
   }
#if 0
   // Original code. Argument to func1 is constant!
   for (i = 0; i < size; i++) sta[i] = 0x8B8B8B8B;
   r = func1(sta[0] ^ sta[mid] ^ sta[size - 1]);
   sta[mid] += r;
   r += NumSeeds;
   sta[mid + lag] += r;
   sta[0] = r;
#else
   // 1. loop: Fill state vector with random numbers from NumSeeds
   const uint32_t factor = 1812433253U;// Multiplication factor
   r = static_cast<uint32_t>(NumSeeds);
   for (i = 0; i < SFMT_N*4; i++) {
      r = factor * (r ^ (r >> 30)) + i;
      sta[i] = r;
   }

#endif

   // 2. loop: Fill state vector with random numbers from seeds[]
   for (i = 1, j = 0; j < count; j++) {
      r = func1(sta[i] ^ sta[(i + mid) % size] ^ sta[(i + size - 1) % size]);
      sta[(i + mid) % size] += r;
      if (j < static_cast<uint32_t>(NumSeeds)) r += static_cast<uint32_t>(seeds[j]);
      r += i;
      sta[(i + mid + lag) % size] += r;
      sta[i] = r;
      i = (i + 1) % size;
   }

   // 3. loop: Randomize some more
   for (j = 0; j < size; j++) {
      r = func2(sta[i] + sta[(i + mid) % size] + sta[(i + size - 1) % size]);
      sta[(i + mid) % size] ^= r;
      r -= i;
      sta[(i + mid + lag) % size] ^= r;
      sta[i] = r;
      i = (i + 1) % size;
   }
   if (UseMother) {
      // 4. loop: Initialize MotherState
      for (j = 0; j < 5; j++) {
         r = func2(r) + j;
         MotherState[j] = r + sta[2*j];
      }
   }
   
   // Further initialization and period certification
   Init2();
}


void CRandomSFMT::Init2() {
   // Various initializations and period certification
   uint32_t i, j, temp;

   // Initialize mask
   static const uint32_t maskinit[4] = {SFMT_MASK};
   mask = _mm_loadu_si128((__m128i*)maskinit);

   // Period certification
   // Define period certification vector
   static const uint32_t parityvec[4] = {SFMT_PARITY};

   // Check if parityvec & state[0] has odd parity
   temp = 0;
   for (i = 0; i < 4; i++) {
      temp ^= parityvec[i] & ((uint32_t*)state)[i];
   }
   for (i = 16; i > 0; i >>= 1) temp ^= temp >> i;
   if (!(temp & 1)) {
      // parity is even. Certification failed
      // Find a nonzero bit in period certification vector
      for (i = 0; i < 4; i++) {
         if (parityvec[i]) {
            for (j = 1; j; j <<= 1) {
               if (parityvec[i] & j) {
                  // Flip the corresponding bit in state[0] to change parity
                  ((uint32_t*)state)[i] ^= j;
                  // Done. Exit i and j loops
                  i = 5;  break;
               }
            }
         }
      }
   }
   // Generate first random numbers and set ix = 0
   Generate();
}


// Subfunction for the sfmt algorithm
static inline __m128i sfmt_recursion(__m128i const &a, __m128i const &b, 
__m128i const &c, __m128i const &d, __m128i const &mask) {
    __m128i a1, b1, c1, d1, z1, z2;
    b1 = _mm_srli_epi32(b, SFMT_SR1);
    a1 = _mm_slli_si128(a, SFMT_SL2);
    c1 = _mm_srli_si128(c, SFMT_SR2);
    d1 = _mm_slli_epi32(d, SFMT_SL1);
    b1 = _mm_and_si128(b1, mask);
    z1 = _mm_xor_si128(a, a1);
    z2 = _mm_xor_si128(b1, d1);
    z1 = _mm_xor_si128(z1, c1);
    z2 = _mm_xor_si128(z1, z2);
    return z2;
}

void CRandomSFMT::Generate() {
   // Fill state array with new random numbers
   int i;
   __m128i r, r1, r2;

   r1 = state[SFMT_N - 2];
   r2 = state[SFMT_N - 1];
   for (i = 0; i < SFMT_N - SFMT_M; i++) {
      r = sfmt_recursion(state[i], state[i + SFMT_M], r1, r2, mask);
      state[i] = r;
      r1 = r2;
      r2 = r;
   }
   for (; i < SFMT_N; i++) {
      r = sfmt_recursion(state[i], state[i + SFMT_M - SFMT_N], r1, r2, mask);
      state[i] = r;
      r1 = r2;
      r2 = r;
   }
   ix = 0;
}

uint32_t CRandomSFMT::BRandom() {
   // Output 32 random bits
   uint32_t y;

   if (ix >= SFMT_N*4) {
      Generate();
   }
   y = ((uint32_t*)state)[ix++];
   if (UseMother) y += MotherBits();
   return y;
}

uint32_t CRandomSFMT::MotherBits() {
   // Get random bits from Mother-Of-All generator
   uint64_t sum;
   sum = 
      static_cast<uint64_t>(2111111111U) * static_cast<uint64_t>(MotherState[3]) +
      static_cast<uint64_t>(1492) * static_cast<uint64_t>(MotherState[2]) +
      static_cast<uint64_t>(1776) * static_cast<uint64_t>(MotherState[1]) +
      static_cast<uint64_t>(5115) * static_cast<uint64_t>(MotherState[0]) +
      static_cast<uint64_t>(MotherState[4]);
   MotherState[3] = MotherState[2];  
   MotherState[2] = MotherState[1];  
   MotherState[1] = MotherState[0];
   MotherState[4] = static_cast<uint32_t>(sum >> 32);       // Carry
   MotherState[0] = static_cast<uint32_t>(sum);               // Low 32 bits of sum
   return MotherState[0];
}

int  CRandomSFMT::IRandom (int min, int max) {
   // Output random integer in the interval min <= x <= max
   // Slightly inaccurate if (max-min+1) is not a power of 2
   if (max <= min) {
      if (max == min) return min; else return 0x80000000;
   }
   // Assume 64 bit integers supported. Use multiply and shift method
   uint32_t interval;                  // Length of interval
   uint64_t longran;                   // Random bits * interval
   uint32_t iran;                      // Longran / 2^32

   interval = static_cast<uint32_t>(max - min + 1);
   longran  = static_cast<uint64_t>(BRandom()) * interval;
   iran = static_cast<uint32_t>(longran >> 32);
   // Convert back to signed and return result
   return static_cast<int32_t>(iran) + min;
}

int  CRandomSFMT::IRandomX (int min, int max) {
   // Output random integer in the interval min <= x <= max
   // Each output value has exactly the same probability.
   // This is obtained by rejecting certain bit values so that the number
   // of possible bit values is divisible by the interval length
   if (max <= min) {
      if (max == min) {
         return min;                   // max == min. Only one possible value
      }
      else {
         return 0x80000000;            // max < min. Error output
      }
   }
   // Assume 64 bit integers supported. Use multiply and shift method
   uint32_t interval;                  // Length of interval
   uint64_t longran;                   // Random bits * interval
   uint32_t iran;                      // Longran / 2^32
   uint32_t remainder;                 // Longran % 2^32

   interval = static_cast<uint32_t>(max - min + 1);
   if (interval != LastInterval) {
      // Interval length has changed. Must calculate rejection limit
      // Reject when remainder = 2^32 / interval * interval
      // RLimit will be 0 if interval is a power of 2. No rejection then.
      RLimit = static_cast<uint32_t>((static_cast<uint64_t>(1) << 32) / interval) * interval - 1;
      LastInterval = interval;
   }
   do { // Rejection loop
      longran  = static_cast<uint64_t>(BRandom()) * interval;
      iran = static_cast<uint32_t>(longran >> 32);
      remainder = static_cast<uint32_t>(longran);
   } while (remainder > RLimit);
   // Convert back to signed and return result
   return static_cast<int32_t>(iran) + min;
}

double CRandomSFMT::Random() {
   // Output random floating point number
   if (ix >= SFMT_N*4-1) {
      // Make sure we have at least two 32-bit numbers
      Generate();
   }
   uint64_t r = *(uint64_t*)((uint32_t*)state+ix);
   ix += 2;
   if (UseMother) {
      // We need 53 bits from Mother-Of-All generator
      // Use the regular 32 bits and the the carry bits rotated
      uint64_t r2 = static_cast<uint64_t>(MotherBits()) << 32;
      r2 |= (MotherState[4] << 16) | (MotherState[4] >> 16);
      r += r2;
   }
   // 53 bits resolution:
   // return (int64_t)(r >> 11) * (1./(67108864.0*134217728.0)); // (r >> 11)*2^(-53)
   // 52 bits resolution for compatibility with assembly version:
   return (int64_t)(r >> 12) * (1./(67108864.0*67108864.0));  // (r >> 12)*2^(-52)
}

} // namespace Foam

// ************************************************************************* //
