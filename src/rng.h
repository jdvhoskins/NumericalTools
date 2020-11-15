// Copyright 2020 Jeffrey D. V. Hoskins
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is furnished
// to do so, subject to the following conditions :
// 
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
// FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE AUTHORS OR
// COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN
// AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// This library contains random number generators for both uniform distribuions
// as well Poisson, Gamma, and Normal (Gaussian) distributions. 

#pragma once

#include <limits>

class rng
{
public:
	rng ();
	rng (long long seed);

	// seed != 0
	void SetState (long long seed);
	void SetRandomState ();
	void GetNext ();


	// Uniform distributions

	bool RandomBool ();

	// Returns a random 8 bit integer in range...
	char RandomInt8 ();  // [CHAR_MIN, CHAR_MAX]
	char RandomInt8 (char cMax);  // [0, cMax)
	char RandomInt8 (char cMin, char cMax);   // [cMin, cMax)

	// Returns a random unsigned 8 bit integer in range...
	unsigned char RandomUInt8 ();  // [0, UCHAR_MAX]
	unsigned char RandomUInt8 (unsigned char ucMax);  // [0, ucMax)
	unsigned char RandomUInt8 (unsigned char ucMin, unsigned char ucMax);   // [ucMin, ucMax)

	// Returns a random 16 bit integer in range...
	short RandomInt16 ();  // [SHRT_MIN, SHRT_MAX]
	short RandomInt16 (short sMax);  // [0, sMax)
	short RandomInt16 (short sMin, short sMax);   // [sMin, sMax)

	// Returns a random unsigned 16 bit integer in range...
	unsigned short RandomUInt16 ();  // [0, USHRT_MAX]
	unsigned short RandomUInt16 (unsigned short usMax);  // [0, usMax)
	unsigned short RandomUInt16 (unsigned short usMin, unsigned short usMax);   // [usMin, usMax)

	// Returns a random 32 bit integer in range...
	int RandomInt32 ();  // [INT_MIN, INT_MAX]
	int RandomInt32 (int iMax);  // [0, iMax)
	int RandomInt32 (int iMin, int iMax);   // [iMin, iMax)

	// Returns a random unsigned 32 bit integer in range... 
	unsigned int RandomUInt32 (); // [0, UINT_MAX]
	unsigned int RandomUInt32 (unsigned int uiMax);  // [0, uiMax)
	unsigned int RandomUInt32 (unsigned int uiMin, unsigned int uiMax);   // [uiMin, uiMax)

	// Returns a random 64 bit integer in range...
	long long RandomInt64 ();  // [LLONG_MIN, LLONG_MAX]
	long long RandomInt64 (long long llMax);  // [0, llMax)
	long long RandomInt64 (long long llMin, long long llMax);   // [lMin, llMax)

	// Returns a random 64 bit unsigned integer in range...
	unsigned long long RandomUInt64 ();  // [0, ULLONG_MAX]
	unsigned long long RandomUInt64 (unsigned long long ullMax);  // [0, ullMax)
	unsigned long long RandomUInt64 (unsigned long long ullMin, unsigned long long ullMax);   // [ullMin, ullMax)


	// Returns a random float in range...
	float RandomFloat ();  // [0, 1)
	float RandomFloat (float fMax);  // [0, fMax)
	float RandomFloat (float fMin, float fMax);  // [fMin, fMax)

	// Returns a random double in range...
	double RandomDouble ();  // [0, 1)
	double RandomDouble (double dMax);  // [0, dMax)
	double RandomDouble (double dMin, double dMax);  // [dMin, dMax)


	// Non-uniform distributions

	// Returns a random integer based on a normal distribution with a mean value of 0.0 and a standard deviation of 1.0
	double RandomNormal ();
	double RandomNormal (double dMean, double dStdv);
	double RandomGamma (int iOrder);
	// Returns a random integer based on the poisson distribution with a mean value of ullMean.
	double RandomPoisson (double dMean);


	// TODO: Move this to a special functions library
	double Factorial (int iX);
	
private:
	long long m_llState;

	// Randomly determined at the time this was writted by rolling, for each set of 4 bits, 1D20 - 1 (reroll 16+)
	static const long long NL_CONST = 0x52ee31c34f014733;
};