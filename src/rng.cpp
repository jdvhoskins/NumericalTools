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

#include "rng.h"
#include <chrono>
#include <cassert>

#ifndef PI

#define PI 3.141592654

#endif // !PI


rng::rng ()
{
	// Randomly determined at the time this was written by rolling, for each set of 4 bits, 1D20 - 1 (reroll 16+)
	m_llState = 0x74a1eb62bf6b40f9; 
}

//--------------------

rng::rng (long long seed)
{
	assert (seed != 0);
	m_llState = seed;

	if (seed == 0)
	{
		// Randomly determined at the time this was written by rolling, for each set of 4 bits, 1D20 - 1 (reroll 16+)
		m_llState = 0x74a1eb62bf6b40f9;
	}
}

//--------------------

void rng::SetState (long long seed)
{
	assert (seed != 0);

	if (seed == 0)
	{
		return;
	}

	m_llState = seed;
}

//--------------------

void rng::SetRandomState ()
{
	m_llState = std::chrono::high_resolution_clock::now ().time_since_epoch ().count ();
}

//--------------------

void rng::GetNext ()
{
	// Shift values were randomly determined at the time this was written by rolling 2D30
	m_llState ^= m_llState << 21;
	m_llState ^= m_llState >> 27;
	m_llState ^= m_llState << 37;
}

//--------------------

bool rng::RandomBool ()
{
	GetNext ();
	return (m_llState % 2) == 0;
}

//--------------------

char rng::RandomInt8 ()
{
	GetNext ();
	return (char) (m_llState * NL_CONST);
}

//--------------------

char rng::RandomInt8 (char cMax)
{
	assert (cMax > 0);
	// Negating the sign bit ensures that the result is positive definite
	return (RandomInt8 () & 0xEF) % cMax;
}

//--------------------

char rng::RandomInt8 (char cMin, char cMax)
{
	assert (cMax > cMin);
	// Negating the sign bit ensures that the result is positive definite
	return cMin + (RandomInt8 () & 0xEF) % (cMax - cMin);
}


//--------------------

unsigned char rng::RandomUInt8 ()
{
	GetNext ();
	return (unsigned char) (m_llState * NL_CONST);
}

//--------------------

unsigned char rng::RandomUInt8 (unsigned char ucMax)
{
	assert (ucMax > 0);
	return RandomUInt8 () % ucMax;
}

//--------------------

unsigned char rng::RandomUInt8 (unsigned char ucMin, unsigned char ucMax)
{
	assert (ucMax > ucMin);
	return ucMin + RandomUInt8 () % (ucMax - ucMin);
}

//--------------------

short rng::RandomInt16 ()
{
	GetNext ();
	return (short) (m_llState * NL_CONST);
}

//--------------------

short rng::RandomInt16 (short sMax)
{
	assert (sMax > 0);
	// Negating the sign bit ensures that the result is positive definite
	return (RandomInt16 () & 0xEFFF) % sMax;
}

//--------------------

short rng::RandomInt16 (short sMin, short sMax)
{
	assert (sMax > sMin);
	// Negating the sign bit ensures that the result is positive definite
	return sMin + (RandomInt16 () & 0xEFFF) % (sMax - sMin);
}

//--------------------

unsigned short rng::RandomUInt16 ()
{
	GetNext ();
	return (unsigned short) (m_llState * NL_CONST);
}

//--------------------

unsigned short rng::RandomUInt16 (unsigned short usMax)
{
	assert (usMax > 0);
	return RandomUInt16 () % usMax;
}

//--------------------

unsigned short rng::RandomUInt16 (unsigned short usMin, unsigned short usMax)
{
	assert (usMax > usMin);
	return usMin + RandomUInt16 () % (usMax - usMin);
}

//--------------------

int rng::RandomInt32 ()
{
	GetNext ();
	return (int) (m_llState * NL_CONST);
}

//--------------------

int rng::RandomInt32 (int iMax)
{
	assert (iMax > 0);
	// Negating the sign bit ensures that the result is positive definite
	return (RandomInt32 () & 0xEFFFFFFF) % iMax;
}

//--------------------

int rng::RandomInt32 (int iMin, int iMax)
{
	assert (iMax > iMin);
	// Negating the sign bit ensures that the result is positive definite
	return iMin + (RandomInt32 () & 0xEFFFFFFF) % (iMax - iMin);
}

//--------------------

unsigned int rng::RandomUInt32 ()
{
	GetNext ();
	return (unsigned int) (m_llState * NL_CONST);
}

//--------------------

unsigned int rng::RandomUInt32 (unsigned int uiMax)
{
	assert (uiMax > 0);
	return RandomUInt32 () % uiMax;
}

//--------------------

unsigned int rng::RandomUInt32 (unsigned int uiMin, unsigned int uiMax)
{
	assert (uiMax > uiMin);
	return uiMin + RandomUInt32 () % (uiMax - uiMin);
}

//--------------------

long long rng::RandomInt64 ()
{
	GetNext ();
	return (long long) (m_llState * NL_CONST);
}

//--------------------

long long rng::RandomInt64 (long long llMax)
{
	assert (llMax > 0);
	// Negating the sign bit ensures that the result is positive definite
	return (RandomInt64 () & 0xEFFFFFFFFFFFFFFF) % llMax;
}

//--------------------

long long rng::RandomInt64 (long long llMin, long long llMax)
{
	assert (llMax > llMin);
	// Negating the sign bit ensures that the result is positive definite
	return llMin + (RandomInt64 () & 0xEFFFFFFFFFFFFFFF) % (llMax - llMin);
}

//--------------------

unsigned long long rng::RandomUInt64 ()
{
	GetNext ();
	return (unsigned long long) (m_llState * NL_CONST);
}

//--------------------

unsigned long long rng::RandomUInt64 (unsigned long long ullMax)
{
	assert (ullMax > 0);
	return RandomUInt64 () % ullMax;
}

//--------------------

unsigned long long rng::RandomUInt64 (unsigned long long ullMin, unsigned long long ullMax)
{
	assert (ullMax > ullMin);
	return ullMin + RandomUInt64 () % (ullMax - ullMin);
}

//--------------------

float rng::RandomFloat ()
{
	// 32-bit floats use 23 bits for the mantissa. As such, 2^23 - 1 is used as the maximum mantissa value to avoid rounding issues.
	unsigned int uiMantissaMax = 0x7FFFFF;
	unsigned int uiRandomInt = RandomUInt32 () & uiMantissaMax;
	return (float) uiRandomInt / (float) uiMantissaMax;
}

//--------------------

float rng::RandomFloat (float fMax)
{
	return fMax * RandomFloat ();
}

//--------------------

float rng::RandomFloat (float fMin, float fMax)
{
	assert (fMax > fMin);
	return fMin + RandomFloat () * (fMax - fMin);
}

//--------------------

double rng::RandomDouble ()
{
	// 64-bit floats use 52 bits for the mantissa. As such, 2^52 - 1 is used as the maximum mantissa value to avoid rounding issues.
	unsigned long long ullMantissaMax = 0xFFFFFFFFFFFFF;
	unsigned long long ullRandomInt = RandomUInt64 () & ullMantissaMax;
	return (double) ullRandomInt / (double) ullMantissaMax;
}

//--------------------

double rng::RandomDouble (double dMax)
{
	return dMax * RandomDouble ();
}

//--------------------

double rng::RandomDouble (double dMin, double dMax)
{
	assert (dMax > dMin);
	return dMin + RandomDouble () * (dMax - dMin);
}

//--------------------

double rng::RandomNormal ()
{
	// This algorithm implemented here is from:
	// W.H. Press, et al., Numerical Recipes in C, Cambridge University Press, 1992.

	static bool bHaveNextValue = false;
	static double dNextValue;
	double dX;
	double dY;
	double dRadiusSquared;
	double dBoxMuller;

	if (!bHaveNextValue) // Compute 2 new values
	{
		do {
			dX = RandomDouble (2.0) - 1.0;
			dY = RandomDouble (2.0) - 1.0;
			dRadiusSquared = dX * dX + dY * dY;
		} while (dRadiusSquared >= 1.0 || dRadiusSquared == 0.0);	// If (dX, dY) is outside the unit circle or is (0, 0), then resample.
		dBoxMuller = sqrt (-2.0 * log (dRadiusSquared) / dRadiusSquared);	// Use the Box-Muller transformation to get two normally distributed values.
		dNextValue = dX * dBoxMuller;  // Save one for next time...
		bHaveNextValue = true;
		return dY * dBoxMuller;  // ...and return the other.
	}
	else // We have the second value from the previous function call
	{
		bHaveNextValue = false;
		return dNextValue;
	}
}

//--------------------

double rng::RandomNormal (double dMean, double dStdv)
{
	return RandomNormal () * dStdv + dMean;
}

//--------------------

double rng::RandomGamma (int iOrder)
{
	// This algorithm implemented here is from:
	// W.H. Press, et al., Numerical Recipes in C, Cambridge University Press, 1992.
	// Comments have been added to fill in the gaps in their derivation.

	// The probability density function for the Gamma distribution (with scale parameter = 1) is given by:
	//
	//		P_a(x) = x^(a-1) * exp(-x) / GammaFn(a)
	//
	// where "a" is the order of the distribution. For a = 1, this yields the exponential distribution:
	// 
	//		P_1(x) = exp(-x)
	//


	double dMu;
	double dRatio;
	double dSigma;
	double dRandX;
	double dRandY;
	double dX;
	double dTangent;
	double dRadiusSquared;

	assert (iOrder >= 1);
	if (iOrder < 1)
	{
		return -1.0;
	}
	
	if (iOrder < 6) //Use direct method.
	{
		dX = 1.0;
		for (int iLoop = 1; iLoop <= iOrder; ++iLoop)
		{
			dX *= RandomDouble (); // Get product of iOrder uniformly distributed values
		}
		dX = -log (dX); // Take the log of that product and return
	}
	else // Use the rejection method by comparing to a lorentzian
	{
		// The equation for a lorentzian is:
		//
		//		L(x) = C / (1 + (x - Mu)^2 / Sigma^2)
		//
		// For a gamma distribution of order a, the optimal values for Mu, Sigma, and C were found by Aherns and Dieter to be:
		//
		//		Mu = a - 1
		//		Sigma = sqrt(2 * a - 1)
		//		C = PI * Sigma * exp(-Mu * (1 - log(Mu))) / GammaFn(a)
		//
		// J.H. Ahrens and U. Dieter. Generating gamma variates by a modified rejectiontechnique. Comm. ACM, 25(1):47–54, January 1982.

		// Thus (via the wonders of algebra) the ratio of the gamma distribution to the optimal lorentzian is given by:
		//
		//		ratio = (1 + x^2) * exp (Mu * log (x / Mu) - Sigma * x)
		//
		// where x = Sigma tan(2 * PI * RandomFloat) + Mu. Thankfully GammaFn(a) is canceled out.
		do
		{
			do
			{
				do { // Compute tan(2 * PI * RandomFloat) = RandY / RandX,  within the unit circle
					dRandX = RandomFloat ();
					dRandY = RandomFloat (2.0) - 1.0;
					dRadiusSquared = dRandX * dRandX + dRandY * dRandY;
				} while (dRadiusSquared > 1.0 || dRadiusSquared == 0.0); // Get a random point in the 1st or 4th quadrant that is within the unit circle
				
				dTangent = dRandY / dRandX;
				dMu = (double) iOrder - 1.0;
				dSigma = sqrt (2.0 * (double) iOrder - 1.0);
				dX = dSigma * dTangent + dMu;
			} while (dX <= 0.0); // If dX is negative, it is outside the gamma distribution and is rejected.

			dRatio = (1.0 + dTangent * dTangent) * exp (dMu * log (dX / dMu) - dSigma * dTangent); // Ratio of gamma distribution to the lorentz distribution.
		} while (RandomFloat () > dRatio); // If a second uniformly distributed value is greater than this ratio, reject dX.
	}

	return dX;
}

//--------------------

double rng::RandomPoisson (double dMean)
{
	// This algorithm is adapted from:
	// W.H. Press, et al., Numerical Recipes in C, Cambridge University Press, 1992.

	static double dSigma = -1.0;
	static double dLogOfMean = -1.0;
	static double dCompare = -1.0;
	static double dOldMean = -1.0; // dOldMean is saved to determine whether or not dMean has changed since last call.
	static double dC = -1;
	double dRetValue, dProbability, dTangent;
	
	if (dMean < 0)
	{
		return -1.0;
	}

	if (dMean < 12.0)	// Use direct method.
	{
		// Compute a new exponential value only if the mean of the distribution has changed since the previous function call.
		if (dMean != dOldMean) 
		{
			dOldMean = dMean;
			dCompare = exp (-dMean);
		}
		dRetValue = -1;
		dProbability = 1.0;
		do {
			// Instead of adding exponential deviates it is equivalent to multiply uniform deviates.
			// We never actually have to take the log, we merely compare the cumulative probability to the precomputed exponential, dCompare.
			++dRetValue;
			dProbability *= RandomDouble ();
		} while (dProbability > dCompare);

	}
	else {
		// Use rejection method.
		if (dMean != dOldMean) // If dMean has changed since the last call, then precompute some functions that occur below.
		{
			dOldMean = dMean;
			dSigma = sqrt (dMean);
			dLogOfMean = log (dMean);
			dC = pow (dMean, dMean) * exp (-dMean) / tgamma (dMean + 1);
		}
		do {
			do {
				// dTangent is a deviate from a Lorentzian comparison function given by:
				//
				//		L(x) = C / (1 + (x - Mu)^2 / Sigma^2)
				//
				// where Mu = dMean and Sigma = sqrt(dMean).
				// The inverse of the integral of L(x) is given by:
				//
				//		y = sigma * tan(x / (C * sigma)) + Mu
				//		y = sqrt(dMean) * tan(Random) + dMean
				//
				// Because C only apears in the argument of the tangent function, it can be subsumed by the uniform random value used there.
				// As such, it is constrained only by the fact that L(x) must be greater than P(x) for all x. A value of C = P_Max for a
				// given dMean will result in parts of L(x) being less than P(x). C = 2 * P_Max yields a ratio of exactly 1 for P(x)/L(x)
				// at x = 0 and dMean = 1, and is less than 1 everywhere else. C = 3* P_Max yields ratios that are everywhere less than 1.

				dTangent = tan (PI * RandomDouble ());
				dRetValue = dSigma * dTangent + dMean;	//dRetValue is the result of the inverse CDF of the above lorentzian.
			} while (dRetValue < 0.0);	// If dRetValue is negative, reject it.

			dRetValue = floor (dRetValue); // Take the floor of dRetValue because we want integers only.
			dCompare = (1.0 / dC) * (1.0 + dTangent * dTangent) * exp (dRetValue * dLogOfMean - dMean - lgamma (dRetValue + 1.0)); // The comparison value is the ratio P(x)/L(x)
			} while (RandomFloat () > dCompare); // If a second uniformly distributed number is greater than dCompare, then reject.
	}
	return dRetValue;
}

//--------------------

// TODO: move this to a special functions library
double rng::Factorial (int iN)
{
	// Factorial is a special case of the Gamma function. For small N a lookup table is used. For large N the Gamma function is used.

	// The largest factorial that can be represented by an unsigned long long is 19!
	// The largest factorial that can be represented by a float is 34!
	// The largest factorial that can be represented by a double is 170!
	// The return value will be a double regardles of whether or not a double is actually needed to hold the value.

	// This algorithm implemented here is from:
	// W.H. Press, et al., Numerical Recipes in C, Cambridge University Press, 1992.


	if (iN > 170 || iN < 0)
	{
		return -1.0;
	}

	if (iN > 34)
	{
		return exp (lgamma(iN + 1));
	}

	int iPrevious;
	static int iLast = 1;
	static double pdFactorials[35] = {1.0, 1.0};
	while (iLast < iN)
	{
		iPrevious = iLast++;
		pdFactorials[iLast] = pdFactorials[iPrevious] * (double) iLast;
	}

	return pdFactorials[iN];
}

//--------------------
