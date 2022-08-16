#pragma once
/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS, * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

 /*
  * This is not the original file distributed by the Apache Software Foundation
  * It has been modified by the Hipparchus project
  */
  //package org.hipparchus.random;

  //import java.io.Serializable;

  //import org.hipparchus.util.FastMath;

  /**
   * This class : a powerful pseudo-random number generator
   * developed by Makoto Matsumoto and Takuji Nishimura during
   * 1996-1997.
   * <p>
   * <b>Caveat:</b> It is recommended to use one of WELL generators rather
   * than the Mersenne_Twister generator (see
   * <a href="http://www.iro.umontreal.ca/~panneton/WELLRNG.html">
   * this paper</a> for more information).
   * <p>
   * This generator features an extremely long period
   * (2<sup>19937</sup>-1) and 623-dimensional equidistribution up to 32
   * bits accuracy. The home page for this generator is located at <a
   * href="http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html">
   * http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html</a>.
   * <p>
   * This generator is described in a paper by Makoto Matsumoto and
   * Takuji Nishimura in 1998: <a
   * href="http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/ARTICLES/mt.pdf">Mersenne
   * Twister: A 623-Dimensionally Equidistributed Uniform Pseudo-Random
   * Number Generator</a>, ACM Transactions on Modeling and Computer
   * Simulation, Vol. 8, No. 1, January 1998, pp 3--30.
   * <p>
   * This class is mainly a Java port of the 2002-01-26 version of
   * the generator written in C by Makoto Matsumoto and Takuji
   * Nishimura. Here is their original copyright:
   * <p>
   * <table border="0" width="80%" cellpadding="10" align="center" bgcolor="#E0E0E0">
   * <tr><td>Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura, *     All rights reserved.</td></tr>
   * <tr><td>Redistribution and use in source and binary forms, with or without
   * modification, are permitted provided that the following conditions
   * are met:
   * <ol>
   *   <li>Redistributions of source code must retain the above copyright
   *       notice, this list of conditions and the following disclaimer.</li>
   *   <li>Redistributions in binary form must reproduce the above copyright
   *       notice, this list of conditions and the following disclaimer in the
   *       documentation and/or other materials provided with the distribution.</li>
   *   <li>The names of its contributors may not be used to endorse or promote
   *       products derived from this software without specific prior written
   *       permission.</li>
   * </ol></td></tr>
   *
   * <tr><td><strong>THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
   * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
   * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
   * DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS
   * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, * OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
   * OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
   * USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
   * DAMAGE.</strong></td></tr>
   * </table>
   */
class Mersenne_Twister extends Int_Random_Generator
{
	20160529L;

	/** Size of the bytes pool. */
	private static const int   N = 624;

	/** Period second parameter. */
	private static const int   M = 397;

	/** X * MATRIX_A for X = {0, 1}. */
	private static const std::vector<int> MAG01 = { 0x0, 0x9908b0df };

	/** Bytes pool. */
	private std::vector<int> mt;

	/** Current index in the bytes pool. */
	private int   mti;

	/**
	 * Creates a random number generator.
	 * <p>
	 * The instance is initialized using the current time plus the
	 * system identity hash code of this instance as the seed.
	 */
	public Mersenne_Twister()
	{
		mt = int[N];
		set_seed(System.current_time_millis() + System.identity_hash_code(this));
	}

	/**
	 * Creates a random number generator using a single int seed.
	 * @param seed the initial seed (32 bits integer)
	 */
	public Mersenne_Twister(const int& seed)
	{
		mt = int[N];
		set_seed(seed);
	}

	/**
	 * Creates a random number generator using an int array seed.
	 * @param seed the initial seed (32 bits integers array), if NULL
	 * the seed of the generator will be related to the current time
	 */
	public Mersenne_Twister(std::vector<int> seed)
	{
		mt = int[N];
		set_seed(seed);
	}

	/**
	 * Creates a random number generator using a single long seed.
	 * @param seed the initial seed (64 bits integer)
	 */
	public Mersenne_Twister(long seed)
	{
		mt = int[N];
		set_seed(seed);
	}

	/**
	 * Reinitialize the generator as if just built with the given int seed.
	 * <p>
	 * The state of the generator is exactly the same as a new
	 * generator built with the same seed.
	 *
	 * @param seed the initial seed (32 bits integer)
	 */
	 //override
	public void set_seed(const int& seed)
	{
		// we use a long masked by 0xffffffffL as a poor man unsigned int
		long long_m_t = seed & 0xffffffffL;
		// NB: unlike original C code, we are working with java longs, // the cast below makes masking unnecessary
		mt[0] = static_cast<int>(long_m_t;
		for (mti = 1; mti < N; ++mti)
		{
			// See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier.
			// initializer from the 2002-01-09 C version by Makoto Matsumoto
			long_m_t = (1812433253l * (long_m_t ^ (long_m_t >> 30)) + mti) & 0xffffffffL;
			mt[mti] = static_cast<int>(long_m_t;
		}

		clear_cache(); // Clear normal deviate cache
	}

	/**
	 * Reinitialize the generator as if just built with the given int array seed.
	 * <p>
	 * The state of the generator is exactly the same as a new
	 * generator built with the same seed.
	 *
	 * @param seed the initial seed (32 bits integers array), if NULL
	 * the seed of the generator will be the current system time plus the
	 * system identity hash code of this instance
	 */
	 //override
	public void set_seed(std::vector<int> seed)
	{
		if (seed == NULL)
		{
			set_seed(System.current_time_millis() + System.identity_hash_code(this));
			return;
		}

		set_seed(19650218);
		int i = 1;
		int j = 0;

		for (int k = std::max(N, seed.size()); k != 0; k--)
		{
			long l0 = (mt[i] & 0x7fffffffl) | ((mt[i] < 0) ? 0x80000000l : 0x0l);
			long l1 = (mt[i - 1] & 0x7fffffffl) | ((mt[i - 1] < 0) ? 0x80000000l : 0x0l);
			long l = (l0 ^ ((l1 ^ (l1 >> 30)) * 1664525l)) + seed[j] + j; // non linear
			mt[i] = static_cast<int>((l & 0xffffffffl);
			i++; j++;
			if (i >= N)
			{
				mt[0] = mt[N - 1];
				i = 1;
			}
			if (j >= seed.size())
			{
				j = 0;
			}
		}

		for (int k = N - 1; k != 0; k--)
		{
			long l0 = (mt[i] & 0x7fffffffl) | ((mt[i] < 0) ? 0x80000000l : 0x0l);
			long l1 = (mt[i - 1] & 0x7fffffffl) | ((mt[i - 1] < 0) ? 0x80000000l : 0x0l);
			long l = (l0 ^ ((l1 ^ (l1 >> 30)) * 1566083941l)) - i; // non linear
			mt[i] = static_cast<int>((l & 0xffffffffL);
			i++;
			if (i >= N)
			{
				mt[0] = mt[N - 1];
				i = 1;
			}
		}

		mt[0] = 0x80000000; // MSB is 1; assuring non-zero initial array

		clear_cache(); // Clear normal deviate cache
	}

	/** {@inherit_doc} */
	//override
	public int next_int()
	{
		int y;

		if (mti >= N) { // generate N words at one time
			int mt_next = mt[0];
			for (int k{}; k < N - M; ++k)
			{
				int mt_curr = mt_next;
				mt_next = mt[k + 1];
				y = (mt_curr & 0x80000000) | (mt_next & 0x7fffffff);
				mt[k] = mt[k + M] ^ (y >> > 1) ^ MAG01[y & 0x1];
			}
			for (int k = N - M; k < N - 1; ++k)
			{
				int mt_curr = mt_next;
				mt_next = mt[k + 1];
				y = (mt_curr & 0x80000000) | (mt_next & 0x7fffffff);
				mt[k] = mt[k + (M - N)] ^ (y >> > 1) ^ MAG01[y & 0x1];
			}
			y = (mt_next & 0x80000000) | (mt[0] & 0x7fffffff);
			mt[N - 1] = mt[M - 1] ^ (y >> > 1) ^ MAG01[y & 0x1];

			mti = 0;
		}

		y = mt[mti++];

		// tempering
		y ^= y >> > 11;
		y ^= (y << 7) & 0x9d2c5680;
		y ^= (y << 15) & 0xefc60000;
		y ^= y >> > 18;

		return y;
	}
}
