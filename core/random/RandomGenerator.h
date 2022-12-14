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
#include <vector>

/**
 * Interface for generators of random number sequences.
 */
class Random_Generator
{
	/**
	 * Sets the seed of the underlying random number generator using an
	 * <code>int</code> seed.
	 * <p>
	 * Sequences of values generated starting with the same seeds
	 * should be identical.
	 *
	 * @param seed the seed value
	 */
	virtual void set_seed(const int& seed) = 0;

	/**
	 * Sets the seed of the underlying random number generator using an
	 * <code>int</code> array seed.
	 * <p>
	 * Sequences of values generated starting with the same seeds
	 * should be identical.
	 *
	 * @param seed the seed value
	 */
	virtual void set_seed(std::vector<int> seed) = 0;

	/**
	 * Sets the seed of the underlying random number generator using a
	 * <code>long</code> seed.
	 * <p>
	 * Sequences of values generated starting with the same seeds
	 * should be identical.
	 *
	 * @param seed the seed value
	 */
	virtual void set_seed(long seed) = 0;

	/**
	 * Generates random bytes and places them into a user-supplied
	 * std::byte array. The number of random bytes produced is equal to
	 * the length of the std::byte array.
	 *
	 * @param bytes the non-null std::byte array in which to put the random bytes
	 */
	virtual void next_bytes(std::vector<std::byte>bytes) = 0;

	/**
	 * Generates random bytes and places them into a user-supplied
	 * std::byte array.
	 *
	 * @param bytes the non-null std::byte array in which to put the random bytes
	 * @param offset the starting index for inserting the generated bytes into
	 * the array
	 * @param len the number of bytes to generate
	 * @org.hipparchus.exception. if {@code offset < 0} or
	 * {@code offset + len >= bytes.size()}
	 */
	virtual void next_bytes(std::vector<std::byte>bytes, int offset, int len) = 0;

	/**
	 * Returns the next pseudorandom, uniformly distributed {@code int}
	 * value from this random number generator's sequence.
	 * <p>
	 * All 2<sup>32</sup> possible {@code int} values should be produced
	 * with (approximately) equal probability.
	 *
	 * @return the next pseudorandom, uniformly distributed {@code int}
	 * value from this random number generator's sequence
	 */
	virtual int next_int() = 0;

	/**
	 * Returns a pseudorandom, uniformly distributed {@code int} value
	 * between 0 (inclusive) and the specified value (exclusive), drawn from
	 * this random number generator's sequence.
	 *
	 * @param n the bound on the random number to be returned. Must be positive.
	 * @return  a pseudorandom, uniformly distributed {@code int}
	 * value between 0 (inclusive) and n (exclusive).
	 * @Illegal_Argument_Exception if n is not positive.
	 */
	virtual int next_int(const int& n) = 0;

	/**
	 * Returns the next pseudorandom, uniformly distributed {@code long}
	 * value from this random number generator's sequence. All 2<sup>64</sup>
	 * possible {@code long} values should be produced with (approximately)
	 * equal probability.
	 *
	 * @return the next pseudorandom, uniformly distributed {@code long}
	 * value from this random number generator's sequence
	 */
	virtual long next_long() = 0;

	/**
	 * Returns a pseudorandom, uniformly distributed {@code int} value
	 * between 0 (inclusive) and the specified value (exclusive), drawn from
	 * this random number generator's sequence.
	 *
	 * @param n the bound on the random number to be returned. Must be positive.
	 * @return  a pseudorandom, uniformly distributed {@code int}
	 * value between 0 (inclusive) and n (exclusive).
	 * @Illegal_Argument_Exception if n is not positive.
	 */
	virtual long next_long(const long& n) = 0;

	/**
	 * Returns the next pseudorandom, uniformly distributed
	 * {@code bool} value from this random number generator's sequence.
	 *
	 * @return  the next pseudorandom, uniformly distributed
	 * <code>bool</code> value from this random number generator's
	 * sequence
	 */
	virtual bool next_boolean() = 0;

	/**
	 * Returns the next pseudorandom, uniformly distributed {@code float}
	 * value between <code>0.0</code> and <code>1.0</code> from this random
	 * number generator's sequence.
	 *
	 * @return  the next pseudorandom, uniformly distributed {@code float}
	 * value between <code>0.0</code> and <code>1.0</code> from this
	 * random number generator's sequence
	 */
	virtual float next_float() = 0;

	/**
	 * Returns the next pseudorandom, uniformly distributed
	 * <code>double</code> value between <code>0.0</code> and
	 * <code>1.0</code> from this random number generator's sequence.
	 *
	 * @return  the next pseudorandom, uniformly distributed
	 *  <code>double</code> value between <code>0.0</code> and
	 *  <code>1.0</code> from this random number generator's sequence
	 */
	virtual double next_double() = 0;

	/**
	 * Returns the next pseudorandom, Gaussian ("normally") distributed
	 * <code>double</code> value with mean <code>0.0</code> and standard
	 * deviation <code>1.0</code> from this random number generator's sequence.
	 *
	 * @return  the next pseudorandom, Gaussian ("normally") distributed
	 * <code>double</code> value with mean <code>0.0</code> and
	 * standard deviation <code>1.0</code> from this random number
	 *  generator's sequence
	 */
	virtual double next_gaussian() = 0;
};