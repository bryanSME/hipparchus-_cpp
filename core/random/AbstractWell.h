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
#include <vector>

  /**
   * This virtual class : the WELL class of pseudo-random number generator
   * from Fran&ccedil;ois Panneton, Pierre L'Ecuyer and Makoto Matsumoto.
   * <p>
   * This generator is described in a paper by Fran&ccedil;ois Panneton, * Pierre L'Ecuyer and Makoto Matsumoto
   * <a href="http://www.iro.umontreal.ca/~lecuyer/myftp/papers/wellrng.pdf">
   * Improved long-_Period Generators Based on Linear Recurrences Modulo 2</a>
   * ACM Transactions on Mathematical Software, 32, 1 (2006). The errata for the paper
   * are in <a href="http://www.iro.umontreal.ca/~lecuyer/myftp/papers/wellrng-errata.txt">
   * wellrng-errata.txt</a>.
   *
   * @see <a href="http://www.iro.umontreal.ca/~panneton/WELLRNG.html">WELL Random number generator</a>
   */
class Abstract_Well : public Int_Random_Generator
{
protected:
	/** Current index in the bytes pool. */
	int index;

	/** Bytes pool. */
	const std::vector<int> v;

	/** Creates a random number generator.
	 * <p>The instance is initialized using the current time plus the
	 * system identity hash code of this instance as the seed.</p>
	 * @param k number of bits in the pool (not necessarily a multiple of 32)
	 */
	Abstract_Well(const int& k)
	{
		this(k, NULL);
	}

	/** Creates a random number generator using a single int seed.
	 * @param k number of bits in the pool (not necessarily a multiple of 32)
	 * @param seed the initial seed (32 bits integer)
	 */
	Abstract_Well(const int& k, const int& seed)
	{
		Abstract_Well(k, std::vector<int> { seed });
	}

	/**
	 * Creates a random number generator using an int array seed.
	 * @param k number of bits in the pool (not necessarily a multiple of 32)
	 * @param seed the initial seed (32 bits integers array), if NULL
	 * the seed of the generator will be related to the current time
	 */
	Abstract_Well(const int& k, const std::vector<int>& seed)
	{
		const int r = calculate_block_count(k);
		this.v = int[r];
		this.index = 0;

		// initialize the pool content
		set_seed(seed);
	}

	/**
	 * Creates a random number generator using a single long seed.
	 * @param k number of bits in the pool (not necessarily a multiple of 32)
	 * @param seed the initial seed (64 bits integer)
	 */
	Abstract_Well(const int& k, const long& seed)
	{
		Abstract_Well(k, std::vector<int> { static_cast<int>((seed >> > 32), static_cast<int>((seed & 0xffffffffl) });
	}

public:
	/**
	 * Reinitialize the generator as if just built with the given int array seed.
	 * <p>
	 * The state of the generator is exactly the same as a new
	 * generator built with the same seed.
	 *
	 * @param seed the initial seed (32 bits integers array). If NULL
	 * the seed of the generator will be the system time plus the system identity
	 * hash code of the instance.
	 */
	 //override
	void set_seed(const std::vector<int>& seed)
	{
		if (seed == NULL)
		{
			set_seed(System.current_time_millis() + System.identity_hash_code(this));
			return;
		}

		System.arraycopy(seed, 0, v, 0, std::min(seed.size(), v.size()));

		if (seed.size() < v.size())
		{
			for (int i = seed.size(); i < v.size(); ++i)
			{
				const long l = v[i - seed.size()];
				v[i] = static_cast<int>(((1812433253l * (l ^ (l >> 30)) + i) & 0xffffffffL);
			}
		}

		index = 0;
		clear_cache(); // Clear normal deviate cache
	}

private:
	/**
	 * Calculate the number of 32-bits blocks.
	 * @param k number of bits in the pool (not necessarily a multiple of 32)
	 * @return the number of 32-bits blocks
	 */
	static int calculate_block_count(const int& k)
	{
		// the bits pool contains k bits, k = r w - p where r is the number
		// of w bits blocks, w is the block size (always 32 in the original paper)
		// and p is the number of unused bits in the last block
		const int w = 32;
		return (k + w - 1) / w;
	}

	/**
	 * Inner class used to store the indirection index table which is fixed
	 * for a given type of WELL class of pseudo-random number generator.
	 */
	static const class Index_Table
	{
	private:
		/**
		 * Index indirection table giving for each index its predecessor
		 * taking table size into account.
		 */
		const std::vector<int> i_rm1;

		/**
		 * Index indirection table giving for each index its second predecessor
		 * taking table size into account.
		 */
		const std::vector<int> i_rm2;

		/**
		 * Index indirection table giving for each index the value index + m1
		 * taking table size into account.
		 */
		const std::vector<int> i1;

		/**
		 * Index indirection table giving for each index the value index + m2
		 * taking table size into account.
		 */
		const std::vector<int> i2;

		/**
		 * Index indirection table giving for each index the value index + m3
		 * taking table size into account.
		 */
		const std::vector<int> i3;

	public:
		/**
		 * Creates a pre-calculated indirection index table.
		 * @param k number of bits in the pool (not necessarily a multiple of 32)
		 * @param m1 first parameter of the algorithm
		 * @param m2 second parameter of the algorithm
		 * @param m3 third parameter of the algorithm
		 */
		Index_Table(const int& k, const int m1, const int m2, const int m3)
		{
			const int r = calculate_block_count(k);

			// precompute indirection index tables. These tables are used for optimizing access
			// they allow saving computations like "(j + r - 2) % r" with costly modulo operations
			i_rm1 = int[r];
			i_rm2 = int[r];
			i1 = int[r];
			i2 = int[r];
			i3 = int[r];
			for (int j{}; j < r; ++j)
			{
				i_rm1[j] = (j + r - 1) % r;
				i_rm2[j] = (j + r - 2) % r;
				i1[j] = (j + m1) % r;
				i2[j] = (j + m2) % r;
				i3[j] = (j + m3) % r;
			}
		}

		/**
		 * Returns the predecessor of the given index modulo the table size.
		 * @param index the index to look at
		 * @return (index - 1) % table size
		 */
		int get_index_pred(const int index)
		{
			return i_rm1[index];
		}

		/**
		 * Returns the second predecessor of the given index modulo the table size.
		 * @param index the index to look at
		 * @return (index - 2) % table size
		 */
		int get_index_pred2(const int index)
		{
			return i_rm2[index];
		}

		/**
		 * Returns index + M1 modulo the table size.
		 * @param index the index to look at
		 * @return (index + M1) % table size
		 */
		int get_index_m1(const int index)
		{
			return i1[index];
		}

		/**
		 * Returns index + M2 modulo the table size.
		 * @param index the index to look at
		 * @return (index + M2) % table size
		 */
		int get_index_m2(const int index)
		{
			return i2[index];
		}

		/**
		 * Returns index + M3 modulo the table size.
		 * @param index the index to look at
		 * @return (index + M3) % table size
		 */
		int get_index_m3(const int& index)
		{
			return i3[index];
		}
	};
};