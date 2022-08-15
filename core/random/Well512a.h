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
 * This class : the WELL512a pseudo-random number generator
 * from Fran&ccedil;ois Panneton, Pierre L'Ecuyer and Makoto Matsumoto.
 * <p>
 * This generator is described in a paper by Fran&ccedil;ois Panneton, * Pierre L'Ecuyer and Makoto Matsumoto <a
 * href="http://www.iro.umontreal.ca/~lecuyer/myftp/papers/wellrng.pdf">Improved
 * long-_Period Generators Based on Linear Recurrences Modulo 2</a> ACM
 * Transactions on Mathematical Software, 32, 1 (2006). The errata for the paper
 * are in <a href="http://www.iro.umontreal.ca/~lecuyer/myftp/papers/wellrng-errata.txt">
 * wellrng-errata.txt</a>.
 *
 * @see <a href="http://www.iro.umontreal.ca/~panneton/WELLRNG.html">WELL Random number generator</a>
 */
class Well512a : Abstract_Well
{
private:

    /** Number of bits in the pool. */
    static constexpr int my_k{ 512 };

    /** First parameter of the algorithm. */
    static constexpr int my_m1{ 13 };

    /** Second parameter of the algorithm. */
    static constexpr int my_m2{ 9 };

    /** Third parameter of the algorithm. */
    static constexpr int my_m3{ 5 };

    /** The indirection index table. */
    static const Index_Table TABLE = Index_Table(my_k, my_m1, my_m2, my_m3);

public:
    /**
     * Creates a random number generator.
     * <p>
     * The instance is initialized using the current time as the seed.
     */
    Well512a()
    {
        super(my_k);
    }

    /**
     * Creates a random number generator using a single int seed.
     * @param seed the initial seed (32 bits integer)
     */
    Well512a(const int& seed)
    {
        super(my_k, seed);
    }

    /**
     * Creates a random number generator using an int array seed.
     * @param seed the initial seed (32 bits integers array), if NULL
     * the seed of the generator will be related to the current time
     */
    Well512a(const std::vector<int>& seed)
    {
        super(my_k, seed);
    }

    /**
     * Creates a random number generator using a single long seed.
     * @param seed the initial seed (64 bits integer)
     */
    Well512a(const long& seed)
    {
        super(my_k, seed);
    }

    /** {@inherit_doc} */
    //override
    int next_int()
    {

        const int index_rm1 = TABLE.get_index_pred(index);

        const int vi = v[index];
        const int vi1 = v[TABLE.get_index_m1(index)];
        const int vi2 = v[TABLE.get_index_m2(index)];
        const int z0 = v[index_rm1];

        // the values below include the errata of the original article
        const int z1 = (vi ^ (vi << 16)) ^ (vi1 ^ (vi1 << 15));
        const int z2 = vi2 ^ (vi2 >> > 11);
        const int z3 = z1 ^ z2;
        const int z4 = (z0 ^ (z0 << 2)) ^ (z1 ^ (z1 << 18)) ^ (z2 << 28) ^ (z3 ^ ((z3 << 5) & 0xda442d24));

        v[index] = z3;
        v[index_rm1] = z4;
        index = index_rm1;

        return z4;
    }

};