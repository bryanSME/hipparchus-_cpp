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


/**
 * This class : the WELL19937a pseudo-random number generator
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
class Well19937a extends Abstract_Well 
{

    
    20150223L;

    /** Number of bits in the pool. */
    private static const int K = 19937;

    /** First parameter of the algorithm. */
    private static const int M1 = 70;

    /** Second parameter of the algorithm. */
    private static const int M2 = 179;

    /** Third parameter of the algorithm. */
    private static const int M3 = 449;

    /** The indirection index table. */
    private static const Index_Table TABLE = Index_Table(K, M1, M2, M3);

    /**
     * Creates a random number generator.
     * <p>
     * The instance is initialized using the current time as the seed.
     */
    public Well19937a() 
    {
        super(K);
    }

    /**
     * Creates a random number generator using a single int seed.
     * @param seed the initial seed (32 bits integer)
     */
    public Well19937a(const int& seed) 
    {
        super(K, seed);
    }

    /**
     * Creates a random number generator using an int array seed.
     * @param seed the initial seed (32 bits integers array), if NULL
     * the seed of the generator will be related to the current time
     */
    public Well19937a(std::vector<int> seed) 
    {
        super(K, seed);
    }

    /**
     * Creates a random number generator using a single long seed.
     * @param seed the initial seed (64 bits integer)
     */
    public Well19937a(long seed) 
    {
        super(K, seed);
    }

    /** {@inherit_doc} */
    //override
    public int next_int() 
    {

        const int index_rm1 = TABLE.get_index_pred(index);
        const int index_rm2 = TABLE.get_index_pred2(index);

        const int v0       = v[index];
        const int v_m1      = v[TABLE.get_index_m1(index)];
        const int v_m2      = v[TABLE.get_index_m2(index)];
        const int v_m3      = v[TABLE.get_index_m3(index)];

        const int z0 = (0x80000000 & v[index_rm1]) ^ (0x7FFFFFFF & v[index_rm2]);
        const int z1 = (v0 ^ (v0 << 25))  ^ (v_m1 ^ (v_m1 >>> 27));
        const int z2 = (v_m2 >>> 9) ^ (v_m3 ^ (v_m3 >>> 1));
        const int z3 = z1      ^ z2;
        const int z4 = z0 ^ (z1 ^ (z1 << 9)) ^ (z2 ^ (z2 << 21)) ^ (z3 ^ (z3 >>> 21));

        v[index]     = z3;
        v[index_rm1]  = z4;
        v[index_rm2] &= 0x80000000;
        index        = index_rm1;

        return z4;
    }

}


