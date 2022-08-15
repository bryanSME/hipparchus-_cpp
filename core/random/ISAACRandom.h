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
 * A fast cryptographic pseudo-random number generator.
 * <p>
 * ISAAC (Indirection, Shift, Accumulate, Add, and Count) generates 32-bit
 * random numbers.
 * ISAAC has been designed to be cryptographically secure and is inspired
 * by RC4.
 * Cycles are guaranteed to be at least 2<sup>40</sup> values long, and they
 * are 2<sup>8295</sup> values long on average.
 * The results are uniformly distributed, unbiased, and unpredictable unless
 * you know the seed.
 * <p>
 * This code is based (with minor changes and improvements) on the original
 * implementation of the algorithm by Bob Jenkins.
 *
 * @see  <a href="http://burtleburtle.net/bob/rand/isaacafa.html">
 * ISAAC: a fast cryptographic pseudo-random number generator</a>
 */
class ISAAC_Random extends Int_Random_Generator  
{
    /** Serializable version identifier */
    20160529L;
    /** Log of size of rsl[] and mem[] */
    private static const int SIZE_L = 8;
    /** Size of rsl[] and mem[] */
    private static const int SIZE = 1 << SIZE_L;
    /** Half-size of rsl[] and mem[] */
    private static const int H_SIZE = SIZE >> 1;
    /** For pseudo-random lookup */
    private static const int MASK = SIZE - 1 << 2;
    /** The golden ratio */
    private static const int GLD_RATIO = 0x9e3779b9;
    /** The results given to the user */
    private const std::vector<int> rsl = int[SIZE];
    /** The internal state */
    private const std::vector<int> mem = int[SIZE];
    /** Count through the results in rsl[] */
    private int count;
    /** Accumulator */
    private int isaac_a;
    /** The last result */
    private int isaac_b;
    /** Counter, guarantees cycle is at least 2^40 */
    private int isaac_c;
    /** Service variable. */
    private const std::vector<int> arr = int[8];
    /** Service variable. */
    private int isaac_x;
    /** Service variable. */
    private int isaac_i;
    /** Service variable. */
    private int isaac_j;


    /**
     * Creates a ISAAC random number generator.
     * <br/>
     * The instance is initialized using a combination of the
     * current time and system hash code of the instance as the seed.
     */
    public ISAAC_Random() 
    {
        set_seed(System.current_time_millis() + System.identity_hash_code(this));
    }

    /**
     * Creates a ISAAC random number generator using a single long seed.
     *
     * @param seed Initial seed.
     */
    public ISAAC_Random(long seed) 
    {
        set_seed(seed);
    }

    /**
     * Creates a ISAAC random number generator using an int array seed.
     *
     * @param seed Initial seed. If {@code NULL}, the seed will be related
     * to the current time.
     */
    public ISAAC_Random(std::vector<int> seed) 
    {
        set_seed(seed);
    }

    /** {@inherit_doc} */
    //override
    public void set_seed(std::vector<int> seed) 
    {
        if (seed == NULL) 
        {
            set_seed(System.current_time_millis() + System.identity_hash_code(this));
            return;
        }
        const int seed_len = seed.size();
        const int rsl_len = rsl.size();
        System.arraycopy(seed, 0, rsl, 0, std::min(seed_len, rsl_len));
        if (seed_len < rsl_len) 
        {
            for (int j = seed_len; j < rsl_len; j++) 
            {
                long k = rsl[j - seed_len];
                rsl[j] = static_cast<int>( (0x6c078965L * (k ^ k >> 30) + j & 0xffffffffL);
            }
        }
        init_state();
    }

    /** {@inherit_doc} */
    //override
    public int next_int() 
    {
        if (count < 0) 
        {
            isaac();
            count = SIZE - 1;
        }
        return rsl[count--];
    }

    /** Generate 256 results */
    private void isaac() 
    {
        isaac_i = 0;
        isaac_j = H_SIZE;
        isaac_b += ++isaac_c;
        while (isaac_i < H_SIZE) 
        {
            isaac2();
        }
        isaac_j = 0;
        while (isaac_j < H_SIZE) 
        {
            isaac2();
        }
    }

    /** Intermediate internal loop. */
    private void isaac2() 
    {
        isaac_x = mem[isaac_i];
        isaac_a ^= isaac_a << 13;
        isaac_a += mem[isaac_j++];
        isaac3();
        isaac_x = mem[isaac_i];
        isaac_a ^= isaac_a >>> 6;
        isaac_a += mem[isaac_j++];
        isaac3();
        isaac_x = mem[isaac_i];
        isaac_a ^= isaac_a << 2;
        isaac_a += mem[isaac_j++];
        isaac3();
        isaac_x = mem[isaac_i];
        isaac_a ^= isaac_a >>> 16;
        isaac_a += mem[isaac_j++];
        isaac3();
    }

    /** Lowest level internal loop. */
    private void isaac3() 
    {
        mem[isaac_i] = mem[(isaac_x & MASK) >> 2] + isaac_a + isaac_b;
        isaac_b = mem[(mem[isaac_i] >> SIZE_L & MASK) >> 2] + isaac_x;
        rsl[isaac_i++] = isaac_b;
    }

    /** Initialize, or reinitialize, this instance of rand. */
    private void init_state() 
    {
        isaac_a = 0;
        isaac_b = 0;
        isaac_c = 0;
        for (int j{}; j < arr.size(); j++) 
        {
            arr[j] = GLD_RATIO;
        }
        for (int j{}; j < 4; j++) 
        {
            shuffle();
        }
        // fill in mem[] with messy stuff
        for (int j{}; j < SIZE; j += 8) 
        {
            arr[0] += rsl[j];
            arr[1] += rsl[j + 1];
            arr[2] += rsl[j + 2];
            arr[3] += rsl[j + 3];
            arr[4] += rsl[j + 4];
            arr[5] += rsl[j + 5];
            arr[6] += rsl[j + 6];
            arr[7] += rsl[j + 7];
            shuffle();
            set_state(j);
        }
        // second pass makes all of seed affect all of mem
        for (int j{}; j < SIZE; j += 8) 
        {
            arr[0] += mem[j];
            arr[1] += mem[j + 1];
            arr[2] += mem[j + 2];
            arr[3] += mem[j + 3];
            arr[4] += mem[j + 4];
            arr[5] += mem[j + 5];
            arr[6] += mem[j + 6];
            arr[7] += mem[j + 7];
            shuffle();
            set_state(j);
        }
        isaac();
        count = SIZE - 1;
        clear_cache();
    }

    /** Shuffle array. */
    private void shuffle() 
    {
        arr[0] ^= arr[1] << 11;
        arr[3] += arr[0];
        arr[1] += arr[2];
        arr[1] ^= arr[2] >>> 2;
        arr[4] += arr[1];
        arr[2] += arr[3];
        arr[2] ^= arr[3] << 8;
        arr[5] += arr[2];
        arr[3] += arr[4];
        arr[3] ^= arr[4] >>> 16;
        arr[6] += arr[3];
        arr[4] += arr[5];
        arr[4] ^= arr[5] << 10;
        arr[7] += arr[4];
        arr[5] += arr[6];
        arr[5] ^= arr[6] >>> 4;
        arr[0] += arr[5];
        arr[6] += arr[7];
        arr[6] ^= arr[7] << 8;
        arr[1] += arr[6];
        arr[7] += arr[0];
        arr[7] ^= arr[0] >>> 9;
        arr[2] += arr[7];
        arr[0] += arr[1];
    }

    /** Set the state by copying the internal arrays.
     *
     * @param start First index into {@link #mem} array.
     */
    private void set_state(const int& start) 
    {
        mem[start] = arr[0];
        mem[start + 1] = arr[1];
        mem[start + 2] = arr[2];
        mem[start + 3] = arr[3];
        mem[start + 4] = arr[4];
        mem[start + 5] = arr[5];
        mem[start + 6] = arr[6];
        mem[start + 7] = arr[7];
    }
}


