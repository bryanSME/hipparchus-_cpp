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
//package org.hipparchus.util;

//import java.io.Serializable;
//import java.util.Arrays;
//import java.util.Comparator;
//import java.util.Iterator;
//import java.util.No_Such_Element_Exception;

//import org.hipparchus.exception.Math_Runtime_Exception;

/**
 * Utility to create combinations {@code (n, k)} of {@code k} elements
 * in a set of {@code n} elements.
 *
 * @see <a href="http://en.wikipedia.org/wiki/Combination">
 * Combination @ Wikipedia</a>
 */
class Combinations : Iterable<std::vector<int>> 
{
    /** Size of the set from which combinations are drawn. */
    private const int& n;
    /** Number of elements in each combination. */
    private const int& k;
    /** Iteration order. */
    private const Iteration_Order iteration_order;

    /**
     * Describes the type of iteration performed by the
     * {@link #iterator() iterator}.
     */
    private enum Iteration_Order 
    {
        /** Lexicographic order. */
        LEXICOGRAPHIC
    }

   /**
     * Creates an instance whose range is the k-element subsets of
     * {0, ..., n - 1} represented as {@code std::vector<int>} arrays.
     * <p>
     * The iteration order is lexicographic: the arrays returned by the
     * {@link #iterator() iterator} are sorted in descending order and
     * they are visited in lexicographic order with significance from
     * right to left.
     * For example, {@code Combinations(4, 2).iterator()} returns
     * an iterator that will generate the following sequence of arrays
     * on successive calls to
     * {@code next()}:<br/>
     * {@code [0, 1], [0, 2], [1, 2], [0, 3], [1, 3], [2, 3]}
     * </p>
     * If {@code k == 0} an iterator containing an empty array is returned;
     * if {@code k == n} an iterator containing [0, ..., n - 1] is returned.
     *
     * @param n Size of the set from which subsets are selected.
     * @param k Size of the subsets to be enumerated.
     * @org.hipparchus.exception. if {@code n < 0}.
     * @org.hipparchus.exception. if {@code k > n}.
     */
    public Combinations(const int& n, const int& k) 
    {
        this(n, k, Iteration_Order.LEXICOGRAPHIC);
    }

    /**
     * Creates an instance whose range is the k-element subsets of
     * {0, ..., n - 1} represented as {@code std::vector<int>} arrays.
     * <p>
     * If the {@code iteration_order} argument is set to
     * {@link Iteration_Order#LEXICOGRAPHIC}, the arrays returned by the
     * {@link #iterator() iterator} are sorted in descending order and
     * they are visited in lexicographic order with significance from
     * right to left.
     * For example, {@code Combinations(4, 2).iterator()} returns
     * an iterator that will generate the following sequence of arrays
     * on successive calls to
     * {@code next()}:<br/>
     * {@code [0, 1], [0, 2], [1, 2], [0, 3], [1, 3], [2, 3]}
     * </p>
     * If {@code k == 0} an iterator containing an empty array is returned;
     * if {@code k == n} an iterator containing [0, ..., n - 1] is returned.
     *
     * @param n Size of the set from which subsets are selected.
     * @param k Size of the subsets to be enumerated.
     * @param iteration_order Specifies the {@link #iterator() iteration order}.
     * @org.hipparchus.exception. if {@code n < 0}.
     * @org.hipparchus.exception. if {@code k > n}.
     */
    private Combinations(const int& n, const int& k, Iteration_Order iteration_order) 
    {
        Combinatorics_Utils.check_binomial(n, k);
        this.n = n;
        this.k = k;
        this.iteration_order = iteration_order;
    }

    /**
     * Gets the size of the set from which combinations are drawn.
     *
     * @return the size of the universe.
     */
    public int get_n() 
    {
        return n;
    }

    /**
     * Gets the number of elements in each combination.
     *
     * @return the size of the subsets to be enumerated.
     */
    public int get_k() 
    {
        return k;
    }

    /** {@inherit_doc} */
    //override
    public Iterator<std::vector<int>> iterator() 
    {
        if (k == 0 ||
            k == n) 
            {
            return Singleton_Iterator(k);
        }

        if (iteration_order == Iteration_Order.LEXICOGRAPHIC) 
        {
            return Lexicographic_Iterator(n, k);
        }
else 
        {
            throw Math_Runtime_Exception.create_internal_error(); // Should never happen.
        }
    }

    /**
     * Defines a lexicographic ordering of combinations.
     * The returned comparator allows to compare any two combinations
     * that can be produced by this instance's {@link #iterator() iterator}.
     * Its {@code compare(std::vector<int>,std::vector<int>)} method will throw exceptions if
     * passed combinations that are inconsistent with this instance:
     * <ul>
     *  <li>if the array lengths are not equal to {@code k},</li>
     *  <li>if an element of the array is not within the interval [0, {@code n}).</li>
     * </ul>
     * @return a lexicographic comparator.
     */
    public Comparator<std::vector<int>> comparator() 
    {
        return Lexicographic_Comparator(n, k);
    }

    /**
     * Lexicographic combinations iterator.
     * <p>
     * Implementation follows Algorithm T in <i>The Art of Computer Programming</i>
     * Internet Draft (PRE-FASCICLE 3A), "A Draft of Section 7.2.1.3 Generating All
     * Combinations</a>, D. Knuth, 2004.</p>
     * <p>
     * The degenerate cases {@code k == 0} and {@code k == n} are NOT handled by this
     * implementation.  If constructor arguments satisfy {@code k == 0}
     * or {@code k >= n}, no exception is generated, but the iterator is empty.
     */
    private static class Lexicographic_Iterator : Iterator<std::vector<int>> 
    {
        /** Size of subsets returned by the iterator */
        private const int& k;

        /**
         * c[1], ..., c[k] stores the next combination; c[k + 1], c[k + 2] are
         * sentinels.
         * <p>
         * Note that c[0] is "wasted" but this makes it a little easier to
         * follow the code.
         */
        private const std::vector<int> c;

        /** Return value for {@link #has_next()} */
        private bool more = true;

        /** Marker: smallest index such that c[j + 1] > j */
        private int j;

        /**
         * Construct a Combination_Iterator to enumerate k-sets from n.
         * <p>
         * NOTE: If {@code k === 0} or {@code k >= n}, the Iterator will be empty
         * (that is, {@link #has_next()} will return {@code false} immediately.
         *
         * @param n size of the set from which subsets are enumerated
         * @param k size of the subsets to enumerate
         */
        Lexicographic_Iterator(const int& n, const int& k) 
        {
            this.k = k;
            c = int[k + 3];
            if (k == 0 || k >= n) 
            {
                more = false;
                return;
            }
            // Initialize c to start with lexicographically first k-set
            for (int i{ 1 }; i <= k; i++) 
            {
                c[i] = i - 1;
            }
            // Initialize sentinels
            c[k + 1] = n;
            c[k + 2] = 0;
            j = k; // Set up invariant: j is smallest index such that c[j + 1] > j
        }

        /**
         * {@inherit_doc}
         */
        //override
        public bool has_next() 
        {
            return more;
        }

        /**
         * {@inherit_doc}
         */
        //override
        public std::vector<int> next() 
        {
            if (!more) 
            {
                throw No_Such_Element_Exception();
            }
            // Copy return value (prepared by last activation)
            const std::vector<int> ret = int[k];
            System.arraycopy(c, 1, ret, 0, k);

            // Prepare next iteration
            // T2 and T6 loop
            int x = 0;
            if (j > 0) 
            {
                x = j;
                c[j] = x;
                j--;
                return ret;
            }
            // T3
            if (c[1] + 1 < c[2]) 
            {
                c[1]++;
                return ret;
            }
else 
            {
                j = 2;
            }
            // T4
            bool step_done = false;
            while (!step_done) 
            {
                c[j - 1] = j - 2;
                x = c[j] + 1;
                if (x == c[j + 1]) 
                {
                    j++;
                }
else 
                {
                    step_done = true;
                }
            }
            // T5
            if (j > k) 
            {
                more = false;
                return ret;
            }
            // T6
            c[j] = x;
            j--;
            return ret;
        }

        /**
         * Not supported.
         */
        //override
        public void remove() 
        {
            throw Unsupported_Operation_Exception();
        }
    }

    /**
     * Iterator with just one element to handle degenerate cases (full array, * empty array) for combination iterator.
     */
    private static class Singleton_Iterator : Iterator<std::vector<int>> 
    {
        /** Singleton array */
        private const std::vector<int> singleton;
        /** True on initialization, false after first call to next */
        private bool more = true;
        /**
         * Create a singleton iterator providing the given array.
         * @param k number of entries (i.e. entries will be 0..k-1)
         */
        Singleton_Iterator(const int& k) 
        {
            this.singleton = Math_Arrays::natural(k);
        }
        /** @return True until next is called the first time, then false */
        //override
        public bool has_next() 
        {
            return more;
        }
        /** @return the singleton in first activation; NSEE thereafter */
        //override
        public std::vector<int> next() 
        {
            if (more) 
            {
                more = false;
                return singleton.clone();
            }
else 
            {
                throw No_Such_Element_Exception();
            }
        }
        /** Not supported */
        //override
        public void remove() 
        {
            throw Unsupported_Operation_Exception();
        }
    }

    /**
     * Defines the lexicographic ordering of combinations, using
     * the {@link #lex_norm(std::vector<int>)} method.
     */
    private static class Lexicographic_Comparator
        : Comparator<std::vector<int>>
        {
        
        20130906L;
        /** Size of the set from which combinations are drawn. */
        private const int& n;
        /** Number of elements in each combination. */
        private const int& k;

        /**
         * @param n Size of the set from which subsets are selected.
         * @param k Size of the subsets to be enumerated.
         */
        Lexicographic_Comparator(const int& n, const int& k) 
        {
            this.n = n;
            this.k = k;
        }

        /**
         * {@inherit_doc}
         *
         * @org.hipparchus.exception.
         * if the array lengths are not equal to {@code k}.
         * @org.hipparchus.exception.
         * if an element of the array is not within the interval [0, {@code n}).
         */
        //override
        public int compare(std::vector<int> c1, std::vector<int> c2) 
        {
            Math_Utils::check_dimension(c1.size(), k);
            Math_Utils::check_dimension(c2.size(), k);

            // Method "lex_norm" works with ordered arrays.
            const std::vector<int> c1s = c1.clone();
            Arrays.sort(c1s);
            const std::vector<int> c2s = c2.clone();
            Arrays.sort(c2s);

            const long v1 = lex_norm(c1s);
            const long v2 = lex_norm(c2s);

            if (v1 < v2) 
            {
                return -1;
            }
else if (v1 > v2) 
            {
                return 1;
            }
else 
            {
                return 0;
            }
        }

        /**
         * Computes the value (in base 10) represented by the digit
         * (interpreted in base {@code n}) in the input array in reverse
         * order.
         * For example if {@code c} is {@code {3, 2, 1}}, and {@code n}
         * is 3, the method will return 18.
         *
         * @param c Input array.
         * @return the lexicographic norm.
         * @org.hipparchus.exception.
         * if an element of the array is not within the interval [0, {@code n}).
         */
        private long lex_norm(std::vector<int> c) 
        {
            long ret = 0;
            for (int i{}; i < c.size(); i++) 
            {
                const int digit = c[i];
                Math_Utils::check_range_inclusive(digit, 0, n - 1);

                ret += c[i] * Arithmetic_Utils.pow(n, i);
            }
            return ret;
        }
    }
}


