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

//import java.util.Iterator;
//import java.util.concurrent.atomic.Atomic_Reference;

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Runtime_Exception;
//import org.hipparchus.special.Gamma;
#include <vector>
#include <cmath>

/**
 * Combinatorial utilities.
 */
class Combinatorics_Utils 
{

    /** All long-representable factorials */
    static const std::vector<long> FACTORIALS = 
    {
                       1l,
                       1l,
                       2l,
                       6l,
                       24l,
                       120l, 
                       720l,          
                       5040l,             
                       40320l,
                       362880l,
                       3628800l,
                       39916800l,
                       479001600l,    
                       6227020800l,
                       87178291200l,
                       1307674368000l,  
                       20922789888000l,  
                       355687428096000l,
                       6402373705728000l,
                       121645100408832000l,
                       2432902008176640000l 
    };

    /** Stirling numbers of the second kind. */
    static const Atomic_Reference<long[][]> STIRLING_S2 = Atomic_Reference<> (null);

    /**
     * Default implementation of {@link #factorial_logstatic_cast<int>(} method:
     * <ul>
     *  <li>No pre-computation</li>
     *  <li>No cache allocation</li>
     * </ul>
     */
    private static const Factorial_Log FACTORIAL_LOG_NO_CACHE = Factorial_Log.create();

    /** Private constructor (class contains only static methods). */
    private Combinatorics_Utils() {}


    /**
     * Returns an exact representation of the <a
     * href="http://mathworld.wolfram.com/BinomialCoefficient.html"> Binomial
     * Coefficient</a>, "{@code n choose k}", the number of
     * {@code k}-element subsets that can be selected from an
     * {@code n}-element set.
     * <p>
     * <Strong>Preconditions</strong>:
     * <ul>
     * <li> {@code 0 <= k <= n } (otherwise
     * {@code } is thrown)</li>
     * <li> The result is small enough to fit into a {@code long}. The
     * largest value of {@code n} for which all coefficients are
     * {@code  < long.MAX_VALUE} is 66. If the computed value exceeds
     * {@code long.MAX_VALUE} a {@code Math_Runtime_Exception} is
     * thrown.</li>
     * </ul></p>
     *
     * @param n the size of the set
     * @param k the size of the subsets to be counted
     * @return {@code n choose k}
     * @ if {@code n < 0}.
     * @ if {@code k > n}.
     * @Math_Runtime_Exception if the result is too large to be
     * represented by a long integer.
     */
    public static long binomial_coefficient(const int& n, const int& k)
    {
        Combinatorics_Utils.check_binomial(n, k);
        if ((n == k) || (k == 0)) 
        {
            return 1;
        }
        if ((k == 1) || (k == n - 1)) 
        {
            return n;
        }
        // Use symmetry for large k
        if (k > n / 2) 
        {
            return binomial_coefficient(n, n - k);
        }

        // We use the formula
        // (n choose k) = n! / (n-k)! / k!
        // (n choose k) == ((n-k+1)*...*n) / (1*...*k)
        // which could be written
        // (n choose k) == (n-1 choose k-1) * n / k
        long result = 1;
        if (n <= 61) 
        {
            // For n <= 61, the naive implementation cannot overflow.
            int i = n - k + 1;
            for (int j{ 1 }; j <= k; j++) 
            {
                result = result * i / j;
                i++;
            }
        }
        else if (n <= 66) 
        {
            // For n > 61 but n <= 66, the result cannot overflow, // but we must take care not to overflow intermediate values.
            int i = n - k + 1;
            for (int j{ 1 }; j <= k; j++) 
            {
                // We know that (result * i) is divisible by j, // but (result * i) may overflow, so we split j:
                // Filter out the gcd, d, so j/d and i/d are integer.
                // result is divisible by (j/d) because (j/d)
                // is relative prime to (i/d) and is a divisor of
                // result * (i/d).
                const long d = Arithmetic_Utils.gcd(i, j);
                result = (result / (j / d)) * (i / d);
                i++;
            }
        }
        else 
        {
            // For n > 66, a result overflow might occur, so we check
            // the multiplication, taking care to not overflow
            // unnecessary.
            int i = n - k + 1;
            for (int j{ 1 }; j <= k; j++) 
            {
                const long d = Arithmetic_Utils.gcd(i, j);
                result = Arithmetic_Utils.mul_and_check(result / (j / d), i / d);
                i++;
            }
        }
        return result;
    }

    /**
     * Returns a {@code double} representation of the <a
     * href="http://mathworld.wolfram.com/BinomialCoefficient.html"> Binomial
     * Coefficient</a>, "{@code n choose k}", the number of
     * {@code k}-element subsets that can be selected from an
     * {@code n}-element set.
     * <p>
     * <Strong>Preconditions</strong>:
     * <ul>
     * <li> {@code 0 <= k <= n } (otherwise
     * {@code Illegal_Argument_Exception} is thrown)</li>
     * <li> The result is small enough to fit into a {@code double}. The
     * largest value of {@code n} for which all coefficients are &lt;
     * Double.MAX_VALUE is 1029. If the computed value exceeds Double.MAX_VALUE, * INFINITY is returned</li>
     * </ul></p>
     *
     * @param n the size of the set
     * @param k the size of the subsets to be counted
     * @return {@code n choose k}
     * @ if {@code n < 0}.
     * @ if {@code k > n}.
     * @Math_Runtime_Exception if the result is too large to be
     * represented by a long integer.
     */
    public static double binomial_coefficient_double(const int& n, const int& k)
    {
        Combinatorics_Utils.check_binomial(n, k);
        if ((n == k) || (k == 0)) 
        {
            return 1;
        }
        if ((k == 1) || (k == n - 1)) 
        {
            return n;
        }
        if (k > n/2) 
        {
            return binomial_coefficient_double(n, n - k);
        }
        if (n < 67) 
        {
            return binomial_coefficient(n,k);
        }

        double result = 1;
        for (int i{ 1 }; i <= k; i++) 
        {
             result *= static_cast<double>((n - k + i) / static_cast<double>(i;
        }

        return std::floor(result + 0.5);
    }

    /**
     * Returns the natural {@code log} of the <a
     * href="http://mathworld.wolfram.com/BinomialCoefficient.html"> Binomial
     * Coefficient</a>, "{@code n choose k}", the number of
     * {@code k}-element subsets that can be selected from an
     * {@code n}-element set.
     * <p>
     * <Strong>Preconditions</strong>:
     * <ul>
     * <li> {@code 0 <= k <= n } (otherwise
     * {@code } is thrown)</li>
     * </ul></p>
     *
     * @param n the size of the set
     * @param k the size of the subsets to be counted
     * @return {@code n choose k}
     * @ if {@code n < 0}.
     * @ if {@code k > n}.
     * @Math_Runtime_Exception if the result is too large to be
     * represented by a long integer.
     */
    public static double binomial_coefficient_log(const int& n, const int& k)
    {
        Combinatorics_Utils::check_binomial(n, k);
        if ((n == k) || (k == 0)) 
        {
            return 0;
        }
        if ((k == 1) || (k == n - 1)) 
        {
            return std::log(n);
        }

        /*
         * For values small enough to do exact integer computation, * return the log of the exact value
         */
        if (n < 67) 
        {
            return std::log(binomial_coefficient(n,k));
        }

        /*
         * Return the log of binomial_coefficient_double for values that will not
         * overflow binomial_coefficient_double
         */
        if (n < 1030) 
        {
            return std::log(binomial_coefficient_double(n, k));
        }

        if (k > n / 2) 
        {
            return binomial_coefficient_log(n, n - k);
        }

        /*
         * Sum logs for values that could overflow
         */
        double log_sum = 0;

        // n!/(n-k)!
        for (int i = n - k + 1; i <= n; i++) 
        {
            log_sum += std::log(i);
        }

        // divide by k!
        for (int i{ 2 }; i <= k; i++) 
        {
            log_sum -= std::log(i);
        }

        return log_sum;
    }

    /**
     * Returns n!. Shorthand for {@code n} <a
     * href="http://mathworld.wolfram.com/Factorial.html"> Factorial</a>, the
     * product of the numbers {@code 1,...,n}.
     * <p>
     * <Strong>Preconditions</strong>:
     * <ul>
     * <li> {@code n >= 0} (otherwise
     * {@code } is thrown)</li>
     * <li> The result is small enough to fit into a {@code long}. The
     * largest value of {@code n} for which {@code n!} does not exceed
     * long.MAX_VALUE} is 20. If the computed value exceeds {@code long.MAX_VALUE}
     * an {@code Math_Runtime_Exception } is thrown.</li>
     * </ul>
     * </p>
     *
     * @param n argument
     * @return {@code n!}
     * @Math_Runtime_Exception if the result is too large to be represented
     * by a {@code long}.
     * @ if {@code n < 0}.
     * @ if {@code n > 20}: The factorial value is too
     * large to fit in a {@code long}.
     */
    public static long factorial(const int& n)  
    {
        if (n < 0) 
        {
            throw std::exception("not implemented");
            //  throw (hipparchus::exception::Localized_Core_Formats_Type::FACTORIAL_NEGATIVE_PARAMETER, n);
        }
        if (n > 20) 
        {
            throw std::exception("not implemented");
            //   throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_LARGE, n, 20);
        }
        return FACTORIALS[n];
    }

    /**
     * Compute n!, the<a href="http://mathworld.wolfram.com/Factorial.html">
     * factorial</a> of {@code n} (the product of the numbers 1 to n), as a
     * {@code double}.
     * The result should be small enough to fit into a {@code double}: The
     * largest {@code n} for which {@code n!} does not exceed
     * {@code Double.MAX_VALUE} is 170. If the computed value exceeds
     * {@code Double.MAX_VALUE}, {@code INFINITY} is returned.
     *
     * @param n Argument.
     * @return {@code n!}
     * @ if {@code n < 0}.
     */
    public static double factorial_double(const int& n)  
    {
        if (n < 0) 
        {
            throw std::exception("not implemented");
            // throw (hipparchus::exception::Localized_Core_Formats_Type::FACTORIAL_NEGATIVE_PARAMETER, n);
        }
        if (n < 21) 
        {
            return FACTORIALS[n];
        }
        return std::floor(std::exp(Combinatorics_Utils.factorial_log(n)) + 0.5);
    }

    /**
     * Compute the natural logarithm of the factorial of {@code n}.
     *
     * @param n Argument.
     * @return {@code log(n!)}
     * @ if {@code n < 0}.
     */
    public static double factorial_log(const int& n)  
    {
        return FACTORIAL_LOG_NO_CACHE.value(n);
    }

    /**
     * Returns the <a
     * href="http://mathworld.wolfram.com/StirlingNumberoftheSecondKind.html">
     * Stirling number of the second kind</a>, "{@code S(n,k)}", the number of
     * ways of partitioning an {@code n}-element set into {@code k} non-empty
     * subsets.
     * <p>
     * The preconditions are {@code 0 <= k <= n } (otherwise
     * {@code } is thrown)
     * </p>
     * @param n the size of the set
     * @param k the number of non-empty subsets
     * @return {@code S(n,k)}
     * @ if {@code k < 0}.
     * @ if {@code k > n}.
     * @Math_Runtime_Exception if some overflow happens, typically for n exceeding 25 and
     * k between 20 and n-2 (S(n,n-1) is handled specifically and does not overflow)
     */
    public static long stirling_s2(const int& n, const int& k)
    {
        if (k < 0) 
        {
            throw std::exception("not implemented");
            //  throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL, k, 0);
        }
        if (k > n) 
        {
            throw std::exception("not implemented");
            //  throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_LARGE, k, n);
        }

        std::vector<std::vector<long>> stirling_s2 = STIRLING_S2.get();

        if (stirling_s2 == NULL) 
        {
            // the cache has never been initialized, compute the first numbers
            // by direct recurrence relation

            // as S(26,9) = 11201516780955125625 is larger than long.MAX_VALUE
            // we must stop computation at row 26
            const int max_index = 26;
            stirling_s2 = std::vector<std::vector<long>>(max_index);
            stirling_s2[0] = std::vector<long>() { 1l };
            for (int i{ 1 }; i < stirling_s2.size(); ++i) 
            {
                stirling_s2[i] = std::vector<long>(i + 1);
                stirling_s2[i][0] = 0;
                stirling_s2[i][1] = 1;
                stirling_s2[i][i] = 1;
                for (int j = 2; j < i; ++j) 
                {
                    stirling_s2[i][j] = j * stirling_s2[i - 1][j] + stirling_s2[i - 1][j - 1];
                }
            }

            // atomically save the cache
            STIRLING_S2.compare_and_set(null, stirling_s2);

        }

        if (n < stirling_s2.size()) 
        {
            // the number is in the small cache
            return stirling_s2[n][k];
        }
        else 
        {
            // use explicit formula to compute the number without caching it
            if (k == 0) 
            {
                return 0;
            }
            else if (k == 1 || k == n) 
            {
                return 1;
            }
            else if (k == 2) 
            {
                return (1l << (n - 1)) - 1l;
            }
            else if (k == n - 1) 
            {
                return binomial_coefficient(n, 2);
            }
            else 
            {
                // definition formula: note that this may trigger some overflow
                long sum = 0;
                long sign = ((k & 0x1) == 0) ? 1 : -1;
                for (int j{ 1 }; j <= k; ++j) 
                {
                    sign = -sign;
                    sum += sign * binomial_coefficient(k, j) * Arithmetic_Utils::pow(j, n);
                    if (sum < 0) 
                    {
                        // there was an overflow somewhere
                        throw std::exception("not implemented");
                        //throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::OUT_OF_RANGE_SIMPLE, n, 0, stirling_s2.size() - 1);
                    }
                }
                return sum / factorial(k);
            }
        }

    }

    /**
     * Returns an iterator whose range is the k-element subsets of {0, ..., n - 1}
     * represented as {@code std::vector<int>} arrays.
     * <p>
     * The arrays returned by the iterator are sorted in descending order and
     * they are visited in lexicographic order with significance from right to
     * left. For example, combinations_iterator(4, 2) returns an Iterator that
     * will generate the following sequence of arrays on successive calls to
     * {@code next()}:</p><p>
     * {@code [0, 1], [0, 2], [1, 2], [0, 3], [1, 3], [2, 3]}
     * </p><p>
     * If {@code k == 0} an Iterator containing an empty array is returned and
     * if {@code k == n} an Iterator containing [0, ..., n -1] is returned.</p>
     *
     * @param n Size of the set from which subsets are selected.
     * @param k Size of the subsets to be enumerated.
     * @return an {@link Iterator iterator} over the k-sets in n.
     * @ if {@code n < 0}.
     * @ if {@code k > n}.
     */
    public static Iterator<std::vector<int>> combinations_iterator(const int& n, const int& k) 
    {
        return Combinations(n, k).iterator();
    }

    /**
     * Check binomial preconditions.
     *
     * @param n Size of the set.
     * @param k Size of the subsets to be counted.
     * @ if {@code n < 0}.
     * @ if {@code k > n}.
     */
    public static void check_binomial(const int& n, const int& k)
    {
        if (n < k) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::BINOMIAL_INVALID_PARAMETERS_ORDER, k, n, true);
        }
        if (n < 0) 
        {
            throw std::exception("not implemented");
            // throw (hipparchus::exception::Localized_Core_Formats_Type::BINOMIAL_NEGATIVE_PARAMETER, n);
        }
    }

    /**
     * Class for computing the natural logarithm of the factorial of {@code n}.
     * It allows to allocate a cache of precomputed values.
     * In case of cache miss, computation is preformed by a call to
     * {@link Gamma#log_gammastatic_cast<double>(}.
     */
    public static const class Factorial_Log 
    {
    private:
        /**
         * Precomputed values of the function:
         * {@code LOG_FACTORIALS[i] = log(i!)}.
         */
        std::vector<double> LOG_FACTORIALS;

        /**
         * Creates an instance, reusing the already computed values if available.
         *
         * @param num_values Number of values of the function to compute.
         * @param cache Existing cache.
         * @throw  if {@code n < 0}.
         */
        Factorial_Log(const int& num_values, const std::vector<double>& cache) 
        {
            if (num_values < 0) 
            {
                throw std::exception("not implemented");
                //throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL, num_values, 0);
            }

            LOG_FACTORIALS = std::vector<double>(num_values];

            const int begin_copy{ 2 };
            const int end_copy = cache == NULL || cache.size() <= begin_copy ?
                begin_copy : cache.size() <= num_values ?
                cache.size() : num_values;

            // Copy available values.
            for (int i{ begin_copy }; i < end_copy; i++)
            {
                LOG_FACTORIALS[i] = cache[i];
            }

            // Precompute.
            for (int i{ end_copy }; i < num_values; i++)
            {
                LOG_FACTORIALS[i] = LOG_FACTORIALS[i - 1] + std::log(i);
            }
        }
    public:
        /**
         * Creates an instance with no precomputed values.
         * @return instance with no precomputed values
         */
        static Factorial_Log create() 
        {
            return Factorial_Log(0, NULL);
        }

        /**
         * Creates an instance with the specified cache size.
         *
         * @param cache_size Number of precomputed values of the function.
         * @return a instance where {@code cache_size} values have been
         * precomputed.
         * @ if {@code n < 0}.
         */
        Factorial_Log with_cache(const int& cache_size) 
        {
            return Factorial_Log(cache_size, LOG_FACTORIALS);
        }

        /**
         * Computes {@code log(n!)}.
         *
         * @param n Argument.
         * @return {@code log(n!)}.
         * @ if {@code n < 0}.
         */
        double value(const int& n) 
        {
            if (n < 0) 
            {
                throw std::exception("not implemented");
                // throw (hipparchus::exception::Localized_Core_Formats_Type::FACTORIAL_NEGATIVE_PARAMETER, n);
            }

            // Use cache of precomputed values.
            if (n < LOG_FACTORIALS.size()) 
            {
                return LOG_FACTORIALS[n];
            }

            // Use cache of precomputed factorial values.
            if (n < FACTORIALS.size()) 
            {
                return std::log(FACTORIALS[n]);
            }

            // Delegate.
            return Gamma::log_gamma(n + 1);
        }
    }
};