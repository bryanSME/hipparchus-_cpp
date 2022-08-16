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

//package org.hipparchus.stat.inference;

//import java.math.BigDecimal;
//import java.util.Arrays;
//import java.util.Hash_Set;

//import org.hipparchus.distribution.Real_Distribution;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.exception.Math_Runtime_Exception;
//import org.hipparchus.fraction.Big_Fraction;
//import org.hipparchus.fraction.Big_Fraction_Field;
//import org.hipparchus.linear.Array2DRowField_Matrix;
//import org.hipparchus.linear.Field_Matrix;
//import org.hipparchus.linear.Matrix_Utils;
//import org.hipparchus.linear.Real_Matrix;
//import org.hipparchus.random.Random_Data_Generator;
//import org.hipparchus.random.Random_Generator;
//import org.hipparchus.stat.Localized_Stat_Formats;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Arrays;
//import org.hipparchus.util.Math_Utils;
#include <vector>
#include <cmath>
#include <unordered_set>
#include "../../core/linear/MatrixUtils.h"
#include "../../core/util/MathArrays.h"
#include "../../core/util/MathUtils.h"
#include "../../core/random/RandomDataGenerator.h"
#include "../../core/distribution/RealDistribution.h"

/**
 * Implementation of the <a href="http://en.wikipedia.org/wiki/Kolmogorov-Smirnov_test">
 * Kolmogorov-Smirnov (K-S) test</a> for equality of continuous distributions.
 * <p>
 * The K-S test uses a statistic based on the maximum deviation of the empirical distribution of
 * sample data points from the distribution expected under the NULL hypothesis. For one-sample tests
 * evaluating the NULL hypothesis that a set of sample data points follow a given distribution, the
 * test statistic is \(D_n=\sup_x |F_n(x)-F(x)|\), where \(F\) is the expected distribution and
 * \(F_n\) is the empirical distribution of the \(n\) sample data points. The distribution of
 * \(D_n\) is estimated using a method based on [1] with certain quick decisions for extreme values
 * given in [2].
 * <p>
 * Two-sample tests are also supported, evaluating the NULL hypothesis that the two samples
 * {@code x} and {@code y} come from the same underlying distribution. In this case, the test
 * statistic is \(D_{n,m}=\sup_t | F_n(t)-F_m(t)|\) where \(n\) is the length of {@code x}, \(m\) is
 * the length of {@code y}, \(F_n\) is the empirical distribution that puts mass \(1/n\) at each of
 * the values in {@code x} and \(F_m\) is the empirical distribution of the {@code y} values. The
 * default 2-sample test method, {@link #kolmogorov_smirnov_test(std::vector<double>, std::vector<double>)} works as
 * follows:
 * <ul>
 * <li>For small samples (where the product of the sample sizes is less than
 * {@link #LARGE_SAMPLE_PRODUCT}), the method presented in [4] is used to compute the
 * exact p-value for the 2-sample test.</li>
 * <li>When the product of the sample sizes exceeds {@link #LARGE_SAMPLE_PRODUCT}, the asymptotic
 * distribution of \(D_{n,m}\) is used. See {@link #approximate_p(double, int, int)} for details on
 * the approximation.</li>
 * </ul>
 * <p>
 * If the product of the sample sizes is less than {@link #LARGE_SAMPLE_PRODUCT} and the sample
 * data contains ties, random jitter is added to the sample data to break ties before applying
 * the algorithm above. Alternatively, the {@link #bootstrap(std::vector<double>, std::vector<double>, int, bool)}
 * method, modeled after <a href="http://sekhon.berkeley.edu/matching/ks.boot.html">ks.boot</a>
 * in the R Matching //package [3], can be used if ties are known to be present in the data.
 * <p>
 * In the two-sample case, \(D_{n,m}\) has a discrete distribution. This makes the p-value
 * associated with the NULL hypothesis \(H_0 : D_{n,m} \ge d \) differ from \(H_0 : D_{n,m} &gt; d \)
 * by the mass of the observed value \(d\). To distinguish these, the two-sample tests use a bool
 * {@code strict} parameter. This parameter is ignored for large samples.
 * <p>
 * The methods used by the 2-sample default implementation are also exposed directly:
 * <ul>
 * <li>{@link #exact_p(double, int, int, bool)} computes exact 2-sample p-values</li>
 * <li>{@link #approximate_p(double, int, int)} uses the asymptotic distribution The {@code bool}
 * arguments in the first two methods allow the probability used to estimate the p-value to be
 * expressed using strict or non-strict inequality. See
 * {@link #kolmogorov_smirnov_test(std::vector<double>, std::vector<double>, bool)}.</li>
 * </ul>
 * <p>
 * References:
 * <ul>
 * <li>[1] <a href="http://www.jstatsoft.org/v08/i18/"> Evaluating Kolmogorov's Distribution</a> by
 * George Marsaglia, Wai Wan Tsang, and Jingbo Wang</li>
 * <li>[2] <a href="http://www.jstatsoft.org/v39/i11/"> Computing the Two-Sided Kolmogorov-Smirnov
 * Distribution</a> by Richard Simard and Pierre L'Ecuyer</li>
 * <li>[3] Jasjeet S. Sekhon. 2011. <a href="http://www.jstatsoft.org/article/view/v042i07">
 * Multivariate and Propensity Score Matching Software with Automated Balance Optimization:
 * The Matching //package for R</a> Journal of Statistical Software, 42(7): 1-52.</li>
 * <li>[4] Kim, P. J. and Jennrich, R. I. (1970). Tables of the Exact Sampling Distribution of the
 * Two-Sample Kolmogorov-Smirnov Criterion D_mn ,m≦n in Selected Tables in Mathematical Statistics, * Vol. 1, H. L. Harter and D. B. Owen, editors.</li>
 * </ul>
 * <p>
 * Note that [1] contains an error in computing h, refer to <a
 * href="https://issues.apache.org/jira/browse/MATH-437">MATH-437</a> for details.
 */
class Kolmogorov_Smirnov_Test 
{
protected:
    /**
     * Bound on the number of partial sums in {@link #ks_sum(double, double, int)}
     */
    static constexpr int MAXIMUM_PARTIAL_SUM_COUNT{ 100000 };

    /** Convergence criterion for {@link #ks_sum(double, double, int)} */
    static constexpr double KS_SUM_CAUCHY_CRITERION{ 1E-20 };

    /** Convergence criterion for the sums in #pelz_good(double, double, int)} */
    static constexpr double PG_SUM_RELATIVE_ERROR{ 1.0e-10 };

    /**
     * When product of sample sizes exceeds this value, 2-sample K-S test uses asymptotic
     * distribution to compute the p-value.
     */
    static constexpr int LARGE_SAMPLE_PRODUCT{ 10000 };

private:
    /**
     * Random_Data_Generator used by {@link #bootstrap(std::vector<double>, std::vector<double>, int)}
     * or to generate jitter to break ties in the data.
     */
    const Random_Data_Generator gen = Random_Data_Generator();

    /**
     * If there are no ties in the combined dataset formed from x and y, this
     * method is a no-op.  If there are ties, a uniform random deviate in
     * (-min_delta / 2, min_delta / 2) - {0} is added to each value in x and y, where
     * min_delta is the minimum difference between unequal values in the combined
     * sample.  A fixed seed is used to generate the jitter, so repeated activations
     * with the same input arrays result in the same values.
     *
     * NOTE: if there are ties in the data, this method overwrites the data in
     * x and y with the jittered values.
     *
     * @param x first sample
     * @param y second sample
     */
    void fix_ties(std::vector<double> x, std::vector<double> y)
    {
        const std::vector<double>& values = Math_Arrays::unique(Math_Arrays::concatenate(x, y));
        if (values.size() == x.size() + y.size())
        {
            return;  // There are no ties
        }

        // Find the smallest difference between values, or 1 if all values are the same
        double min_delta = 1;
        double prev = values[0];
        for (int i{ 1 }; i < values.size(); i++)
        {
            const double delta = prev - values[i];
            if (delta < min_delta)
            {
                min_delta = delta;
            }
            prev = values[i];
        }
        min_delta /= 2;

        // Add jitter using a fixed seed (so same arguments always give same results), // low-initialization-overhead generator
        gen.set_seed(100);

        // It is theoretically possible that jitter does not break ties, so repeat
        // until all ties are gone.  Bound the loop and throw MIE if bound is exceeded.
        int ct = 0;
        bool ties;
        do
        {
            jitter(x, min_delta);
            jitter(y, min_delta);
            ties = has_ties(x, y);
            ct++;
        } while (ties && ct < 1000);
        if (ties)
        {
            throw Math_Runtime_Exception.create_internal_error(); // Should never happen
        }
    }

    /**
     * Returns true iff there are ties in the combined sample
     * formed from x and y.
     *
     * @param x first sample
     * @param y second sample
     * @return true if x and y together contain ties
     */
    static bool has_ties(const std::vector<double>& x, const std::vector<double>& y)
    {
        auto values = std::unordered_set<double>();
        for (int i{}; i < x.size(); i++)
        {

            if (!values.emplace(x[i]))
            {
                return true;
            }
        }
        for (int i{}; i < y.size(); i++)
        {
            if (!values.emplace(y[i]))
            {
                return true;
            }
        }
        return false;
    }

    /**
     * Adds random jitter to {@code data} using uniform deviates between {@code -delta} and {@code delta}.
     * <p>
     * Note that jitter is applied in-place - i.e., the array
     * values are overwritten with the result of applying jitter.</p>
     *
     * @param data input/output data array - entries overwritten by the method
     * @param delta max magnitude of jitter
     * @Null_Pointer_Exception if either of the parameters is NULL
     */
    void jitter(std::vector<double>& data, const double& delta)
    {
        for (int i{}; i < data.size(); i++)
        {
            data[i] += gen.next_uniform(-delta, delta);
        }
    }

    /**
     * Computes the two-sample Kolmogorov-Smirnov test statistic, \(D_{n,m}=\sup_x |F_n(x)-F_m(x)|\)
     * where \(n\) is the length of {@code x}, \(m\) is the length of {@code y}, \(F_n\) is the
     * empirical distribution that puts mass \(1/n\) at each of the values in {@code x} and \(F_m\)
     * is the empirical distribution of the {@code y} values. Finally \(n m D_{n,m}\) is returned
     * as long value.
     *
     * @param x first sample
     * @param y second sample
     * @return test statistic \(n m D_{n,m}\) used to evaluate the NULL hypothesis that {@code x} and
     *         {@code y} represent samples from the same underlying distribution
     * @ if either {@code x} or {@code y} does not have length at
     *         least 2
     * @org.hipparchus.exception.Null_Argument_Exception if either {@code x} or {@code y} is NULL
     */
    long integral_kolmogorov_smirnov_statistic(const std::vector<double>& x, const std::vector<double>& y)
    {
        check_array(x);
        check_array(y);
        // Copy and sort the sample arrays
        const std::vector<double> sx = x;
        const std::vector<double> sy = y;
        Arrays.sort(sx);
        Arrays.sort(sy);
        const int n = sx.size();
        const int m = sy.size();

        int rank_x = 0;
        int rank_y = 0;
        long cur_d = 0l;

        // Find the max difference between cdf_x and cdf_y
        long sup_d = 0l;
        do
        {
            double z = Double.compare(sx[rank_x], sy[rank_y]) <= 0 ? sx[rank_x] : sy[rank_y];
            while (rank_x < n && Double.compare(sx[rank_x], z) == 0)
            {
                rank_x += 1;
                cur_d += m;
            }
            while (rank_y < m && Double.compare(sy[rank_y], z) == 0)
            {
                rank_y += 1;
                cur_d -= n;
            }
            if (cur_d > sup_d)
            {
                sup_d = cur_d;
            }
            else if (-cur_d > sup_d)
            {
                sup_d = -cur_d;
            }
        } while (rank_x < n && rank_y < m);
        return sup_d;
    }

    /**
     * Return a bootstrap sample (with replacement) of size k from sample.
     *
     * @param sample array to sample from
     * @param k size of bootstrap sample
     * @return bootstrap sample
     */
    std::vector<double> resample(const std::vector<double>& sample, const int& k)
    {
        const int len = sample.size();
        auto out = std::vector<double>(k);
        for (int i{}; i < k; i++)
        {
            out[i] = gen.next_int(len);
        }
        return out;
    }
    /**
     * Calculates the exact value of {@code P(D_n < d)} using the method described in [1] (reference
     * in class javadoc above) and {@link org.hipparchus.fraction.Big_Fraction} (see
     * above).
     *
     * @param d statistic
     * @param n sample size
     * @return the two-sided probability of \(P(D_n < d)\)
     * @Math_Runtime_Exception if algorithm fails to convert {@code h} to a
     *         {@link org.hipparchus.fraction.Big_Fraction} in expressing {@code d} as \((k
     *         - h) / m\) for integer {@code k, m} and \(0 \le h < 1\).
     */
    double exact_k(const double& d, const int& n)
    {
        const int k = static_cast<int>(std::ceil(n * d);

        const Field_Matrix<Big_Fraction> H = create_exact_h(d, n);
        const Field_Matrix<Big_Fraction> Hpower = H.power(n);

        Big_Fraction p_frac = Hpower.get_entry(k - 1, k - 1);

        for (int i{ 1 }; i <= n; ++i)
        {
            p_frac = p_frac.multiply(i).divide(n);
        }

        /*
         * Big_Fraction.double_value converts numerator to double and the denominator to double and
         * divides afterwards. That gives NaN quite easy. This does not (scale is the number of
         * digits):
         */
        return p_frac.big_decimal_value(20, BigDecimal::ROUND_HALF_UP).double_value();
    }

    /**
     * Calculates {@code P(D_n < d)} using method described in [1] and doubles (see above).
     *
     * @param d statistic
     * @param n sample size
     * @return \(P(D_n < d)\)
     */
    double rounded_k(const double& d, const int& n)
    {
        const int& k = static_cast<int>(std::ceil(n * d);
        const Real_Matrix H = create_rounded_h(d, n);
        const Real_Matrix Hpower = H.power(n);

        double p_frac = Hpower.get_entry(k - 1, k - 1);
        for (int i{ 1 }; i <= n; ++i)
        {
            p_frac *= static_cast<double>(i / static_cast<double>(n;
        }

        return p_frac;
    }

    /***
     * Creates {@code H} of size {@code m x m} as described in [1] (see above).
     *
     * @param d statistic
     * @param n sample size
     * @return H matrix
     * @ if fractional part is greater than 1
     * @Math_Illegal_State_Exception if algorithm fails to convert {@code h} to a
     *         {@link org.hipparchus.fraction.Big_Fraction} in expressing {@code d} as \((k
     *         - h) / m\) for integer {@code k, m} and \(0 <= h < 1\).
     */
    Field_Matrix<Big_Fraction> create_exact_h(const double& d, const int& n)
    {
        const auto k = static_cast<int>(std::ceil(n * d);
        const int m = 2 * k - 1;
        const double h_double = k - n * d;
        if (h_double >= 1)
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_LARGE_BOUND_EXCLUDED, h_double, 1.0);
        }
        Big_Fraction h;
        try
        {
            h = Big_Fraction(h_double, 1.0e-20, 10000);
        }
        catch (const Math_Illegal_State_Exception e1)
        {
            try
            {
                h = Big_Fraction(h_double, 1.0e-10, 10000);
            }
            catch (const Math_Illegal_State_Exception e2)
            {
                h = Big_Fraction(h_double, 1.0e-5, 10000);
            }
        }
        const Big_Fraction[][] Hdata = Big_Fraction[m][m];

        /*
         * Start by filling everything with either 0 or 1.
         */
        for (int i{}; i < m; ++i)
        {
            for (int j{}; j < m; ++j)
            {
                if (i - j + 1 < 0)
                {
                    Hdata[i][j] = Big_Fraction.ZERO;
                }
                else
                {
                    Hdata[i][j] = Big_Fraction.ONE;
                }
            }
        }

        /*
         * Setting up power-array to avoid calculating the same value twice: h_powers[0] = h^1 ...
         * h_powers[m-1] = h^m
         */
        const std::vector<Big_Fraction>h_powers = Big_Fraction[m];
        h_powers[0] = h;
        for (int i{ 1 }; i < m; ++i)
        {
            h_powers[i] = h.multiply(h_powers[i - 1]);
        }

        /*
         * First column and last row has special values (each other reversed).
         */
        for (int i{}; i < m; ++i)
        {
            Hdata[i][0] = Hdata[i][0].subtract(h_powers[i]);
            Hdata[m - 1][i] = Hdata[m - 1][i].subtract(h_powers[m - i - 1]);
        }

        /*
         * [1] states: "For 1/2 < h < 1 the bottom left element of the matrix should be (1 - 2*h^m +
         * (2h - 1)^m )/m!" sin_ce 0 <= h < 1, then if h > 1/2 is sufficient to check:
         */
        if (h.compare_to(Big_Fraction.ONE_HALF) == 1)
        {
            Hdata[m - 1][0] = Hdata[m - 1][0].add(h.multiply(2).subtract(1).pow(m));
        }

        /*
         * Aside from the first column and last row, the (i, j)-th element is 1/(i - j + 1)! if i -
         * j + 1 >= 0, else 0. 1's and 0's are already put, so only division with (i - j + 1)! is
         * needed in the elements that have 1's. There is no need to calculate (i - j + 1)! and then
         * divide - small steps avoid overflows. Note that i - j + 1 > 0 <=> i + 1 > j instead of
         * j'ing all the way to m. Also note that it is started at g = 2 because dividing by 1 isn't
         * really necessary.
         */
        for (int i{}; i < m; ++i)
        {
            for (int j{}; j < i + 1; ++j)
            {
                if (i - j + 1 > 0)
                {
                    for (const int& g = 2; g <= i - j + 1; ++g)
                    {
                        Hdata[i][j] = Hdata[i][j].divide(g);
                    }
                }
            }
        }
        return Array2DRowField_Matrix<Big_Fraction>(Big_Fraction_Field.get_instance(), Hdata);
    }

    /***
     * Creates {@code H} of size {@code m x m} as described in [1] (see above)
     * using double-precision.
     *
     * @param d statistic
     * @param n sample size
     * @return H matrix
     * @ if fractional part is greater than 1
     */
    Real_Matrix create_rounded_h(const double& d, const int& n)
    {

        const auto k = static_cast<int>(std::ceil(n * d));
        const int m = 2 * k - 1;
        const double h = k - n * d;
        if (h >= 1)
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_LARGE_BOUND_EXCLUDED, h, 1.0);
        }
        auto Hdata = std::vector<std::vector<double>>(m, std::vector<double>(m));

        /*
         * Start by filling everything with either 0 or 1.
         */
        for (int i{}; i < m; ++i)
        {
            for (int j{}; j < m; ++j)
            {
                Hdata[i][j] = (i - j + 1 < 0)
                    ? 0
                    : 1;
            }
        }

        /*
         * Setting up power-array to avoid calculating the same value twice: h_powers[0] = h^1 ...
         * h_powers[m-1] = h^m
         */
        auto h_powers = std::vector<double>(m);
        h_powers[0] = h;
        for (int i{ 1 }; i < m; ++i)
        {
            h_powers[i] = h * h_powers[i - 1];
        }

        /*
         * First column and last row has special values (each other reversed).
         */
        for (int i{}; i < m; ++i)
        {
            Hdata[i][0] = Hdata[i][0] - h_powers[i];
            Hdata[m - 1][i] -= h_powers[m - i - 1];
        }

        /*
         * [1] states: "For 1/2 < h < 1 the bottom left element of the matrix should be (1 - 2*h^m +
         * (2h - 1)^m )/m!" sin_ce 0 <= h < 1, then if h > 1/2 is sufficient to check:
         */
        if (Double.compare(h, 0.5) > 0)
        {
            Hdata[m - 1][0] += std::pow(2 * h - 1, m);
        }

        /*
         * Aside from the first column and last row, the (i, j)-th element is 1/(i - j + 1)! if i -
         * j + 1 >= 0, else 0. 1's and 0's are already put, so only division with (i - j + 1)! is
         * needed in the elements that have 1's. There is no need to calculate (i - j + 1)! and then
         * divide - small steps avoid overflows. Note that i - j + 1 > 0 <=> i + 1 > j instead of
         * j'ing all the way to m. Also note that it is started at g = 2 because dividing by 1 isn't
         * really necessary.
         */
        for (int i{}; i < m; ++i)
        {
            for (int j{}; j < i + 1; ++j)
            {
                if (i - j + 1 > 0)
                {
                    for (const int& g = 2; g <= i - j + 1; ++g)
                    {
                        Hdata[i][j] /= g;
                    }
                }
            }
        }
        return Matrix_Utils::create_real_matrix(Hdata);
    }

    /**
     * Verifies that {@code array} has length at least 2.
     *
     * @param array array to test
     * @org.hipparchus.exception.Null_Argument_Exception if array is NULL
     * @ if array is too short
     */
    void check_array(std::vector<double>& arr)
    {
        //Math_Utils::check_not_null(array);
        if (arr.size() < 2)
        {
            throw std::exception("not implemented");
            // throw (hipparchus::exception::Localized_Core_Formats_Type::INSUFFICIENT_OBSERVED_POINTS_IN_SAMPLE, arr.size(), 2);
        }
    }

    /**
     * Normalizes a value to an integral multiple of 1/mn between 0 and 1.
     * If d < 1/mn, 0 is returned; if d > 1, 1 is returned; if d is very close
     * to an integral multiple of 1/mn, that value is returned; otherwise the
     * returned value is the smallest multiple of 1/mn less than or equal to d.
     *
     * @param d d value
     * @param n first sample size
     * @param m second sample size
     * @return d value suitable for input to exact_p_at_meshpoint(d, m, n)
     */
    double normalize_d(double d, int n, int m)
    {
        const double resolution = 1 / (static_cast<double>(n) * m);
        const double tol = 1e-12;

        // If d is smaller that the first mesh point, return 0
        // If greater than 1, return 1
        if (d < resolution)
        {
            return 0;
        }
        if (d > 1)
        {
            return 1;
        }

        // Normalize d to the smallest mesh point less than or equal to d;
        // except if d is less than tol less than the next mesh point, bump it up
        const double resolutions = d / resolution;
        const double ceil = std::ceil(resolutions);
        return (ceil - resolutions < tol)
            ? ceil * resolution
            : std::floor(resolutions) * resolution;
    }

    /**
     * Computes \(P(D_{n,m} &gt; d)\) where \(D_{n,m}\) is the 2-sample Kolmogorov-Smirnov statistic. See
     * {@link #kolmogorov_smirnov_statistic(std::vector<double>, std::vector<double>)} for the definition of \(D_{n,m}\).
     * <p>
     * The returned probability is exact, implemented by unwinding the recursive function
     * definitions presented in [4].
     *
     * @param d D-statistic value (must be a "meshpoint" - i.e., a possible actual value of D(m,n)).
     * @param n first sample size
     * @param m second sample size
     * @return probability that a randomly selected m-n partition of m + n generates \(D_{n,m}\)
     *         greater than (resp. greater than or equal to) {@code d}
     */
    double exact_p_at_meshpoint(const double& d, const int& n, const int& m)
    {
        const int& nn = std::max(n, m);
        const int mm = std::min(n, m);
        auto u = std::vector<double>(nn + 2);
        const double k = mm * nn * d + 0.5;
        u[1] = 1;
        for (int j{ 1 }; j < nn + 1; j++)
        {
            u[j + 1] = 1;
            if (mm * j > k)
            {
                u[j + 1] = 0;
            }
        }
        for (int i{ 1 }; i < mm + 1; i++)
        {
            const auto w = static_cast<double>(i / static_cast<double>((i + nn)));
            u[1] = w * u[1];
            if (nn * i > k)
            {
                u[1] = 0;
            }
            for (int j{ 1 }; j < nn + 1; j++)
            {
                u[j + 1] = u[j] + u[j + 1] * w;
                if (std::abs(nn * i - mm * j) > k)
                {
                    u[j + 1] = 0;
                }
            }
        }
        return 1 - u[nn + 1];
    }

public:
    /**
     * Construct a Kolmogorov_Smirnov_Test instance.
     */
    Kolmogorov_Smirnov_Test() 
    {
        super();
    }

    /**
     * Construct a Kolmogorov_Smirnov_Test instance providing a seed for the PRNG
     * used by the {@link #bootstrap(std::vector<double>, std::vector<double>, int)} method.
     *
     * @param seed the seed for the PRNG
     */
    Kolmogorov_Smirnov_Test(const long& seed) 
    {
        super();
        gen.set_seed(seed);
    }

    /**
     * Computes the <i>p-value</i>, or <i>observed significance level</i>, of a one-sample <a
     * href="http://en.wikipedia.org/wiki/Kolmogorov-Smirnov_test"> Kolmogorov-Smirnov test</a>
     * evaluating the NULL hypothesis that {@code data} conforms to {@code distribution}. If
     * {@code exact} is true, the distribution used to compute the p-value is computed using
     * extended precision. See {@link #cdf_exact(double, int)}.
     *
     * @param distribution reference distribution
     * @param data sample being being evaluated
     * @param exact whether or not to force exact computation of the p-value
     * @return the p-value associated with the NULL hypothesis that {@code data} is a sample from
     *         {@code distribution}
     * @ if {@code data} does not have length at least 2
     * @org.hipparchus.exception.Null_Argument_Exception if {@code data} is NULL
     */
    double kolmogorov_smirnov_test(const Real_Distribution& distribution, const std::vector<double>& data, bool exact) 
    {
        return 1.0- cdf(kolmogorov_smirnov_statistic(distribution, data), data.size(), exact);
    }

    /**
     * Computes the one-sample Kolmogorov-Smirnov test statistic, \(D_n=\sup_x |F_n(x)-F(x)|\) where
     * \(F\) is the distribution (cdf) function associated with {@code distribution}, \(n\) is the
     * length of {@code data} and \(F_n\) is the empirical distribution that puts mass \(1/n\) at
     * each of the values in {@code data}.
     *
     * @param distribution reference distribution
     * @param data sample being evaluated
     * @return Kolmogorov-Smirnov statistic \(D_n\)
     * @ if {@code data} does not have length at least 2
     * @org.hipparchus.exception.Null_Argument_Exception if {@code data} is NULL
     */
    double kolmogorov_smirnov_statistic(const Real_Distribution& distribution, const std::vector<double>& data) 
    {
        check_array(data);
        const int n = data.size();
        const double nd = n;
        const std::vector<double> data_copy = std::vector<double>(n];
        System.arraycopy(data, 0, data_copy, 0, n);
        Arrays.sort(data_copy);
        double d{};
        for (int i{ 1 }; i <= n; i++) 
        {
            const double yi = distribution.cumulative_probability(data_copy[i - 1]);
            const double curr_d = std::max(yi - (i - 1) / nd, i / nd - yi);
            if (curr_d > d) 
            {
                d = curr_d;
            }
        }
        return d;
    }

    /**
     * Computes the <i>p-value</i>, or <i>observed significance level</i>, of a two-sample <a
     * href="http://en.wikipedia.org/wiki/Kolmogorov-Smirnov_test"> Kolmogorov-Smirnov test</a>
     * evaluating the NULL hypothesis that {@code x} and {@code y} are samples drawn from the same
     * probability distribution. Specifically, what is returned is an estimate of the probability
     * that the {@link #kolmogorov_smirnov_statistic(std::vector<double>, std::vector<double>)} associated with a randomly
     * selected partition of the combined sample into subsamples of sizes {@code x.size()} and
     * {@code y.size()} will strictly exceed (if {@code strict} is {@code true}) or be at least as
     * large as {@code strict = false}) as {@code kolmogorov_smirnov_statistic(x, y)}.
     * <ul>
     * <li>For small samples (where the product of the sample sizes is less than
     * {@link #LARGE_SAMPLE_PRODUCT}), the exact p-value is computed using the method presented
     * in [4], implemented in {@link #exact_p(double, int, int, bool)}. </li>
     * <li>When the product of the sample sizes exceeds {@link #LARGE_SAMPLE_PRODUCT}, the
     * asymptotic distribution of \(D_{n,m}\) is used. See {@link #approximate_p(double, int, int)}
     * for details on the approximation.</li>
     * </ul><p>
     * If {@code x.size() * y.size()} &lt; {@link #LARGE_SAMPLE_PRODUCT} and the combined set of values in
     * {@code x} and {@code y} contains ties, random jitter is added to {@code x} and {@code y} to
     * break ties before computing \(D_{n,m}\) and the p-value. The jitter is uniformly distributed
     * on (-min_delta / 2, min_delta / 2) where min_delta is the smallest pairwise difference between
     * values in the combined sample.</p>
     * <p>
     * If ties are known to be present in the data, {@link #bootstrap(std::vector<double>, std::vector<double>, int, bool)}
     * may be used as an alternative method for estimating the p-value.</p>
     *
     * @param x first sample dataset
     * @param y second sample dataset
     * @param strict whether or not the probability to compute is expressed as a strict inequality
     *        (ignored for large samples)
     * @return p-value associated with the NULL hypothesis that {@code x} and {@code y} represent
     *         samples from the same distribution
     * @ if either {@code x} or {@code y} does not have length at
     *         least 2
     * @org.hipparchus.exception.Null_Argument_Exception if either {@code x} or {@code y} is NULL
     * @see #bootstrap(std::vector<double>, std::vector<double>, int, bool)
     */
    double kolmogorov_smirnov_test(const std::vector<double>& x, const std::vector<double>& y, bool strict) 
    {
        const long length_product = static_cast<long>( x.size() * y.size());
        std::vector<double> xa;
        std::vector<double> ya;
        if (length_product < LARGE_SAMPLE_PRODUCT && has_ties(x,y)) 
        {
            xa = x.clone();
            ya = y.clone();
            fix_ties(xa, ya);
        }
        else 
        {
            xa = x;
            ya = y;
        }
        if (length_product < LARGE_SAMPLE_PRODUCT) 
        {
            return exact_p(kolmogorov_smirnov_statistic(xa, ya), x.size(), y.size(), strict);
        }
        return approximate_p(kolmogorov_smirnov_statistic(x, y), x.size(), y.size());
    }

    /**
     * Computes the <i>p-value</i>, or <i>observed significance level</i>, of a two-sample <a
     * href="http://en.wikipedia.org/wiki/Kolmogorov-Smirnov_test"> Kolmogorov-Smirnov test</a>
     * evaluating the NULL hypothesis that {@code x} and {@code y} are samples drawn from the same
     * probability distribution. Assumes the strict form of the inequality used to compute the
     * p-value. See {@link #kolmogorov_smirnov_test(Real_Distribution, std::vector<double>, bool)}.
     *
     * @param x first sample dataset
     * @param y second sample dataset
     * @return p-value associated with the NULL hypothesis that {@code x} and {@code y} represent
     *         samples from the same distribution
     * @ if either {@code x} or {@code y} does not have length at
     *         least 2
     * @org.hipparchus.exception.Null_Argument_Exception if either {@code x} or {@code y} is NULL
     */
    double kolmogorov_smirnov_test(const std::vector<double>& x, const std::vector<double>& y) 
    {
        return kolmogorov_smirnov_test(x, y, true);
    }

    /**
     * Computes the two-sample Kolmogorov-Smirnov test statistic, \(D_{n,m}=\sup_x |F_n(x)-F_m(x)|\)
     * where \(n\) is the length of {@code x}, \(m\) is the length of {@code y}, \(F_n\) is the
     * empirical distribution that puts mass \(1/n\) at each of the values in {@code x} and \(F_m\)
     * is the empirical distribution of the {@code y} values.
     *
     * @param x first sample
     * @param y second sample
     * @return test statistic \(D_{n,m}\) used to evaluate the NULL hypothesis that {@code x} and
     *         {@code y} represent samples from the same underlying distribution
     * @ if either {@code x} or {@code y} does not have length at
     *         least 2
     * @org.hipparchus.exception.Null_Argument_Exception if either {@code x} or {@code y} is NULL
     */
    double kolmogorov_smirnov_statistic(const std::vector<double>& x, const std::vector<double>& y) 
    {
        return integral_kolmogorov_smirnov_statistic(x, y)/(static_cast<double>((x.size() * static_cast<long>(y.size()));
    }

    /**
     * Computes the <i>p-value</i>, or <i>observed significance level</i>, of a one-sample <a
     * href="http://en.wikipedia.org/wiki/Kolmogorov-Smirnov_test"> Kolmogorov-Smirnov test</a>
     * evaluating the NULL hypothesis that {@code data} conforms to {@code distribution}.
     *
     * @param distribution reference distribution
     * @param data sample being being evaluated
     * @return the p-value associated with the NULL hypothesis that {@code data} is a sample from
     *         {@code distribution}
     * @ if {@code data} does not have length at least 2
     * @org.hipparchus.exception.Null_Argument_Exception if {@code data} is NULL
     */
    double kolmogorov_smirnov_test(const Real_Distribution& distribution, const std::vector<double>& data) 
    {
        return kolmogorov_smirnov_test(distribution, data, false);
    }

    /**
     * Performs a <a href="http://en.wikipedia.org/wiki/Kolmogorov-Smirnov_test"> Kolmogorov-Smirnov
     * test</a> evaluating the NULL hypothesis that {@code data} conforms to {@code distribution}.
     *
     * @param distribution reference distribution
     * @param data sample being being evaluated
     * @param alpha significance level of the test
     * @return true iff the NULL hypothesis that {@code data} is a sample from {@code distribution}
     *         can be rejected with confidence 1 - {@code alpha}
     * @ if {@code data} does not have length at least 2
     * @org.hipparchus.exception.Null_Argument_Exception if {@code data} is NULL
     */
    bool kolmogorov_smirnov_test(const Real_Distribution& distribution, const std::vector<double>& data, const double& alpha) 
    {
        if ((alpha <= 0) || (alpha > 0.5)) 
        {
            throw std::exception("not implemented");
            //throw (Localized_Stat_Formats.OUT_OF_BOUND_SIGNIFICANCE_LEVEL, alpha, 0, 0.5);
        }
        return kolmogorov_smirnov_test(distribution, data) < alpha;
    }

    /**
     * Estimates the <i>p-value</i> of a two-sample
     * <a href="http://en.wikipedia.org/wiki/Kolmogorov-Smirnov_test"> Kolmogorov-Smirnov test</a>
     * evaluating the NULL hypothesis that {@code x} and {@code y} are samples drawn from the same
     * probability distribution. This method estimates the p-value by repeatedly sampling sets of size
     * {@code x.size()} and {@code y.size()} from the empirical distribution of the combined sample.
     * When {@code strict} is true, this is equivalent to the algorithm implemented in the R function
     * {@code ks.boot}, described in <pre>
     * Jasjeet S. Sekhon. 2011. 'Multivariate and Propensity Score Matching
     * Software with Automated Balance Optimization: The Matching //package for R.'
     * Journal of Statistical Software, 42(7): 1-52.
     * </pre>
     * @param x first sample
     * @param y second sample
     * @param iterations number of bootstrap resampling iterations
     * @param strict whether or not the NULL hypothesis is expressed as a strict inequality
     * @return estimated p-value
     */
    double bootstrap(const std::vector<double>& x, const std::vector<double>& y, const int& iterations, bool strict)
    {
        const int x_length = x.size();
        const int y_length = y.size();
        const auto combined = std::vector<double>(x_length + y_length);
        System.arraycopy(x, 0, combined, 0, x_length);
        System.arraycopy(y, 0, combined, x_length, y_length);
        const long d = integral_kolmogorov_smirnov_statistic(x, y);
        int greater_count = 0;
        int equal_count = 0;
        std::vector<double> cur_x;
        std::vector<double> cur_y;
        long cur_d;
        for (int i{}; i < iterations; i++) 
        {
            cur_x = resample(combined, x_length);
            cur_y = resample(combined, y_length);
            cur_d = integral_kolmogorov_smirnov_statistic(cur_x, cur_y);
            if (cur_d > d) 
            {
                greater_count++;
            }
            else if (cur_d == d) 
            {
                equal_count++;
            }
        }
        return strict ? greater_count / static_cast<double>( iterations :
            (greater_count + equal_count) / static_cast<double>( iterations;
    }

    /**
     * Computes {@code bootstrap(x, y, iterations, true)}.
     * This is equivalent to ks.boot(x,y, nboots=iterations) using the R Matching
     * //package function. See #bootstrap(std::vector<double>, std::vector<double>, int, bool).
     *
     * @param x first sample
     * @param y second sample
     * @param iterations number of bootstrap resampling iterations
     * @return estimated p-value
     */
    double bootstrap(const std::vector<double>& x, const std::vector<double>& y, const int& iterations)
    {
        return bootstrap(x, y, iterations, true);
    }

    /**
     * Calculates {@code P(D_n < d)} using the method described in [1] with quick decisions for extreme
     * values given in [2] (see above). The result is not exact as with
     * {@link #cdf_exact(double, int)} because calculations are based on
     * {@code double} rather than {@link org.hipparchus.fraction.Big_Fraction}.
     *
     * @param d statistic
     * @param n sample size
     * @return \(P(D_n &lt; d)\)
     * @Math_Runtime_Exception if algorithm fails to convert {@code h} to a
     *         {@link org.hipparchus.fraction.Big_Fraction} in expressing {@code d} as \((k
     *         - h) / m\) for integer {@code k, m} and \(0 &lt;= h &lt; 1\)
     */
    double cdf(const double& d, const int& n)
    {
        return cdf(d, n, false);
    }

    /**
     * Calculates {@code P(D_n < d)}. The result is exact in the sense that Big_Fraction/BigReal is
     * used everywhere at the expense of very slow execution time. Almost never choose this in real
     * applications unless you are very sure; this is almost solely for verification purposes.
     * Normally, you would choose {@link #cdf(double, int)}. See the class
     * javadoc for definitions and algorithm description.
     *
     * @param d statistic
     * @param n sample size
     * @return \(P(D_n &lt; d)\)
     * @Math_Runtime_Exception if the algorithm fails to convert {@code h} to a
     *         {@link org.hipparchus.fraction.Big_Fraction} in expressing {@code d} as \((k
     *         - h) / m\) for integer {@code k, m} and \(0 &lt;= h &lt; 1\)
     */
    double cdf_exact(const double& d, const int& n)
    {
        return cdf(d, n, true);
    }

    /**
     * Calculates {@code P(D_n < d)} using method described in [1] with quick decisions for extreme
     * values given in [2] (see above).
     *
     * @param d statistic
     * @param n sample size
     * @param exact whether the probability should be calculated exact using
     *        {@link org.hipparchus.fraction.Big_Fraction} everywhere at the expense of
     *        very slow execution time, or if {@code double} should be used convenient places to
     *        gain speed. Almost never choose {@code true} in real applications unless you are very
     *        sure; {@code true} is almost solely for verification purposes.
     * @return \(P(D_n &lt; d)\)
     * @Math_Runtime_Exception if algorithm fails to convert {@code h} to a
     *         {@link org.hipparchus.fraction.Big_Fraction} in expressing {@code d} as \((k
     *         - h) / m\) for integer {@code k, m} and \(0 \lt;= h &lt; 1\).
     */
    double cdf(const double& d, const int& n, bool exact)
    {
        const double ninv{ 1 / (static_cast<double>(n)) };
        const double ninvhalf{ 0.5 * ninv };

        if (d <= ninvhalf) 
        {
            return 0;
        }
        if (ninvhalf < d && d <= ninv) 
        {
            double res{ 1 };
            const double f{ 2 * d - ninv };
            // n! f^n = n*f * (n-1)*f * ... * 1*x
            for (int i{ 1 }; i <= n; ++i) 
            {
                res *= i * f;
            }
            return res;
        }
        if (1 - ninv <= d && d < 1) 
        {
            return 1 - 2 * std::pow(1 - d, n);
        }
        if (1 <= d) 
        {
            return 1;
        }
        if (exact) 
        {
            return exact_k(d, n);
        }
        if (n <= 140) 
        {
            return rounded_k(d, n);
        }
        return pelz_good(d, n);
    }

    /**
     * Computes the Pelz-Good approximation for \(P(D_n &lt; d)\) as described in [2] in the class javadoc.
     *
     * @param d value of d-statistic (x in [2])
     * @param n sample size
     * @return \(P(D_n &lt; d)\)
     */
    double pelz_good(const double& d, const int& n)
    {
        // Change the variable since approximation is for the distribution evaluated at d / sqrt(n)
        const double sqrt_n = std::sqrt(n);
        const double z = d * sqrt_n;
        const double z2 = d * d * n;
        const double z4 = z2 * z2;
        const double z6 = z4 * z2;
        const double z8 = z4 * z4;

        // Compute K_0(z)
        double sum{};
        double z2_term = Math_Utils::PI_SQUARED / (8 * z2);
        int k = 1;
        for (; k < MAXIMUM_PARTIAL_SUM_COUNT; k++) 
        {
            const double k_term = 2 * k - 1;
            const double increment = std::exp(-z2_term * k_term * k_term);
            sum += increment;
            if (increment <= PG_SUM_RELATIVE_ERROR * sum) 
            {
                break;
            }
        }
        if (k == MAXIMUM_PARTIAL_SUM_COUNT) 
        {
            throw std::exception("not implemented");
            //throw Math_Illegal_State_Exception(hipparchus::exception::Localized_Core_Formats_Type::MAX_COUNT_EXCEEDED, MAXIMUM_PARTIAL_SUM_COUNT);
        }
        double ret = sum * std::sqrt(2 * std::numbers::pi) / z;

        // K_1(z)
        // Sum is -inf to inf, but k term is always (k + 1/2) ^ 2, so really have
        // twice the sum from k = 0 to inf (k = -1 is same as 0, -2 same as 1, ...)
        const double two_z2 = 2 * z2;
        sum = 0;
        for (k = 0; k < MAXIMUM_PARTIAL_SUM_COUNT; k++) 
        {
            const double k_term = k + 0.5;
            const double k_term2 = k_term * k_term;
            const double increment = (Math_Utils::PI_SQUARED * k_term2 - z2) * std::exp(-Math_Utils::PI_SQUARED * k_term2 / two_z2);
            sum += increment;
            if (std::abs(increment) < PG_SUM_RELATIVE_ERROR * std::abs(sum)) 
            {
                break;
            }
        }
        if (k == MAXIMUM_PARTIAL_SUM_COUNT) 
        {
            throw std::exception("not implemented");
            //throw Math_Illegal_State_Exception(hipparchus::exception::Localized_Core_Formats_Type::MAX_COUNT_EXCEEDED, MAXIMUM_PARTIAL_SUM_COUNT);
        }
        const double sqrt_half_pi = std::sqrt(std::numbers::pi / 2);
        // Instead of doubling sum, divide by 3 instead of 6
        ret += sum * sqrt_half_pi / (3 * z4 * sqrt_n);

        // K_2(z)
        // Same drill as K_1, but with two doubly infinite sums, all k terms are even powers.
        const double z4_term{ 2 * z4 };
        const double z6_term{ 6 * z6 };
        z2_term = 5 * z2;
        const double pi4 = Math_Utils::PI_SQUARED * Math_Utils::PI_SQUARED;
        sum = 0;
        for (k = 0; k < MAXIMUM_PARTIAL_SUM_COUNT; k++) 
        {
            const double k_term = k + 0.5;
            const double k_term2 = k_term * k_term;
            const double increment =  (z6_term + z4_term + Math_Utils::PI_SQUARED * (z4_term - z2_term) * k_term2 +
                    pi4 * (1 - two_z2) * k_term2 * k_term2) * std::exp(-Math_Utils::PI_SQUARED * k_term2 / two_z2);
            sum += increment;
            if (std::abs(increment) < PG_SUM_RELATIVE_ERROR * std::abs(sum)) 
            {
                break;
            }
        }
        if (k == MAXIMUM_PARTIAL_SUM_COUNT) 
        {
            throw std::exception("not implemented");
            //throw Math_Illegal_State_Exception(hipparchus::exception::Localized_Core_Formats_Type::MAX_COUNT_EXCEEDED, MAXIMUM_PARTIAL_SUM_COUNT);
        }
        double sum2{};
        for (k = 1; k < MAXIMUM_PARTIAL_SUM_COUNT; k++) 
        {
            const double k_term2 = k * k;
            const double increment = Math_Utils::PI_SQUARED * k_term2 * std::exp(-Math_Utils::PI_SQUARED * k_term2 / two_z2);
            sum2 += increment;
            if (std::abs(increment) < PG_SUM_RELATIVE_ERROR * std::abs(sum2)) 
            {
                break;
            }
        }
        if (k == MAXIMUM_PARTIAL_SUM_COUNT) 
        {
            throw std::exception("not implemented");
            //throw Math_Illegal_State_Exception(hipparchus::exception::Localized_Core_Formats_Type::MAX_COUNT_EXCEEDED, MAXIMUM_PARTIAL_SUM_COUNT);
        }
        // Again, adjust coefficients instead of doubling sum, sum2
        ret += (sqrt_half_pi / n) * (sum / (36 * z2 * z2 * z2 * z) - sum2 / (18 * z2 * z));

        // K_3(z) One more time with feeling - two doubly infinite sums, all k powers even.
        // Multiply coefficient denominators by 2, so omit doubling sums.
        const double pi6{ pi4 * Math_Utils::PI_SQUARED };
        sum = 0;
        for (k = 0; k < MAXIMUM_PARTIAL_SUM_COUNT; k++) 
        {
            const double k_term = k + 0.5;
            const double k_term2 = k_term * k_term;
            const double k_term4 = k_term2 * k_term2;
            const double k_term6 = k_term4 * k_term2;
            const double increment = (pi6 * k_term6 * (5 - 30 * z2) + pi4 * k_term4 * (-60 * z2 + 212 * z4) +
                            Math_Utils::PI_SQUARED * k_term2 * (135 * z4 - 96 * z6) - 30 * z6 - 90 * z8) *
                    std::exp(-Math_Utils::PI_SQUARED * k_term2 / two_z2);
            sum += increment;
            if (std::abs(increment) < PG_SUM_RELATIVE_ERROR * std::abs(sum)) 
            {
                break;
            }
        }
        if (k == MAXIMUM_PARTIAL_SUM_COUNT) 
        {
            throw std::exception("not implemented");
            //throw Math_Illegal_State_Exception(hipparchus::exception::Localized_Core_Formats_Type::MAX_COUNT_EXCEEDED, MAXIMUM_PARTIAL_SUM_COUNT);
        }
        sum2 = 0;
        for (k = 1; k < MAXIMUM_PARTIAL_SUM_COUNT; k++) 
        {
            const double k_term2{ k * k };
            const double k_term4{ k_term2 * k_term2 };
            const double increment = (-pi4 * k_term4 + 3 * Math_Utils::PI_SQUARED * k_term2 * z2) *
                    std::exp(-Math_Utils::PI_SQUARED * k_term2 / two_z2);
            sum2 += increment;
            if (std::abs(increment) < PG_SUM_RELATIVE_ERROR * std::abs(sum2)) 
            {
                break;
            }
        }
        if (k == MAXIMUM_PARTIAL_SUM_COUNT) 
        {
            throw std::exception("not implemented");
            //throw Math_Illegal_State_Exception(hipparchus::exception::Localized_Core_Formats_Type::MAX_COUNT_EXCEEDED, MAXIMUM_PARTIAL_SUM_COUNT);
        }
        return ret + (sqrt_half_pi / (sqrt_n * n)) * (sum / (3240 * z6 * z4) +
                + sum2 / (108 * z6));

    }

    /**
     * Computes \( 1 + 2 \sum_{i=1}^\infty (-1)^i e^{-2 i^2 t^2} \) stopping when successive partial
     * sums are within {@code tolerance} of one another, or when {@code max_iterations} partial sums
     * have been computed. If the sum does not converge before {@code max_iterations} iterations a
     * {@link Math_Illegal_State_Exception} is thrown.
     *
     * @param t argument
     * @param tolerance Cauchy criterion for partial sums
     * @param max_iterations maximum number of partial sums to compute
     * @return Kolmogorov sum evaluated at t
     * @Math_Illegal_State_Exception if the series does not converge
     */
    double ks_sum(const double& t, const double& tolerance, const int& max_iterations)
    {
        if (t == 0.0) 
        {
            return 0;
        }

        // TODO: for small t (say less than 1), the alternative expansion in part 3 of [1]
        // from class javadoc should be used.

        const double x{ -2 * t * t };
        int sign{ -1 };
        long i{ 1 };
        double partial_sum{ 0.5 };
        double delta{ 1 };
        while (delta > tolerance && i < max_iterations) 
        {
            delta = std::exp(x * i * i);
            partial_sum += sign * delta;
            sign *= -1;
            i++;
        }
        if (i == max_iterations) 
        {
            throw std::exception("not implemented");
            //throw Math_Illegal_State_Exception(hipparchus::exception::Localized_Core_Formats_Type::MAX_COUNT_EXCEEDED, max_iterations);
        }
        return partial_sum * 2;
    }

    /**
     * Computes \(P(D_{n,m} &gt; d)\) if {@code strict} is {@code true}; otherwise \(P(D_{n,m} \ge
     * d)\), where \(D_{n,m}\) is the 2-sample Kolmogorov-Smirnov statistic. See
     * {@link #kolmogorov_smirnov_statistic(std::vector<double>, std::vector<double>)} for the definition of \(D_{n,m}\).
     * <p>
     * The returned probability is exact, implemented by unwinding the recursive function
     * definitions presented in [4] from the class javadoc.
     * </p>
     *
     * @param d D-statistic value
     * @param n first sample size
     * @param m second sample size
     * @param strict whether or not the probability to compute is expressed as a strict inequality
     * @return probability that a randomly selected m-n partition of m + n generates \(D_{n,m}\)
     *         greater than (resp. greater than or equal to) {@code d}
     */
    double exact_p(const double& d, const int& n, const int& m, bool strict)
    {
        if (d < 1 / static_cast<double>(( m * n)))
        {
            return 1.0;
        }
        else if (d >= 1) 
        {
            return 0;
        }
        double normalize_d = normalize_d(d, n, m);
        if (!strict) 
        {
            normalize_d -= 1 / (static_cast<double>(n * m);
        }
        return exact_p_at_meshpoint(normalize_d, n, m);
    }

    /**
     * Uses the Kolmogorov-Smirnov distribution to approximate \(P(D_{n,m} &gt; d)\) where \(D_{n,m}\)
     * is the 2-sample Kolmogorov-Smirnov statistic. See
     * {@link #kolmogorov_smirnov_statistic(std::vector<double>, std::vector<double>)} for the definition of \(D_{n,m}\).
     * <p>
     * Specifically, what is returned is \(1 - k(d \sqrt{mn / (m + n)})\) where \(k(t) = 1 + 2
     * \sum_{i=1}^\infty (-1)^i e^{-2 i^2 t^2}\). See {@link #ks_sum(double, double, int)} for
     * details on how convergence of the sum is determined. This implementation passes {@code ks_sum}
     * {@link #KS_SUM_CAUCHY_CRITERION} as {@code tolerance} and
     * {@link #MAXIMUM_PARTIAL_SUM_COUNT} as {@code max_iterations}.
     * </p>
     *
     * @param d D-statistic value
     * @param n first sample size
     * @param m second sample size
     * @return approximate probability that a randomly selected m-n partition of m + n generates
     *         \(D_{n,m}\) greater than {@code d}
     */
    double approximate_p(const double& d, const int& n, const int& m)
    {
        return 1 - ks_sum(d * std::sqrt((m * n) / (m + n)), KS_SUM_CAUCHY_CRITERION, MAXIMUM_PARTIAL_SUM_COUNT);
    }

    /**
     * Fills a bool array randomly with a fixed number of {@code true} values.
     * The method uses a simplified version of the Fisher-Yates shuffle algorithm.
     * By processing first the {@code true} values followed by the remaining {@code false} values
     * less random numbers need to be generated. The method is optimized for the case
     * that the number of {@code true} values is larger than or equal to the number of
     * {@code false} values.
     *
     * @param b bool array
     * @param number_of_true_values number of {@code true} values the bool array should constly have
     * @param rng random data generator
     */
    static void fill_boolean_array_randomly_with_fixed_number_true_values(std::vector<bool> b, const int& number_of_true_values, const Random_Generator& rng) 
    {
        Arrays.fill(b, true);
        for (int k{ number_of_true_values }; k < b.size(); k++)
        {
            const int r = rng.next_int(k + 1);
            b[b[r] ? r : k] = false;
        }
    }
};