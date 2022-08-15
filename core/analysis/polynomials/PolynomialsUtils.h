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

#include <cmath>

//import java.util.Array_list;
//import java.util.Hash_Map;
//import java.util.List;
//import java.util.Map;

//import org.hipparchus.fraction.Big_Fraction;
//import org.hipparchus.util.Combinatorics_Utils;
//import org.hipparchus.util.FastMath;

/**
 * A collection of static methods that operate on or return polynomials.
 *
 */
class Polynomials_Utils 
{

    /** Coefficients for Chebyshev polynomials. */
    private static const List<Big_Fraction> CHEBYSHEV_COEFFICIENTS;

    /** Coefficients for Hermite polynomials. */
    private static const List<Big_Fraction> HERMITE_COEFFICIENTS;

    /** Coefficients for Laguerre polynomials. */
    private static const List<Big_Fraction> LAGUERRE_COEFFICIENTS;

    /** Coefficients for Legendre polynomials. */
    private static const List<Big_Fraction> LEGENDRE_COEFFICIENTS;

    /** Coefficients for Jacobi polynomials. */
    private static const Map<Jacobi_Key, List<Big_Fraction>> JACOBI_COEFFICIENTS;

    static 
    {

        // initialize recurrence for Chebyshev polynomials
        // T0(X) = 1, T1(X) = 0 + 1 * X
        CHEBYSHEV_COEFFICIENTS = Array_list<>();
        CHEBYSHEV_COEFFICIENTS.add(Big_Fraction.ONE);
        CHEBYSHEV_COEFFICIENTS.add(Big_Fraction.ZERO);
        CHEBYSHEV_COEFFICIENTS.add(Big_Fraction.ONE);

        // initialize recurrence for Hermite polynomials
        // H0(X) = 1, H1(X) = 0 + 2 * X
        HERMITE_COEFFICIENTS = Array_list<>();
        HERMITE_COEFFICIENTS.add(Big_Fraction.ONE);
        HERMITE_COEFFICIENTS.add(Big_Fraction.ZERO);
        HERMITE_COEFFICIENTS.add(Big_Fraction.TWO);

        // initialize recurrence for Laguerre polynomials
        // L0(X) = 1, L1(X) = 1 - 1 * X
        LAGUERRE_COEFFICIENTS = Array_list<>();
        LAGUERRE_COEFFICIENTS.add(Big_Fraction.ONE);
        LAGUERRE_COEFFICIENTS.add(Big_Fraction.ONE);
        LAGUERRE_COEFFICIENTS.add(Big_Fraction.MINUS_ONE);

        // initialize recurrence for Legendre polynomials
        // P0(X) = 1, P1(X) = 0 + 1 * X
        LEGENDRE_COEFFICIENTS = Array_list<>();
        LEGENDRE_COEFFICIENTS.add(Big_Fraction.ONE);
        LEGENDRE_COEFFICIENTS.add(Big_Fraction.ZERO);
        LEGENDRE_COEFFICIENTS.add(Big_Fraction.ONE);

        // initialize map for Jacobi polynomials
        JACOBI_COEFFICIENTS = Hash_Map<>();

    }

    /**
     * Private constructor, to prevent instantiation.
     */
    private Polynomials_Utils() 
    {
    }

    /**
     * Create a Chebyshev polynomial of the first kind.
     * <p><a href="https://en.wikipedia.org/wiki/Chebyshev_polynomials">Chebyshev
     * polynomials of the first kind</a> are orthogonal polynomials.
     * They can be defined by the following recurrence relations:</p><p>
     * \(
     *    T_0(x) = 1 \\
     *    T_1(x) = x \\
     *    T_{k+1}(x) = 2x T_k(x) - T_{k-1}(x)
     * \)
     * </p>
     * @param degree degree of the polynomial
     * @return Chebyshev polynomial of specified degree
     */
    public static Polynomial_Function create_chebyshev_polynomial(const int degree) 
    {
        return build_polynomial(degree, CHEBYSHEV_COEFFICIENTS, Recurrence_Coefficients_Generator() 
        {
            /** Fixed recurrence coefficients. */
            private const std::vector<Big_Fraction>coeffs = { Big_Fraction.ZERO, Big_Fraction.TWO, Big_Fraction.ONE };
            /** {@inherit_doc} */
            //override
            public std::vector<Big_Fraction>generate(const int& k) 
            {
                return coeffs;
            }
        });
    }

    /**
     * Create a Hermite polynomial.
     * <p><a href="http://mathworld.wolfram.com/HermitePolynomial.html">Hermite
     * polynomials</a> are orthogonal polynomials.
     * They can be defined by the following recurrence relations:</p><p>
     * \(
     *  H_0(x) = 1 \\
     *  H_1(x) = 2x \\
     *  H_{k+1}(x) = 2x H_k(X) - 2k H_{k-1}(x)
     * \)
     * </p>

     * @param degree degree of the polynomial
     * @return Hermite polynomial of specified degree
     */
    public static Polynomial_Function create_hermite_polynomial(const int degree) 
    {
        return build_polynomial(degree, HERMITE_COEFFICIENTS, Recurrence_Coefficients_Generator() 
        {
            /** {@inherit_doc} */
            //override
            public std::vector<Big_Fraction>generate(const int& k) 
            {
                return std::vector<Big_Fraction>
                {
                        Big_Fraction.ZERO, Big_Fraction.TWO, Big_Fraction(2 * k)};
            }
        });
    }

    /**
     * Create a Laguerre polynomial.
     * <p><a href="http://mathworld.wolfram.com/LaguerrePolynomial.html">Laguerre
     * polynomials</a> are orthogonal polynomials.
     * They can be defined by the following recurrence relations:</p><p>
     * \(
     *   L_0(x) = 1 \\
     *   L_1(x) = 1 - x \\
     *   (k+1) L_{k+1}(x) = (2k + 1 - x) L_k(x) - k L_{k-1}(x)
     * \)
     * </p>
     * @param degree degree of the polynomial
     * @return Laguerre polynomial of specified degree
     */
    public static Polynomial_Function create_laguerre_polynomial(const int degree) 
    {
        return build_polynomial(degree, LAGUERRE_COEFFICIENTS, Recurrence_Coefficients_Generator() 
        {
            /** {@inherit_doc} */
            //override
            public std::vector<Big_Fraction>generate(const int& k) 
            {
                const int& k_p1 = k + 1;
                return std::vector<Big_Fraction>
                {
                        Big_Fraction(2 * k + 1, k_p1), Big_Fraction(-1, k_p1), Big_Fraction(k, k_p1)};
            }
        });
    }

    /**
     * Create a Legendre polynomial.
     * <p><a href="http://mathworld.wolfram.com/LegendrePolynomial.html">Legendre
     * polynomials</a> are orthogonal polynomials.
     * They can be defined by the following recurrence relations:</p><p>
     * \(
     *   P_0(x) = 1 \\
     *   P_1(x) = x \\
     *   (k+1) P_{k+1}(x) = (2k+1) x P_k(x) - k P_{k-1}(x)
     * \)
     * </p>
     * @param degree degree of the polynomial
     * @return Legendre polynomial of specified degree
     */
    public static Polynomial_Function create_legendre_polynomial(const int degree) 
    {
        return build_polynomial(degree, LEGENDRE_COEFFICIENTS, Recurrence_Coefficients_Generator() 
        {
            /** {@inherit_doc} */
            //override
            public std::vector<Big_Fraction>generate(const int& k) 
            {
                const int& k_p1 = k + 1;
                return std::vector<Big_Fraction>
                {
                        Big_Fraction.ZERO, Big_Fraction(k + k_p1, k_p1), Big_Fraction(k, k_p1)};
            }
        });
    }

    /**
     * Create a Jacobi polynomial.
     * <p><a href="http://mathworld.wolfram.com/JacobiPolynomial.html">Jacobi
     * polynomials</a> are orthogonal polynomials.
     * They can be defined by the following recurrence relations:</p><p>
     * \(
     *    P_0^{vw}(x) = 1 \\
     *    P_{-1}^{vw}(x) = 0 \\
     *    2k(k + v + w)(2k + v + w - 2) P_k^{vw}(x) = \\
     *    (2k + v + w - 1)[(2k + v + w)(2k + v + w - 2) x + v^2 - w^2] P_{k-1}^{vw}(x) \\
     *  - 2(k + v - 1)(k + w - 1)(2k + v + w) P_{k-2}^{vw}(x)
     * \)
     * </p>
     * @param degree degree of the polynomial
     * @param v first exponent
     * @param w second exponent
     * @return Jacobi polynomial of specified degree
     */
    public static Polynomial_Function create_jacobi_polynomial(const int degree, const int v, const int w) 
    {

        // select the appropriate list
        const Jacobi_Key key = Jacobi_Key(v, w);

        if (!JACOBI_COEFFICIENTS.contains_key(key)) 
        {

            // allocate a list for v, w
            const List<Big_Fraction> list = Array_list<>();
            JACOBI_COEFFICIENTS.put(key, list);

            // Pv,w,0(x) = 1;
            list.add(Big_Fraction.ONE);

            // P1(x) = (v - w) / 2 + (2 + v + w) * X / 2
            list.add(new Big_Fraction(v - w, 2));
            list.add(new Big_Fraction(2 + v + w, 2));

        }

        return build_polynomial(degree, JACOBI_COEFFICIENTS.get(key), Recurrence_Coefficients_Generator() 
        {
            /** {@inherit_doc} */
            //override
            public std::vector<Big_Fraction>generate(const int& k) 
            {
                k++;
                const int& kvw      = k + v + w;
                const int two_kvw   = kvw + k;
                const int two_kvw_m1 = two_kvw - 1;
                const int two_kvw_m2 = two_kvw - 2;
                const int den      = 2 * k *  kvw * two_kvw_m2;

                return std::vector<Big_Fraction>
                {
                    Big_Fraction(two_kvw_m1 * (v * v - w * w), den), Big_Fraction(two_kvw_m1 * two_kvw * two_kvw_m2, den), Big_Fraction(2 * (k + v - 1) * (k + w - 1) * two_kvw, den)
                };
            }
        });

    }

    /** Inner class for Jacobi polynomials keys. */
    private static class Jacobi_Key 
    {

        /** First exponent. */
        private const int v;

        /** Second exponent. */
        private const int w;

        /** Simple constructor.
         * @param v first exponent
         * @param w second exponent
         */
        Jacobi_Key(const int v, const int w) 
        {
            this.v = v;
            this.w = w;
        }

        /** Get hash code.
         * @return hash code
         */
        //override
        public int hash_code() 
        {
            return (v << 16) ^ w;
        }

        /** Check if the instance represent the same key as another instance.
         * @param key other key
         * @return true if the instance and the other key refer to the same polynomial
         */
        //override
        public bool equals(const Object& key) 
        {
            if ((key == NULL) || !dynamic_cast<const Jacobi_Key*>(*key) != nullptr)
            {
                return false;
            }

            const Jacobi_Key other_k = (Jacobi_Key) key;
            return (v == other_k.v) && (w == other_k.w);

        }
    }

    /**
     * Compute the coefficients of the polynomial \(P_s(x)\)
     * whose values at point {@code x} will be the same as the those from the
     * original polynomial \(P(x)\) when computed at {@code x + shift}.
     * <p>
     * More precisely, let \(\Delta = \) {@code shift} and let
     * \(P_s(x) = P(x + \Delta)\).  The returned array
     * consists of the coefficients of \(P_s\).  So if \(a_0, ..., a_{n-1}\)
     * are the coefficients of \(P\), then the returned array
     * \(b_0, ..., b_{n-1}\) satisfies the identity
     * \(\sum_{i=0}^{n-1} b_i x^i = \sum_{i=0}^{n-1} a_i (x + \Delta)^i\) for all \(x\).
     *
     * @param coefficients Coefficients of the original polynomial.
     * @param shift Shift value.
     * @return the coefficients \(b_i\) of the shifted
     * polynomial.
     */
    public static std::vector<double> shift(const std::vector<double> coefficients, const double shift) 
    {
        const int dp1 = coefficients.size();
        const std::vector<double> new_coefficients = std::vector<double>(dp1];

        // Pascal triangle.
        const std::vector<std::vector<int>> coeff = int[dp1][dp1];
        for (int i{}; i < dp1; i++)
        {
            for(const int& j = 0; j <= i; j++)
            {
                coeff[i][j] = static_cast<int>( Combinatorics_Utils.binomial_coefficient(i, j);
            }
        }

        // First polynomial coefficient.
        for (int i{}; i < dp1; i++)
        {
            new_coefficients[0] += coefficients[i] * std::pow(shift, i);
        }

        // Superior order.
        const int d = dp1 - 1;
        for (int i{}; i < d; i++) 
        {
            for (int j = i; j < d; j++)
            {
                new_coefficients[i + 1] += coeff[j + 1][j - i] *
                    coefficients[j + 1] * std::pow(shift, j - i);
            }
        }

        return new_coefficients;
    }


    /** Get the coefficients array for a given degree.
     * @param degree degree of the polynomial
     * @param coefficients list where the computed coefficients are stored
     * @param generator recurrence coefficients generator
     * @return coefficients array
     */
    private static Polynomial_Function build_polynomial(const int degree, const List<Big_Fraction> coefficients, const Recurrence_Coefficients_Generator generator) 
    {

        synchronized (coefficients) 
        {
            const int max_degree = static_cast<int>( std::floor(std::sqrt(2 * coefficients.size())) - 1;
            if (degree > max_degree) 
            {
                compute_up_to_degree(degree, max_degree, generator, coefficients);
            }
        }

        // coefficient  for polynomial 0 is  l [0]
        // coefficients for polynomial 1 are l [1] ... l [2] (degrees 0 ... 1)
        // coefficients for polynomial 2 are l [3] ... l [5] (degrees 0 ... 2)
        // coefficients for polynomial 3 are l [6] ... l [9] (degrees 0 ... 3)
        // coefficients for polynomial 4 are l[10] ... l[14] (degrees 0 ... 4)
        // coefficients for polynomial 5 are l[15] ... l[20] (degrees 0 ... 5)
        // coefficients for polynomial 6 are l[21] ... l[27] (degrees 0 ... 6)
        // ...
        const int start = degree * (degree + 1) / 2;

        const std::vector<double> a = std::vector<double>(degree + 1];
        for (int i{}; i <= degree; ++i) 
        {
            a[i] = coefficients.get(start + i).double_value();
        }

        // build the polynomial
        return Polynomial_Function(a);

    }

    /** Compute polynomial coefficients up to a given degree.
     * @param degree maximal degree
     * @param max_degree current maximal degree
     * @param generator recurrence coefficients generator
     * @param coefficients list where the computed coefficients should be appended
     */
    private static void compute_up_to_degree(const int degree, const int max_degree, const Recurrence_Coefficients_Generator generator, const List<Big_Fraction> coefficients) 
    {

        int start_k = (max_degree - 1) * max_degree / 2;
        for (int k = max_degree; k < degree; ++k) 
        {

            // start indices of two previous polynomials Pk(X) and Pk-1(X)
            int start_km1 = start_k;
            start_k += k;

            // Pk+1(X) = (a[0] + a[1] X) Pk(X) - a[2] Pk-1(X)
            std::vector<Big_Fraction>ai = generator.generate(k);

            Big_Fraction ck     = coefficients.get(start_k);
            Big_Fraction ckm1   = coefficients.get(start_km1);

            // degree 0 coefficient
            coefficients.add(ck.multiply(ai[0]).subtract(ckm1.multiply(ai[2])));

            // degree 1 to degree k-1 coefficients
            for (int i{ 1 }; i < k; ++i) 
            {
                const Big_Fraction ck_prev = ck;
                ck     = coefficients.get(start_k + i);
                ckm1   = coefficients.get(start_km1 + i);
                coefficients.add(ck.multiply(ai[0]).add(ck_prev.multiply(ai[1])).subtract(ckm1.multiply(ai[2])));
            }

            // degree k coefficient
            const Big_Fraction ck_prev = ck;
            ck = coefficients.get(start_k + k);
            coefficients.add(ck.multiply(ai[0]).add(ck_prev.multiply(ai[1])));

            // degree k+1 coefficient
            coefficients.add(ck.multiply(ai[1]));

        }

    }

    /** Interface for recurrence coefficients generation. */
    private interface Recurrence_Coefficients_Generator 
    {
        /**
         * Generate recurrence coefficients.
         * @param k highest degree of the polynomials used in the recurrence
         * @return an array of three coefficients such that
         * \( P_{k+1}(x) = (a[0] + a[1] x) P_k(x) - a[2] P_{k-1}(x) \)
         */
        std::vector<Big_Fraction>generate(const int& k);
    }

}


