#pragma once
/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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
//package org.hipparchus.special.elliptic.jacobi;

//import org.hipparchus.complex.std::complex<double>;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Field_Sin_Cos;
//import org.hipparchus.util.Precision;
#include <cmath>

/** Algorithm computing Jacobi theta functions.
 * @since 2.0
 */
class Jacobi_Theta
{
private:
    /** Maximum number of terms in the Fourier series. */
    static constexpr int N_MAX{ 100 };

    /** Nome. */
    const double my_q;

    /** q². */
    const double my_q_square;

    /** ∜q. */
    const double my_q_fourth;

public:
    /** Simple constructor.
     * <p>
     * The nome {@code q} can be computed using ratios of complete elliptic integrals
     * ({@link org.hipparchus.special.elliptic.legendre.Legendre_Elliptic_Integral#nomestatic_cast<double>(
     * Legendre_Elliptic_Integral.nome(m)} which are themselves defined in term of parameter m, * where m=k² and k is the elliptic modulus.
     * </p>
     * @param q nome
     */
    Jacobi_Theta(const double q) : my_q{ q }, my_q_square{ q * q }, my_q_fourth{ std::sqrt(std::sqrt(q)) } {};

    /** Get the nome.
     * @return nome
     */
    double get_q() const
    {
        return my_q;
    }

    /** Evaluate the Jacobi theta functions.
     * @param z argument of the functions
     * @return container for the four Jacobi theta functions θ₁(z|τ), θ₂(z|τ), θ₃(z|τ), and θ₄(z|τ)
     */
    public Theta values(const std::complex<double> z)
    {

        // the computation is based on Fourier series, // see Digital Library of Mathematical Functions section 20.2
        // https://dlmf.nist.gov/20.2

        // base angle for Fourier Series
        const Field_Sin_Cos<std::complex<double>> sc1 = Sin_Cos(z);

        // recursion rules initialization
        double               sgn = 1.0;
        double               qNN = 1.0;
        double               qTwoN = 1.0;
        double               qNNp1 = 1.0;
        Field_Sin_Cos<std::complex<double>> sc2n1 = sc1;

        // Fourier series
        std::complex<double> sum1 = sc1.sin();
        std::complex<double> sum2 = sc1.cos();
        std::complex<double> sum3 = std::complex<double>.ZERO;
        std::complex<double> sum4 = std::complex<double>.ZERO;
        for (const int n{ 1 }; n < N_MAX; ++n)
        {

            sgn = -sgn;              // (-1)ⁿ⁻¹     ← (-1)ⁿ
            qNN = qNN * qTwoN * q; // q⁽ⁿ⁻¹⁾⁽ⁿ⁻¹⁾ ← qⁿⁿ
            qTwoN = qTwoN * my_q_Square;   // q²⁽ⁿ⁻¹⁾     ← q²ⁿ
            qNNp1 = qNNp1 * qTwoN;     // q⁽ⁿ⁻¹⁾ⁿ     ← qⁿ⁽ⁿ⁺¹⁾

            sc2n1 = Field_Sin_Cos.sum(sc2n1, sc1); // {sin|cos}([2n-1] z) ← {sin|cos}(2n z)
            sum3 = sum3.add(sc2n1.cos().multiply(qNN));
            sum4 = sum4.add(sc2n1.cos().multiply(sgn * qNN));

            sc2n1 = Field_Sin_Cos.sum(sc2n1, sc1); // {sin|cos}(2n z) ← {sin|cos}([2n+1] z)
            sum1 = sum1.add(sc2n1.sin().multiply(sgn * qNNp1));
            sum2 = sum2.add(sc2n1.cos().multiply(qNNp1));

            if (std::abs(qNNp1) <= Precision.EPSILON)
            {
                // we have reach convergence
                break;
            }

        }

        return Theta(
            sum1.multiply(2 * my_q_fourth),
            sum2.multiply(2 * my_q_fourth),
            sum3.multiply(2).add(1),
            sum4.multiply(2).add(1)
        );
    }
};