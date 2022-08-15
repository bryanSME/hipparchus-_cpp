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

//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Field_Sin_Cos;
#include <type_traits>
#include "../../../CalculusFieldElement.hpp"

/** Algorithm computing Jacobi theta functions.
 * @param <T> the type of the field elements
 * @since 2.0
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
class Field_Jacobi_Theta
{
private:
    /** Maximum number of terms in the Fourier series. */
    static constexpr int N_MAX{ 100 };

    /** Nome. */
    const T my_q;

    /** q\xc2\xb2. */
    const T my_q_Square;

    /** \xe2\x88\x9cq. */
    const T my_q_fourth;

public:
    /** Simple constructor.
     * <p>
     * The nome {@code q} can be computed using ratios of complete elliptic integrals
     * ({@link org.hipparchus.special.elliptic.legendre.Legendre_Elliptic_Integral#nome(Calculus_Field_Element)
     * Legendre_Elliptic_Integral.nome(m)} which are themselves defined in term of parameter m, * where m=k\xc2\xb2 and k is the elliptic modulus.
     * </p>
     * @param q nome
     */
    Field_Jacobi_Theta(const T q)
    {
        this.q = q;
        this.my_q_Square = q.multiply(q);
        this.my_q_fourth = std::sqrt(std::sqrt(q));
    }

    /** Get the nome.
     * @return nome
     */
    T get_q()
    {
        return q;
    }

    /** Evaluate the Jacobi theta functions.
     * @param z argument of the functions
     * @return container for the four Jacobi theta functions \xce\xb8\xe2\x82\x81(z|\xcf\x84), \xce\xb8\xe2\x82\x82(z|\xcf\x84), \xce\xb8\xe2\x82\x83(z|\xcf\x84), and \xce\xb8\xe2\x82\x84(z|\xcf\x84)
     */
    Field_Theta<T> values(const T z)
    {

        // the computation is based on Fourier series, // see Digital Library of Mathematical Functions section 20.2
        // https://dlmf.nist.gov/20.2
        const T zero = q.get_field().get_zero();
        const T one = q.get_field().get_one();

        // base angle for Fourier Series
        const Field_Sin_Cos<T> sc1 = Sin_Cos(z);

        // recursion rules initialization
        double         sgn = 1.0;
        T              qNN = one;
        T              qTwoN = one;
        T              qNNp1 = one;
        Field_Sin_Cos<T> sc2n1 = sc1;
        const double   eps = FastMath.ulp(one).get_real();

        // Fourier series
        T sum1 = sc1.sin();
        T sum2 = sc1.cos();
        T sum3 = zero;
        T sum4 = zero;
        for (const int n{ 1 }; n < N_MAX; ++n)
        {

            sgn = -sgn;                            // (-1)\xe2\x81\xbf\xe2\x81\xbb\xc2\xb9     \xe2\x86\x90 (-1)\xe2\x81\xbf
            qNN = qNN.multiply(qTwoN).multiply(q); // q\xe2\x81\xbd\xe2\x81\xbf\xe2\x81\xbb\xc2\xb9\xe2\x81\xbe\xe2\x81\xbd\xe2\x81\xbf\xe2\x81\xbb\xc2\xb9\xe2\x81\xbe \xe2\x86\x90 q\xe2\x81\xbf\xe2\x81\xbf
            qTwoN = qTwoN.multiply(my_q_Square);         // q\xc2\xb2\xe2\x81\xbd\xe2\x81\xbf\xe2\x81\xbb\xc2\xb9\xe2\x81\xbe     \xe2\x86\x90 q\xc2\xb2\xe2\x81\xbf
            qNNp1 = qNNp1.multiply(qTwoN);           // q\xe2\x81\xbd\xe2\x81\xbf\xe2\x81\xbb\xc2\xb9\xe2\x81\xbe\xe2\x81\xbf     \xe2\x86\x90 q\xe2\x81\xbf\xe2\x81\xbd\xe2\x81\xbf\xe2\x81\xba\xc2\xb9\xe2\x81\xbe

            sc2n1 = Field_Sin_Cos.sum(sc2n1, sc1); // {sin|cos}([2n-1] z) \xe2\x86\x90 {sin|cos}(2n z)
            sum3 = sum3.add(sc2n1.cos().multiply(qNN));
            sum4 = sum4.add(sc2n1.cos().multiply(qNN.multiply(sgn)));

            sc2n1 = Field_Sin_Cos.sum(sc2n1, sc1); // {sin|cos}(2n z) \xe2\x86\x90 {sin|cos}([2n+1] z)
            sum1 = sum1.add(sc2n1.sin().multiply(qNNp1.multiply(sgn)));
            sum2 = sum2.add(sc2n1.cos().multiply(qNNp1));

            if (qNNp1.norm() <= eps)
            {
                // we have reach convergence
                break;
            }

        }

        return Field_Theta<>(sum1.multiply(my_q_fourth.multiply(2)), sum2.multiply(my_q_fourth.multiply(2)), sum3.multiply(2).add(1), sum4.multiply(2).add(1));

    }

};