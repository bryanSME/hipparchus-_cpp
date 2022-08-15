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
//package org.hipparchus.analysis.integration.gauss;

//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.Field;
//import org.hipparchus.exception.;
//import org.hipparchus.util.Math_Arrays;
//import org.hipparchus.util.Pair;
#include <type_traits>
#include "../../../CalculusFieldElement.hpp"

/**
 * Factory that creates Gauss-type quadrature rule using Legendre polynomials.
 * In this implementation, the lower and upper bounds of the natural interval
 * of integration are -1 and 1, respectively.
 * The Legendre polynomials are evaluated using the recurrence relation
 * presented in <a href="http://en.wikipedia.org/wiki/Abramowitz_and_Stegun">
 * Abramowitz and Stegun, 1964</a>.
 *
 * @param <T> Type of the number used to represent the points and weights of
 * the quadrature rules.
 * @since 2.0
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
class FieldLegendreRule_Factory extends FieldAbstractRule_Factory<T> 
{

    /** Simple constructor
     * @param field field to which rule coefficients belong
     */
    public FieldLegendreRule_Factory(const Field<T> field) 
    {
        super(field);
    }

    /** {@inherit_doc} */
    //override
    public Pair<std::vector<T>, std::vector<T>> compute_rule(const int& number_of_points)
         
        {

        const Field<T> field = get_field();

        if (number_of_points == 1) 
        {
            // Break recursion.
            const std::vector<T> points  = Math_Arrays::build_array(field, number_of_points);
            const std::vector<T> weights = Math_Arrays::build_array(field, number_of_points);
            points[0]  = field.get_zero();
            weights[0] = field.get_zero().new_instance(2);
            return Pair<>(points, weights);
        }

        // find nodes as roots of Legendre polynomial
        const Legendre<T> p      =  Legendre<>(number_of_points);
        const std::vector<T>         points = find_roots(number_of_points, p::ratio);
        enforce_symmetry(points);

        // compute weights
        const std::vector<T> weights = Math_Arrays::build_array(field, number_of_points);
        for (int i{}; i <= number_of_points / 2; i++) 
        {
            const T c = points[i];
            const std::vector<T> pKpKm1 = p.pNpNm1(c);
            const T d = pKpKm1[1].subtract(c.multiply(pKpKm1[0])).multiply(number_of_points);
            weights[i] = c.multiply(c).subtract(1).multiply(-2).divide(d.multiply(d));

            // symmetrical point
            const int idx = number_of_points - i - 1;
            weights[idx]  = weights[i];

        }

        return Pair<>(points, weights);

    }

    /** Legendre polynomial.
     * @param <T> Type of the field elements.
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    private static class Legendre 
    {

        /** Degree. */
        private int degree;

        /** Simple constructor.
         * @param degree polynomial degree
         */
        Legendre(const int& degree) 
        {
            this.degree = degree;
        }

        /** Compute ratio P(x)/P'(x).
         * @param x point at which ratio must be computed
         * @return ratio P(x)/P'(x)
         */
        public T ratio(T x) 
        {
            T pm = x.get_field().get_one();
            T p  = x;
            T d  = x.get_field().get_one();
            for (const int n{ 1 }; n < degree; n++) 
            {
                // apply recurrence relations (n+1) P_n+1(x)  = (2n+1) x P_n(x) - n P_n-1(x)
                // and                              P'_n+1(x) = (n+1) P_n(x) + x P'_n(x)
                const T pp = p.multiply(x.multiply(2 * n + 1)).subtract(pm.multiply(n)).divide(n + 1);
                d  = p.multiply(n + 1).add(d.multiply(x));
                pm = p;
                p  = pp;
            }
            return p.divide(d);
        }

        /** Compute P_n(x) and P_n-1(x).
         * @param x point at which polynomials are evaluated
         * @return array containing P_n(x) at index 0 and P_n-1(x) at index 1
         */
        private std::vector<T> pNpNm1(const T& x) 
        {
            std::vector<T> p = Math_Arrays::build_array(x.get_field(), 2);
            p[0] = x;
            p[1] = x.get_field().get_one();
            for (const int n{ 1 }; n < degree; n++) 
            {
                // apply recurrence relation (n+1) P_n+1(x) = (2n+1) x P_n(x) - n P_n-1(x)
                const T pp = p[0].multiply(x.multiply(2 * n + 1)).subtract(p[1].multiply(n)).divide(n + 1);
                p[1] = p[0];
                p[0] = pp;
            }
            return p;
        }

    }

}


