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
//package org.hipparchus.analysis.integration.gauss;

//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.analysis.Calculus_Field_Univariate_Function;
//import org.hipparchus.exception.;
//import org.hipparchus.util.Pair;
#include <type_traits>
#include "../../../CalculusFieldElement.hpp"

/**
 * This class's : {@link #integrate(Calculus_Field_Univariate_Function) integrate}
 * method assuming that the integral is symmetric about 0.
 * This allows to reduce numerical errors.
 *
 * @param <T> Type of the field elements.
 * @since 2.0
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
class SymmetricFieldGauss_Integrator extends FieldGauss_Integrator<T> 
{
    /**
     * Creates an integrator from the given {@code points} and {@code weights}.
     * The integration interval is defined by the first and last value of
     * {@code points} which must be sorted in increasing order.
     *
     * @param points Integration points.
     * @param weights Weights of the corresponding integration nodes.
     * @ if the {@code points} are not
     * sorted in increasing order.
     * @ if points and weights don't have the same length
     */
    public SymmetricFieldGauss_Integrator(std::vector<T> points, std::vector<T> weights)
         
        {
        super(points, weights);
    }

    /**
     * Creates an integrator from the given pair of points (first element of
     * the pair) and weights (second element of the pair.
     *
     * @param points_and_weights Integration points and corresponding weights.
     * @ if the {@code points} are not
     * sorted in increasing order.
     *
     * @see #SymmetricFieldGauss_Integrator(Calculus_Field_Element[], Calculus_Field_Element[])
     */
    public SymmetricFieldGauss_Integrator(Pair<std::vector<T>, std::vector<T>> points_and_weights)
         
        {
        this(points_and_weights.get_first(), points_and_weights.get_second());
    }

    /**
     * {@inherit_doc}
     */
    //override
    public T integrate(Calculus_Field_Univariate_Function<T> f) 
    {
        const int rule_length = get_number_of_points();

        const T zero = get_point(0).get_field().get_zero();
        if (rule_length == 1) 
        {
            return get_weight(0).multiply(f.value(zero));
        }

        const int i_max = rule_length / 2;
        T s = zero;
        T c = zero;
        for (int i{}; i < i_max; i++) 
        {
            const T p = get_point(i);
            const T w = get_weight(i);

            const T f1 = f.value(p);
            const T f2 = f.value(p.negate());

            const T y = w.multiply(f1.add(f2)).subtract(c);
            const T t = s.add(y);

            c = t.subtract(s).subtract(y);
            s = t;
        }

        if (rule_length % 2 != 0) 
        {
            const T w = get_weight(i_max);

            const T y = w.multiply(f.value(zero)).subtract(c);
            const T t = s.add(y);

            s = t;
        }

        return s;
    }
}


