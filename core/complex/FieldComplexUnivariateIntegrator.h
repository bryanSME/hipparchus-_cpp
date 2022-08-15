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
//package org.hipparchus.complex;

//import java.util.function.Function;

//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.analysis.Calculus_Field_Univariate_Function;
//import org.hipparchus.analysis.integration.Field_Univariate_Integrator;
#include <type_traits>
#include "../CalculusFieldElement.hpp"

/**
 * Wrapper to perform univariate complex integration using an underlying real integration algorithms.
 * @param <T> the type of the field elements
 * @since 2.0
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
class Field_Complex_Univariate_Integrator  
{

    /** Underlying real integrator. */
    private Field_Univariate_Integrator<T> integrator;

    /** Crate a complex integrator from a real integrator.
     * @param integrator underlying real integrator to use
     */
    public Field_Complex<double>_Univariate_Integrator(const Field_Univariate_Integrator<T> integrator) 
    {
        this.integrator = integrator;
    }

    /**
     * Integrate a function along a straight path between points.
     *
     * @param max_eval maximum number of evaluations (real and imaginary
     * parts are evaluated separately, so up to twice this number may be used)
     * @param f the integrand function
     * @param start start point of the integration path
     * @param end end point of the integration path
     * @return the value of integral along the straight path
     */
    public Field_Complex<T> integrate(const int max_eval, const Calculus_Field_Univariate_Function<Field_Complex<T>> f, const Field_Complex<T> start, const Field_Complex<T> end) 
    {

        // linear mapping from real interval [0; 1] to function value along complex straight path from start to end
        const Field_Complex<T>              rate   = end.subtract(start);
        const Function<T, Field_Complex<T>> mapped = t -> f.value(start.add(rate.multiply(t)));

        const T zero = start.get_real_part().get_field().get_zero();
        const T one  = start.get_real_part().get_field().get_one();

        // integrate real and imaginary parts separately
        const T real      = integrator.integrate(max_eval, t -> mapped.apply(t).get_real_part(),      zero, one);
        const T imaginary = integrator.integrate(max_eval, t -> mapped.apply(t).get_imaginary_part(), zero, one);

        // combine integrals
        return Field_Complex<double><>(real, imaginary).multiply(rate);

    }

    /**
     * Integrate a function along a polyline path between any number of points.
     *
     * @param max_eval maximum number of evaluations (real and imaginary
     * parts are evaluated separately and each path segments are also evaluated
     * separately, so up to 2n times this number may be used for n segments)
     * @param f the integrand function
     * @param start start point of the integration path
     * @param path successive points defining the path vertices
     * @return the value of integral along the polyline path
     */
    public Field_Complex<T> integrate(const int max_eval, const Calculus_Field_Univariate_Function<Field_Complex<T>> f, const Field_Complex<T> start, ////@Suppress_Warnings("unchecked") const Field_Complex<T>...path) 
    {
        Field_Complex<T> sum      = start.new_instance(0);
        Field_Complex<T> previous = start;
        for (const Field_Complex<T> current : path) 
        {
            sum = sum.add(integrate(max_eval, f, previous, current));
            previous = current;
        }
        return sum;
    }

}


