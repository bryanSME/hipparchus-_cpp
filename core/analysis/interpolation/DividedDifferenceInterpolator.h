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
//package org.hipparchus.analysis.interpolation;

//import java.io.Serializable;

//import org.hipparchus.analysis.polynomials.Polynomial_Function_Lagrange_Form;
//import org.hipparchus.analysis.polynomials.Polynomial_FunctionNewtonForm;
//import org.hipparchus.exception.;

/**
 * Implements the <a href=
 * "http://mathworld.wolfram.com/NewtonsDividedDifferenceInterpolationFormula.html">
 * Divided Difference Algorithm</a> for interpolation of real univariate
 * functions. For reference, see <b>Introduction to Numerical Analysis</b>, * ISBN 038795452X, chapter 2.
 * <p>
 * The actual code of Neville's evaluation is in Polynomial_Function_Lagrange_Form, * this class provides an easy-to-use interface to it.</p>
 *
 */
class Divided_Difference_Interpolator
    : Univariate_Interpolator
    {
    /** serializable version identifier */
    107049519551235069L;

    /**
     * Compute an interpolating function for the dataset.
     *
     * @param x Interpolating points array.
     * @param y Interpolating values array.
     * @return a function which interpolates the dataset.
     * @ if the array lengths are different.
     * @ if the number of points is less than 2.
     * @ if {@code x} is not sorted in
     * strictly increasing order.
     */
    //override
    public Polynomial_FunctionNewtonForm interpolate(const std::vector<double>& x, const std::vector<double>& y)
         
        {
        /**
         * a[] and c[] are defined in the general formula of Newton form:
         * p(x) = a[0] + a[1](x-c[0]) + a[2](x-c[0])(x-c[1]) + ... +
         *        a[n](x-c[0])(x-c[1])...(x-c[n-1])
         */
        Polynomial_Function_Lagrange_Form.verify_interpolation_array(x, y, true);

        /**
         * When used for interpolation, the Newton form formula becomes
         * p(x) = f[x0] + f[x0,x1](x-x0) + f[x0,x1,x2](x-x0)(x-x1) + ... +
         *        f[x0,x1,...,x[n-1]](x-x0)(x-x1)...(x-x[n-2])
         * Therefore, a[k] = f[x0,x1,...,xk], c[k] = x[k].
         * <p>
         * Note x[], y[], a[] have the same length but c[]'s size is one less.</p>
         */
        const std::vector<double> c = std::vector<double>(x.size()-1];
        System.arraycopy(x, 0, c, 0, c.size());

        const std::vector<double> a = compute_divided_difference(x, y);
        return Polynomial_FunctionNewtonForm(a, c);
    }

    /**
     * Return a copy of the divided difference array.
     * <p>
     * The divided difference array is defined recursively by <pre>
     * f[x0] = f(x0)
     * f[x0,x1,...,xk] = (f[x1,...,xk] - f[x0,...,x[k-1]]) / (xk - x0)
     * </pre>
     * <p>
     * The computational complexity is \(O(n^2)\) where \(n\) is the common
     * length of {@code x} and {@code y}.</p>
     *
     * @param x Interpolating points array.
     * @param y Interpolating values array.
     * @return a fresh copy of the divided difference array.
     * @ if the array lengths are different.
     * @ if the number of points is less than 2.
     * @
     * if {@code x} is not sorted in strictly increasing order.
     */
    protected static std::vector<double> compute_divided_difference(const std::vector<double>& x, const std::vector<double>& y)
         
        {
        Polynomial_Function_Lagrange_Form.verify_interpolation_array(x, y, true);

        const std::vector<double> divdiff = y.clone(); // initialization

        const int n = x.size();
        const std::vector<double> a = double [n];
        a[0] = divdiff[0];
        for (int i{ 1 }; i < n; i++) 
        {
            for (int j{}; j < n-i; j++) 
            {
                const double denominator = x[j+i] - x[j];
                divdiff[j] = (divdiff[j+1] - divdiff[j]) / denominator;
            }
            a[i] = divdiff[0];
        }

        return a;
    }
}


