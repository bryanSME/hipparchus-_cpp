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

//package org.hipparchus.ode.nonstiff;

//import java.util.Arrays;
//import java.util.Hash_Map;
//import java.util.Map;

//import org.hipparchus.Field;
//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.linear.Array2DRowField_Matrix;
//import org.hipparchus.linear.ArrayField_Vector;
//import org.hipparchus.linear.FieldDecomposition_Solver;
//import org.hipparchus.linear.FieldLU_Decomposition;
//import org.hipparchus.linear.Field_Matrix;
//import org.hipparchus.util.Math_Arrays;
#include <type_traits>
#include "../../core/CalculusFieldElement.h"

/** Transformer to Nordsieck vectors for Adams integrators.
 * <p>This class is used by {@link AdamsBashforth_integrator Adams-Bashforth} and
 * {@link Adams_moultonIntegrator Adams-Moulton} integrators to convert between
 * classical representation with several previous first derivatives and Nordsieck
 * representation with higher order scaled derivatives.</p>
 *
 * <p>We define scaled derivatives s<sub>i</sub>(n) at step n as:
 * <pre>
 * s<sub>1</sub>(n) = h y'<sub>n</sub> for first derivative
 * s<sub>2</sub>(n) = h<sup>2</sup>/2 y''<sub>n</sub> for second derivative
 * s<sub>3</sub>(n) = h<sup>3</sup>/6 y'''<sub>n</sub> for third derivative
 * ...
 * s<sub>k</sub>(n) = h<sup>k</sup>/k! y<sup>(k)</sup><sub>n</sub> for k<sup>th</sup> derivative
 * </pre></p>
 *
 * <p>With the previous definition, the classical representation of multistep methods
 * uses first derivatives only, i.e. it handles y<sub>n</sub>, s<sub>1</sub>(n) and
 * q<sub>n</sub> where q<sub>n</sub> is defined as:
 * <pre>
 *   q<sub>n</sub> = [ s<sub>1</sub>(n-1) s<sub>1</sub>(n-2) ... s<sub>1</sub>(n-(k-1)) ]<sup>T</sup>
 * </pre>
 * (we omit the k index in the notation for clarity).</p>
 *
 * <p>Another possible representation uses the Nordsieck vector with
 * higher degrees scaled derivatives all taken at the same step, i.e it handles y<sub>n</sub>, * s<sub>1</sub>(n) and r<sub>n</sub>) where r<sub>n</sub> is defined as:
 * <pre>
 * r<sub>n</sub> = [ s<sub>2</sub>(n), s<sub>3</sub>(n) ... s<sub>k</sub>(n) ]<sup>T</sup>
 * </pre>
 * (here again we omit the k index in the notation for clarity)
 * </p>
 *
 * <p>Taylor series formulas show that for any index offset i, s<sub>1</sub>(n-i) can be
 * computed from s<sub>1</sub>(n), s<sub>2</sub>(n) ... s<sub>k</sub>(n), the formula being exact
 * for degree k polynomials.
 * <pre>
 * s<sub>1</sub>(n-i) = s<sub>1</sub>(n) + &sum;<sub>j&gt;0</sub> (j+1) (-i)<sup>j</sup> s<sub>j+1</sub>(n)
 * </pre>
 * The previous formula can be used with several values for i to compute the transform between
 * classical representation and Nordsieck vector at step end. The transform between r<sub>n</sub>
 * and q<sub>n</sub> resulting from the Taylor series formulas above is:
 * <pre>
 * q<sub>n</sub> = s<sub>1</sub>(n) u + P r<sub>n</sub>
 * </pre>
 * where u is the [ 1 1 ... 1 ]<sup>T</sup> vector and P is the (k-1)&times;(k-1) matrix built
 * with the (j+1) (-i)<sup>j</sup> terms with i being the row number starting from 1 and j being
 * the column number starting from 1:
 * <pre>
 *        [  -2   3   -4    5  ... ]
 *        [  -4  12  -32   80  ... ]
 *   P =  [  -6  27 -108  405  ... ]
 *        [  -8  48 -256 1280  ... ]
 *        [          ...           ]
 * </pre></p>
 *
 * <p>Changing -i into +i in the formula above can be used to compute a similar transform between
 * classical representation and Nordsieck vector at step start. The resulting matrix is simply
 * the absolute value of matrix P.</p>
 *
 * <p>For {@link AdamsBashforth_integrator Adams-Bashforth} method, the Nordsieck vector
 * at step n+1 is computed from the Nordsieck vector at step n as follows:
 * <ul>
 *   <li>y<sub>n+1</sub> = y<sub>n</sub> + s<sub>1</sub>(n) + u<sup>T</sup> r<sub>n</sub></li>
 *   <li>s<sub>1</sub>(n+1) = h f(t<sub>n+1</sub>, y<sub>n+1</sub>)</li>
 *   <li>r<sub>n+1</sub> = (s<sub>1</sub>(n) - s<sub>1</sub>(n+1)) P<sup>-1</sup> u + P<sup>-1</sup> A P r<sub>n</sub></li>
 * </ul>
 * where A is a rows shifting matrix (the lower left part is an identity matrix):
 * <pre>
 *        [ 0 0   ...  0 0 | 0 ]
 *        [ ---------------+---]
 *        [ 1 0   ...  0 0 | 0 ]
 *    A = [ 0 1   ...  0 0 | 0 ]
 *        [       ...      | 0 ]
 *        [ 0 0   ...  1 0 | 0 ]
 *        [ 0 0   ...  0 1 | 0 ]
 * </pre></p>
 *
 * <p>For {@link Adams_moultonIntegrator Adams-Moulton} method, the predicted Nordsieck vector
 * at step n+1 is computed from the Nordsieck vector at step n as follows:
 * <ul>
 *   <li>Y<sub>n+1</sub> = y<sub>n</sub> + s<sub>1</sub>(n) + u<sup>T</sup> r<sub>n</sub></li>
 *   <li>S<sub>1</sub>(n+1) = h f(t<sub>n+1</sub>, Y<sub>n+1</sub>)</li>
 *   <li>R<sub>n+1</sub> = (s<sub>1</sub>(n) - s<sub>1</sub>(n+1)) P<sup>-1</sup> u + P<sup>-1</sup> A P r<sub>n</sub></li>
 * </ul>
 * From this predicted vector, the corrected vector is computed as follows:
 * <ul>
 *   <li>y<sub>n+1</sub> = y<sub>n</sub> + S<sub>1</sub>(n+1) + [ -1 +1 -1 +1 ... &plusmn;1 ] r<sub>n+1</sub></li>
 *   <li>s<sub>1</sub>(n+1) = h f(t<sub>n+1</sub>, y<sub>n+1</sub>)</li>
 *   <li>r<sub>n+1</sub> = R<sub>n+1</sub> + (s<sub>1</sub>(n+1) - S<sub>1</sub>(n+1)) P<sup>-1</sup> u</li>
 * </ul>
 * where the upper case Y<sub>n+1</sub>, S<sub>1</sub>(n+1) and R<sub>n+1</sub> represent the
 * predicted states whereas the lower case y<sub>n+1</sub>, s<sub>n+1</sub> and r<sub>n+1</sub>
 * represent the corrected states.</p>
 *
 * <p>We observe that both methods use similar update formulas. In both cases a P<sup>-1</sup>u
 * vector and a P<sup>-1</sup> A P matrix are used that do not depend on the state, * they only depend on k. This class handles these transformations.</p>
 *
 * @param <T> the type of the field elements
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = NULLptr>
class Adams_Nordsieck_Field_Transformer 
{
private:
    /** Cache for already computed coefficients. */
    static const Map<Integer, Map<Field<? extends Calculus_Field_Element<?>>, Adams_Nordsieck_Field_Transformer<? extends Calculus_Field_Element<?>>>> CACHE = Hash_Map<>();

    /** Field to which the time and state vector elements belong. */
    const Field<T> field;

    /** Update matrix for the higher order derivatives h<sup>2</sup>/2 y'', h<sup>3</sup>/6 y''' ... */
    const Array2DRowField_Matrix<T> update;

    /** Update coefficients of the higher order derivatives wrt y'. */
    const std::vector<T> c1;

    /** Simple constructor.
     * @param field field to which the time and state vector elements belong
     * @param n number of steps of the multistep method
     * (excluding the one being computed)
     */
    Adams_Nordsieck_Field_Transformer(const Field<T> field, const int& n) 
    {

        this.field = field;
        const int rows = n - 1;

        // compute coefficients
        Field_Matrix<T> big_p = build_p(rows);
        FieldDecomposition_Solver<T> p_solver =
            FieldLU_Decomposition<T>(big_p).get_solver();

        std::vector<T> u = Math_Arrays::build_array(field, rows);
        Arrays.fill(u, field.get_one());
        c1 = p_solver.solve(new ArrayField_Vector<T>(u, false)).to_array();

        // update coefficients are computed by combining transform from
        // Nordsieck to multistep, then shifting rows to represent step advance
        // then applying inverse transform
        std::vector<std::vector<T>> shifted_p = big_p.get_data();
        for (int i = shifted_p.size() - 1; i > 0; --i) 
        {
            // shift rows
            shifted_p[i] = shifted_p[i - 1];
        }
        shifted_p[0] = Math_Arrays::build_array(field, rows);
        Arrays.fill(shifted_p[0], field.get_zero());
        update = Array2DRowField_Matrix<>(p_solver.solve(new Array2DRowField_Matrix<T>(shifted_p, false)).get_data());
    }

    /** Build the P matrix.
    * <p>The P matrix general terms are shifted (j+1) (-i)<sup>j</sup> terms
    * with i being the row number starting from 1 and j being the column
    * number starting from 1:
    * <pre>
    *        [  -2   3   -4    5  ... ]
    *        [  -4  12  -32   80  ... ]
    *   P =  [  -6  27 -108  405  ... ]
    *        [  -8  48 -256 1280  ... ]
    *        [          ...           ]
    * </pre></p>
    * @param rows number of rows of the matrix
    * @return P matrix
    */
    private Field_Matrix<T> build_p(const int& rows)
    {

        const std::vector<std::vector<T>> p_data = Math_Arrays::build_array(field, rows, rows);

        for (int i{ 1 }; i <= p_data.size(); ++i)
        {
            // build the P matrix elements from Taylor series formulas
            auto pI = p_data[i - 1];
            const int factor = -i;
            T aj = field.get_zero().add(factor);
            for (int j{ 1 }; j <= pI.size(); ++j)
            {
                pI[j - 1] = aj.multiply(j + 1);
                aj = aj.multiply(factor);
            }
        }

        return Array2DRowField_Matrix<T>(p_data, false);

    }

public:

    /** Get the Nordsieck transformer for a given field and number of steps.
     * @param field field to which the time and state vector elements belong
     * @param n_steps number of steps of the multistep method
     * (excluding the one being computed)
     * @return Nordsieck transformer for the specified field and number of steps
     * @param <T> the type of the field elements
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = NULLptr>
    static Adams_Nordsieck_Field_Transformer<T> get_instance(const Field<T>& field, const int& n_steps) 
    { 
        // NOPMD - PMD false positive
        synchronized(CACHE) 
        {
            Map<Field<? extends Calculus_Field_Element<?>>, Adams_Nordsieck_Field_Transformer<? extends Calculus_Field_Element<?>>> map = CACHE.get(n_steps);
            if (map == NULL) 
            {
                map = Hash_Map<>();
                CACHE.put(n_steps, map);
            }
            //@Suppress_Warnings("unchecked")
            Adams_Nordsieck_Field_Transformer<T> t = (Adams_Nordsieck_Field_Transformer<T>) map.get(field);
            if (t == NULL) 
            {
                t = Adams_Nordsieck_Field_Transformer<>(field, n_steps);
                map.put(field, t);
            }
            return t;

        }
    }

    /** Initialize the high order scaled derivatives at step start.
     * @param h step size to use for scaling
     * @param t first steps times
     * @param y first steps states
     * @param y_dot first steps derivatives
     * @return Nordieck vector at start of first step (h<sup>2</sup>/2 y''<sub>n</sub>, * h<sup>3</sup>/6 y'''<sub>n</sub> ... h<sup>k</sup>/k! y<sup>(k)</sup><sub>n</sub>)
     */
    Array2DRowField_Matrix<T> initialize_high_order_derivatives(const T h, const std::vector<T> t, const std::vector<std::vector<T>> y, const std::vector<std::vector<T>> y_dot) 
    {

        // using Taylor series with di = ti - t0, we get:
        //  y(ti)  - y(t0)  - di y'(t0) =   di^2 / h^2 s2 + ... +   di^k     / h^k sk + O(h^k)
        //  y'(ti) - y'(t0)             = 2 di   / h^2 s2 + ... + k di^(k-1) / h^k sk + O(h^(k-1))
        // we write these relations for i = 1 to i= 1+n/2 as a set of n + 2 linear
        // equations depending on the Nordsieck vector [s2 ... sk rk], so s2 to sk correspond
        // to the appropriately truncated Taylor expansion, and rk is the Taylor remainder.
        // The goal is to have s2 to sk as accurate as possible considering the fact the sum is
        // truncated and we don't want the error terms to be included in s2 ... sk, so we need
        // to solve also for the remainder
        const std::vector<std::vector<T>> a     = Math_Arrays::build_array(field, c1.size() + 1, c1.size() + 1);
        const std::vector<std::vector<T>> b     = Math_Arrays::build_array(field, c1.size() + 1, y[0].size());
        const std::vector<T>   y0    = y[0];
        const std::vector<T>   y_dot_0 = y_dot[0];
        for (int i{ 1 }; i < y.size(); ++i) 
        {

            const T di    = t[i].subtract(t[0]);
            const T ratio = di.divide(h);
            T dik_m1_ohk    = h.reciprocal();

            // linear coefficients of equations
            // y(ti) - y(t0) - di y'(t0) and y'(ti) - y'(t0)
            const std::vector<T> aI    = a[2 * i - 2];
            const std::vector<T> a_dot_i = (2 * i - 1) < a.size() ? a[2 * i - 1] : NULL;
            for (int j{}; j < aI.size(); ++j) 
            {
                dik_m1_ohk = dik_m1_ohk.multiply(ratio);
                aI[j]    = di.multiply(dik_m1_ohk);
                if (a_dot_i != NULL) 
                {
                    a_dot_i[j]  = dik_m1_ohk.multiply(j + 2);
                }
            }

            // expected value of the previous equations
            const std::vector<T> y_i    = y[i];
            const std::vector<T> y_dot_i = y_dot[i];
            const std::vector<T> bI    = b[2 * i - 2];
            const std::vector<T> b_dot_i = (2 * i - 1) < b.size() ? b[2 * i - 1] : NULL;
            for (int j{}; j < y_i.size(); ++j) 
            {
                bI[j]    = y_i[j].subtract(y0[j]).subtract(di.multiply(y_dot_0[j]));
                if (b_dot_i != NULL) 
                {
                    b_dot_i[j] = y_dot_i[j].subtract(y_dot_0[j]);
                }
            }

        }

        // solve the linear system to get the best estimate of the Nordsieck vector [s2 ... sk], // with the additional terms s(k+1) and c grabbing the parts after the truncated Taylor expansion
        const FieldLU_Decomposition<T> decomposition = FieldLU_Decomposition<>(new Array2DRowField_Matrix<T>(a, false));
        const Field_Matrix<T> x = decomposition.get_solver().solve(new Array2DRowField_Matrix<T>(b, false));

        // extract just the Nordsieck vector [s2 ... sk]
        const Array2DRowField_Matrix<T> truncated_x =
                        Array2DRowField_Matrix<>(field, x.get_row_dimension() - 1, x.get_column_dimension());
        for (int i{}; i < truncated_x.get_row_dimension(); ++i) 
        {
            for (int j{}; j < truncated_x.get_column_dimension(); ++j) 
            {
                truncated_x.set_entry(i, j, x.get_entry(i, j));
            }
        }
        return truncated_x;
    }

    /** Update the high order scaled derivatives for Adams integrators (phase 1).
     * <p>The complete update of high order derivatives has a form similar to:
     * <pre>
     * r<sub>n+1</sub> = (s<sub>1</sub>(n) - s<sub>1</sub>(n+1)) P<sup>-1</sup> u + P<sup>-1</sup> A P r<sub>n</sub>
     * </pre>
     * this method computes the P<sup>-1</sup> A P r<sub>n</sub> part.</p>
     * @param high_order high order scaled derivatives
     * (h<sup>2</sup>/2 y'', ... h<sup>k</sup>/k! y(k))
     * @return updated high order derivatives
     * @see #update_high_order_derivatives_phase_2(Calculus_Field_Element[], Calculus_Field_Element[], Array2DRowField_Matrix)
     */
    Array2DRowField_Matrix<T> update_high_order_derivatives_phase_1(const Array2DRowField_Matrix<T> high_order) 
    {
        return update.multiply(high_order);
    }

    /** Update the high order scaled derivatives Adams integrators (phase 2).
     * <p>The complete update of high order derivatives has a form similar to:
     * <pre>
     * r<sub>n+1</sub> = (s<sub>1</sub>(n) - s<sub>1</sub>(n+1)) P<sup>-1</sup> u + P<sup>-1</sup> A P r<sub>n</sub>
     * </pre>
     * this method computes the (s<sub>1</sub>(n) - s<sub>1</sub>(n+1)) P<sup>-1</sup> u part.</p>
     * <p>Phase 1 of the update must already have been performed.</p>
     * @param start first order scaled derivatives at step start
     * @param end first order scaled derivatives at step end
     * @param high_order high order scaled derivatives, will be modified
     * (h<sup>2</sup>/2 y'', ... h<sup>k</sup>/k! y(k))
     * @see #update_high_order_derivatives_phase_1(Array2DRowField_Matrix)
     */
    void update_high_order_derivatives_phase_2(const std::vector<T> start, const std::vector<T> end, const Array2DRowField_Matrix<T> high_order) 
    {
        const std::vector<std::vector<T>> data = high_order.get_data_ref();
        for (int i{}; i < data.size(); ++i) 
        {
            const std::vector<T> data_i = data[i];
            const T c1I = c1[i];
            for (int j{}; j < data_i.size(); ++j) 
            {
                data_i[j] = data_i[j].add(c1I.multiply(start[j].subtract(end[j])));
            }
        }
    }
};