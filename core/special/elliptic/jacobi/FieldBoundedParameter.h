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
//import org.hipparchus.special.elliptic.legendre.Legendre_Elliptic_Integral;
#include <type_traits>
#include "../../../CalculusFieldElement.hpp"

/** Algorithm for computing the principal Jacobi functions for parameter m in [0; 1].
 * @param <T> the type of the field elements
 * @since 2.0
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
class FieldBounded_Parameter : public Field_Jacobi_Elliptic<T>
{
    /** Jacobi \xce\xb8 functions. */
    private const Field_Jacobi_Theta<T> jacobi_theta;

    /** Value of Jacobi \xce\xb8 functions at origin. */
    private const Field_Theta<T> t0;

    /** Scaling factor. */
    private const T scaling;

    /** Simple constructor.
     * @param m parameter of the Jacobi elliptic function
     */
    FieldBounded_Parameter(const T m)
    {

        super(m);

        // compute nome
        const T q = Legendre_Elliptic_Integral.nome(m);

        // prepare underlying Jacobi \xce\xb8 functions
        this.jacobi_theta = Field_Jacobi_Theta<>(q);
        this.t0 = jacobi_theta.values(m.get_field().get_zero());
        this.scaling = Legendre_Elliptic_Integral.big_k(m).reciprocal().multiply(m.get_pi().multiply(0.5));

    }

    /** {@inherit_doc}
     * <p>
     * The algorithm for evaluating the functions is based on {@link Field_Jacobi_Theta
     * Jacobi theta functions}.
     * </p>
     */
     //override
    public Field_Copolar_N<T> values_n(T u)
    {

        // evaluate Jacobi \xce\xb8 functions at argument
        const Field_Theta<T> t_z = jacobi_theta.values(u.multiply(scaling));

        // convert to Jacobi elliptic functions
        const T sn = t0.theta3().multiply(t_z.theta1()).divide(t0.theta2().multiply(t_z.theta4()));
        const T cn = t0.theta4().multiply(t_z.theta2()).divide(t0.theta2().multiply(t_z.theta4()));
        const T dn = t0.theta4().multiply(t_z.theta3()).divide(t0.theta3().multiply(t_z.theta4()));

        return Field_Copolar_N<>(sn, cn, dn);

    }

};