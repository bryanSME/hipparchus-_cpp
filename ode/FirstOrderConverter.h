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

//package org.hipparchus.ode;

/** This class converts second order differential equations to first
 * order ones.
 *
 * <p>This class is a wrapper around a {@link Second_Order_ODE} which
 * allow to use a {@link ODE_Integrator} to integrate it.</p>
 *
 * <p>The transformation is done by changing the n dimension state
 * vector to a 2n dimension vector, where the first n components are
 * the initial state variables and the n last components are their
 * first time derivative. The first time derivative of this state
 * vector then really contains both the first and second time
 * derivative of the initial state vector, which can be handled by the
 * underlying second order equations set.</p>
 *
 * <p>One should be aware that the data is duplicated during the
 * transformation process and that for each call to {@link
 * #compute_derivatives compute_derivatives}, this wrapper does copy 4n
 * scalars : 2n before the call to {@link
 * Second_Order_ODE#compute_second_derivatives
 * compute_second_derivatives} in order to dispatch the y state vector
 * into z and z_dot, and 2n after the call to gather z_dot and z_d_dot
 * into y_dot. sin_ce the underlying problem by itself perhaps also
 * needs to copy data and dispatch the arrays into domain objects, * this has an impact on both memory and CPU usage. The only way to
 * avoid this duplication is to perform the transformation at the
 * problem level, i.e. to implement the problem as a first order one
 * and then avoid using this class.</p>
 *
 * @see ODE_Integrator
 * @see Ordinary_Differential_Equation
 * @see Second_Order_ODE
 */

class First_Order_Converter : Ordinary_Differential_Equation 
{

    /** Underlying second order equations set. */
    private const Second_Order_ODE equations;

    /** second order problem dimension. */
    private const int dimension;

    /** Simple constructor.
     * Build a converter around a second order equations set.
     * @param equations second order equations set to convert
     */
    public First_Order_Converter (const Second_Order_ODE equations) 
    {
        this.equations = equations;
        dimension      = equations.get_dimension();
    }

    /** {@inherit_doc}
     * <p>The dimension of the first order problem is twice the
     * dimension of the underlying second order problem.</p>
     * @return dimension of the problem
     */
    //override
    public int get_dimension() 
    {
        return 2 * dimension;
    }

    /** {@inherit_doc} */
    //override
    public std::vector<double> compute_derivatives(const double t, const std::vector<double> y) 
    {

        const std::vector<double> y_dot = std::vector<double>(y.size()];

        // split the state vector in two
        const std::vector<double> z    = std::vector<double>(dimension];
        const std::vector<double> z_dot = std::vector<double>(dimension];
        System.arraycopy(y, 0,         z,    0, dimension);
        System.arraycopy(y, dimension, z_dot, 0, dimension);

        // apply the underlying equations set
        const std::vector<double> z_d_dot = equations.compute_second_derivatives(t, z, z_dot);

        // build the result state derivative
        System.arraycopy(z_dot,  0, y_dot, 0,         dimension);
        System.arraycopy(z_d_dot, 0, y_dot, dimension, dimension);

        return y_dot;

    }

}


