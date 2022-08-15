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

#include <vector>

//import org.hipparchus.complex.std::complex<double>;

/** This interface represents a first order differential equations set for {@link std::complex<double> complex state}.
 *
 * @see Ordinary_Differential_Equation
 * @see std::complex<double>_ODE_Converter
 * @since 1.4
 *
 */
class std::complex<double>Ordinary_Differential_Equation
{

public:
    /** Get the dimension of the problem.
     * @return dimension of the problem
     */
    int get_dimension();

    /** Initialize equations at the start of an ODE integration.
     * <p>
     * This method is called once at the start of the integration. It
     * may be used by the equations to initialize some internal data
     * if needed.
     * </p>
     * <p>
     * The default implementation does nothing.
     * </p>
     * @param t0 value of the independent <I>time</I> variable at integration start
     * @param y0 array containing the value of the state vector at integration start
     * @param const_time target time for the integration
     */
    void init(const double& t0, std::vector<std::complex<double>>& y0, const double& const_time)
    {
        // do nothing by default
    }

    /** Get the current time derivative of the state vector.
     * @param t current value of the independent <I>time</I> variable
     * @param y array containing the current value of the state vector
     * @return time derivative of the state vector
     */
    std::vector<std::complex<double>>compute_derivatives(const double& t, const std::vector<std::complex<double>>& y);
};