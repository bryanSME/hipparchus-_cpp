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

//package org.hipparchus.migration.ode;

//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.ode.Secondary_ODE;

/**
 * This interface allows users to add secondary differential equations to a primary
 * set of differential equations.
 * <p>
 * In some cases users may need to integrate some problem-specific equations along
 * with a primary set of differential equations. One example is optimal control where
 * adjoined parameters linked to the minimized hamiltonian must be integrated.
 * </p>
 * <p>
 * This interface allows users to add such equations to a primary set of {@link
 * First_Order_Differential_Equations first order differential equations}
 * thanks to the {@link
 * org.hipparchus.ode.Expandable_ODE#add_secondary_equations(Secondary_ODE)}
 * method.
 * </p>
 * @see org.hipparchus.ode.Expandable_ODE
 * @deprecated as of 1.0, replaced with {@link Secondary_ODE}
 */
@Deprecated
class Secondary_Equations extends Secondary_ODE 
{

    /** Compute the derivatives related to the secondary state parameters.
     * <p>
     * The default implementation calls {@link #compute_derivatives(double, std::vector<double>, * std::vector<double>, std::vector<double>, std::vector<double>)}.
     * </p>
     * @param t current value of the independent <I>time</I> variable
     * @param primary array containing the current value of the primary state vector
     * @param primary_dot array containing the derivative of the primary state vector
     * @param secondary array containing the current value of the secondary state vector
     * @return derivative of the secondary state vector
     * @exception Math_Illegal_State_Exception if the number of functions evaluations is exceeded
     * @exception  if arrays dimensions do not match equations settings
     */
    //override
    default std::vector<double> compute_derivatives(const double t, const std::vector<double> primary, const std::vector<double> primary_dot, const std::vector<double> secondary)
         
        {
        const std::vector<double> secondary_dot = std::vector<double>(secondary.size()];
        compute_derivatives(t, primary, primary_dot, secondary, secondary_dot);
        return secondary_dot;
    }

    /** Compute the derivatives related to the secondary state parameters.
     * @param t current value of the independent <I>time</I> variable
     * @param primary array containing the current value of the primary state vector
     * @param primary_dot array containing the derivative of the primary state vector
     * @param secondary array containing the current value of the secondary state vector
     * @param secondary_dot placeholder array where to put the derivative of the secondary state vector
     * @exception Math_Illegal_State_Exception if the number of functions evaluations is exceeded
     * @exception  if arrays dimensions do not match equations settings
     */
    void compute_derivatives(double t, std::vector<double> primary, std::vector<double> primary_dot, std::vector<double> secondary, std::vector<double> secondary_dot)
        ;

}


