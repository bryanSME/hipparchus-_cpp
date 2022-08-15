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
//package org.hipparchus.ode;

//import java.util.Collections;
//import java.util.List;

//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
#include <vector>
#include <string>
#include "OrdinaryDifferentialEquation.h"
#include "NamedParameterJacobianProvider.h"

/** Interface expanding {@link Ordinary_Differential_Equation first order
 *  differential equations} in order to compute exactly the Jacobian
 *  matrices for {@link Variational_Equation partial derivatives equations}.
 */
class ODE_Jacobians_Provider : public Ordinary_Differential_Equation, public Named_Parameter_Jacobian_Provider 
{

    /** Compute the Jacobian matrix of ODE with respect to state.
     * @param t current value of the independent <I>time</I> variable
     * @param y array containing the current value of the main state vector
     * @param y_dot array containing the current value of the time derivative of the main state vector
     * @return Jacobian matrix of the ODE w.r.t. the main state vector
     * @exception Math_Illegal_State_Exception if the number of functions evaluations is exceeded
     * @exception  if arrays dimensions do not match equations settings
     */
    std::vector<std::vector<double>> compute_main_state_jacobian(double t, std::vector<double> y, std::vector<double> y_dot);

    /** {@inherit_doc}
     * <p>
     * The default implementation has no parameters at all.
     * </p>
     */
    //override
    std::vector<std::string> get_parameters_names() 
    {
        throw std::exception("Not implemented");
        //return Collections.empty_list();
    }

    /** {@inherit_doc}
     * <p>
     * The default implementation supports no parameters at all.
     * </p>
     */
    //override
    bool is_supported(std::string name) 
    {
        return false;
    }

    /** {@inherit_doc}
     * <p>
     * The default implementation supports no parameters at all.
     * </p>
     */
    //override
    std::vector<double> compute_parameter_jacobian(double t, std::vector<double> y, std::vector<double> y_dot, std::string param_name)
    {
        throw std::exception("not implemented");
        //throw (Localized_ODE_Formats.UNKNOWN_PARAMETER, param_name);
    }

}


