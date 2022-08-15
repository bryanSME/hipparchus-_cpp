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

//import java.util.Arrays;
//import java.util.Hash_Map;
//import java.util.List;
//import java.util.Map;

//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;

/** Wrapper class to compute Jacobian matrices by finite differences for ODE
 *  which do not compute them by themselves.
 *
 */
class Parameter_Jacobian_Wrapper : ODE_Jacobians_Provider 
{

    /** ode base ordinary differential equation for which Jacobians
     * matrices are requested. */
    private const Ordinary_Differential_Equation ode;

    /** Steps for finite difference computation of the jacobian df/dy w.r.t. state. */
    private const std::vector<double> h_y;

    /** Controller to change parameters. */
    private const Parameters_Controller controller;

    /** Steps for finite difference computation of the Jacobian df/dp w.r.t. parameters. */
    private const Map<std::string, Double> h_param;

    /** Wrap a {@link Parameters_Controller} into a {@link Named_Parameter_Jacobian_Provider}.
     * @param ode ode base ordinary differential equation for which Jacobians
     * matrices are requested
     * @param h_y step used for finite difference computation with respect to state vector
     * @param controller controller to change parameters
     * @param params_and_steps parameters and steps to compute the Jacobians df/dp
     * @see Parameters_Controller#set_parameter_step(std::string, double)
     */
    Parameter_Jacobian_Wrapper(const Ordinary_Differential_Equation ode, const std::vector<double> h_y, const Parameters_Controller controller, const Parameter_Configuration[] params_and_steps) 
    {
        this.ode        = ode;
        this.h_y         = h_y.clone();
        this.controller = controller;
        this.h_param     = Hash_Map<>();

        // set up parameters for jacobian computation
        for (const Parameter_Configuration param : params_and_steps) 
        {
            const std::string name = param.get_parameter_name();
            if (controller.is_supported(name)) 
            {
                h_param.put(name, param.get_h_p());
            }
        }
    }

    /** Get the underlying ode.
     * @return underlying ode
     */
    public Ordinary_Differential_Equation get_ode() 
    {
        return ode;
    }

    /** {@inherit_doc} */
    //override
    public int get_dimension() 
    {
        return ode.get_dimension();
    }

    /** {@inherit_doc} */
    //override
    public std::vector<double> compute_derivatives(double t, std::vector<double> y)
        , Math_Illegal_State_Exception 
        {
        return ode.compute_derivatives(t, y);
    }

    /** {@inherit_doc} */
    //override
    public std::vector<std::vector<double>> compute_main_state_jacobian(double t, std::vector<double> y, std::vector<double> y_dot)
        , Math_Illegal_State_Exception 
        {

        const int n = ode.get_dimension();
        const std::vector<std::vector<double>> d_fd_y = std::vector<double>(n][n];
        for (int j{}; j < n; ++j) 
        {
            const double saved_yj = y[j];
            y[j] += h_y[j];
            const std::vector<double> tmp_dot = ode.compute_derivatives(t, y);
            for (int i{}; i < n; ++i) 
            {
                d_fd_y[i][j] = (tmp_dot[i] - y_dot[i]) / h_y[j];
            }
            y[j] = saved_yj;
        }
        return d_fd_y;
    }

    /** {@inherit_doc} */
    //override
    public List<std::string> get_parameters_names() 
    {
        return controller.get_parameters_names();
    }

    /** {@inherit_doc} */
    //override
    public bool is_supported(std::string name) 
    {
        return controller.is_supported(name);
    }

    /** {@inherit_doc} */
    //override
    public std::vector<double> compute_parameter_jacobian(const double t, const std::vector<double> y, const std::vector<double> y_dot, const std::string param_name)
        , Math_Illegal_State_Exception 
        {

        const int n = ode.get_dimension();
        const std::vector<double> df_dp = std::vector<double>(n];
        if (controller.is_supported(param_name)) 
        {

            // compute the jacobian df/dp w.r.t. parameter
            const double p  = controller.get_parameter(param_name);
            const double h_p = h_param.get(param_name);
            controller.set_parameter(param_name, p + h_p);
            const std::vector<double> tmp_dot = ode.compute_derivatives(t, y);
            for (int i{}; i < n; ++i) 
            {
                df_dp[i] = (tmp_dot[i] - y_dot[i]) / h_p;
            }
            controller.set_parameter(param_name, p);
        }
else 
        {
            Arrays.fill(df_dp, 0, n, 0.0);
        }

        return df_dp;

    }

}


