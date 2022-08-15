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

//import java.lang.reflect.Array;

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
#include <string>
#include <vector>
#include "ODEJacobiansProvider.h"


/**
 * This class defines a set of {@link Secondary_ODE secondary equations} to
 * compute the global Jacobian matrices with respect to the initial state
 * vector and, if any, to some parameters of the primary ODE set.
 * <p>
 * The primary set of ODE for which Jaobian matrices are requested may be:
 * <ul>
 * <li>a full-fledged {@link ODE_Jacobians_Provider} that computes by itself
 * both the ODE and its local partial derivatives,</li>
 * <li>a simple {@link Ordinary_Differential_Equation} which must therefore
 * be completed with a finite differences configuration to compute local
 * partial derivatives (so-called internal differentiation).</li>
 * </ul>
 * </p>
 * <p>
 * As the variational equation automatically inserts {@link
 * Expandable_ODE#add_secondary_equations(Secondary_ODE) secondary differential
 * equations}, in the {@link Expandable_ODE expandable ODE}, data for
 * initial state must also be inserted before integration and matrices
 * result must be extracted after integration. This implies a precise
 * scheduling of the calls to the various methods of this class. The
 * proper scheduling is the following one:
 * </p>
 * <pre>
 *   // set up equations
 *   ODE_Jacobians_Provider jode       = My_ODE(...);
 *   Expandable_ODE        expandable = Expandable(jode);
 *   Variational_Equation  ve         = Variational_Equation(expandable, jode);
 *
 *   // set up initial state
 *   ODE_State init_without_derivatives = ODE_State(t0, y0);
 *   ve.set_initial_main_state_jacobian(dYdY0); // only needed if the default identity matrix is not suitable
 *   ve.set_initial_parameter_jacobian(name, dYdP); // only needed if the default zero matrix is not suitable
 *   ODE_State init_with_derivatives = ve.set_up_initial_state(init_without_derivatives);
 *
 *   // perform integration on the expanded equations with the expanded initial state
 *   ODE_State_And_Derivative const_state = integrator.integrate(expandable, init_with_derivatives, const_t);
 *
 *   // extract Jacobian matrices
 *   dYdY0 = ve.extract_main_set_jacobian(const_state);
 *   dYdP  = ve.extract_parameter_jacobian(const_state, name);
 * </pre>
 * <p>
 * The most important part is to not forget to call {@link #set_up_initial_state(ODE_State)} to add
 * the secondary state with the initial matrices to the {@link ODE_State} used in the
 * {@link ODE_Integrator#integrate(Expandable_ODE, ODE_State, double) integrate} method.
 * Forgetting to do this and passing only a {@link ODE_State} without the secondary state
 * set up will trigger an error as the state vector will not have the correct dimension.
 * </p>
 *
 * @see Expandable_ODE
 * @see ODE_Jacobians_Provider
 * @see Ordinary_Differential_Equation
 * @see Named_Parameter_Jacobian_Provider
 * @see Parameters_Controller
 *
 */
class Variational_Equation 
{
private:
    /** ODE with Jacobian computation skill. */
    const ODE_Jacobians_Provider my_jode;

    /** Expandable first order differential equation. */
    const Expandable_ODE my_expandable;

    /** Index of the instance in the expandable set. */
    const int index;

    /** State and parameters Jacobian matrices in a row. */
    std::vector<double> my_matrices_data;

    /** Check array dimensions.
     * @param expected expected dimension
     * @param array (may be NULL if expected is 0)
     * @ if the array dimension does not match the expected one
     */
    void check_dimension(const int expected, const Object array)
    {
        int array_dimension = (array == NULL) ? 0 : Array.get_length(array);
        if (array_dimension != expected)
        {
            throw (Localized_Core_Formats.DIMENSIONS_MISMATCH, array_dimension, expected);
        }
    }

    /** Local implementation of secondary equations.
     * <p>
     * This class is an inner class to ensure proper scheduling of calls
     * by forcing the use of {@link Variational_Equation#register_variational_equations(Expandable_ODE)}.
     * </p>
     */
    class JacobiansSecondary_ODE : public Secondary_ODE
    {
    public:
        /** {@inherit_doc} */
        //override
        int get_dimension()
        {
            return my_jode.get_dimension() * (my_jode.get_dimension() + my_jode.get_parameters_names().size());
        }

        /** {@inherit_doc} */
        //override
        std::vector<double> compute_derivatives(const double t, const std::vector<double> y, const std::vector<double> y_dot, const std::vector<double> z)
        {

            const std::vector<double> z_dot = std::vector<double>(z.size()];

            // variational equations:
            // from d[dy/dt]/dy0 and d[dy/dt]/dp to d[dy/dy0]/dt and d[dy/dp]/dt

            // compute Jacobian matrix with respect to primary state
            std::vector<std::vector<double>> d_fd_y = my_jode.compute_main_state_jacobian(t, y, y_dot);

            // Dispatch Jacobian matrix in the compound secondary state vector
            for (int i{}; i < my_jode.get_dimension(); ++i)
            {
                const std::vector<double> d_fd_yi = d_fd_y[i];
                for (int j{}; j < my_jode.get_dimension(); ++j)
                {
                    double s = 0;
                    const int start_index = j;
                    int z_index = start_index;
                    for (const int& l = 0; l < my_jode.get_dimension(); ++l)
                    {
                        s += d_fd_yi[l] * z[z_index];
                        z_index += my_jode.get_dimension();
                    }
                    z_dot[start_index + i * my_jode.get_dimension()] = s;
                }
            }

            // compute Jacobian matrices with respect to parameters
            int start_index = my_jode.get_dimension() * my_jode.get_dimension();
            for (const std::string name : my_jode.get_parameters_names())
            {
                const std::vector<double> df_dp = my_jode.compute_parameter_jacobian(t, y, y_dot, name);
                for (int i{}; i < my_jode.get_dimension(); ++i)
                {
                    const std::vector<double> d_fd_yi = d_fd_y[i];
                    int z_index = start_index;
                    double s = df_dp[i];
                    for (const int& l = 0; l < my_jode.get_dimension(); ++l)
                    {
                        s += d_fd_yi[l] * z[z_index];
                        z_index++;
                    }
                    z_dot[start_index + i] = s;
                }
                start_index += my_jode.get_dimension();
            }

            return z_dot;

            }
    }


public:
    /** Build variational equation using finite differences for local
     * partial derivatives.
     * @param expandable expandable set into which variational equations should be registered
     * @param ode base ordinary differential equation for which Jacobians
     * matrices are requested
     * @param h_y step used for finite difference computation with respect to state vector
     * @param controller controller to change parameters
     * @param params_and_steps parameters and steps to compute the Jacobians df/dp
     * @exception Mismatched_Equations if the primary set of the expandable set does
     * not match the {@code ode}
     */
    Variational_Equation(const Expandable_ODE& expandable, const Ordinary_Differential_Equation& ode, const std::vector<double>& h_y, const Parameters_Controller& controller, const Parameter_Configuration ... params_and_steps) 
    {
        Variational_Equation(expandable, Parameter_Jacobian_Wrapper(ode, h_y, controller, params_and_steps));
    }

    /** Build variational equation using analytical local partial derivatives.
     * <p>
     * Parameters must belong to the supported ones given by {@link
     * Parameterizable#get_parameters_names()}, so the primary set of differential
     * equations must be {@link Parameterizable}.
     * </p>
     * <p>Note that each selection clears the previous selected parameters.</p>
     *
     * @param expandable expandable set into which variational equations should be registered
     * @param jode the primary first order differential equations set to extend
     * @exception Mismatched_Equations if the primary set of the expandable set does
     * not match the {@code ode}
     */
    Variational_Equation(const Expandable_ODE& expandable, const ODE_Jacobians_Provider& jode)
    {
        // safety checks
        Ordinary_Differential_Equation ode;
        if (jode instanceof Parameter_Jacobian_Wrapper)
        {
            ode = ((Parameter_Jacobian_Wrapper)jode).get_ode();
        }
        else
        {
            ode = jode;
        }
        if (expandable.get_primary() != ode)
        {
            throw Mismatched_Equations();
        }

        this.jode = jode;
        this.expandable = expandable;
        this.index = expandable.add_secondary_equations(new JacobiansSecondary_ODE());

        // set the default initial state Jacobian to the identity
        // and the default initial parameters Jacobian to the NULL matrix
        matrices_data = std::vector<double>((jode.get_dimension() + jode.get_parameters_names().size()) * jode.get_dimension()];
        for (int i{}; i < jode.get_dimension(); ++i)
        {
            matrices_data[i * (jode.get_dimension() + 1)] = 1.0;
        }

    };

    /** Set the initial value of the Jacobian matrix with respect to state.
     * <p>
     * If this method is not called, the initial value of the Jacobian
     * matrix with respect to state is set to identity.
     * </p>
     * <p>
     * This method must be called <em>before {@link #set_up_initial_state(ODE_State)}</em>
     * </p>
     * @param dYdY0 initial Jacobian matrix w.r.t. state
     * @exception  if matrix dimensions are incorrect
     */
    void set_initial_main_state_jacobian(const std::vector<std::vector<double>>& dYdY0)     
    {
        // Check dimensions
        check_dimension(my_jode.get_dimension(), dYdY0);
        check_dimension(my_jode.get_dimension(), dYdY0[0]);

        // store the matrix in row major order as a single dimension array
        int i{};
        for (const std::vector<double> row : dYdY0) 
        {
            System.arraycopy(row, 0, my_matrices_data, i, my_jode.get_dimension());
            i += my_jode.get_dimension();
        }

    }

    /** Set the initial value of a column of the Jacobian matrix with respect to one parameter.
     * <p>
     * If this method is not called for some parameter, the initial value of
     * the column of the Jacobian matrix with respect to this parameter is set to zero.
     * </p>
     * <p>
     * This method must be called <em>before {@link #set_up_initial_state(ODE_State)}</em>
     * </p>
     * @param p_name parameter name
     * @param dYdP initial Jacobian column vector with respect to the parameter
     * @exception  if a parameter is not supported
     * @ if the column vector does not match state dimension
     */
    void set_initial_parameter_jacobian(const std::string& p_name, const std::vector<double>& dYdP)
    {
        // Check dimensions
        check_dimension(my_jode.get_dimension(), dYdP);

        // store the column in a global single dimension array
        int i = my_jode.get_dimension() * my_jode.get_dimension();
        for (const std::string known_parameter : my_jode.get_parameters_names())
        {
            if (p_name.equals(known_parameter)) 
            {
                System.arraycopy(dYdP, 0, matrices_data, i, my_jode.get_dimension());
                return;
            }
            i += my_jode.get_dimension();
        }

        throw (Localized_ODE_Formats.UNKNOWN_PARAMETER, p_name);

    }

    /** Set up initial state.
     * <p>
     * This method inserts the initial Jacobian matrices data into
     * an {@link ODE_State ODE state} by overriding the additional
     * state components corresponding to the instance. It must be
     * called prior to integrate the equations.
     * </p>
     * <p>
     * This method must be called <em>after</em>
     * {@link #set_initial_main_state_jacobian(std::vector<std::vector<double>>)} and
     * {@link #set_initial_parameter_jacobian(std::string, std::vector<double>)}.
     * </p>
     * @param initial_state initial state, without the initial Jacobians
     * matrices
     * @return a instance of initial state, with the initial Jacobians
     * matrices properly initialized
     */
    ODE_State set_up_initial_state(const ODE_State initial_state)
    { // NOPMD - PMD false positive

        // insert the matrices data into secondary states
        auto secondary = std::vector<std::vector<double>>(my_expandable.get_mapper().get_number_of_equations() - 1);
        for (int i{}; i < initial_state.get_number_of_secondary_states(); ++i) 
        {
            if (i + 1 != index) 
            {
                secondary[i] = initial_state.get_secondary_state(i + 1);
            }
        }
        secondary[index - 1] = matrices_data;

        // create an updated initial state
        return ODE_State(initial_state.get_time(), initial_state.get_primary_state(), secondary);

    }

    /** Extract the Jacobian matrix with respect to state.
     * @param state state from which to extract Jacobian matrix
     * @return Jacobian matrix dY/dY0 with respect to state.
     */
    std::vector<std::vector<double>> extract_main_set_jacobian(const ODE_State state) 
    {

        // get current state for this set of equations from the expandable fode
        const std::vector<double> p = state.get_secondary_state(index);

        const std::vector<std::vector<double>> dYdY0 = std::vector<double>(my_jode.get_dimension()][my_jode.get_dimension()];
        int j = 0;
        for (int i{}; i < my_jode.get_dimension(); i++)
        {
            System.arraycopy(p, j, dYdY0[i], 0, my_jode.get_dimension());
            j += my_jode.get_dimension();
        }

        return dYdY0;

    }

    /** Extract the Jacobian matrix with respect to one parameter.
     * @param state state from which to extract Jacobian matrix
     * @param p_name name of the parameter for the computed Jacobian matrix
     * @return Jacobian matrix dY/dP with respect to the named parameter
     */
    std::vector<double> extract_parameter_jacobian(const ODE_State state, const std::string p_name) 
    {

        // get current state for this set of equations from the expandable fode
        const std::vector<double> p = state.get_secondary_state(index);

        const std::vector<double> dYdP = std::vector<double>(my_jode.get_dimension()];
        int i = my_jode.get_dimension() * my_jode.get_dimension();
        for (const std::string known_parameter : my_jode.get_parameters_names())
        {
            if (p_name.equals(known_parameter)) 
            {
                System.arraycopy(p, i, dYdP, 0, my_jode.get_dimension());
                break;
            }
            i += my_jode.get_dimension();
        }

        return dYdP;

    }


    /**
     * Special exception for equations mismatch.
     */
    static class Mismatched_Equations extends  
    {

        /** Serializable UID. */
        private static const long serial_version_uid = 20120902L;

        /** Simple constructor. */
        public Mismatched_Equations() 
        {
            super(Localized_ODE_Formats.UNMATCHED_ODE_IN_EXPANDED_SET);
        }

    }

}



