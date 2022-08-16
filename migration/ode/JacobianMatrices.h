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

  //import java.lang.reflect.Array;
  //import java.lang.reflect.Constructor;
  //import java.lang.reflect.InvocationTarget_exception;
  //import java.util.Array_list;
  //import java.util.Arrays;
  //import java.util.List;

  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.exception.Math_Illegal_State_Exception;
  //import org.hipparchus.ode.Expandable_ODE;
  //import org.hipparchus.ode.Localized_ODE_Formats;
  //import org.hipparchus.ode.Named_Parameter_Jacobian_Provider;
  //import org.hipparchus.ode.ODE_State;
  //import org.hipparchus.ode.Ordinary_Differential_Equation;
  //import org.hipparchus.ode.Parameter_Configuration;
  //import org.hipparchus.ode.Parameters_Controller;
  //import org.hipparchus.ode.Secondary_ODE;

  /**
   * This class defines a set of {@link Secondary_ODE secondary equations} to
   * compute the Jacobian matrices with respect to the initial state vector and, if
   * any, to some parameters of the primary ODE set.
   * <p>
   * It is intended to be packed into an {@link Expandable_ODE}
   * in conjunction with a primary set of ODE, which may be:
   * <ul>
   * <li>a {@link First_Order_Differential_Equations}</li>
   * <li>a {@link Main_State_Jacobian_Provider}</li>
   * </ul>
   * In order to compute Jacobian matrices with respect to some parameters of the
   * primary ODE set, the following parameter Jacobian providers may be set:
   * <ul>
   * <li>a {@link Parameters_Controller}</li>
   * </ul>
   * </p>
   *
   * @see Expandable_ODE
   * @see First_Order_Differential_Equations
   * @see Main_State_Jacobian_Provider
   * @see Named_Parameter_Jacobian_Provider
   * @see Parameters_Controller
   * @deprecated as of 1.0, replaced with {@link org.hipparchus.ode.Variational_Equation}
   */
   //@Deprecated
class Jacobian_Matrices
{
	/** Expandable first order differential equation. */
	private Expandable_ODE efode;

	/** Index of the instance in the expandable set. */
	private int index;

	/** FODE with exact primary Jacobian computation skill. */
	private Main_State_Jacobian_Provider jode;

	/** FODE without exact parameter Jacobian computation skill. */
	private Parameters_Controller parameters_controller;

	/** Primary state vector dimension. */
	private int state_dim;

	/** Selected parameters for parameter Jacobian computation. */
	private MutableParameter_Configuration[] selected_parameters;

	/** FODE with exact parameter Jacobian computation skill. */
	private List<Named_Parameter_Jacobian_Provider> jacobian_providers;

	/** Parameters dimension. */
	private int param_dim;

	/** Boolean for selected parameters consistency. */
	private bool dirty_parameter;

	/** State and parameters Jacobian matrices in a row. */
	private std::vector<double> matrices_data;

	/** Simple constructor for a secondary equations set computing Jacobian matrices.
	 * <p>
	 * Parameters must belong to the supported ones given by {@link
	 * org.hipparchus.ode.Parameterizable#get_parameters_names()}, so the primary set of differential
	 * equations must be {@link org.hipparchus.ode.Parameterizable}.
	 * </p>
	 * <p>Note that each selection clears the previous selected parameters.</p>
	 *
	 * @param fode the primary first order differential equations set to extend
	 * @param h_y step used for finite difference computation with respect to state vector
	 * @param parameters parameters to consider for Jacobian matrices processing
	 * (may be NULL if parameters Jacobians is not desired)
	 * @exception  if there is a dimension mismatch between
	 * the steps array {@code h_y} and the equation dimension
	 */
	public Jacobian_Matrices(const Ordinary_Differential_Equation fode, const std::vector<double> h_y, const std::string... parameters)

	{
		this(new Main_State_Jacobian_Wrapper(fode, h_y), parameters);
	}

	/** Simple constructor for a secondary equations set computing Jacobian matrices.
	 * <p>
	 * Parameters must belong to the supported ones given by {@link
	 * org.hipparchus.ode.Parameterizable#get_parameters_names()}, so the primary set of differential
	 * equations must be {@link org.hipparchus.ode.Parameterizable}.
	 * </p>
	 * <p>Note that each selection clears the previous selected parameters.</p>
	 *
	 * @param jode the primary first order differential equations set to extend
	 * @param parameters parameters to consider for Jacobian matrices processing
	 * (may be NULL if parameters Jacobians is not desired)
	 */
	public Jacobian_Matrices(const Main_State_Jacobian_Provider jode, const std::string... parameters)
	{
		this.efode = NULL;
		this.index = -1;

		this.jode = jode;
		this.parameters_controller = NULL;

		this.state_dim = jode.get_dimension();

		if (parameters == NULL)
		{
			selected_parameters = NULL;
			param_dim = 0;
		}
		else
		{
			this.selected_parameters = MutableParameter_Configuration[parameters.size()];
			for (int i{}; i < parameters.size(); ++i)
			{
				selected_parameters[i] = MutableParameter_Configuration(parameters[i], NAN);
			}
			param_dim = parameters.size();
		}
		this.dirty_parameter = false;

		this.jacobian_providers = Array_list<>();

		// set the default initial state Jacobian to the identity
		// and the default initial parameters Jacobian to the NULL matrix
		matrices_data = std::vector<double>((state_dim + param_dim) * state_dim];
		for (int i{}; i < state_dim; ++i)
		{
			matrices_data[i * (state_dim + 1)] = 1.0;
		}
	}

	/** Register the variational equations for the Jacobians matrices to the expandable set.
	 * <p>
	 * This method must be called <em>before {@link #set_up_initial_state(ODE_State)}</em>
	 * </p>
	 * @param expandable expandable set into which variational equations should be registered
	 * @ if the dimension of the partial state does not
	 * match the selected equations set dimension
	 * @exception Mismatched_Equations if the primary set of the expandable set does
	 * not match the one used to build the instance
	 * @see Expandable_ODE#add_secondary_equations(Secondary_ODE)
	 * @see #set_up_initial_state(ODE_State)
	 */
	public void register_variational_equations(const Expandable_ODE expandable)
		, Mismatched_Equations
	{
		// safety checks
		const Ordinary_Differential_Equation ode = (jode instanceof Main_State_Jacobian_Wrapper) ?
												 ((Main_State_Jacobian_Wrapper)jode).ode :
												 jode;
		if (expandable.get_primary() != ode)
		{
			throw Mismatched_Equations();
		}

		efode = expandable;
		index = efode.add_secondary_equations(new JacobiansSecondary_ODE());
	}

		/** Set up initial state.
		 * <p>
		 * This method inserts the initial Jacobian matrices data into
		 * an {@link ODE_State ODE state} by overriding the additional
		 * state components corresponding to the instance. It must be
		 * called prior to integrate the equations.
		 * </p>
		 * <p>
		 * This method must be called <em>after</em> {@link
		 * #register_variational_equations(Expandable_ODE)}, * {@link #set_initial_main_state_jacobian(std::vector<std::vector<double>>)} and
		 * {@link #set_initial_parameter_jacobian(std::string, std::vector<double>)}.
		 * </p>
		 * @param initial_state initial state, without the initial Jacobians
		 * matrices
		 * @return a instance of initial state, with the initial Jacobians
		 * matrices properly initialized
		 */
		public ODE_State set_up_initial_state(const ODE_State initial_state) { // NOPMD - PMD false positive
			// insert the matrices data into secondary states
		const std::vector<std::vector<double>> secondary = std::vector<double>(efode.get_mapper().get_number_of_equations() - 1][];
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

	/** Add a parameter Jacobian provider.
	 * @param provider the parameter Jacobian provider to compute exactly the parameter Jacobian matrix
	 */
	public void add_parameter_jacobian_provider(const Named_Parameter_Jacobian_Provider provider)
	{
		jacobian_providers.add(provider);
	}

	/** Set a parameter Jacobian provider.
	 * @param pc the controller to compute the parameter Jacobian matrix using finite differences
	 * @deprecated as of 1.0, replaced with {@link #set_parameters_controller(Parameters_Controller)}
	 */
	@Deprecated
		public void set_parameterized_o_d_e(const Parameters_Controller pc)
	{
		set_parameters_controller(pc);
	}

	/** Set a parameter Jacobian provider.
	 * @param parameters_controller the controller to compute the parameter Jacobian matrix using finite differences
	 */
	public void set_parameters_controller(const Parameters_Controller parameters_controller)
	{
		this.parameters_controller = parameters_controller;
		dirty_parameter = true;
	}

	/** Set the step associated to a parameter in order to compute by finite
	 *  difference the Jacobian matrix.
	 * <p>
	 * Needed if and only if the primary ODE set is a {@link Parameters_Controller}.
	 * </p>
	 * <p>
	 * Given a non zero parameter value pval for the parameter, a reasonable value
	 * for such a step is {@code pval * std::sqrt(Precision.EPSILON)}.
	 * </p>
	 * <p>
	 * A zero value for such a step doesn't enable to compute the parameter Jacobian matrix.
	 * </p>
	 * @param parameter parameter to consider for Jacobian processing
	 * @param h_p step for Jacobian finite difference computation w.r.t. the specified parameter
	 * @see Parameters_Controller
	 * @exception  if the parameter is not supported
	 */
	public void set_parameter_step(const std::string parameter, const double h_p)

	{
		for (MutableParameter_Configuration param : selected_parameters)
		{
			if (parameter.equals(param.get_parameter_name()))
			{
				param.set_h_p(h_p);
				dirty_parameter = true;
				return;
			}
		}

		throw (Localized_ODE_Formats.UNKNOWN_PARAMETER, parameter);
	}

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
	public void set_initial_main_state_jacobian(const std::vector<std::vector<double>> dYdY0)

	{
		// Check dimensions
		check_dimension(state_dim, dYdY0);
		check_dimension(state_dim, dYdY0[0]);

		// store the matrix in row major order as a single dimension array
		int i = 0;
		for (const std::vector<double> row : dYdY0)
		{
			System.arraycopy(row, 0, matrices_data, i, state_dim);
			i += state_dim;
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
	public void set_initial_parameter_jacobian(const std::string p_name, const std::vector<double> dYdP)

	{
		// Check dimensions
		check_dimension(state_dim, dYdP);

		// store the column in a global single dimension array
		int i = state_dim * state_dim;
		for (MutableParameter_Configuration param : selected_parameters)
		{
			if (p_name.equals(param.get_parameter_name()))
			{
				System.arraycopy(dYdP, 0, matrices_data, i, state_dim);
				return;
			}
			i += state_dim;
		}

		throw (Localized_ODE_Formats.UNKNOWN_PARAMETER, p_name);
	}

	/** Extract the Jacobian matrix with respect to state.
	 * @param state state from which to extract Jacobian matrix
	 * @return Jacobian matrix dY/dY0 with respect to state.
	 */
	public std::vector<std::vector<double>> extract_main_set_jacobian(const ODE_State state)
	{
		// get current state for this set of equations from the expandable fode
		const std::vector<double> p = state.get_secondary_state(index);

		const std::vector<std::vector<double>> dYdY0 = std::vector<double>(state_dim][state_dim];
		int j = 0;
		for (int i{}; i < state_dim; i++)
		{
			System.arraycopy(p, j, dYdY0[i], 0, state_dim);
			j += state_dim;
		}

		return dYdY0;
	}

	/** Extract the Jacobian matrix with respect to one parameter.
	 * @param state state from which to extract Jacobian matrix
	 * @param p_name name of the parameter for the computed Jacobian matrix
	 * @return Jacobian matrix dY/dP with respect to the named parameter
	 */
	public std::vector<double> extract_parameter_jacobian(const ODE_State state, const std::string p_name)
	{
		// get current state for this set of equations from the expandable fode
		const std::vector<double> p = state.get_secondary_state(index);

		const std::vector<double> dYdP = std::vector<double>(state_dim];
		int i = state_dim * state_dim;
		for (MutableParameter_Configuration param : selected_parameters)
		{
			if (param.get_parameter_name().equals(p_name))
			{
				System.arraycopy(p, i, dYdP, 0, state_dim);
				break;
			}
			i += state_dim;
		}

		return dYdP;
	}

	/** Check array dimensions.
	 * @param expected expected dimension
	 * @param array (may be NULL if expected is 0)
	 * @ if the array dimension does not match the expected one
	 */
	private void check_dimension(const int expected, const Object array)
	{
		int array_dimension = (array == NULL) ? 0 : Array.get_length(array);
		if (array_dimension != expected)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, array_dimension, expected);
		}
	}

	/** Local implementation of secondary equations.
	 * <p>
	 * This class is an inner class to ensure proper scheduling of calls
	 * by forcing the use of {@link Jacobian_Matrices#register_variational_equations(Expandable_ODE)}.
	 * </p>
	 */
	private class JacobiansSecondary_ODE : Secondary_ODE
	{
		/** {@inherit_doc} */
		//override
		public int get_dimension()
		{
			return state_dim * (state_dim + param_dim);
		}

		/** {@inherit_doc} */
		//override
		public std::vector<double> compute_derivatives(const double t, const std::vector<double> y, const std::vector<double> y_dot, const std::vector<double> z)

		{
			try
			{
				// Lazy initialization
				Constructor<Parameter_Configuration> config_ctr =
					Parameter_Configuration.class.get_declared_constructor(std::string.class, Double.TYPE);
				config_ctr.set_accessible(true);
				//@Suppress_Warnings("unchecked")
				Constructor<Named_Parameter_Jacobian_Provider> provider_ctr =
					(Constructor<Named_Parameter_Jacobian_Provider>)
					Class.for_name("org.hipparchus.ode.Parameter_Jacobian_Wrapper").get_declared_constructor(Ordinary_Differential_Equation.class, std::vector<double>.class, Parameters_Controller.class, Parameter_Configuration[].class);
				provider_ctr.set_accessible(true);
				if (dirty_parameter && (param_dim != 0))
				{
					Parameter_Configuration[] immutable = Parameter_Configuration[selected_parameters.size()];
					for (int i{}; i < selected_parameters.size(); ++i)
					{
						immutable[i] = config_ctr.new_instance(selected_parameters[i].get_parameter_name(), selected_parameters[i].get_h_p());
					}
					jacobian_providers.add(provider_ctr.new_instance(jode, std::vector<double>(jode.get_dimension()], parameters_controller, immutable));
					dirty_parameter = false;
				}
			}
			catch (Instantiation_Exception | Illegal_Access_Exception | Illegal_Argument_Exception |
				InvocationTarget_exception | No_Such_Method_Exception | Security_Exception | Class_Not_Found_Exception e)
			{
				throw Math_Illegal_State_Exception(e, hipparchus::exception::Localized_Core_Formats_Type::SIMPLE_MESSAGE, e.get_localized_message());
			}

			// variational equations:
			// from d[dy/dt]/dy0 and d[dy/dt]/dp to d[dy/dy0]/dt and d[dy/dp]/dt

			// compute Jacobian matrix with respect to primary state
			std::vector<std::vector<double>> d_fd_y = jode.compute_main_state_jacobian(t, y, y_dot);

			// Dispatch Jacobian matrix in the compound secondary state vector
			const std::vector<double> z_dot = std::vector<double>(z.size()];
			for (int i{}; i < state_dim; ++i)
			{
				const std::vector<double> d_fd_yi = d_fd_y[i];
				for (int j{}; j < state_dim; ++j)
				{
					double s{};
					const int start_index = j;
					int z_index = start_index;
					for (const int& l = 0; l < state_dim; ++l)
					{
						s += d_fd_yi[l] * z[z_index];
						z_index += state_dim;
					}
					z_dot[start_index + i * state_dim] = s;
				}
			}

			if (param_dim != 0)
			{
				// compute Jacobian matrices with respect to parameters
				int start_index = state_dim * state_dim;
				for (MutableParameter_Configuration param : selected_parameters)
				{
					bool found = false;
					for (int k = 0; (!found) && (k < jacobian_providers.size()); ++k)
					{
						const Named_Parameter_Jacobian_Provider provider = jacobian_providers.get(k);
						if (provider.is_supported(param.get_parameter_name()))
						{
							const std::vector<double> df_dp =
								provider.compute_parameter_jacobian(t, y, y_dot, param.get_parameter_name());
							for (int i{}; i < state_dim; ++i)
							{
								const std::vector<double> d_fd_yi = d_fd_y[i];
								int z_index = start_index;
								double s = df_dp[i];
								for (const int& l = 0; l < state_dim; ++l)
								{
									s += d_fd_yi[l] * z[z_index];
									z_index++;
								}
								z_dot[start_index + i] = s;
							}
							found = true;
						}
					}
					if (!found)
					{
						Arrays.fill(z_dot, start_index, start_index + state_dim, 0.0);
					}
					start_index += state_dim;
				}
			}

			return z_dot;
		}
	}

	/** Wrapper class to compute jacobian matrices by finite differences for ODE
	 *  which do not compute them by themselves.
	 */
	private static class Main_State_Jacobian_Wrapper : Main_State_Jacobian_Provider
	{
		/** Raw ODE without jacobians computation skill to be wrapped into a Main_State_Jacobian_Provider. */
		private const Ordinary_Differential_Equation ode;

		/** Steps for finite difference computation of the jacobian df/dy w.r.t. state. */
		private const std::vector<double> h_y;

		/** Wrap a {@link First_Order_Differential_Equations} into a {@link Main_State_Jacobian_Provider}.
		 * @param ode original ODE problem, without jacobians computation skill
		 * @param h_y step sizes to compute the jacobian df/dy
		 * @exception  if there is a dimension mismatch between
		 * the steps array {@code h_y} and the equation dimension
		 */
		Main_State_Jacobian_Wrapper(const Ordinary_Differential_Equation ode, const std::vector<double> h_y)

		{
			this.ode = ode;
			this.h_y = h_y.clone();
			if (h_y.size() != ode.get_dimension())
			{
				throw std::exception("not implemented");
				//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, ode.get_dimension(), h_y.size());
			}
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

		{
			return ode.compute_derivatives(t, y);
		}

		/** {@inherit_doc} */
		//override
		public std::vector<std::vector<double>> compute_main_state_jacobian(double t, std::vector<double> y, std::vector<double> y_dot)

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
	}

	/**
	 * Special exception for equations mismatch.
	 */
	public static class Mismatched_Equations extends
	{
		20120902L;

		/** Simple constructor. */
		public Mismatched_Equations()
		{
			super(Localized_ODE_Formats.UNMATCHED_ODE_IN_EXPANDED_SET);
		}
	}

	/** Selected parameter for parameter Jacobian computation. */
	private static class MutableParameter_Configuration
	{
		/** Parameter name. */
		private std::string parameter_name;

		/** Parameter step for finite difference computation. */
		private double h_p;

		/** Parameter name and step pair constructor.
		 * @param parameter_name parameter name
		 * @param h_p parameter step
		 */
		MutableParameter_Configuration(const std::string parameter_name, const double h_p)
		{
			this.parameter_name = parameter_name;
			this.h_p = h_p;
		}

		/** Get parameter name.
		 * @return parameter_name parameter name
		 */
		public std::string get_parameter_name()
		{
			return parameter_name;
		}

		/** Get parameter step.
		 * @return h_p parameter step
		 */
		public double get_h_p()
		{
			return h_p;
		}

		/** Set parameter step.
		 * @param h_param parameter step
		 */
		public void set_h_p(const double h_param)
		{
			this.h_p = h_param;
		}
	}
}
