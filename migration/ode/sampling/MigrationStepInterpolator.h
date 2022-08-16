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

  //package org.hipparchus.migration.ode.sampling;

  //import java.io.Serializable;

  //import org.hipparchus.exception.Math_Illegal_State_Exception;
  //import org.hipparchus.ode.ODE_State_And_Derivative;
  //import org.hipparchus.ode.sampling.ODE_StateInterpolator;

  /** This interface represents an interpolator over the last step
   * during an ODE integration.
   *
   * <p>The various ODE integrators provide objects implementing this
   * interface to the step handlers. These objects are often custom
   * objects tightly bound to the integrator internal algorithms. The
   * handlers can use these objects to retrieve the state vector at
   * intermediate times between the previous and the current grid points
   * (this feature is often called dense output).</p>
   * </p>
   *
   * @see org.hipparchus.ode.First_Order_Integrator
   * @see Step_Handler
   * @deprecated as of 1.0, this class is a temporary wrapper between
   * {@link ODE_StateInterpolator} and {@link Migration_Step_Interpolator}
   */
   //@Deprecated
class Migration_Step_Interpolator : public Step_Interpolator
{
	/** Underlying interpolator. */
	private const ODE_StateInterpolator interpolator;

	/** Interpolated state. */
	private ODE_State_And_Derivative interpolated;

	/** Simple constructor.
	 * @param interpolator underlying interpolator
	 */
	Migration_Step_Interpolator(const ODE_StateInterpolator interpolator)
	{
		this.interpolator = interpolator;
		this.interpolated = interpolator.get_current_state();
	}

	/** {@inherit_doc} */
	//override
	@Deprecated
		public double get_previous_time()
	{
		return get_previous_state().get_time();
	}

	/** {@inherit_doc} */
	//override
	@Deprecated
		public double get_current_time() const
	{
		return get_current_state().get_time();
	}

	/** {@inherit_doc} */
	//override
	@Deprecated
		public double get_interpolated_time()
	{
		return interpolated.get_time();
	}

	/** {@inherit_doc} */
	//override
	@Deprecated
		public void set_interpolated_time(const double time)
	{
		interpolated = get_interpolated_state(time);
	}

	/** {@inherit_doc} */
	//override
	@Deprecated
		public std::vector<double> get_interpolated_state() Math_Illegal_State_Exception
	{
		return interpolated.get_primary_state();
	}

	/** {@inherit_doc} */
	//override
	@Deprecated
		public std::vector<double> get_interpolated_derivatives() Math_Illegal_State_Exception
	{
		return interpolated.get_primary_derivative();
	}

	/** {@inherit_doc} */
	//override
	@Deprecated
		public std::vector<double> get_interpolated_secondary_state(const int index) Math_Illegal_State_Exception
	{
		return interpolated.get_secondary_state(index);
	}

	/** {@inherit_doc} */
	//override
	@Deprecated
		public std::vector<double> get_interpolated_secondary_derivatives(const int index) Math_Illegal_State_Exception
	{
		return interpolated.get_secondary_derivative(index);
	}

	/** {@inherit_doc} */
	//override
	public bool is_forward()
	{
		return interpolator.is_forward();
	}

	/** {@inherit_doc} */
	//override
	public Migration_Step_Interpolator copy() Math_Illegal_State_Exception
	{
		return Migration_Step_Interpolator(interpolator);
	}

	/** {@inherit_doc} */
	//override
	public ODE_State_And_Derivative get_previous_state()
	{
		return interpolator.get_previous_state();
	}

	/** {@inherit_doc} */
	//override
	public bool is_previous_state_interpolated()
	{
		return interpolator.is_previous_state_interpolated();
	}

	/** {@inherit_doc} */
	//override
	public ODE_State_And_Derivative get_current_state()
	{
		return interpolator.get_current_state();
	}

	/** {@inherit_doc} */
	//override
	public bool is_current_state_interpolated()
	{
		return interpolator.is_current_state_interpolated();
	}

	/** {@inherit_doc} */
	//override
	public ODE_State_And_Derivative get_interpolated_state(const double time)
	{
		return interpolator.get_interpolated_state(time);
	}

	/**
	 * Replace the instance with a data transfer object for serialization.
	 * @return data transfer object that will be serialized
	 */
	private Object write_replace()
	{
		return Data_Transfer_Object(interpolator, interpolated.get_time());
	}

	/** Internal class used only for serialization. */
	private static class Data_Transfer_Object
	{
		20160328L;

		/** Underlying interpolator.
		 * @serial
		 */
		private const ODE_StateInterpolator interpolator;

		/** Interpolation time.
		 * @serial
		 */
		private const double time;

		/** Simple constructor.
		 * @param interpolator underlying interpolator
		 * @param time interpolation time
		 */
		Data_Transfer_Object(const ODE_StateInterpolator interpolator, const double time)
		{
			this.interpolator = interpolator;
			this.time = time;
		}

		/** Replace the deserialized data transfer object with a {@link Migration_Step_Interpolator}.
		 * @return replacement {@link Migration_Step_Interpolator}
		 */
		private Object read_resolve()
		{
			const Migration_Step_Interpolator msi = Migration_Step_Interpolator(interpolator);
			msi.set_interpolated_time(time);
			return msi;
		}
	}
}
