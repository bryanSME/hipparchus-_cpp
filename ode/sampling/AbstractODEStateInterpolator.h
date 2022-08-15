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

//package org.hipparchus.ode.sampling;

//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.ode.Equations_mapper;
//import org.hipparchus.ode.ODE_State_And_Derivative;
//import org.hipparchus.util.FastMath;

/** This virtual class represents an interpolator over the last step
 * during an ODE integration.
 *
 * <p>The various ODE integrators provide objects extending this class
 * to the step handlers. The handlers can use these objects to
 * retrieve the state vector at intermediate times between the
 * previous and the current grid points (dense output).</p>
 *
 * @see org.hipparchus.ode.ODE_Integrator
 * @see ODE_Step_Handler
 */

class AbstractODE_StateInterpolator
    : ODE_StateInterpolator 
    {

    /** Serializable UID. */
    private static const long serial_version_uid = 20160328L;

    /** Global previous state. */
    private const ODE_State_And_Derivative global_previous_state;

    /** Global current state. */
    private const ODE_State_And_Derivative global_current_state;

    /** Soft previous state. */
    private const ODE_State_And_Derivative soft_previous_state;

    /** Soft current state. */
    private const ODE_State_And_Derivative soft_current_state;

    /** integration direction. */
    private const bool forward;

    /** Mapper for ODE equations primary and secondary components. */
    private Equations_mapper mapper;

    /** Simple constructor.
     * @param is_forward integration direction indicator
     * @param global_previous_state start of the global step
     * @param global_current_state end of the global step
     * @param soft_previous_state start of the restricted step
     * @param soft_current_state end of the restricted step
     * @param equations_mapper mapper for ODE equations primary and secondary components
     */
    protected AbstractODE_StateInterpolator(const bool is_forward, const ODE_State_And_Derivative global_previous_state, const ODE_State_And_Derivative global_current_state, const ODE_State_And_Derivative soft_previous_state, const ODE_State_And_Derivative soft_current_state, const Equations_mapper equations_mapper) 
    {
        this.forward             = is_forward;
        this.global_previous_state = global_previous_state;
        this.global_current_state  = global_current_state;
        this.soft_previous_state   = soft_previous_state;
        this.soft_current_state    = soft_current_state;
        this.mapper              = equations_mapper;
    }

    /** Create a restricted version of the instance.
     * <p>
     * The instance is not changed at all.
     * </p>
     * @param previous_state start of the restricted step
     * @param current_state end of the restricted step
     * @return restricted version of the instance
     * @see #get_previous_state()
     * @see #get_current_state()
     */
    public AbstractODE_StateInterpolator restrict_step(const ODE_State_And_Derivative previous_state, const ODE_State_And_Derivative current_state) 
    {
        return create(forward, global_previous_state, global_current_state, previous_state, current_state, mapper);
    }

    /** Create a instance.
     * @param new_forward integration direction indicator
     * @param new_global_previous_state start of the global step
     * @param new_global_current_state end of the global step
     * @param new_soft_previous_state start of the restricted step
     * @param new_soft_current_state end of the restricted step
     * @param new_mapper equations mapper for the all equations
     * @return a instance
     */
    protected virtual AbstractODE_StateInterpolator create(bool new_forward, ODE_State_And_Derivative new_global_previous_state, ODE_State_And_Derivative new_global_current_state, ODE_State_And_Derivative new_soft_previous_state, ODE_State_And_Derivative new_soft_current_state, Equations_mapper new_mapper);

    /**
     * Get the previous global grid point state.
     * @return previous global grid point state
     */
    public ODE_State_And_Derivative get_global_previous_state() 
    {
        return global_previous_state;
    }

    /**
     * Get the current global grid point state.
     * @return current global grid point state
     */
    public ODE_State_And_Derivative get_global_current_state() 
    {
        return global_current_state;
    }

    /** {@inherit_doc} */
    //override
    public ODE_State_And_Derivative get_previous_state() 
    {
        return soft_previous_state;
    }

    /** {@inherit_doc} */
    //override
    public bool is_previous_state_interpolated() 
    {
        return soft_previous_state != global_previous_state;
    }

    /** {@inherit_doc} */
    //override
    public ODE_State_And_Derivative get_current_state() 
    {
        return soft_current_state;
    }

    /** {@inherit_doc} */
    //override
    public bool is_current_state_interpolated() 
    {
        return soft_current_state != global_current_state;
    }

    /** {@inherit_doc} */
    //override
    public ODE_State_And_Derivative get_interpolated_state(const double time) 
    {
        if (std::abs(global_current_state.get_time() - global_previous_state.get_time()) <=
                FastMath.ulp(global_current_state.get_time())) 
                {
            return global_current_state;
        }
        const double theta_h         = time - global_previous_state.get_time();
        const double one_minus_theta_h = global_current_state.get_time() - time;
        const double theta          = theta_h / (global_current_state.get_time() - global_previous_state.get_time());
        return compute_interpolated_state_and_derivatives(mapper, time, theta, theta_h, one_minus_theta_h);
    }

    /** {@inherit_doc} */
    //override
    public bool is_forward() 
    {
        return forward;
    }

    /** Get the mapper for ODE equations primary and secondary components.
     * @return mapper for ODE equations primary and secondary components
     */
    protected Equations_mapper get_mapper() 
    {
        return mapper;
    }

    /** Compute the state and derivatives at the interpolated time.
     * This is the main processing method that should be implemented by
     * the derived classes to perform the interpolation.
     * @param equations_mapper mapper for ODE equations primary and secondary components
     * @param time interpolation time
     * @param theta normalized interpolation abscissa within the step
     * (theta is zero at the previous time step and one at the current time step)
     * @param theta_h time gap between the previous time and the interpolated time
     * @param one_minus_theta_h time gap between the interpolated time and
     * the current time
     * @return interpolated state and derivatives
     * @exception Math_Illegal_State_Exception if the number of functions evaluations is exceeded
     */
    protected virtual ODE_State_And_Derivative compute_interpolated_state_and_derivatives(Equations_mapper equations_mapper, double time, double theta, double theta_h, double one_minus_theta_h)
        Math_Illegal_State_Exception;

}


