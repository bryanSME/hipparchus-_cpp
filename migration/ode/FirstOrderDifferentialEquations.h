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
//import org.hipparchus.ode.Ordinary_Differential_Equation;



/** This interface represents a first order differential equations set.
 *
 * <p>This interface should be implemented by all real first order
 * differential equation problems before they can be handled by the
 * integrators {@link org.hipparchus.ode.ODE_Integrator#integrate} method.</p>
 *
 * <p>A first order differential equations problem, as seen by an
 * integrator is the time derivative <code>dY/dt</code> of a state
 * vector <code>Y</code>, both being one dimensional arrays. From the
 * integrator point of view, this derivative depends only on the
 * current time <code>t</code> and on the state vector
 * <code>Y</code>.</p>
 *
 * <p>For real problems, the derivative depends also on parameters
 * that do not belong to the state vector (dynamical model constants
 * for example). These constants are completely outside of the scope
 * of this interface, the classes that implement it are allowed to
 * handle them as they want.</p>
 *
 * @see org.hipparchus.ode.ODE_Integrator
 * @see org.hipparchus.ode.First_Order_Converter
 * @see Second_Order_Differential_Equations
 * @deprecated as of 1.0, replaced with {@link Ordinary_Differential_Equation}
 */
@Deprecated
class First_Order_Differential_Equations extends Ordinary_Differential_Equation 
{

    /** {@inherit_doc}
     * <p>
     * The default implementation calls {@link #compute_derivatives(double, std::vector<double>, std::vector<double>)}.
     * </p>
     */
    //override
    default std::vector<double> compute_derivatives(double t, std::vector<double> y) 
    {
        const std::vector<double> y_dot = std::vector<double>(y.size()];
        compute_derivatives(t, y, y_dot);
        return y_dot;
    }

    /** Get the current time derivative of the state vector.
     * @param t current value of the independent <I>time</I> variable
     * @param y array containing the current value of the state vector
     * @param y_dot placeholder array where to put the time derivative of the state vector
     * @exception Math_Illegal_State_Exception if the number of functions evaluations is exceeded
     * @exception  if arrays dimensions do not match equations settings
     */
    void compute_derivatives(double t, std::vector<double> y, std::vector<double> y_dot)
        , Math_Illegal_State_Exception;

}


