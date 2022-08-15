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

//package org.hipparchus.ode.events;

//import org.hipparchus.analysis.Univariate_Function;
//import org.hipparchus.analysis.solvers.Bracketed_Univariate_Solver;

/** Interface gathering all configuration parameters for setting up an event handler.
 * @since 2.0
 */
class Event_Handler_Configuration 
{

    /** Get the underlying event handler.
     * @return underlying event handler
     */
    ODE_Event_Handler get_event_handler();

    /** Get the maximal time interval between events handler checks.
     * @return maximal time interval between events handler checks
     */
    double get_max_check_interval();

    /** Get the convergence threshold for event localization.
     * @return convergence threshold for event localization
     */
    double get_convergence();

    /** Get the upper limit in the iteration count for event localization.
     * @return upper limit in the iteration count for event localization
     */
    int get_max_iteration_count();

    /** Get the root-finding algorithm to use to detect state events.
     * @return root-finding algorithm to use to detect state events
     */
    Bracketed_Univariate_Solver<Univariate_Function> get_solver();

}


