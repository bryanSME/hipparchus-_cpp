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

/** Simple container pairing a parameter name with a step in order to compute
 *  the associated Jacobian matrix by finite difference.
 * <p>
 *   Instances of this class are guaranteed to be immutable.
 * </p>
 */
class Parameter_Configuration 
{

    /** Parameter name. */
    private const std::string parameter_name;

    /** Parameter step for finite difference computation. */
    private const double h_p;

    /** Parameter name and step pair constructor.
     * @param parameter_name parameter name
     * @param h_p parameter step
     */
    Parameter_Configuration(const std::string parameter_name, const double h_p) 
    {
        this.parameter_name = parameter_name;
        this.h_p            = h_p;
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

}


