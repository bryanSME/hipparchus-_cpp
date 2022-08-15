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

//import org.hipparchus.ode.sampling.ODE_StateInterpolator;
//import org.hipparchus.ode.sampling.ODE_Step_Handler;

/**
 * This class is a step handler that does nothing.

 * <p>This class is provided as a convenience for users who are only
 * interested in the const state of an integration and not in the
 * intermediate steps. Its handle_step method does nothing.</p>
 *
 * <p>sin_ce this class has no internal state, it is implemented using
 * the Singleton design pattern. This means that only one instance is
 * ever created, which can be retrieved using the get_instance
 * method. This explains why there is no public constructor.</p>
 *
 * @deprecated as of 1.0, this class is not used anymore
 */
@Deprecated
class DummyStep_Handler : ODE_Step_Handler 
{

    /** Private constructor.
     * The constructor is private to prevent users from creating
     * instances (Singleton design-pattern).
     */
    private DummyStep_Handler() 
    {
    }

    /** Get the only instance.
     * @return the only instance
     */
    public static DummyStep_Handler get_instance() 
    {
        return Lazy_Holder.INSTANCE;
    }

    /** {@inherit_doc} */
    //override
    public void handle_step(const ODE_StateInterpolator interpolator) 
    {
    }

    // CHECKSTYLE: stop Hide_Utility_Class_Constructor
    /** Holder for the instance.
     * <p>We use here the Initialization On Demand Holder Idiom.</p>
     */
    private static class Lazy_Holder 
    {
        /** Cached field instance. */
        private static const DummyStep_Handler INSTANCE = DummyStep_Handler();
    }
    // CHECKSTYLE: resume Hide_Utility_Class_Constructor

    /** Handle deserialization of the singleton.
     * @return the singleton instance
     */
    private Object read_resolve() 
    {
        // return the singleton instance
        return Lazy_Holder.INSTANCE;
    }

}


