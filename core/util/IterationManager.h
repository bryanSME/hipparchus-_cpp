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
//package org.hipparchus.util;

//import java.util.Collection;
//import java.util.concurrent.Copy_On_Write_Array_list;
#include <vector>
#include "IterationListener.h"
#include "Incrementor.hpp"


//import org.hipparchus.exception.Math_Illegal_State_Exception;

/**
 * This virtual class provides a general framework for managing iterative
 * algorithms. The maximum number of iterations can be set, and methods are
 * provided to monitor the current iteration count. A lightweight event
 * framework is also provided.
 */
class Iteration_Manager 
{
private:
    /** Keeps a count of the number of iterations. */
    Incrementor my_iterations;

    /** The collection of all listeners attached to this iterative algorithm. */
    std::vector<Iteration_Listener> my_listeners;

public:
    /**
     * Creates a instance of this class.
     *
     * @param max_iterations the maximum number of iterations
     */
    Iteration_Manager(const int max_iterations) : my_iterations{ Incrementor(max_iterations) }, my_listeners{ Copy_On_Write_Array_list<>() } {};

    /**
     * Creates a instance of this class.
     *
     * @param max_iterations the maximum number of iterations
     * @param call_back the function to be called when the maximum number of
     * iterations has been reached
     * @org.hipparchus.exception.Null_Argument_Exception if {@code call_back} is {@code NULL}
     */
    Iteration_Manager(const int max_iterations, const Incrementor::Max_Count_Exceeded_Callback call_back) : my_iterations{ Incrementor(max_iterations, call_back) }, my_listeners{ Copy_On_Write_Array_list<>() } {};

    /**
     * Attaches a listener to this manager.
     *
     * @param listener A {@code Iteration_Listener} object.
     */
    void add_iteration_listener(const Iteration_Listener& listener) 
    {
        my_listeners.push_back(listener);
    }

    /**
     * Informs all registered listeners that the initial phase (prior to the
     * main iteration loop) has been completed.
     *
     * @param e The {@link Iteration_Event} object.
     */
    void fire_initialization_event(const Iteration_Event& e) 
    {
        for (auto& l : my_listeners) 
        {
            l.initialization_performed(e);
        }
    }

    /**
     * Informs all registered listeners that a iteration (in the main
     * iteration loop) has been performed.
     *
     * @param e The {@link Iteration_Event} object.
     */
    void fire_iteration_performed_event(const Iteration_Event& e) 
    {
        for (auto& l : my_listeners)
        {
            l.iteration_performed(e);
        }
    }

    /**
     * Informs all registered listeners that a iteration (in the main
     * iteration loop) has been started.
     *
     * @param e The {@link Iteration_Event} object.
     */
    void fire_iteration_started_event(const Iteration_Event& e) 
    {
        for (auto& l : my_listeners)
        {
            l.iteration_started(e);
        }
    }

    /**
     * Informs all registered listeners that the const phase (post-iterations)
     * has been completed.
     *
     * @param e The {@link Iteration_Event} object.
     */
    void fire_termination_event(const Iteration_Event& e) 
    {
        for (auto& l : my_listeners)
        {
            l.termination_performed(e);
        }
    }

    /**
     * Returns the number of iterations of this solver, 0 if no iterations has
     * been performed yet.
     *
     * @return the number of iterations.
     */
    int get_iterations() const
    {
        return my_iterations.get_count();
    }

    /**
     * Returns the maximum number of iterations.
     *
     * @return the maximum number of iterations.
     */
    int get_max_iterations() const
    {
        return my_iterations.get_maximal_count();
    }

    /**
     * Increments the iteration count by one, and an exception if the
     * maximum number of iterations is reached. This method should be called at
     * the beginning of a iteration.
     *
     * @Math_Illegal_State_Exception if the maximum number of iterations is
     * reached.
     */
    void increment_iteration_count()
    {
        my_iterations.increment();
    }

    /**
     * Removes the specified iteration listener from the list of listeners
     * currently attached to {@code this} object. Attempting to remove a
     * listener which was <em>not</em> previously registered does not cause any
     * error.
     *
     * @param listener The {@link Iteration_Listener} to be removed.
     */
    void remove_iteration_listener(const Iteration_Listener& listener) 
    {
        my_listeners.erase(std::find(my_listeners.begin(), my_listeners.end(), listener));
    }

    /**
     * Sets the iteration count to 0. This method must be called during the
     * initial phase.
     */
    void reset_iteration_count() 
    {
        my_iterations.reset();
    }
};