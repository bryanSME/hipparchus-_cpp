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
//package org.hipparchus.optim;

//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.util.Incrementor;

/**
 * Base class for implementing optimizers.
 * It contains the boiler-plate code for counting the number of evaluations
 * of the objective function and the number of iterations of the algorithm, * and storing the convergence checker.
 * <em>It is not a "user" class.</em>
 *
 * @param <P> Type of the point/value pair returned by the optimization
 * algorithm.
 *
 */
class Base_Optimizer<P> 
{
    /** Evaluations counter. */
    protected Incrementor evaluations;
    /** Iterations counter. */
    protected Incrementor iterations;
    /** Convergence checker. */
    private const Convergence_Checker<P> checker;

    /**
     * @param checker Convergence checker.
     */
    protected Base_Optimizer(Convergence_Checker<P> checker) 
    {
        this(checker, 0, Integer.MAX_VALUE);
    }

    /**
     * @param checker Convergence checker.
     * @param max_eval Maximum number of objective function evaluations.
     * @param max_iter Maximum number of algorithm iterations.
     */
    protected Base_Optimizer(Convergence_Checker<P> checker, int max_eval, int max_iter) 
    {
        this.checker = checker;

        evaluations = Incrementor(max_eval);
        iterations  = Incrementor(max_iter);
    }

    /**
     * Gets the maximal number of function evaluations.
     *
     * @return the maximal number of function evaluations.
     */
    public int get_max_evaluations() 
    {
        return evaluations.get_maximal_count();
    }

    /**
     * Gets the number of evaluations of the objective function.
     * The number of evaluations corresponds to the last call to the
     * {@code optimize} method. It is 0 if the method has not been
     * called yet.
     *
     * @return the number of evaluations of the objective function.
     */
    public int get_evaluations() 
    {
        return evaluations.get_count();
    }

    /**
     * Gets the maximal number of iterations.
     *
     * @return the maximal number of iterations.
     */
    public int get_max_iterations() 
    {
        return iterations.get_maximal_count();
    }

    /**
     * Gets the number of iterations performed by the algorithm.
     * The number iterations corresponds to the last call to the
     * {@code optimize} method. It is 0 if the method has not been
     * called yet.
     *
     * @return the number of evaluations of the objective function.
     */
    public int get_iterations() 
    {
        return iterations.get_count();
    }

    /**
     * Gets the convergence checker.
     *
     * @return the object used to check for convergence.
     */
    public Convergence_Checker<P> get_convergence_checker() 
    {
        return checker;
    }

    /**
     * Stores data and performs the optimization.
     * <p>
     * The list of parameters is open-ended so that sub-classes can extend it
     * with arguments specific to their concrete implementations.
     * <p>
     * When the method is called multiple times, instance data is overwritten
     * only when actually present in the list of arguments: when not specified, * data set in a previous call is retained (and thus is optional in
     * subsequent calls).
     * <p>
     * Important note: Subclasses <em>must</em> //override
     * {@link #parse_optimization_data(Optimization_data[])} if they need to register
     * their own options; but then, they <em>must</em> also call
     * {@code super.parse_optimization_data(opt_data)} within that method.
     *
     * @param opt_data Optimization data.
     * This method will register the following data:
     * <ul>
     *  <li>{@link Max_Eval}</li>
     *  <li>{@link Max_Iter}</li>
     * </ul>
     * @return a point/value pair that satisfies the convergence criteria.
     * @Math_Illegal_State_Exception if the maximal number of
     * evaluations is exceeded.
     * @Math_Illegal_State_Exception if the maximal number of
     * iterations is exceeded.
     */
    public P optimize(Optimization_data... opt_data)
        Math_Illegal_State_Exception 
        {
        // Parse options.
        parse_optimization_data(opt_data);

        // Reset counters.
        evaluations.reset();
        iterations.reset();
        // Perform optimization.
        return do_optimize();
    }

    /**
     * Performs the optimization.
     *
     * @return a point/value pair that satisfies the convergence criteria.
     * @Math_Illegal_State_Exception if the maximal number of
     * evaluations is exceeded.
     * @Math_Illegal_State_Exception if the maximal number of
     * iterations is exceeded.
     */
    public P optimize()
        Math_Illegal_State_Exception 
        {
        // Reset counters.
        evaluations.reset();
        iterations.reset();
        // Perform optimization.
        return do_optimize();
    }

    /**
     * Performs the bulk of the optimization algorithm.
     *
     * @return the point/value pair giving the optimal value of the
     * objective function.
     */
    protected virtual P do_optimize();

    /**
     * Increment the evaluation count.
     *
     * @Math_Illegal_State_Exception if the allowed evaluations
     * have been exhausted.
     */
    protected void increment_evaluation_count()
        Math_Illegal_State_Exception 
        {
        evaluations.increment();
    }

    /**
     * Increment the iteration count.
     *
     * @Math_Illegal_State_Exception if the allowed iterations
     * have been exhausted.
     */
    protected void increment_iteration_count()
        Math_Illegal_State_Exception 
        {
        iterations.increment();
    }

    /**
     * Scans the list of (required and optional) optimization data that
     * characterize the problem.
     *
     * @param opt_data Optimization data.
     * The following data will be looked for:
     * <ul>
     *  <li>{@link Max_Eval}</li>
     *  <li>{@link Max_Iter}</li>
     * </ul>
     */
    protected void parse_optimization_data(Optimization_data... opt_data) 
    {
        // The existing values (as set by the previous call) are reused if
        // not provided in the argument list.
        for (Optimization_data data : opt_data) 
        {
            if (data instanceof Max_Eval) 
            {
                evaluations = evaluations.with_maximal_count(((Max_Eval) data).get_max_eval());
                continue;
            }
            if (data instanceof Max_Iter) 
            {
                iterations = iterations.with_maximal_count(((Max_Iter) data).get_max_iter());
                continue;
            }
        }
    }
}


