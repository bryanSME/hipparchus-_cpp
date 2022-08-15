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
//package org.hipparchus.optim.nonlinear.vector.leastsquares;

//import org.hipparchus.analysis.Multivariate_Matrix_Function;
//import org.hipparchus.analysis.MultivariateVector_function;
//import org.hipparchus.linear.Array_Real_Vector;
//import org.hipparchus.linear.Real_Matrix;
//import org.hipparchus.linear.Real_Vector;
//import org.hipparchus.optim.Convergence_Checker;
//import org.hipparchus.optim.Point_Vector_Value_Pair;
//import org.hipparchus.optim.nonlinear.vector.leastsquares.Least_Squares_Problem.Evaluation;

/**
 * A mutable builder for {@link Least_Squares_Problem}s.
 *
 * @see Least_Squares_Factory
 */
class Least_Squares_Builder 
{

    /** max evaluations */
    private int max_evaluations;
    /** max iterations */
    private int max_iterations;
    /** convergence checker */
    private Convergence_Checker<Evaluation> checker;
    /** model function */
    private Multivariate_Jacobian_Function model;
    /** observed values */
    private Real_Vector target;
    /** initial guess */
    private Real_Vector start;
    /** weight matrix */
    private Real_Matrix weight;
    /**
     * Lazy evaluation.
     *
     */
    private bool lazy_evaluation;
    /** Validator.
     *
     */
    private Parameter_Validator param_validator;


    /**
     * Construct a {@link Least_Squares_Problem} from the data in this builder.
     *
     * @return a {@link Least_Squares_Problem}.
     */
    public Least_Squares_Problem build() 
    {
        return Least_Squares_Factory.create(model, target, start, weight, checker, max_evaluations, max_iterations, lazy_evaluation, param_validator);
    }

    /**
     * Configure the max evaluations.
     *
     * @param new_max_evaluations the maximum number of evaluations permitted.
     * @return this
     */
    public Least_Squares_Builder max_evaluations(const int& new_max_evaluations) 
    {
        this.max_evaluations = new_max_evaluations;
        return this;
    }

    /**
     * Configure the max iterations.
     *
     * @param new_max_iterations the maximum number of iterations permitted.
     * @return this
     */
    public Least_Squares_Builder max_iterations(const int& new_max_iterations) 
    {
        this.max_iterations = new_max_iterations;
        return this;
    }

    /**
     * Configure the convergence checker.
     *
     * @param new_checker the convergence checker.
     * @return this
     */
    public Least_Squares_Builder checker(const Convergence_Checker<Evaluation> new_checker) 
    {
        this.checker = new_checker;
        return this;
    }

    /**
     * Configure the convergence checker.
     * <p/>
     * This function is an overloaded version of {@link #checker(Convergence_Checker)}.
     *
     * @param new_checker the convergence checker.
     * @return this
     */
    public Least_Squares_Builder checker_pair(const Convergence_Checker<Point_Vector_Value_Pair> new_checker) 
    {
        return this.checker(Least_Squares_Factory.evaluation_checker(new_checker));
    }

    /**
     * Configure the model function.
     *
     * @param value the model function value
     * @param jacobian the Jacobian of {@code value}
     * @return this
     */
    public Least_Squares_Builder model(const MultivariateVector_function value, const Multivariate_Matrix_Function jacobian) 
    {
        return model(Least_Squares_Factory.model(value, jacobian));
    }

    /**
     * Configure the model function.
     *
     * @param new_model the model function value and Jacobian
     * @return this
     */
    public Least_Squares_Builder model(const Multivariate_Jacobian_Function new_model) 
    {
        this.model = new_model;
        return this;
    }

    /**
     * Configure the observed data.
     *
     * @param new_target the observed data.
     * @return this
     */
    public Least_Squares_Builder target(const Real_Vector new_target) 
    {
        this.target = new_target;
        return this;
    }

    /**
     * Configure the observed data.
     *
     * @param new_target the observed data.
     * @return this
     */
    public Least_Squares_Builder target(const std::vector<double> new_target) 
    {
        return target(new Array_Real_Vector(new_target, false));
    }

    /**
     * Configure the initial guess.
     *
     * @param new_start the initial guess.
     * @return this
     */
    public Least_Squares_Builder start(const Real_Vector new_start) 
    {
        this.start = new_start;
        return this;
    }

    /**
     * Configure the initial guess.
     *
     * @param new_start the initial guess.
     * @return this
     */
    public Least_Squares_Builder start(const std::vector<double> new_start) 
    {
        return start(new Array_Real_Vector(new_start, false));
    }

    /**
     * Configure the weight matrix.
     *
     * @param new_weight the weight matrix
     * @return this
     */
    public Least_Squares_Builder weight(const Real_Matrix new_weight) 
    {
        this.weight = new_weight;
        return this;
    }

    /**
     * Configure whether evaluation will be lazy or not.
     *
     * @param new_value Whether to perform lazy evaluation.
     * @return this object.
     *
     */
    public Least_Squares_Builder lazy_evaluation(const bool new_value) 
    {
        lazy_evaluation = new_value;
        return this;
    }

    /**
     * Configure the validator of the model parameters.
     *
     * @param new_validator Parameter validator.
     * @return this object.
     *
     */
    public Least_Squares_Builder parameter_validator(const Parameter_Validator new_validator) 
    {
        param_validator = new_validator;
        return this;
    }
}


