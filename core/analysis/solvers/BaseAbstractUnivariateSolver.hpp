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

//package org.hipparchus.analysis.solvers;

//import org.hipparchus.analysis.Univariate_Function;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.exception.Null_Argument_Exception;
//import org.hipparchus.util.Incrementor;
//import org.hipparchus.util.Math_Utils;
#include <type_traits>
#include "UnivariateSolverUtils.hpp"
#include "../UnivariateFunction.h"
#include "../../util/MathUtils.h"
#include "../../util/Incrementor.h"

/**
 * Provide a default implementation for several functions useful to generic
 * solvers.
 * The default values for relative and function tolerances are 1e-14
 * and 1e-15, respectively. It is however highly recommended to not
 * rely on the default, but rather carefully consider values that match
 * user's expectations, as well as the specifics of each implementation.
 *
 * @param <F> Type of function to solve.
 *
 */
template<typename F, typename std::enable_if<std::is_base_of<Univariate_Function, F>::value>::type* = nullptr>
class Base_Abstract_Univariate_Solver : public Base_Univariate_Solver<F> 
{
private:
    /** Default relative accuracy. */
    static constexpr double DEFAULT_RELATIVE_ACCURACY{ 1e-14 };
    /** Default function value accuracy. */
    static constexpr double DEFAULT_FUNCTION_VALUE_ACCURACY{ 1e-15 };
    /** Function value accuracy. */
    const double my_function_value_accuracy;
    /** Absolute accuracy. */
    const double my_absolute_accuracy;
    /** Relative accuracy. */
    const double my_relative_accuracy;
    /** Evaluations counter. */
    Incrementor my_evaluations = Incrementor();
    /** Lower end of search interval. */
    double my_search_min;
    /** Higher end of search interval. */
    double my_search_max;
    /** Initial guess. */
    double my_search_start;
    /** Function to solve. */
    F my_function;

protected:
    /**
     * Construct a solver with given absolute accuracy.
     *
     * @param absolute_accuracy Maximum absolute error.
     */
    Base_Abstract_Univariate_Solver(const double& absolute_accuracy) 
    {
        this(DEFAULT_RELATIVE_ACCURACY, absolute_accuracy, DEFAULT_FUNCTION_VALUE_ACCURACY);
    }

    /**
     * Construct a solver with given accuracies.
     *
     * @param relative_accuracy Maximum relative error.
     * @param absolute_accuracy Maximum absolute error.
     */
    Base_Abstract_Univariate_Solver(const double& relative_accuracy, const double& absolute_accuracy) 
    {
        Base_Abstract_Univariate_Solver(relative_accuracy, absolute_accuracy, DEFAULT_FUNCTION_VALUE_ACCURACY);
    }

    /**
     * Construct a solver with given accuracies.
     *
     * @param relative_accuracy Maximum relative error.
     * @param absolute_accuracy Maximum absolute error.
     * @param function_value_accuracy Maximum function value error.
     */
    Base_Abstract_Univariate_Solver(const double& relative_accuracy, const double& absolute_accuracy, const double& function_value_accuracy) 
        :
        my_absolute_accuracy{ absolute_accuracy },
        my_relative_accuracy{ relative_accuracy },
        my_function_value_accuracy{ function_value_accuracy }
    {
    }

    /**
     * Compute the objective function value.
     *
     * @param point Point at which the objective function must be evaluated.
     * @return the objective function value at specified point.
     * @Math_Illegal_State_Exception if the maximal number of evaluations
     * is exceeded.
     */
    double compute_objective_value(const double& point)
    {
        increment_evaluation_count();
        return my_function.value(point);
    }

    /**
     * Prepare for computation.
     * Subclasses must call this method if they //override any of the
     * {@code solve} methods.
     *
     * @param f Function to solve.
     * @param min Lower bound for the interval.
     * @param max Upper bound for the interval.
     * @param start_value Start value to use.
     * @param max_eval Maximum number of evaluations.
     * @exception Null_Argument_Exception if f is NULL
     */
    void setup(const int& max_eval, const F& f, const double& min, const double& max, const double& start_value)
    {
        // Checks.
        //Math_Utils::check_not_null(f);

        // Reset.
        my_search_min = min;
        my_search_max = max;
        my_search_start = start_value;
        my_function = f;
        my_evaluations = my_evaluations.with_maximal_count(max_eval);
    }

    /**
     * Method for implementing actual optimization algorithms in derived
     * classes.
     *
     * @return the root.
     * @Math_Illegal_State_Exception if the maximal number of evaluations
     * is exceeded.
     * @ if the initial search interval does not bracket
     * a root and the solver requires it.
     */
    virtual double do_solve()

    /**
     * Check whether the function takes opposite signs at the endpoints.
     *
     * @param lower Lower endpoint.
     * @param upper Upper endpoint.
     * @return {@code true} if the function values have opposite signs at the
     * given points.
     */
    bool is_bracketing(const double& lower, const double& upper)
    {
        return Univariate_Solver_Utils::is_bracketing(my_function, lower, upper);
    }

    /**
     * Check whether the arguments form a (strictly) increasing sequence.
     *
     * @param start First number.
     * @param mid Second number.
     * @param end Third number.
     * @return {@code true} if the arguments form an increasing sequence.
     */
    bool is_sequence(const double& start, const double& mid, const double& end)
    {
        return Univariate_Solver_Utils::is_sequence(start, mid, end);
    }

    /**
     * Check that the endpoints specify an interval.
     *
     * @param lower Lower endpoint.
     * @param upper Upper endpoint.
     * @ if {@code lower >= upper}.
     */
    void verify_interval(const double& lower, const double& upper)
    {
        Univariate_Solver_Utils::verify_interval(lower, upper);
    }

    /**
     * Check that {@code lower < initial < upper}.
     *
     * @param lower Lower endpoint.
     * @param initial Initial value.
     * @param upper Upper endpoint.
     * @ if {@code lower >= initial} or
     * {@code initial >= upper}.
     */
    void verify_sequence(const double& lower, const double& initial, const double& upper)
    {
        Univariate_Solver_Utils::verify_sequence(lower, initial, upper);
    }

    /**
     * Check that the endpoints specify an interval and the function takes
     * opposite signs at the endpoints.
     *
     * @param lower Lower endpoint.
     * @param upper Upper endpoint.
     * @Null_Argument_Exception if the function has not been set.
     * @ if the function has the same sign at
     * the endpoints.
     */
    void verify_bracketing(const double& lower, const double& upper)
    {
        Univariate_Solver_Utils::verify_bracketing(function, lower, upper);
    }

    /**
     * Increment the evaluation count by one.
     * Method {@link #compute_objective_valuestatic_cast<double>(} calls this method internally.
     * It is provided for subclasses that do not exclusively use
     * {@code compute_objective_value} to solve the function.
     * See e.g. {@link AbstractUnivariate_Differentiable_Solver}.
     *
     * @Math_Illegal_State_Exception when the allowed number of function
     * evaluations has been exhausted.
     */
    void increment_evaluation_count()
    {
        my_evaluations.increment();
    }

public:
    /** {@inherit_doc} */
    //override
    int get_max_evaluations() const
    {
        return my_evaluations.get_maximal_count();
    }

    /** {@inherit_doc} */
    //override
    int get_evaluations() const
    {
        return my_evaluations.get_count();
    }

    /**
     * @return the lower end of the search interval.
     */
    double get_min() const
    {
        return my_search_min;
    }

    /**
     * @return the higher end of the search interval.
     */
    double get_max() const
    {
        return my_search_max;
    }

    /**
     * @return the initial guess.
     */
    double get_start_value() 
    {
        return my_search_start;
    }

    /**
     * {@inherit_doc}
     */
    //override
    double get_absolute_accuracy() const
    {
        return my_absolute_accuracy;
    }

    /**
     * {@inherit_doc}
     */
    //override
    double get_relative_accuracy() const
    {
        return my_relative_accuracy;
    }

    /**
     * {@inherit_doc}
     */
    //override
    double get_function_value_accuracy() const
    {
        return my_function_value_accuracy;
    }

    /** {@inherit_doc} */
    //override
    double solve(const int& max_eval, F f, const double& min,  const double& max,  const double& start_value)
    {
        // Initialization.
        setup(max_eval, f, min, max, start_value);

        // Perform computation.
        return do_solve();
    }

    /** {@inherit_doc} */
    //override
    double solve(const int& max_eval, const F& f, const double& min, const double& max) 
    {
        return solve(max_eval, f, min, max, min + 0.5 * (max - min));
    }

    /** {@inherit_doc} */
    //override
    double solve(const int& max_eval, const F& f, const double& start_value)
    {
        return solve(max_eval, f,NAN,NAN, start_value);
    }
};