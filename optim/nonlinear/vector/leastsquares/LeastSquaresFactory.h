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
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.linear.Array_2D_Row_Real_Matrix;
//import org.hipparchus.linear.Array_Real_Vector;
//import org.hipparchus.linear.Diagonal_Matrix;
//import org.hipparchus.linear.Eigen_Decomposition;
//import org.hipparchus.linear.Real_Matrix;
//import org.hipparchus.linear.Real_Vector;
//import org.hipparchus.optim.AbstractOptimization_Problem;
//import org.hipparchus.optim.Convergence_Checker;
//import org.hipparchus.optim.Localized_Optim_Formats;
//import org.hipparchus.optim.Point_Vector_Value_Pair;
//import org.hipparchus.optim.nonlinear.vector.leastsquares.Least_Squares_Problem.Evaluation;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Incrementor;
//import org.hipparchus.util.Pair;

/**
 * A Factory for creating {@link Least_Squares_Problem}s.
 *
 */
class Least_Squares_Factory 
{

    /** Prevent instantiation. */
    private Least_Squares_Factory() {}

    /**
     * Create a {@link org.hipparchus.optim.nonlinear.vector.leastsquares.Least_Squares_Problem}
     * from the given elements. There will be no weights applied (unit weights).
     *
     * @param model          the model function. Produces the computed values.
     * @param observed       the observed (target) values
     * @param start          the initial guess.
     * @param weight         the weight matrix
     * @param checker        convergence checker
     * @param max_evaluations the maximum number of times to evaluate the model
     * @param max_iterations  the maximum number to times to iterate in the algorithm
     * @param lazy_evaluation Whether the call to {@link Evaluation#evaluate(Real_Vector)}
     * will defer the evaluation until access to the value is requested.
     * @param param_validator Model parameters validator.
     * @return the specified General Least Squares problem.
     *
     */
    public static Least_Squares_Problem create(const Multivariate_Jacobian_Function model, const Real_Vector observed, const Real_Vector start, const Real_Matrix weight, const Convergence_Checker<Evaluation> checker, const int max_evaluations, const int max_iterations, const bool lazy_evaluation, const Parameter_Validator param_validator) 
    {
        const Least_Squares_Problem p = LocalLeast_Squares_Problem(model, observed, start, checker, max_evaluations, max_iterations, lazy_evaluation, param_validator);
        if (weight != NULL) 
        {
            return weight_matrix(p, weight);
        }
else 
        {
            return p;
        }
    }

    /**
     * Create a {@link org.hipparchus.optim.nonlinear.vector.leastsquares.Least_Squares_Problem}
     * from the given elements. There will be no weights applied (unit weights).
     *
     * @param model          the model function. Produces the computed values.
     * @param observed       the observed (target) values
     * @param start          the initial guess.
     * @param checker        convergence checker
     * @param max_evaluations the maximum number of times to evaluate the model
     * @param max_iterations  the maximum number to times to iterate in the algorithm
     * @return the specified General Least Squares problem.
     */
    public static Least_Squares_Problem create(const Multivariate_Jacobian_Function model, const Real_Vector observed, const Real_Vector start, const Convergence_Checker<Evaluation> checker, const int max_evaluations, const int max_iterations) 
    {
        return create(model, observed, start, NULL, checker, max_evaluations, max_iterations, false, NULL);
    }

    /**
     * Create a {@link org.hipparchus.optim.nonlinear.vector.leastsquares.Least_Squares_Problem}
     * from the given elements.
     *
     * @param model          the model function. Produces the computed values.
     * @param observed       the observed (target) values
     * @param start          the initial guess.
     * @param weight         the weight matrix
     * @param checker        convergence checker
     * @param max_evaluations the maximum number of times to evaluate the model
     * @param max_iterations  the maximum number to times to iterate in the algorithm
     * @return the specified General Least Squares problem.
     */
    public static Least_Squares_Problem create(const Multivariate_Jacobian_Function model, const Real_Vector observed, const Real_Vector start, const Real_Matrix weight, const Convergence_Checker<Evaluation> checker, const int max_evaluations, const int max_iterations) 
    {
        return weight_matrix(create(model, observed, start, checker, max_evaluations, max_iterations), weight);
    }

    /**
     * Create a {@link org.hipparchus.optim.nonlinear.vector.leastsquares.Least_Squares_Problem}
     * from the given elements.
     * <p>
     * This factory method is provided for continuity with previous interfaces. Newer
     * applications should use {@link #create(Multivariate_Jacobian_Function, Real_Vector, * Real_Vector, Convergence_Checker, int, int)}, or {@link #create(Multivariate_Jacobian_Function, * Real_Vector, Real_Vector, Real_Matrix, Convergence_Checker, int, int)}.
     *
     * @param model          the model function. Produces the computed values.
     * @param jacobian       the jacobian of the model with respect to the parameters
     * @param observed       the observed (target) values
     * @param start          the initial guess.
     * @param weight         the weight matrix
     * @param checker        convergence checker
     * @param max_evaluations the maximum number of times to evaluate the model
     * @param max_iterations  the maximum number to times to iterate in the algorithm
     * @return the specified General Least Squares problem.
     */
    public static Least_Squares_Problem create(const MultivariateVector_function model, const Multivariate_Matrix_Function jacobian, const std::vector<double> observed, const std::vector<double> start, const Real_Matrix weight, const Convergence_Checker<Evaluation> checker, const int max_evaluations, const int max_iterations) 
    {
        return create(model(model, jacobian), Array_Real_Vector(observed, false), Array_Real_Vector(start, false), weight, checker, max_evaluations, max_iterations);
    }

    /**
     * Apply a dense weight matrix to the {@link Least_Squares_Problem}.
     *
     * @param problem the unweighted problem
     * @param weights the matrix of weights
     * @return a {@link Least_Squares_Problem} with the weights applied. The original
     *         {@code problem} is not modified.
     */
    public static Least_Squares_Problem weight_matrix(const Least_Squares_Problem problem, const Real_Matrix weights) 
    {
        const Real_Matrix weight_square_root = square_root(weights);
        return Least_Squares_Adapter(problem) 
        {
            /** {@inherit_doc} */
            //override
            public Evaluation evaluate(const Real_Vector point) 
            {
                return DenseWeighted_Evaluation(super.evaluate(point), weight_square_root);
            }
        };
    }

    /**
     * Apply a diagonal weight matrix to the {@link Least_Squares_Problem}.
     *
     * @param problem the unweighted problem
     * @param weights the diagonal of the weight matrix
     * @return a {@link Least_Squares_Problem} with the weights applied. The original
     *         {@code problem} is not modified.
     */
    public static Least_Squares_Problem weight_diagonal(const Least_Squares_Problem problem, const Real_Vector weights) 
    {
        // TODO more efficient implementation
        return weight_matrix(problem, Diagonal_Matrix(weights.to_array()));
    }

    /**
     * Count the evaluations of a particular problem. The {@code counter} will be
     * incremented every time {@link Least_Squares_Problem#evaluate(Real_Vector)} is called on
     * the <em>returned</em> problem.
     *
     * @param problem the problem to track.
     * @param counter the counter to increment.
     * @return a least squares problem that tracks evaluations
     */
    public static Least_Squares_Problem count_evaluations(const Least_Squares_Problem problem, const Incrementor counter) 
    {
        return Least_Squares_Adapter(problem) 
        {

            /** {@inherit_doc} */
            //override
            public Evaluation evaluate(const Real_Vector point) 
            {
                counter.increment();
                return super.evaluate(point);
            }

            // Delegate the rest.
        };
    }

    /**
     * View a convergence checker specified for a {@link Point_Vector_Value_Pair} as one
     * specified for an {@link Evaluation}.
     *
     * @param checker the convergence checker to adapt.
     * @return a convergence checker that delegates to {@code checker}.
     */
    public static Convergence_Checker<Evaluation> evaluation_checker(const Convergence_Checker<Point_Vector_Value_Pair> checker) 
    {
        return Convergence_Checker<Evaluation>() 
        {
            /** {@inherit_doc} */
            //override
            public bool converged(const int iteration, const Evaluation previous, const Evaluation current) 
            {
                return checker.converged(
                        iteration, Point_Vector_Value_Pair(
                                previous.get_point().to_array(), previous.get_residuals().to_array(), false), Point_Vector_Value_Pair(
                                current.get_point().to_array(), current.get_residuals().to_array(), false)
                );
            }
        };
    }

    /**
     * Computes the square-root of the weight matrix.
     *
     * @param m Symmetric, positive-definite (weight) matrix.
     * @return the square-root of the weight matrix.
     */
    private static Real_Matrix square_root(const Real_Matrix& m) 
    {
        if (m instanceof Diagonal_Matrix) 
        {
            const int dim = m.get_row_dimension();
            const Real_Matrix sqrt_m = Diagonal_Matrix(dim);
            for (int i{}; i < dim; i++) 
            {
                sqrt_m.set_entry(i, i, std::sqrt(m.get_entry(i, i)));
            }
            return sqrt_m;
        }
else 
        {
            const Eigen_Decomposition dec = Eigen_Decomposition(m);
            return dec.get_square_root();
        }
    }

    /**
     * Combine a {@link MultivariateVector_function} with a {@link
     * Multivariate_Matrix_Function} to produce a {@link Multivariate_Jacobian_Function}.
     *
     * @param value    the vector value function
     * @param jacobian the Jacobian function
     * @return a function that computes both at the same time
     */
    public static Multivariate_Jacobian_Function model(const MultivariateVector_function value, const Multivariate_Matrix_Function jacobian) 
    {
        return LocalValue_And_Jacobian_Function(value, jacobian);
    }

    /**
     * Combine a {@link MultivariateVector_function} with a {@link
     * Multivariate_Matrix_Function} to produce a {@link Multivariate_Jacobian_Function}.
     */
    private static class LocalValue_And_Jacobian_Function
        : Value_And_Jacobian_Function 
        {
        /** Model. */
        private const MultivariateVector_function value;
        /** Model's Jacobian. */
        private const Multivariate_Matrix_Function jacobian;

        /**
         * @param value Model function.
         * @param jacobian Model's Jacobian function.
         */
        LocalValue_And_Jacobian_Function(const MultivariateVector_function value, const Multivariate_Matrix_Function jacobian) 
        {
            this.value = value;
            this.jacobian = jacobian;
        }

        /** {@inherit_doc} */
        //override
        public Pair<Real_Vector, Real_Matrix> value(const Real_Vector point) 
        {
            //TODO get array from Real_Vector without copying?
            const std::vector<double> p = point.to_array();

            // Evaluate.
            return Pair<Real_Vector, Real_Matrix>(compute_value(p), compute_jacobian(p));
        }

        /** {@inherit_doc} */
        //override
        public Real_Vector compute_value(const std::vector<double> params) 
        {
            return Array_Real_Vector(value.value(params), false);
        }

        /** {@inherit_doc} */
        //override
        public Real_Matrix compute_jacobian(const std::vector<double> params) 
        {
            return Array_2D_Row_Real_Matrix(jacobian.value(params), false);
        }
    }


    /**
     * A private, "field" immutable (not "real" immutable) implementation of {@link
     * Least_Squares_Problem}.
     */
    private static class LocalLeast_Squares_Problem
            extends AbstractOptimization_Problem<Evaluation>
            : Least_Squares_Problem 
            {

        /** Target values for the model function at optimum. */
        private const Real_Vector target;
        /** Model function. */
        private const Multivariate_Jacobian_Function model;
        /** Initial guess. */
        private const Real_Vector start;
        /** Whether to use lazy evaluation. */
        private const bool lazy_evaluation;
        /** Model parameters validator. */
        private const Parameter_Validator param_validator;

        /**
         * Create a {@link Least_Squares_Problem} from the given data.
         *
         * @param model          the model function
         * @param target         the observed data
         * @param start          the initial guess
         * @param checker        the convergence checker
         * @param max_evaluations the allowed evaluations
         * @param max_iterations  the allowed iterations
         * @param lazy_evaluation Whether the call to {@link Evaluation#evaluate(Real_Vector)}
         * will defer the evaluation until access to the value is requested.
         * @param param_validator Model parameters validator.
         */
        LocalLeast_Squares_Problem(const Multivariate_Jacobian_Function model, const Real_Vector target, const Real_Vector start, const Convergence_Checker<Evaluation> checker, const int max_evaluations, const int max_iterations, const bool lazy_evaluation, const Parameter_Validator param_validator) 
        {
            super(max_evaluations, max_iterations, checker);
            this.target = target;
            this.model = model;
            this.start = start;
            this.lazy_evaluation = lazy_evaluation;
            this.param_validator = param_validator;

            if (lazy_evaluation &&
                !(model instanceof Value_And_Jacobian_Function)) 
                {
                // Lazy evaluation requires that value and Jacobian
                // can be computed separately.
                throw Math_Illegal_State_Exception(Localized_Optim_Formats.INVALID_IMPLEMENTATION, model.get_class().get_name());
            }
        }

        /** {@inherit_doc} */
        //override
        public int get_observation_size() 
        {
            return target.get_dimension();
        }

        /** {@inherit_doc} */
        //override
        public int get_parameter_size() 
        {
            return start.get_dimension();
        }

        /** {@inherit_doc} */
        //override
        public Real_Vector get_start() 
        {
            return start == NULL ? NULL : start.copy();
        }

        /** {@inherit_doc} */
        //override
        public Evaluation evaluate(const Real_Vector point) 
        {
            // Copy so optimizer can change point without changing our instance.
            const Real_Vector p = param_validator == NULL ?
                point.copy() :
                param_validator.validate(point.copy());

            if (lazy_evaluation) 
            {
                return Lazy_Unweighted_Evaluation((Value_And_Jacobian_Function) model, target, p);
            }
else 
            {
                // Evaluate value and jacobian in one function call.
                const Pair<Real_Vector, Real_Matrix> value = model.value(p);
                return Unweighted_Evaluation(value.get_first(), value.get_second(), target, p);
            }
        }

        /**
         * Container with the model evaluation at a particular point.
         */
        private static class Unweighted_Evaluation extends Abstract_Evaluation 
        {
            /** Point of evaluation. */
            private const Real_Vector point;
            /** Derivative at point. */
            private const Real_Matrix jacobian;
            /** Computed residuals. */
            private const Real_Vector residuals;

            /**
             * Create an {@link Evaluation} with no weights.
             *
             * @param values   the computed function values
             * @param jacobian the computed function Jacobian
             * @param target   the observed values
             * @param point    the abscissa
             */
            private Unweighted_Evaluation(const Real_Vector values, const Real_Matrix jacobian, const Real_Vector target, const Real_Vector point) 
            {
                super(target.get_dimension());
                this.jacobian = jacobian;
                this.point = point;
                this.residuals = target.subtract(values);
            }

            /** {@inherit_doc} */
            //override
            public Real_Matrix get_jacobian() 
            {
                return jacobian;
            }

            /** {@inherit_doc} */
            //override
            public Real_Vector get_point() 
            {
                return point;
            }

            /** {@inherit_doc} */
            //override
            public Real_Vector get_residuals() 
            {
                return residuals;
            }
        }

        /**
         * Container with the model <em>lazy</em> evaluation at a particular point.
         */
        private static class Lazy_Unweighted_Evaluation extends Abstract_Evaluation 
        {
            /** Point of evaluation. */
            private const Real_Vector point;
            /** Model and Jacobian functions. */
            private const Value_And_Jacobian_Function model;
            /** Target values for the model function at optimum. */
            private const Real_Vector target;

            /**
             * Create an {@link Evaluation} with no weights.
             *
             * @param model  the model function
             * @param target the observed values
             * @param point  the abscissa
             */
            private Lazy_Unweighted_Evaluation(const Value_And_Jacobian_Function model, const Real_Vector target, const Real_Vector point) 
            {
                super(target.get_dimension());
                // Safe to cast as long as we control usage of this class.
                this.model = model;
                this.point = point;
                this.target = target;
            }

            /** {@inherit_doc} */
            //override
            public Real_Matrix get_jacobian() 
            {
                return model.compute_jacobian(point.to_array());
            }

            /** {@inherit_doc} */
            //override
            public Real_Vector get_point() 
            {
                return point;
            }

            /** {@inherit_doc} */
            //override
            public Real_Vector get_residuals() 
            {
                return target.subtract(model.compute_value(point.to_array()));
            }
        }
    }
}



