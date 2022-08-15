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
//package org.hipparchus.optim.nonlinear.scalar;

//import org.hipparchus.analysis.Multivariate_Function;
//import org.hipparchus.analysis.Univariate_Function;
//import org.hipparchus.analysis.function.Logit;
//import org.hipparchus.analysis.function.Sigmoid;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Utils;

/**
 * <p>Adapter for mapping bounded {@link Multivariate_Function} to unbounded ones.</p>
 *
 * <p>
 * This adapter can be used to wrap functions subject to simple bounds on
 * parameters so they can be used by optimizers that do <em>not</em> directly
 * support simple bounds.
 * </p>
 * <p>
 * The principle is that the user function that will be wrapped will see its
 * parameters bounded as required, i.e when its {@code value} method is called
 * with argument array {@code point}, the elements array will fulfill requirement
 * {@code lower[i] <= point[i] <= upper[i]} for all i. Some of the components
 * may be unbounded or bounded only on one side if the corresponding bound is
 * set to an infinite value. The optimizer will not manage the user function by
 * itself, but it will handle this adapter and it is this adapter that will take
 * care the bounds are fulfilled. The adapter {@link #value(std::vector<double>)} method will
 * be called by the optimizer with unbound parameters, and the adapter will map
 * the unbounded value to the bounded range using appropriate functions like
 * {@link Sigmoid} for double bounded elements for example.
 * </p>
 * <p>
 * As the optimizer sees only unbounded parameters, it should be noted that the
 * start point or simplex expected by the optimizer should be unbounded, so the
 * user is responsible for converting his bounded point to unbounded by calling
 * {@link #bounded_to_unbounded(std::vector<double>)} before providing them to the optimizer.
 * For the same reason, the point returned by the {@link
 * org.hipparchus.optim.BaseMultivariate_Optimizer#optimize(Optimization_data[])}
 * method is unbounded. So to convert this point to bounded, users must call
 * {@link #unbounded_to_bounded(std::vector<double>)} by themselves!</p>
 * <p>
 * This adapter is only a poor man solution to simple bounds optimization constraints
 * that can be used with simple optimizers like
 * {@link org.hipparchus.optim.nonlinear.scalar.noderiv.Simplex_Optimizer
 * Simplex_Optimizer}.
 * A better solution is to use an optimizer that directly supports simple bounds like
 * {@link org.hipparchus.optim.nonlinear.scalar.noderiv.CMAES_Optimizer
 * CMAES_Optimizer} or
 * {@link org.hipparchus.optim.nonlinear.scalar.noderiv.BOBYQA_Optimizer
 * BOBYQA_Optimizer}.
 * One caveat of this poor-man's solution is that behavior near the bounds may be
 * numerically unstable as bounds are mapped from infinite values.
 * Another caveat is that convergence values are evaluated by the optimizer with
 * respect to unbounded variables, so there will be scales differences when
 * converted to bounded variables.
 * </p>
 *
 * @see Multivariate_FunctionPenalty_adapter
 *
 */
class Multivariate_FunctionMappingAdapter
    : Multivariate_Function 
    {
    /** Underlying bounded function. */
    private const Multivariate_Function bounded;
    /** Mapping functions. */
    private const Mapper[] mappers;

    /** Simple constructor.
     * @param bounded bounded function
     * @param lower lower bounds for each element of the input parameters array
     * (some elements may be set to {@code -INFINITY} for
     * unbounded values)
     * @param upper upper bounds for each element of the input parameters array
     * (some elements may be set to {@code INFINITY} for
     * unbounded values)
     * @exception  if lower and upper bounds are not
     * consistent, either according to dimension or to values
     */
    public Multivariate_FunctionMappingAdapter(const Multivariate_Function bounded, const std::vector<double> lower, const std::vector<double> upper) 
    {
        // safety checks
        //Math_Utils::check_not_null(lower);
        //Math_Utils::check_not_null(upper);
        if (lower.size() != upper.size()) 
        {
            throw (Localized_Core_Formats.DIMENSIONS_MISMATCH, lower.size(), upper.size());
        }
        for (int i{}; i < lower.size(); ++i) 
        {
            if (!(upper[i] >= lower[i])) { // NOPMD - the test is written this way so it also fails for NaN
                throw (Localized_Core_Formats.NUMBER_TOO_SMALL, upper[i], lower[i]);
            }
        }

        this.bounded = bounded;
        this.mappers = Mapper[lower.size()];
        for (int i{}; i < mappers.size(); ++i) 
        {
            if (Double.std::isinfinite(lower[i])) 
            {
                if (Double.std::isinfinite(upper[i])) 
                {
                    // element is unbounded, no transformation is needed
                    mappers[i] = NoBounds_mapper();
                }
else 
                {
                    // element is simple-bounded on the upper side
                    mappers[i] = Upper_Bound_Mapper(upper[i]);
                }
            }
else 
            {
                if (Double.std::isinfinite(upper[i])) 
                {
                    // element is simple-bounded on the lower side
                    mappers[i] = Lower_Bound_Mapper(lower[i]);
                }
else 
                {
                    // element is double-bounded
                    mappers[i] = LowerUpper_Bound_Mapper(lower[i], upper[i]);
                }
            }
        }
    }

    /**
     * Maps an array from unbounded to bounded.
     *
     * @param point Unbounded values.
     * @return the bounded values.
     */
    public std::vector<double> unbounded_to_bounded(std::vector<double> point) 
    {
        // Map unbounded input point to bounded point.
        const std::vector<double> mapped = std::vector<double>(mappers.size()];
        for (int i{}; i < mappers.size(); ++i) 
        {
            mapped[i] = mappers[i].unbounded_to_bounded(point[i]);
        }

        return mapped;
    }

    /**
     * Maps an array from bounded to unbounded.
     *
     * @param point Bounded values.
     * @return the unbounded values.
     */
    public std::vector<double> bounded_to_unbounded(std::vector<double> point) 
    {
        // Map bounded input point to unbounded point.
        const std::vector<double> mapped = std::vector<double>(mappers.size()];
        for (int i{}; i < mappers.size(); ++i) 
        {
            mapped[i] = mappers[i].bounded_to_unbounded(point[i]);
        }

        return mapped;
    }

    /**
     * Compute the underlying function value from an unbounded point.
     * <p>
     * This method simply bounds the unbounded point using the mappings
     * set up at construction and calls the underlying function using
     * the bounded point.
     * </p>
     * @param point unbounded value
     * @return underlying function value
     * @see #unbounded_to_bounded(std::vector<double>)
     */
    //override
    public double value(std::vector<double> point) 
    {
        return bounded.value(unbounded_to_bounded(point));
    }

    /** Mapping interface. */
    private interface Mapper 
    {
        /**
         * Maps a value from unbounded to bounded.
         *
         * @param y Unbounded value.
         * @return the bounded value.
         */
        double unbounded_to_bounded(double y);

        /**
         * Maps a value from bounded to unbounded.
         *
         * @param x Bounded value.
         * @return the unbounded value.
         */
        double bounded_to_unbounded(double x);
    }

    /** Local class for no bounds mapping. */
    private static class NoBounds_mapper : Mapper 
    {
        /** {@inherit_doc} */
        //override
        public double unbounded_to_bounded(const double y) 
        {
            return y;
        }

        /** {@inherit_doc} */
        //override
        public double bounded_to_unbounded(const double& x) 
        {
            return x;
        }
    }

    /** Local class for lower bounds mapping. */
    private static class Lower_Bound_Mapper : Mapper 
    {
        /** Low bound. */
        private const double lower;

        /**
         * Simple constructor.
         *
         * @param lower lower bound
         */
        Lower_Bound_Mapper(const double lower) 
        {
            this.lower = lower;
        }

        /** {@inherit_doc} */
        //override
        public double unbounded_to_bounded(const double y) 
        {
            return lower + std::exp(y);
        }

        /** {@inherit_doc} */
        //override
        public double bounded_to_unbounded(const double& x) 
        {
            return std::log(x - lower);
        }

    }

    /** Local class for upper bounds mapping. */
    private static class Upper_Bound_Mapper : Mapper 
    {

        /** Upper bound. */
        private const double upper;

        /** Simple constructor.
         * @param upper upper bound
         */
        Upper_Bound_Mapper(const double upper) 
        {
            this.upper = upper;
        }

        /** {@inherit_doc} */
        //override
        public double unbounded_to_bounded(const double y) 
        {
            return upper - std::exp(-y);
        }

        /** {@inherit_doc} */
        //override
        public double bounded_to_unbounded(const double& x) 
        {
            return -std::log(upper - x);
        }

    }

    /** Local class for lower and bounds mapping. */
    private static class LowerUpper_Bound_Mapper : Mapper 
    {
        /** Function from unbounded to bounded. */
        private const Univariate_Function bounding_function;
        /** Function from bounded to unbounded. */
        private const Univariate_Function unbounding_function;

        /**
         * Simple constructor.
         *
         * @param lower lower bound
         * @param upper upper bound
         */
        LowerUpper_Bound_Mapper(const double lower, const double upper) 
        {
            bounding_function   = Sigmoid(lower, upper);
            unbounding_function = Logit(lower, upper);
        }

        /** {@inherit_doc} */
        //override
        public double unbounded_to_bounded(const double y) 
        {
            return bounding_function.value(y);
        }

        /** {@inherit_doc} */
        //override
        public double bounded_to_unbounded(const double& x) 
        {
            return unbounding_function.value(x);
        }
    }
}


