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

//package org.hipparchus.analysis;

//import org.hipparchus.analysis.differentiation.DS_Factory;
//import org.hipparchus.analysis.differentiation.Derivative;
//import org.hipparchus.analysis.differentiation.Derivative_Structure;
//import org.hipparchus.analysis.differentiation.Multivariate_Differentiable_Function;
//import org.hipparchus.analysis.differentiation.Univariate_Differentiable_Function;
//import org.hipparchus.analysis.function.Identity;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.util.Math_Arrays;
//import org.hipparchus.util.Math_Utils;
#include <vector>

/**
 * Utilities for manipulating function objects.
 *
 */
class Function_Utils 
{
private:
    /**
     * Class only contains static methods.
     */
    Function_Utils() {}

public:
    /**
     * Composes functions.
     * <p>
     * The functions in the argument list are composed sequentially, in the
     * given order.  For example, compose(f1,f2,f3) acts like f1(f2(f3(x))).</p>
     *
     * @param f List of functions.
     * @return the composite function.
     */
    static Univariate_Function compose(const Univariate_Function ... f) 
    {
        return Univariate_Function() 
        {
            /** {@inherit_doc} */
            override
            public double value(double x) 
            {
                double r = x;
                for (int i = f.length - 1; i >= 0; i--) 
                {
                    r = f[i].value(r);
                }
                return r;
            }
        };
    }

    /**
     * Composes functions.
     * <p>
     * The functions in the argument list are composed sequentially, in the
     * given order.  For example, compose(f1,f2,f3) acts like f1(f2(f3(x))).</p>
     *
     * @param f List of functions.
     * @return the composite function.
     */
    static Univariate_Differentiable_Function compose(const Univariate_Differentiable_Function ... f) 
    {
        return Univariate_Differentiable_Function() 
        {

            /** {@inherit_doc} */
            override
            public double value(const double t) 
            {
                double r = t;
                for (int i = f.length - 1; i >= 0; i--) 
                {
                    r = f[i].value(r);
                }
                return r;
            }

            /** {@inherit_doc} */
            override
            public <T extends Derivative<T>> T value(const T t) 
            {
                T r = t;
                for (int i = f.length - 1; i >= 0; i--) 
                {
                    r = f[i].value(r);
                }
                return r;
            }

        };
    }

    /**
     * Adds functions.
     *
     * @param f List of functions.
     * @return a function that computes the sum of the functions.
     */
    public static Univariate_Function add(const Univariate_Function ... f) 
    {
        return Univariate_Function() 
        {
            /** {@inherit_doc} */
            override
            public double value(double x) 
            {
                double r = f[0].value(x);
                for (int i{ 1 }; i < f.length; i++) 
                {
                    r += f[i].value(x);
                }
                return r;
            }
        };
    }

    /**
     * Adds functions.
     *
     * @param f List of functions.
     * @return a function that computes the sum of the functions.
     */
    public static Univariate_Differentiable_Function add(const Univariate_Differentiable_Function ... f) 
    {
        return Univariate_Differentiable_Function() 
        {

            /** {@inherit_doc} */
            override
            public double value(const double t) 
            {
                double r = f[0].value(t);
                for (int i{ 1 }; i < f.length; i++) 
                {
                    r += f[i].value(t);
                }
                return r;
            }

            /** {@inherit_doc}
             * @ if functions are not consistent with each other
             */
            override
            public <T extends Derivative<T>> T value(const T t)
                 
                {
                T r = f[0].value(t);
                for (int i{ 1 }; i < f.length; i++) 
                {
                    r = r.add(f[i].value(t));
                }
                return r;
            }

        };
    }

    /**
     * Multiplies functions.
     *
     * @param f List of functions.
     * @return a function that computes the product of the functions.
     */
    public static Univariate_Function multiply(const Univariate_Function ... f) 
    {
        return Univariate_Function() 
        {
            /** {@inherit_doc} */
            override
            public double value(double x) 
            {
                double r = f[0].value(x);
                for (int i{ 1 }; i < f.length; i++) 
                {
                    r *= f[i].value(x);
                }
                return r;
            }
        };
    }

    /**
     * Multiplies functions.
     *
     * @param f List of functions.
     * @return a function that computes the product of the functions.
     */
    public static Univariate_Differentiable_Function multiply(const Univariate_Differentiable_Function ... f) 
    {
        return Univariate_Differentiable_Function() 
        {

            /** {@inherit_doc} */
            override
            public double value(const double t) 
            {
                double r = f[0].value(t);
                for (int i{ 1 }; i < f.length; i++) 
                {
                    r  *= f[i].value(t);
                }
                return r;
            }

            /** {@inherit_doc} */
            override
            public <T extends Derivative<T>> T value(const T t) 
            {
                T r = f[0].value(t);
                for (int i{ 1 }; i < f.length; i++) 
                {
                    r = r.multiply(f[i].value(t));
                }
                return r;
            }

        };
    }

    /**
     * Returns the univariate function
     * {@code h(x) = combiner(f(x), g(x)).}
     *
     * @param combiner Combiner function.
     * @param f Function.
     * @param g Function.
     * @return the composite function.
     */
    public static Univariate_Function combine(const Bivariate_Function combiner, const Univariate_Function& f, const Univariate_Function g) 
    {
        return Univariate_Function() 
        {
            /** {@inherit_doc} */
            override
            public double value(double x) 
            {
                return combiner.value(f.value(x), g.value(x));
            }
        };
    }

    /**
     * Returns a Multivariate_Function h(x[]) defined by <pre> <code>
     * h(x[]) = combiner(...combiner(combiner(initial_value,f(x[0])),f(x[1]))...),f(x[x.length-1]))
     * </code></pre>
     *
     * @param combiner Combiner function.
     * @param f Function.
     * @param initial_value Initial value.
     * @return a collector function.
     */
    public static Multivariate_Function collector(const Bivariate_Function combiner, const Univariate_Function& f, const double initial_value) 
    {
        return Multivariate_Function() 
        {
            /** {@inherit_doc} */
            override
            public double value(std::vector<double> point) 
            {
                double result = combiner.value(initial_value, f.value(point[0]));
                for (int i{ 1 }; i < point.length; i++) 
                {
                    result = combiner.value(result, f.value(point[i]));
                }
                return result;
            }
        };
    }

    /**
     * Returns a Multivariate_Function h(x[]) defined by <pre> <code>
     * h(x[]) = combiner(...combiner(combiner(initial_value,x[0]),x[1])...),x[x.length-1])
     * </code></pre>
     *
     * @param combiner Combiner function.
     * @param initial_value Initial value.
     * @return a collector function.
     */
    public static Multivariate_Function collector(const Bivariate_Function combiner, const double initial_value) 
    {
        return collector(combiner, Identity(), initial_value);
    }

    /**
     * Creates a unary function by fixing the first argument of a binary function.
     *
     * @param f Binary function.
     * @param fixed value to which the first argument of {@code f} is set.
     * @return the unary function h(x) = f(fixed, x)
     */
    public static Univariate_Function fix_1st_argument(const Bivariate_Function f, const double fixed) 
    {
        return Univariate_Function() 
        {
            /** {@inherit_doc} */
            override
            public double value(double x) 
            {
                return f.value(fixed, x);
            }
        };
    }
    /**
     * Creates a unary function by fixing the second argument of a binary function.
     *
     * @param f Binary function.
     * @param fixed value to which the second argument of {@code f} is set.
     * @return the unary function h(x) = f(x, fixed)
     */
    public static Univariate_Function fix2nd_argument(const Bivariate_Function f, const double fixed) 
    {
        return Univariate_Function() 
        {
            /** {@inherit_doc} */
            override
            public double value(double x) 
            {
                return f.value(x, fixed);
            }
        };
    }

    /**
     * Samples the specified univariate real function on the specified interval.
     * <p>
     * The interval is divided equally into {@code n} sections and sample points
     * are taken from {@code min} to {@code max - (max - min) / n}; therefore
     * {@code f} is not sampled at the upper bound {@code max}.</p>
     *
     * @param f Function to be sampled
     * @param min Lower bound of the interval (included).
     * @param max Upper bound of the interval (excluded).
     * @param n Number of sample points.
     * @return the array of samples.
     * @ if the lower bound {@code min} is
     * greater than, or equal to the upper bound {@code max}.
     * @ if the number of sample points
     * {@code n} is negative.
     */
    public static std::vector<double> sample(const Univariate_Function& f, const double& min, const double& max, const int& n) 
       {

        if (n <= 0) 
        {
            throw (
                    hipparchus::exception::Localized_Core_Formats_Type::NOT_POSITIVE_NUMBER_OF_SAMPLES, Integer.value_of(n));
        }
        if (min >= max) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_LARGE_BOUND_EXCLUDED, min, max);
        }

        const std::vector<double> s = std::vector<double>(n];
        const double h = (max - min) / n;
        for (int i{}; i < n; i++) 
        {
            s[i] = f.value(min + i * h);
        }
        return s;
    }

    /** Convert regular functions to {@link Univariate_Differentiable_Function}.
     * <p>
     * This method handle the case with one free parameter and several derivatives.
     * For the case with several free parameters and only first order derivatives, * see {@link #to_differentiable(Multivariate_Function, Multivariate_Vector_function)}.
     * There are no direct support for intermediate cases, with several free parameters
     * and order 2 or more derivatives, as is would be difficult to specify all the
     * cross derivatives.
     * </p>
     * <p>
     * Note that the derivatives are expected to be computed only with respect to the
     * raw parameter x of the base function, i.e. they are df/dx, df<sup>2</sup>/dx<sup>2</sup>, ...
     * Even if the built function is later used in a composition like f(sin(t)), the provided
     * derivatives should <em>not</em> apply the composition with sine and its derivatives by
     * themselves. The composition will be done automatically here and the result will properly
     * contain f(sin(t)), df(sin(t))/dt, df<sup>2</sup>(sin(t))/dt<sup>2</sup> despite the
     * provided derivatives functions know nothing about the sine function.
     * </p>
     * @param f base function f(x)
     * @param derivatives derivatives of the base function, in increasing differentiation order
     * @return a differentiable function with value and all specified derivatives
     * @see #to_differentiable(Multivariate_Function, Multivariate_Vector_function)
     * @see #derivative(Univariate_Differentiable_Function, int)
     */
    public static Univariate_Differentiable_Function to_differentiable(const Univariate_Function& f, const Univariate_Function ... derivatives) 
    {

        return Univariate_Differentiable_Function() 
        {

            /** {@inherit_doc} */
            //override
            public double value(const double& x) 
            {
                return f.value(x);
            }

            /** {@inherit_doc} */
            //override
            public <T extends Derivative<T>> T value(const T& x) 
            {
                if (x.get_order() > derivatives.length) 
                {
                    throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_LARGE, x.get_order(), derivatives.length);
                }
                const std::vector<double> packed = std::vector<double>(x.get_order() + 1];
                packed[0] = f.value(x.get_value());
                for (int i{}; i < x.get_order(); ++i) 
                {
                    packed[i + 1] = derivatives[i].value(x.get_value());
                }
                return x.compose(packed);
            }

        };

    }

    /** Convert regular functions to {@link Multivariate_Differentiable_Function}.
     * <p>
     * This method handle the case with several free parameters and only first order derivatives.
     * For the case with one free parameter and several derivatives, * see {@link #to_differentiable(Univariate_Function, Univariate_Function...)}.
     * There are no direct support for intermediate cases, with several free parameters
     * and order 2 or more derivatives, as is would be difficult to specify all the
     * cross derivatives.
     * </p>
     * <p>
     * Note that the gradient is expected to be computed only with respect to the
     * raw parameter x of the base function, i.e. it is df/dx<sub>1</sub>, df/dx<sub>2</sub>, ...
     * Even if the built function is later used in a composition like f(sin(t), cos(t)), the provided
     * gradient should <em>not</em> apply the composition with sine or cosine and their derivative by
     * itself. The composition will be done automatically here and the result will properly
     * contain f(sin(t), cos(t)), df(sin(t), cos(t))/dt despite the provided derivatives functions
     * know nothing about the sine or cosine functions.
     * </p>
     * @param f base function f(x)
     * @param gradient gradient of the base function
     * @return a differentiable function with value and gradient
     * @see #to_differentiable(Univariate_Function, Univariate_Function...)
     * @see #derivative(Multivariate_Differentiable_Function, std::vector<int>)
     */
    public static Multivariate_Differentiable_Function to_differentiable(const Multivariate_Function f, const Multivariate_Vector_function gradient) 
    {

        return Multivariate_Differentiable_Function() 
        {

            /** {@inherit_doc} */
            override
            public double value(const std::vector<double> point) 
            {
                return f.value(point);
            }

            /** {@inherit_doc} */
            override
            public Derivative_Structure value(const Derivative_Structure[] point) 
            {

                // set up the input parameters
                const std::vector<double> d_point = std::vector<double>(point.length];
                for (int i{}; i < point.length; ++i) 
                {
                    d_point[i] = point[i].get_value();
                    if (point[i].get_order() > 1) 
                    {
                        throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_LARGE, point[i].get_order(), 1);
                    }
                }

                // evaluate regular functions
                const double    v = f.value(d_point);
                const std::vector<double> dv = gradient.value(d_point);
                Math_Utils::check_dimension(dv.length, point.length);

                // build the combined derivative
                const int parameters = point[0].get_free_parameters();
                const std::vector<double> partials = std::vector<double>(point.length];
                const std::vector<double> packed = std::vector<double>(parameters + 1];
                packed[0] = v;
                const int orders[] = int[parameters];
                for (int i{}; i < parameters; ++i) 
                {

                    // we differentiate once with respect to parameter i
                    orders[i] = 1;
                    for (int j{}; j < point.length; ++j) 
                    {
                        partials[j] = point[j].get_partial_derivative(orders);
                    }
                    orders[i] = 0;

                    // compose partial derivatives
                    packed[i + 1] = Math_Arrays::linear_combination(dv, partials);

                }

                return point[0].get_factory().build(packed);

            }

        };

    }

    /** Convert an {@link Univariate_Differentiable_Function} to an
     * {@link Univariate_Function} computing n<sup>th</sup> order derivative.
     * <p>
     * This converter is only a convenience method. Beware computing only one derivative does
     * not save any computation as the original function will really be called under the hood.
     * The derivative will be extracted from the full {@link Derivative_Structure} result.
     * </p>
     * @param f original function, with value and all its derivatives
     * @param order of the derivative to extract
     * @return function computing the derivative at required order
     * @see #derivative(Multivariate_Differentiable_Function, std::vector<int>)
     * @see #to_differentiable(Univariate_Function, Univariate_Function...)
     */
    public static Univariate_Function derivative(const Univariate_Differentiable_Function f, const int order) 
    {

        const DS_Factory factory = DS_Factory(1, order);

        return Univariate_Function() 
        {

            /** {@inherit_doc} */
            override
            public double value(const double& x) 
            {
                const Derivative_Structure ds_x = factory.variable(0, x);
                return f.value(ds_x).get_partial_derivative(order);
            }

        };
    }

    /** Convert an {@link Multivariate_Differentiable_Function} to an
     * {@link Multivariate_Function} computing n<sup>th</sup> order derivative.
     * <p>
     * This converter is only a convenience method. Beware computing only one derivative does
     * not save any computation as the original function will really be called under the hood.
     * The derivative will be extracted from the full {@link Derivative_Structure} result.
     * </p>
     * @param f original function, with value and all its derivatives
     * @param orders of the derivative to extract, for each free parameters
     * @return function computing the derivative at required order
     * @see #derivative(Univariate_Differentiable_Function, int)
     * @see #to_differentiable(Multivariate_Function, Multivariate_Vector_function)
     */
    public static Multivariate_Function derivative(const Multivariate_Differentiable_Function f, const std::vector<int> orders) 
    {

        // the maximum differentiation order is the sum of all orders
        int sum{};
        for (const auto& order : orders) 
        {
            sum += order;
        }
        const int sum_orders = sum;

        return Multivariate_Function() 
        {

            /** Factory used for building derivatives. */
            private DS_Factory factory;

            /** {@inherit_doc} */
            //override
            public double value(const std::vector<double>& point) 
            {

                if (factory == NULL || point.length != factory.get_compiler().get_free_parameters()) 
                {
                    // rebuild the factory in case of mismatch
                    factory = DS_Factory(point.length, sum_orders);
                }

                // set up the input parameters
                auto ds_point = std::vector<Derivative_Structure>(point.length);
                for (int i{}; i < point.length; ++i) 
                {
                    ds_point[i] = factory.variable(i, point[i]);
                }

                return f.value(ds_point).get_partial_derivative(orders);

            }

        };
    }

};