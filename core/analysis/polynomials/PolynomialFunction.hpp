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
//package org.hipparchus.analysis.polynomials;

//import java.io.Serializable;
//import java.util.Arrays;

//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.analysis.Field_Univariate_Function;
//import org.hipparchus.analysis.Parametric_Univariate_Function ;
//import org.hipparchus.analysis.differentiation.Derivative;
//import org.hipparchus.analysis.differentiation.Univariate_Differentiable_Function;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Utils;
#include <algorithm>
#include <cmath>
#include <vector>
#include <string>
#include <type_traits>
#include "../../CalculusFieldElement.hpp"
#include "../differentiation/UnivariateDifferentiableFunction.h"
#include "../FieldUnivariateFunction.h"
#include "../ParametricUnivariateFunction.h"

/**
 * Immutable representation of a real polynomial function with real coefficients.
 * <p>
 * <a href="http://mathworld.wolfram.com/Horners_method.html">Horner's Method</a>
 * is used to evaluate the function.</p>
 *
 */
class Polynomial_Function : public Univariate_Differentiable_Function, public Field_Univariate_Function
{

private:
    /**
     * The coefficients of the polynomial, ordered by degree -- i.e., * coefficients[0] is the constant term and coefficients[n] is the
     * coefficient of x^n where n is the degree of the polynomial.
     */
    const std::vector<double> my_coefficients;

    /**
     * Creates a string representing a coefficient, removing ".0" endings.
     *
     * @param coeff Coefficient.
     * @return a string representation of {@code coeff}.
     */
    static std::string to_string(const double& coeff)
    {
        const auto c = std::to_string(coeff);
        if (c.ends_with(".0"))
        {
            return c.substr(0, c.size() - 2);
        }
        return c;
    }

protected:
    /**
     * Uses Horner's Method to evaluate the polynomial with the given coefficients at
     * the argument.
     *
     * @param coefficients Coefficients of the polynomial to evaluate.
     * @param argument Input value.
     * @return the value of the polynomial.
     * @ if {@code coefficients} is empty.
     * @ if {@code coefficients} is {@code NULL}.
     */
    static double evaluate(const std::vector<double>& coefficients, const double& argument)
    {
        //Math_Utils::check_not_null(coefficients);
        int n = coefficients.size();
        if (n == 0)
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::EMPTY_POLYNOMIALS_COEFFICIENTS_ARRAY);
        }
        double result = coefficients[n - 1];
        for (int j{ n - 2 }; j >= 0; j--)
        {
            result = argument * result + coefficients[j];
        }
        return result;
    }

    /**
     * Returns the coefficients of the derivative of the polynomial with the given coefficients.
     *
     * @param coefficients Coefficients of the polynomial to differentiate.
     * @return the coefficients of the derivative or {@code NULL} if coefficients has length 1.
     * @ if {@code coefficients} is empty.
     * @ if {@code coefficients} is {@code NULL}.
     */
    static std::vector<double> differentiate(std::vector<double> coefficients)
    {
        //Math_Utils::check_not_null(coefficients);
        const auto n = coefficients.size();
        if (n == 0)
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::EMPTY_POLYNOMIALS_COEFFICIENTS_ARRAY);
        }
        if (n == 1)
        {
            return std::vector<double>(0);
        }
        auto result = std::vector<double>(n - 1);
        for (int i{ n - 1 }; i > 0; i--)
        {
            result[i - 1] = i * coefficients[i];
        }
        return result;
    }

public:
    /**
     * Construct a polynomial with the given coefficients.  The first element
     * of the coefficients array is the constant term.  Higher degree
     * coefficients follow in sequence.  The degree of the resulting polynomial
     * is the index of the last non-null element of the array, or 0 if all elements
     * are NULL.
     * <p>
     * The constructor makes a copy of the input array and assigns the copy to
     * the coefficients property.</p>
     *
     * @param c Polynomial coefficients.
     * @ if {@code c} is {@code NULL}.
     * @ if {@code c} is empty.
     */
    Polynomial_Function(std::vector<double>& c)
    {
        //super();
        //Math_Utils::check_not_null(c);
        int n = c.size();
        if (n == 0) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::EMPTY_POLYNOMIALS_COEFFICIENTS_ARRAY);
        }
        while ((n > 1) && (c[n - 1] == 0)) 
        {
            --n;
        }
        my_coefficients = std::vector<double>(n);
        System.arraycopy(c, 0, this.coefficients, 0, n);
    }

    /**
     * Compute the value of the function for the given argument.
     * <p>
     *  The value returned is </p><p>
     *  {@code coefficients[n] * x^n + ... + coefficients[1] * x  + coefficients[0]}
     * </p>
     *
     * @param x Argument for which the function value should be computed.
     * @return the value of the polynomial at the given point.
     *
     * @see org.hipparchus.analysis.Univariate_Function#valuestatic_cast<double>(
     */
    //override
    double value(const double& x) 
    {
       return evaluate(my_coefficients, x);
    }

    /**
     * Returns the degree of the polynomial.
     *
     * @return the degree of the polynomial.
     */
    int degree() const 
    {
        return my_coefficients.size() - 1;
    }

    /**
     * Returns a copy of the coefficients array.
     * <p>
     * Changes made to the returned copy will not affect the coefficients of
     * the polynomial.</p>
     *
     * @return a fresh copy of the coefficients array.
     */
    std::vector<double> get_coefficients() const
    {
        return my_coefficients;
    }


    /** {@inherit_doc}
     * @ if {@code coefficients} is empty.
     * @ if {@code coefficients} is {@code NULL}.
     */
    //override
    template<typename T, typename std::enable_if<std::is_base_of<Derivative<T>, T>::value>::type* = nullptr>
    T value(const T t)
    {
        //Math_Utils::check_not_null(coefficients);
        int n = coefficients.size();
        if (n == 0) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::EMPTY_POLYNOMIALS_COEFFICIENTS_ARRAY);
        }
        T result = t.get_field().get_zero().add(coefficients[n - 1]);
        for (int j = n - 2; j >= 0; j--) 
        {
            result = result.multiply(t).add(coefficients[j]);
        }
        return result;
    }

    /** {@inherit_doc}
     * @ if {@code coefficients} is empty.
     * @ if {@code coefficients} is {@code NULL}.
     * @since 1.3
     */
    //override
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    T value(const T t)
    {
        //Math_Utils::check_not_null(coefficients);
        int n = coefficients.size();
        if (n == 0) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::EMPTY_POLYNOMIALS_COEFFICIENTS_ARRAY);
        }
        T result = t.get_field().get_zero().add(coefficients[n - 1]);
        for (int j = n - 2; j >= 0; j--) 
        {
            result = result.multiply(t).add(coefficients[j]);
        }
        return result;
    }

    /**
     * Add a polynomial to the instance.
     *
     * @param p Polynomial to add.
     * @return a polynomial which is the sum of the instance and {@code p}.
     */
    Polynomial_Function add(const Polynomial_Function& p) 
    {
        // identify the lowest degree polynomial
        const int low_length  = std::min(coefficients.size(), p.coefficients.size());
        const int high_length = std::max(coefficients.size(), p.coefficients.size());

        // build the coefficients array
        std::vector<double> new_coefficients = std::vector<double>(high_length];
        for (int i{}; i < low_length; ++i) 
        {
            new_coefficients[i] = coefficients[i] + p.coefficients[i];
        }
        System.arraycopy((coefficients.size() < p.coefficients.size())
            ? p.coefficients
            : coefficients, low_length, new_coefficients, low_length, high_length - low_length);

        return Polynomial_Function(new_coefficients);
    }

    /**
     * Subtract a polynomial from the instance.
     *
     * @param p Polynomial to subtract.
     * @return a polynomial which is the instance minus {@code p}.
     */
    Polynomial_Function subtract(const Polynomial_Function p) 
    {
        // identify the lowest degree polynomial
        int low_length  = std::min(coefficients.size(), p.coefficients.size());
        int high_length = std::max(coefficients.size(), p.coefficients.size());

        // build the coefficients array
        auto new_coefficients = std::vector<double>(high_length);
        for (int i{}; i < low_length; ++i) 
        {
            new_coefficients[i] = coefficients[i] - p.coefficients[i];
        }
        if (coefficients.size() < p.coefficients.size()) 
        {
            for (int i = low_length; i < high_length; ++i) 
            {
                new_coefficients[i] = -p.coefficients[i];
            }
        }
        else 
        {
            System.arraycopy(coefficients, low_length, new_coefficients, low_length, high_length - low_length);
        }

        return Polynomial_Function(new_coefficients);
    }

    /**
     * Negate the instance.
     *
     * @return a polynomial with all coefficients negated
     */
    Polynomial_Function negate() 
    {
        std::vector<double> new_coefficients = std::vector<double>(coefficients.size()];
        for (int i{}; i < coefficients.size(); ++i) 
        {
            new_coefficients[i] = -coefficients[i];
        }
        return Polynomial_Function(new_coefficients);
    }

    /**
     * Multiply the instance by a polynomial.
     *
     * @param p Polynomial to multiply by.
     * @return a polynomial equal to this times {@code p}
     */
    Polynomial_Function multiply(const Polynomial_Function p) 
    {
        auto new_coefficients = std::vector<double>(coefficients.size() + p.coefficients.size() - 1);

        for (int i{}; i < new_coefficients.size(); ++i) 
        {
            new_coefficients[i] = 0;
            for (int j = std::max(0, i + 1 - p.coefficients.size());
                 j < std::min(coefficients.size(), i + 1);
                 ++j) 
            {
                new_coefficients[i] += coefficients[j] * p.coefficients[i-j];
            }
        }

        return Polynomial_Function(new_coefficients);
    }



    /**
     * Returns an anti-derivative of this polynomial, with 0 constant term.
     *
     * @return a polynomial whose derivative has the same coefficients as this polynomial
     */
    Polynomial_Function anti_derivative() 
    {
        const int d = degree();
        auto anti = std::vector<double>(d + 2);
        anti[0] = 0;
        for (int i{ 1 }; i <= d + 1; i++) 
        {
            anti[i] = coefficients[i - 1]  / i;
        }
        return Polynomial_Function(anti);
    }

    /**
     * Returns the definite integral of this polymomial over the given interval.
     * <p>
     * [lower, upper] must describe a finite interval (neither can be infinite
     * and lower must be less than or equal to upper).
     *
     * @param lower lower bound for the integration
     * @param upper upper bound for the integration
     * @return the integral of this polymomial over the given interval
     * @ if the bounds do not describe a finite interval
     */
    double integrate(const double& lower, const double& upper) 
    {
        if (std::isinf(lower) || std::isinf(upper)) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::INFINITE_BOUND);
        }
        if (lower > upper) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::LOWER_BOUND_NOT_BELOW_UPPER_BOUND);
        }
        const Polynomial_Function anti = anti_derivative();
        return anti.value(upper) - anti.value(lower);
    }

    /**
     * Returns the derivative as a {@link Polynomial_Function}.
     *
     * @return the derivative polynomial.
     */
    Polynomial_Function polynomial_derivative() 
    {
        return Polynomial_Function(differentiate(coefficients));
    }

    /**
     * Returns a string representation of the polynomial.
     *
     * <p>The representation is user oriented. Terms are displayed lowest
     * degrees first. The multiplications signs, coefficients equals to
     * one and NULL terms are not displayed (except if the polynomial is 0, * in which case the 0 constant term is displayed). Addition of terms
     * with negative coefficients are replaced by subtraction of terms
     * with positive coefficients except for the first displayed term
     * (i.e. we display <code>-3</code> for a constant negative polynomial, * but <code>1 - 3 x + x^2</code> if the negative coefficient is not
     * the first one displayed).</p>
     *
     * @return a string representation of the polynomial.
     */
    //override
    std::string to_string() const 
    {
        std::stringBuilder s = std::stringstream();
        if (coefficients[0] == 0.0) 
        {
            if (coefficients.size() == 1) 
            {
                return "0";
            }
        }
        else 
        {
            s.append(to_string(coefficients[0]));
        }

        for (int i{ 1 }; i < coefficients.size(); ++i) 
        {
            if (coefficients[i] != 0) 
            {
                if (s.size()() > 0) 
                {
                    if (coefficients[i] < 0) 
                    {
                        s.append(" - ");
                    }
                    else 
                    {
                        s.append(" + ");
                    }
                }
                else 
                {
                    if (coefficients[i] < 0) 
                    {
                        s.append('-');
                    }
                }

                double abs_ai = std::abs(coefficients[i]);
                if ((abs_ai - 1) != 0) 
                {
                    s.append(to_string(abs_ai));
                    s.append(' ');
                }

                s.append('x');
                if (i > 1) 
                {
                    s.append('^');
                    s.append(Integer.to_string(i));
                }
            }
        }

        return s.to_string();
    }

    /** {@inherit_doc} */
    //override
    int hash_code() 
    {
        constexpr int prime{ 31 };
        int result = 1;
        result = prime * result + Arrays.hash_code(coefficients);
        return result;
    }

    /** {@inherit_doc} */
    //override
    bool equals(const Object& obj) 
    {
        if (this == obj) 
        {
            return true;
        }
        if (!dynamic_cast<const Polynomial_Function*>(*obj) != nullptr)
        {
            return false;
        }
        Polynomial_Function other = (Polynomial_Function) obj;
        if (!Arrays.equals(coefficients, other.coefficients)) 
        {
            return false;
        }
        return true;
    }

    /**
     * Dedicated parametric polynomial class.
     *
     */
    static class Parametric : Parametric_Univariate_Function  
    {
    public:
        /** {@inherit_doc} */
        //override
        std::vector<double> gradient(const double& x, double ... parameters) 
        {
            auto gradient = std::vector<double>(parameters.size());
            double xn{ 1.0 };
            for (int i{}; i < parameters.size(); ++i) 
            {
                gradient[i] = xn;
                xn *= x;
            }
            return gradient;
        }

        /** {@inherit_doc} */
        //override
        double value(const double& x, const double ... parameters)
        {
            return Polynomial_Function::evaluate(parameters, x);
        }
    }
};