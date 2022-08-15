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
//package org.hipparchus.analysis.differentiation;
#include <cmath>
#include <vector>
#include "FieldUnivariateDerivative2.h"

//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.Field;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Field_Sin_Cos;
//import org.hipparchus.util.Field_Sinh_Cosh;
//import org.hipparchus.util.Math_Arrays;
//import org.hipparchus.util.Math_Utils;
#include <type_traits>
#include "../../CalculusFieldElement.hpp"

/** Class representing both the value and the differentials of a function.
 * <p>This class is a stripped-down version of {@link Field_Derivative_Structure}
 * with only one {@link Field_Derivative_Structure#get_free_parameters() free parameter}
 * and {@link Field_Derivative_Structure#get_order() derivation order} also limited to one.
 * It should have less overhead than {@link Field_Derivative_Structure} in its domain.</p>
 * <p>This class is an implementation of Rall's numbers. Rall's numbers are an
 * extension to the real numbers used throughout mathematical expressions; they hold
 * the derivative together with the value of a function.</p>
 * <p>{@link Field_Univariate_Derivative_1} instances can be used directly thanks to
 * the arithmetic operators to the mathematical functions provided as
 * methods by this class (+, -, *, /, %, sin, cos ...).</p>
 * <p>Implementing complex expressions by hand using these classes is
 * a tedious and error-prone task but has the advantage of having no limitation
 * on the derivation order despite not requiring users to compute the derivatives by
 * themselves.</p>
 * <p>Instances of this class are guaranteed to be immutable.</p>
 * @param <T> the type of the function parameters and value
 * @see Derivative_Structure
 * @see Univariate_Derivative_1
 * @see Univariate_Derivative_2
 * @see Gradient
 * @see Field_Derivative_Structure
 * @see Field_Univariate_Derivative_2
 * @see Field_Gradient
 * @since 1.7
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
class Field_Univariate_Derivative_1 : public Field_Univariate_Derivative<T, Field_Univariate_Derivative_1<T>> 
{
private:
    /** Value of the function. */
    const T my_f0;

    /** First derivative of the function. */
    const T my_f1;

public:
    /** Build an instance with values and derivative.
     * @param f0 value of the function
     * @param f1 first derivative of the function
     */
    Field_Univariate_Derivative_1(const T& f0, const T& f1) : my_f0{ f0 }, my_f1{ f1 } {};

    /** Build an instance from a {@link Derivative_Structure}.
     * @param ds derivative structure
     * @exception  if either {@code ds} parameters
     * is not 1 or {@code ds} order is not 1
     */
    Field_Univariate_Derivative_1(const Field_Derivative_Structure<T>& ds)  
    {
        Math_Utils::check_dimension(ds.get_free_parameters(), 1);
        Math_Utils::check_dimension(ds.get_order(), 1);
        my_f0 = ds.get_value();
        my_f1 = ds.get_partial_derivative(1);
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> new_instance(const double value) 
    {
        const T zero = f0.get_field().get_zero();
        return Field_Univariate_Derivative_1<>(zero.new_instance(value), zero);
    }

    /** {@inherit_doc} */
    //override
    double get_real() const
    {
        return get_value().get_real();
    }

    /** Get the value part of the univariate derivative.
     * @return value part of the univariate derivative
     */
    //override
    T get_value() const
    {
        return my_f0;
    }

    /** Get a derivative from the univariate derivative.
     * @param n derivation order (must be between 0 and {@link #get_order()}, both inclusive)
     * @return n<sup>th</sup> derivative, or {@code NaN} if n is
     * either negative or strictly larger than {@link #get_order()}
     */
    //override
    T get_derivative(const int& n) 
    {
        switch (n) 
        {
            case 0 :
                return my_f0;
            case 1 :
                return my_f1;
            default :
                throw (hipparchus::exception::Localized_Core_Formats_Type::DERIVATION_ORDER_NOT_ALLOWED, n);
        }
    }

    /** Get the derivation order.
     * @return derivation order
     */
    //override
    int get_order() const
    {
        return 1;
    }

    /** Get the first derivative.
     * @return first derivative
     * @see #get_value()
     */
    T get_first_derivative() const
    {
        return my_f1;
    }

    /** Get the {@link Field} the value and parameters of the function belongs to.
     * @return {@link Field} the value and parameters of the function belongs to
     */
    Field<T> get_value_field() const
    {
        return my_f0.get_field();
    }

    /** Convert the instance to a {@link Field_Derivative_Structure}.
     * @return derivative structure with same value and derivative as the instance
     */
    //override
    Field_Derivative_Structure<T> to_derivative_structure() 
    {
        return get_field().get_conversion_factory().build(my_f0, my_f1);
    }

    /** '+' operator.
     * @param a right hand side parameter of the operator
     * @return this+a
     */
    Field_Univariate_Derivative_1<T> add(const T& a) const
    {
        return Field_Univariate_Derivative_1<>(f0.add(a), my_f1);
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> add(const double& a) const
    {
        return Field_Univariate_Derivative_1<>(my_f0.add(a), my_f1);
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> add(const Field_Univariate_Derivative_1<T>& a) 
    {
        return Field_Univariate_Derivative_1<>(my_f0.add(a.f0), f1.add(a.f1));
    }

    /** '-' operator.
     * @param a right hand side parameter of the operator
     * @return this-a
     */
    Field_Univariate_Derivative_1<T> subtract(const T& a) 
    {
        return Field_Univariate_Derivative_1<>(f0.subtract(a), my_f1);
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> subtract(const double& a)
    {
        return Field_Univariate_Derivative_1<>(f0.subtract(a), my_f1);
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> subtract(const Field_Univariate_Derivative_1<T>& a) 
    {
        return Field_Univariate_Derivative_1<>(my_f0.subtract(a.f0), my_f1.subtract(a.f1));
    }

    /** '&times;' operator.
     * @param a right hand side parameter of the operator
     * @return this&times;a
     */
    Field_Univariate_Derivative_1<T> multiply(const T& a) 
    {
        return Field_Univariate_Derivative_1<>(my_f0.multiply(a), my_f1.multiply(a));
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> multiply(const int& n) 
    {
        return Field_Univariate_Derivative_1<>(my_f0.multiply(n), my_f1.multiply(n));
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> multiply(const double& a) 
    {
        return Field_Univariate_Derivative_1<>(my_f0.multiply(a), my_f1.multiply(a));
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> multiply(const Field_Univariate_Derivative_1<T>& a) 
    {
        return Field_Univariate_Derivative_1<>(my_f0.multiply(a.f0), a.f0.linear_combination(f1, a.f0, my_f0, a.f1));
    }

    /** '&divide;' operator.
     * @param a right hand side parameter of the operator
     * @return this&divide;a
     */
    Field_Univariate_Derivative_1<T> divide(const T& a) 
    {
        const T inv1 = a.reciprocal();
        return Field_Univariate_Derivative_1<>(my_f0.multiply(inv1), my_f1.multiply(inv1));
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> divide(const double& a) 
    {
        const double inv1 = 1.0 / a;
        return Field_Univariate_Derivative_1<>(my_f0.multiply(inv1), my_f1.multiply(inv1));
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> divide(const Field_Univariate_Derivative_1<T>& a) 
    {
        const T inv1 = a.f0.reciprocal();
        const T inv2 = inv1.multiply(inv1);
        return Field_Univariate_Derivative_1<>(my_f0.multiply(inv1), a.f0.linear_combination(my_f1, a.f0, my_f0.negate(), a.f1).multiply(inv2));
    }

    /** IEEE remainder operator.
     * @param a right hand side parameter of the operator
     * @return this - n &times; a where n is the closest integer to this/a
     * (the even integer is chosen for n if this/a is halfway between two integers)
     */
    Field_Univariate_Derivative_1<T> remainder(const T& a) 
    {
        return Field_Univariate_Derivative_1<>(std::remainder(my_f0, a), my_f1);
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> remainder(const double& a) 
    {
        return Field_Univariate_Derivative_1<>(std::remainder(my_f0, a), my_f1);
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> remainder(const Field_Univariate_Derivative_1<T>& a) 
    {

        // compute k such that lhs % rhs = lhs - k rhs
        const T rem = std::remainder(my_f0, a.f0);
        const T k   = std::rint(my_f0.subtract(rem).divide(a.f0));

        return Field_Univariate_Derivative_1<>(rem, my_f1.subtract(k.multiply(a.f1)));

    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> negate() 
    {
        return Field_Univariate_Derivative_1<>(my_f0.negate(), my_f1.negate());
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> abs() 
    {
        if (Double.double_to_long_bits(my_f0.get_real()) < 0)
        {
            // we use the bits representation to also handle -0.0
            return negate();
        }
        return *this;
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> ceil() 
    {
        return Field_Univariate_Derivative_1<>(std::ceil(my_f0), my_f0.get_field().get_zero());
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> floor() 
    {
        return Field_Univariate_Derivative_1<>(std::floor(my_f0), my_f0.get_field().get_zero());
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> rint() 
    {
        return Field_Univariate_Derivative_1<>(std::rint(my_f0), my_f0.get_field().get_zero());
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> sign() 
    {
        return Field_Univariate_Derivative_1<>(FastMath.sign(my_f0), my_f0.get_field().get_zero());
    }

    /**
     * Returns the instance with the sign of the argument.
     * A NaN {@code sign} argument is treated as positive.
     *
     * @param sign the sign for the returned value
     * @return the instance with the same sign as the {@code sign} argument
     */
    Field_Univariate_Derivative_1<T> copy_sign(const T& sign) 
    {
        const long m = Double.double_to_long_bits(my_f0.get_real());
        const long s = Double.double_to_long_bits(sign.get_real());
        if ((m >= 0 && s >= 0) || (m < 0 && s < 0))
        {
            // Sign is currently OK
            return this;
        }
        return negate(); // flip sign
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> copy_sign(const Field_Univariate_Derivative_1<T>& sign) 
    {
        const long m = Double.double_to_long_bits(my_f0.get_real());
        const long s = Double.double_to_long_bits(sign.f0.get_real());
        if ((m >= 0 && s >= 0) || (m < 0 && s < 0))
        {
            // Sign is currently OK
            return this;
        }
        return negate(); // flip sign
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> copy_sign(const double& sign) 
    {
        const long m = Double.double_to_long_bits(my_f0.get_real());
        const long s = Double.double_to_long_bits(sign);
        if ((m >= 0 && s >= 0) || (m < 0 && s < 0))
        {
            // Sign is currently OK
            return *this;
        }
        return negate(); // flip sign
    }

    /** {@inherit_doc} */
    //override
    int get_exponent() 
    {
        return FastMath.get_exponent(my_f0.get_real());
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> scalb(const int& n) 
    {
        return Field_Univariate_Derivative_1<>(std::scalbn(my_f0, n), std::scalbn(my_f1, n));
    }

    /** {@inherit_doc}
     * <p>
     * The {@code ulp} function is a step function, hence all its derivatives are 0.
     * </p>
     * @since 2.0
     */
    //override
    Field_Univariate_Derivative_1<T> ulp() 
    {
        return Field_Univariate_Derivative_1<>(FastMath.ulp(my_f0), get_value_field().get_zero());
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> hypot(const Field_Univariate_Derivative_1<T>& y) 
    {

        if (std::isinf(my_f0.get_real()) || std::isinf(y.f0.get_real()))
        {
            return Field_Univariate_Derivative_1<>(my_f0.new_instance(INFINITY), my_f0.get_field().get_zero());
        }
        if (std::isnan(my_f0.get_real()) || std::isnan(y.f0.get_real()))
        {
            return Field_Univariate_Derivative_1<>(f0.new_instance(Double.NaN), my_f0.get_field().get_zero());
        }

        const int exp_x = get_exponent();
        const int exp_y = y.get_exponent();
        if (exp_x > exp_y + 27) 
        {
            // y is neglectible with respect to x
            return abs();
        }
        if (exp_y > exp_x + 27) 
        {
            // x is neglectible with respect to y
            return y.abs();
        }


        // find an intermediate scale to avoid both overflow and underflow
        const int middle_exp = (exp_x + exp_y) / 2;

        // scale parameters without losing precision
        const Field_Univariate_Derivative_1<T> scaled_x = scalb(-middle_exp);
        const Field_Univariate_Derivative_1<T> scaled_y = y.scalb(-middle_exp);

        // compute scaled hypotenuse
        const Field_Univariate_Derivative_1<T> scaled_h =
                scaled_x.multiply(scaled_x).add(scaled_y.multiply(scaled_y)).sqrt();

        // remove scaling
        return scaled_h.scalb(middle_exp);        
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> reciprocal() 
    {
        const T inv1 = my_f0.reciprocal();
        const T inv2 = inv1.multiply(inv1);
        return Field_Univariate_Derivative_1<>(inv1, f1.negate().multiply(inv2));
    }

    /** Compute composition of the instance by a function.
     * @param g0 value of the function at the current point (i.e. at {@code g(get_value())})
     * @param g1 first derivative of the function at the current point (i.e. at {@code g'(get_value())})
     * @return g(this)
     */
    Field_Univariate_Derivative_1<T> compose(const T& g0, const T& g1) 
    {
        return Field_Univariate_Derivative_1<>(g0, g1.multiply(my_f1));
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> sqrt() 
    {
        const T s = std::sqrt(my_f0);
        return compose(s, s.add(s).reciprocal());
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> cbrt() 
    {
        const T c = std::cbrt(my_f0);
        return compose(c, c.multiply(c).multiply(3).reciprocal());
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> root_n(const int& n) 
    {
        if (n == 2) 
        {
            return sqrt();
        }
        if (n == 3) 
        {
            return cbrt();
        }
        const T r = std::pow(my_f0, 1.0 / n);
        return compose(r, std::pow(r, n - 1).multiply(n).reciprocal());
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1Field<T> get_field() 
    {
        return Field_Univariate_Derivative_1Field.get_univariate_derivative1_field(my_f0.get_field());
    }

    /** Compute a<sup>x</sup> where a is a double and x a {@link Field_Univariate_Derivative_1}
     * @param a number to exponentiate
     * @param x power to apply
     * @param <T> the type of the function parameters and value
     * @return a<sup>x</sup>
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  Field_Univariate_Derivative_1<T> pow(const double& a, const Field_Univariate_Derivative_1<T>& x) 
    {
        if (a == 0) 
        {
            return x.get_field().get_zero();
        }
        const T a_x = std::pow(x.f0.new_instance(a), x.f0);
        return Field_Univariate_Derivative_1<>(a_x, a_x.multiply(std::log(a)).multiply(x.f1));
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> pow(const double& p) 
    {
        if (p == 0) 
        {
            return get_field().get_one();
        }
        const T f0_pm1 = std::pow(my_f0, p - 1);
        return compose(f0_pm1.multiply(my_f0), f0_pm1.multiply(p));
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> pow(const int& n) 
    {
        if (n == 0) 
        {
            return get_field().get_one();
        }
        const T f0_nm1 = std::pow(my_f0, n - 1);
        return compose(f0_nm1.multiply(my_f0), f0_nm1.multiply(n));
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> pow(const Field_Univariate_Derivative_1<T>& e) 
    {
        return log().multiply(e).exp();
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> exp() 
    {
        const T exp = std::exp(my_f0);
        return compose(exp, exp);
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> expm1() 
    {
        const T exp   = std::exp(my_f0);
        const T exp_m1 = std::expm1(my_f0);
        return compose(exp_m1, exp);
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> log() 
    {
        return compose(std::log(my_f0), my_f0.reciprocal());
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> log1p() 
    {
        return compose(std::log1p(my_f0), my_f0.add(1).reciprocal());
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> log10() 
    {
        return compose(std::log10(my_f0), my_f0.multiply(std::log(10.0)).reciprocal());
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> cos() 
    {
        const Field_Sin_Cos<T> sin_cos = Sin_Cos(my_f0);
        return compose(sin_cos.cos(), sin_cos.sin().negate());
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> sin() 
    {
        const Field_Sin_Cos<T> sin_cos = Sin_Cos(my_f0);
        return compose(sin_cos.sin(), sin_cos.cos());
    }

    /** {@inherit_doc} */
    //override
    Field_Sin_Cos<Field_Univariate_Derivative_1<T>> sin_cos() 
    {
        const Field_Sin_Cos<T> sin_cos = Sin_Cos(my_f0);
        return Field_Sin_Cos<>(Field_Univariate_Derivative_1<>(sin_cos.sin(), my_f1.multiply(sin_cos.cos())), Field_Univariate_Derivative_1<>(sin_cos.cos(), my_f1.multiply(sin_cos.sin()).negate()));
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> tan() 
    {
        const T tan = std::tan(my_f0);
        return compose(tan, tan.multiply(tan).add(1));
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> acos() 
    {
        return compose(std::acos(my_f0), my_f0.multiply(my_f0).negate().add(1).sqrt().reciprocal().negate());
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> asin() 
    {
        return compose(std::asin(my_f0), my_f0.multiply(my_f0).negate().add(1).sqrt().reciprocal());
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> atan() 
    {
        return compose(std::atan(my_f0), my_f0.multiply(my_f0).add(1).reciprocal());
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> atan2(const Field_Univariate_Derivative_1<T>& x) 
    {
        const T inv = my_f0.multiply(my_f0).add(x.f0.multiply(x.f0)).reciprocal();
        return Field_Univariate_Derivative_1<>(std::atan2(my_f0, x.f0), my_f0.linear_combination(x.f0, my_f1, x.f1.negate(), my_f0).multiply(inv));
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> cosh() 
    {
        return compose(std::cosh(my_f0), std::sinh(my_f0));
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> sinh() 
    {
        return compose(std::sinh(my_f0), std::cosh(my_f0));
    }

    /** {@inherit_doc} */
    //override
    Field_Sinh_Cosh<Field_Univariate_Derivative_1<T>> sinh_cosh() 
    {
        const Field_Sinh_Cosh<T> sinh_cosh = std::sinh_cosh(my_f0);
        return Field_Sinh_Cosh<>(new Field_Univariate_Derivative_1<>(sinh_cosh.sinh(), my_f1.multiply(sinh_cosh.cosh())), Field_Univariate_Derivative_1<>(sinh_cosh.cosh(), my_f1.multiply(sinh_cosh.sinh())));
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> tanh() 
    {
        const T tanh = std::tanh(my_f0);
        return compose(tanh, tanh.multiply(tanh).negate().add(1));
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> acosh() 
    {
        return compose(std::acosh(my_f0), my_f0.multiply(my_f0).subtract(1).sqrt().reciprocal());
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> asinh() 
    {
        return compose(std::asinh(my_f0), my_f0.multiply(my_f0).add(1).sqrt().reciprocal());
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> atanh() 
    {
        return compose(std::atanh(my_f0), my_f0.multiply(my_f0).negate().add(1).reciprocal());
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> to_degrees() 
    {
        return Field_Univariate_Derivative_1<>(FastMath.to_degrees(my_f0), FastMath.to_degrees(my_f1));
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> to_radians() 
    {
        return Field_Univariate_Derivative_1<>(FastMath.to_radians(my_f0), FastMath.to_radians(my_f1));
    }

    /** Evaluate Taylor expansion of a univariate derivative.
     * @param delta parameter offset \xce\x94x
     * @return value of the Taylor expansion at x + \xce\x94x
     */
    T taylor(const double delta) 
    {
        return my_f0.add(my_f1.multiply(delta));
    }

    /** Evaluate Taylor expansion of a univariate derivative.
     * @param delta parameter offset \xce\x94x
     * @return value of the Taylor expansion at x + \xce\x94x
     */
    T taylor(const T delta) 
    {
        return my_f0.add(f1.multiply(delta));
    }

    /**
     * Compute a linear combination.
     * @param a Factors.
     * @param b Factors.
     * @return <code>&Sigma;<sub>i</sub> a<sub>i</sub> b<sub>i</sub></code>.
     * @ if arrays dimensions don't match
     */
    Field_Univariate_Derivative_1<T> linear_combination(const std::vector<T>& a, const std::vector<Field_Univariate_Derivative_1<T>>& b) 
    {
        // extract values and first derivatives
        const Field<T> field = b[0].f0.get_field();
        const int      n  = b.size();
        std::vector<T> b0 = Math_Arrays::build_array(field, n);
        std::vector<T> b1 = Math_Arrays::build_array(field, n);
        for (int i{}; i < n; ++i) 
        {
            b0[i] = b[i].f0;
            b1[i] = b[i].f1;
        }

        return Field_Univariate_Derivative_1<>(b[0].f0.linear_combination(a, b0), b[0].f0.linear_combination(a, b1));

    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> linear_combination(const std::vector<Field_Univariate_Derivative_1<T>> a, const std::vector<Field_Univariate_Derivative_1<T>> b) 
    {

        // extract values and first derivatives
        const Field<T> field = a[0].f0.get_field();
        const int n  = a.size();
        const std::vector<T> a0 = Math_Arrays::build_array(field, n);
        const std::vector<T> b0 = Math_Arrays::build_array(field, n);
        const std::vector<T> a1 = Math_Arrays::build_array(field, 2 * n);
        const std::vector<T> b1 = Math_Arrays::build_array(field, 2 * n);
        for (int i{}; i < n; ++i) 
        {
            const Field_Univariate_Derivative_1<T> ai = a[i];
            const Field_Univariate_Derivative_1<T> bi = b[i];
            a0[i]         = ai.f0;
            b0[i]         = bi.f0;
            a1[2 * i]     = ai.f0;
            a1[2 * i + 1] = ai.f1;
            b1[2 * i]     = bi.f1;
            b1[2 * i + 1] = bi.f0;
        }

        return Field_Univariate_Derivative_1<>(a[0].f0.linear_combination(a0, b0), a[0].f0.linear_combination(a1, b1));

    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> linear_combination(const std::vector<double> a, const std::vector<Field_Univariate_Derivative_1<T>> b) 
    {

        // extract values and first derivatives
        const Field<T> field = b[0].f0.get_field();
        const int      n  = b.size();
        const std::vector<T> b0 = Math_Arrays::build_array(field, n);
        const std::vector<T> b1 = Math_Arrays::build_array(field, n);
        for (int i{}; i < n; ++i) 
        {
            b0[i] = b[i].f0;
            b1[i] = b[i].f1;
        }

        return Field_Univariate_Derivative_1<>(b[0].f0.linear_combination(a, b0), b[0].f0.linear_combination(a, b1));

    }

    /** {@inherit_doc} */
    //override
    public Field_Univariate_Derivative_1<T> linear_combination(const Field_Univariate_Derivative_1<T>& a1, const Field_Univariate_Derivative_1<T> b1, const Field_Univariate_Derivative_1<T> a2, const Field_Univariate_Derivative_1<T>& b2) 
    {
        return Field_Univariate_Derivative_1<>(a1.f0.linear_combination(a1.f0, b1.f0, a2.f0, b2.f0), a1.f0.linear_combination(a1.f0, b1.f1, a1.f1, b1.f0, a2.f0, b2.f1, a2.f1, b2.f0));
    }

    /** {@inherit_doc} */
    //override
    public Field_Univariate_Derivative_1<T> linear_combination(const double& a1, const Field_Univariate_Derivative_1<T> b1, const double& a2, const Field_Univariate_Derivative_1<T> b2) 
    {
        return Field_Univariate_Derivative_1<>(b1.f0.linear_combination(a1, b1.f0, a2, b2.f0), b1.f0.linear_combination(a1, b1.f1, a2, b2.f1));
    }

    /** {@inherit_doc} */
    //override
    public Field_Univariate_Derivative_1<T> linear_combination(const Field_Univariate_Derivative_1<T> a1, const Field_Univariate_Derivative_1<T> b1, const Field_Univariate_Derivative_1<T>& a2, const Field_Univariate_Derivative_1<T>& b2, const Field_Univariate_Derivative_1<T> a3, const Field_Univariate_Derivative_1<T> b3) 
    {
        const Field<T> field = a1.f0.get_field();
        const std::vector<T> a = Math_Arrays::build_array(field, 6);
        const std::vector<T> b = Math_Arrays::build_array(field, 6);
        a[0] = a1.f0;
        a[1] = a1.f1;
        a[2] = a2.f0;
        a[3] = a2.f1;
        a[4] = a3.f0;
        a[5] = a3.f1;
        b[0] = b1.f1;
        b[1] = b1.f0;
        b[2] = b2.f1;
        b[3] = b2.f0;
        b[4] = b3.f1;
        b[5] = b3.f0;
        return Field_Univariate_Derivative_1<>(a1.f0.linear_combination(a1.f0, b1.f0, a2.f0, b2.f0, a3.f0, b3.f0), a1.f0.linear_combination(a, b));
    }

    /**
     * Compute a linear combination.
     * @param a1 first factor of the first term
     * @param b1 second factor of the first term
     * @param a2 first factor of the second term
     * @param b2 second factor of the second term
     * @param a3 first factor of the third term
     * @param b3 second factor of the third term
     * @return a<sub>1</sub>&times;b<sub>1</sub> +
     * a<sub>2</sub>&times;b<sub>2</sub> + a<sub>3</sub>&times;b<sub>3</sub>
     * @see #linear_combination(double, Field_Univariate_Derivative_1, double, Field_Univariate_Derivative_1)
     * @see #linear_combination(double, Field_Univariate_Derivative_1, double, Field_Univariate_Derivative_1, double, Field_Univariate_Derivative_1, double, Field_Univariate_Derivative_1)
     * @exception  if number of free parameters or orders are inconsistent
     */
    Field_Univariate_Derivative_1<T> linear_combination(const T a1, const Field_Univariate_Derivative_1<T>& b1, const T a2, const Field_Univariate_Derivative_1<T>& b2, const T& a3, const Field_Univariate_Derivative_1<T>& b3) 
    {
        return Field_Univariate_Derivative_1<>(b1.f0.linear_combination(a1, b1.f0, a2, b2.f0, a3, b3.f0), b1.f0.linear_combination(a1, b1.f1, a2, b2.f1, a3, b3.f1));
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> linear_combination(const double& a1, const Field_Univariate_Derivative_1<T>& b1, const double& a2, const Field_Univariate_Derivative_1<T>& b2, const double& a3, const Field_Univariate_Derivative_1<T>& b3) 
    {
        return Field_Univariate_Derivative_1<>(b1.f0.linear_combination(a1, b1.f0, a2, b2.f0, a3, b3.f0), b1.f0.linear_combination(a1, b1.f1, a2, b2.f1, a3, b3.f1));
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> linear_combination(const Field_Univariate_Derivative_1<T>& a1, const Field_Univariate_Derivative_1<T>& b1, const Field_Univariate_Derivative_1<T>& a2, const Field_Univariate_Derivative_1<T> b2, const Field_Univariate_Derivative_1<T>& a3, const Field_Univariate_Derivative_1<T>& b3, const Field_Univariate_Derivative_1<T>& a4, const Field_Univariate_Derivative_1<T>& b4) 
    {
        const Field<T> field = a1.f0.get_field();
        const std::vector<T> a = Math_Arrays::build_array(field, 8);
        const std::vector<T> b = Math_Arrays::build_array(field, 8);
        a[0] = a1.f0;
        a[1] = a1.f1;
        a[2] = a2.f0;
        a[3] = a2.f1;
        a[4] = a3.f0;
        a[5] = a3.f1;
        a[6] = a4.f0;
        a[7] = a4.f1;
        b[0] = b1.f1;
        b[1] = b1.f0;
        b[2] = b2.f1;
        b[3] = b2.f0;
        b[4] = b3.f1;
        b[5] = b3.f0;
        b[6] = b4.f1;
        b[7] = b4.f0;
        return Field_Univariate_Derivative_1<>(a1.f0.linear_combination(a1.f0, b1.f0, a2.f0, b2.f0, a3.f0, b3.f0, a4.f0, b4.f0), a1.f0.linear_combination(a, b));
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> linear_combination(const double& a1, const Field_Univariate_Derivative_1<T>& b1, const double& a2, const Field_Univariate_Derivative_1<T>& b2, const double& a3, const Field_Univariate_Derivative_1<T>& b3, const double& a4, const Field_Univariate_Derivative_1<T>& b4) 
    {
        return Field_Univariate_Derivative_1<>(b1.f0.linear_combination(a1, b1.f0, a2, b2.f0, a3, b3.f0, a4, b4.f0), b1.f0.linear_combination(a1, b1.f1, a2, b2.f1, a3, b3.f1, a4, b4.f1));
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_1<T> get_pi() 
    {
        const T zero = get_value_field().get_zero();
        return Field_Univariate_Derivative_1<>(zero.get_pi(), zero);
    }

    /** Test for the equality of two univariate derivatives.
     * <p>
     * univariate derivatives are considered equal if they have the same derivatives.
     * </p>
     * @param other Object to test for equality to this
     * @return true if two univariate derivatives are equal
     */
    //override
    bool equals(Object other) 
    {

        if (this == other) 
        {
            return true;
        }

        if (dynamic_cast<const Field_Univariate_Derivative_1*>(*other) != nullptr)
        {
            //@Suppress_Warnings("unchecked")
            const Field_Univariate_Derivative_1<T> rhs = (Field_Univariate_Derivative_1<T>) other;
            return f0.equals(rhs.f0) && f1.equals(rhs.f1);
        }

        return false;

    }

    /** Get a hash_code for the univariate derivative.
     * @return a hash code value for this object
     */
    //override
    int hash_code() 
    {
        return 453 - 19 * f0.hash_code() + 37 * f1.hash_code();
    }

}


