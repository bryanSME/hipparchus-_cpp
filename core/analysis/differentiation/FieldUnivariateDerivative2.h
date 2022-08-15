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
 * and {@link Field_Derivative_Structure#get_order() derivation order} limited to two.
 * It should have less overhead than {@link Field_Derivative_Structure} in its domain.</p>
 * <p>This class is an implementation of Rall's numbers. Rall's numbers are an
 * extension to the real numbers used throughout mathematical expressions; they hold
 * the derivative together with the value of a function.</p>
 * <p>{@link Field_Univariate_Derivative_2} instances can be used directly thanks to
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
 * @see Field_Univariate_Derivative_1
 * @see Field_Gradient
 * @since 1.7
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
class Field_Univariate_Derivative_2 : public Field_Univariate_Derivative<T, Field_Univariate_Derivative_2<T>> 
{
private:
    /** Value of the function. */
    const T my_f0;

    /** First derivative of the function. */
    const T my_f1;

    /** Second derivative of the function. */
    const T my_f2;

public:
    /** Build an instance with values and derivative.
     * @param my_f0 value of the function
     * @param my_f1 first derivative of the function
     * @param my_f2 second derivative of the function
     */
    Field_Univariate_Derivative_2(const T& my_f0, const T& my_f1, const T& my_f2) : my_f0{ my_f0 }, my_f1{ m1 }, my_f2{ m2 } {};

    /** Build an instance from a {@link Derivative_Structure}.
     * @param ds derivative structure
     * @exception  if either {@code ds} parameters
     * is not 1 or {@code ds} order is not 2
     */
    Field_Univariate_Derivative_2(const Field_Derivative_Structure<T>& ds)  
    {
        Math_Utils::check_dimension(ds.get_free_parameters(), 1);
        Math_Utils::check_dimension(ds.get_order(), 2);
       my_f0 = ds.get_value();
        my_f1 = ds.get_partial_derivative(1);
        my_f2 = ds.get_partial_derivative(2);
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> new_instance(const double& value) 
    {
        const T zero = my_f0.get_field().get_zero();
        return Field_Univariate_Derivative_2<>(zero.new_instance(value), zero, zero);
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
    T get_derivative(const int& n) const
    {
        switch (n) 
        {
            case 0 :
                return my_f0;
            case 1 :
                return my_f1;
            case 2 :
                return my_f2;
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
        return 2;
    }

    /** Get the first derivative.
     * @return first derivative
     * @see #get_value()
     */
    T get_first_derivative() const
    {
        return my_f1;
    }

    /** Get the second derivative.
     * @return second derivative
     * @see #get_value()
     * @see #get_first_derivative()
     */
    T get_second_derivative() const
    {
        return my_f2;
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
        return get_field().get_conversion_factory().build(my_f0, my_f1, my_f2);
    }

    /** '+' operator.
     * @param a right hand side parameter of the operator
     * @return this+a
     */
    Field_Univariate_Derivative_2<T> add(const T& a) 
    {
        return Field_Univariate_Derivative_2<>(my_f0.add(a), my_f1, my_f2);
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> add(const double& a) 
    {
        return Field_Univariate_Derivative_2<>(my_f0.add(a), my_f1, my_f2);
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> add(const Field_Univariate_Derivative_2<T>& a) 
    {
        return Field_Univariate_Derivative_2<>(my_f0.add(a.f0), my_f1.add(a.f1), my_f2.add(a.f2));
    }

    /** '-' operator.
     * @param a right hand side parameter of the operator
     * @return this-a
     */
    Field_Univariate_Derivative_2<T> subtract(const T& a) 
    {
        return Field_Univariate_Derivative_2<>(my_f0.subtract(a), my_f1, my_f2);
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> subtract(const double& a) 
    {
        return Field_Univariate_Derivative_2<>(my_f0.subtract(a), my_f1, my_f2);
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> subtract(const Field_Univariate_Derivative_2<T>& a) 
    {
        return Field_Univariate_Derivative_2<>(my_f0.subtract(a.f0), my_f1.subtract(a.f1), my_f2.subtract(a.f2));
    }

    /** '&times;' operator.
     * @param a right hand side parameter of the operator
     * @return this&times;a
     */
    Field_Univariate_Derivative_2<T> multiply(const T& a) 
    {
        return Field_Univariate_Derivative_2<>(my_f0.multiply(a), my_f1.multiply(a), my_f2.multiply(a));
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> multiply(const int& n) 
    {
        return Field_Univariate_Derivative_2<>(my_f0.multiply(n), my_f1.multiply(n), my_f2.multiply(n));
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> multiply(const double& a) 
    {
        return Field_Univariate_Derivative_2<>(my_f0.multiply(a), my_f1.multiply(a), my_f2.multiply(a));
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> multiply(const Field_Univariate_Derivative_2<T>& a) 
    {
        return Field_Univariate_Derivative_2<>(my_f0.multiply(a.f0), a.f0.linear_combination(my_f1, a.f0, my_f0, a.f1), a.f0.linear_combination(f2, a.f0, my_f1.add(my_f1), a.f1, my_f0, a.f2));
    }

    /** '&divide;' operator.
     * @param a right hand side parameter of the operator
     * @return this&divide;a
     */
    Field_Univariate_Derivative_2<T> divide(const T& a) 
    {
        const T inv1 = a.reciprocal();
        return Field_Univariate_Derivative_2<>(my_f0.multiply(inv1), my_f1.multiply(inv1), my_f2.multiply(inv1));
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> divide(const double& a) 
    {
        const double inv1 = 1.0 / a;
        return Field_Univariate_Derivative_2<>(my_f0.multiply(inv1), my_f1.multiply(inv1), my_f2.multiply(inv1));
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> divide(const Field_Univariate_Derivative_2<T>& a) 
    {
        const T inv1 = a.f0.reciprocal();
        const T inv2 = inv1.multiply(inv1);
        const T inv3 = inv1.multiply(inv2);
        return Field_Univariate_Derivative_2<>(my_f0.multiply(inv1), a.f0.linear_combination(my_f1, a.f0, my_f0.negate(), a.f1).multiply(inv2), a.f0.linear_combination(f2, a.f0.multiply(a.f0), my_f1.multiply(-2), a.f0.multiply(a.f1), my_f0.add(my_f0), a.f1.multiply(a.f1), my_f0.negate(), a.f0.multiply(a.f2)).multiply(inv3));
    }

    /** IEEE remainder operator.
     * @param a right hand side parameter of the operator
     * @return this - n &times; a where n is the closest integer to this/a
     * (the even integer is chosen for n if this/a is halfway between two integers)
     */
    Field_Univariate_Derivative_2<T> remainder(const T& a) 
    {
        return Field_Univariate_Derivative_2<>(std::remainder(my_f0, a), my_f1, my_f2);
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> remainder(const double& a) 
    {
        return Field_Univariate_Derivative_2<>(std::remainder(my_f0, a), my_f1, my_f2);
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> remainder(const Field_Univariate_Derivative_2<T>& a) 
    {

        // compute k such that lhs % rhs = lhs - k rhs
        const T rem = std::remainder(my_f0, a.f0);
        const T k   = std::rint(my_f0.subtract(rem).divide(a.f0));

        return Field_Univariate_Derivative_2<>(rem, my_f1.subtract(k.multiply(a.f1)), my_f2.subtract(k.multiply(a.f2)));

    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> negate() 
    {
        return Field_Univariate_Derivative_2<>(my_f0.negate(), my_f1.negate(), my_f2.negate());
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> abs() 
    {
        if (Double.double_to_long_bits(my_f0.get_real()) < 0) 
        {
            // we use the bits representation to also handle -0.0
            return negate();
        }
        return this;
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> ceil() 
    {
        const T zero = my_f0.get_field().get_zero();
        return Field_Univariate_Derivative_2<>(std::ceil(my_f0), zero, zero);
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> floor() 
    {
        const T zero = my_f0.get_field().get_zero();
        return Field_Univariate_Derivative_2<>(std::floor(my_f0), zero, zero);
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> rint() 
    {
        const T zero = my_f0.get_field().get_zero();
        return Field_Univariate_Derivative_2<>(std::rint(my_f0), zero, zero);
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> sign() 
    {
        const T zero = my_f0.get_field().get_zero();
        return Field_Univariate_Derivative_2<>(FastMath.sign(my_f0), zero, zero);
    }

    /**
     * Returns the instance with the sign of the argument.
     * A NaN {@code sign} argument is treated as positive.
     *
     * @param sign the sign for the returned value
     * @return the instance with the same sign as the {@code sign} argument
     */
    Field_Univariate_Derivative_2<T> copy_sign(const T& sign) 
    {
        long m = Double.double_to_long_bits(my_f0.get_real());
        long s = Double.double_to_long_bits(sign.get_real());
        if ((m >= 0 && s >= 0) || (m < 0 && s < 0)) { // Sign is currently OK
            return this;
        }
        return negate(); // flip sign
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> copy_sign(const Field_Univariate_Derivative_2<T>& sign) 
    {
        long m = Double.double_to_long_bits(my_f0.get_real());
        long s = Double.double_to_long_bits(sign.f0.get_real());
        if ((m >= 0 && s >= 0) || (m < 0 && s < 0)) { // Sign is currently OK
            return this;
        }
        return negate(); // flip sign
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> copy_sign(const double sign) 
    {
        long m = Double.double_to_long_bits(my_f0.get_real());
        long s = Double.double_to_long_bits(sign);
        if ((m >= 0 && s >= 0) || (m < 0 && s < 0)) { // Sign is currently OK
            return this;
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
    Field_Univariate_Derivative_2<T> scalb(const int& n) 
    {
        return Field_Univariate_Derivative_2<>(std::scalbn(my_f0, n), std::scalbn(my_f1, n), std::scalbn(f2, n));
    }

    /** {@inherit_doc}
     * <p>
     * The {@code ulp} function is a step function, hence all its derivatives are 0.
     * </p>
     * @since 2.0
     */
    //override
    Field_Univariate_Derivative_2<T> ulp() 
    {
        const T zero = get_value_field().get_zero();
        return Field_Univariate_Derivative_2<>(FastMath.ulp(my_f0), zero, zero);
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> hypot(const Field_Univariate_Derivative_2<T>& y) 
    {

        if (std::isinf(my_f0.get_real()) || std::isinf(y.f0.get_real())) 
        {
            const T zero = my_f0.get_field().get_zero();
            return Field_Univariate_Derivative_2<>(my_f0.new_instance(INFINITY), zero, zero);
        }
        if (std::isnan(my_f0.get_real()) || std::isnan(y.f0.get_real())) 
        {
            const T zero = my_f0.get_field().get_zero();
            return Field_Univariate_Derivative_2<>(my_f0.new_instance(Double.NaN), zero, zero);
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
        const Field_Univariate_Derivative_2<T>& scaled_x = scalb(-middle_exp);
        const Field_Univariate_Derivative_2<T>& scaled_y = y.scalb(-middle_exp);

        // compute scaled hypotenuse
        const Field_Univariate_Derivative_2<T>& scaled_h = scaled_x.multiply(scaled_x).add(scaled_y.multiply(scaled_y)).sqrt();

        // remove scaling
        return scaled_h.scalb(middle_exp);
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> reciprocal() 
    {
        const T inv1 = my_f0.reciprocal();
        const T inv2 = inv1.multiply(inv1);
        const T inv3 = inv1.multiply(inv2);
        return Field_Univariate_Derivative_2<>(inv1, my_f1.negate().multiply(inv2), my_f0.linear_combination(my_f1.add(my_f1), my_f1, my_f0.negate(), my_f2).multiply(inv3));
    }

    /** Compute composition of the instance by a function.
     * @param g0 value of the function at the current point (i.e. at {@code g(get_value())})
     * @param g1 first derivative of the function at the current point (i.e. at {@code g'(get_value())})
     * @param g2 second derivative of the function at the current point (i.e. at {@code g''(get_value())})
     * @return g(this)
     */
    Field_Univariate_Derivative_2<T> compose(const T& g0, const T& g1, const T& g2) 
    {
        return Field_Univariate_Derivative_2<>(
            g0, g1.multiply(my_f1),
            my_f0.linear_combination(g1, my_f2, g2, my_f1.multiply(my_f1))
        );
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> sqrt() 
    {
        const T s = std::sqrt(my_f0);
        return compose(s, s.add(s).reciprocal(), s.multiply(-4).multiply(my_f0).reciprocal());
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> cbrt() 
    {
        const T c  = std::cbrt(my_f0);
        const T c2 = c.multiply(c);
        return compose(c, c2.multiply(3).reciprocal(), c2.multiply(4.5).multiply(my_f0).reciprocal());
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> root_n(const int& n) 
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
        const T z = std::pow(r, n - 1).multiply(n);
        return compose(r, z.reciprocal(), z.multiply(z).multiply(r).reciprocal().multiply(1 -n));
       
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2_Field<T> get_field() 
    {
        return Field_Univariate_Derivative_2_Field.get_univariate__derivative_2__field(my_f0.get_field());
    }

    /** Compute a<sup>x</sup> where a is a double and x a {@link Field_Univariate_Derivative_2}
     * @param a number to exponentiate
     * @param x power to apply
     * @param <T> the type of the function parameters and value
     * @return a<sup>x</sup>
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  Field_Univariate_Derivative_2<T> pow(const double& a, const Field_Univariate_Derivative_2<T>& x) 
    {
        if (a == 0) 
        {
            return x.get_field().get_zero();
        }

        const T      a_x    = std::pow(x.f0.new_instance(a), x.f0);
        const double ln_a   = std::log(a);
        const T      a_xln_a = a_x.multiply(ln_a);
        return Field_Univariate_Derivative_2<>(a_x, a_xln_a.multiply(x.f1), a_xln_a.multiply(x.f1.multiply(x.f1).multiply(ln_a).add(x.f2)));
        
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> pow(const double& p) 
    {
        if (p == 0) 
        {
            return get_field().get_one();
        }

        const T my_f0_pm2 = std::pow(my_f0, p - 2);
        const T my_f0_pm1 = my_f0_pm2.multiply(my_f0);
        const T my_f0P   = my_f0_pm1.multiply(my_f0);
        return compose(my_f0P, my_f0_pm1.multiply(p), my_f0_pm2.multiply(p * (p - 1)));
        
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> pow(const int& n) 
    {
        if (n == 0) 
        {
            return get_field().get_one();
        }

        const T my_f0_nm2 = std::pow(my_f0, n - 2);
        const T my_f0_nm1 = my_f0_nm2.multiply(my_f0);
        const T my_f0N   = my_f0_nm1.multiply(my_f0);
        return compose(my_f0N, my_f0_nm1.multiply(n), my_f0_nm2.multiply(n * (n - 1)));
        
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> pow(const Field_Univariate_Derivative_2<T>& e) 
    {
        return log().multiply(e).exp();
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> exp() 
    {
        const T exp = std::exp(my_f0);
        return compose(exp, exp, exp);
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> expm1() 
    {
        const T exp   = std::exp(my_f0);
        const T exp_m1 = std::expm1(my_f0);
        return compose(exp_m1, exp, exp);
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> log() 
    {
        const T inv = my_f0.reciprocal();
        return compose(std::log(my_f0), inv, inv.multiply(inv).negate());
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> log1p() 
    {
        const T inv = my_f0.add(1).reciprocal();
        return compose(std::log1p(my_f0), inv, inv.multiply(inv).negate());
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> log10() 
    {
        const T inv_f0 = my_f0.reciprocal();
        const T inv = inv_f0.divide(std::log(10.0));
        return compose(std::log10(my_f0), inv, inv.multiply(inv_f0).negate());
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> cos() 
    {
        const Field_Sin_Cos<T> sin_cos = Sin_Cos(my_f0);
        return compose(sin_cos.cos(), sin_cos.sin().negate(), sin_cos.cos().negate());
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> sin() 
    {
        const Field_Sin_Cos<T> sin_cos = Sin_Cos(my_f0);
        return compose(sin_cos.sin(), sin_cos.cos(), sin_cos.sin().negate());
    }

    /** {@inherit_doc} */
    //override
    Field_Sin_Cos<Field_Univariate_Derivative_2<T>> sin_cos() 
    {
        const Field_Sin_Cos<T> sin_cos = Sin_Cos(my_f0);
        const T m_sin = sin_cos.sin().negate();
        const T m_cos = sin_cos.cos().negate();
        return Field_Sin_Cos<>(compose(sin_cos.sin(), sin_cos.cos(), m_sin), compose(sin_cos.cos(), m_sin, m_cos));
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> tan() 
    {
        const T tan  = std::tan(my_f0);
        const T sec2 = tan.multiply(tan).add(1);
        return compose(tan, sec2, sec2.add(sec2).multiply(tan));
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> acos() 
    {
        const T inv = my_f0.multiply(my_f0).negate().add(1).reciprocal();
        const T mS  = inv.sqrt().negate();
        return compose(std::acos(my_f0), mS, mS.multiply(my_f0).multiply(inv));
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> asin() 
    {
        const T inv = my_f0.multiply(my_f0).negate().add(1).reciprocal();
        const T s   = inv.sqrt();
        return compose(std::asin(my_f0), s, s.multiply(my_f0).multiply(inv));
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> atan() 
    {
        const T inv = my_f0.multiply(my_f0).add(1).reciprocal();
        return compose(std::atan(my_f0), inv, my_f0.multiply(-2).multiply(inv).multiply(inv));
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> atan2(const Field_Univariate_Derivative_2<T>& x) 
    {
        const T x2    = x.f0.multiply(x.f0);
        const T my_f02   = my_f0.add(my_f0);
        const T inv   = my_f0.multiply(my_f0).add(x2).reciprocal();
        const T atan0 = std::atan2(my_f0, x.f0);
        const T atan1 = my_f0.linear_combination(x.f0, my_f1, x.f1.negate(), my_f0).multiply(inv);
        const T c     = my_f0.linear_combination(f2, x2, my_f1.multiply(-2), x.f0.multiply(x.f1), my_f02, x.f1.multiply(x.f1), my_f0.negate(), x.f0.multiply(x.f2)).multiply(inv);
        return Field_Univariate_Derivative_2<>(atan0, atan1, c.subtract(my_f02.multiply(atan1).multiply(atan1)).divide(x.f0));
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> cosh() 
    {
        const T c = std::cosh(my_f0);
        const T s = std::sinh(my_f0);
        return compose(c, s, c);
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> sinh() 
    {
        const T c = std::cosh(my_f0);
        const T s = std::sinh(my_f0);
        return compose(s, c, s);
    }

    /** {@inherit_doc} */
    //override
    Field_Sinh_Cosh<Field_Univariate_Derivative_2<T>> sinh_cosh() 
    {
        const Field_Sinh_Cosh<T> sinh_cosh = std::sinh_cosh(my_f0);
        return Field_Sinh_Cosh<>(compose(sinh_cosh.sinh(), sinh_cosh.cosh(), sinh_cosh.sinh()), compose(sinh_cosh.cosh(), sinh_cosh.sinh(), sinh_cosh.cosh()));
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> tanh() 
    {
        const T tanh  = std::tanh(my_f0);
        const T sech2 = tanh.multiply(tanh).negate().add(1);
        return compose(tanh, sech2, sech2.multiply(-2).multiply(tanh));
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> acosh() 
    {
        const T inv = my_f0.multiply(my_f0).subtract(1).reciprocal();
        const T s   = inv.sqrt();
        return compose(std::acosh(my_f0), s, my_f0.negate().multiply(s).multiply(inv));
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> asinh() 
    {
        const T inv = my_f0.multiply(my_f0).add(1).reciprocal();
        const T s   = inv.sqrt();
        return compose(std::asinh(my_f0), s, my_f0.negate().multiply(s).multiply(inv));
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> atanh() 
    {
        const T inv = my_f0.multiply(my_f0).negate().add(1).reciprocal();
        return compose(std::atanh(my_f0), inv, my_f0.add(my_f0).multiply(inv).multiply(inv));
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> to_degrees() 
    {
        return Field_Univariate_Derivative_2<>(FastMath.to_degrees(my_f0), FastMath.to_degrees(my_f1), FastMath.to_degrees(f2));
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> to_radians() 
    {
        return Field_Univariate_Derivative_2<>(FastMath.to_radians(my_f0), FastMath.to_radians(my_f1), FastMath.to_radians(f2));
    }

    /** Evaluate Taylor expansion a univariate derivative.
     * @param delta parameter offset \xce\x94x
     * @return value of the Taylor expansion at x + \xce\x94x
     */
    T taylor(const double delta) 
    {
        return my_f0.add(my_f1.add(f2.multiply(0.5 * delta)).multiply(delta));
    }

    /** Evaluate Taylor expansion a univariate derivative.
     * @param delta parameter offset \xce\x94x
     * @return value of the Taylor expansion at x + \xce\x94x
     */
    T taylor(const T& delta) 
    {
        return my_f0.add(my_f1.add(f2.multiply(delta.multiply(0.5))).multiply(delta));
    }

    /**
     * Compute a linear combination.
     * @param a Factors.
     * @param b Factors.
     * @return <code>&Sigma;<sub>i</sub> a<sub>i</sub> b<sub>i</sub></code>.
     * @ if arrays dimensions don't match
     */
    Field_Univariate_Derivative_2<T> linear_combination(const std::vector<T>& a, const std::vector<Field_Univariate_Derivative_2<T>&b) 
    {
        // extract values and derivatives
        const Field<T> field = b[0].f0.get_field();
        const auto n  = b.size();
        const std::vector<T> b0 = Math_Arrays::build_array(field, n);
        const std::vector<T> b1 = Math_Arrays::build_array(field, n);
        const std::vector<T> b2 = Math_Arrays::build_array(field, n);
        for (int i{}; i < n; ++i) 
        {
            b0[i] = b[i].f0;
            b1[i] = b[i].f1;
            b2[i] = b[i].f2;
        }

        return Field_Univariate_Derivative_2<>(b[0].f0.linear_combination(a, b0), b[0].f0.linear_combination(a, b1), b[0].f0.linear_combination(a, b2));

    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> linear_combination(const std::vector<Field_Univariate_Derivative_2<T>& a, const std::vector<Field_Univariate_Derivative_2<T>& b) 
    {
        // extract values and derivatives
        const Field<T> field = a[0].f0.get_field();
        const int& n  = a.size();
        const std::vector<T> a0 = Math_Arrays::build_array(field, n);
        const std::vector<T> b0 = Math_Arrays::build_array(field, n);
        const std::vector<T> a1 = Math_Arrays::build_array(field, 2 * n);
        const std::vector<T> b1 = Math_Arrays::build_array(field, 2 * n);
        std::vector<T> a2 = Math_Arrays::build_array(field, 3 * n);
        std::vector<T> b2 = Math_Arrays::build_array(field, 3 * n);
        for (int i{}; i < n; ++i) 
        {
            Field_Univariate_Derivative_2<T> ai = a[i];
            Field_Univariate_Derivative_2<T> bi = b[i];
            a0[i]         = ai.f0;
            b0[i]         = bi.f0;
            a1[2 * i]     = ai.f0;
            a1[2 * i + 1] = ai.f1;
            b1[2 * i]     = bi.f1;
            b1[2 * i + 1] = bi.f0;
            a2[3 * i]     = ai.f0;
            a2[3 * i + 1] = ai.f1.add(ai.f1);
            a2[3 * i + 2] = ai.f2;
            b2[3 * i]     = bi.f2;
            b2[3 * i + 1] = bi.f1;
            b2[3 * i + 2] = bi.f0;
        }

        return Field_Univariate_Derivative_2<>(a[0].f0.linear_combination(a0, b0), a[0].f0.linear_combination(a1, b1), a[0].f0.linear_combination(a2, b2));

    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> linear_combination(const std::vector<double>& a, const std::vector<Field_Univariate_Derivative_2<T>& b) 
    {

        // extract values and derivatives
        const Field<T> field = b[0].f0.get_field();
        const int      n  = b.size();
        const std::vector<T> b0 = Math_Arrays::build_array(field, n);
        const std::vector<T> b1 = Math_Arrays::build_array(field, n);
        const std::vector<T> b2 = Math_Arrays::build_array(field, n);
        for (int i{}; i < n; ++i) 
        {
            b0[i] = b[i].f0;
            b1[i] = b[i].f1;
            b2[i] = b[i].f2;
        }

        return Field_Univariate_Derivative_2<>(b[0].f0.linear_combination(a, b0), b[0].f0.linear_combination(a, b1), b[0].f0.linear_combination(a, b2));

    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> linear_combination(const Field_Univariate_Derivative_2<T>& a1, const Field_Univariate_Derivative_2<T>& b1, const Field_Univariate_Derivative_2<T>& a2, const Field_Univariate_Derivative_2<T>& b2) 
    {
        const Field<T> field = a1.f0.get_field();
        const std::vector<T>      u2    = Math_Arrays::build_array(field, 6);
        const std::vector<T>      v2    = Math_Arrays::build_array(field, 6);
        u2[0] = a1.f0;
        u2[1] = a1.f1.add(a1.f1);
        u2[2] = a1.f2;
        u2[3] = a2.f0;
        u2[4] = a2.f1.add(a2.f1);
        u2[5] = a2.f2;
        v2[0] = b1.f2;
        v2[1] = b1.f1;
        v2[2] = b1.f0;
        v2[3] = b2.f2;
        v2[4] = b2.f1;
        v2[5] = b2.f0;
        return Field_Univariate_Derivative_2<>(a1.f0.linear_combination(a1.f0, b1.f0, a2.f0, b2.f0), a1.f0.linear_combination(a1.f0, b1.f1, a1.f1, b1.f0, a2.f0, b2.f1, a2.f1, b2.f0), a1.f0.linear_combination(u2, v2));
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> linear_combination(const double& a1, const Field_Univariate_Derivative_2<T>& b1, const double& a2, const Field_Univariate_Derivative_2<T>& b2) 
    {
        return Field_Univariate_Derivative_2<>(b1.f0.linear_combination(a1, b1.f0, a2, b2.f0), b1.f0.linear_combination(a1, b1.f1, a2, b2.f1), b1.f0.linear_combination(a1, b1.f2, a2, b2.f2));
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> linear_combination(const Field_Univariate_Derivative_2<T>& a1, const Field_Univariate_Derivative_2<T>& b1, const Field_Univariate_Derivative_2<T>& a2, const Field_Univariate_Derivative_2<T>& b2, const Field_Univariate_Derivative_2<T>& a3, const Field_Univariate_Derivative_2<T>& b3) 
    {
        const Field<T> field = a1.f0.get_field();
        const std::vector<T>      u1     = Math_Arrays::build_array(field, 6);
        const std::vector<T>      v1     = Math_Arrays::build_array(field, 6);
        u1[0] = a1.f0;
        u1[1] = a1.f1;
        u1[2] = a2.f0;
        u1[3] = a2.f1;
        u1[4] = a3.f0;
        u1[5] = a3.f1;
        v1[0] = b1.f1;
        v1[1] = b1.f0;
        v1[2] = b2.f1;
        v1[3] = b2.f0;
        v1[4] = b3.f1;
        v1[5] = b3.f0;
        const std::vector<T>      u2     = Math_Arrays::build_array(field, 9);
        const std::vector<T>      v2     = Math_Arrays::build_array(field, 9);
        u2[0] = a1.f0;
        u2[1] = a1.f1.add(a1.f1);
        u2[2] = a1.f2;
        u2[3] = a2.f0;
        u2[4] = a2.f1.add(a2.f1);
        u2[5] = a2.f2;
        u2[6] = a3.f0;
        u2[7] = a3.f1.add(a3.f1);
        u2[8] = a3.f2;
        v2[0] = b1.f2;
        v2[1] = b1.f1;
        v2[2] = b1.f0;
        v2[3] = b2.f2;
        v2[4] = b2.f1;
        v2[5] = b2.f0;
        v2[6] = b3.f2;
        v2[7] = b3.f1;
        v2[8] = b3.f0;
        return Field_Univariate_Derivative_2<>(a1.f0.linear_combination(a1.f0, b1.f0, a2.f0, b2.f0, a3.f0, b3.f0), a1.f0.linear_combination(u1, v1), a1.f0.linear_combination(u2, v2));
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
     * @see #linear_combination(double, Field_Univariate_Derivative_2, double, Field_Univariate_Derivative_2)
     * @see #linear_combination(double, Field_Univariate_Derivative_2, double, Field_Univariate_Derivative_2, double, Field_Univariate_Derivative_2, double, Field_Univariate_Derivative_2)
     * @exception  if number of free parameters or orders are inconsistent
     */
    Field_Univariate_Derivative_2<T> linear_combination(const T& a1, const Field_Univariate_Derivative_2<T>& b1, const T a2, const Field_Univariate_Derivative_2<T>& b2, const T a3, const Field_Univariate_Derivative_2<T>& b3) 
    {
        return Field_Univariate_Derivative_2<>(b1.f0.linear_combination(a1, b1.f0, a2, b2.f0, a3, b3.f0), b1.f0.linear_combination(a1, b1.f1, a2, b2.f1, a3, b3.f1), b1.f0.linear_combination(a1, b1.f2, a2, b2.f2, a3, b3.f2));
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> linear_combination(const double& a1, const Field_Univariate_Derivative_2<T>& b1, const double& a2, const Field_Univariate_Derivative_2<T>& b2, const double& a3, const Field_Univariate_Derivative_2<T>& b3) 
    {
        return Field_Univariate_Derivative_2<>(b1.f0.linear_combination(a1, b1.f0, a2, b2.f0, a3, b3.f0), b1.f0.linear_combination(a1, b1.f1, a2, b2.f1, a3, b3.f1), b1.f0.linear_combination(a1, b1.f2, a2, b2.f2, a3, b3.f2));
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> linear_combination(const Field_Univariate_Derivative_2<T>& a1, const Field_Univariate_Derivative_2<T>& b1, const Field_Univariate_Derivative_2<T>& a2, const Field_Univariate_Derivative_2<T>& b2, const Field_Univariate_Derivative_2<T>& a3, const Field_Univariate_Derivative_2<T>& b3, const Field_Univariate_Derivative_2<T>& a4, const Field_Univariate_Derivative_2<T>& b4) 
    {
        const Field<T> field = a1.f0.get_field();
        const std::vector<T> u1 = Math_Arrays::build_array(field, 8);
        const std::vector<T> v1 = Math_Arrays::build_array(field, 8);
        u1[0] = a1.f0;
        u1[1] = a1.f1;
        u1[2] = a2.f0;
        u1[3] = a2.f1;
        u1[4] = a3.f0;
        u1[5] = a3.f1;
        u1[6] = a4.f0;
        u1[7] = a4.f1;
        v1[0] = b1.f1;
        v1[1] = b1.f0;
        v1[2] = b2.f1;
        v1[3] = b2.f0;
        v1[4] = b3.f1;
        v1[5] = b3.f0;
        v1[6] = b4.f1;
        v1[7] = b4.f0;
        const std::vector<T> u2 = Math_Arrays::build_array(field, 12);
        const std::vector<T> v2 = Math_Arrays::build_array(field, 12);
        u2[ 0] = a1.f0;
        u2[ 1] = a1.f1.add(a1.f1);
        u2[ 2] = a1.f2;
        u2[ 3] = a2.f0;
        u2[ 4] = a2.f1.add(a2.f1);
        u2[ 5] = a2.f2;
        u2[ 6] = a3.f0;
        u2[ 7] = a3.f1.add(a3.f1);
        u2[ 8] = a3.f2;
        u2[ 9] = a4.f0;
        u2[10] = a4.f1.add(a4.f1);
        u2[11] = a4.f2;
        v2[ 0] = b1.f2;
        v2[ 1] = b1.f1;
        v2[ 2] = b1.f0;
        v2[ 3] = b2.f2;
        v2[ 4] = b2.f1;
        v2[ 5] = b2.f0;
        v2[ 6] = b3.f2;
        v2[ 7] = b3.f1;
        v2[ 8] = b3.f0;
        v2[ 9] = b4.f2;
        v2[10] = b4.f1;
        v2[11] = b4.f0;
        return Field_Univariate_Derivative_2<>(a1.f0.linear_combination(a1.f0, b1.f0, a2.f0, b2.f0, a3.f0, b3.f0, a4.f0, b4.f0), a1.f0.linear_combination(u1, v1), a1.f0.linear_combination(u2, v2));
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> linear_combination(const double& a1, const Field_Univariate_Derivative_2<T>& b1, const double& a2, const Field_Univariate_Derivative_2<T>& b2, const double& a3, const Field_Univariate_Derivative_2<T>& b3, const double& a4, const Field_Univariate_Derivative_2<T>& b4) 
    {
        return Field_Univariate_Derivative_2<>(b1.f0.linear_combination(a1, b1.f0, a2, b2.f0, a3, b3.f0, a4, b4.f0), b1.f0.linear_combination(a1, b1.f1, a2, b2.f1, a3, b3.f1, a4, b4.f1), b1.f0.linear_combination(a1, b1.f2, a2, b2.f2, a3, b3.f2, a4, b4.f2));
    }

    /** {@inherit_doc} */
    //override
    Field_Univariate_Derivative_2<T> get_pi() 
    {
        const T zero = get_value_field().get_zero();
        return Field_Univariate_Derivative_2<>(zero.get_pi(), zero, zero);
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

        if (dynamic_cast<const Field_Univariate_Derivative_2*>(*other) != nullptr) 
        {
            //@Suppress_Warnings("unchecked")
            const Field_Univariate_Derivative_2<T>& rhs = (Field_Univariate_Derivative_2<T>) other;
            return my_f0.equals(rhs.f0) && my_f1.equals(rhs.f1) && my_f2.equals(rhs.f2);
        }

        return false;

    }

    /** Get a hash_code for the univariate derivative.
     * @return a hash code value for this object
     */
    //override
    int hash_code() 
    {
        return 317 - 41 * my_f0.hash_code() + 57 * my_f1.hash_code() - 103 * my_f2.hash_code();
    }
};