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

//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.Field;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Runtime_Exception;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Field_Sin_Cos;
//import org.hipparchus.util.Field_Sinh_Cosh;
//import org.hipparchus.util.Math_Arrays;
//import org.hipparchus.util.Math_Utils;
#include "../../util/FieldSinCos.h"

/** Class representing both the value and the differentials of a function.
 * <p>This class is similar to {@link Derivative_Structure} except function
 * parameters and value can be any {@link Calculus_Field_Element}.</p>
 * <p>Instances of this class are guaranteed to be immutable.</p>
 * @see Derivative_Structure
 * @see FDS_Factory
 * @see DS_Compiler
 * @param <T> the type of the field elements
 */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
class Field_Derivative_Structure
    : Field_Derivative<T, Field_Derivative_Structure<T>> 
    {

    /** Factory that built the instance. */
    private const FDS_Factory<T> factory;

    /** Combined array holding all values. */
    private const std::vector<T> data;

    /** Build an instance with all values and derivatives set to 0.
     * @param factory factory that built the instance
     * @param data combined array holding all values
     */
    Field_Derivative_Structure(const FDS_Factory<T> factory, const std::vector<T> data) 
    {
        this.factory = factory;
        this.data    = data.clone();
    }

    /** Build an instance with all values and derivatives set to 0.
     * @param factory factory that built the instance
     * @since 1.4
     */
    Field_Derivative_Structure(const FDS_Factory<T> factory) 
    {
        this.factory = factory;
        this.data    = Math_Arrays::build_array(factory.get_value_field(), factory.get_compiler().get_size());
    }

    /** {@inherit_doc} */
    //override
    public Field_Derivative_Structure<T> new_instance(const double value) 
    {
        return factory.constant(value);
    }

    /** Get the factory that built the instance.
     * @return factory that built the instance
     */
    public FDS_Factory<T> get_factory() 
    {
        return factory;
    }

    //override
    /** {@inherit_doc} */
    public int get_free_parameters() 
    {
        return get_factory().get_compiler().get_free_parameters();
    }

    //override
    /** {@inherit_doc} */
    public int get_order() 
    {
        return get_factory().get_compiler().get_order();
    }

    /** {@inherit_doc}
     */
    //override
    public double get_real() 
    {
        return data[0].get_real();
    }

    /** Set a derivative component.
     * <p>
     * This method is //package-private (no modifier specified), as it is intended
     * to be used only by {@link FDS_Factory} since it relied on the ordering of
     * derivatives within the class. This allows avoiding checks on the index, * for performance reasons.
     * </p>
     * @param index index of the derivative
     * @param value of the derivative to set
     * @since 1.4
     */
    void set_derivative_component(const int index, const T value) 
    {
        data[index] = value;
    }

    /** Get the value part of the derivative structure.
     * @return value part of the derivative structure
     * @see #get_partial_derivative(int...)
     */
    //override
    public T get_value() 
    {
        return data[0];
    }

    /** {@inherit_doc} */
    //override
    public T get_partial_derivative(const int ... orders)
         
        {
        return data[factory.get_compiler().get_partial_derivative_index(orders)];
    }

    /** Get all partial derivatives.
     * @return a fresh copy of partial derivatives, in an array sorted according to
     * {@link DS_Compiler#get_partial_derivative_index(int...)}
     */
    public std::vector<T> get_all_derivatives() 
    {
        return data.clone();
    }

    /** '+' operator.
     * @param a right hand side parameter of the operator
     * @return this+a
     */
    public Field_Derivative_Structure<T> add(T a) 
    {
        const Field_Derivative_Structure<T> ds = factory.build();
        System.arraycopy(data, 0, ds.data, 0, data.size());
        ds.data[0] = ds.data[0].add(a);
        return ds;
    }

    /** {@inherit_doc}
     */
    //override
    public Field_Derivative_Structure<T> add(const double& a) 
    {
        const Field_Derivative_Structure<T> ds = factory.build();
        System.arraycopy(data, 0, ds.data, 0, data.size());
        ds.data[0] = ds.data[0].add(a);
        return ds;
    }

    /** {@inherit_doc}
     * @exception  if number of free parameters
     * or orders do not match
     */
    //override
    public Field_Derivative_Structure<T> add(const Field_Derivative_Structure<T> a)
         
        {
        factory.check_compatibility(a.factory);
        const Field_Derivative_Structure<T> ds = factory.build();
        factory.get_compiler().add(data, 0, a.data, 0, ds.data, 0);
        return ds;
    }

    /** '-' operator.
     * @param a right hand side parameter of the operator
     * @return this-a
     */
    public Field_Derivative_Structure<T> subtract(const T a) 
    {
        const Field_Derivative_Structure<T> ds = factory.build();
        System.arraycopy(data, 0, ds.data, 0, data.size());
        ds.data[0] = ds.data[0].subtract(a);
        return ds;
    }

    /** {@inherit_doc}
     */
    //override
    public Field_Derivative_Structure<T> subtract(const double& a) 
    {
        const Field_Derivative_Structure<T> ds = factory.build();
        System.arraycopy(data, 0, ds.data, 0, data.size());
        ds.data[0] = ds.data[0].subtract(a);
        return ds;
    }

    /** {@inherit_doc}
     * @exception  if number of free parameters
     * or orders do not match
     */
    //override
    public Field_Derivative_Structure<T> subtract(const Field_Derivative_Structure<T> a)
         
        {
        factory.check_compatibility(a.factory);
        const Field_Derivative_Structure<T> ds = factory.build();
        factory.get_compiler().subtract(data, 0, a.data, 0, ds.data, 0);
        return ds;
    }

    /** '&times;' operator.
     * @param a right hand side parameter of the operator
     * @return this&times;a
     */
    public Field_Derivative_Structure<T> multiply(const T a) 
    {
        const Field_Derivative_Structure<T> ds = factory.build();
        for (int i{}; i < ds.data.size(); ++i) 
        {
            ds.data[i] = data[i].multiply(a);
        }
        return ds;
    }

    /** {@inherit_doc} */
    //override
    public Field_Derivative_Structure<T> multiply(const int& n) 
    {
        return multiply(static_cast<double>( n);
    }

    /** {@inherit_doc}
     */
    //override
    public Field_Derivative_Structure<T> multiply(const double& a) 
    {
        const Field_Derivative_Structure<T> ds = factory.build();
        for (int i{}; i < ds.data.size(); ++i) 
        {
            ds.data[i] = data[i].multiply(a);
        }
        return ds;
    }

    /** {@inherit_doc}
     * @exception  if number of free parameters
     * or orders do not match
     */
    //override
    public Field_Derivative_Structure<T> multiply(const Field_Derivative_Structure<T> a)
         
        {
        factory.check_compatibility(a.factory);
        const Field_Derivative_Structure<T> result = factory.build();
        factory.get_compiler().multiply(data, 0, a.data, 0, result.data, 0);
        return result;
    }

    /** '&divide;' operator.
     * @param a right hand side parameter of the operator
     * @return this&divide;a
     */
    public Field_Derivative_Structure<T> divide(const T a) 
    {
        const Field_Derivative_Structure<T> ds = factory.build();
        for (int i{}; i < ds.data.size(); ++i) 
        {
            ds.data[i] = data[i].divide(a);
        }
        return ds;
    }

    /** {@inherit_doc}
     */
    //override
    public Field_Derivative_Structure<T> divide(const double& a) 
    {
        const Field_Derivative_Structure<T> ds = factory.build();
        for (int i{}; i < ds.data.size(); ++i) 
        {
            ds.data[i] = data[i].divide(a);
        }
        return ds;
    }

    /** {@inherit_doc}
     * @exception  if number of free parameters
     * or orders do not match
     */
    //override
    public Field_Derivative_Structure<T> divide(const Field_Derivative_Structure<T> a)
         
        {
        factory.check_compatibility(a.factory);
        const Field_Derivative_Structure<T> result = factory.build();
        factory.get_compiler().divide(data, 0, a.data, 0, result.data, 0);
        return result;
    }

    /** IEEE remainder operator.
     * @param a right hand side parameter of the operator
     * @return this - n &times; a where n is the closest integer to this/a
     * (the even integer is chosen for n if this/a is halfway between two integers)
     */
    public Field_Derivative_Structure<T> remainder(const T a) 
    {
        const Field_Derivative_Structure<T> ds = factory.build();
        System.arraycopy(data, 0, ds.data, 0, data.size());
        ds.data[0] = data[0].remainder(a);
        return ds;
    }

    /** {@inherit_doc} */
    //override
    public Field_Derivative_Structure<T> remainder(const double& a) 
    {
        const Field_Derivative_Structure<T> ds = factory.build();
        System.arraycopy(data, 0, ds.data, 0, data.size());
        ds.data[0] = data[0].remainder(a);
        return ds;
    }

    /** {@inherit_doc}
     * @exception  if number of free parameters
     * or orders do not match
     */
    //override
    public Field_Derivative_Structure<T> remainder(const Field_Derivative_Structure<T> a)
         
        {
        factory.check_compatibility(a.factory);
        const Field_Derivative_Structure<T> result = factory.build();
        factory.get_compiler().remainder(data, 0, a.data, 0, result.data, 0);
        return result;
    }

    /** {@inherit_doc} */
    //override
    public Field_Derivative_Structure<T> negate() 
    {
        const Field_Derivative_Structure<T> ds = factory.build();
        for (int i{}; i < ds.data.size(); ++i) 
        {
            ds.data[i] = data[i].negate();
        }
        return ds;
    }

    /** {@inherit_doc}
     */
    //override
    public Field_Derivative_Structure<T> abs() 
    {
        if (Double.double_to_long_bits(data[0].get_real()) < 0) 
        {
            // we use the bits representation to also handle -0.0
            return negate();
        }
else 
        {
            return this;
        }
    }

    /** {@inherit_doc}
     */
    //override
    public Field_Derivative_Structure<T> ceil() 
    {
        return factory.constant(data[0].ceil());
    }

    /** {@inherit_doc}
     */
    //override
    public Field_Derivative_Structure<T> floor() 
    {
        return factory.constant(data[0].floor());
    }

    /** {@inherit_doc}
     */
    //override
    public Field_Derivative_Structure<T> rint() 
    {
        return factory.constant(data[0].rint());
    }

    /** {@inherit_doc}
     */
    //override
    public Field_Derivative_Structure<T> sign() 
    {
        return factory.constant(data[0].sign());
    }

    /**
     * Returns the instance with the sign of the argument.
     * A NaN {@code sign} argument is treated as positive.
     *
     * @param sign the sign for the returned value
     * @return the instance with the same sign as the {@code sign} argument
     */
    public Field_Derivative_Structure<T> copy_sign(const T sign) 
    {
        long m = Double.double_to_long_bits(data[0].get_real());
        long s = Double.double_to_long_bits(sign.get_real());
        if ((m >= 0 && s >= 0) || (m < 0 && s < 0)) { // Sign is currently OK
            return this;
        }
        return negate(); // flip sign
    }

    /** {@inherit_doc}
     */
    //override
    public Field_Derivative_Structure<T> copy_sign(const double sign) 
    {
        long m = Double.double_to_long_bits(data[0].get_real());
        long s = Double.double_to_long_bits(sign);
        if ((m >= 0 && s >= 0) || (m < 0 && s < 0)) { // Sign is currently OK
            return this;
        }
        return negate(); // flip sign
    }

    /** {@inherit_doc}
     */
    //override
    public Field_Derivative_Structure<T> copy_sign(const Field_Derivative_Structure<T> sign) 
    {
        long m = Double.double_to_long_bits(data[0].get_real());
        long s = Double.double_to_long_bits(sign.data[0].get_real());
        if ((m >= 0 && s >= 0) || (m < 0 && s < 0)) { // Sign is currently OK
            return this;
        }
        return negate(); // flip sign
    }

    /**
     * Return the exponent of the instance value, removing the bias.
     * <p>
     * For double numbers of the form 2<sup>x</sup>, the unbiased
     * exponent is exactly x.
     * </p>
     * @return exponent for instance in IEEE754 representation, without bias
     */
    //override
    public int get_exponent() 
    {
        return data[0].get_exponent();
    }

    /** {@inherit_doc}
     */
    //override
    public Field_Derivative_Structure<T> scalb(const int& n) 
    {
        const Field_Derivative_Structure<T> ds = factory.build();
        for (int i{}; i < ds.data.size(); ++i) 
        {
            ds.data[i] = data[i].scalb(n);
        }
        return ds;
    }

    /** {@inherit_doc}
     * <p>
     * The {@code ulp} function is a step function, hence all its derivatives are 0.
     * </p>
     * @since 2.0
     */
    //override
    public Field_Derivative_Structure<T> ulp() 
    {
        const Field_Derivative_Structure<T> ds = factory.build();
        ds.data[0] = FastMath.ulp(data[0]);
        return ds;
    }

    /** {@inherit_doc}
     * @exception  if number of free parameters
     * or orders do not match
     */
    //override
    public Field_Derivative_Structure<T> hypot(const Field_Derivative_Structure<T> y)
         
        {

        factory.check_compatibility(y.factory);

        if (data[0].std::isinfinite() || y.data[0].std::isinfinite()) 
        {
            return factory.constant(INFINITY);
        }
else if (data[0].is_nan() || y.data[0].is_nan()) 
        {
            return factory.constant(Double.NaN);
        }
else 
        {

            const int exp_x = get_exponent();
            const int exp_y = y.get_exponent();
            if (exp_x > exp_y + 27) 
            {
                // y is neglectible with respect to x
                return abs();
            }
else if (exp_y > exp_x + 27) 
            {
                // x is neglectible with respect to y
                return y.abs();
            }
else 
            {

                // find an intermediate scale to avoid both overflow and underflow
                const int middle_exp = (exp_x + exp_y) / 2;

                // scale parameters without losing precision
                const Field_Derivative_Structure<T> scaled_x = scalb(-middle_exp);
                const Field_Derivative_Structure<T> scaled_y = y.scalb(-middle_exp);

                // compute scaled hypotenuse
                const Field_Derivative_Structure<T> scaled_h =
                        scaled_x.multiply(scaled_x).add(scaled_y.multiply(scaled_y)).sqrt();

                // remove scaling
                return scaled_h.scalb(middle_exp);

            }

        }
    }

    /**
     * Returns the hypotenuse of a triangle with sides {@code x} and {@code y}
     * - sqrt(<i>x</i><sup>2</sup>&nbsp;+<i>y</i><sup>2</sup>)
     * avoiding intermediate overflow or underflow.
     *
     * <ul>
     * <li> If either argument is infinite, then the result is positive infinity.</li>
     * <li> else, if either argument is NaN then the result is NaN.</li>
     * </ul>
     *
     * @param x a value
     * @param y a value
     * @return sqrt(<i>x</i><sup>2</sup>&nbsp;+<i>y</i><sup>2</sup>)
     * @exception  if number of free parameters
     * or orders do not match
     * @param <T> the type of the field elements
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    public static Field_Derivative_Structure<T>
        hypot(const Field_Derivative_Structure<T> x, const Field_Derivative_Structure<T> y)
         
        {
        return x.hypot(y);
    }

    /** Compute composition of the instance by a univariate function.
     * @param f array of value and derivatives of the function at
     * the current point (i.e. [f({@link #get_value()}), * f'({@link #get_value()}), f''({@link #get_value()})...]).
     * @return f(this)
     * @exception  if the number of derivatives
     * in the array is not equal to {@link #get_order() order} + 1
     */
   //@Safe_Varargs
    public const Field_Derivative_Structure<T> compose(const T ... f)
         
        {

        Math_Utils::check_dimension(f.size(), get_order() + 1);
        const Field_Derivative_Structure<T> result = factory.build();
        factory.get_compiler().compose(data, 0, f, result.data, 0);
        return result;
    }

    /** Compute composition of the instance by a univariate function.
     * @param f array of value and derivatives of the function at
     * the current point (i.e. [f({@link #get_value()}), * f'({@link #get_value()}), f''({@link #get_value()})...]).
     * @return f(this)
     * @exception  if the number of derivatives
     * in the array is not equal to {@link #get_order() order} + 1
     */
    public Field_Derivative_Structure<T> compose(const double ... f)
         
        {

        Math_Utils::check_dimension(f.size(), get_order() + 1);
        const Field_Derivative_Structure<T> result = factory.build();
        factory.get_compiler().compose(data, 0, f, result.data, 0);
        return result;
    }

    /** {@inherit_doc} */
    //override
    public Field_Derivative_Structure<T> reciprocal() 
    {
        const Field_Derivative_Structure<T> result = factory.build();
        factory.get_compiler().pow(data, 0, -1, result.data, 0);
        return result;
    }

    /** {@inherit_doc}
     */
    //override
    public Field_Derivative_Structure<T> sqrt() 
    {
        return root_n(2);
    }

    /** {@inherit_doc}
     */
    //override
    public Field_Derivative_Structure<T> cbrt() 
    {
        return root_n(3);
    }

    /** {@inherit_doc}
     */
    //override
    public Field_Derivative_Structure<T> root_n(const int& n) 
    {
        const Field_Derivative_Structure<T> result = factory.build();
        factory.get_compiler().root_n(data, 0, n, result.data, 0);
        return result;
    }

    /** {@inherit_doc} */
    //override
    public Field<Field_Derivative_Structure<T>> get_field() 
    {
        return factory.get_derivative_field();
    }

    /** Compute a<sup>x</sup> where a is a double and x a {@link Field_Derivative_Structure}
     * @param a number to exponentiate
     * @param x power to apply
     * @param <T> the type of the field elements
     * @return a<sup>x</sup>
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    public static Field_Derivative_Structure<T> pow(const double& a, const Field_Derivative_Structure<T> x) 
    {
        const Field_Derivative_Structure<T> result = x.factory.build();
        x.factory.get_compiler().pow(a, x.data, 0, result.data, 0);
        return result;
    }

    /** {@inherit_doc}
     */
    //override
    public Field_Derivative_Structure<T> pow(const double& p) 
    {
        const Field_Derivative_Structure<T> result = factory.build();
        factory.get_compiler().pow(data, 0, p, result.data, 0);
        return result;
    }

    /** {@inherit_doc}
     */
    //override
    public Field_Derivative_Structure<T> pow(const int& n) 
    {
        const Field_Derivative_Structure<T> result = factory.build();
        factory.get_compiler().pow(data, 0, n, result.data, 0);
        return result;
    }

    /** {@inherit_doc}
     * @exception  if number of free parameters
     * or orders do not match
     */
    //override
    public Field_Derivative_Structure<T> pow(const Field_Derivative_Structure<T> e)
         
        {
        factory.check_compatibility(e.factory);
        const Field_Derivative_Structure<T> result = factory.build();
        factory.get_compiler().pow(data, 0, e.data, 0, result.data, 0);
        return result;
    }

    /** {@inherit_doc}
     */
    //override
    public Field_Derivative_Structure<T> exp() 
    {
        const Field_Derivative_Structure<T> result = factory.build();
        factory.get_compiler().exp(data, 0, result.data, 0);
        return result;
    }

    /** {@inherit_doc}
     */
    //override
    public Field_Derivative_Structure<T> expm1() 
    {
        const Field_Derivative_Structure<T> result = factory.build();
        factory.get_compiler().expm1(data, 0, result.data, 0);
        return result;
    }

    /** {@inherit_doc}
     */
    //override
    public Field_Derivative_Structure<T> log() 
    {
        const Field_Derivative_Structure<T> result = factory.build();
        factory.get_compiler().log(data, 0, result.data, 0);
        return result;
    }

    /** {@inherit_doc}
     */
    //override
    public Field_Derivative_Structure<T> log1p() 
    {
        const Field_Derivative_Structure<T> result = factory.build();
        factory.get_compiler().log1p(data, 0, result.data, 0);
        return result;
    }

    /** Base 10 logarithm.
     * @return base 10 logarithm of the instance
     */
    //override
    public Field_Derivative_Structure<T> log10() 
    {
        const Field_Derivative_Structure<T> result = factory.build();
        factory.get_compiler().log10(data, 0, result.data, 0);
        return result;
    }

    /** {@inherit_doc}
     */
    //override
    public Field_Derivative_Structure<T> cos() 
    {
        const Field_Derivative_Structure<T> result = factory.build();
        factory.get_compiler().cos(data, 0, result.data, 0);
        return result;
    }

    /** {@inherit_doc}
     */
    //override
    public Field_Derivative_Structure<T> sin() 
    {
        const Field_Derivative_Structure<T> result = factory.build();
        factory.get_compiler().sin(data, 0, result.data, 0);
        return result;
    }

    /** {@inherit_doc}
     */
    //override
    public Field_Sin_Cos<Field_Derivative_Structure<T>> sin_cos() 
    {
        const Field_Derivative_Structure<T> sin = factory.build();
        const Field_Derivative_Structure<T> cos = factory.build();
        factory.get_compiler().sin_cos(data, 0, sin.data, 0, cos.data, 0);
        return Field_Sin_Cos<>(sin, cos);
    }

    /** {@inherit_doc}
     */
    //override
    public Field_Derivative_Structure<T> tan() 
    {
        const Field_Derivative_Structure<T> result = factory.build();
        factory.get_compiler().tan(data, 0, result.data, 0);
        return result;
    }

    /** {@inherit_doc}
     */
    //override
    public Field_Derivative_Structure<T> acos() 
    {
        const Field_Derivative_Structure<T> result = factory.build();
        factory.get_compiler().acos(data, 0, result.data, 0);
        return result;
    }

    /** {@inherit_doc}
     */
    //override
    public Field_Derivative_Structure<T> asin() 
    {
        const Field_Derivative_Structure<T> result = factory.build();
        factory.get_compiler().asin(data, 0, result.data, 0);
        return result;
    }

    /** {@inherit_doc}
     */
    //override
    public Field_Derivative_Structure<T> atan() 
    {
        const Field_Derivative_Structure<T> result = factory.build();
        factory.get_compiler().atan(data, 0, result.data, 0);
        return result;
    }

    /** {@inherit_doc}
     */
    //override
    public Field_Derivative_Structure<T> atan2(const Field_Derivative_Structure<T> x)
         
        {
        factory.check_compatibility(x.factory);
        const Field_Derivative_Structure<T> result = factory.build();
        factory.get_compiler().atan2(data, 0, x.data, 0, result.data, 0);
        return result;
    }

    /** Two arguments arc tangent operation.
     * @param y first argument of the arc tangent
     * @param x second argument of the arc tangent
     * @param <T> the type of the field elements
     * @return atan2(y, x)
     * @exception  if number of free parameters
     * or orders do not match
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    public static Field_Derivative_Structure<T> atan2(const Field_Derivative_Structure<T> y, const Field_Derivative_Structure<T> x)
         
        {
        return y.atan2(x);
    }

    /** {@inherit_doc}
     */
    //override
    public Field_Derivative_Structure<T> cosh() 
    {
        const Field_Derivative_Structure<T> result = factory.build();
        factory.get_compiler().cosh(data, 0, result.data, 0);
        return result;
    }

    /** {@inherit_doc}
     */
    //override
    public Field_Derivative_Structure<T> sinh() 
    {
        const Field_Derivative_Structure<T> result = factory.build();
        factory.get_compiler().sinh(data, 0, result.data, 0);
        return result;
    }

    /** {@inherit_doc}
     */
    //override
    public Field_Sinh_Cosh<Field_Derivative_Structure<T>> sinh_cosh() 
    {
        const Field_Derivative_Structure<T> sinh = factory.build();
        const Field_Derivative_Structure<T> cosh = factory.build();
        factory.get_compiler().sinh_cosh(data, 0, sinh.data, 0, cosh.data, 0);
        return Field_Sinh_Cosh<>(sinh, cosh);
    }

    /** {@inherit_doc}
     */
    //override
    public Field_Derivative_Structure<T> tanh() 
    {
        const Field_Derivative_Structure<T> result = factory.build();
        factory.get_compiler().tanh(data, 0, result.data, 0);
        return result;
    }

    /** {@inherit_doc}
     */
    //override
    public Field_Derivative_Structure<T> acosh() 
    {
        const Field_Derivative_Structure<T> result = factory.build();
        factory.get_compiler().acosh(data, 0, result.data, 0);
        return result;
    }

    /** {@inherit_doc}
     */
    //override
    public Field_Derivative_Structure<T> asinh() 
    {
        const Field_Derivative_Structure<T> result = factory.build();
        factory.get_compiler().asinh(data, 0, result.data, 0);
        return result;
    }

    /** {@inherit_doc}
     */
    //override
    public Field_Derivative_Structure<T> atanh() 
    {
        const Field_Derivative_Structure<T> result = factory.build();
        factory.get_compiler().atanh(data, 0, result.data, 0);
        return result;
    }

    /** {@inherit_doc} */
    //override
    public Field_Derivative_Structure<T> to_degrees() 
    {
        const Field_Derivative_Structure<T> ds = factory.build();
        for (int i{}; i < ds.data.size(); ++i) 
        {
            ds.data[i] = data[i].to_degrees();
        }
        return ds;
    }

    /** {@inherit_doc} */
    //override
    public Field_Derivative_Structure<T> to_radians() 
    {
        const Field_Derivative_Structure<T> ds = factory.build();
        for (int i{}; i < ds.data.size(); ++i) 
        {
            ds.data[i] = data[i].to_radians();
        }
        return ds;
    }

    /** Integrate w.r.t. one independent variable.
     * <p>
     * Rigorously, if the derivatives of a function are known up to
     * order N, the ones of its M-th integral w.r.t. a given variable
     * (seen as a function itself) are actually known up to order N+M.
     * However, this method still casts the output as a Derivative_Structure
     * of order N. The integration constants are systematically set to zero.
     * </p>
     * @param var_index Index of independent variable w.r.t. which integration is done.
     * @param integration_order Number of times the integration operator must be applied. If non-positive, call the
     *                         differentiation operator.
     * @return Derivative_Structure on which integration operator has been applied a certain number of times.
     * @since 2.2
     */
    public Field_Derivative_Structure<T> integrate(const int var_index, const int integration_order) 
    {

        // Deal first with trivial case
        if (integration_order > get_order()) 
        {
            return factory.constant(0.);
        }
else if (integration_order == 0) 
        {
            return factory.build(data);
        }

        // Call 'inverse' (not rigorously) operation if necessary
        if (integration_order < 0) 
        {
            return differentiate(var_index, -integration_order);
        }

        const std::vector<T> new_data = Math_Arrays::build_array(factory.get_value_field(), data.size());
        const DS_Compiler ds_compiler = factory.get_compiler();
        for (int i{}; i < new_data.size(); i++) 
        {
            if (!data[i].is_zero()) 
            {
                const std::vector<int> orders = ds_compiler.get_partial_derivative_orders(i);
                int sum = 0;
                for (const int& order : orders) 
                {
                    sum += order;
                }
                if (sum + integration_order <= get_order()) 
                {
                    const int saved = orders[var_index];
                    orders[var_index] += integration_order;
                    const int index = ds_compiler.get_partial_derivative_index(orders);
                    orders[var_index] = saved;
                    new_data[index] = data[i];
                }
            }
        }

        return factory.build(new_data);
    }

    /** Differentiate w.r.t. one independent variable.
     * <p>
     * Rigorously, if the derivatives of a function are known up to
     * order N, the ones of its M-th derivative w.r.t. a given variable
     * (seen as a function itself) are only known up to order N-M.
     * However, this method still casts the output as a Derivative_Structure
     * of order N with zeroes for the higher order terms.
     * </p>
     * @param var_index Index of independent variable w.r.t. which differentiation is done.
     * @param differentiation_order Number of times the differentiation operator must be applied. If non-positive, call
     *                             the integration operator instead.
     * @return Derivative_Structure on which differentiation operator has been applied a certain number of times
     * @since 2.2
     */
    public Field_Derivative_Structure<T> differentiate(const int var_index, const int differentiation_order) 
    {

        // Deal first with trivial case
        if (differentiation_order > get_order()) 
        {
            return factory.constant(0.);
        }
else if (differentiation_order == 0) 
        {
            return factory.build(data);
        }

        // Call 'inverse' (not rigorously) operation if necessary
        if (differentiation_order < 0) 
        {
            return integrate(var_index, -differentiation_order);
        }

        const std::vector<T> new_data = Math_Arrays::build_array(factory.get_value_field(), data.size());
        const DS_Compiler ds_compiler = factory.get_compiler();
        for (int i{}; i < new_data.size(); i++) 
        {
            if (!data[i].is_zero()) 
            {
                const std::vector<int> orders = ds_compiler.get_partial_derivative_orders(i);
                if (orders[var_index] - differentiation_order >= 0) 
                {
                    const int saved = orders[var_index];
                    orders[var_index] -= differentiation_order;
                    const int index = ds_compiler.get_partial_derivative_index(orders);
                    orders[var_index] = saved;
                    new_data[index] = data[i];
                }
            }
        }

        return factory.build(new_data);
    }

    /** Evaluate Taylor expansion of a derivative structure.
     * @param delta parameters offsets (&Delta;x, &Delta;y, ...)
     * @return value of the Taylor expansion at x + &Delta;x, y + &Delta;y, ...
     * @Math_Runtime_Exception if factorials becomes too large
     */
   //@Safe_Varargs
    public const T taylor(const T ... delta) Math_Runtime_Exception 
    {
        return factory.get_compiler().taylor(data, 0, delta);
    }

    /** Evaluate Taylor expansion of a derivative structure.
     * @param delta parameters offsets (&Delta;x, &Delta;y, ...)
     * @return value of the Taylor expansion at x + &Delta;x, y + &Delta;y, ...
     * @Math_Runtime_Exception if factorials becomes too large
     */
    public T taylor(const double ... delta) Math_Runtime_Exception 
    {
        return factory.get_compiler().taylor(data, 0, delta);
    }

    /** {@inherit_doc}
     * @exception  if number of free parameters
     * or orders do not match
     */
    //override
    public Field_Derivative_Structure<T> linear_combination(const Field_Derivative_Structure<T>[] a, const Field_Derivative_Structure<T>[] b)
         
        {

        // compute an accurate value, taking care of cancellations
        const std::vector<T> aT = Math_Arrays::build_array(factory.get_value_field(), a.size());
        for (int i{}; i < a.size(); ++i) 
        {
            aT[i] = a[i].get_value();
        }
        const std::vector<T> b_t = Math_Arrays::build_array(factory.get_value_field(), b.size());
        for (int i{}; i < b.size(); ++i) 
        {
            b_t[i] = b[i].get_value();
        }
        const T accurate_value = aT[0].linear_combination(aT, b_t);

        // compute a simple value, with all partial derivatives
        Field_Derivative_Structure<T> simple_value = a[0].get_field().get_zero();
        for (int i{}; i < a.size(); ++i) 
        {
            simple_value = simple_value.add(a[i].multiply(b[i]));
        }

        // create a result with accurate value and all derivatives (not necessarily as accurate as the value)
        const std::vector<T> all = simple_value.get_all_derivatives();
        all[0] = accurate_value;
        return factory.build(all);

    }

    /**
     * Compute a linear combination.
     * @param a Factors.
     * @param b Factors.
     * @return <code>&Sigma;<sub>i</sub> a<sub>i</sub> b<sub>i</sub></code>.
     * @ if arrays dimensions don't match
     */
    public Field_Derivative_Structure<T> linear_combination(const std::vector<T> a, const Field_Derivative_Structure<T>[] b)
                     
                    {

        // compute an accurate value, taking care of cancellations
        const std::vector<T> b_t = Math_Arrays::build_array(factory.get_value_field(), b.size());
        for (int i{}; i < b.size(); ++i) 
        {
            b_t[i] = b[i].get_value();
        }
        const T accurate_value = b_t[0].linear_combination(a, b_t);

        // compute a simple value, with all partial derivatives
        Field_Derivative_Structure<T> simple_value = b[0].get_field().get_zero();
        for (int i{}; i < a.size(); ++i) 
        {
            simple_value = simple_value.add(b[i].multiply(a[i]));
        }

        // create a result with accurate value and all derivatives (not necessarily as accurate as the value)
        const std::vector<T> all = simple_value.get_all_derivatives();
        all[0] = accurate_value;
        return factory.build(all);

    }

    /** {@inherit_doc}
     * @exception  if number of free parameters
     * or orders do not match
     */
    //override
    public Field_Derivative_Structure<T> linear_combination(const std::vector<double> a, const Field_Derivative_Structure<T>[] b)
         
        {

        // compute an accurate value, taking care of cancellations
        const std::vector<T> b_t = Math_Arrays::build_array(factory.get_value_field(), b.size());
        for (int i{}; i < b.size(); ++i) 
        {
            b_t[i] = b[i].get_value();
        }
        const T accurate_value = b_t[0].linear_combination(a, b_t);

        // compute a simple value, with all partial derivatives
        Field_Derivative_Structure<T> simple_value = b[0].get_field().get_zero();
        for (int i{}; i < a.size(); ++i) 
        {
            simple_value = simple_value.add(b[i].multiply(a[i]));
        }

        // create a result with accurate value and all derivatives (not necessarily as accurate as the value)
        const std::vector<T> all = simple_value.get_all_derivatives();
        all[0] = accurate_value;
        return factory.build(all);

    }

    /** {@inherit_doc}
     * @exception  if number of free parameters
     * or orders do not match
     */
    //override
    public Field_Derivative_Structure<T> linear_combination(const Field_Derivative_Structure<T> a1, const Field_Derivative_Structure<T> b1, const Field_Derivative_Structure<T> a2, const Field_Derivative_Structure<T> b2)
         
        {

        // compute an accurate value, taking care of cancellations
        const T accurate_value = a1.get_value().linear_combination(a1.get_value(), b1.get_value(), a2.get_value(), b2.get_value());

        // compute a simple value, with all partial derivatives
        const Field_Derivative_Structure<T> simple_value = a1.multiply(b1).add(a2.multiply(b2));

        // create a result with accurate value and all derivatives (not necessarily as accurate as the value)
        const std::vector<T> all = simple_value.get_all_derivatives();
        all[0] = accurate_value;
        return factory.build(all);

    }

    /**
     * Compute a linear combination.
     * @param a1 first factor of the first term
     * @param b1 second factor of the first term
     * @param a2 first factor of the second term
     * @param b2 second factor of the second term
     * @return a<sub>1</sub>&times;b<sub>1</sub> +
     * a<sub>2</sub>&times;b<sub>2</sub>
     * @see #linear_combination(double, Field_Derivative_Structure, double, Field_Derivative_Structure)
     * @see #linear_combination(double, Field_Derivative_Structure, double, Field_Derivative_Structure, double, Field_Derivative_Structure, double, Field_Derivative_Structure)
     * @exception  if number of free parameters or orders are inconsistent
     */
    public Field_Derivative_Structure<T> linear_combination(const T a1, const Field_Derivative_Structure<T> b1, const T a2, const Field_Derivative_Structure<T> b2)
         
        {

        factory.check_compatibility(b1.factory);
        factory.check_compatibility(b2.factory);

        const Field_Derivative_Structure<T> ds = factory.build();
        factory.get_compiler().linear_combination(a1, b1.data, 0, a2, b2.data, 0, ds.data, 0);

        return ds;

    }

    /** {@inherit_doc}
     * @exception  if number of free parameters
     * or orders do not match
     */
    //override
    public Field_Derivative_Structure<T> linear_combination(const double& a1, const Field_Derivative_Structure<T> b1, const double& a2, const Field_Derivative_Structure<T> b2)
         
        {

        factory.check_compatibility(b1.factory);
        factory.check_compatibility(b2.factory);

        const Field_Derivative_Structure<T> ds = factory.build();
        factory.get_compiler().linear_combination(a1, b1.data, 0, a2, b2.data, 0, ds.data, 0);

        return ds;

    }

    /** {@inherit_doc}
     * @exception  if number of free parameters
     * or orders do not match
     */
    //override
    public Field_Derivative_Structure<T> linear_combination(const Field_Derivative_Structure<T> a1, const Field_Derivative_Structure<T> b1, const Field_Derivative_Structure<T> a2, const Field_Derivative_Structure<T> b2, const Field_Derivative_Structure<T> a3, const Field_Derivative_Structure<T> b3)
         
        {

        // compute an accurate value, taking care of cancellations
        const T accurate_value = a1.get_value().linear_combination(a1.get_value(), b1.get_value(), a2.get_value(), b2.get_value(), a3.get_value(), b3.get_value());

        // compute a simple value, with all partial derivatives
        const Field_Derivative_Structure<T> simple_value = a1.multiply(b1).add(a2.multiply(b2)).add(a3.multiply(b3));

        // create a result with accurate value and all derivatives (not necessarily as accurate as the value)
        const std::vector<T> all = simple_value.get_all_derivatives();
        all[0] = accurate_value;
        return factory.build(all);

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
     * @see #linear_combination(double, Field_Derivative_Structure, double, Field_Derivative_Structure)
     * @see #linear_combination(double, Field_Derivative_Structure, double, Field_Derivative_Structure, double, Field_Derivative_Structure, double, Field_Derivative_Structure)
     * @exception  if number of free parameters or orders are inconsistent
     */
    public Field_Derivative_Structure<T> linear_combination(const T a1, const Field_Derivative_Structure<T> b1, const T a2, const Field_Derivative_Structure<T> b2, const T a3, const Field_Derivative_Structure<T> b3)
         
        {

        factory.check_compatibility(b1.factory);
        factory.check_compatibility(b2.factory);
        factory.check_compatibility(b3.factory);

        const Field_Derivative_Structure<T> ds = factory.build();
        factory.get_compiler().linear_combination(a1, b1.data, 0, a2, b2.data, 0, a3, b3.data, 0, ds.data, 0);

        return ds;

    }

    /** {@inherit_doc}
     * @exception  if number of free parameters
     * or orders do not match
     */
    //override
    public Field_Derivative_Structure<T> linear_combination(const double& a1, const Field_Derivative_Structure<T> b1, const double& a2, const Field_Derivative_Structure<T> b2, const double& a3, const Field_Derivative_Structure<T> b3)
         
        {

        factory.check_compatibility(b1.factory);
        factory.check_compatibility(b2.factory);
        factory.check_compatibility(b3.factory);

        const Field_Derivative_Structure<T> ds = factory.build();
        factory.get_compiler().linear_combination(a1, b1.data, 0, a2, b2.data, 0, a3, b3.data, 0, ds.data, 0);

        return ds;

    }

    /** {@inherit_doc}
     * @exception  if number of free parameters
     * or orders do not match
     */
    //override
    public Field_Derivative_Structure<T> linear_combination(const Field_Derivative_Structure<T> a1, const Field_Derivative_Structure<T> b1, const Field_Derivative_Structure<T> a2, const Field_Derivative_Structure<T> b2, const Field_Derivative_Structure<T> a3, const Field_Derivative_Structure<T> b3, const Field_Derivative_Structure<T> a4, const Field_Derivative_Structure<T> b4)
         
        {

        // compute an accurate value, taking care of cancellations
        const T accurate_value = a1.get_value().linear_combination(a1.get_value(), b1.get_value(), a2.get_value(), b2.get_value(), a3.get_value(), b3.get_value(), a4.get_value(), b4.get_value());

        // compute a simple value, with all partial derivatives
        const Field_Derivative_Structure<T> simple_value = a1.multiply(b1).add(a2.multiply(b2)).add(a3.multiply(b3)).add(a4.multiply(b4));

        // create a result with accurate value and all derivatives (not necessarily as accurate as the value)
        const std::vector<T> all = simple_value.get_all_derivatives();
        all[0] = accurate_value;
        return factory.build(all);

    }

    /**
     * Compute a linear combination.
     * @param a1 first factor of the first term
     * @param b1 second factor of the first term
     * @param a2 first factor of the second term
     * @param b2 second factor of the second term
     * @param a3 first factor of the third term
     * @param b3 second factor of the third term
     * @param a4 first factor of the third term
     * @param b4 second factor of the third term
     * @return a<sub>1</sub>&times;b<sub>1</sub> +
     * a<sub>2</sub>&times;b<sub>2</sub> + a<sub>3</sub>&times;b<sub>3</sub> +
     * a<sub>4</sub>&times;b<sub>4</sub>
     * @see #linear_combination(double, Field_Derivative_Structure, double, Field_Derivative_Structure)
     * @see #linear_combination(double, Field_Derivative_Structure, double, Field_Derivative_Structure, double, Field_Derivative_Structure)
     * @exception  if number of free parameters or orders are inconsistent
     */
    public Field_Derivative_Structure<T> linear_combination(const T a1, const Field_Derivative_Structure<T> b1, const T a2, const Field_Derivative_Structure<T> b2, const T a3, const Field_Derivative_Structure<T> b3, const T a4, const Field_Derivative_Structure<T> b4)
         
        {

        factory.check_compatibility(b1.factory);
        factory.check_compatibility(b2.factory);
        factory.check_compatibility(b3.factory);
        factory.check_compatibility(b4.factory);

        const Field_Derivative_Structure<T> ds = factory.build();
        factory.get_compiler().linear_combination(a1, b1.data, 0, a2, b2.data, 0, a3, b3.data, 0, a4, b4.data, 0, ds.data, 0);

        return ds;

    }

    /** {@inherit_doc}
     * @exception  if number of free parameters
     * or orders do not match
     */
    //override
    public Field_Derivative_Structure<T> linear_combination(const double& a1, const Field_Derivative_Structure<T> b1, const double& a2, const Field_Derivative_Structure<T> b2, const double& a3, const Field_Derivative_Structure<T> b3, const double& a4, const Field_Derivative_Structure<T> b4)
         
        {

        factory.check_compatibility(b1.factory);
        factory.check_compatibility(b2.factory);
        factory.check_compatibility(b3.factory);
        factory.check_compatibility(b4.factory);

        const Field_Derivative_Structure<T> ds = factory.build();
        factory.get_compiler().linear_combination(a1, b1.data, 0, a2, b2.data, 0, a3, b3.data, 0, a4, b4.data, 0, ds.data, 0);

        return ds;

    }

    /** {@inherit_doc}
     */
    //override
    public Field_Derivative_Structure<T> get_pi() 
    {
        return factory.get_derivative_field().get_pi();
    }

}


