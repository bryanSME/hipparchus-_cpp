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
#include "../../util/FieldSinCos.h"
//import java.util.Arrays;

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
 * with {@link Field_Derivative_Structure#get_order() derivation order} limited to one.
 * It should have less overhead than {@link Field_Derivative_Structure} in its domain.</p>
 * <p>This class is an implementation of Rall's numbers. Rall's numbers are an
 * extension to the real numbers used throughout mathematical expressions; they hold
 * the derivative together with the value of a function.</p>
 * <p>{@link Field_Gradient} instances can be used directly thanks to
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
 * @see Field_Univariate_Derivative_2
 * @since 1.7
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
class Field_Gradient : Field_Derivative<T, Field_Gradient<T>>
{
private:
	/** Value of the function. */
	const T my_value;

	/** Gradient of the function. */
	const std::vector<T> my_grad;

	/** Build an instance with values and unitialized derivatives array.
	 * @param value value of the function
	 * @param free_parameters number of free parameters
	 */
	Field_Gradient(const T& value, const int& free_parameters) : my_value{ value }
	{
		this.grad = Math_Arrays::build_array(value.get_field(), free_parameters);
	}

public:
	/** Build an instance with values and derivative.
	 * @param value value of the function
	 * @param my_gradient my_gradient of the function
	 */
	 /*@Safe_Varargs*/
	Field_Gradient(const T value, const T... my_gradient)
	{
		this(value, my_gradient.size());
		System.arraycopy(gradient, 0, my_grad, 0, my_grad.size());
	}

	/** Build an instance from a {@link Derivative_Structure}.
	 * @param ds derivative structure
	 * @exception  if {@code ds} order
	 * is not 1
	 */
	Field_Gradient(const Field_Derivative_Structure<T> ds)
	{
		this(ds.get_value(), ds.get_free_parameters());
		Math_Utils::check_dimension(ds.get_order(), 1);
		System.arraycopy(ds.get_all_derivatives(), 1, my_grad, 0, my_grad.size());
	}

	/** Build an instance corresponding to a constant value.
	 * @param free_parameters number of free parameters (i.e. dimension of the my_gradient)
	 * @param value constant value of the function
	 * @param <T> the type of the function parameters and value
	 * @return a {@code Field_Gradient} with a constant value and all derivatives set to 0.0
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	static  Field_Gradient<T> constant(const int free_parameters, const T value)
	{
		const Field_Gradient<T> g = Field_Gradient<>(value, free_parameters);
		Arrays.fill(g.grad, value.get_field().get_zero());
		return g;
	}

	/** Build a {@code Gradient} representing a variable.
	 * <p>Instances built using this method are considered
	 * to be the free variables with respect to which differentials
	 * are computed. As such, their differential with respect to
	 * themselves is +1.</p>
	 * @param free_parameters number of free parameters (i.e. dimension of the my_gradient)
	 * @param index index of the variable (from 0 to {@link #get_free_parameters() get_free_parameters()} - 1)
	 * @param value value of the variable
	 * @param <T> the type of the function parameters and value
	 * @return a {@code Field_Gradient} with a constant value and all derivatives set to 0.0 except the
	 * one at {@code index} which will be set to 1.0
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	static  Field_Gradient<T> variable(const int free_parameters, const int index, const T& value)
	{
		const Field_Gradient<T> g = Field_Gradient<>(value, free_parameters);
		const Field<T> field = value.get_field();
		Arrays.fill(g.grad, field.get_zero());
		g.grad[index] = field.get_one();
		return g;
	}

	/** Get the {@link Field} the value and parameters of the function belongs to.
	 * @return {@link Field} the value and parameters of the function belongs to
	 */
	Field<T> get_value_field() const
	{
		return my_value.get_field();
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> new_instance(const double& c)
	{
		return new_instance(get_value_field().get_zero().new_instance(c));
	}

	/** Create an instance corresponding to a constant real value.
	 * <p>
	 * The default implementation creates the instance by adding
	 * the value to {@code get_field().get_zero()}. This is not optimal
	 * and does not work when called with a negative zero as the
	 * sign of zero is lost with the addition. The default implementation
	 * should therefore be overridden in concrete classes. The default
	 * implementation will be removed at the next major version.
	 * </p>
	 * @param c constant real value
	 * @return instance corresponding to a constant real value
	 */
	Field_Gradient<T> new_instance(const T c)
	{
		return Field_Gradient<>(c, Math_Arrays::build_array(value.get_field(), my_grad.size()));
	}

	/** {@inherit_doc} */
	//override
	double get_real()
	{
		return get_value().get_real();
	}

	/** Get the value part of the function.
	 * @return value part of the value of the function
	 */
	 //override
	T get_value() const
	{
		return my_value;
	}

	/** Get the my_gradient part of the function.
	 * @return my_gradient part of the value of the function
	 */
	std::vector<T> get_gradient() const
	{
		return my_grad;
	}

	/** Get the number of free parameters.
	 * @return number of free parameters
	 */
	 //override
	int get_free_parameters() const
	{
		return my_grad.size();
	}

	/** {@inherit_doc} */
	//override
	int get_order() const
	{
		return 1;
	}

	/** {@inherit_doc} */
	//override
	T get_partial_derivative(const int ... orders)
	{
		// check the number of components
		if (orders.size() != my_grad.size())
		{
			throw std::exception("not implmented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, orders.size(), my_grad.size());
		}

		// check that either all derivation orders are set to 0, // or that only one is set to 1 and all other ones are set to 0
		int selected = -1;
		for (int i{}; i < orders.size(); ++i)
		{
			if (orders[i] != 0)
			{
				if (selected >= 0 || orders[i] != 1)
				{
					throw std::exception("not implmented");
					//throw (hipparchus::exception::Localized_Core_Formats_Type::DERIVATION_ORDER_NOT_ALLOWED, orders[i]);
				}
				// found the component set to derivation order 1
				selected = i;
			}
		}

		return (selected < 0)
			? value
			: my_grad[selected];
	}

	/** Get the partial derivative with respect to one parameter.
	 * @param n index of the parameter (counting from 0)
	 * @return partial derivative with respect to the n<sup>th</sup> parameter
	 * @exception  if n is either negative or larger
	 * or equal to {@link #get_free_parameters()}
	 */
	T get_partial_derivative(const int& n)
	{
		if (n < 0 || n >= my_grad.size())
		{
			throw std::exception("not implmented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::OUT_OF_RANGE_SIMPLE, n, 0, my_grad.size() - 1);
		}
		return my_grad[n];
	}

	/** Convert the instance to a {@link Field_Derivative_Structure}.
	 * @return derivative structure with same value and derivative as the instance
	 */
	Field_Derivative_Structure<T> to_derivative_structure()
	{
		const std::vector<T> derivatives = Math_Arrays::build_array(get_value_field(), 1 + my_grad.size());
		derivatives[0] = value;
		System.arraycopy(my_grad, 0, derivatives, 1, my_grad.size());
		return get_field().get_conversion_factory().build(derivatives);
	}

	/** '+' operator.
	 * @param a right hand side parameter of the operator
	 * @return this+a
	 */
	Field_Gradient<T> add(const T a)
	{
		return Field_Gradient<>(value.add(a), my_grad);
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> add(const double& a)
	{
		return Field_Gradient<>(value.add(a), my_grad);
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> add(const Field_Gradient<T>& a)
	{
		const Field_Gradient<T> result = new_instance(value.add(a.value));
		for (int i{}; i < my_grad.size(); ++i)
		{
			result.grad[i] = my_grad[i].add(a.grad[i]);
		}
		return result;
	}

	/** '-' operator.
	 * @param a right hand side parameter of the operator
	 * @return this-a
	 */
	Field_Gradient<T> subtract(const T a)
	{
		return Field_Gradient<>(value.subtract(a), my_grad);
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> subtract(const double& a)
	{
		return Field_Gradient<>(value.subtract(a), my_grad);
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> subtract(const Field_Gradient<T>& a)
	{
		const Field_Gradient<T> result = new_instance(my_value.subtract(a.value));
		for (int i{}; i < my_grad.size(); ++i)
		{
			result.grad[i] = my_grad[i].subtract(a.grad[i]);
		}
		return result;
	}

	/** '&times;' operator.
	 * @param n right hand side parameter of the operator
	 * @return this&times;n
	 */
	Field_Gradient<T> multiply(const T n)
	{
		const Field_Gradient<T> result = new_instance(value.multiply(n));
		for (int i{}; i < my_grad.size(); ++i)
		{
			result.grad[i] = my_grad[i].multiply(n);
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> multiply(const int& n)
	{
		const Field_Gradient<T> result = new_instance(value.multiply(n));
		for (int i{}; i < my_grad.size(); ++i)
		{
			result.grad[i] = my_grad[i].multiply(n);
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> multiply(const double& a)
	{
		const Field_Gradient<T> result = new_instance(value.multiply(a));
		for (int i{}; i < my_grad.size(); ++i)
		{
			result.grad[i] = my_grad[i].multiply(a);
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> multiply(const Field_Gradient<T>& a)
	{
		const Field_Gradient<T> result = new_instance(value.multiply(a.value));
		for (int i{}; i < my_grad.size(); ++i)
		{
			result.grad[i] = my_grad[i].multiply(a.value).add(value.multiply(a.grad[i]));
		}
		return result;
	}

	/** '&divide;' operator.
	 * @param a right hand side parameter of the operator
	 * @return this&divide;a
	 */
	Field_Gradient<T> divide(const T a)
	{
		const Field_Gradient<T> result = new_instance(value.divide(a));
		for (int i{}; i < my_grad.size(); ++i)
		{
			result.grad[i] = my_grad[i].divide(a);
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> divide(const double& a)
	{
		const Field_Gradient<T> result = new_instance(value.divide(a));
		for (int i{}; i < my_grad.size(); ++i)
		{
			result.grad[i] = my_grad[i].divide(a);
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> divide(const Field_Gradient<T>& a)
	{
		const T inv1 = a.value.reciprocal();
		const T inv2 = inv1.multiply(inv1);
		const Field_Gradient<T> result = new_instance(value.multiply(inv1));
		for (int i{}; i < my_grad.size(); ++i)
		{
			result.grad[i] = my_grad[i].multiply(a.value).subtract(value.multiply(a.grad[i])).multiply(inv2);
		}
		return result;
	}

	/** IEEE remainder operator.
	 * @param a right hand side parameter of the operator
	 * @return this - n &times; a where n is the closest integer to this/a
	 * (the even integer is chosen for n if this/a is halfway between two integers)
	 */
	Field_Gradient<T> remainder(const T a)
	{
		return Field_Gradient<>(std::remainder(value, a), my_grad);
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> remainder(const double& a)
	{
		return Field_Gradient<>(std::remainder(value, a), my_grad);
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> remainder(const Field_Gradient<T>& a)
	{
		// compute k such that lhs % rhs = lhs - k rhs
		const T rem = std::remainder(value, a.value);
		const T k = std::rint(value.subtract(rem).divide(a.value));

		const Field_Gradient<T> result = new_instance(rem);
		for (int i{}; i < my_grad.size(); ++i)
		{
			result.grad[i] = my_grad[i].subtract(k.multiply(a.grad[i]));
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> negate()
	{
		const Field_Gradient<T> result = new_instance(value.negate());
		for (int i{}; i < my_grad.size(); ++i)
		{
			result.grad[i] = my_grad[i].negate();
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> abs()
	{
		if (Double.double_to_long_bits(value.get_real()) < 0)
		{
			// we use the bits representation to also handle -0.0
			return negate();
		}
		else
		{
			return this;
		}
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> ceil()
	{
		return new_instance(std::ceil(value));
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> floor()
	{
		return new_instance(std::floor(value));
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> rint()
	{
		return new_instance(std::rint(value));
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> sign()
	{
		return new_instance(FastMath.sign(value));
	}

	/**
	 * Returns the instance with the sign of the argument.
	 * A NaN {@code sign} argument is treated as positive.
	 *
	 * @param sign the sign for the returned value
	 * @return the instance with the same sign as the {@code sign} argument
	 */
	Field_Gradient<T> copy_sign(const T sign)
	{
		long m = Double.double_to_long_bits(value.get_real());
		long s = Double.double_to_long_bits(sign.get_real());
		if ((m >= 0 && s >= 0) || (m < 0 && s < 0)) { // Sign is currently OK
			return this;
		}
		return negate(); // flip sign
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> copy_sign(const Field_Gradient<T> sign)
	{
		long m = Double.double_to_long_bits(value.get_real());
		long s = Double.double_to_long_bits(sign.value.get_real());
		if ((m >= 0 && s >= 0) || (m < 0 && s < 0)) { // Sign is currently OK
			return this;
		}
		return negate(); // flip sign
	};

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> copy_sign(const double sign)
	{
		long m = Double.double_to_long_bits(value.get_real());
		long s = Double.double_to_long_bits(sign);
		if ((m >= 0 && s >= 0) || (m < 0 && s < 0)) { // Sign is currently OK
			return this;
		}
		return negate(); // flip sign
	};

	/** {@inherit_doc} */
	//override
	int get_exponent()
	{
		return FastMath.get_exponent(my_value.get_real());
	};

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> scalb(const int& n)
	{
		const Field_Gradient<T> result = new_instance(std::scalbn(value, n));
		for (int i{}; i < my_grad.size(); ++i)
		{
			result.grad[i] = std::scalbn(grad[i], n);
		}
		return result;
	}

	/** {@inherit_doc}
	 * <p>
	 * The {@code ulp} function is a step function, hence all its derivatives are 0.
	 * </p>
	 * @since 2.0
	 */
	 //override
	Field_Gradient<T> ulp()
	{
		return new_instance(FastMath.ulp(value));
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> hypot(const Field_Gradient<T>& y)
	{
		if (std::isinf(value.get_real()) || std::isinf(y.value.get_real()))
		{
			return new_instance(INFINITY);
		}
		if (std::isnan(value.get_real()) || std::isnan(y.value.get_real()))
		{
			return new_instance(Double.NaN);
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
		const Field_Gradient<T> scaled_x = scalb(-middle_exp);
		const Field_Gradient<T> scaled_y = y.scalb(-middle_exp);

		// compute scaled hypotenuse
		const Field_Gradient<T> scaled_h =
			scaled_x.multiply(scaled_x).add(scaled_y.multiply(scaled_y)).sqrt();

		// remove scaling
		return scaled_h.scalb(middle_exp);
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> reciprocal()
	{
		const T inv1 = value.reciprocal();
		const T m_inv2 = inv1.multiply(inv1).negate();
		const Field_Gradient<T> result = new_instance(inv1);
		for (int i{}; i < my_grad.size(); ++i)
		{
			result.grad[i] = my_grad[i].multiply(m_inv2);
		}
		return result;
	}

	/** Compute composition of the instance by a function.
	 * @param g0 value of the function at the current point (i.e. at {@code g(get_value())})
	 * @param g1 first derivative of the function at the current point (i.e. at {@code g'(get_value())})
	 * @return g(this)
	 */
	Field_Gradient<T> compose(const T& g0, const T& g1)
	{
		const Field_Gradient<T> result = new_instance(g0);
		for (int i{}; i < my_grad.size(); ++i)
		{
			result.grad[i] = g1.multiply(grad[i]);
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> sqrt()
	{
		const T s = std::sqrt(value);
		return compose(s, s.add(s).reciprocal());
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> cbrt()
	{
		const T c = std::cbrt(value);
		return compose(c, c.multiply(c).multiply(3).reciprocal());
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> root_n(const int& n)
	{
		if (n == 2)
		{
			return sqrt();
		}
		if (n == 3)
		{
			return cbrt();
		}
		const T r = std::pow(value, 1.0 / n);
		return compose(r, std::pow(r, n - 1).multiply(n).reciprocal());
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient_Field<T> get_field()
	{
		return Field_Gradient_Field.get_field(get_value_field(), get_free_parameters());
	}

	/** Compute a<sup>x</sup> where a is a double and x a {@link Field_Gradient}
	 * @param a number to exponentiate
	 * @param x power to apply
	 * @param <T> the type of the function parameters and value
	 * @return a<sup>x</sup>
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	static  Field_Gradient<T> pow(const double& a, const Field_Gradient<T>& x)
	{
		if (a == 0)
		{
			return x.get_field().get_zero();
		}

		const T a_x = std::pow(x.value.new_instance(a), x.value);
		const T a_xln_a = a_x.multiply(std::log(a));
		const Field_Gradient<T> result = x.new_instance(a_x);
		for (int i{}; i < x.grad.size(); ++i)
		{
			result.grad[i] = a_xln_a.multiply(x.grad[i]);
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> pow(const double& p)
	{
		if (p == 0)
		{
			return get_field().get_one();
		}
		const T f0_pm1 = std::pow(value, p - 1);
		return compose(f0_pm1.multiply(value), f0_pm1.multiply(p));
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> pow(const int& n)
	{
		if (n == 0)
		{
			return get_field().get_one();
		}
		const T f0_nm1 = std::pow(value, n - 1);
		return compose(f0_nm1.multiply(value), f0_nm1.multiply(n));
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> pow(const Field_Gradient<T>& e)
	{
		return log().multiply(e).exp();
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> exp()
	{
		const T exp = std::exp(value);
		return compose(exp, exp);
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> expm1()
	{
		const T exp = std::exp(value);
		const T exp_m1 = std::expm1(value);
		return compose(exp_m1, exp);
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> log()
	{
		return compose(std::log(value), value.reciprocal());
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> log1p()
	{
		return compose(std::log1p(value), value.add(1).reciprocal());
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> log10()
	{
		return compose(std::log10(value), value.multiply(std::log(10.0)).reciprocal());
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> cos()
	{
		const Field_Sin_Cos<T> sin_cos = Sin_Cos(value);
		return compose(sin_cos.cos(), sin_cos.sin().negate());
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> sin()
	{
		const Field_Sin_Cos<T> sin_cos = Sin_Cos(value);
		return compose(sin_cos.sin(), sin_cos.cos());
	}

	/** {@inherit_doc} */
	//override
	Field_Sin_Cos<Field_Gradient<T>> sin_cos()
	{
		const Field_Sin_Cos<T> sin_cos = Sin_Cos(value);
		const Field_Gradient<T> sin = new_instance(sin_cos.sin());
		const Field_Gradient<T> cos = new_instance(sin_cos.cos());
		const T m_sin = sin_cos.sin().negate();
		for (int i{}; i < my_grad.size(); ++i)
		{
			sin.grad[i] = my_grad[i].multiply(sin_cos.cos());
			cos.grad[i] = my_grad[i].multiply(m_sin);
		}
		return Field_Sin_Cos<>(sin, cos);
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> tan()
	{
		const T tan = std::tan(value);
		return compose(tan, tan.multiply(tan).add(1));
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> acos()
	{
		return compose(std::acos(value), value.multiply(value).negate().add(1).sqrt().reciprocal().negate());
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> asin()
	{
		return compose(std::asin(value), value.multiply(value).negate().add(1).sqrt().reciprocal());
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> atan()
	{
		return compose(std::atan(value), value.multiply(value).add(1).reciprocal());
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> atan2(const Field_Gradient<T>& x)
	{
		const T inv = value.multiply(value).add(x.value.multiply(x.value)).reciprocal();
		const Field_Gradient<T> result = new_instance(std::atan2(value, x.value));
		const T x_value_inv = x.value.multiply(inv);
		const T m_value_inv = value.negate().multiply(inv);
		for (int i{}; i < my_grad.size(); ++i)
		{
			result.grad[i] = x_value_inv.multiply(grad[i]).add(x.grad[i].multiply(m_value_inv));
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> cosh()
	{
		return compose(std::cosh(value), std::sinh(value));
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> sinh()
	{
		return compose(std::sinh(value), std::cosh(value));
	}

	/** {@inherit_doc} */
	//override
	Field_Sinh_Cosh<Field_Gradient<T>> sinh_cosh()
	{
		const Field_Sinh_Cosh<T> sinh_cosh = std::sinh_cosh(value);
		const Field_Gradient<T> sinh = new_instance(sinh_cosh.sinh());
		const Field_Gradient<T> cosh = new_instance(sinh_cosh.cosh());
		for (int i{}; i < my_grad.size(); ++i)
		{
			sinh.grad[i] = my_grad[i].multiply(sinh_cosh.cosh());
			cosh.grad[i] = my_grad[i].multiply(sinh_cosh.sinh());
		}
		return Field_Sinh_Cosh<>(sinh, cosh);
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> tanh()
	{
		const T tanh = std::tanh(value);
		return compose(tanh, tanh.multiply(tanh).negate().add(1));
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> acosh()
	{
		return compose(std::acosh(value), value.multiply(value).subtract(1).sqrt().reciprocal());
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> asinh()
	{
		return compose(std::asinh(value), value.multiply(value).add(1).sqrt().reciprocal());
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> atanh()
	{
		return compose(std::atanh(value), value.multiply(value).negate().add(1).reciprocal());
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> to_degrees()
	{
		const Field_Gradient<T> result = new_instance(FastMath.to_degrees(value));
		for (int i{}; i < my_grad.size(); ++i)
		{
			result.grad[i] = FastMath.to_degrees(grad[i]);
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Gradient<T> to_radians()
	{
		const Field_Gradient<T> result = new_instance(FastMath.to_radians(value));
		for (int i{}; i < my_grad.size(); ++i)
		{
			result.grad[i] = FastMath.to_radians(grad[i]);
		}
		return result;
	}

	/** Evaluate Taylor expansion of a my_gradient.
	 * @param delta parameters offsets (&Delta;x, &Delta;y, ...)
	 * @return value of the Taylor expansion at x + &Delta;x, y + &Delta;y, ...
	 */
	T taylor(const double... delta)
	{
		T result = value;
		for (int i{}; i < my_grad.size(); ++i)
		{
			result = result.add(grad[i].multiply(delta[i]));
		}
		return result;
	}

	/** Evaluate Taylor expansion of a my_gradient.
	 * @param delta parameters offsets (&Delta;x, &Delta;y, ...)
	 * @return value of the Taylor expansion at x + &Delta;x, y + &Delta;y, ...
	 */
	T taylor(////@Suppress_Warnings("unchecked") const T... delta)
		{
			T result = value;
			for (int i{}; i < my_grad.size(); ++i)
			{
				result = result.add(grad[i].multiply(delta[i]));
			}
			return result;
		}

		/** {@inherit_doc} */
		//override
			Field_Gradient<T> linear_combination(const std::vector<Field_Gradient<T>>& a, const std::vector<Field_Gradient<T>>& b)
		{
			// extract values and first derivatives
			const Field<T> field = a[0].value.get_field();
			const int n = a.size();
			const std::vector<T> a0 = Math_Arrays::build_array(field, n);
			const std::vector<T> b0 = Math_Arrays::build_array(field, n);
			const std::vector<T> a1 = Math_Arrays::build_array(field, 2 * n);
			const std::vector<T> b1 = Math_Arrays::build_array(field, 2 * n);
			for (int i{}; i < n; ++i)
			{
				const Field_Gradient<T>& ai = a[i];
				const Field_Gradient<T>& bi = b[i];
				a0[i] = ai.value;
				b0[i] = bi.value;
				a1[2 * i] = ai.value;
				b1[2 * i + 1] = bi.value;
			}

			const Field_Gradient<T> result = new_instance(a[0].value.linear_combination(a0, b0));
			for (int k{}; k < my_grad.size(); ++k)
			{
				for (int i{}; i < n; ++i)
				{
					a1[2 * i + 1] = a[i].grad[k];
					b1[2 * i] = b[i].grad[k];
				}
				result.grad[k] = a[0].value.linear_combination(a1, b1);
			}
			return result;
		}

		/**
		 * Compute a linear combination.
		 * @param a Factors.
		 * @param b Factors.
		 * @return <code>&Sigma;<sub>i</sub> a<sub>i</sub> b<sub>i</sub></code>.
		 * @ if arrays dimensions don't match
		 */
		Field_Gradient<T> linear_combination(const std::vector<T> a, const std::vector<Field_Gradient<T>>& b)
		{
			// extract values and first derivatives
			const Field<T> field = b[0].value.get_field();
			const int      n = b.size();
			const std::vector<T> b0 = Math_Arrays::build_array(field, n);
			const std::vector<T> b1 = Math_Arrays::build_array(field, n);
			for (int i{}; i < n; ++i)
			{
				b0[i] = b[i].value;
			}

			const Field_Gradient<T> result = new_instance(b[0].value.linear_combination(a, b0));
			for (int k{}; k < my_grad.size(); ++k)
			{
				for (int i{}; i < n; ++i)
				{
					b1[i] = b[i].grad[k];
				}
				result.grad[k] = b[0].value.linear_combination(a, b1);
			}
			return result;
		}

		/** {@inherit_doc} */
		//override
		Field_Gradient<T> linear_combination(const std::vector<double>& a, const std::vector<Field_Gradient<T>>& b)
		{
			// extract values and first derivatives
			const Field<T> field = b[0].value.get_field();
			const auto n = b.size();
			const std::vector<T> b0 = Math_Arrays::build_array(field, n);
			const std::vector<T> b1 = Math_Arrays::build_array(field, n);
			for (int i{}; i < n; ++i)
			{
				b0[i] = b[i].value;
			}

			const Field_Gradient<T> result = new_instance(b[0].value.linear_combination(a, b0));
			for (int k{}; k < my_grad.size(); ++k)
			{
				for (int i{}; i < n; ++i)
				{
					b1[i] = b[i].grad[k];
				}
				result.grad[k] = b[0].value.linear_combination(a, b1);
			}
			return result;
		}

		/** {@inherit_doc} */
		//override
		Field_Gradient<T> linear_combination(const Field_Gradient<T>& a1, const Field_Gradient<T>& b1, const Field_Gradient<T>& a2, const Field_Gradient<T>& b2)
		{
			const Field_Gradient<T> result = new_instance(a1.value.linear_combination(a1.value, b1.value, a2.value, b2.value));
			for (int i{}; i < b1.grad.size(); ++i)
			{
				result.grad[i] = a1.value.linear_combination(
					a1.value,
					b1.grad[i],
					a1.grad[i],
					b1.value,
					a2.value,
					b2.grad[i],
					a2.grad[i],
					b2.value
				);
			}
			return result;
		}

		/** {@inherit_doc} */
		//override
		Field_Gradient<T> linear_combination(const double& a1, const Field_Gradient<T>& b1, const double& a2, const Field_Gradient<T>& b2)
		{
			const Field_Gradient<T> result = new_instance(b1.value.linear_combination(a1, b1.value, a2, b2.value));
			for (int i{}; i < b1.grad.size(); ++i)
			{
				result.grad[i] = b1.value.linear_combination(a1, b1.grad[i], a2, b2.grad[i]);
			}
			return result;
		}

		/** {@inherit_doc} */
		//override
		Field_Gradient<T> linear_combination(const Field_Gradient<T>& a1, const Field_Gradient<T>& b1, const Field_Gradient<T>& a2, const Field_Gradient<T>& b2, const Field_Gradient<T>& a3, const Field_Gradient<T>& b3)
		{
			const Field<T> field = a1.value.get_field();
			const std::vector<T> a = Math_Arrays::build_array(field, 6);
			const std::vector<T> b = Math_Arrays::build_array(field, 6);
			a[0] = a1.value;
			a[2] = a2.value;
			a[4] = a3.value;
			b[1] = b1.value;
			b[3] = b2.value;
			b[5] = b3.value;
			const Field_Gradient<T> result = new_instance(a1.value.linear_combination(a1.value, b1.value, a2.value, b2.value, a3.value, b3.value));
			for (int i{}; i < b1.grad.size(); ++i)
			{
				a[1] = a1.grad[i];
				a[3] = a2.grad[i];
				a[5] = a3.grad[i];
				b[0] = b1.grad[i];
				b[2] = b2.grad[i];
				b[4] = b3.grad[i];
				result.grad[i] = a1.value.linear_combination(a, b);
			}
			return result;
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
		 * @see #linear_combination(double, Field_Gradient, double, Field_Gradient)
		 * @see #linear_combination(double, Field_Gradient, double, Field_Gradient, double, Field_Gradient, double, Field_Gradient)
		 * @exception  if number of free parameters or orders are inconsistent
		 */
		Field_Gradient<T> linear_combination(const T a1, const Field_Gradient<T>& b1, const T a2, const Field_Gradient<T>& b2, const T a3, const Field_Gradient<T>& b3)
		{
			const Field_Gradient<T> result = new_instance(b1.value.linear_combination(a1, b1.value, a2, b2.value, a3, b3.value));
			for (int i{}; i < b1.grad.size(); ++i)
			{
				result.grad[i] = b1.value.linear_combination(a1, b1.grad[i], a2, b2.grad[i], a3, b3.grad[i]);
			}
			return result;
		}

		/** {@inherit_doc} */
		//override
		Field_Gradient<T> linear_combination(const double& a1, const Field_Gradient<T>& b1, const double& a2, const Field_Gradient<T>& b2, const double& a3, const Field_Gradient<T>& b3)
		{
			const Field_Gradient<T> result = new_instance(b1.value.linear_combination(a1, b1.value, a2, b2.value, a3, b3.value));
			for (int i{}; i < b1.grad.size(); ++i)
			{
				result.grad[i] = b1.value.linear_combination(a1, b1.grad[i], a2, b2.grad[i], a3, b3.grad[i]);
			}
			return result;
		}

		/** {@inherit_doc} */
		//override
		Field_Gradient<T> linear_combination(const Field_Gradient<T>& a1, const Field_Gradient<T>& b1, const Field_Gradient<T>& a2, const Field_Gradient<T>& b2, const Field_Gradient<T>& a3, const Field_Gradient<T>& b3, const Field_Gradient<T>& a4, const Field_Gradient<T>& b4)
		{
			const Field<T> field = a1.value.get_field();
			const std::vector<T> a = Math_Arrays::build_array(field, 8);
			const std::vector<T> b = Math_Arrays::build_array(field, 8);
			a[0] = a1.value;
			a[2] = a2.value;
			a[4] = a3.value;
			a[6] = a4.value;
			b[1] = b1.value;
			b[3] = b2.value;
			b[5] = b3.value;
			b[7] = b4.value;
			const Field_Gradient<T> result = new_instance(a1.value.linear_combination(a1.value, b1.value, a2.value, b2.value, a3.value, b3.value, a4.value, b4.value));
			for (int i{}; i < b1.grad.size(); ++i)
			{
				a[1] = a1.grad[i];
				a[3] = a2.grad[i];
				a[5] = a3.grad[i];
				a[7] = a4.grad[i];
				b[0] = b1.grad[i];
				b[2] = b2.grad[i];
				b[4] = b3.grad[i];
				b[6] = b4.grad[i];
				result.grad[i] = a1.value.linear_combination(a, b);
			}
			return result;
		}

		/** {@inherit_doc} */
		//override
		Field_Gradient<T> linear_combination(const double& a1, const Field_Gradient<T>& b1, const double& a2, const Field_Gradient<T>& b2, const double& a3, const Field_Gradient<T>& b3, const double& a4, const Field_Gradient<T>& b4)
		{
			const Field_Gradient<T> result = new_instance(b1.value.linear_combination(a1, b1.value, a2, b2.value, a3, b3.value, a4, b4.value));
			for (int i{}; i < b1.grad.size(); ++i)
			{
				result.grad[i] = b1.value.linear_combination(a1, b1.grad[i], a2, b2.grad[i], a3, b3.grad[i], a4, b4.grad[i]);
			}
			return result;
		}

		/** {@inherit_doc} */
		//override
		Field_Gradient<T> get_pi()
		{
			return Field_Gradient<>(get_value_field().get_zero().get_pi(), get_free_parameters());
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

			if (dynamic_cast<const Field_Gradient*>(*other) != nullptr)
			{
				//@Suppress_Warnings("unchecked")
				const Field_Gradient<T> rhs = (Field_Gradient<T>) other;
				if (!value.equals(rhs.value) || my_grad.size() != rhs.grad.size())
				{
					return false;
				}
				for (int i{}; i < my_grad.size(); ++i)
				{
					if (!grad[i].equals(rhs.grad[i]))
					{
						return false;
					}
				}
				return true;
			}

			return false;
		}

		/** Get a hash_code for the univariate derivative.
		 * @return a hash code value for this object
		 */
		 //override
		int hash_code() const
		{
			return 129 + 7 * my_value.hash_code() - 15 * Arrays.hash_code(my_grad);
		}
};