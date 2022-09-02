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
  //package org.hipparchus.analysis.differentiation;

  //import java.io.Serializable;

  //import org.hipparchus.Field;
  //import org.hipparchus.exception.;
  //import org.hipparchus.exception.Math_Runtime_Exception;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Field_Sin_Cos;
  //import org.hipparchus.util.Field_Sinh_Cosh;
  //import org.hipparchus.util.Math_Arrays;
  //import org.hipparchus.util.Math_Utils;
#include "../../util/FieldSinCos.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include "../../util/MathUtils.h"
#include "../../util/MathArrays.h"
#include "Derivative.h"

/** Class representing both the value and the differentials of a function.
 * <p>This class is the workhorse of the differentiation //package.</p>
 * <p>This class is an implementation of the extension to Rall's
 * numbers described in Dan Kalman's paper <a
 * href="http://www.dankalman.net/AUhome/pdffiles/mmgautodiff.pdf">Doubly
 * Recursive Multivariate Automatic Differentiation</a>, Mathematics Magazine, vol. 75, * no. 3, June 2002. Rall's numbers are an extension to the real numbers used
 * throughout mathematical expressions; they hold the derivative together with the
 * value of a function. Dan Kalman's derivative structures hold all partial derivatives
 * up to any specified order, with respect to any number of free parameters. Rall's
 * numbers therefore can be seen as derivative structures for order one derivative and
 * one free parameter, and real numbers can be seen as derivative structures with zero
 * order derivative and no free parameters.</p>
 * <p>{@link Derivative_Structure} instances can be used directly thanks to
 * the arithmetic operators to the mathematical functions provided as
 * methods by this class (+, -, *, /, %, sin, cos ...).</p>
 * <p>Implementing complex expressions by hand using these classes is
 * a tedious and error-prone task but has the advantage of having no limitation
 * on the derivation order despite not requiring users to compute the derivatives by
 * themselves. Implementing complex expression can also be done by developing computation
 * code using standard primitive double values and to use {@link
 * Univariate_Function_differentiator differentiators} to create the {@link
 * Derivative_Structure}-based instances. This method is simpler but may be limited in
 * the accuracy and derivation orders and may be computationally intensive (this is
 * typically the case for {@link Finite_Differences_Differentiator finite differences
 * differentiator}.</p>
 * <p>Instances of this class are guaranteed to be immutable.</p>
 * @see DS_Compiler
 * @see Field_Derivative_Structure
 */
class Derivative_Structure : Derivative<Derivative_Structure>
{
private:
	/** Factory that built the instance. */
	const DS_Factory my_factory;

	/** Combined array holding all values. */
	std::vector<double> my_data;

public:
	/** Build an instance with all values and derivatives set to 0.
	 * @param factory factory that built the instance
	 * @param data combined array holding all values
	 */
	Derivative_Structure(const DS_Factory& factory, const std::vector<double>& data) : my_factory{ factory }, my_data{ data } {};

	/** Build an instance with all values and derivatives set to 0.
	 * @param factory factory that built the instance
	 * @since 1.4
	 */
	Derivative_Structure(const DS_Factory factory) : my_factory{ factory }, my_data{ std::vector<double>(factory.get_compiler().get_size()) } {};

	/** {@inherit_doc} */
	//override
	Derivative_Structure new_instance(const double& value)
	{
		return my_factory.constant(value);
	}

	/** Get the factory that built the instance.
	 * @return factory that built the instance
	 */
	DS_Factory get_factory() const
	{
		return my_factory;
	}

	//override
	/** {@inherit_doc} */
	int get_free_parameters()
	{
		return get_factory().get_compiler().get_free_parameters();
	}

	//override
	/** {@inherit_doc} */
	int get_order()
	{
		return get_factory().get_compiler().get_order();
	}

	/** Set a derivative component.
	 * <p>
	 * This method is //package-private (no modifier specified), as it is intended
	 * to be used only by {@link DS_Factory} since it relied on the ordering of
	 * derivatives within the class. This allows avoiding checks on the index, * for performance reasons.
	 * </p>
	 * @param index index of the derivative
	 * @param value of the derivative to set
	 * @since 1.4
	 */
	void set_derivative_component(const int& index, const double& value)
	{
		my_data[index] = value;
	}

	/** {@inherit_doc}
	 */
	 //override
	double get_real() const
	{
		return my_data[0];
	}

	/** Get the value part of the derivative structure.
	 * @return value part of the derivative structure
	 * @see #get_partial_derivative(int...)
	 */
	 //override
	double get_value() const
	{
		return my_data[0];
	}

	/** {@inherit_doc} */
	//override
	double get_partial_derivative(const int ... orders)
	{
		return my_data[get_factory().get_compiler().get_partial_derivative_index(orders)];
	}

	/** Get all partial derivatives.
	 * @return a fresh copy of partial derivatives, in an array sorted according to
	 * {@link DS_Compiler#get_partial_derivative_index(int...)}
	 */
	std::vector<double> get_all_derivatives() const
	{
		return my_data;
	}

	/** {@inherit_doc}
	 */
	 //override
	Derivative_Structure add(const double& a)
	{
		const Derivative_Structure ds = my_factory.build();
		System.arraycopy(my_data, 0, ds.data, 0, my_data.size());
		ds.data[0] += a;
		return ds;
	}

	/** {@inherit_doc}
	 * @exception  if number of free parameters
	 * or orders do not match
	 */
	 //override
	Derivative_Structure add(const Derivative_Structure& a)
	{
		my_factory.check_compatibility(a.factory);
		const Derivative_Structure ds = my_factory.build();
		my_factory.get_compiler().add(my_data, 0, a.data, 0, ds.data, 0);
		return ds;
	}

	/** {@inherit_doc}
	 */
	 //override
	Derivative_Structure subtract(const double& a)
	{
		return add(-a);
	}

	/** {@inherit_doc}
	 * @exception  if number of free parameters
	 * or orders do not match
	 */
	 //override
	Derivative_Structure subtract(const Derivative_Structure& a)
	{
		my_factory.check_compatibility(a.factory);
		const Derivative_Structure ds = my_factory.build();
		my_factory.get_compiler().subtract(my_data, 0, a.data, 0, ds.data, 0);
		return ds;
	}

	/** {@inherit_doc} */
	//override
	Derivative_Structure multiply(const int& n)
	{
		return multiply(static_cast<double>(n));
	}

	/** {@inherit_doc}
	 */
	 //override
	Derivative_Structure multiply(const double& a)
	{
		const Derivative_Structure ds = my_factory.build();
		for (int i{}; i < ds.data.size(); ++i)
		{
			ds.data[i] = data[i] * a;
		}
		return ds;
	}

	/** {@inherit_doc}
	 * @exception  if number of free parameters
	 * or orders do not match
	 */
	 //override
	Derivative_Structure multiply(const Derivative_Structure& a)
	{
		my_factory.check_compatibility(a.factory);
		const Derivative_Structure result = my_factory.build();
		my_factory.get_compiler().multiply(my_data, 0, a.data, 0, result.data, 0);
		return result;
	}

	/** {@inherit_doc}
	 */
	 //override
	Derivative_Structure divide(const double& a)
	{
		const Derivative_Structure ds = my_factory.build();
		const double inv = 1.0 / a;
		for (int i{}; i < ds.data.size(); ++i)
		{
			ds.data[i] = data[i] * inv;
		}
		return ds;
	}

	/** {@inherit_doc}
	 * @exception  if number of free parameters
	 * or orders do not match
	 */
	 //override
	Derivative_Structure divide(const Derivative_Structure& a)
	{
		my_factory.check_compatibility(a.factory);
		const Derivative_Structure result = my_factory.build();
		my_factory.get_compiler().divide(my_data, 0, a.data, 0, result.data, 0);
		return result;
	}

	/** {@inherit_doc} */
	//override
	Derivative_Structure remainder(const double& a)
	{
		const Derivative_Structure ds = my_factory.build();
		System.arraycopy(my_data, 0, ds.data, 0, my_data.size());
		ds.data[0] = std::remainder(ds.data[0], a);
		return ds;
	}

	/** {@inherit_doc}
	 * @exception  if number of free parameters
	 * or orders do not match
	 */
	 //override
	Derivative_Structure remainder(const Derivative_Structure& a)
	{
		my_factory.check_compatibility(a.factory);
		const Derivative_Structure result = my_factory.build();
		my_factory.get_compiler().remainder(my_data, 0, a.data, 0, result.data, 0);
		return result;
	}

	/** {@inherit_doc} */
	//override
	Derivative_Structure negate()
	{
		const Derivative_Structure ds = my_factory.build();
		for (int i{}; i < ds.data.size(); ++i)
		{
			ds.data[i] = -my_data[i];
		}
		return ds;
	}

	/** {@inherit_doc}
	 */
	 //override
	Derivative_Structure abs()
	{
		if (Double.double_to_long_bits(my_data[0]) < 0)
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
	Derivative_Structure ceil()
	{
		return my_factory.constant(std::ceil(my_data[0]));
	}

	/** {@inherit_doc}
	 */
	 //override
	Derivative_Structure floor()
	{
		return my_factory.constant(std::floor(my_data[0]));
	}

	/** {@inherit_doc}
	 */
	 //override
	Derivative_Structure rint()
	{
		return my_factory.constant(std::rint(my_data[0]));
	}

	/** {@inherit_doc}
	 */
	 //override
	Derivative_Structure sign()
	{
		return my_factory.constant(std::signum(my_data[0]));
	}

	/** {@inherit_doc}
	 */
	 //override
	Derivative_Structure copy_sign(const Derivative_Structure& sign)
	{
		long m = Double.double_to_long_bits(my_data[0]);
		long s = Double.double_to_long_bits(sign.data[0]);
		if ((m >= 0 && s >= 0) || (m < 0 && s < 0)) { // Sign is currently OK
			return this;
		}
		return negate(); // flip sign
	}

	/** {@inherit_doc}
	 */
	 //override
	Derivative_Structure copy_sign(const double sign)
	{
		long m = Double.double_to_long_bits(my_data[0]);
		long s = Double.double_to_long_bits(sign);
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
	int get_exponent()
	{
		return FastMath.get_exponent(my_data[0]);
	}

	/** {@inherit_doc}
	 */
	 //override
	Derivative_Structure scalb(const int& n)
	{
		const Derivative_Structure ds = my_factory.build();
		for (int i{}; i < ds.data.size(); ++i)
		{
			ds.data[i] = std::scalbn(my_data[i], n);
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
	Derivative_Structure ulp()
	{
		const Derivative_Structure ds = my_factory.build();
		ds.data[0] = FastMath.ulp(my_data[0]);
		return ds;
	}

	/** {@inherit_doc}
	 * @exception  if number of free parameters
	 * or orders do not match
	 */
	 //override
	Derivative_Structure hypot(const Derivative_Structure& y)
	{
		my_factory.check_compatibility(y.factory);

		if (std::isinf(my_data[0]) || std::isinf(y.data[0]))
		{
			return my_factory.constant(INFINITY);
		}
		if (std::isnan(my_data[0]) || std::isnan(y.data[0]))
		{
			return my_factory.constant(NAN);
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
		const Derivative_Structure scaled_x = scalb(-middle_exp);
		const Derivative_Structure scaled_y = y.scalb(-middle_exp);

		// compute scaled hypotenuse
		const Derivative_Structure scaled_h = scaled_x.multiply(scaled_x).add(scaled_y.multiply(scaled_y)).sqrt();

		// remove scaling
		return scaled_h.scalb(middle_exp);
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
	 */
	static Derivative_Structure hypot(const Derivative_Structure& x, const Derivative_Structure& y)
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
	 //override
	Derivative_Structure compose(const double ... f)
	{
		Math_Utils::check_dimension(f.size(), get_order() + 1);
		const Derivative_Structure result = my_factory.build();
		my_factory.get_compiler().compose(my_data, 0, f, result.data, 0);
		return result;
	}

	/** {@inherit_doc} */
	//override
	Derivative_Structure reciprocal()
	{
		const Derivative_Structure result = my_factory.build();
		my_factory.get_compiler().pow(my_data, 0, -1, result.data, 0);
		return result;
	}

	/** {@inherit_doc}
	 */
	 //override
	Derivative_Structure sqrt()
	{
		return root_n(2);
	}

	/** {@inherit_doc}
	 */
	 //override
	Derivative_Structure cbrt()
	{
		return root_n(3);
	}

	/** {@inherit_doc}
	 */
	 //override
	Derivative_Structure root_n(const int& n)
	{
		const Derivative_Structure result = my_factory.build();
		my_factory.get_compiler().root_n(my_data, 0, n, result.data, 0);
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field<Derivative_Structure> get_field()
	{
		return my_factory.get_derivative_field();
	}

	/** Compute a<sup>x</sup> where a is a double and x a {@link Derivative_Structure}
	 * @param a number to exponentiate
	 * @param x power to apply
	 * @return a<sup>x</sup>
	 */
	static Derivative_Structure pow(const double& a, const Derivative_Structure x)
	{
		const Derivative_Structure result = x.factory.build();
		x.factory.get_compiler().pow(a, x.data, 0, result.data, 0);
		return result;
	}

	/** {@inherit_doc}
	 */
	 //override
	Derivative_Structure pow(const double& p)
	{
		const Derivative_Structure result = my_factory.build();
		my_factory.get_compiler().pow(my_data, 0, p, result.data, 0);
		return result;
	}

	/** {@inherit_doc}
	 */
	 //override
	Derivative_Structure pow(const int& n)
	{
		const Derivative_Structure result = my_factory.build();
		my_factory.get_compiler().pow(my_data, 0, n, result.data, 0);
		return result;
	}

	/** {@inherit_doc}
	 * @exception  if number of free parameters
	 * or orders do not match
	 */
	 //override
	Derivative_Structure pow(const Derivative_Structure e)

	{
		my_factory.check_compatibility(e.factory);
		const Derivative_Structure result = my_factory.build();
		my_factory.get_compiler().pow(my_data, 0, e.data, 0, result.data, 0);
		return result;
	}

	/** {@inherit_doc}
	 */
	 //override
	Derivative_Structure exp()
	{
		const Derivative_Structure result = my_factory.build();
		my_factory.get_compiler().exp(my_data, 0, result.data, 0);
		return result;
	}

	/** {@inherit_doc}
	 */
	 //override
	Derivative_Structure expm1()
	{
		const Derivative_Structure result = my_factory.build();
		my_factory.get_compiler().expm1(my_data, 0, result.data, 0);
		return result;
	}

	/** {@inherit_doc}
	 */
	 //override
	Derivative_Structure log()
	{
		const Derivative_Structure result = my_factory.build();
		my_factory.get_compiler().log(my_data, 0, result.data, 0);
		return result;
	}

	/** {@inherit_doc}
	 */
	 //override
	Derivative_Structure log1p()
	{
		const Derivative_Structure result = my_factory.build();
		my_factory.get_compiler().log1p(my_data, 0, result.data, 0);
		return result;
	}

	/** Base 10 logarithm.
	 * @return base 10 logarithm of the instance
	 */
	 //override
	Derivative_Structure log10()
	{
		const Derivative_Structure result = my_factory.build();
		my_factory.get_compiler().log10(my_data, 0, result.data, 0);
		return result;
	}

	/** {@inherit_doc}
	 */
	 //override
	Derivative_Structure cos()
	{
		const Derivative_Structure result = my_factory.build();
		my_factory.get_compiler().cos(my_data, 0, result.data, 0);
		return result;
	}

	/** {@inherit_doc}
	 */
	 //override
	Derivative_Structure sin()
	{
		const Derivative_Structure result = my_factory.build();
		my_factory.get_compiler().sin(my_data, 0, result.data, 0);
		return result;
	}

	/** {@inherit_doc}
	 */
	 //override
	Field_Sin_Cos<Derivative_Structure> sin_cos()
	{
		const Derivative_Structure sin = my_factory.build();
		const Derivative_Structure cos = my_factory.build();
		my_factory.get_compiler().sin_cos(my_data, 0, sin.data, 0, cos.data, 0);
		return Field_Sin_Cos<>(sin, cos);
	}

	/** {@inherit_doc}
	 */
	 //override
	Derivative_Structure tan()
	{
		const Derivative_Structure result = my_factory.build();
		my_factory.get_compiler().tan(my_data, 0, result.data, 0);
		return result;
	}

	/** {@inherit_doc}
	 */
	 //override
	Derivative_Structure acos()
	{
		const Derivative_Structure result = my_factory.build();
		my_factory.get_compiler().acos(my_data, 0, result.data, 0);
		return result;
	}

	/** {@inherit_doc}
	 */
	 //override
	Derivative_Structure asin()
	{
		const Derivative_Structure result = my_factory.build();
		my_factory.get_compiler().asin(my_data, 0, result.data, 0);
		return result;
	}

	/** {@inherit_doc}
	 */
	 //override
	Derivative_Structure atan()
	{
		const Derivative_Structure result = my_factory.build();
		my_factory.get_compiler().atan(my_data, 0, result.data, 0);
		return result;
	}

	/** {@inherit_doc}
	 */
	 //override
	Derivative_Structure atan2(const Derivative_Structure x)

	{
		my_factory.check_compatibility(x.factory);
		const Derivative_Structure result = my_factory.build();
		my_factory.get_compiler().atan2(my_data, 0, x.data, 0, result.data, 0);
		return result;
	}

	/** Two arguments arc tangent operation.
	 * @param y first argument of the arc tangent
	 * @param x second argument of the arc tangent
	 * @return atan2(y, x)
	 * @exception  if number of free parameters
	 * or orders do not match
	 */
	static Derivative_Structure atan2(const Derivative_Structure y, const Derivative_Structure x)

	{
		return y.atan2(x);
	}

	/** {@inherit_doc}
	 */
	 //override
	Derivative_Structure cosh()
	{
		const Derivative_Structure result = my_factory.build();
		my_factory.get_compiler().cosh(my_data, 0, result.data, 0);
		return result;
	}

	/** {@inherit_doc}
	 */
	 //override
	Derivative_Structure sinh()
	{
		const Derivative_Structure result = my_factory.build();
		my_factory.get_compiler().sinh(my_data, 0, result.data, 0);
		return result;
	}

	/** {@inherit_doc}
	 */
	 //override
	Field_Sinh_Cosh<Derivative_Structure> sinh_cosh()
	{
		const Derivative_Structure sinh = my_factory.build();
		const Derivative_Structure cosh = my_factory.build();
		my_factory.get_compiler().sinh_cosh(my_data, 0, sinh.data, 0, cosh.data, 0);
		return Field_Sinh_Cosh<>(sinh, cosh);
	}

	/** {@inherit_doc}
	 */
	 //override
	Derivative_Structure tanh()
	{
		const Derivative_Structure result = my_factory.build();
		my_factory.get_compiler().tanh(my_data, 0, result.data, 0);
		return result;
	}

	/** {@inherit_doc}
	 */
	 //override
	Derivative_Structure acosh()
	{
		const Derivative_Structure result = my_factory.build();
		my_factory.get_compiler().acosh(my_data, 0, result.data, 0);
		return result;
	}

	/** {@inherit_doc}
	 */
	 //override
	Derivative_Structure asinh()
	{
		const Derivative_Structure result = my_factory.build();
		my_factory.get_compiler().asinh(my_data, 0, result.data, 0);
		return result;
	}

	/** {@inherit_doc}
	 */
	 //override
	Derivative_Structure atanh()
	{
		const Derivative_Structure result = my_factory.build();
		my_factory.get_compiler().atanh(my_data, 0, result.data, 0);
		return result;
	}

	/** {@inherit_doc} */
	//override
	Derivative_Structure to_degrees()
	{
		const Derivative_Structure ds = my_factory.build();
		for (int i{}; i < ds.data.size(); ++i)
		{
			ds.data[i] = FastMath.to_degrees(data[i]);
		}
		return ds;
	}

	/** {@inherit_doc} */
	//override
	Derivative_Structure to_radians()
	{
		const Derivative_Structure ds = my_factory.build();
		for (int i{}; i < ds.data.size(); ++i)
		{
			ds.data[i] = FastMath.to_radians(data[i]);
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
	Derivative_Structure integrate(const int var_index, const int integration_order)
	{
		// Deal first with trivial case
		if (integration_order > get_order())
		{
			return my_factory.constant(0.);
		}
		else if (integration_order == 0)
		{
			return my_factory.build(my_data);
		}

		// Call 'inverse' (not rigorously) operation if necessary
		if (integration_order < 0)
		{
			return differentiate(var_index, -integration_order);
		}

		const std::vector<double> new_data = std::vector<double>(my_data.size());
		const DS_Compiler ds_compiler = my_factory.get_compiler();
		for (int i{}; i < new_data.size(); i++)
		{
			if (my_data[i] != 0.)
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
					new_data[index] = my_data[i];
				}
			}
		}

		return my_factory.build(new_data);
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
	Derivative_Structure differentiate(const int var_index, const int differentiation_order)
	{
		// Deal first with trivial case
		if (differentiation_order > get_order())
		{
			return my_factory.constant(0.);
		}
		if (differentiation_order == 0)
		{
			return my_factory.build(my_data);
		}

		// Call 'inverse' (not rigorously) operation if necessary
		if (differentiation_order < 0)
		{
			return integrate(var_index, -differentiation_order);
		}

		auto new_data = std::vector<double>(my_data.size());
		const DS_Compiler ds_compiler = my_factory.get_compiler();
		for (int i{}; i < new_data.size(); i++)
		{
			if (my_data[i] != 0.)
			{
				std::vector<int> orders = ds_compiler.get_partial_derivative_orders(i);
				if (orders[var_index] - differentiation_order >= 0)
				{
					const int saved = orders[var_index];
					orders[var_index] -= differentiation_order;
					const int index = ds_compiler.get_partial_derivative_index(orders);
					orders[var_index] = saved;
					new_data[index] = my_data[i];
				}
			}
		}

		return my_factory.build(new_data);
	}

	/** Evaluate Taylor expansion a derivative structure.
	 * @param delta parameters offsets (&Delta;x, &Delta;y, ...)
	 * @return value of the Taylor expansion at x + &Delta;x, y + &Delta;y, ...
	 * @Math_Runtime_Exception if factorials becomes too large
	 */
	double taylor(const double ... delta)
	{
		return my_factory.get_compiler().taylor(my_data, 0, delta);
	}

	/** {@inherit_doc}
	 * @exception  if number of free parameters
	 * or orders do not match
	 */
	 //override
	Derivative_Structure linear_combination(const std::vector<Derivative_Structure>& a, const std::vector<Derivative_Structure>& b)
	{
		// compute an accurate value, taking care of cancellations
		auto a_double = std::vector<double>(a.size());
		for (int i{}; i < a.size(); ++i)
		{
			a_double[i] = a[i].get_value();
		}
		auto b_double = std::vector<double>(b.size());
		for (int i{}; i < b.size(); ++i)
		{
			b_double[i] = b[i].get_value();
		}
		const double accurate_value = Math_Arrays::linear_combination(a_double, b_double);

		// compute a simple value, with all partial derivatives
		Derivative_Structure simple_value = a[0].get_field().get_zero();
		for (int i{}; i < a.size(); ++i)
		{
			simple_value = simple_value.add(a[i].multiply(b[i]));
		}

		// create a result with accurate value and all derivatives (not necessarily as accurate as the value)
		auto all = simple_value.get_all_derivatives();
		all[0] = accurate_value;
		return my_factory.build(all);
	}

	/** {@inherit_doc}
	 * @exception  if number of free parameters
	 * or orders do not match
	 */
	 //override
	Derivative_Structure linear_combination(const std::vector<double> a, const Derivative_Structure[] b)

	{
		// compute an accurate value, taking care of cancellations
		auto b_double = std::vector<double>(b.size());
		for (int i{}; i < b.size(); ++i)
		{
			b_double[i] = b[i].get_value();
		}
		const double& accurate_value = Math_Arrays::linear_combination(a, b_double);

		// compute a simple value, with all partial derivatives
		Derivative_Structure simple_value = b[0].get_field().get_zero();
		for (int i{}; i < a.size(); ++i)
		{
			simple_value = simple_value.add(b[i].multiply(a[i]));
		}

		// create a result with accurate value and all derivatives (not necessarily as accurate as the value)
		auto all = simple_value.get_all_derivatives();
		all[0] = accurate_value;
		return my_factory.build(all);
	}

	/** {@inherit_doc}
	 * @exception  if number of free parameters
	 * or orders do not match
	 */
	 //override
	Derivative_Structure linear_combination(const Derivative_Structure a1, const Derivative_Structure b1, const Derivative_Structure a2, const Derivative_Structure b2)

	{
		// compute an accurate value, taking care of cancellations
		const double& accurate_value = Math_Arrays::linear_combination(a1.get_value(), b1.get_value(), a2.get_value(), b2.get_value());

		// compute a simple value, with all partial derivatives
		const Derivative_Structure simple_value = a1.multiply(b1).add(a2.multiply(b2));

		// create a result with accurate value and all derivatives (not necessarily as accurate as the value)
		auto all = simple_value.get_all_derivatives();
		all[0] = accurate_value;
		return my_factory.build(all);
	}

	/** {@inherit_doc}
	 * @exception  if number of free parameters
	 * or orders do not match
	 */
	 //override
	Derivative_Structure linear_combination(const double& a1, const Derivative_Structure b1, const double& a2, const Derivative_Structure b2)

	{
		my_factory.check_compatibility(b1.factory);
		my_factory.check_compatibility(b2.factory);

		const Derivative_Structure ds = my_factory.build();
		my_factory.get_compiler().linear_combination(a1, b1.data, 0, a2, b2.data, 0, ds.data, 0);

		return ds;
	}

	/** {@inherit_doc}
	 * @exception  if number of free parameters
	 * or orders do not match
	 */
	 //override
	Derivative_Structure linear_combination(const Derivative_Structure a1, const Derivative_Structure b1, const Derivative_Structure a2, const Derivative_Structure b2, const Derivative_Structure a3, const Derivative_Structure b3)

	{
		// compute an accurate value, taking care of cancellations
		const double& accurate_value = Math_Arrays::linear_combination(a1.get_value(), b1.get_value(), a2.get_value(), b2.get_value(), a3.get_value(), b3.get_value());

		// compute a simple value, with all partial derivatives
		const Derivative_Structure simple_value = a1.multiply(b1).add(a2.multiply(b2)).add(a3.multiply(b3));

		// create a result with accurate value and all derivatives (not necessarily as accurate as the value)
		auto all = simple_value.get_all_derivatives();
		all[0] = accurate_value;
		return my_factory.build(all);
	}

	/** {@inherit_doc}
	 * @exception  if number of free parameters
	 * or orders do not match
	 */
	 //override
	Derivative_Structure linear_combination(const double& a1, const Derivative_Structure b1, const double& a2, const Derivative_Structure b2, const double& a3, const Derivative_Structure b3)

	{
		my_factory.check_compatibility(b1.factory);
		my_factory.check_compatibility(b2.factory);
		my_factory.check_compatibility(b3.factory);

		const Derivative_Structure ds = my_factory.build();
		my_factory.get_compiler().linear_combination(a1, b1.data, 0, a2, b2.data, 0, a3, b3.data, 0, ds.data, 0);

		return ds;
	}

	/** {@inherit_doc}
	 * @exception  if number of free parameters
	 * or orders do not match
	 */
	 //override
	Derivative_Structure linear_combination(const Derivative_Structure a1, const Derivative_Structure b1, const Derivative_Structure a2, const Derivative_Structure b2, const Derivative_Structure a3, const Derivative_Structure b3, const Derivative_Structure a4, const Derivative_Structure b4)

	{
		// compute an accurate value, taking care of cancellations
		const double& accurate_value = Math_Arrays::linear_combination(a1.get_value(), b1.get_value(), a2.get_value(), b2.get_value(), a3.get_value(), b3.get_value(), a4.get_value(), b4.get_value());

		// compute a simple value, with all partial derivatives
		const Derivative_Structure simple_value = a1.multiply(b1).add(a2.multiply(b2)).add(a3.multiply(b3)).add(a4.multiply(b4));

		// create a result with accurate value and all derivatives (not necessarily as accurate as the value)
		auto all = simple_value.get_all_derivatives();
		all[0] = accurate_value;
		return my_factory.build(all);
	}

	/** {@inherit_doc}
	 * @exception  if number of free parameters
	 * or orders do not match
	 */
	 //override
	Derivative_Structure linear_combination(const double& a1, const Derivative_Structure b1, const double& a2, const Derivative_Structure b2, const double& a3, const Derivative_Structure b3, const double& a4, const Derivative_Structure b4)

	{
		my_factory.check_compatibility(b1.factory);
		my_factory.check_compatibility(b2.factory);
		my_factory.check_compatibility(b3.factory);
		my_factory.check_compatibility(b4.factory);

		const Derivative_Structure ds = my_factory.build();
		my_factory.get_compiler().linear_combination(a1, b1.data, 0, a2, b2.data, 0, a3, b3.data, 0, a4, b4.data, 0, ds.data, 0);

		return ds;
	}

	/** {@inherit_doc}
	 */
	 //override
	Derivative_Structure get_pi()
	{
		return my_factory.get_derivative_field().get_pi();
	}

	/**
	 * Test for the equality of two derivative structures.
	 * <p>
	 * Derivative structures are considered equal if they have the same number
	 * of free parameters, the same derivation order, and the same derivatives.
	 * </p>
	 * @param other Object to test for equality to this
	 * @return true if two derivative structures are equal
	 */
	 //override
	bool equals(Object other)
	{
		if (this == other)
		{
			return true;
		}

		if (dynamic_cast<const Derivative_Structure*>(*other) != nullptr)
		{
			const Derivative_Structure rhs = (Derivative_Structure)other;
			return (get_free_parameters() == rhs.get_free_parameters()) &&
				(get_order() == rhs.get_order()) &&
				Math_Arrays::equals(my_data, rhs.data);
		}

		return false;
	}

	/**
	 * Get a hash_code for the derivative structure.
	 * @return a hash code value for this object
	 */
	 //override
	int hash_code()
	{
		return 227 + 229 * get_free_parameters() + 233 * get_order() + 239 * Math_Utils::hash(my_data);
	}
};