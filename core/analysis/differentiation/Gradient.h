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
#include "../../CalculusFieldElement.hpp"
#include "Derivative.h"
#include <vector>
#include <algorithm>
#include <numbers>
#include "../../util/MathArrays.h"
#include "../../util/SinCos.h"
//import java.io.Serializable;
//import java.util.Arrays;

//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Field_Sin_Cos;
//import org.hipparchus.util.Field_Sinh_Cosh;
//import org.hipparchus.util.Math_Arrays;
//import org.hipparchus.util.Math_Utils;
//import org.hipparchus.util.Sin_Cos;
//import org.hipparchus.util.Sinh_Cosh;

/** Class representing both the value and the differentials of a function.
 * <p>This class is a stripped-down version of {@link Derivative_Structure}
 * with {@link Derivative_Structure#get_order() derivation order} limited to one.
 * It should have less overhead than {@link Derivative_Structure} in its domain.</p>
 * <p>This class is an implementation of Rall's numbers. Rall's numbers are an
 * extension to the real numbers used throughout mathematical expressions; they hold
 * the derivative together with the value of a function.</p>
 * <p>{@link Gradient} instances can be used directly thanks to
 * the arithmetic operators to the mathematical functions provided as
 * methods by this class (+, -, *, /, %, sin, cos ...).</p>
 * <p>Implementing complex expressions by hand using these classes is
 * a tedious and error-prone task but has the advantage of having no limitation
 * on the derivation order despite not requiring users to compute the derivatives by
 * themselves.</p>
 * <p>Instances of this class are guaranteed to be immutable.</p>
 * @see Derivative_Structure
 * @see Univariate_Derivative_1
 * @see Univariate_Derivative_2
 * @see Field_Derivative_Structure
 * @see Field_Univariate_Derivative_1
 * @see Field_Univariate_Derivative_2
 * @see Field_Gradient
 * @since 1.7
 */
class Gradient : Derivative<Gradient>, Calculus_Field_Element<Gradient>
{
private:

	/** Value of the function. */
	const double my_value;

	/** Gradient of the function. */
	const std::vector<double> my_grad;

	/** Build an instance with values and unitialized derivatives array.
	 * @param value value of the function
	 * @param free_parameters number of free parameters
	 */
	Gradient(const double& value, int free_parameters) : my_value{ value }, my_grad{ std::vector<double>(free_parameters) } {};

	template<typename T>
	static int sgn(const T& val)
	{
		return (T(0) < val) - (val < T(0));
	}

public:
	/** Build an instance with values and derivative.
	 * @param value value of the function
	 * @param gradient gradient of the function
	 */
	Gradient(const double& value, const std::vector<double>& gradient)
	{
		*(my_value, gradient.size());
		System.arraycopy(gradient, 0, my_grad, 0, my_grad.size());
	}

	/** Build an instance from a {@link Derivative_Structure}.
	 * @param ds derivative structure
	 * @exception  if {@code ds} order
	 * is not 1
	 */
	Gradient(const Derivative_Structure& ds)
	{
		this(ds.get_value(), ds.get_free_parameters());
		Math_Utils::check_dimension(ds.get_order(), 1);
		System.arraycopy(ds.get_all_derivatives(), 1, my_grad, 0, my_grad.size());
	}

	/** Build an instance corresponding to a constant value.
	 * @param free_parameters number of free parameters (i.e. dimension of the gradient)
	 * @param value constant value of the function
	 * @return a {@code Gradient} with a constant value and all derivatives set to 0.0
	 */
	static Gradient constant(const int& free_parameters, const double& value)
	{
		return Gradient(value, free_parameters);
	}

	/** Build a {@code Gradient} representing a variable.
	 * <p>Instances built using this method are considered
	 * to be the free variables with respect to which differentials
	 * are computed. As such, their differential with respect to
	 * themselves is +1.</p>
	 * @param free_parameters number of free parameters (i.e. dimension of the gradient)
	 * @param index index of the variable (from 0 to {@link #get_free_parameters() get_free_parameters()} - 1)
	 * @param value value of the variable
	 * @return a {@code Gradient} with a constant value and all derivatives set to 0.0 except the
	 * one at {@code index} which will be set to 1.0
	 */
	static Gradient variable(const int& free_parameters, const int& index, const double& value)
	{
		const Gradient g = Gradient(value, free_parameters);
		g.get_gradient()[index] = 1.0;
		return g;
	}

	/** {@inherit_doc} */
	//override
	Gradient new_instance(const double& c)
	{
		return Gradient(c, std::vector<double>(my_grad.size()));
	}

	/** {@inherit_doc} */
	//override
	double get_real() const
	{
		return get_value();
	}

	/** Get the value part of the function.
	 * @return value part of the value of the function
	 */
	 //override
	double get_value() const
	{
		return my_value;
	}

	/** Get the gradient part of the function.
	 * @return gradient part of the value of the function
	 * @see #get_partial_derivativestatic_cast<int>(
	 */
	std::vector<double> get_gradient() const
	{
		return my_grad;
	}

	/** {@inherit_doc} */
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
	double get_partial_derivative(const std::vector<int>& orders)
	{
		// check the number of components
		if (orders.size() != my_grad.size())
		{
			throw std::exception("not implmented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, orders.size(), my_grad.size());
		}

		// check that either all derivation orders are set to 0, // or that only one is set to 1 and all other ones are set to 0
		int selected{ -1 };
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

		return selected < 0
			? my_value
			: my_grad[selected];
	}

	/** Get the partial derivative with respect to one parameter.
	 * @param n index of the parameter (counting from 0)
	 * @return partial derivative with respect to the n<sup>th</sup> parameter
	 * @exception  if n is either negative or larger
	 * or equal to {@link #get_free_parameters()}
	 */
	double get_partial_derivative(const int& n)
	{
		if (n < 0 || n >= my_grad.size())
		{
			throw std::exception("not implmented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::OUT_OF_RANGE_SIMPLE, n, 0, my_grad.size() - 1);
		}
		return my_grad[n];
	}

	/** Convert the instance to a {@link Derivative_Structure}.
	 * @return derivative structure with same value and derivative as the instance
	 */
	Derivative_Structure to_derivative_structure()
	{
		auto derivatives = std::vector<double>(1 + my_grad.size());
		derivatives[0] = my_value;
		System.arraycopy(my_grad, 0, derivatives, 1, my_grad.size());
		return get_field().get_conversion_factory().build(derivatives);
	}

	/** {@inherit_doc} */
	//override
	Gradient add(const double& a)
	{
		return Gradient(my_value + a, my_grad);
	}

	/** {@inherit_doc} */
	//override
	Gradient add(const Gradient& a)
	{
		auto result = new_instance(my_value + a.get_value());
		for (int i{}; i < my_grad.size(); ++i)
		{
			result.get_gradient()[i] = my_grad[i] + a.get_gradient()[i];
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Gradient subtract(const double& a) const
	{
		return Gradient(my_value - a, my_grad);
	}

	/** {@inherit_doc} */
	//override
	Gradient subtract(const Gradient& a)
	{
		auto result = new_instance(my_value - a.get_value());
		for (int i{}; i < my_grad.size(); ++i)
		{
			result.get_gradient()[i] = my_grad[i] - a.get_gradient()[i];
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Gradient multiply(const int& n)
	{
		auto result = new_instance(my_value * n);
		for (int i{}; i < my_grad.size(); ++i)
		{
			result.get_gradient()[i] = my_grad[i] * n;
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Gradient multiply(const double& a)
	{
		auto result = new_instance(my_value * a);
		for (int i{}; i < my_grad.size(); ++i)
		{
			result.get_gradient()[i] = my_grad[i] * a;
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Gradient multiply(const Gradient& a)
	{
		auto result = new_instance(my_value * a.get_value());
		for (int i{}; i < my_grad.size(); ++i)
		{
			result.get_gradient()[i] = my_grad[i] * a.get_value() + my_value * a.get_gradient()[i];
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Gradient divide(const double& a)
	{
		auto result = new_instance(my_value / a);
		for (int i{}; i < my_grad.size(); ++i)
		{
			result.get_gradient()[i] = my_grad[i] / a;
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Gradient divide(const Gradient& a)
	{
		const double inv1 = 1.0 / a.get_value();
		const double inv2 = inv1 * inv1;
		auto result = new_instance(my_value * inv1);
		for (int i{}; i < my_grad.size(); ++i)
		{
			result.get_gradient()[i] = (my_grad[i] * a.get_value() - my_value * a.get_gradient()[i]) * inv2;
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Gradient remainder(const double& a)
	{
		return Gradient(std::remainder(my_value, a), my_grad);
	}

	/** {@inherit_doc} */
	//override
	Gradient remainder(const Gradient& a)
	{
		// compute k such that lhs % rhs = lhs - k rhs
		const double rem = std::remainder(my_value, a.get_value());
		const double k = std::rint((my_value - rem) / a.get_value());

		auto result = new_instance(rem);
		for (int i{}; i < my_grad.size(); ++i)
		{
			result.get_gradient()[i] = my_grad[i] - k * a.get_gradient()[i];
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Gradient negate()
	{
		auto result = new_instance(-value);
		for (int i{}; i < my_grad.size(); ++i)
		{
			result.get_gradient()[i] = -grad[i];
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Gradient abs()
	{
		if (Double.double_to_long_bits(my_value) < 0)
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
	Gradient ceil()
	{
		return new_instance(std::ceil(my_value));
	}

	/** {@inherit_doc} */
	//override
	Gradient floor()
	{
		return new_instance(std::floor(my_value));
	}

	/** {@inherit_doc} */
	//override
	Gradient rint()
	{
		return new_instance(std::rint(my_value));
	}

	/** {@inherit_doc} */
	//override
	Gradient sign()
	{
		return new_instance(sgn(my_value));
	}

	/** {@inherit_doc} */
	//override
	Gradient copy_sign(const Gradient& sign)
	{
		long m = Double.double_to_long_bits(my_value);
		long s = Double.double_to_long_bits(sign.get_value());
		if ((m >= 0 && s >= 0) || (m < 0 && s < 0)) { // Sign is currently OK
			return this;
		}
		return negate(); // flip sign
	}

	/** {@inherit_doc} */
	//override
	Gradient copy_sign(const double sign)
	{
		long m = Double.double_to_long_bits(my_value);
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
		return FastMath.get_exponent(my_value);
	}

	/** {@inherit_doc} */
	//override
	Gradient scalb(const int& n)
	{
		const auto result = new_instance(std::scalbn(my_value, n));
		for (int i{}; i < my_grad.size(); ++i)
		{
			result.get_gradient()[i] = std::scalbn(my_grad[i], n);
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
	Gradient ulp()
	{
		return new_instance(std::ulp(my_value));
	}

	/** {@inherit_doc} */
	//override
	Gradient hypot(const Gradient& y)
	{
		if (std::isinf(my_value) || std::isinf(y.get_value()))
		{
			return new_instance(INFINITY);
		}
		if (std::isnan(my_value) || std::isnan(y.get_value()))
		{
			return new_instance(NAN);
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
		const auto scaled_x = scalb(-middle_exp);
		const auto scaled_y = y.scalb(-middle_exp);

		// compute scaled hypotenuse
		const Gradient scaled_h =
			scaled_x.multiply(scaled_x).add(scaled_y.multiply(scaled_y)).sqrt();

		// remove scaling
		return scaled_h.scalb(middle_exp);
	}

	/** {@inherit_doc} */
	//override
	Gradient reciprocal()
	{
		const double inv1 = 1.0 / my_value;
		const double inv2 = inv1 * inv1;
		auto result = new_instance(inv1);
		for (int i{}; i < my_grad.size(); ++i)
		{
			result.get_gradient()[i] = -my_grad[i] * inv2;
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Gradient compose(const std::vector<double>& f)
	{
		Math_Utils::check_dimension(f.size(), get_order() + 1);
		auto result = new_instance(f[0]);
		for (int i{}; i < my_grad.size(); ++i)
		{
			result.get_gradient()[i] = f[1] * my_grad[i];
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Gradient sqrt()
	{
		const double s = std::sqrt(my_value);
		return compose(s, 1 / (2 * s));
	}

	/** {@inherit_doc} */
	//override
	Gradient cbrt()
	{
		const auto c = std::cbrt(my_value);
		return compose(c, 1 / (3 * c * c));
	}

	/** {@inherit_doc} */
	//override
	Gradient root_n(const int& n)
	{
		if (n == 2)
		{
			return sqrt();
		}
		if (n == 3)
		{
			return cbrt();
		}

		const double r = std::pow(my_value, 1.0 / n);
		return compose(r, 1 / (n * std::pow(r, n - 1)));
	}

	/** {@inherit_doc} */
	//override
	Gradient_Field get_field()
	{
		return Gradient_Field.get_field(get_free_parameters());
	}

	/** Compute a<sup>x</sup> where a is a double and x a {@link Gradient}
	 * @param a number to exponentiate
	 * @param x power to apply
	 * @return a<sup>x</sup>
	 */
	static Gradient pow(const double& a, const Gradient& x)
	{
		if (a == 0)
		{
			return x.get_field().get_zero();
		}

		const double& a_x = std::pow(a, x.get_value());
		const double& a_xln_a = a_x * std::log(a);
		auto result = x.new_instance(a_x);
		for (int i{}; i < x.get_gradient().size(); ++i)
		{
			result.get_gradient()[i] = a_xln_a * x.get_gradient()[i];
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Gradient pow(const double& p)
	{
		if (p == 0)
		{
			return get_field().get_one();
		}

		const auto value_pm1 = std::pow(my_value, p - 1);
		return compose(value_pm1 * my_value, p * value_pm1);
	}

	/** {@inherit_doc} */
	//override
	Gradient pow(const int& n)
	{
		if (n == 0)
		{
			return get_field().get_one();
		}

		const double value_nm1 = std::pow(my_value, n - 1);
		return compose(value_nm1 * my_value, n * value_nm1);
	}

	/** {@inherit_doc} */
	//override
	Gradient pow(const Gradient& e)
	{
		return log().multiply(e).exp();
	}

	/** {@inherit_doc} */
	//override
	Gradient exp()
	{
		const double exp = std::exp(my_value);
		return compose(exp, exp);
	}

	/** {@inherit_doc} */
	//override
	Gradient expm1()
	{
		const double exp = std::exp(my_value);
		const double exp_m1 = std::expm1(my_value);
		return compose(exp_m1, exp);
	}

	/** {@inherit_doc} */
	//override
	Gradient log()
	{
		return compose(std::log(my_value), 1 / my_value);
	}

	/** {@inherit_doc} */
	//override
	Gradient log1p()
	{
		return compose(std::log1p(my_value), 1 / (1 + my_value));
	}

	/** {@inherit_doc} */
	//override
	Gradient log10()
	{
		return compose(std::log10(my_value), 1 / (my_value * std::log(10.0)));
	}

	/** {@inherit_doc} */
	//override
	Gradient cos()
	{
		const Sin_Cos sin_cos = Sin_Cos(my_value);
		return compose(sin_cos.cos(), -sin_cos.sin());
	}

	/** {@inherit_doc} */
	//override
	Gradient sin()
	{
		const Sin_Cos sin_cos = Sin_Cos(my_value);
		return compose(sin_cos.sin(), sin_cos.cos());
	}

	/** {@inherit_doc} */
	//override
	Field_Sin_Cos<Gradient> sin_cos()
	{
		const Sin_Cos sin_cos = Sin_Cos(my_value);
		auto sin = new_instance(sin_cos.sin());
		auto cos = new_instance(sin_cos.cos());
		for (int i{}; i < my_grad.size(); ++i)
		{
			sin.get_gradient()[i] = +my_grad[i] * sin_cos.cos();
			cos.get_gradient()[i] = -my_grad[i] * sin_cos.sin();
		}
		return Field_Sin_Cos<>(sin, cos);
	}

	/** {@inherit_doc} */
	//override
	Gradient tan()
	{
		const double tan = std::tan(my_value);
		return compose(tan, 1 + tan * tan);
	}

	/** {@inherit_doc} */
	//override
	Gradient acos()
	{
		return compose(std::acos(my_value), -1 / std::sqrt(1 - my_value * my_value));
	}

	/** {@inherit_doc} */
	//override
	Gradient asin()
	{
		return compose(std::asin(my_value), 1 / std::sqrt(1 - my_value * my_value));
	}

	/** {@inherit_doc} */
	//override
	Gradient atan()
	{
		return compose(std::atan(my_value), 1 / (1 + my_value * my_value));
	}

	/** {@inherit_doc} */
	//override
	Gradient atan2(const Gradient& x)
	{
		const double inv = 1.0 / (my_value * my_value + x.get_value() * x.get_value());
		auto result = new_instance(std::atan2(my_value, x.get_value()));
		for (int i{}; i < my_grad.size(); ++i)
		{
			result.get_gradient()[i] = (x.get_value() * my_grad[i] - x.get_gradient()[i] * my_value) * inv;
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Gradient cosh()
	{
		return compose(std::cosh(my_value), std::sinh(my_value));
	}

	/** {@inherit_doc} */
	//override
	Gradient sinh()
	{
		return compose(std::sinh(my_value), std::cosh(my_value));
	}

	/** {@inherit_doc} */
	//override
	Field_Sinh_Cosh<Gradient> sinh_cosh()
	{
		const Sinh_Cosh sinh_cosh = std::sinh_cosh(my_value);
		const Gradient sinh = new_instance(sinh_cosh.sinh());
		const Gradient cosh = new_instance(sinh_cosh.cosh());
		for (int i{}; i < my_grad.size(); ++i)
		{
			sinh.get_gradient()[i] = my_grad[i] * sinh_cosh.cosh();
			cosh.get_gradient()[i] = my_grad[i] * sinh_cosh.sinh();
		}
		return Field_Sinh_Cosh<>(sinh, cosh);
	}

	/** {@inherit_doc} */
	//override
	Gradient tanh()
	{
		const double tanh = std::tanh(my_value);
		return compose(tanh, 1 - tanh * tanh);
	}

	/** {@inherit_doc} */
	//override
	Gradient acosh()
	{
		return compose(std::acosh(my_value), 1 / std::sqrt(my_value * my_value - 1));
	}

	/** {@inherit_doc} */
	//override
	Gradient asinh()
	{
		return compose(std::asinh(my_value), 1 / std::sqrt(my_value * my_value + 1));
	}

	/** {@inherit_doc} */
	//override
	Gradient atanh()
	{
		return compose(std::atanh(my_value), 1 / (1 - my_value * my_value));
	}

	/** {@inherit_doc} */
	//override
	Gradient to_degrees()
	{
		auto result = new_instance(FastMath.to_degrees(my_value));
		for (int i{}; i < my_grad.size(); ++i)
		{
			result.get_gradient()[i] = FastMath.to_degrees(my_grad[i]);
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Gradient to_radians()
	{
		auto result = new_instance(FastMath.to_radians(my_value));
		for (int i{}; i < my_grad.size(); ++i)
		{
			result.get_gradient()[i] = FastMath.to_radians(my_grad[i]);
		}
		return result;
	}

	/** Evaluate Taylor expansion a derivative structure.
	 * @param delta parameters offsets (&Delta;x, &Delta;y, ...)
	 * @return value of the Taylor expansion at x + &Delta;x, y + &Delta;y, ...
	 */
	double taylor(const std::vector<double>& delta)
	{
		auto result = my_value;
		for (int i{}; i < my_grad.size(); ++i)
		{
			result += delta[i] * my_grad[i];
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Gradient linear_combination(const std::vector<Gradient>& a, const std::vector<Gradient>& b)
	{
		// extract values and first derivatives
		const auto n = a.size();
		auto a0 = std::vector<double>(n);
		auto b0 = std::vector<double>(n);
		auto a1 = std::vector<double>(2 * n);
		auto b1 = std::vector<double>(2 * n);
		for (int i{}; i < n; ++i)
		{
			const auto ai = a[i];
			const auto bi = b[i];
			a0[i] = ai.get_value();
			b0[i] = bi.get_value();
			a1[2 * i] = ai.get_value();
			b1[2 * i + 1] = bi.get_value();
		}

		auto result = new_instance(Math_Arrays::linear_combination(a0, b0));
		for (int k{}; k < my_grad.size(); ++k)
		{
			for (int i{}; i < n; ++i)
			{
				a1[2 * i + 1] = a[i].get_gradient()[k];
				b1[2 * i] = b[i].get_gradient()[k];
			}
			result.get_gradient()[k] = Math_Arrays::linear_combination(a1, b1);
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Gradient linear_combination(const std::vector<double> a, const std::vector<Gradient>& b)
	{
		// extract values and first derivatives
		const auto n = b.size();
		auto b0 = std::vector<double>(n);
		auto b1 = std::vector<double>(n);
		for (int i{}; i < n; ++i)
		{
			b0[i] = b[i].get_value();
		}

		auto result = new_instance(Math_Arrays::linear_combination(a, b0));
		for (int k{}; k < my_grad.size(); ++k)
		{
			for (int i{}; i < n; ++i)
			{
				b1[i] = b[i].get_gradient()[k];
			}
			result.get_gradient()[k] = Math_Arrays::linear_combination(a, b1);
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Gradient linear_combination(const Gradient& a1, const Gradient& b1, const Gradient& a2, const Gradient& b2)
	{
		auto result = new_instance(Math_Arrays::linear_combination(a1.get_value(), b1.get_value(), a2.get_value(), b2.get_value()));
		for (int i{}; i < b1.get_gradient().size(); ++i)
		{
			result.get_gradient()[i] = Math_Arrays::linear_combination(a1.get_value(), b1.get_gradient()[i], a1.get_gradient()[i], b1.get_value(), a2.get_value(), b2.get_gradient()[i], a2.get_gradient()[i], b2.get_value());
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Gradient linear_combination(const double& a1, const Gradient b1, const double& a2, const Gradient b2)
	{
		auto result = new_instance(Math_Arrays::linear_combination(a1, b1.get_value(), a2, b2.get_value()));
		for (int i{}; i < b1.get_gradient().size(); ++i)
		{
			result.get_gradient()[i] = Math_Arrays::linear_combination(a1, b1.get_gradient()[i], a2, b2.get_gradient()[i]);
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Gradient linear_combination(const Gradient& a1, const Gradient b1, const Gradient a2, const Gradient b2, const Gradient a3, const Gradient b3)
	{
		const std::vector<double> a =
		{
			a1.get_value(), 0, a2.get_value(), 0, a3.get_value(), 0
		};
		const std::vector<double> b =
		{
			0, b1.get_value(), 0, b2.get_value(), 0, b3.get_value()
		};
		auto result = new_instance(Math_Arrays::linear_combination(a1.get_value(), b1.get_value(), a2.get_value(), b2.get_value(), a3.get_value(), b3.get_value()));
		for (int i{}; i < b1.get_gradient().size(); ++i)
		{
			a[1] = a1.get_gradient()[i];
			a[3] = a2.get_gradient()[i];
			a[5] = a3.get_gradient()[i];
			b[0] = b1.get_gradient()[i];
			b[2] = b2.get_gradient()[i];
			b[4] = b3.get_gradient()[i];
			result.get_gradient()[i] = Math_Arrays::linear_combination(a, b);
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Gradient linear_combination(const double& a1, const Gradient b1, const double& a2, const Gradient b2, const double& a3, const Gradient b3)
	{
		auto result = new_instance(Math_Arrays::linear_combination(a1, b1.get_value(), a2, b2.get_value(), a3, b3.get_value()));
		for (int i{}; i < b1.get_gradient().size(); ++i)
		{
			result.get_gradient()[i] = Math_Arrays::linear_combination(a1, b1.get_gradient()[i], a2, b2.get_gradient()[i], a3, b3.get_gradient()[i]);
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Gradient linear_combination(const Gradient& a1, const Gradient b1, const Gradient a2, const Gradient b2, const Gradient a3, const Gradient b3, const Gradient a4, const Gradient b4)
	{
		std::vector<double> a =
		{
			a1.get_value(), 0, a2.get_value(), 0, a3.get_value(), 0, a4.get_value(), 0
		};
		std::vector<double> b =
		{
			0, b1.get_value(), 0, b2.get_value(), 0, b3.get_value(), 0, b4.get_value()
		};
		auto result = new_instance(Math_Arrays::linear_combination(a1.get_value(), b1.get_value(), a2.get_value(), b2.get_value(), a3.get_value(), b3.get_value(), a4.get_value(), b4.get_value()));
		for (int i{}; i < b1.get_gradient().size(); ++i)
		{
			a[1] = a1.get_gradient()[i];
			a[3] = a2.get_gradient()[i];
			a[5] = a3.get_gradient()[i];
			a[7] = a4.get_gradient()[i];
			b[0] = b1.get_gradient()[i];
			b[2] = b2.get_gradient()[i];
			b[4] = b3.get_gradient()[i];
			b[6] = b4.get_gradient()[i];
			result.get_gradient()[i] = Math_Arrays::linear_combination(a, b);
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Gradient linear_combination(const double& a1, const Gradient b1, const double& a2, const Gradient b2, const double& a3, const Gradient b3, const double& a4, const Gradient b4)
	{
		auto result = new_instance(Math_Arrays::linear_combination(a1, b1.get_value(), a2, b2.get_value(), a3, b3.get_value(), a4, b4.get_value()));
		for (int i{}; i < b1.get_gradient().size(); ++i)
		{
			result.get_gradient()[i] = Math_Arrays::linear_combination(a1, b1.get_gradient()[i], a2, b2.get_gradient()[i], a3, b3.get_gradient()[i], a4, b4.get_gradient()[i]);
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Gradient get_pi()
	{
		return Gradient(std::numbers::pi, get_free_parameters());
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

		if (dynamic_cast<const Gradient*>(*other) != nullptr)
		{
			const auto rhs = static_cast<Gradient>(other);
			return my_value == rhs.get_value() && Math_Arrays::equals(my_grad, rhs.grad);
		}

		return false;
	}

	/** Get a hash_code for the univariate derivative.
	 * @return a hash code value for this object
	 */
	 //override
	int hash_code()
	{
		return 129 + 7 * Double.hash_code(my_value) - 15 * Arrays.hash_code(my_grad);
	}
};