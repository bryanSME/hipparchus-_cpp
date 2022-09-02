#pragma once
/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *	  http://www.apache.org/licenses/LICENSE-2.0
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
  //package org.hipparchus.util;

  //import org.hipparchus.Calculus_Field_Element;
  //import org.hipparchus.Field;
  //import org.hipparchus.exception.;
#include <vector>
#include <numbers>
#include "../CalculusFieldElement.hpp"

/**
 * This class wraps a {@code double} value in an object. It is similar to the
 * standard class {@link Double}, while also implementing the
 * {@link Calculus_Field_Element} interface.
 */
class Decimal64 extends Number : public Calculus_Field_Element<Decimal64>, public Comparable<Decimal64>
{
private:
	/** The primitive {@code double} value of this object. */
	const double my_value;

public:
	/** The constant value of {@code 0} as a {@code Decimal64}. */
	static const Decimal64 ZERO{ 0 };

	/** The constant value of {@code 1d} as a {@code Decimal64}. */
	static const Decimal64 ONE{ 1 };

	/** The constant value of π as a {@code Decimal64}. */
	static const Decimal64 PI{ std::numbers::pi };

	/**
	 * The constant value of {@link Double#NEGATIVE_INFINITY} as a
	 * {@code Decimal64}.
	 */
	static const Decimal64 NEGATIVE_INFINITY{ - std::numeric_limits<double>::infinity() };

	/**
	 * The constant value of {@link Double#POSITIVE_INFINITY} as a
	 * {@code Decimal64}.
	 */
	static const Decimal64 POSITIVE_INFINITY{ std::numeric_limits<double>::infinity() };

	/** The constant value of {@link Double#NaN} as a {@code Decimal64}. */
	static const Decimal64{ std::numeric_limits<double>::quiet_NaN() };

	/**
	 * Creates a instance of this class.
	 *
	 * @param x the primitive {@code double} value of the object to be created
	 */
	Decimal64(const double& x) : my_value{ x }
	{
	}

	/*
	 * Methods from the Field_Element interface.
	 */

	 /** {@inherit_doc} */
	 //override
	 Decimal64 new_instance(const double& v)
	 {
		 return Decimal64(v);
	 }

	 /** {@inherit_doc} */
	 //override
	 Field<Decimal64> get_field()
	 {
		 return Decimal64_Field.get_instance();
	 }

	 /**
	  * {@inherit_doc}
	  *
	  * The current implementation strictly enforces
	  * {@code this.add(a).equals(new Decimal64(this.double_value()
	  * + a.double_value()))}.
	  */
	 //override
	 Decimal64 add(const Decimal64& a)
	 {
		return Decimal64(my_value + a.value);
	 }

	/**
	 * {@inherit_doc}
	 *
	 * The current implementation strictly enforces
	 * {@code this.subtract(a).equals(new Decimal64(this.double_value()
	 * - a.double_value()))}.
	 */
	//override
	Decimal64 subtract(const Decimal64& a)
	{
		return Decimal64(my_value - a.value);
	}

	/**
	 * {@inherit_doc}
	 *
	 * The current implementation strictly enforces
	 * {@code this.negate().equals(new Decimal64(-this.double_value()))}.
	 */
	//override
	Decimal64 negate()
	{
		return Decimal64(-my_value);
	}

	/**
	 * {@inherit_doc}
	 *
	 * The current implementation strictly enforces
	 * {@code this.multiply(a).equals(new Decimal64(this.double_value()
	 * * a.double_value()))}.
	 */
	//override
	Decimal64 multiply(const Decimal64& a)
	{
		return Decimal64(my_value * a.value);
	}

	/**
	 * {@inherit_doc}
	 *
	 * The current implementation strictly enforces
	 * {@code this.multiply(n).equals(new Decimal64(n * this.double_value()))}.
	 */
	//override
	Decimal64 multiply(const int& n)
	{
		return Decimal64(n * my_value);
	}

	/**
	 * {@inherit_doc}
	 *
	 * The current implementation strictly enforces
	 * {@code this.divide(a).equals(new Decimal64(this.double_value()
	 * / a.double_value()))}.
	 *
	 */
	//override
	Decimal64 divide(const Decimal64& a)
	{
		return Decimal64(my_value / a.value);
	}

	/**
	 * {@inherit_doc}
	 *
	 * The current implementation strictly enforces
	 * {@code this.reciprocal().equals(new Decimal64(1.0
	 * / this.double_value()))}.
	 */
	//override
	Decimal64 reciprocal()
	{
		return Decimal64(1.0 / my_value);
	}

	/*
	 * Methods from the Number virtual class
	 */

	/**
	 * {@inherit_doc}
	 *
	 * The current implementation performs casting to a {@code byte}.
	 */
	//override
	std::byte byte_value() const
	{
		return (byte)value;
	}

	/**
	 * {@inherit_doc}
	 *
	 * The current implementation performs casting to a {@code short}.
	 */
	//override
	short short_value() const
	{
		return static_cast<short>(my_value);
	}

	/**
	 * {@inherit_doc}
	 *
	 * The current implementation performs casting to a {@code int}.
	 */
	//override
	int int_value() const
	{
		return static_cast<int>(my_value);
	}

	/**
	 * {@inherit_doc}
	 *
	 * The current implementation performs casting to a {@code long}.
	 */
	//override
	long long_value() const
	{
		return static_cast<long>(my_value);
	}

	/**
	 * {@inherit_doc}
	 *
	 * The current implementation performs casting to a {@code float}.
	 */
	//override
	float float_value() const
	{
		return static_cast<float>(my_value);
	}

	/** {@inherit_doc} */
	//override
	double double_value() const
	{
		return my_value;
	}

	/*
	 * Methods from the Comparable interface.
	 */

	/**
	 * {@inherit_doc}
	 *
	 * The current implementation returns the same value as
	 * <center> {@code Double(this.double_value()).compare_to(new
	 * Double(o.double_value()))} </center>
	 *
	 * @see Double#compare_tostatic_cast<double>(
	 */
	//override
	int compare_to(const Decimal64& o)
	{
		return Double.compare(my_value, o.value);
	}

	/*
	 * Methods from the Object virtual class.
	 */

	/** {@inherit_doc} */
	//override
	bool equals(const Object& obj)
	{
		if (dynamic_cast<const Decimal64*>(*obj) != nullptr)
		{
			const Decimal64 that = (Decimal64)obj;
			return Double.double_to_long_bits(my_value) == Double.double_to_long_bits(that.value);
		}
		return false;
	}

	/** {@inherit_doc}
	 * <p>
	 * This implementation considers +0.0 and -0.0 to be equal.
	 * </p>
	 * @since 1.8
	 */
	//override
	bool is_zero() const
	{
		return my_value == 0.0;
	}

	/**
	 * {@inherit_doc}
	 *
	 * The current implementation returns the same value as
	 * {@code Double(this.double_value()).hash_code()}
	 *
	 * @see Double#hash_code()
	 */
	//override
	int hash_code()
	{
		long v = Double.double_to_long_bits(my_value);
		return static_cast<int>((v ^ (v >> > 32));
	}

	/**
	 * {@inherit_doc}
	 *
	 * The returned {@code std::string} is equal to
	 * {@code std::to_string(this.double_value())}
	 *
	 * @see Double#to_stringstatic_cast<double>(
	 */
	//override
	std::string to_string() const
	{
		return std::to_string(my_value);
	}

	/*
	 * Methods inspired by the Double class.
	 */

	/**
	 * Returns {@code true} if {@code this} double precision number is infinite
	 * ({@link Double#POSITIVE_INFINITY} or {@link Double#NEGATIVE_INFINITY}).
	 *
	 * @return {@code true} if {@code this} number is infinite
	 */
	//override
	bool is_infinite() const
	{
		return std::isinf(my_value);
	}

	/**
	 * Returns {@code true} if {@code this} double precision number is
	 * Not-a-Number ({@code NaN}), false otherwise.
	 *
	 * @return {@code true} if {@code this} is {@code NaN}
	 */
	//override
	bool is_nan() const
	{
		return std::isnan(my_value);
	}

	/** {@inherit_doc} */
	//override
	double get_real() const
	{
		return my_value;
	}

	/** {@inherit_doc} */
	//override
	Decimal64 add(const double& a)
	{
		return Decimal64(my_value + a);
	}

	/** {@inherit_doc} */
	//override
	Decimal64 subtract(const double& a)
	{
		return Decimal64(my_value - a);
	}

	/** {@inherit_doc} */
	//override
	Decimal64 multiply(const double& a)
	{
		return Decimal64(my_value * a);
	}

	/** {@inherit_doc} */
	//override
	Decimal64 divide(const double& a)
	{
		return Decimal64(my_value / a);
	}

	/** {@inherit_doc} */
	//override
	Decimal64 remainder(const double& a) const
	{
		return Decimal64(std::remainder(my_value, a));
	}

	/** {@inherit_doc} */
	//override
	Decimal64 remainder(const Decimal64& a) const
	{
		return Decimal64(std::remainder(my_value, a.value));
	}

	/** {@inherit_doc} */
	//override
	Decimal64 abs() const
	{
		return Decimal64(std::abs(my_value));
	}

	/** {@inherit_doc} */
	//override
	Decimal64 ceil() const
	{
		return Decimal64(std::ceil(my_value));
	}

	/** {@inherit_doc} */
	//override
	Decimal64 floor() const
	{
		return Decimal64(std::floor(my_value));
	}

	/** {@inherit_doc} */
	//override
	Decimal64 rint() const
	{
		return Decimal64(std::rint(my_value));
	}

	/** {@inherit_doc} */
	//override
	Decimal64 sign() const
	{
		return Decimal64(FastMath.signum(my_value));
	}

	/** {@inherit_doc} */
	//override
	Decimal64 copy_sign(const Decimal64& sign) const
	{
		return Decimal64(std::copysign(my_value, sign.value));
	}

	/** {@inherit_doc} */
	//override
	Decimal64 copy_sign(const double& sign) const
	{
		return Decimal64(std::copysign(my_value, sign));
	}

	/** {@inherit_doc} */
	//override
	Decimal64 scalb(const int& n) const
	{
		return Decimal64(std::scalbn(my_value, n));
	}

	/** {@inherit_doc} */
	//override
	Decimal64 ulp() const
	{
		return Decimal64(FastMath.ulp(my_value));
	}

	/** {@inherit_doc} */
	//override
	Decimal64 hypot(const Decimal64& y) const
	{
		return Decimal64(std::hypot(my_value, y.value));
	}

	/** {@inherit_doc} */
	//override
	Decimal64 sqrt() const
	{
		return Decimal64(std::sqrt(my_value));
	}

	/** {@inherit_doc} */
	//override
	Decimal64 cbrt() const
	{
		return Decimal64(std::cbrt(my_value));
	}

	/** {@inherit_doc} */
	//override
	Decimal64 root_n(const int& n) const
	{
		if (my_value < 0)
		{
			return (n % 2 == 0)
				? Decimal64.NAN
				: Decimal64(-std::pow(-value, 1.0 / n));
		}
		return Decimal64(std::pow(my_value, 1.0 / n));
	}

	/** {@inherit_doc} */
	//override
	Decimal64 pow(const double& p) const
	{
		return Decimal64(std::pow(my_value, p));
	}

	/** {@inherit_doc} */
	//override
	Decimal64 pow(const int& n) const
	{
		return Decimal64(std::pow(my_value, n));
	}

	/** {@inherit_doc} */
	//override
	Decimal64 pow(const Decimal64& e) const
	{
		return Decimal64(std::pow(my_value, e.value));
	}

	/** {@inherit_doc} */
	//override
	Decimal64 exp() const
	{
		return Decimal64(std::exp(my_value));
	}

	/** {@inherit_doc} */
	//override
	Decimal64 expm1() const
	{
		return Decimal64(std::expm1(my_value));
	}

	/** {@inherit_doc} */
	//override
	Decimal64 log() const
	{
		return Decimal64(std::log(my_value));
	}

	/** {@inherit_doc} */
	//override
	Decimal64 log1p() const
	{
		return Decimal64(std::log1p(my_value));
	}

	/** Base 10 logarithm.
	* @return base 10 logarithm of the instance
	*/
	//override
	Decimal64 log10() const
	{
		return Decimal64(std::log10(my_value));
	}

	/** {@inherit_doc} */
	//override
	Decimal64 cos() const
	{
		return Decimal64(std::cos(my_value));
	}

	/** {@inherit_doc} */
	//override
	Decimal64 sin() const
	{
		return Decimal64(std::sin(my_value));
	}

	/** {@inherit_doc} */
	//override
	Field_Sin_Cos<Decimal64> sin_cos() const
	{
		const auto sc = Sin_Cos(my_value);
		return Field_Sin_Cos<>(Decimal64(sc.sin()), Decimal64(sc.cos()));
	}

	/** {@inherit_doc} */
	//override
	Decimal64 tan() const
	{
		return Decimal64(std::tan(my_value));
	}

	/** {@inherit_doc} */
	//override
	Decimal64 acos() const
	{
		return Decimal64(std::acos(my_value));
	}

	/** {@inherit_doc} */
	//override
	Decimal64 asin() const
	{
		return Decimal64(std::asin(my_value));
	}

	/** {@inherit_doc} */
	//override
	Decimal64 atan() const
	{
		return Decimal64(std::atan(my_value));
	}

	/** {@inherit_doc} */
	//override
	Decimal64 atan2(const Decimal64& x) const
	{
		return Decimal64(std::atan2(my_value, x.value));
	}

	/** {@inherit_doc} */
	//override
	Decimal64 cosh() const
	{
		return Decimal64(std::cosh(my_value));
	}

	/** {@inherit_doc} */
	//override
	Decimal64 sinh() const
	{
		return Decimal64(std::sinh(my_value));
	}

	/** {@inherit_doc} */
	//override
	Field_Sinh_Cosh<Decimal64> sinh_cosh() const
	{
		const Sinh_Cosh sch = std::sinh_cosh(my_value);
		return Field_Sinh_Cosh<>(Decimal64(sch.sinh()), Decimal64(sch.cosh()));
	}

	/** {@inherit_doc} */
	//override
	Decimal64 tanh() const
	{
		return Decimal64(std::tanh(my_value));
	}

	/** {@inherit_doc} */
	//override
	Decimal64 acosh() const
	{
		return Decimal64(std::acosh(my_value));
	}

	/** {@inherit_doc} */
	//override
	Decimal64 asinh() const
	{
		return Decimal64(std::asinh(my_value));
	}

	/** {@inherit_doc} */
	//override
	Decimal64 atanh() const
	{
		return Decimal64(std::atanh(my_value));
	}

	/** {@inherit_doc} */
	//override
	Decimal64 to_degrees() const
	{
		return Decimal64(FastMath.to_degrees(my_value));
	}

	/** {@inherit_doc} */
	//override
	Decimal64 to_radians() const
	{
		return Decimal64(FastMath.to_radians(my_value));
	}

	/** {@inherit_doc} */
	//override
	Decimal64 linear_combination(const std::vector<Decimal64>& a, const std::vector<Decimal64>& b) const
	{
		Math_Utils::check_dimension(a.size(), b.size());
		auto a_double = std::vector<double>(a.size());
		auto b_double = std::vector<double>(b.size());
		for (int i{}; i < a.size(); ++i)
		{
			a_double[i] = a[i].value;
			b_double[i] = b[i].value;
		}
		return Decimal64(Math_Arrays::linear_combination(a_double, b_double));
	}

	/** {@inherit_doc} */
	//override
	Decimal64 linear_combination(const std::vector<double>& a, const std::vector<Decimal64>& b) const
	{
		Math_Utils::check_dimension(a.size(), b.size());
		auto b_double = std::vector<double>(b.size());
		for (int i{}; i < a.size(); ++i)
		{
			b_double[i] = b[i].value;
		}
		return Decimal64(Math_Arrays::linear_combination(a, b_double));
	}

	/** {@inherit_doc} */
	//override
	Decimal64 linear_combination(const Decimal64& a1, const Decimal64& b1, const Decimal64& a2, const Decimal64& b2) const
	{
		return Decimal64(Math_Arrays::linear_combination(a1.value, b1.value, a2.value, b2.value));
	}

	/** {@inherit_doc} */
	//override
	Decimal64 linear_combination(const double& a1, const Decimal64& b1, const double& a2, const Decimal64& b2) const
	{
		return Decimal64(Math_Arrays::linear_combination(a1, b1.value, a2, b2.value));
	}

	/** {@inherit_doc} */
	//override
	Decimal64 linear_combination(const Decimal64& a1, const Decimal64& b1, const Decimal64& a2, const Decimal64& b2, const Decimal64& a3, const Decimal64& b3) const
	{
		return Decimal64(Math_Arrays::linear_combination(a1.value, b1.value, a2.value, b2.value, a3.value, b3.value));
	}

	/** {@inherit_doc} */
	//override
	Decimal64 linear_combination(const double& a1, const Decimal64& b1, const double& a2, const Decimal64& b2, const double& a3, const Decimal64& b3) const
	{
		return Decimal64(Math_Arrays::linear_combination(a1, b1.value, a2, b2.value, a3, b3.value));
	}

	/** {@inherit_doc} */
	//override
	Decimal64 linear_combination(const Decimal64& a1, const Decimal64& b1, const Decimal64& a2, const Decimal64& b2, const Decimal64& a3, const Decimal64& b3, const Decimal64& a4, const Decimal64& b4) const
	{
		return Decimal64(Math_Arrays::linear_combination(a1.value, b1.value, a2.value, b2.value, a3.value, b3.value, a4.value, b4.value));
	}

	/** {@inherit_doc} */
	//override
	Decimal64 linear_combination(const double& a1, const Decimal64& b1, const double& a2, const Decimal64& b2, const double& a3, const Decimal64& b3, const double& a4, const Decimal64& b4) const
	{
		return Decimal64(Math_Arrays::linear_combination(a1, b1.value, a2, b2.value, a3, b3.value, a4, b4.value));
	}

	/** {@inherit_doc} */
	//override
	Decimal64 get_pi() const
	{
		return PI;
	}
};