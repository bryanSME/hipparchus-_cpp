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

  //package org.hipparchus.complex;

  //import org.hipparchus.Calculus_Field_Element;
  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Field_Sin_Cos;
  //import org.hipparchus.util.Sin_Cos;
#include <vector>
#include <type_traits>
#include <complex>
#include "../CalculusFieldElement.hpp"
#include "../util/SinCos.h"
#include "../complex/FieldComplex.h"

/**
 * Static implementations of common {@link std::complex<double>} utilities functions.
 */
template<typename T>
class Complex_Utils
{
private:
	/**
	 * Default constructor.
	 */
	Complex_Utils() {}

public:

	/**
	 * Creates a complex number from the given polar representation.
	 * <p>
	 * The value returned is <code>r&middot;e<sup>i&middot;theta</sup></code>, * computed as <code>r&middot;cos(theta) + r&middot;sin(theta)i</code></p>
	 * <p>
	 * If either <code>r</code> or <code>theta</code> is NaN, or
	 * <code>theta</code> is infinite, {@link std::complex<double>#NaN} is returned.</p>
	 * <p>
	 * If <code>r</code> is infinite and <code>theta</code> is finite, * infinite or NaN values may be returned in parts of the result, following
	 * the rules for double arithmetic.<pre>
	 * Examples:
	 * <code>
	 * polar_2_complex(INFINITY, &pi;/4) = INFINITY + INFINITY i
	 * polar_2_complex(INFINITY, 0) = INFINITY + NaN i
	 * polar_2_complex(INFINITY, -&pi;/4) = INFINITY - INFINITY i
	 * polar_2_complex(INFINITY, 5&pi;/4) = -INFINITY - INFINITY i </code></pre></p>
	 *
	 * @param r the modulus of the complex number to create
	 * @param theta  the argument of the complex number to create
	 * @return <code>r&middot;e<sup>i&middot;theta</sup></code>
	 * @ if {@code r} is negative.
	 */
	static std::complex<double> polar_2_complex(const double& r, const double& theta)
	{
		if (r < 0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NEGATIVE_COMPLEX_MODULE, r);
		}
		const auto sc = Sin_Cos(theta);
		return std::complex<double>(r * sc.cos(), r * sc.sin());
	}

	/**
	 * Creates a complex number from the given polar representation.
	 * <p>
	 * The value returned is <code>r&middot;e<sup>i&middot;theta</sup></code>, * computed as <code>r&middot;cos(theta) + r&middot;sin(theta)i</code></p>
	 * <p>
	 * If either <code>r</code> or <code>theta</code> is NaN, or
	 * <code>theta</code> is infinite, {@link std::complex<double>#NaN} is returned.</p>
	 * <p>
	 * If <code>r</code> is infinite and <code>theta</code> is finite, * infinite or NaN values may be returned in parts of the result, following
	 * the rules for double arithmetic.<pre>
	 * Examples:
	 * <code>
	 * polar_2_complex(INFINITY, &pi;/4) = INFINITY + INFINITY i
	 * polar_2_complex(INFINITY, 0) = INFINITY + NaN i
	 * polar_2_complex(INFINITY, -&pi;/4) = INFINITY - INFINITY i
	 * polar_2_complex(INFINITY, 5&pi;/4) = -INFINITY - INFINITY i </code></pre></p>
	 *
	 * @param r the modulus of the complex number to create
	 * @param theta  the argument of the complex number to create
	 * @param <T> type of the field elements
	 * @return <code>r&middot;e<sup>i&middot;theta</sup></code>
	 * @ if {@code r} is negative.
	 * @since 2.0
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	static Field_Complex<T> polar_2_complex(T r, T theta)
	{
		if (r.get_real() < 0)
		{
			throw std::exception("not implemented")
				//throw (hipparchus::exception::Localized_Core_Formats_Type::NEGATIVE_COMPLEX_MODULE, r);
		}
		const Field_Sin_Cos<T> sc = Sin_Cos(theta);
		return Field_Complex<double><>(r.multiply(sc.cos()), r.multiply(sc.sin()));
	}

	/**
	 * Convert an array of primitive doubles to an array of {@code std::complex<double>} objects.
	 *
	 * @param real Array of numbers to be converted to their {@code std::complex<double>} equivalent.
	 * @return an array of {@code std::complex<double>} objects.
	 */
	static std::vector<std::complex<double>>convert_to_complex(const std::vector<double>& real)
	{
		const c = std::vector<std::complex<double>>(real.size());
		for (int i{}; i < real.size(); i++)
		{
			c[i] = std::complex<double>(real[i], 0);
		}

		return c;
	}
};