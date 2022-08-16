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

#include <type_traits>
#include <complex>
#include "../../../CalculusFieldElement.hpp"
#include "NegativeParameter.h"
#include "BigParameter.h"
#include "NearZeroParameter.h"
#include "NearOneParameter.h"
#include "JacobiElliptic.h"
#include "FieldJacobiElliptic.h"

 /** Builder for algorithms compmuting Jacobi elliptic functions.
  * <p>
  * The Jacobi elliptic functions are related to elliptic integrals.
  * </p>
  * <p>
  * There are different conventions to interpret the arguments of
  * Jacobi elliptic functions. The first argument may be  the amplitude \xcf\x86, * but is more often the variable u (with sn(u) = sin(\xcf\x86) and cn(u) = cos(\xcf\x86)).
  * The second argument  is either the modulus k or the parameter m with m = k\xc2\xb2.
  * In Hipparchus, we adopted the convention to use u and m.
  * </p>
  * @since 2.0
  */
class Jacobi_Elliptic_Builder
{
private:
	/** Threshold near 0 for using specialized algorithm. */
	static constexpr double NEAR_ZERO{ 1.0e-9 };

	/** Threshold near 1 for using specialized algorithm. */
	static constexpr double NEAR_ONE{ 1.0 - NEAR_ZERO };

	/**
	 * Private constructor for utility class.
	 * nothing to do
	 */
	Jacobi_Elliptic_Builder() = default;

public:
	/** Build an algorithm for computing Jacobi elliptic functions.
	 * @param m parameter of the Jacobi elliptic function
	 * @return selected algorithm
	 */
	static Jacobi_Elliptic build(const double& m)
	{
		if (m < 0)
		{
			return Negative_Parameter(m);
		}
		if (m > 1)
		{
			return Big_Parameter(m);
		}
		if (m < NEAR_ZERO)
		{
			return Near_Zero_Parameter(m);
		}
		if (m > NEAR_ONE)
		{
			return Near_One_Parameter(m);
		}
		return Bounded_Parameter(m);
	}

	/** Build an algorithm for computing Jacobi elliptic functions.
	 * @param m parameter of the Jacobi elliptic function
	 * @param <T> type of the field elements
	 * @return selected algorithm
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	static Field_Jacobi_Elliptic<T> build(const T& m)
	{
		if (m.get_real() < 0)
		{
			return Field_Negative_Parameter<>(m);
		}
		if (m.get_real() > 1)
		{
			return FieldBig_Parameter<>(m);
		}
		if (m.get_real() < NEAR_ZERO)
		{
			return Field_Near_Zero_Parameter<>(m);
		}
		if (m.get_real() > NEAR_ONE)
		{
			return Field_Near_One_Parameter<>(m);
		}
		return FieldBounded_Parameter<>(m);
	}

	/** Build an algorithm for computing Jacobi elliptic functions.
	 * @param m parameter of the Jacobi elliptic function
	 * @return selected algorithm
	 */
	static Field_Jacobi_Elliptic<std::complex<double>> build(const std::complex<double>& m)
	{
		return Complex_Parameter(m);
	}

	/** Build an algorithm for computing Jacobi elliptic functions.
	 * @param m parameter of the Jacobi elliptic function
	 * @param <T> type of the field elements
	 * @return selected algorithm
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	static Field_Jacobi_Elliptic<Field_Complex<T>> build(const Field_Complex<T>& m)
	{
		return Field_Complex_Parameter<>(m);
	}
};