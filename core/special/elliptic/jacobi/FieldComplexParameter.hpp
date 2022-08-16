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
 //package org.hipparchus.special.elliptic.jacobi;

 //import org.hipparchus.Calculus_Field_Element;
 //import org.hipparchus.complex.Field_Complex<double>;
 //import org.hipparchus.special.elliptic.legendre.Legendre_Elliptic_Integral;
 //import org.hipparchus.util.FastMath;
#include <type_traits>
#include "../../../CalculusFieldElement.hpp"

/** Algorithm for computing the principal Jacobi functions for complex parameter m.
 * @param <T> the type of the field elements
 * @since 2.0
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
class Field_Complex_Parameter : public Field_Jacobi_Elliptic<Field_Complex<T>>
{
private:

	/** Jacobi \xce\xb8 functions. */
	const Field_Jacobi_Theta<Field_Complex<T>> my_jacobi_theta;

	/** Quarter period K. */
	const Field_Complex<T> my_big_k;

	/** Quarter period iK'. */
	const Field_Complex<T> my_i_big_k_prime;

	/** Real periodic factor for K. */
	const T my_rK;

	/** Imaginary periodic factor for K. */
	const T my_iK;

	/** Real periodic factor for iK'. */
	const T my_r_k_prime;

	/** Imaginary periodic factor for iK'. */
	const T my_i_k_prime;

	/** Value of Jacobi \xce\xb8 functions at origin. */
	const Field_Theta<Field_Complex<T>> my_t0;

	/** Scaling factor. */
	const Field_Complex<T> my_scaling;

public:
	/** Simple constructor.
	 * @param m parameter of the Jacobi elliptic function
	 */
	Field_Complex_Parameter(const Field_Complex<T> m)
	{
		super(m);

		// compute nome
		const Field_Complex<T> q = Legendre_Elliptic_Integral.nome(m);

		// compute periodic factors such that
		// z = 4K [rK Re(z) + iK Im(z)] + 4K' [rK' Re(z) + iK' Im(z)]
		my_big_k = Legendre_Elliptic_Integral.big_k(m);
		my_i_big_k_prime = Legendre_Elliptic_Integral.big_k_prime(m).multiply_plus_i();
		const T inverse = big_k.get_real_part().multiply(i_big_k_prime.get_imaginary_part()).
			subtract(big_k.get_imaginary_part().multiply(i_big_k_prime.get_real_part())).
			multiply(4).reciprocal();
		my_rK = i_big_k_prime.get_imaginary_part().multiply(inverse);
		my_iK = i_big_k_prime.get_real_part().multiply(inverse).negate();
		my_r_k_prime = big_k.get_imaginary_part().multiply(inverse).negate();
		my_i_k_prime = big_k.get_real_part().multiply(inverse);

		// prepare underlying Jacobi \xce\xb8 functions
		my_jacobi_theta = Field_Jacobi_Theta<>(q);
		my_t0 = jacobi_theta.values(m.get_field().get_zero());
		my_scaling = big_k.reciprocal().multiply(m.get_pi().multiply(0.5));
	}

	/** {@inherit_doc}
	 * <p>
	 * The algorithm for evaluating the functions is based on {@link Field_Jacobi_Theta
	 * Jacobi theta functions}.
	 * </p>
	 */
	 //override
	Field_Copolar_N<Field_Complex<T>> values_n(Field_Complex<T> u)
	{
		// perform argument reduction
		const T cK = rK.multiply(u.get_real_part()).add(iK.multiply(u.get_imaginary_part()));
		const T c_k_prime = r_k_prime.multiply(u.get_real_part()).add(i_k_prime.multiply(u.get_imaginary_part()));
		const Field_Complex<T> reduced_u = u.linear_combination(1.0, u, -4 * std::rint(cK.get_real()), big_k, -4 * std::rint(c_k_prime.get_real()), i_big_k_prime);

		// evaluate Jacobi \xce\xb8 functions at argument
		const Field_Theta<Field_Complex<T>> t_z = jacobi_theta.values(reduced_u.multiply(scaling));

		// convert to Jacobi elliptic functions
		const Field_Complex<T> sn = t0.theta3().multiply(t_z.theta1()).divide(t0.theta2().multiply(t_z.theta4()));
		const Field_Complex<T> cn = t0.theta4().multiply(t_z.theta2()).divide(t0.theta2().multiply(t_z.theta4()));
		const Field_Complex<T> dn = t0.theta4().multiply(t_z.theta3()).divide(t0.theta3().multiply(t_z.theta4()));

		return Field_Copolar_N<>(sn, cn, dn);
	}
};