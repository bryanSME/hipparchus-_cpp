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

 //import org.hipparchus.complex.std::complex<double>;
 //import org.hipparchus.special.elliptic.legendre.Legendre_Elliptic_Integral;
 //import org.hipparchus.util.FastMath;
 //import org.hipparchus.util.Math_Utils;

 /** Algorithm for computing the principal Jacobi functions for complex parameter m.
  * @since 2.0
  */
class Complex_Parameter : public Field_Jacobi_Elliptic<std::complex<double>>
{
private:
	/** Jacobi θ functions. */
	const Field_Jacobi_Theta<std::complex<double>> my_jacobi_theta;

	/** Quarter period K. */
	const std::complex<double> my_big_k;

	/** Quarter period iK'. */
	const std::complex<double> my_i_big_k_prime;

	/** Real periodic factor for K. */
	const double my_rK;

	/** Imaginary periodic factor for K. */
	const double my_iK;

	/** Real periodic factor for iK'. */
	const double my_r_k_prime;

	/** Imaginary periodic factor for iK'. */
	const double my_i_k_prime;

	/** Value of Jacobi θ functions at origin. */
	const Field_Theta<std::complex<double>> my_t0;

	/** Scaling factor. */
	const std::complex<double> my_scaling;

public:
	/** Simple constructor.
	 * @param m parameter of the Jacobi elliptic function
	 */
	Complex_Parameter(const std::complex<double>& m)
	{
		super(m);

		// compute nome
		const std::complex<double> q = Legendre_Elliptic_Integral.nome(m);

		// compute periodic factors such that
		// z = 4 K [rK Re(z) + iK Im(z)] + 4i K' [rK' Re(z) + iK' Im(z)]
		my_big_k = Legendre_Elliptic_Integral.big_k(m);
		my_i_big_k_prime = Legendre_Elliptic_Integral.big_k_prime(m).multiply_plus_i();
		const double inverse = 0.25 /
			(big_k.get_real_part() * i_big_k_prime.get_imaginary_part() -
				big_k.get_imaginary_part() * i_big_k_prime.get_real_part());
		my_rK = i_big_k_prime.get_imaginary_part() * inverse;
		my_iK = i_big_k_prime.get_real_part() * -inverse;
		my_r_k_prime = big_k.get_imaginary_part() * -inverse;
		my_i_k_prime = big_k.get_real_part() * inverse;

		// prepare underlying Jacobi θ functions
		my_jacobi_theta = Field_Jacobi_Theta<>(q);
		my_t0 = jacobi_theta.values(m.get_field().get_zero());
		my_scaling = big_k.reciprocal().multiply(Math_Utils::SEMI_PI);
	}

	/** {@inherit_doc}
	 * <p>
	 * The algorithm for evaluating the functions is based on {@link Field_Jacobi_Theta
	 * Jacobi theta functions}.
	 * </p>
	 */
	 //override
	Field_Copolar_N<std::complex<double>> values_n(const std::complex<double>& u)
	{
		// perform argument reduction
		const double cK = rK * u.get_real_part() + iK * u.get_imaginary_part();
		const double c_k_prime = r_k_prime * u.get_real_part() + i_k_prime * u.get_imaginary_part();
		const std::complex<double> reduced_u = u.linear_combination(1.0, u, -4 * std::rint(cK), big_k, -4 * std::rint(c_k_prime), i_big_k_prime);

		// evaluate Jacobi θ functions at argument
		const Field_Theta<std::complex<double>> t_z = jacobi_theta.values(reduced_u.multiply(scaling));

		// convert to Jacobi elliptic functions
		const std::complex<double> sn = t0.theta3().multiply(t_z.theta1()).divide(t0.theta2().multiply(t_z.theta4()));
		const std::complex<double> cn = t0.theta4().multiply(t_z.theta2()).divide(t0.theta2().multiply(t_z.theta4()));
		const std::complex<double> dn = t0.theta4().multiply(t_z.theta3()).divide(t0.theta3().multiply(t_z.theta4()));

		return Field_Copolar_N<>(sn, cn, dn);
	}
};