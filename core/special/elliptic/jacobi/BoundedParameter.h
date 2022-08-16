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

#include "JacobiElliptic.h"
#include "JacobiTheta.h"
#include "Theta.h"
#include "../../../util/MathUtils.h"

 /** Algorithm for computing the principal Jacobi functions for parameter m in [0; 1].
  * @since 2.0
  */
class Bounded_Parameter : public Jacobi_Elliptic
{
private:
	/** Jacobi θ functions. */
	const Jacobi_Theta my_jacobi_theta;

	/** Value of Jacobi θ functions at origin. */
	const Theta my_t0;

	/** Scaling factor. */
	const double my_scaling;

public:
	/** Simple constructor.
	 * @param m parameter of the Jacobi elliptic function
	 */
	Bounded_Parameter(const double& m)
	{
		super(m);

		// compute nome
		const double q = Legendre_Elliptic_Integral.nome(m);

		// prepare underlying Jacobi θ functions
		my_jacobi_theta = Jacobi_Theta(q);
		my_t0 = my_jacobi_theta.values(std::complex<double>.ZERO);
		my_scaling = Math_Utils::SEMI_PI / Legendre_Elliptic_Integral.big_k(m);
	}

	/** {@inherit_doc}
	 * <p>
	 * The algorithm for evaluating the functions is based on {@link Jacobi_Theta
	 * Jacobi theta functions}.
	 * </p>
	 */
	 //override
	Copolar_N values_n(const double& u)
	{
		// evaluate Jacobi θ functions at argument
		const Theta t_z = jacobi_theta.values(new std::complex<double>(u * scaling));

		// convert to Jacobi elliptic functions
		const double sn = my_t0.theta3().multiply(t_z.theta1()).divide(my_t0.theta2().multiply(t_z.theta4())).get_real_part();
		const double cn = my_t0.theta4().multiply(t_z.theta2()).divide(my_t0.theta2().multiply(t_z.theta4())).get_real_part();
		const double dn = my_t0.theta4().multiply(t_z.theta3()).divide(my_t0.theta3().multiply(t_z.theta4())).get_real_part();

		return Copolar_N(sn, cn, dn);
	}
};