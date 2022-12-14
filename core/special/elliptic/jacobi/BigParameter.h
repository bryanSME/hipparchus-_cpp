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

#include <cmath>
#include "JacobiElliptic.h"
#include "JacobiEllipticBuilder.hpp"
#include "CopolarN.h"

 //import org.hipparchus.util.FastMath;

 /** Algorithm for computing the principal Jacobi functions for parameter m greater than 1.
  * <p>
  * The rules for reciprocal parameter change are given in Abramowitz and Stegun, * sections 16.11 and 17.4.15.
  * </p>
  * @since 2.0
  */
class Big_Parameter : public Jacobi_Elliptic
{
private:
	/** Algorithm to use for the positive parameter. */
	const Jacobi_Elliptic my_algorithm;

	/** Input scaling factor. */
	const double my_input_scale;

	/** output scaling factor. */
	const double my_output_scale;

public:
	/** Simple constructor.
	 * @param m parameter of the Jacobi elliptic function (must be greater than 1 here)
	 */
	Big_Parameter(const double& m)
		:
		my_algorithm{ Jacobi_Elliptic_Builder::build(1.0 / m) },
		my_input_scale{ std::sqrt(m) },
		my_output_scale{ 1.0 / my_input_scale }
	{
		//super(m);
	};

	/** {@inherit_doc} */
	//override
	Copolar_N values_n(const double& u)
	{
		const Copolar_N trio_n = my_algorithm.values_n(u * my_input_scale);
		return Copolar_N(my_output_scale * trio_n.sn(), trio_n.dn(), trio_n.cn());
	}
};