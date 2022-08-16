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
 //import org.hipparchus.util.FastMath;
#include <type_traits>
#include "../../../CalculusFieldElement.hpp"

/** Algorithm for computing the principal Jacobi functions for negative parameter m.
 * <p>
 * The rules for negative parameter change are given in Abramowitz and Stegun, section 16.10.
 * </p>
 * @param <T> the type of the field elements
 * @since 2.0
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
class Field_Negative_Parameter : public Field_Jacobi_Elliptic<T>
{
	/** Algorithm to use for the positive parameter. */
	private const Field_Jacobi_Elliptic<T> algorithm;

	/** Input scaling factor. */
	private const T input_scale;

	/** output scaling factor. */
	private const T output_scale;

	/** Simple constructor.
	 * @param m parameter of the Jacobi elliptic function (must be negative here)
	 */
	Field_Negative_Parameter(const T m)
	{
		super(m);
		const T om_m = m.get_field().get_one().subtract(m);
		algorithm = Jacobi_Elliptic_Builder.build(m.negate().divide(om_m));
		input_scale = std::sqrt(om_m);
		output_scale = input_scale.reciprocal();
	}

	/** {@inherit_doc} */
	//override
	public Field_Copolar_N<T> values_n(const T u)
	{
		const FieldCopolar_D<T> trio_d = FieldCopolar_D<>(algorithm.values_n(u.multiply(input_scale)));
		return Field_Copolar_N<>(output_scale.multiply(trio_d.sd()), trio_d.cd(), trio_d.nd());
	}
};