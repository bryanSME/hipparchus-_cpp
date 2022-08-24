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

 //import org.hipparchus.Calculus_Field_Element;
 //import org.hipparchus.util.Pair;
#include <type_traits>
#include "../../../FieldElement.h"


 /**
  * Factory converting {@link Calculus_Field_Element field-based} {@link FieldRule_Factory} into {@link Rule_Factory}.
  * @param <T> Type of the number used to represent the points and weights of
  * the quadrature rules.
  * @since 2.0
  */
template<typename T, typename std::enable_if<std::is_base_of<Field_Element<T>, T>::value>::type* = nullptr>
class ConvertingRule_Factory : public AbstractRule_Factory
{
private:
	/** Underlying field-based factory. */
	const FieldRule_Factory<T> field_factory;

public:
	/** Simple constructor.
	 * @param field_factory field-based factory to convert
	 */
	ConvertingRule_Factory(const FieldRule_Factory<T>& field_factory)
		:
		my_field_factory{ field_factory }
	{

	}

protected:
	/** {@inherit_doc} */
	//override
	Pair<std::vector<double>, std::vector<double>> compute_rule(const int& number_of_points)
	{
		// get the field-based rule
		Pair<std::vector<T>, std::vector<T>> rule = field_factory.get_rule(number_of_points);

		// convert the nodes and weights
		const std::vector<T> pT = rule.get_first();
		const std::vector<T> wT = rule.get_second();

		const int len = pT.size();
		const std::vector<double> pD = std::vector<double>(len];
		const std::vector<double> wD = std::vector<double>(len];

		for (int i{}; i < len; i++)
		{
			pD[i] = pT[i].get_real();
			wD[i] = wT[i].get_real();
		}

		return Pair<>(pD, wD);
	}
};