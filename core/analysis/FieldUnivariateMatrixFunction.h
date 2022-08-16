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

 //import org.hipparchus.Field;
 //import org.hipparchus.Calculus_Field_Element;
#include <type_traits>
#include "../CalculusFieldElement.hpp"

/**
 * An interface representing a univariate matrix function for any field type.
 *
 * @since 1.3
 */
class FieldUnivariate_Matrix_Function
{
	/** Convert to a {@link CalculusFieldUnivariate_Matrix_Function} with a specific type.
	 * @param <T> the type of the field elements
	 * @param field field for the argument and value
	 * @return converted function
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	CalculusFieldUnivariate_Matrix_Function<T> to_calculus_field_univariate__matrix__function(Field<T> field)
	{
		return this::value;
	}

	/**
	 * Compute the value for the function.
	 * @param <T> the type of the field elements
	 * @param x the point for which the function value should be computed
	 * @return the value
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	virtual  std::vector<std::vector<T>> value(T x) = 0;
};