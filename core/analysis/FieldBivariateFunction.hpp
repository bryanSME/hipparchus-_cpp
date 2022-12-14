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
#include "../CalculusFieldElement.hpp"
#include "../Field.h"
#include "CalculusFieldBivariateFunction.hpp"

/**
 * An interface representing a bivariate field function.
 * @since 1.5
 */
class Field_Bivariate_Function
{
	/** Convert to a {@link Calculus_Field_Bivariate_Function} with a specific type.
	 * @param <T> the type of the field elements
	 * @param field field for the argument and value
	 * @return converted function
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	Calculus_Field_Bivariate_Function<T> to_calculus_field_bivariate_function(const Field<T>& field)
	{
		return this::value;
	}

	/**
	 * Compute the value for the function.
	 *
	 * @param x Abscissa for which the function value should be computed.
	 * @param y Ordinate for which the function value should be computed.
	 * @param <T> type of the field elements
	 * @return the value.
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	virtual T value(const T& x, const T& y) = 0;
};