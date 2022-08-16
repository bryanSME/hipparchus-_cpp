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
 //import org.hipparchus.analysis.Calculus_Field_Univariate_Function;
#include <type_traits>
#include "../../CalculusFieldElement.hpp"
#include "../CalculusFieldUnivariateFunction.hpp"

/**
 * Interface representing a univariate field interpolating function.
 * @since 1.5
 */
class Field_Univariate_Interpolator
{
	/**
	 * Compute an interpolating function for the dataset.
	 *
	 * @param xval Arguments for the interpolation points.
	 * @param yval Values for the interpolation points.
	 * @param <T> the type of the field elements
	 * @return a function which interpolates the dataset.
	 * @
	 * if the arguments violate assumptions made by the interpolation
	 * algorithm.
	 * @ if arrays lengthes do not match
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	Calculus_Field_Univariate_Function<T> interpolate(T xval[], T yval[]);
};