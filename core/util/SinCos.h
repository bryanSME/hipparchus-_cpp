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

#include <exception>

 /** Holder for both sine and cosine values.
  * <p>
  * This class is a simple container, it does not provide any computational method.
  * </p>
  * @see FastMath#sin_cosstatic_cast<double>(
  * @since 1.3
  */
class Sin_Cos
{
private:
	/** Value of the sine. */
	const double my_sin;

	/** Value of the cosine. */
	const double my_cos;

public:
	/** Simple constructor.
	 * @param sin value of the sine
	 * @param cos value of the cosine
	 */
	Sin_Cos(const double sin, const double cos) : my_sin{ sin }, my_cos{ cos } {};

	/** Get the value of the sine.
	 * @return value of the sine
	 */
	double sin() const
	{
		return my_sin;
	}

	/** Get the value of the cosine.
	 * @return value of the cosine
	 */
	double cos() const
	{
		return my_cos;
	}

	/** Compute sine and cosine of angles sum.
	 * @param sc_alpha \((\sin \alpha, \cos \alpha)\)
	 * @param sc_beta \((\sin \beta, \cos \beta)\)
	 * @return \((\sin \alpha+\beta, \cos \alpha+\beta)\)
	 * @since 1.8
	 */
	static Sin_Cos sum(const Sin_Cos& sc_alpha, const Sin_Cos& sc_beta)
	{
		throw std::exception("not implemented");
		//return Sin_Cos(Math_Arrays::linear_combination(sc_alpha.sin(), sc_beta.cos(), sc_alpha.cos(), sc_beta.sin()), Math_Arrays::linear_combination(sc_alpha.cos(), sc_beta.cos(), -sc_alpha.sin(), sc_beta.sin()));
	}

	/** Compute sine and cosine of angles difference.
	 * @param sc_alpha \((\sin \alpha, \cos \alpha)\)
	 * @param sc_beta \((\sin \beta, \cos \beta)\)
	 * @return \((\sin \alpha+\beta, \cos \alpha-\beta)\)
	 * @since 1.8
	 */
	static Sin_Cos difference(const Sin_Cos& sc_alpha, const Sin_Cos& sc_beta)
	{
		throw std::exception("not implemented");
		//return Sin_Cos(Math_Arrays::linear_combination(sc_alpha.sin(), sc_beta.cos(), -sc_alpha.cos(), sc_beta.sin()), Math_Arrays::linear_combination(sc_alpha.cos(), sc_beta.cos(), sc_alpha.sin(), sc_beta.sin()));
	}
};