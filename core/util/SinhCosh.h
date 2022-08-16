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
 //package org.hipparchus.util;

 /** Holder for both hyperbolic sine and hyperbolic cosine values.
  * <p>
  * This class is a simple container, it does not provide any computational method.
  * </p>
  * @see FastMath#sinh_coshstatic_cast<double>(
  * @since 2.0
  */
class Sinh_Cosh
{
	/** Value of the hyperbolic sine. */
	private const double sinh;

	/** Value of the hyperbolic cosine. */
	private const double cosh;

	/** Simple constructor.
	 * @param sinh value of the hyperbolic sine
	 * @param cosh value of the hyperbolic cosine
	 */
	Sinh_Cosh(const double sinh, const double cosh)
	{
		this.sinh = sinh;
		this.cosh = cosh;
	}

	/** Get the value of the hyperbolic sine.
	 * @return value of the hyperbolic sine
	 */
	public double sinh()
	{
		return sinh;
	}

	/** Get the value of the hyperbolic cosine.
	 * @return value of the hyperbolic cosine
	 */
	public double cosh()
	{
		return cosh;
	}

	/** Compute hyperbolic sine and hyperbolic cosine of angles sum.
	 * @param sch_alpha \((\sinh \alpha, \cosh \alpha)\)
	 * @param sch_beta \((\sinh \beta, \cosh \beta)\)
	 * @return \((\sinh \alpha+\beta, \cosh \alpha+\beta)\)
	 */
	public static Sinh_Cosh sum(const Sinh_Cosh sch_alpha, const Sinh_Cosh sch_beta)
	{
		return Sinh_Cosh(Math_Arrays::linear_combination(sch_alpha.sinh, sch_beta.cosh, sch_alpha.cosh, sch_beta.sinh), Math_Arrays::linear_combination(sch_alpha.cosh, sch_beta.cosh, sch_alpha.sinh, sch_beta.sinh));
	}

	/** Compute hyperbolic sine and hyperbolic cosine of angles difference.
	 * @param sch_alpha \((\sinh \alpha, \cosh \alpha)\)
	 * @param sch_beta \((\sinh \beta, \cosh \beta)\)
	 * @return \((\sinh \alpha+\beta, \cosh \alpha-\beta)\)
	 */
	public static Sinh_Cosh difference(const Sinh_Cosh sch_alpha, const Sinh_Cosh sch_beta)
	{
		return Sinh_Cosh(Math_Arrays::linear_combination(sch_alpha.sinh, sch_beta.cosh, sch_alpha.cosh, -sch_beta.sinh), Math_Arrays::linear_combination(sch_alpha.cosh, sch_beta.cosh, sch_alpha.sinh, -sch_beta.sinh));
	}
}
