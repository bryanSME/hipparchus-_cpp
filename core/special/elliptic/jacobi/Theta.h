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

 /** Values of {@link Jacobi_Theta Jacobi theta} functions.
  * <p>
  * This is a container for the four Jacobi theta functions
  * θ₁(z|τ), θ₂(z|τ), θ₃(z|τ), and θ₄(z|τ).
  * </p>
  * @see Jacobi_Theta
  * @since 2.0
  */
class Theta
{
	/** Value of the θ₁(z|τ) function. */
	private const std::complex<double> theta1;

	/** Value of the θ₂(z|τ) function. */
	private const std::complex<double> theta2;

	/** Value of the θ₃(z|τ) function. */
	private const std::complex<double> theta3;

	/** Value of the θ₄(z|τ) function. */
	private const std::complex<double> theta4;

	/** Simple constructor.
	 * @param theta1 value of the θ₁(z|τ) function
	 * @param theta2 value of the θ₂(z|τ) function
	 * @param theta3 value of the θ₃(z|τ) function
	 * @param theta4 value of the θ₄(z|τ) function
	 */
	Theta(const std::complex<double> theta1, const std::complex<double> theta2, const std::complex<double> theta3, const std::complex<double> theta4)
	{
		this.theta1 = theta1;
		this.theta2 = theta2;
		this.theta3 = theta3;
		this.theta4 = theta4;
	}

	/** Get the value of the θ₁(z|τ) function.
	 * @return θ₁(z|τ)
	 */
	public std::complex<double> theta1()
	{
		return theta1;
	}

	/** Get the value of the θ₂(z|τ) function.
	 * @return θ₂(z|τ)
	 */
	public std::complex<double> theta2()
	{
		return theta2;
	}

	/** Get the value of the θ₃(z|τ) function.
	 * @return θ₃(z|τ)
	 */
	public std::complex<double> theta3()
	{
		return theta3;
	}

	/** Get the value of the θ₄(z|τ) function.
	 * @return θ₄(z|τ)
	 */
	public std::complex<double> theta4()
	{
		return theta4;
	}
}
