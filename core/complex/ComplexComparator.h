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

#include "../util/Comparable.h"
#include <complex>

 /**
  * Comparator for std::complex<double> Numbers.
  *
  */
class Complex_Comparator : public Comparable<std::complex<double>>
{
public:

	/** Compare two complex numbers, using real ordering as the primary sort order and
	 * imaginary ordering as the secondary sort order.
	 * @param o1 first complex number
	 * @param o2 second complex number
	 * @return a negative value if o1 real part is less than o2 real part
	 * or if real parts are equal and o1 imaginary part is less than o2 imaginary part
	 */
	 //override
	double compare(const std::complex<double>& o1, const std::complex<double>& o2) override
	{
		auto cR = o1.real() - o2.real();
		return cR == 0
			? o1.imag() - o2.imag()
			: cR;
	}
};