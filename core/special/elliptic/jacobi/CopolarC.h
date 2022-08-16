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
#include "CopolarN.h"

 /** Copolar trio with pole at point c in Glaisher’s Notation.
  * <p>
  * This is a container for the three subsidiary Jacobi elliptic functions
  * {@code dc(u|m)}, {@code nc(u|m)}, and {@code sc(u|m)}.
  * </p>
  * @since 2.0
  */
class Copolar_C
{
private:

	/** Value of the dc function. */
	const double my_dc;

	/** Value of the nc function. */
	const double my_nc;

	/** Value of the sc function. */
	const double my_sc;

public:
	/** Simple constructor.
	 * @param trio_n copolar trio with pole at point n in Glaisher’s Notation
	 */
	Copolar_C(const Copolar_N& trio_n)
		:
		my_nc{ 1.0 / trio_n.cn() },
		my_sc{ my_nc * trio_n.sn() },
		my_dc{ my_nc * trio_n.dn() }
	{};

	/** Get the value of the dc function.
	 * @return dc(u|m)
	 */
	double dc() const
	{
		return my_dc;
	}

	/** Get the value of the nc function.
	 * @return nc(u|m)
	 */
	double nc() const
	{
		return my_nc;
	}

	/** Get the value of the sc function.
	 * @return sc(u|m)
	 */
	double sc() const
	{
		return my_sc;
	}
};