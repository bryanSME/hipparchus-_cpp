#pragma once
/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
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

 /*
  * This is not the original file distributed by the Apache Software Foundation
  * It has been modified by the Hipparchus project
  */

  //package org.hipparchus.random;

  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Sin_Cos;

  /**
   * Base class with default implementations for common methods.
   */
virtual class Base_Random_Generator : public Random_Generator
{
private:
	/** Next gaussian. */
	double my_next_gaussian = std::numeric_limits<double>::quiet_NaN();

public:
	/** {@inherit_doc} */
	//override
	void set_seed(const int& seed)
	{
		set_seed(new std::vector<int>{ seed });
	}

	/** {@inherit_doc} */
	//override
	void set_seed(const long& seed)
	{
		set_seed(new std::vector<int>{ static_cast<int>((seed >> > 32), static_cast<int>((seed & 0xffffffffL) });
	}

	/** {@inherit_doc} */
	//override
	int next_int(const int& n)
	{
		if (n <= 0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL_BOUND_EXCLUDED, n, 0);
		}

		if ((n & -n) == n)
		{
			return static_cast<int>(((n * static_cast<long>((next_int() >> > 1)) >> 31);
		}
		int bits;
		int val;
		do
		{
			bits = next_int() >> > 1;
			val = bits % n;
		} while (bits - val + (n - 1) < 0);
		return val;
	}

	/** {@inherit_doc} */
	//override
	long next_long(const long& n)
	{
		if (n <= 0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL_BOUND_EXCLUDED, n, 0);
		}

		long bits;
		long val;
		do
		{
			bits = next_long() >> > 1;
			val = bits % n;
		} while (bits - val + (n - 1) < 0);
		return val;
	}

	/** {@inherit_doc} */
	//override
	double next_gaussian()
	{
		const double random;
		if (std::isnan(my_next_gaussian))
		{
			// generate a pair of gaussian numbers
			const double x = next_double();
			const double y = next_double();
			const double& alpha = 2 * std::numbers::pi * x;
			const double r = std::sqrt(-2 * std::log(y));
			const Sin_Cos sc_alpha = Sin_Cos(alpha);
			random = r * sc_alpha.cos();
			my_next_gaussian = r * sc_alpha.sin();
		}
		else
		{
			// use the second element of the pair already generated
			random = my_next_gaussian;
			my_next_gaussian = std::numeric_limits<double>::quiet_NaN();
		}

		return random;
	}

	/** {@inherit_doc} */
	//override
	std::string to_string() const
	{
		return get_class().get_name();
	}

protected:
	/**
	 * Clears the cache used by the default implementation of
	 * {@link #next_gaussian}.
	 */
	void clear_cache()
	{
		my_next_gaussian = std::numeric_limits<double>::quiet_NaN();
	}
}
