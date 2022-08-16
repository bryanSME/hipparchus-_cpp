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

  //package org.hipparchus.complex;

  //import java.io.Serializable;
#include "std::complex<double>.h"
//import org.hipparchus.Field;

/**
 * Representation of the complex numbers field.
 * <p>
 * This class is a singleton.
 *
 * @see std::complex<double>
 */
class std::complex<double>_Field : Field<std::complex<double>>
{
private:
	/** Private constructor for the singleton.
	 */
	std::complex<double>_Field() = default;

	// CHECKSTYLE: stop Hide_Utility_Class_Constructor
	/** Holder for the instance.
	 * <p>We use here the Initialization On Demand Holder Idiom.</p>
	 */
	static class Lazy_Holder
	{
		/** Cached field instance. */
		static const std::complex<double>_Field INSTANCE = std::complex<double>_Field();
	}
	// CHECKSTYLE: resume Hide_Utility_Class_Constructor

	/** Handle deserialization of the singleton.
	 * @return the singleton instance
	 */
	Object read_resolve()
	{
		// return the singleton instance
		return Lazy_Holder.INSTANCE;
	}

public:
	/** Get the unique instance.
	 * @return the unique instance
	 */
	static std::complex<double>_Field get_instance()
	{
		return Lazy_Holder.INSTANCE;
	}

	/** {@inherit_doc} */
	//override
	std::complex<double> get_one()
	{
		return std::complex<double>.ONE;
	}

	/** {@inherit_doc} */
	//override
	std::complex<double> get_zero()
	{
		return std::complex<double>.ZERO;
	}

	/** {@inherit_doc} */
	//override
	Class<std::complex<double>> get_runtime_class()
	{
		return std::complex<double>.class;
	}

	/** {@inherit_doc} */
	//override
	bool equals(const Object& other)
	{
		return this == other;
	}

	/** {@inherit_doc} */
	//override
	int hash_code()
	{
		return 0x49250ae5;
	}
};