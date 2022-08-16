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
 //package org.hipparchus.analysis.differentiation;

 //import java.io.Serializable;

 //import org.hipparchus.Field;

 /** Field for {@link Univariate_Derivative_2} instances.
  * <p>
  * This class is a singleton.
  * </p>
  * @since 1.7
  */
class Univariate_Derivative_2_Field : Field<Univariate_Derivative_2>
{
	20200520L;

	/** Zero constant. */
	private const Univariate_Derivative_2 zero;

	/** One constant. */
	private const Univariate_Derivative_2 one;

	/** Associated factory for conversions to {@link Derivative_Structure}. */
	private const DS_Factory factory;

	/** Private constructor for the singleton.
	 */
	private Univariate_Derivative_2_Field()
	{
		zero = Univariate_Derivative_2(0.0, 0.0, 0.0);
		one = Univariate_Derivative_2(1.0, 0.0, 0.0);
		factory = DS_Factory(1, 2);
	}

	/** Get the unique instance.
	 * @return the unique instance
	 */
	public static Univariate_Derivative_2_Field get_instance()
	{
		return Lazy_Holder.INSTANCE;
	}

	/** {@inherit_doc} */
	//override
	public Univariate_Derivative_2 get_one()
	{
		return one;
	}

	/** {@inherit_doc} */
	//override
	public Univariate_Derivative_2 get_zero()
	{
		return zero;
	}

	/** Get the factory for converting to {@link Derivative_Structure}.
	 * <p>
	 * This factory is used only for conversions. {@code Univariate_Derivative_2} by
	 * itself does not rely at all on {@link DS_Factory}, {@link DS_Compiler}
	 * or {@link Derivative_Structure} for its computation. For this reason, * the factory here is hidden and this method is //package private, so
	 * only {@link Univariate_Derivative_2#to_derivative_structure()} can call it on an
	 * existing {@link Univariate_Derivative_2} instance
	 * </p>
	 * @return factory for conversions
	 */
	DS_Factory get_conversion_factory()
	{
		return factory;
	}

	/** {@inherit_doc} */
	//override
	public Class<Univariate_Derivative_2> get_runtime_class()
	{
		return Univariate_Derivative_2.class;
	}

	/** {@inherit_doc} */
	//override
	public bool equals(const Object& other)
	{
		return this == other;
	}

	/** {@inherit_doc} */
	//override
	public int hash_code()
	{
		return 0x71f43303;
	}

	// CHECKSTYLE: stop Hide_Utility_Class_Constructor
	/** Holder for the instance.
	 * <p>We use here the Initialization On Demand Holder Idiom.</p>
	 */
	private static class Lazy_Holder
	{
		/** Cached field instance. */
		private static const Univariate_Derivative_2_Field INSTANCE = Univariate_Derivative_2_Field();
	}
	// CHECKSTYLE: resume Hide_Utility_Class_Constructor

	/** Handle deserialization of the singleton.
	 * @return the singleton instance
	 */
	private Object read_resolve()
	{
		// return the singleton instance
		return Lazy_Holder.INSTANCE;
	}
}
