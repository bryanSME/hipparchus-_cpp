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

 //import java.util.Hash_Map;
 //import java.util.Map;
#include <unordered_map>
#include <type_traits>
#include "../../CalculusFieldElement.hpp"

//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.Field;

/** Field for {@link Field_Univariate_Derivative_1} instances.
 * @param <T> the type of the function parameters and value
 * @since 1.7
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
class Field_Univariate_Derivative_1Field : public Field<Field_Univariate_Derivative_1<T>>
{
private:
	/** Cached fields. */
	static const std::map<Field<>, Field_Univariate_Derivative_1Field<>> CACHE = Hash_Map<>();

	/** Zero constant. */
	const Field_Univariate_Derivative_1<T> my_zero;

	/** One constant. */
	const Field_Univariate_Derivative_1<T> my_one;

	/** Associated factory for conversions to {@link Field_Derivative_Structure}. */
	const FDS_Factory<T> factory;

	/** Private constructor for populating the cache.
	 * @param value_field field for the function parameters and value
	 */
	Field_Univariate_Derivative_1Field(const Field<T> value_field)
		:
		my_zero{ Field_Univariate_Derivative_1<>(value_field.get_zero(), value_field.get_zero()) },
		my_one{ Field_Univariate_Derivative_1<>(value_field.get_one(), value_field.get_zero()) },
		my_factory{ FDS_Factory<>(value_field, 1, 1) }
	{
	}

public:

	/** Get the univariate derivative field corresponding to a value field.
	 * @param value_field field for the function parameters and value
	 * @param <T> the type of the function parameters and value
	 * @return univariate derivative field
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	static Field_Univariate_Derivative_1Field<T> get_univariate_derivative1_field(const Field<T> value_field)
	{
		synchronized(CACHE)
		{
			Field_Univariate_Derivative_1Field<> cached = CACHE.get(value_field);
			if (cached == NULL)
			{
				cached = Field_Univariate_Derivative_1Field<>(value_field);
				CACHE.put(value_field, cached);
			}
			//@Suppress_Warnings("unchecked")
			const Field_Univariate_Derivative_1Field<T> t_cached = (Field_Univariate_Derivative_1Field<T>) cached;
			return t_cached;
		}
	}

	/** {@inherit_doc} */
	//override
	Field_Univariate_Derivative_1<T> get_one() const
	{
		return my_one;
	}

	/** {@inherit_doc} */
	//override
	Field_Univariate_Derivative_1<T> get_zero() const
	{
		return my_zero;
	}

	/** Get the factory for converting to {@link Derivative_Structure}.
	 * <p>
	 * This factory is used only for conversions. {@code Univariate_Derivative_1} by
	 * itself does not rely at all on {@link DS_Factory}, {@link DS_Compiler}
	 * or {@link Derivative_Structure} for its computation. For this reason, * the factory here is hidden and this method is //package private, so
	 * only {@link Univariate_Derivative_1#to_derivative_structure()} can call it on an
	 * existing {@link Univariate_Derivative_1} instance
	 * </p>
	 * @return factory for conversion
	 */
	FDS_Factory<T> get_conversion_factory() const
	{
		return my_factory;
	}

	/** {@inherit_doc} */
	//@Suppress_Warnings("unchecked")
	//override
	Class<Field_Univariate_Derivative_1<T>> get_runtime_class()
	{
		return (Class<Field_Univariate_Derivative_1<T>>) zero.get_class();
	}

	/** {@inherit_doc} */
	//override
	bool equals(const Object& other) const
	{
		return this == other;
	}

	/** {@inherit_doc} */
	//override
	int hash_code() const
	{
		return 0x712c0fc7;
	}
};