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
  //package org.hipparchus.linear;

  //import java.io.Serializable;
  //import java.util.Arrays;

  //import org.hipparchus.Field;
  //import org.hipparchus.Field_Element;
  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.exception.Math_Runtime_Exception;
  //import org.hipparchus.exception.;
  //import org.hipparchus.util.Math_Arrays;
  //import org.hipparchus.util.Math_Utils;
#include <type_traits>
#include <exception>
#include <vector>
#include "../FieldElement.h"
#include "../Field.h"

  /**
   * This class : the {@link Field_Vector} interface with a {@link Field_Element} array.
   * @param <T> the type of the field elements
   */
template<typename T, typename std::enable_if<std::is_base_of<Field_Element<T>, T>::value>::type* = nullptr>
class Array_Field_Vector : public Field_Vector<T>
{
private:
	/** Entries of the vector. */
	std::vector<T> my_data;

	/** Field to which the elements belong. */
	const Field<T> my_field;

	/**
	 * Check if an index is valid.
	 *
	 * @param index Index to check.
	 * @exception  if the index is not valid.
	 */
	void check_index(const int& index)
	{
		if (index < 0 || index >= get_dimension())
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::INDEX, index, 0, get_dimension() - 1);
		}
	}

	/**
	 * Checks that the indices of a subvector are valid.
	 *
	 * @param start the index of the first entry of the subvector
	 * @param end the index of the last entry of the subvector (inclusive)
	 * @ if {@code start} of {@code end} are not valid
	 * @ if {@code end < start}
	 */
	void check_indices(const int& start, const int& end)
	{
		const int dim = get_dimension();
		if ((start < 0) || (start >= dim))
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::INDEX, start, 0, dim - 1);
		}
		if ((end < 0) || (end >= dim))
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::INDEX, end, 0, dim - 1);
		}
		if (end < start)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::INITIAL_ROW_AFTER_FINAL_ROW, end, start, false);
		}
	}

protected:
	/**
	 * Check if instance and specified vectors have the same dimension.
	 * @param v vector to compare instance with
	 * @exception  if the vectors do not
	 * have the same dimensions
	 */
	void check_vector_dimensions(const Field_Vector<T>& v)
	{
		check_vector_dimensions(v.get_bimension());
	}

	/**
	 * Check if instance dimension is equal to some expected value.
	 *
	 * @param n Expected dimension.
	 * @ if the dimension is not equal to the
	 * size of {@code this} vector.
	 */
	void check_vector_dimensions(const int& n)
	{
		if (data.size() != n)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, data.size(), n);
		}
	}

public:
	/**
	 * Build a 0-length vector.
	 * Zero-length vectors may be used to initialize construction of vectors
	 * by data gathering. We start with zero-length and use either the {@link
	 * #Array_Field_Vector(Field_Vector, Field_Vector)} constructor
	 * or one of the {@code append} methods ({@link #add(Field_Vector)} or
	 * {@link #append(Array_Field_Vector)}) to gather data into this vector.
	 *
	 * @param field field to which the elements belong
	 */
	Array_Field_Vector(const Field<T>& field)
	{
		Array_Field_Vector(field, 0);
	}

	/**
	 * Construct a vector of zeroes.
	 *
	 * @param field Field to which the elements belong.
	 * @param size Size of the vector.
	 */
	Array_Field_Vector(const Field<T>& field, const int& size)
		:
		my_field{ field },
		my_data{ Math_Arrays::build_array(field, size) }
	{
	}

	/**
	 * Construct a vector with preset values.
	 *
	 * @param size Size of the vector.
	 * @param preset All entries will be set with this value.
	 */
	Array_Field_Vector(const int& size, const T& preset)
	{
		Array_Field_Vector(preset.get_field(), size);
		Arrays.fill(data, preset);
	}

	/**
	 * Construct a vector from an array, copying the input array.
	 * This constructor needs a non-empty {@code d} array to retrieve
	 * the field from its first element. This implies it cannot build
	 * 0 length vectors. To build vectors from any size, one should
	 * use the {@link #Array_Field_Vector(Field, Field_Element[])} constructor.
	 *
	 * @param d Array.
	 * @ if {@code d} is {@code NULL}.
	 * @ if {@code d} is empty.
	 * @see #Array_Field_Vector(Field, Field_Element[])
	 */
	Array_Field_Vector(const std::vector<T>& d)
	{
		//Math_Utils::check_not_null(d);
		try
		{
			my_field = d[0].get_field();
			my_data = d.clone();
		}
		catch (Array_indexOutOfboundsException e)
		{
			throw std::exception("not implemented");
			//throw (e, hipparchus::exception::Localized_Core_Formats_Type::VECTOR_MUST_HAVE_AT_LEAST_ONE_ELEMENT);
		}
	}

	/**
	 * Construct a vector from an array, copying the input array.
	 *
	 * @param field Field to which the elements belong.
	 * @param d Array.
	 * @ if {@code d} is {@code NULL}.
	 * @see #Array_Field_Vector(Field_Element[])
	 */
	Array_Field_Vector(const Field<T>& field, const std::vector<T>& d)
		:
		my_field{ field },
		my_data{ d }
	{
		//Math_Utils::check_not_null(d);
	}

	/**
	 * Create a Array_Field_Vector using the input array as the underlying
	 * data array.
	 * If an array is built specially in order to be embedded in a
	 * Array_Field_Vector and not used directly, the {@code copy_array} may be
	 * set to {@code false}. This will prevent the copying and improve
	 * performance as no array will be built and no data will be copied.
	 * This constructor needs a non-empty {@code d} array to retrieve
	 * the field from its first element. This implies it cannot build
	 * 0 length vectors. To build vectors from any size, one should
	 * use the {@link #Array_Field_Vector(Field, Field_Element[], bool)}
	 * constructor.
	 *
	 * @param d Data for the vector.
	 * @param copy_array If {@code true}, the input array will be copied, * otherwise it will be referenced.
	 * @ if {@code d} is {@code NULL}.
	 * @ if {@code d} is empty.
	 * @see #Array_Field_Vector(Field_Element[])
	 * @see #Array_Field_Vector(Field, Field_Element[], bool)
	 */
	Array_Field_Vector(const std::vector<T>& d, bool copy_array)
	{
		//Math_Utils::check_not_null(d);
		if (d.size() == 0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::VECTOR_MUST_HAVE_AT_LEAST_ONE_ELEMENT);
		}
		my_field = d[0].get_field();
		my_data = copy_array
			? d.clone() 
			: d;
	}

	/**
	 * Create a Array_Field_Vector using the input array as the underlying
	 * data array.
	 * If an array is built specially in order to be embedded in a
	 * Array_Field_Vector and not used directly, the {@code copy_array} may be
	 * set to {@code false}. This will prevent the copying and improve
	 * performance as no array will be built and no data will be copied.
	 *
	 * @param field Field to which the elements belong.
	 * @param d Data for the vector.
	 * @param copy_array If {@code true}, the input array will be copied, * otherwise it will be referenced.
	 * @ if {@code d} is {@code NULL}.
	 * @see #Array_Field_Vector(Field_Element[], bool)
	 */
	Array_Field_Vector(const Field<T>& field, const std::vector<T>& d, bool copy_array)
	{
		//Math_Utils::check_not_null(d);
		my_field = field;
		data = copy_array
			? d.clone()
			: d;
	}

	/**
	 * Construct a vector from part of a array.
	 *
	 * @param d Array.
	 * @param pos Position of the first entry.
	 * @param size Number of entries to copy.
	 * @ if {@code d} is {@code NULL}.
	 * @ if the size of {@code d} is less
	 * than {@code pos + size}.
	 */
	Array_Field_Vector(const std::vector<T>& d, int pos, int size)

	{
		//Math_Utils::check_not_null(d);
		if (d.size() < pos + size)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_LARGE, pos + size, d.size());
		}
		field = d[0].get_field();
		data = Math_Arrays::build_array(field, size);
		System.arraycopy(d, pos, data, 0, size);
	}

	/**
	 * Construct a vector from part of a array.
	 *
	 * @param field Field to which the elements belong.
	 * @param d Array.
	 * @param pos Position of the first entry.
	 * @param size Number of entries to copy.
	 * @ if {@code d} is {@code NULL}.
	 * @ if the size of {@code d} is less
	 * than {@code pos + size}.
	 */
	Array_Field_Vector(Field<T> field, const std::vector<T>& d, int pos, int size)
	{
		//Math_Utils::check_not_null(d);
		if (d.size() < pos + size)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_LARGE, pos + size, d.size());
		}
		this.field = field;
		data = Math_Arrays::build_array(field, size);
		System.arraycopy(d, pos, data, 0, size);
	}

	/**
	 * Construct a vector from another vector, using a deep copy.
	 *
	 * @param v Vector to copy.
	 * @ if {@code v} is {@code NULL}.
	 */
	Array_Field_Vector(const Field_Vector<T>& v)

	{
		//Math_Utils::check_not_null(v);
		field = v.get_field();
		data = Math_Arrays::build_array(field, v.get_dimension());
		for (int i{}; i < data.size(); ++i)
		{
			data[i] = v.get_entry(i);
		}
	}

	/**
	 * Construct a vector from another vector, using a deep copy.
	 *
	 * @param v Vector to copy.
	 * @ if {@code v} is {@code NULL}.
	 */
	Array_Field_Vector(const Array_Field_Vector<T>& v)
		:
		my_field{ v.get_field() },
		my_data{ v.data.clone() }
	{
		//Math_Utils::check_not_null(v);
	}

	/**
	 * Construct a vector from another vector.
	 *
	 * @param v Vector to copy.
	 * @param deep If {@code true} perform a deep copy, otherwise perform
	 * a shallow copy
	 * @ if {@code v} is {@code NULL}.
	 */
	Array_Field_Vector(Array_Field_Vector<T> v, bool deep)

	{
		//Math_Utils::check_not_null(v);
		field = v.get_field();
		data = deep ? v.data.clone() : v.data;
	}

	/**
	 * Construct a vector by appending one vector to another vector.
	 *
	 * @param v1 First vector (will be put in front of the vector).
	 * @param v2 Second vector (will be put at back of the vector).
	 * @ if {@code v1} or {@code v2} is
	 * {@code NULL}.
	 */
	Array_Field_Vector(Field_Vector<T> v1, Field_Vector<T> v2)

	{
		//Math_Utils::check_not_null(v1);
		//Math_Utils::check_not_null(v2);
		field = v1.get_field();
		const std::vector<T> v1_data =
			(v1 instanceof Array_Field_Vector) ? ((Array_Field_Vector<T>) v1).data : v1.to_array();
		const std::vector<T> v2_data =
			(v2 instanceof Array_Field_Vector) ? ((Array_Field_Vector<T>) v2).data : v2.to_array();
		data = Math_Arrays::build_array(field, v1_data.size() + v2_data.size());
		System.arraycopy(v1_data, 0, data, 0, v1_data.size());
		System.arraycopy(v2_data, 0, data, v1_data.size(), v2_data.size());
	}

	/**
	 * Construct a vector by appending one vector to another vector.
	 *
	 * @param v1 First vector (will be put in front of the vector).
	 * @param v2 Second vector (will be put at back of the vector).
	 * @ if {@code v1} or {@code v2} is
	 * {@code NULL}.
	 */
	Array_Field_Vector(Field_Vector<T> v1, std::vector<T> v2)

	{
		//Math_Utils::check_not_null(v1);
		//Math_Utils::check_not_null(v2);
		field = v1.get_field();
		const std::vector<T> v1_data =
			(v1 instanceof Array_Field_Vector) ? ((Array_Field_Vector<T>) v1).data : v1.to_array();
		data = Math_Arrays::build_array(field, v1_data.size() + v2.size());
		System.arraycopy(v1_data, 0, data, 0, v1_data.size());
		System.arraycopy(v2, 0, data, v1_data.size(), v2.size());
	}

	/**
	 * Construct a vector by appending one vector to another vector.
	 *
	 * @param v1 First vector (will be put in front of the vector).
	 * @param v2 Second vector (will be put at back of the vector).
	 * @ if {@code v1} or {@code v2} is
	 * {@code NULL}.
	 */
	Array_Field_Vector(std::vector<T> v1, Field_Vector<T> v2)

	{
		//Math_Utils::check_not_null(v1);
		//Math_Utils::check_not_null(v2);
		field = v2.get_field();
		const std::vector<T> v2_data =
			(v2 instanceof Array_Field_Vector) ? ((Array_Field_Vector<T>) v2).data : v2.to_array();
		data = Math_Arrays::build_array(field, v1.size() + v2_data.size());
		System.arraycopy(v1, 0, data, 0, v1.size());
		System.arraycopy(v2_data, 0, data, v1.size(), v2_data.size());
	}

	/**
	 * Construct a vector by appending one vector to another vector.
	 * This constructor needs at least one non-empty array to retrieve
	 * the field from its first element. This implies it cannot build
	 * 0 length vectors. To build vectors from any size, one should
	 * use the {@link #Array_Field_Vector(Field, Field_Element[], Field_Element[])}
	 * constructor.
	 *
	 * @param v1 First vector (will be put in front of the vector).
	 * @param v2 Second vector (will be put at back of the vector).
	 * @ if {@code v1} or {@code v2} is
	 * {@code NULL}.
	 * @ if both arrays are empty.
	 * @see #Array_Field_Vector(Field, Field_Element[], Field_Element[])
	 */
	Array_Field_Vector(std::vector<T> v1, std::vector<T> v2)

	{
		//Math_Utils::check_not_null(v1);
		//Math_Utils::check_not_null(v2);
		if (v1.size() + v2.size() == 0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::VECTOR_MUST_HAVE_AT_LEAST_ONE_ELEMENT);
		}
		data = Math_Arrays::build_array(v1[0].get_field(), v1.size() + v2.size());
		System.arraycopy(v1, 0, data, 0, v1.size());
		System.arraycopy(v2, 0, data, v1.size(), v2.size());
		field = data[0].get_field();
	}

	/**
	 * Construct a vector by appending one vector to another vector.
	 *
	 * @param field Field to which the elements belong.
	 * @param v1 First vector (will be put in front of the vector).
	 * @param v2 Second vector (will be put at back of the vector).
	 * @ if {@code v1} or {@code v2} is
	 * {@code NULL}.
	 * @ if both arrays are empty.
	 * @see #Array_Field_Vector(Field_Element[], Field_Element[])
	 */
	Array_Field_Vector(Field<T> field, std::vector<T> v1, std::vector<T> v2)

	{
		//Math_Utils::check_not_null(v1);
		//Math_Utils::check_not_null(v2);
		if (v1.size() + v2.size() == 0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::VECTOR_MUST_HAVE_AT_LEAST_ONE_ELEMENT);
		}
		data = Math_Arrays::build_array(field, v1.size() + v2.size());
		System.arraycopy(v1, 0, data, 0, v1.size());
		System.arraycopy(v2, 0, data, v1.size(), v2.size());
		this.field = field;
	}

	/** {@inherit_doc} */
	//override
	Field<T> get_field()
	{
		return field;
	}

	/** {@inherit_doc} */
	//override
	Field_Vector<T> copy()
	{
		return Array_Field_Vector<T>(this, true);
	}

	/** {@inherit_doc} */
	//override
	Field_Vector<T> add(const Field_Vector<T>& v)
	{
		if (dynamic_cast<const Array_Field_Vector*>(*v) != nullptr)
		{
			return add((Array_Field_Vector<T>) v);
		}
		check_vector_dimensions(v);
		auto out = Math_Arrays::build_array(field, data.size());
		for (int i{}; i < data.size(); i++)
		{
			out[i] = data[i].add(v.get_entry(i));
		}
		return Array_Field_Vector<T>(field, out, false);
	}

	/**
	 * Compute the sum of {@code this} and {@code v}.
	 * @param v vector to be added
	 * @return {@code this + v}
	 * @ if {@code v} is not the same size as
	 * {@code this}
	 */
	Array_Field_Vector<T> add(const Array_Field_Vector<T>& v)

	{
		check_vector_dimensions(v.data.size());
		auto out = Math_Arrays::build_array(field, data.size());
		for (int i{}; i < data.size(); i++)
		{
			out[i] = data[i].add(v.data[i]);
		}
		return Array_Field_Vector<T>(field, out, false);
	}

	/** {@inherit_doc} */
	//override
	Field_Vector<T> subtract(const Field_Vector<T>& v)
	{
		if (dynamic_cast<const Array_Field_Vector*>(*v) != nullptr)
		{
			return subtract((Array_Field_Vector<T>) v);
		}
		check_vector_dimensions(v);
		auto out = Math_Arrays::build_array(field, data.size());
		for (int i{}; i < data.size(); i++)
		{
			out[i] = data[i].subtract(v.get_entry(i));
		}
		return Array_Field_Vector<T>(field, out, false);
	}

	/**
	 * Compute {@code this} minus {@code v}.
	 * @param v vector to be subtracted
	 * @return {@code this - v}
	 * @ if {@code v} is not the same size as
	 * {@code this}
	 */
	Array_Field_Vector<T> subtract(const Array_Field_Vector<T>& v)

	{
		check_vector_dimensions(v.data.size());
		auto out = Math_Arrays::build_array(field, data.size());
		for (int i{}; i < data.size(); i++)
		{
			out[i] = data[i].subtract(v.data[i]);
		}
		return Array_Field_Vector<T>(field, out, false);
	}

	/** {@inherit_doc} */
	//override
	Field_Vector<T> map_add(const T& d)
	{
		auto out = Math_Arrays::build_array(field, data.size());
		for (int i{}; i < data.size(); i++)
		{
			out[i] = data[i].add(d);
		}
		return Array_Field_Vector<T>(field, out, false);
	}

	/** {@inherit_doc} */
	//override
	Field_Vector<T> map_add_to_self(const T& d)
	{
		for (int i{}; i < data.size(); i++)
		{
			data[i] = data[i].add(d);
		}
		return *this;
	}

	/** {@inherit_doc} */
	//override
	Field_Vector<T> map_subtract(const T& d)
	{
		auto out = Math_Arrays::build_array(field, data.size());
		for (int i{}; i < data.size(); i++)
		{
			out[i] = data[i].subtract(d);
		}
		return Array_Field_Vector<T>(field, out, false);
	}

	/** {@inherit_doc} */
	//override
	Field_Vector<T> map_subtract_to_self(const T& d)
	{
		for (int i{}; i < data.size(); i++)
		{
			data[i] = data[i].subtract(d);
		}
		return *this;
	}

	/** {@inherit_doc} */
	//override
	Field_Vector<T> map_multiply(const T& d)
	{
		auto out = Math_Arrays::build_array(field, data.size());
		for (int i{}; i < data.size(); i++)
		{
			out[i] = data[i].multiply(d);
		}
		return Array_Field_Vector<T>(field, out, false);
	}

	/** {@inherit_doc} */
	//override
	Field_Vector<T> map_multiply_to_self(const T& d)
	{
		for (int i{}; i < data.size(); i++)
		{
			data[i] = data[i].multiply(d);
		}
		return *this;
	}

	/** {@inherit_doc} */
	//override
	Field_Vector<T> map_divide(const T& d)
		
	{
		//Math_Utils::check_not_null(d);
		auto out = Math_Arrays::build_array(field, data.size());
		for (int i{}; i < data.size(); i++)
		{
			out[i] = data[i].divide(d);
		}
		return Array_Field_Vector<T>(field, out, false);
	}

		/** {@inherit_doc} */
		//override
		Field_Vector<T> map_divide_to_self(const T& d)
		
	{
		//Math_Utils::check_not_null(d);
		for (int i{}; i < data.size(); i++)
		{
			data[i] = data[i].divide(d);
		}
		return *this;
	}

		/** {@inherit_doc} */
		//override
		Field_Vector<T> map_inv()
	{
		auto out = Math_Arrays::build_array(field, data.size());
		const T one = field.get_one();
		for (int i{}; i < data.size(); i++)
		{
			try
			{
				out[i] = one.divide(data[i]);
			}
			catch (const Math_Runtime_Exception e)
			{
				throw std::exception("not implemented");
				//throw Math_Runtime_Exception(e, hipparchus::exception::Localized_Core_Formats_Type::INDEX, i);
			}
		}
		return Array_Field_Vector<T>(field, out, false);
	}

	/** {@inherit_doc} */
	//override
	Field_Vector<T> map_inv_to_self()
	{
		const T one = field.get_one();
		for (int i{}; i < data.size(); i++)
		{
			try
			{
				data[i] = one.divide(data[i]);
			}
			catch (const Math_Runtime_Exception e)
			{
				throw std::exception("not implemented");
				//throw Math_Runtime_Exception(e, hipparchus::exception::Localized_Core_Formats_Type::INDEX, i);
			}
		}
		return *this;
	}

	/** {@inherit_doc} */
	//override
	Field_Vector<T> ebe_multiply(const Field_Vector<T>& v)
	{
		if (dynamic_cast<const Array_Field_Vector*>(*v) != nullptr)
		{
			return ebe_multiply((Array_Field_Vector<T>) v);
		}
		check_vector_dimensions(v);
		auto out = Math_Arrays::build_array(field, data.size());
		for (int i{}; i < data.size(); i++)
		{
			out[i] = data[i].multiply(v.get_entry(i));
		}
		return Array_Field_Vector<T>(field, out, false);
	}

	/**
	 * Element-by-element multiplication.
	 * @param v vector by which instance elements must be multiplied
	 * @return a vector containing {@code this[i] * v[i]} for all {@code i}
	 * @ if {@code v} is not the same size as
	 * {@code this}
	 */
	Array_Field_Vector<T> ebe_multiply(const Array_Field_Vector<T>& v)

	{
		check_vector_dimensions(v.data.size());
		auto out = Math_Arrays::build_array(field, data.size());
		for (int i{}; i < data.size(); i++)
		{
			out[i] = data[i].multiply(v.data[i]);
		}
		return Array_Field_Vector<T>(field, out, false);
	}

	/** {@inherit_doc} */
	//override
	Field_Vector<T> ebe_divide(const Field_Vector<T>& v)
	{
		if (dynamic_cast<const Array_Field_Vector*>(*v) != nullptr)
		{
			return ebe_divide((Array_Field_Vector<T>) v);
		}
		check_vector_dimensions(v);
		auto out = Math_Arrays::build_array(field, data.size());
		for (int i{}; i < data.size(); i++)
		{
			try
			{
				out[i] = data[i].divide(v.get_entry(i));
			}
			catch (const Math_Runtime_Exception e)
			{
				throw std::exception("not implemented");
				//throw Math_Runtime_Exception(e, hipparchus::exception::Localized_Core_Formats_Type::INDEX, i);
			}
		}
		return Array_Field_Vector<T>(field, out, false);
	}

	/**
	 * Element-by-element division.
	 * @param v vector by which instance elements must be divided
	 * @return a vector containing {@code this[i] / v[i]} for all {@code i}
	 * @ if {@code v} is not the same size as
	 * {@code this}
	 * @Math_Runtime_Exception if one entry of {@code v} is zero.
	 */
	Array_Field_Vector<T> ebe_divide(const Array_Field_Vector<T>& v)
	{
		check_vector_dimensions(v.data.size());
		auto out = Math_Arrays::build_array(field, data.size());
		for (int i{}; i < data.size(); i++)
		{
			try
			{
				out[i] = data[i].divide(v.data[i]);
			}
			catch (const Math_Runtime_Exception e)
			{
				throw std::exception("not implemented");
					//throw Math_Runtime_Exception(e, hipparchus::exception::Localized_Core_Formats_Type::INDEX, i);
			}
		}
			return Array_Field_Vector<T>(field, out, false);
	}

	/**
		* Returns a reference to the underlying data array.
		* <p>Does not make a fresh copy of the underlying data.</p>
		* @return array of entries
		*/
	std::vector<T> get_data_ref() const
	{
		return my_data; // NOPMD - returning an internal array is intentional and documented here
	}

	/** {@inherit_doc} */
	//override
	T dot_product(const Field_Vector<T>& v)
	{
		if (dynamic_cast<const Array_Field_Vector*>(*v) != nullptr)
		{
			return dot_product((Array_Field_Vector<T>) v);
		}
		check_vector_dimensions(v);
		T dot = field.get_zero();
		for (int i{}; i < data.size(); i++)
		{
			dot = dot.add(data[i].multiply(v.get_entry(i)));
		}
		return dot;
	}

	/**
	 * Compute the dot product.
	 * @param v vector with which dot product should be computed
	 * @return the scalar dot product of {@code this} and {@code v}
	 * @ if {@code v} is not the same size as
	 * {@code this}
	 */
	T dot_product(const Array_Field_Vector<T>& v)
	{
		check_vector_dimensions(v.data.size());
		T dot = field.get_zero();
		for (int i{}; i < data.size(); i++)
		{
			dot = dot.add(data[i].multiply(v.data[i]));
		}
		return dot;
	}

	/** {@inherit_doc} */
	//override
	Field_Vector<T> projection(const Field_Vector<T>& v)
	{
		return v.map_multiply(dot_product(v).divide(v.dot_product(v)));
	}

	/** Find the orthogonal projection of this vector onto another vector.
		* @param v vector onto which {@code this} must be projected
		* @return projection of {@code this} onto {@code v}
		* @ if {@code v} is not the same size as
		* {@code this}
		* @Math_Runtime_Exception if {@code v} is the NULL vector.
		*/
	Array_Field_Vector<T> projection(const Array_Field_Vector<T>& v)
	{
		return (Array_Field_Vector<T>) v.map_multiply(dot_product(v).divide(v.dot_product(v)));
	}

	/** {@inherit_doc} */
	//override
	Field_Matrix<T> outer_product(const Field_Vector<T>& v)
	{
		if (dynamic_cast<const Array_Field_Vector*>(*v) != nullptr)
		{
			return outer_product((Array_Field_Vector<T>) v);
		}
		const int m = data.size();
		const int n = v.get_dimension();
		auto out = Array2DRowField_Matrix<>(field, m, n);
		for (int i{}; i < m; i++)
		{
			for (int j{}; j < n; j++)
			{
				out.set_entry(i, j, data[i].multiply(v.get_entry(j)));
			}
		}
		return out;
	}

	/**
	 * Compute the outer product.
	 * @param v vector with which outer product should be computed
	 * @return the matrix outer product between instance and v
	 */
	Field_Matrix<T> outer_product(const Array_Field_Vector<T>& v)
	{
		const int m = my_data.size();
		const int n = v.data.size();
		auto out = Array2DRowField_Matrix<>(field, m, n);
		for (int i{}; i < m; i++)
		{
			for (int j{}; j < n; j++)
			{
				out.set_entry(i, j, data[i].multiply(v.data[j]));
			}
		}
		return out;
	}

	/** {@inherit_doc} */
	//override
	T get_entry(const int& index) const
	{
		return my_data[index];
	}

	/** {@inherit_doc} */
	//override
	int get_dimension() const
	{
		return my_data.size();
	}

	/** {@inherit_doc} */
	//override
	Field_Vector<T> append(const Field_Vector<T>& v)
	{
		if (dynamic_cast<const Array_Field_Vector*>(*v) != nullptr)
		{
			return append((Array_Field_Vector<T>) v);
		}
		return Array_Field_Vector<T>(*this, Array_Field_Vector<T>(v));
	}

	/**
	 * Construct a vector by appending a vector to this vector.
	 * @param v vector to append to this one.
	 * @return a vector
	 */
	Array_Field_Vector<T> append(const Array_Field_Vector<T>& v)
	{
		return Array_Field_Vector<T>(*this, v);
	}

	/** {@inherit_doc} */
	//override
	Field_Vector<T> append(T in)
	{
		auto out = Math_Arrays::build_array(field, data.size() + 1);
		System.arraycopy(data, 0, out, 0, data.size());
		out[data.size()] = in;
		return Array_Field_Vector<T>(field, out, false);
	}

	/** {@inherit_doc} */
	//override
	Field_Vector<T> get_sub_vector(const int& index, const int& n)
	{
		if (n < 0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_OF_ELEMENTS_SHOULD_BE_POSITIVE, n);
		}
		auto out = Array_Field_Vector<>(field, n);
		try
		{
			System.arraycopy(data, index, out.data, 0, n);
		}
		catch (Index_Out_Of_Bounds_Exception e)
		{
			check_index(index);
			check_index(index + n - 1);
		}
		return out;
	}

	/** {@inherit_doc} */
	//override
	void set_entry(const int& index, const T& value)
	{
		try
		{
			data[index] = value;
		}
		catch (Index_Out_Of_Bounds_Exception e)
		{
			check_index(index);
		}
	}

	/** {@inherit_doc} */
	//override
	void set_sub_vector(const int& index, const Field_Vector<T>& v)
	{
		try
		{
			if (dynamic_cast<const Array_Field_Vector*>(*v) != nullptr)
			{
				set(index, (Array_Field_Vector<T>) v);
			}
			else
			{
				for (int i = index; i < index + v.get_dimension(); ++i)
				{
					data[i] = v.get_entry(i - index);
				}
			}
		}
		catch (Index_Out_Of_Bounds_Exception e)
		{
			check_index(index);
			check_index(index + v.get_dimension() - 1);
		}
	}

	/**
	 * Set a set of consecutive elements.
	 *
	 * @param index index of first element to be set.
	 * @param v vector containing the values to set.
	 * @ if the index is invalid.
	 */
	void set(const int& index, const Array_Field_Vector<T>& v)
	{
		try
		{
			System.arraycopy(v.data, 0, data, index, v.data.size());
		}
		catch (Index_Out_Of_Bounds_Exception e)
		{
			check_index(index);
			check_index(index + v.data.size() - 1);
		}
	}

	/** {@inherit_doc} */
	//override
	void set(const T& value)
	{
		Arrays.fill(my_data, value);
	}

	/** {@inherit_doc} */
	//override
	std::vector<T> to_array()
	{
		return my_data.clone();
	}

	/**
	 * Visits (but does not alter) all entries of this vector in default order
	 * (increasing index).
	 *
	 * @param visitor the visitor to be used to process the entries of this
	 * vector
	 * @return the value returned by {@link Field_VectorPreservingVisitor#end()}
	 * at the end of the walk
	 */
	T walk_in_default_order(const Field_VectorPreservingVisitor<T>& visitor)
	{
		const int dim = get_dimension();
		visitor.start(dim, 0, dim - 1);
		for (int i{}; i < dim; i++)
		{
			visitor.visit(i, get_entry(i));
		}
		return visitor.end();
	}

	/**
	 * Visits (but does not alter) some entries of this vector in default order
	 * (increasing index).
	 *
	 * @param visitor visitor to be used to process the entries of this vector
	 * @param start the index of the first entry to be visited
	 * @param end the index of the last entry to be visited (inclusive)
	 * @return the value returned by {@link Field_VectorPreservingVisitor#end()}
	 * at the end of the walk
	 * @ if {@code end < start}.
	 * @ if the indices are not valid.
	 */
	T walk_in_default_order(const Field_VectorPreservingVisitor<T>& visitor, const int& start, const int& end)
	{
		check_indices(start, end);
		visitor.start(get_dimension(), start, end);
		for (int i{ start }; i <= end; i++)
		{
			visitor.visit(i, get_entry(i));
		}
		return visitor.end();
	}

	/**
	 * Visits (but does not alter) all entries of this vector in optimized
	 * order. The order in which the entries are visited is selected so as to
	 * lead to the most efficient implementation; it might depend on the
	 * concrete implementation of this virtual class.
	 *
	 * @param visitor the visitor to be used to process the entries of this
	 * vector
	 * @return the value returned by {@link Field_VectorPreservingVisitor#end()}
	 * at the end of the walk
	 */
	T walk_in_optimized_order(const Field_VectorPreservingVisitor<T>& visitor)
	{
		return walk_in_default_order(visitor);
	}

	/**
	 * Visits (but does not alter) some entries of this vector in optimized
	 * order. The order in which the entries are visited is selected so as to
	 * lead to the most efficient implementation; it might depend on the
	 * concrete implementation of this virtual class.
	 *
	 * @param visitor visitor to be used to process the entries of this vector
	 * @param start the index of the first entry to be visited
	 * @param end the index of the last entry to be visited (inclusive)
	 * @return the value returned by {@link Field_VectorPreservingVisitor#end()}
	 * at the end of the walk
	 * @ if {@code end < start}.
	 * @ if the indices are not valid.
	 */
	T walk_in_optimized_order(const Field_VectorPreservingVisitor<T>& visitor, const int& start, const int& end)
	{
		return walk_in_default_order(visitor, start, end);
	}

	/**
	 * Visits (and possibly alters) all entries of this vector in default order
	 * (increasing index).
	 *
	 * @param visitor the visitor to be used to process and modify the entries
	 * of this vector
	 * @return the value returned by {@link Field_VectorChangingVisitor#end()}
	 * at the end of the walk
	 */
	T walk_in_default_order(const Field_VectorChangingVisitor<T>& visitor)
	{
		const int dim = get_dimension();
		visitor.start(dim, 0, dim - 1);
		for (int i{}; i < dim; i++)
		{
			set_entry(i, visitor.visit(i, get_entry(i)));
		}
		return visitor.end();
	}

	/**
	 * Visits (and possibly alters) some entries of this vector in default order
	 * (increasing index).
	 *
	 * @param visitor visitor to be used to process the entries of this vector
	 * @param start the index of the first entry to be visited
	 * @param end the index of the last entry to be visited (inclusive)
	 * @return the value returned by {@link Field_VectorChangingVisitor#end()}
	 * at the end of the walk
	 * @ if {@code end < start}.
	 * @ if the indices are not valid.
	 */
	T walk_in_default_order(const Field_VectorChangingVisitor<T>& visitor, const int& start, const int& end)
	{
		check_indices(start, end);
		visitor.start(get_dimension(), start, end);
		for (int i = start; i <= end; i++)
		{
			set_entry(i, visitor.visit(i, get_entry(i)));
		}
		return visitor.end();
	}

	/**
	 * Visits (and possibly alters) all entries of this vector in optimized
	 * order. The order in which the entries are visited is selected so as to
	 * lead to the most efficient implementation; it might depend on the
	 * concrete implementation of this virtual class.
	 *
	 * @param visitor the visitor to be used to process the entries of this
	 * vector
	 * @return the value returned by {@link Field_VectorChangingVisitor#end()}
	 * at the end of the walk
	 */
	T walk_in_optimized_order(const Field_VectorChangingVisitor<T>& visitor)
	{
		return walk_in_default_order(visitor);
	}

	/**
	 * Visits (and possibly change) some entries of this vector in optimized
	 * order. The order in which the entries are visited is selected so as to
	 * lead to the most efficient implementation; it might depend on the
	 * concrete implementation of this virtual class.
	 *
	 * @param visitor visitor to be used to process the entries of this vector
	 * @param start the index of the first entry to be visited
	 * @param end the index of the last entry to be visited (inclusive)
	 * @return the value returned by {@link Field_VectorChangingVisitor#end()}
	 * at the end of the walk
	 * @ if {@code end < start}.
	 * @ if the indices are not valid.
	 */
	T walk_in_optimized_order(const Field_VectorChangingVisitor<T>& visitor, const int& start, const int& end)
	{
		return walk_in_default_order(visitor, start, end);
	}
	/** {@inherit_doc}
	 * @since 2.0
	 */
	 //override
	std::string to_string() const
	{
		auto builder = std::stringstream();
		builder.append('{');
		for (int i{}; i < data.size(); ++i)
		{
			if (i > 0)
			{
				builder.append("; ");
			}
			builder.append(data[i].to_string());
		}
		builder.append('}');
		return builder.to_string();
	}

	/**
	 * Test for the equality of two vectors.
	 *
	 * @param other Object to test for equality.
	 * @return {@code true} if two vector objects are equal, {@code false}
	 * otherwise.
	 */
	 //override
	bool equals(cosnt Object& other)
	{
		if (*this == other)
		{
			return true;
		}
		if (other == NULL)
		{
			return false;
		}

		try
		{
			////@Suppress_Warnings("unchecked") // May fail, but we ignore Class_Cast_Exception
			Field_Vector<T> rhs = (Field_Vector<T>) other;
			if (data.size() != rhs.get_dimension())
			{
				return false;
			}

			for (int i{}; i < data.size(); ++i)
			{
				if (!data[i].equals(rhs.get_entry(i)))
				{
					return false;
				}
			}
			return true;
		}
		catch (Class_Cast_Exception ex)
		{
			// ignore exception
			return false;
		}
	}

	/**
	 * Get a hash_code for the real vector.
	 * <p>All NaN values have the same hash code.</p>
	 * @return a hash code value for this object
	 */
	 //override
	int hash_code()
	{
		int h{ 3542 };
		for (const T a : my_data)
		{
			h ^= a.hash_code();
		}
		return h;
	}
};