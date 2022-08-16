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

  /**
   * This class : the {@link Field_Vector} interface with a {@link Field_Element} array.
   * @param <T> the type of the field elements
   */
class ArrayField_Vector<T extends Field_Element<T>> : Field_Vector<T>
{
	7648186910365927050L;

	/** Entries of the vector. */
	private std::vector<T> data;

	/** Field to which the elements belong. */
	private const Field<T> field;

	/**
	 * Build a 0-length vector.
	 * Zero-length vectors may be used to initialize construction of vectors
	 * by data gathering. We start with zero-length and use either the {@link
	 * #ArrayField_Vector(Field_Vector, Field_Vector)} constructor
	 * or one of the {@code append} methods ({@link #add(Field_Vector)} or
	 * {@link #append(ArrayField_Vector)}) to gather data into this vector.
	 *
	 * @param field field to which the elements belong
	 */
	public ArrayField_Vector(const Field<T> field)
	{
		this(field, 0);
	}

	/**
	 * Construct a vector of zeroes.
	 *
	 * @param field Field to which the elements belong.
	 * @param size Size of the vector.
	 */
	public ArrayField_Vector(Field<T> field, int size)
	{
		this.field = field;
		this.data = Math_Arrays::build_array(field, size);
	}

	/**
	 * Construct a vector with preset values.
	 *
	 * @param size Size of the vector.
	 * @param preset All entries will be set with this value.
	 */
	public ArrayField_Vector(const int& size, T preset)
	{
		this(preset.get_field(), size);
		Arrays.fill(data, preset);
	}

	/**
	 * Construct a vector from an array, copying the input array.
	 * This constructor needs a non-empty {@code d} array to retrieve
	 * the field from its first element. This implies it cannot build
	 * 0 length vectors. To build vectors from any size, one should
	 * use the {@link #ArrayField_Vector(Field, Field_Element[])} constructor.
	 *
	 * @param d Array.
	 * @ if {@code d} is {@code NULL}.
	 * @ if {@code d} is empty.
	 * @see #ArrayField_Vector(Field, Field_Element[])
	 */
	public ArrayField_Vector(std::vector<T> d)

	{
		//Math_Utils::check_not_null(d);
		try
		{
			field = d[0].get_field();
			data = d.clone();
		}
		catch (Array_indexOutOfboundsException e)
		{
			throw (e, hipparchus::exception::Localized_Core_Formats_Type::VECTOR_MUST_HAVE_AT_LEAST_ONE_ELEMENT);
		}
	}

	/**
	 * Construct a vector from an array, copying the input array.
	 *
	 * @param field Field to which the elements belong.
	 * @param d Array.
	 * @ if {@code d} is {@code NULL}.
	 * @see #ArrayField_Vector(Field_Element[])
	 */
	public ArrayField_Vector(Field<T> field, std::vector<T> d)

	{
		//Math_Utils::check_not_null(d);
		this.field = field;
		data = d.clone();
	}

	/**
	 * Create a ArrayField_Vector using the input array as the underlying
	 * data array.
	 * If an array is built specially in order to be embedded in a
	 * ArrayField_Vector and not used directly, the {@code copy_array} may be
	 * set to {@code false}. This will prevent the copying and improve
	 * performance as no array will be built and no data will be copied.
	 * This constructor needs a non-empty {@code d} array to retrieve
	 * the field from its first element. This implies it cannot build
	 * 0 length vectors. To build vectors from any size, one should
	 * use the {@link #ArrayField_Vector(Field, Field_Element[], bool)}
	 * constructor.
	 *
	 * @param d Data for the vector.
	 * @param copy_array If {@code true}, the input array will be copied, * otherwise it will be referenced.
	 * @ if {@code d} is {@code NULL}.
	 * @ if {@code d} is empty.
	 * @see #ArrayField_Vector(Field_Element[])
	 * @see #ArrayField_Vector(Field, Field_Element[], bool)
	 */
	public ArrayField_Vector(std::vector<T> d, bool copy_array)

	{
		//Math_Utils::check_not_null(d);
		if (d.size() == 0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::VECTOR_MUST_HAVE_AT_LEAST_ONE_ELEMENT);
		}
		field = d[0].get_field();
		data = copy_array ? d.clone() : d;
	}

	/**
	 * Create a ArrayField_Vector using the input array as the underlying
	 * data array.
	 * If an array is built specially in order to be embedded in a
	 * ArrayField_Vector and not used directly, the {@code copy_array} may be
	 * set to {@code false}. This will prevent the copying and improve
	 * performance as no array will be built and no data will be copied.
	 *
	 * @param field Field to which the elements belong.
	 * @param d Data for the vector.
	 * @param copy_array If {@code true}, the input array will be copied, * otherwise it will be referenced.
	 * @ if {@code d} is {@code NULL}.
	 * @see #ArrayField_Vector(Field_Element[], bool)
	 */
	public ArrayField_Vector(Field<T> field, std::vector<T> d, bool copy_array)

	{
		//Math_Utils::check_not_null(d);
		this.field = field;
		data = copy_array ? d.clone() : d;
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
	public ArrayField_Vector(std::vector<T> d, int pos, int size)

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
	public ArrayField_Vector(Field<T> field, std::vector<T> d, int pos, int size)
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
	public ArrayField_Vector(Field_Vector<T> v)

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
	public ArrayField_Vector(ArrayField_Vector<T> v)

	{
		//Math_Utils::check_not_null(v);
		field = v.get_field();
		data = v.data.clone();
	}

	/**
	 * Construct a vector from another vector.
	 *
	 * @param v Vector to copy.
	 * @param deep If {@code true} perform a deep copy, otherwise perform
	 * a shallow copy
	 * @ if {@code v} is {@code NULL}.
	 */
	public ArrayField_Vector(ArrayField_Vector<T> v, bool deep)

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
	public ArrayField_Vector(Field_Vector<T> v1, Field_Vector<T> v2)

	{
		//Math_Utils::check_not_null(v1);
		//Math_Utils::check_not_null(v2);
		field = v1.get_field();
		const std::vector<T> v1_data =
			(v1 instanceof ArrayField_Vector) ? ((ArrayField_Vector<T>) v1).data : v1.to_array();
		const std::vector<T> v2_data =
			(v2 instanceof ArrayField_Vector) ? ((ArrayField_Vector<T>) v2).data : v2.to_array();
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
	public ArrayField_Vector(Field_Vector<T> v1, std::vector<T> v2)

	{
		//Math_Utils::check_not_null(v1);
		//Math_Utils::check_not_null(v2);
		field = v1.get_field();
		const std::vector<T> v1_data =
			(v1 instanceof ArrayField_Vector) ? ((ArrayField_Vector<T>) v1).data : v1.to_array();
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
	public ArrayField_Vector(std::vector<T> v1, Field_Vector<T> v2)

	{
		//Math_Utils::check_not_null(v1);
		//Math_Utils::check_not_null(v2);
		field = v2.get_field();
		const std::vector<T> v2_data =
			(v2 instanceof ArrayField_Vector) ? ((ArrayField_Vector<T>) v2).data : v2.to_array();
		data = Math_Arrays::build_array(field, v1.size() + v2_data.size());
		System.arraycopy(v1, 0, data, 0, v1.size());
		System.arraycopy(v2_data, 0, data, v1.size(), v2_data.size());
	}

	/**
	 * Construct a vector by appending one vector to another vector.
	 * This constructor needs at least one non-empty array to retrieve
	 * the field from its first element. This implies it cannot build
	 * 0 length vectors. To build vectors from any size, one should
	 * use the {@link #ArrayField_Vector(Field, Field_Element[], Field_Element[])}
	 * constructor.
	 *
	 * @param v1 First vector (will be put in front of the vector).
	 * @param v2 Second vector (will be put at back of the vector).
	 * @ if {@code v1} or {@code v2} is
	 * {@code NULL}.
	 * @ if both arrays are empty.
	 * @see #ArrayField_Vector(Field, Field_Element[], Field_Element[])
	 */
	public ArrayField_Vector(std::vector<T> v1, std::vector<T> v2)

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
	 * @see #ArrayField_Vector(Field_Element[], Field_Element[])
	 */
	public ArrayField_Vector(Field<T> field, std::vector<T> v1, std::vector<T> v2)

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
	public Field<T> get_field()
	{
		return field;
	}

	/** {@inherit_doc} */
	//override
	public Field_Vector<T> copy()
	{
		return ArrayField_Vector<T>(this, true);
	}

	/** {@inherit_doc} */
	//override
	public Field_Vector<T> add(const Field_Vector<T>& v)
	{
		if (dynamic_cast<const ArrayField_Vector*>(*v) != nullptr)
		{
			return add((ArrayField_Vector<T>) v);
		}
		check_vector_dimensions(v);
		auto out = Math_Arrays::build_array(field, data.size());
		for (int i{}; i < data.size(); i++)
		{
			out[i] = data[i].add(v.get_entry(i));
		}
		return ArrayField_Vector<T>(field, out, false);
	}

	/**
	 * Compute the sum of {@code this} and {@code v}.
	 * @param v vector to be added
	 * @return {@code this + v}
	 * @ if {@code v} is not the same size as
	 * {@code this}
	 */
	public ArrayField_Vector<T> add(ArrayField_Vector<T> v)

	{
		check_vector_dimensions(v.data.size());
		std::vector<T> out = Math_Arrays::build_array(field, data.size());
		for (int i{}; i < data.size(); i++)
		{
			out[i] = data[i].add(v.data[i]);
		}
		return ArrayField_Vector<T>(field, out, false);
	}

	/** {@inherit_doc} */
	//override
	public Field_Vector<T> subtract(Field_Vector<T> v)
	{
		if (dynamic_cast<const ArrayField_Vector*>(*v) != nullptr)
		{
			return subtract((ArrayField_Vector<T>) v);
		}
		check_vector_dimensions(v);
		std::vector<T> out = Math_Arrays::build_array(field, data.size());
		for (int i{}; i < data.size(); i++)
		{
			out[i] = data[i].subtract(v.get_entry(i));
		}
		return ArrayField_Vector<T>(field, out, false);
	}

	/**
	 * Compute {@code this} minus {@code v}.
	 * @param v vector to be subtracted
	 * @return {@code this - v}
	 * @ if {@code v} is not the same size as
	 * {@code this}
	 */
	public ArrayField_Vector<T> subtract(ArrayField_Vector<T> v)

	{
		check_vector_dimensions(v.data.size());
		std::vector<T> out = Math_Arrays::build_array(field, data.size());
		for (int i{}; i < data.size(); i++)
		{
			out[i] = data[i].subtract(v.data[i]);
		}
		return ArrayField_Vector<T>(field, out, false);
	}

	/** {@inherit_doc} */
	//override
	public Field_Vector<T> map_add(T d)
	{
		std::vector<T> out = Math_Arrays::build_array(field, data.size());
		for (int i{}; i < data.size(); i++)
		{
			out[i] = data[i].add(d);
		}
		return ArrayField_Vector<T>(field, out, false);
	}

	/** {@inherit_doc} */
	//override
	public Field_Vector<T> map_add_to_self(T d)
	{
		for (int i{}; i < data.size(); i++)
		{
			data[i] = data[i].add(d);
		}
		return this;
	}

	/** {@inherit_doc} */
	//override
	public Field_Vector<T> map_subtract(T d)
	{
		std::vector<T> out = Math_Arrays::build_array(field, data.size());
		for (int i{}; i < data.size(); i++)
		{
			out[i] = data[i].subtract(d);
		}
		return ArrayField_Vector<T>(field, out, false);
	}

	/** {@inherit_doc} */
	//override
	public Field_Vector<T> map_subtract_to_self(T d)
	{
		for (int i{}; i < data.size(); i++)
		{
			data[i] = data[i].subtract(d);
		}
		return this;
	}

	/** {@inherit_doc} */
	//override
	public Field_Vector<T> map_multiply(T d)
	{
		std::vector<T> out = Math_Arrays::build_array(field, data.size());
		for (int i{}; i < data.size(); i++)
		{
			out[i] = data[i].multiply(d);
		}
		return ArrayField_Vector<T>(field, out, false);
	}

	/** {@inherit_doc} */
	//override
	public Field_Vector<T> map_multiply_to_self(T d)
	{
		for (int i{}; i < data.size(); i++)
		{
			data[i] = data[i].multiply(d);
		}
		return this;
	}

	/** {@inherit_doc} */
	//override
	public Field_Vector<T> map_divide(T d)
		, Math_Runtime_Exception
	{
		//Math_Utils::check_not_null(d);
		std::vector<T> out = Math_Arrays::build_array(field, data.size());
		for (int i{}; i < data.size(); i++)
		{
			out[i] = data[i].divide(d);
		}
		return ArrayField_Vector<T>(field, out, false);
	}

		/** {@inherit_doc} */
		//override
		public Field_Vector<T> map_divide_to_self(T d)
		, Math_Runtime_Exception
	{
		//Math_Utils::check_not_null(d);
		for (int i{}; i < data.size(); i++)
		{
			data[i] = data[i].divide(d);
		}
		return this;
	}

		/** {@inherit_doc} */
		//override
		public Field_Vector<T> map_inv() Math_Runtime_Exception
	{
		std::vector<T> out = Math_Arrays::build_array(field, data.size());
		const T one = field.get_one();
		for (int i{}; i < data.size(); i++)
		{
			try
			{
				out[i] = one.divide(data[i]);
			}
			catch (const Math_Runtime_Exception e)
			{
				throw Math_Runtime_Exception(e, hipparchus::exception::Localized_Core_Formats_Type::INDEX, i);
			}
		}
		return ArrayField_Vector<T>(field, out, false);
	}

	/** {@inherit_doc} */
	//override
	public Field_Vector<T> map_inv_to_self() Math_Runtime_Exception
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
				throw Math_Runtime_Exception(e, hipparchus::exception::Localized_Core_Formats_Type::INDEX, i);
			}
		}
		return this;
	}

	/** {@inherit_doc} */
	//override
	public Field_Vector<T> ebe_multiply(Field_Vector<T> v)
	{
		if (dynamic_cast<const ArrayField_Vector*>(*v) != nullptr)
		{
			return ebe_multiply((ArrayField_Vector<T>) v);
		}
		check_vector_dimensions(v);
		std::vector<T> out = Math_Arrays::build_array(field, data.size());
		for (int i{}; i < data.size(); i++)
		{
			out[i] = data[i].multiply(v.get_entry(i));
		}
		return ArrayField_Vector<T>(field, out, false);
	}

	/**
	 * Element-by-element multiplication.
	 * @param v vector by which instance elements must be multiplied
	 * @return a vector containing {@code this[i] * v[i]} for all {@code i}
	 * @ if {@code v} is not the same size as
	 * {@code this}
	 */
	public ArrayField_Vector<T> ebe_multiply(ArrayField_Vector<T> v)

	{
		check_vector_dimensions(v.data.size());
		std::vector<T> out = Math_Arrays::build_array(field, data.size());
		for (int i{}; i < data.size(); i++)
		{
			out[i] = data[i].multiply(v.data[i]);
		}
		return ArrayField_Vector<T>(field, out, false);
	}

	/** {@inherit_doc} */
	//override
	public Field_Vector<T> ebe_divide(Field_Vector<T> v)
	{
		if (dynamic_cast<const ArrayField_Vector*>(*v) != nullptr)
		{
			return ebe_divide((ArrayField_Vector<T>) v);
		}
		check_vector_dimensions(v);
		std::vector<T> out = Math_Arrays::build_array(field, data.size());
		for (int i{}; i < data.size(); i++)
		{
			try
			{
				out[i] = data[i].divide(v.get_entry(i));
			}
			catch (const Math_Runtime_Exception e)
			{
				throw Math_Runtime_Exception(e, hipparchus::exception::Localized_Core_Formats_Type::INDEX, i);
			}
		}
		return ArrayField_Vector<T>(field, out, false);
	}

	/**
	 * Element-by-element division.
	 * @param v vector by which instance elements must be divided
	 * @return a vector containing {@code this[i] / v[i]} for all {@code i}
	 * @ if {@code v} is not the same size as
	 * {@code this}
	 * @Math_Runtime_Exception if one entry of {@code v} is zero.
	 */
	public ArrayField_Vector<T> ebe_divide(ArrayField_Vector<T> v)
		, Math_Runtime_Exception
	{
	check_vector_dimensions(v.data.size());
	std::vector<T> out = Math_Arrays::build_array(field, data.size());
	for (int i{}; i < data.size(); i++)
	{
		try
		{
			out[i] = data[i].divide(v.data[i]);
		}
catch (const Math_Runtime_Exception e)
			{
				throw Math_Runtime_Exception(e, hipparchus::exception::Localized_Core_Formats_Type::INDEX, i);
			}
		}
		return ArrayField_Vector<T>(field, out, false);
	}

		/**
		 * Returns a reference to the underlying data array.
		 * <p>Does not make a fresh copy of the underlying data.</p>
		 * @return array of entries
		 */
		public std::vector<T> get_data_ref()
	{
		return data; // NOPMD - returning an internal array is intentional and documented here
	}

	/** {@inherit_doc} */
	//override
	public T dot_product(Field_Vector<T> v)
	{
		if (dynamic_cast<const ArrayField_Vector*>(*v) != nullptr)
		{
			return dot_product((ArrayField_Vector<T>) v);
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
	public T dot_product(ArrayField_Vector<T> v)

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
	public Field_Vector<T> projection(Field_Vector<T> v)
		, Math_Runtime_Exception
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
		public ArrayField_Vector<T> projection(ArrayField_Vector<T> v)
		, Math_Runtime_Exception
	{
	return (ArrayField_Vector<T>) v.map_multiply(dot_product(v).divide(v.dot_product(v)));
	}

		/** {@inherit_doc} */
		//override
		public Field_Matrix<T> outer_product(Field_Vector<T> v)
	{
		if (dynamic_cast<const ArrayField_Vector*>(*v) != nullptr)
		{
			return outer_product((ArrayField_Vector<T>) v);
		}
		const int m = data.size();
		const int n = v.get_dimension();
		const Field_Matrix<T> out = Array2DRowField_Matrix<>(field, m, n);
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
	public Field_Matrix<T> outer_product(ArrayField_Vector<T> v)
	{
		const int m = data.size();
		const int n = v.data.size();
		const Field_Matrix<T> out = Array2DRowField_Matrix<>(field, m, n);
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
	public T get_entry(const int& index)
	{
		return data[index];
	}

	/** {@inherit_doc} */
	//override
	public int get_dimension()
	{
		return data.size();
	}

	/** {@inherit_doc} */
	//override
	public Field_Vector<T> append(const Field_Vector<T>& v)
	{
		if (dynamic_cast<const ArrayField_Vector*>(*v) != nullptr)
		{
			return append((ArrayField_Vector<T>) v);
		}
		return ArrayField_Vector<T>(this, ArrayField_Vector<T>(v));
	}

	/**
	 * Construct a vector by appending a vector to this vector.
	 * @param v vector to append to this one.
	 * @return a vector
	 */
	public ArrayField_Vector<T> append(const ArrayField_Vector<T>& v)
	{
		return ArrayField_Vector<T>(*this, v);
	}

	/** {@inherit_doc} */
	//override
	public Field_Vector<T> append(T in)
	{
		const std::vector<T> out = Math_Arrays::build_array(field, data.size() + 1);
		System.arraycopy(data, 0, out, 0, data.size());
		out[data.size()] = in;
		return ArrayField_Vector<T>(field, out, false);
	}

	/** {@inherit_doc} */
	//override
	public Field_Vector<T> get_sub_vector(const int& index, int n)
	{
		if (n < 0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_OF_ELEMENTS_SHOULD_BE_POSITIVE, n);
		}
		ArrayField_Vector<T> out = ArrayField_Vector<>(field, n);
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
	public void set_entry(const int& index, T value)
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
	public void set_sub_vector(const int& index, Field_Vector<T> v)
	{
		try
		{
			if (dynamic_cast<const ArrayField_Vector*>(*v) != nullptr)
			{
				set(index, (ArrayField_Vector<T>) v);
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
	public void set(const int& index, ArrayField_Vector<T> v)
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
	public void set(T value)
	{
		Arrays.fill(data, value);
	}

	/** {@inherit_doc} */
	//override
	public std::vector<T> to_array()
	{
		return data.clone();
	}

	/**
	 * Check if instance and specified vectors have the same dimension.
	 * @param v vector to compare instance with
	 * @exception  if the vectors do not
	 * have the same dimensions
	 */
	protected void check_vector_dimensions(Field_Vector<T> v)

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
	protected void check_vector_dimensions(const int& n)
	{
		if (data.size() != n)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, data.size(), n);
		}
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
	public T walk_in_default_order(const Field_VectorPreservingVisitor<T> visitor)
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
	public T walk_in_default_order(const Field_VectorPreservingVisitor<T> visitor, const int start, const int end)

	{
		check_indices(start, end);
		visitor.start(get_dimension(), start, end);
		for (int i = start; i <= end; i++)
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
	public T walk_in_optimized_order(const Field_VectorPreservingVisitor<T> visitor)
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
	public T walk_in_optimized_order(const Field_VectorPreservingVisitor<T> visitor, const int start, const int end)

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
	public T walk_in_default_order(const Field_VectorChangingVisitor<T> visitor)
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
	public T walk_in_default_order(const Field_VectorChangingVisitor<T> visitor, const int start, const int end)

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
	public T walk_in_optimized_order(const Field_VectorChangingVisitor<T> visitor)
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
	public T walk_in_optimized_order(const Field_VectorChangingVisitor<T> visitor, const int start, const int end)

	{
		return walk_in_default_order(visitor, start, end);
	}
	/** {@inherit_doc}
	 * @since 2.0
	 */
	 //override
	public std::string to_string() const
	{
		const std::stringBuilder builder = std::stringstream();
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
	public bool equals(Object other)
	{
		if (this == other)
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
	public int hash_code()
	{
		int h = 3542;
		for (const T a : data)
		{
			h ^= a.hash_code();
		}
		return h;
	}

	/**
	 * Check if an index is valid.
	 *
	 * @param index Index to check.
	 * @exception  if the index is not valid.
	 */
	private void check_index(const int index)
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
	private void check_indices(const int start, const int end)
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
}
