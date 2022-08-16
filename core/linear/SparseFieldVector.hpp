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

  //import org.hipparchus.Field;
  //import org.hipparchus.Field_Element;
  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.exception.Math_Runtime_Exception;
  //import org.hipparchus.exception.;
  //import org.hipparchus.util.Math_Arrays;
  //import org.hipparchus.util.Math_Utils;
  //import org.hipparchus.util.Open_Int_To_Field_Hash_Map;
#include <type_traits>
#include <vector>
#include "../FieldElement.h"

/**
 * This class : the {@link Field_Vector} interface with a {@link Open_Int_To_Field_Hash_Map} backing store.
 * <p>
 *  Caveat: This implementation assumes that, for any {@code x}, *  the equality {@code x * 0 == 0} holds. But it is is not true for
 *  {@code NaN}. Moreover, zero entries will lose their sign.
 *  Some operations (that involve {@code NaN} and/or infinities) may
 *  thus give incorrect results.
 * </p>
 * @param <T> the type of the field elements
 */
template<typename T, typename std::enable_if<std::is_base_of<Field_Element<T>, T>::value>::type* = nullptr>
class Sparse_Field_Vector : Field_Vector<T>
{
private:
	/** Field to which the elements belong. */
	const Field<T> my_field;
	/** Entries of the vector. */
	const Open_Int_To_Field_Hash_Map<T> my_entries;
	/** Dimension of the vector. */
	const int my_virtual_size;

	/**
	 * Get the entries of this instance.
	 *
	 * @return the entries of this instance
	 */
	Open_Int_To_Field_Hash_Map<T> get_entries() const
	{
		return my_entries;
	}

	/**
	 * Check whether an index is valid.
	 *
	 * @param index Index to check.
	 * @ if the index is not valid.
	 */
	void check_index(const int& index)
	{
		Math_Utils::check_range_inclusive(index, 0, get_dimension() - 1);
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
		throw std::exception("not implemented");
		//if ((start < 0) || (start >= dim))
		//{
		//    throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::INDEX, start, 0, dim - 1);
		//}
		//if ((end < 0) || (end >= dim))
		//{
		//    throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::INDEX, end, 0, dim - 1);
		//}
		//if (end < start)
		//{
		//    throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::INITIAL_ROW_AFTER_FINAL_ROW, end, start, false);
		//}
	}

public:
	/**
	 * Build a 0-length vector.
	 * Zero-length vectors may be used to initialize construction of vectors
	 * by data gathering. We start with zero-length and use either the {@link
	 * #Sparse_Field_Vector(Sparse_Field_Vector, int)} constructor
	 * or one of the {@code append} method ({@link #append(Field_Vector)} or
	 * {@link #append(Sparse_Field_Vector)}) to gather data into this vector.
	 *
	 * @param field Field to which the elements belong.
	 */
	Sparse_Field_Vector(const Field<T>& field)
	{
		Sparse_Field_Vector(field, 0);
	}

	/**
	 * Construct a vector of zeroes.
	 *
	 * @param field Field to which the elements belong.
	 * @param dimension Size of the vector.
	 */
	Sparse_Field_Vector(const Field<T>& field, const int& dimension)
		:
		my_field{ field },
		my_virtual_size{ dimension },
		my_entries{ Open_Int_To_Field_Hash_Map<>(field) }
	{
	}

protected:
	/**
	 * Build a resized vector, for use with append.
	 *
	 * @param v Original vector
	 * @param resize Amount to add.
	 */
	Sparse_Field_Vector(const Sparse_Field_Vector<T>& v, const int& resize)
		:
		my_field{ v.field },
		my_virtual_size{ v.get_dimension() + resize },
		my_entries{ Open_Int_To_Field_Hash_Map<>(v.entries) }
	{
	}

	/**
	 * Check if instance dimension is equal to some expected value.
	 *
	 * @param n Expected dimension.
	 * @ if the dimensions do not match.
	 */
	void check_vector_dimensions(const int& n)
	{
		if (get_dimension() != n)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, get_dimension(), n);
		}
	}

public:
	/**
	 * Build a vector with known the sparseness (for advanced use only).
	 *
	 * @param field Field to which the elements belong.
	 * @param dimension Size of the vector.
	 * @param expected_size Expected number of non-zero entries.
	 */
	Sparse_Field_Vector(const Field<T>& field, const int& dimension, const int& expected_size)
		:
		my_field{ field },
		my_virtual_size{ dimension },
		my_entries{ Open_Int_To_Field_Hash_Map<>(field,expected_size) }
	{
	}

	/**
	 * Create from a Field array.
	 * Only non-zero entries will be stored.
	 *
	 * @param field Field to which the elements belong.
	 * @param values Set of values to create from.
	 * @exceptionif values is NULL
	 */
	Sparse_Field_Vector(const Field<T>& field, const std::vector<T>& values)
		:
		my_field{ field },
		my_virtual_size{ values.size() },
		my_entries{ Open_Int_To_Field_Hash_Map<>(field) }
	{
		//Math_Utils::check_not_null(values);
		for (int key{}; key < values.size(); key++)
		{
			T value = values[key];
			entries.put(key, value);
		}
	}

	/**
	 * Copy constructor.
	 *
	 * @param v Instance to copy.
	 */
	Sparse_Field_Vector(const Sparse_Field_Vector<T>& v)
		:
		my_field{ v.my_field },
		my_virtual_size{ v.get_dimension() },
		my_entries{ Open_Int_To_Field_Hash_Map<>(v.get_entries()) }
	{
	}

	/**
	 * Optimized method to add sparse vectors.
	 *
	 * @param v Vector to add.
	 * @return {@code this + v}.
	 * @ if {@code v} is not the same size as
	 * {@code this}.
	 */
	Field_Vector<T> add(const Sparse_Field_Vector<T>& v)
	{
		check_vector_dimensions(v.get_bimension());
		Sparse_Field_Vector<T> res = (Sparse_Field_Vector<T>)copy();
		Open_Int_To_Field_Hash_Map<T>.Iterator iter = v.get_entries().iterator();
		while (iter.has_next())
		{
			iter.advance();
			int key = iter.key();
			T value = iter.value();
			if (entries.contains_key(key))
			{
				res.set_entry(key, entries.get(key).add(value));
			}
			else
			{
				res.set_entry(key, value);
			}
		}
		return res;
	}

	/**
	 * Construct a vector by appending a vector to this vector.
	 *
	 * @param v Vector to append to this one.
	 * @return a vector.
	 */
	Field_Vector<T> append(Sparse_Field_Vector<T> v)
	{
		Sparse_Field_Vector<T> res = Sparse_Field_Vector<>(this, v.get_dimension());
		Open_Int_To_Field_Hash_Map<T>.Iterator iter = v.entries.iterator();
		while (iter.has_next())
		{
			iter.advance();
			res.set_entry(iter.key() + virtual_size, iter.value());
		}
		return res;
	}

	/** {@inherit_doc} */
	//override
	Field_Vector<T> append(Field_Vector<T> v)
	{
		if (dynamic_cast<const Sparse_Field_Vector*>(*v) != nullptr)
		{
			return append((Sparse_Field_Vector<T>) v);
		}
		const int n = v.get_dimension();
		Field_Vector<T> res = Sparse_Field_Vector<>(this, n);
		for (int i{}; i < n; i++)
		{
			res.set_entry(i + virtual_size, v.get_entry(i));
		}
		return res;
	}

	/** {@inherit_doc}
	 * @exceptionif d is NULL
	 */
	 //override
	Field_Vector<T> append(const T& d)
	{
		//Math_Utils::check_not_null(d);
		Field_Vector<T> res = Sparse_Field_Vector<>(*this, 1);
		res.set_entry(virtual_size, d);
		return res;
	}

	/** {@inherit_doc} */
	//override
	Field_Vector<T> copy()
	{
		return Sparse_Field_Vector<T>(*this);
	}

	/** {@inherit_doc} */
	//override
	T dot_product(const Field_Vector<T>& v)
	{
		check_vector_dimensions(v.get_bimension());
		T res = field.get_zero();
		Open_Int_To_Field_Hash_Map<T>.Iterator iter = entries.iterator();
		while (iter.has_next())
		{
			iter.advance();
			res = res.add(v.get_entry(iter.key()).multiply(iter.value()));
		}
		return res;
	}

	/** {@inherit_doc} */
	//override
	Field_Vector<T> ebe_divide(const Field_Vector<T>& v)
	{
		check_vector_dimensions(v.get_bimension());
		Sparse_Field_Vector<T> res = Sparse_Field_Vector<>(this);
		Open_Int_To_Field_Hash_Map<T>.Iterator iter = res.entries.iterator();
		while (iter.has_next())
		{
			iter.advance();
			res.set_entry(iter.key(), iter.value().divide(v.get_entry(iter.key())));
		}
		return res;
	}

	/** {@inherit_doc} */
	//override
	Field_Vector<T> ebe_multiply(const Field_Vector<T>& v)
	{
		check_vector_dimensions(v.get_bimension());
		Sparse_Field_Vector<T> res = Sparse_Field_Vector<>(this);
		Open_Int_To_Field_Hash_Map<T>.Iterator iter = res.entries.iterator();
		while (iter.has_next())
		{
			iter.advance();
			res.set_entry(iter.key(), iter.value().multiply(v.get_entry(iter.key())));
		}
		return res;
	}

	/** {@inherit_doc} */
	//override
	int get_dimension() const
	{
		return my_virtual_size;
	}

	/** {@inherit_doc} */
	//override
	T get_entry(const int& index)
	{
		check_index(index);
		return entries.get(index);
	}

	/** {@inherit_doc} */
	//override
	Field<T> get_field() const
	{
		return my_field;
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
		check_index(index);
		check_index(index + n - 1);
		Sparse_Field_Vector<T> res = Sparse_Field_Vector<>(field, n);
		int end{ index + n };
		Open_Int_To_Field_Hash_Map<T>.Iterator iter = entries.iterator();
		while (iter.has_next())
		{
			iter.advance();
			int key{ iter.key() };
			if (key >= index && key < end)
			{
				res.set_entry(key - index, iter.value());
			}
		}
		return res;
	}

	/** {@inherit_doc} */
	//override
	Field_Vector<T> map_add(const T& d)
	{
		return copy().map_add_to_self(d);
	}

	/** {@inherit_doc} */
	//override
	Field_Vector<T> map_add_to_self(const T& d)
	{
		for (int i{}; i < virtual_size; i++)
		{
			set_entry(i, get_entry(i).add(d));
		}
		return *this;
	}

	/** {@inherit_doc} */
	//override
	Field_Vector<T> map_divide(const T& d)
	{
		return copy().map_divide_to_self(d);
	}

	/** {@inherit_doc} */
	//override
	Field_Vector<T> map_divide_to_self(cosnt T& d)
	{
		Open_Int_To_Field_Hash_Map<T>.Iterator iter = entries.iterator();
		while (iter.has_next())
		{
			iter.advance();
			entries.put(iter.key(), iter.value().divide(d));
		}
		return *this;
	}

	/** {@inherit_doc} */
	//override
	Field_Vector<T> map_inv()
	{
		return copy().map_inv_to_self();
	}

	/** {@inherit_doc} */
	//override
	Field_Vector<T> map_inv_to_self()
	{
		for (int i{}; i < virtual_size; i++)
		{
			set_entry(i, field.get_one().divide(get_entry(i)));
		}
		return *this;
	}

	/** {@inherit_doc} */
	//override
	Field_Vector<T> map_multiply(const T& d)
	{
		return copy().map_multiply_to_self(d);
	}

	/** {@inherit_doc} */
	//override
	Field_Vector<T> map_multiply_to_self(const T& d)
	{
		Open_Int_To_Field_Hash_Map<T>.Iterator iter = entries.iterator();
		while (iter.has_next())
		{
			iter.advance();
			entries.put(iter.key(), iter.value().multiply(d));
		}
		return *this;
	}

	/** {@inherit_doc} */
	//override
	Field_Vector<T> map_subtract(T d)
	{
		return copy().map_subtract_to_self(d);
	}

	/** {@inherit_doc} */
	//override
	Field_Vector<T> map_subtract_to_self(const T& d)
	{
		return map_add_to_self(field.get_zero().subtract(d));
	}

	/**
	 * Optimized method to compute outer product when both vectors are sparse.
	 * @param v vector with which outer product should be computed
	 * @return the matrix outer product between instance and v
	 */
	Field_Matrix<T> outer_product(const Sparse_Field_Vector<T>& v)
	{
		const int n{ v.get_dimension() };
		SparseField_Matrix<T> res = SparseField_Matrix<>(field, virtual_size, n);
		Open_Int_To_Field_Hash_Map<T>.Iterator iter = entries.iterator();
		while (iter.has_next())
		{
			iter.advance();
			Open_Int_To_Field_Hash_Map<T>.Iterator iter2 = v.entries.iterator();
			while (iter2.has_next())
			{
				iter2.advance();
				res.set_entry(iter.key(), iter2.key(), iter.value().multiply(iter2.value()));
			}
		}
		return res;
	}

	/** {@inherit_doc} */
	//override
	Field_Matrix<T> outer_product(const Field_Vector<T>& v)
	{
		if (dynamic_cast<const Sparse_Field_Vector<T>*>(*v) != nullptr)
		{
			return outer_product((Sparse_Field_Vector<T>)v);
		}
		const int n = v.get_dimension();
		Field_Matrix<T> res = SparseField_Matrix<>(field, virtual_size, n);
		Open_Int_To_Field_Hash_Map<T>.Iterator iter = entries.iterator();
		while (iter.has_next())
		{
			iter.advance();
			int row{ iter.key() };
			Field_Element<T> value = iter.value();
			for (int col{}; col < n; col++)
			{
				res.set_entry(row, col, value.multiply(v.get_entry(col)));
			}
		}
		return res;
	}

	/** {@inherit_doc} */
	//override
	Field_Vector<T> projection(const Field_Vector<T>& v)
	{
		check_vector_dimensions(v.get_bimension());
		return v.map_multiply(dot_product(v).divide(v.dot_product(v)));
	}

	/** {@inherit_doc}
	 * @exceptionif value is NULL
	 */
	 //override
	void set(const T& value)
	{
		//Math_Utils::check_not_null(value);
		for (int i{}; i < virtual_size; i++)
		{
			set_entry(i, value);
		}
	}

	/** {@inherit_doc}
	 * @exceptionif value is NULL
	 */
	 //override
	void set_entry(const int& index, const T& value),
	{
		//Math_Utils::check_not_null(value);
		check_index(index);
		entries.put(index, value);
	}

		/** {@inherit_doc} */
		//override
		void set_sub_vector(const int& index, const Field_Vector<T>& v)
	{
		check_index(index);
		check_index(index + v.get_dimension() - 1);
		const int n{ v.get_dimension() };
		for (int i{}; i < n; i++)
		{
			set_entry(i + index, v.get_entry(i));
		}
	}

	/**
	 * Optimized method to compute {@code this} minus {@code v}.
	 * @param v vector to be subtracted
	 * @return {@code this - v}
	 * @ if {@code v} is not the same size as
	 * {@code this}.
	 */
	Sparse_Field_Vector<T> subtract(const Sparse_Field_Vector<T>& v)
	{
		check_vector_dimensions(v.get_bimension());
		Sparse_Field_Vector<T> res = (Sparse_Field_Vector<T>)copy();
		Open_Int_To_Field_Hash_Map<T>.Iterator iter = v.get_entries().iterator();
		while (iter.has_next())
		{
			iter.advance();
			int key = iter.key();
			if (entries.contains_key(key))
			{
				res.set_entry(key, entries.get(key).subtract(iter.value()));
			}
			else
			{
				res.set_entry(key, field.get_zero().subtract(iter.value()));
			}
		}
		return res;
	}

	/** {@inherit_doc} */
	//override
	Field_Vector<T> subtract(const Field_Vector<T>& v)
	{
		if (dynamic_cast<const Sparse_Field_Vector<T>*>(*v) != nullptr)
		{
			return subtract((Sparse_Field_Vector<T>)v);
		}
		const int n{ v.get_dimension() };
		check_vector_dimensions(n);
		Sparse_Field_Vector<T> res = Sparse_Field_Vector<>(*this);
		for (int i{}; i < n; i++)
		{
			if (entries.contains_key(i))
			{
				res.set_entry(i, entries.get(i).subtract(v.get_entry(i)));
			}
			else
			{
				res.set_entry(i, field.get_zero().subtract(v.get_entry(i)));
			}
		}
		return res;
	}

	/** {@inherit_doc} */
	//override
	std::vector<T> to_array()
	{
		std::vector<T> res = Math_Arrays::build_array(field, virtual_size);
		Open_Int_To_Field_Hash_Map<T>.Iterator iter = entries.iterator();
		while (iter.has_next())
		{
			iter.advance();
			res[iter.key()] = iter.value();
		}
		return res;
	}

	/** {@inherit_doc} */
	//override
	Field_Vector<T> add(Field_Vector<T> v)
	{
		if (dynamic_cast<const Sparse_Field_Vector<T>*>(*v) != nullptr)
		{
			return add((Sparse_Field_Vector<T>) v);
		}
		const int n{ v.get_dimension() };
		check_vector_dimensions(n);
		Sparse_Field_Vector<T> res = Sparse_Field_Vector<>(field, get_dimension());
		for (int i{}; i < n; i++)
		{
			res.set_entry(i, v.get_entry(i).add(get_entry(i)));
		}
		return res;
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
	T walk_in_default_order(const Field_VectorPreservingVisitor<T> visitor)
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
		const int dim{ get_dimension() };
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
		for (int i{ start }; i <= end; i++)
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

	/** {@inherit_doc} */
	//override
	int hash_code()
	{
		constexpr int prime{ 31 };
		int result{ 1 };
		result = prime * result + ((field == NULL) ? 0 : field.hash_code());
		result = prime * result + virtual_size;
		Open_Int_To_Field_Hash_Map<T>.Iterator iter = entries.iterator();
		while (iter.has_next())
		{
			iter.advance();
			int temp = iter.value().hash_code();
			result = prime * result + temp;
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	bool equals(const Object& obj)
	{
		if (*this == obj)
		{
			return true;
		}

		if (!(obj instanceof Sparse_Field_Vector< ? >))
		{
			return false;
		}

		////@Suppress_Warnings("unchecked") // OK, because "else if" check below ensures that
									   // other must be the same type as this
		Sparse_Field_Vector<T> other = (Sparse_Field_Vector<T>) obj;
		if (field == NULL)
		{
			if (other.field != NULL)
			{
				return false;
			}
		}
		else if (!field.equals(other.field))
		{
			return false;
		}
		if (virtual_size != other.virtual_size)
		{
			return false;
		}

		Open_Int_To_Field_Hash_Map<T>.Iterator iter = entries.iterator();
		while (iter.has_next())
		{
			iter.advance();
			T test = other.get_entry(iter.key());
			if (!test.equals(iter.value()))
			{
				return false;
			}
		}
		iter = other.get_entries().iterator();
		while (iter.has_next())
		{
			iter.advance();
			T test = iter.value();
			if (!test.equals(get_entry(iter.key())))
			{
				return false;
			}
		}
		return true;
	}
};