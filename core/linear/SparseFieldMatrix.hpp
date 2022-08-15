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
#include <type_traits>
#include "MatrixUtils.h"
#include "../FieldElement.h"
#include "FieldMatrix.h"
//import org.hipparchus.Field;
//import org.hipparchus.Field_Element;
//import org.hipparchus.exception.;
//import org.hipparchus.util.Open_Int_To_Field_Hash_Map;

/**
 * Sparse matrix implementation based on an open addressed map.
 *
 * <p>
 *  Caveat: This implementation assumes that, for any {@code x}, *  the equality {@code x * 0 == 0} holds. But it is is not true for
 *  {@code NaN}. Moreover, zero entries will lose their sign.
 *  Some operations (that involve {@code NaN} and/or infinities) may
 *  thus give incorrect results.
 * </p>
 * @param <T> the type of the field elements
 */
template<typename T, typename std::enable_if<std::is_base_of<Field_Element<T>, T>::value>::type* = nullptr>
class SparseField_Matrix : public Abstract_Field_Matrix<T>
{
private:
	/** Storage for (sparse) matrix elements. */
	const Open_Int_To_Field_Hash_Map<T> my_entries;
	/** Row dimension. */
	const int my_rows;
	/** Column dimension. */
	const int my_columns;

	/**
	 * Compute the key to access a matrix element.
	 *
	 * @param row Row index of the matrix element.
	 * @param column Column index of the matrix element.
	 * @return the key within the map to access the matrix element.
	 */
	int compute_key(const int& row, const int& column)
	{
		return row * my_columns + column;
	}

public:
	/**
	 * Create a matrix with no data.
	 *
	 * @param field Field to which the elements belong.
	 */
	SparseField_Matrix(const Field<T>& field)
		:
		my_rows{ 0 },
		my_columns{ 0 },
		my_entries{ Open_Int_To_Field_Hash_Map<>(field) }
	{
		super(field);
	};

	/**
	 * Create a SparseField_Matrix<T> with the supplied row and column
	 * dimensions.
	 *
	 * @param field Field to which the elements belong.
	 * @param row_dimension Number of rows in the matrix.
	 * @param column_dimension Number of columns in the matrix.
	 * @org.hipparchus.exception.
	 * if row or column dimension is not positive.
	 */
	SparseField_Matrix(const Field<T>& field, const int& row_dimension, const int& column_dimension)
		:
		my_rows{ row_dimension },
		my_columns{ column_dimension },
		my_entries{ Open_Int_To_Field_Hash_Map<>(field) }
	{
		super(field, row_dimension, column_dimension);
	};

	/**
	 * Copy constructor.
	 *
	 * @param other Instance to copy.
	 */
	SparseField_Matrix(const SparseField_Matrix<T>& other)
		:
		my_rows{ other.get_row_dimension() },
		my_columns{ other.get_column_dimension() },
		my_entries{ Open_Int_To_Field_Hash_Map<>(other.entries) }
	{
		super(other.get_field(), other.get_row_dimension(), other.get_column_dimension());
	};

	/**
	 * Generic copy constructor.
	 *
	 * @param other Instance to copy.
	 */
	SparseField_Matrix(const Field_Matrix<T>& other)
	{
		super(other.get_field(), other.get_row_dimension(), other.get_column_dimension());
		rows = other.get_row_dimension();
		columns = other.get_column_dimension();
		entries = Open_Int_To_Field_Hash_Map<>(get_field());
		for (int i{}; i < rows; i++)
		{
			for (int j{}; j < columns; j++)
			{
				set_entry(i, j, other.get_entry(i, j));
			}
		}
	}

	/** {@inherit_doc} */
	//override
	void add_to_entry(const int& row, const int& column, const T& increment)
	{
		check_row_index(row);
		check_column_index(column);
		const int& key = compute_key(row, column);
		const T value = entries.get(key).add(increment);
		if (get_field().get_zero().equals(value))
		{
			entries.remove(key);
		}
		else
		{
			entries.put(key, value);
		}
	}

	/** {@inherit_doc} */
	//override
	Field_Matrix<T> copy()
	{
		return SparseField_Matrix<T>(*this);
	}

	/** {@inherit_doc} */
	//override
	Field_Matrix<T> create_matrix(const int& row_dimension, const int& column_dimension)
	{
		return SparseField_Matrix<T>(get_field(), row_dimension, column_dimension);
	}

	/** {@inherit_doc} */
	//override
	int get_column_dimension() const
	{
		return my_columns;
	}

	/** {@inherit_doc} */
	//override
	T get_entry(const int& row, const int& column)
	{
		check_row_index(row);
		check_column_index(column);
		return my_entries.get(compute_key(row, column));
	}

	/** {@inherit_doc} */
	//override
	int get_row_dimension() const
	{
		return my_rows;
	}

	/** {@inherit_doc} */
	//override
	void multiply_entry(const int& row, const int& column, const T& factor)
	{
		check_row_index(row);
		check_column_index(column);
		const int& key = compute_key(row, column);
		const T value = my_entries.get(key).multiply(factor);
		if (get_field().get_zero().equals(value))
		{
			my_entries.remove(key);
		}
		else
		{
			my_entries.put(key, value);
		}
	}

	/** {@inherit_doc} */
	//override
	void set_entry(const int& row, const int& column, const T& value)
	{
		check_row_index(row);
		check_column_index(column);
		if (get_field().get_zero().equals(value))
		{
			my_entries.remove(compute_key(row, column));
		}
		else
		{
			my_entries.put(compute_key(row, column), value);
		}
	}

	/**
	 * {@inherit_doc}
	 *
	 * @ if {@code m} is an
	 * {@code Open_Map_Real_Matrix}, and the total number of entries of the product
	 * is larger than {@code std::numeric_limits<int>::max()}.
	 */
	 //override
	Field_Matrix<T> multiply_transposed(const Field_Matrix<T>& m)
	{
		Matrix_Utils::check_same_column_dimension(*this, m);

		const int out_cols = m.get_row_dimension();
		Field_Matrix<T> out = m.create_matrix(my_rows, out_cols);
		for (Open_Int_To_Field_Hash_Map<T>.Iterator iterator = my_entries.iterator(); iterator.has_next();)
		{
			iterator.advance();
			const T value = iterator.value();
			const int key = iterator.key();
			const int i = key / my_columns;
			const int k = key % my_columns;
			for (int j{}; j < out_cols; ++j)
			{
				out.add_to_entry(i, j, value.multiply(m.get_entry(j, k)));
			}
		}

		return out;
	}

	/**
	 * {@inherit_doc}
	 *
	 * @ if {@code m} is an
	 * {@code Open_Map_Real_Matrix}, and the total number of entries of the product
	 * is larger than {@code std::numeric_limits<int>::max()}.
	 */
	 //override
	Field_Matrix<T> transpose_multiply(const Field_Matrix<T>& m)
	{
		Matrix_Utils::check_same_row_dimension(*this, m);

		const int out_cols = m.get_column_dimension();
		Field_Matrix<T> out = m.create_matrix(my_columns, out_cols);
		for (Open_Int_To_Field_Hash_Map<T>.Iterator iterator = my_entries.iterator(); iterator.has_next();)
		{
			iterator.advance();
			const T value = iterator.value();
			const int key = iterator.key();
			const int k = key / my_columns;
			const int i = key % my_columns;
			for (int j{}; j < out_cols; ++j)
			{
				out.add_to_entry(i, j, value.multiply(m.get_entry(k, j)));
			}
		}

		return out;
	}
};