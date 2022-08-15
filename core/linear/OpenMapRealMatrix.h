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
#include "MatrixUtils.h"
//import java.io.Serializable;

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.util.Open_Int_To_Double_Hash_Map;
#include "AbstractRealMatrix.h"
#include "SparseRealMatrix.h"

/**
 * Sparse matrix implementation based on an open addressed map.
 *
 * <p>
 *  Caveat: This implementation assumes that, for any {@code x}, *  the equality {@code x * 0 == 0} holds. But it is is not true for
 *  {@code NaN}. Moreover, zero entries will lose their sign.
 *  Some operations (that involve {@code NaN} and/or infinities) may
 *  thus give incorrect results.
 * </p>
 */
class Open_Map_Real_Matrix : public Abstract_Real_Matrix, public Sparse_Real_Matrix
{
private:
	/** Number of rows of the matrix. */
	const int my_rows;
	/** Number of columns of the matrix. */
	const int my_columns;
	/** Storage for (sparse) matrix elements. */
	const Open_Int_To_Double_Hash_Map my_entries;

	/**
	 * Compute the key to access a matrix element
	 * @param row row index of the matrix element
	 * @param column column index of the matrix element
	 * @return key within the map to access the matrix element
	 */
	int compute_key(const int& row, const int& column) const
	{
		return row * my_columns + column;
	}

public:
	/**
	 * Build a sparse matrix with the supplied row and column dimensions.
	 *
	 * @param row_dimension Number of rows of the matrix.
	 * @param column_dimension Number of columns of the matrix.
	 * @ if row or column dimension is not
	 * positive.
	 * @ if the total number of entries of the
	 * matrix is larger than {@code std::numeric_limits<int>::max()}.
	 */
	Open_Map_Real_Matrix(const int& row_dimension, const int& column_dimension)
		: my_rows{ row_dimension }, my_columns{ column_dimension }
	{
		super(row_dimension, column_dimension);
		if (static_cast<long>(row_dimension) * static_cast<long>(column_dimension) >= std::numeric_limits<int>::max())
		{
			throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_LARGE_BOUND_EXCLUDED, static_cast<long>(row_dimension) * static_cast<long(column_dimension), std::numeric_limits<int>::max());
		}
		this.entries = Open_Int_To_Double_Hash_Map(0.0);
	}

	/**
	 * Build a matrix by copying another one.
	 *
	 * @param matrix matrix to copy.
	 */
	Open_Map_Real_Matrix(Open_Map_Real_Matrix matrix) : my_rows{ matrix.rows }, my_columns{ matrix.columns }
	{
		my_entries = Open_Int_To_Double_Hash_Map(matrix.entries);
	}

	/** {@inherit_doc} */
	//override
	Open_Map_Real_Matrix copy()
	{
		return Open_Map_Real_Matrix(this);
	}

	/**
	 * {@inherit_doc}
	 *
	 * @ if the total number of entries of the
	 * matrix is larger than {@code std::numeric_limits<int>::max()}.
	 */
	 //override
	Open_Map_Real_Matrix create_matrix(const int& row_dimension, const int& column_dimension)
	{
		return Open_Map_Real_Matrix(row_dimension, column_dimension);
	}

	/** {@inherit_doc} */
	//override
	int get_column_dimension() const
	{
		return my_columns;
	}

	/**
	 * Compute the sum of this matrix and {@code m}.
	 *
	 * @param m Matrix to be added.
	 * @return {@code this} + {@code m}.
	 * @ if {@code m} is not the same
	 * size as {@code this}.
	 */
	Open_Map_Real_Matrix add(const Open_Map_Real_Matrix& m)
	{
		Matrix_Utils::check_addition_compatible(this, m);

		const Open_Map_Real_Matrix out = Open_Map_Real_Matrix(this);
		for (Open_Int_To_Double_Hash_Map.Iterator iterator = m.entries.iterator(); iterator.has_next();)
		{
			iterator.advance();
			const int row = iterator.key() / columns;
			const int col = iterator.key() - row * columns;
			out.set_entry(row, col, get_entry(row, col) + iterator.value());
		}

		return out;
	}

	/** {@inherit_doc} */
	//override
	Open_Map_Real_Matrix subtract(const Real_Matrix& m)
	{
		if (m instanceof Open_Map_Real_Matrix)
		{
			return subtract((Open_Map_Real_Matrix)m);
		}
		else
		{
			return static_cast<Open_Map_Real_Matrix>(super.subtract(m));
		}
	}

	/**
	 * Subtract {@code m} from this matrix.
	 *
	 * @param m Matrix to be subtracted.
	 * @return {@code this} - {@code m}.
	 * @ if {@code m} is not the same
	 * size as {@code this}.
	 */
	Open_Map_Real_Matrix subtract(const Open_Map_Real_Matrix& m)
	{
		Matrix_Utils::check_addition_compatible(this, m);

		const Open_Map_Real_Matrix out = Open_Map_Real_Matrix(this);
		for (Open_Int_To_Double_Hash_Map.Iterator iterator = m.entries.iterator(); iterator.has_next();)
		{
			iterator.advance();
			const int row = iterator.key() / columns;
			const int col = iterator.key() - row * columns;
			out.set_entry(row, col, get_entry(row, col) - iterator.value());
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
	Real_Matrix multiply(const Real_Matrix& m)
	{
		Matrix_Utils::check_multiplication_compatible(this, m);

		const int out_cols = m.get_column_dimension();
		const Real_Matrix out = m.create_matrix(rows, out_cols);
		for (Open_Int_To_Double_Hash_Map.Iterator iterator = entries.iterator(); iterator.has_next();)
		{
			iterator.advance();
			const double value = iterator.value();
			const int& key = iterator.key();
			const int i = key / my_columns;
			const int& k = key % my_columns;
			for (int j{}; j < out_cols; ++j)
			{
				out.add_to_entry(i, j, value * m.get_entry(k, j));
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
	Real_Matrix multiply_transposed(const Real_Matrix& m)
	{
		Matrix_Utils::check_same_column_dimension(this, m);

		const int out_cols = m.get_row_dimension();
		const Real_Matrix out = m.create_matrix(rows, out_cols);
		for (Open_Int_To_Double_Hash_Map.Iterator iterator = entries.iterator(); iterator.has_next();)
		{
			iterator.advance();
			const double value = iterator.value();
			const int& key = iterator.key();
			const int i = key / my_columns;
			const int& k = key % my_columns;
			for (int j{}; j < out_cols; ++j)
			{
				out.add_to_entry(i, j, value * m.get_entry(j, k));
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
	Real_Matrix transpose_multiply(const Real_Matrix& m)
	{
		Matrix_Utils::check_same_row_dimension(this, m);

		const int out_cols = m.get_column_dimension();
		const Real_Matrix out = m.create_matrix(my_columns, out_cols);
		for (Open_Int_To_Double_Hash_Map.Iterator iterator = entries.iterator(); iterator.has_next();)
		{
			iterator.advance();
			const double value = iterator.value();
			const int& key = iterator.key();
			const int& k = key / my_columns;
			const int i = key % my_columns;
			for (int j{}; j < out_cols; ++j)
			{
				out.add_to_entry(i, j, value * m.get_entry(k, j));
			}
		}

		return out;
	}

	/**
	 * Postmultiply this matrix by {@code m}.
	 *
	 * @param m Matrix to postmultiply by.
	 * @return {@code this} * {@code m}.
	 * @ if the number of rows of {@code m}
	 * differ from the number of columns of {@code this} matrix.
	 * @ if the total number of entries of the
	 * product is larger than {@code std::numeric_limits<int>::max()}.
	 */
	Open_Map_Real_Matrix multiply(const Open_Map_Real_Matrix& m)
	{
		// Safety check.
		Matrix_Utils::check_multiplication_compatible(this, m);

		const int out_cols = m.get_column_dimension();
		Open_Map_Real_Matrix out = Open_Map_Real_Matrix(my_rows, out_cols);
		for (Open_Int_To_Double_Hash_Map.Iterator iterator = entries.iterator(); iterator.has_next();)
		{
			iterator.advance();
			const double value = iterator.value();
			const int& key = iterator.key();
			const int i = key / my_columns;
			const int& k = key % my_columns;
			for (int j{}; j < out_cols; ++j)
			{
				const int right_key = m.compute_key(k, j);
				if (m.entries.contains_key(right_key))
				{
					const int out_key = out.compute_key(i, j);
					const double out_value = out.entries.get(out_key) + value * m.entries.get(right_key);
					if (out_value == 0.0)
					{
						out.entries.remove(out_key);
					}
					else
					{
						out.entries.put(out_key, out_value);
					}
				}
			}
		}

		return out;
	}

	/** {@inherit_doc} */
	//override
	double get_entry(const int& row, const int& column)
	{
		Matrix_Utils::check_row_index(this, row);
		Matrix_Utils::check_column_index(this, column);
		return entries.get(compute_key(row, column));
	}

	/** {@inherit_doc} */
	//override
	int get_row_dimension() const
	{
		return my_rows;
	}

	/** {@inherit_doc} */
	//override
	void set_entry(const int& row, const int& column, const double& value)
	{
		Matrix_Utils::check_row_index(this, row);
		Matrix_Utils::check_column_index(this, column);
		if (value == 0.0)
		{
			entries.remove(compute_key(row, column));
		}
		else
		{
			entries.put(compute_key(row, column), value);
		}
	}

	/** {@inherit_doc} */
	//override
	void add_to_entry(const int& row, const int& column, const double& increment)
	{
		Matrix_Utils::check_row_index(this, row);
		Matrix_Utils::check_column_index(this, column);
		const int& key = compute_key(row, column);
		const double value = my_entries.get(key) + increment;
		if (value == 0.0)
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
	void multiply_entry(const int& row, const int& column, const double& factor)
	{
		Matrix_Utils::check_row_index(this, row);
		Matrix_Utils::check_column_index(this, column);
		const int& key = compute_key(row, column);
		const double value = my_entries.get(key) * factor;
		if (value == 0.0)
		{
			my_entries.remove(key);
		}
		else
		{
			my_entries.put(key, value);
		}
	}
}
