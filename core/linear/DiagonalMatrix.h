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

  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.exception.;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Math_Utils;
  //import org.hipparchus.util.Precision;
#include "AbstractRealMatrix.h"
#include "MatrixUtils.h"

/**
 * Implementation of a diagonal matrix.
 *
 */
class Diagonal_Matrix : public Abstract_Real_Matrix
{
private:
	/** Entries of the diagonal. */
	std::vector<double> my_data;

	/** Ensure a value is zero.
	 * @param value value to check
	 * @exception  if value is not zero
	 */
	void ensure_zero(const double value)
	{
		/*if (!Precision::equals(0.0, value, 1))
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_LARGE, std::abs(value), 0);
		}*/
	}

public:
	/**
	 * Creates a matrix with the supplied dimension.
	 *
	 * @param dimension Number of rows and columns in the matrix.
	 * @ if the dimension is
	 * not positive.
	 */
	Diagonal_Matrix(const int& dimension)

	{
		//super(dimension, dimension);
		my_data = std::vector<double>(dimension);
	}

	/**
	 * Creates a matrix using the input array as the underlying data.
	 * <br/>
	 * The input array is copied, not referenced.
	 *
	 * @param d Data for the matrix.
	 */
	Diagonal_Matrix(const std::vector<double>& d)
	{
		Diagonal_Matrix(d, true);
	}

	/**
	 * Creates a matrix using the input array as the underlying data.
	 * <br/>
	 * If an array is created specially in order to be embedded in a
	 * this instance and not used directly, the {@code copy_array} may be
	 * set to {@code false}.
	 * This will prevent the copying and improve performance as no new
	 * array will be built and no data will be copied.
	 *
	 * @param d Data for matrix.
	 * @param copy_array if {@code true}, the input array will be copied, * otherwise it will be referenced.
	 * @exception  if d is NULL
	 */
	Diagonal_Matrix(const std::vector<double>& d, const bool& copy_array)
	{
		//Math_Utils::check_not_null(d);
		my_data = copy_array
			? d.clone()
			: d;
	}

	/**
	 * {@inherit_doc}
	 *
	 * @ if the requested dimensions are not equal.
	 */
	 //override
	Real_Matrix create_matrix(const int& row_dimension, const int& column_dimension)
	{
		/*if (row_dimension != column_dimension)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, row_dimension, column_dimension);
		}*/

		return Diagonal_Matrix(row_dimension);
	}

	/** {@inherit_doc} */
	//override
	Real_Matrix copy()
	{
		return Diagonal_Matrix(my_data);
	}

	/**
	 * Compute the sum of {@code this} and {@code m}.
	 *
	 * @param m Matrix to be added.
	 * @return {@code this + m}.
	 * @ if {@code m} is not the same
	 * size as {@code this}.
	 */
	Diagonal_Matrix add(const Diagonal_Matrix& m)

	{
		// Safety check.
		//Matrix_Utils::check_addition_compatible(this, m);

		const int dim = get_row_dimension();
		auto out_data = std::vector<double>(dim);
		for (int i{}; i < dim; i++)
		{
			out_data[i] = my_data[i] + m.get_data_ref()[i];
		}

		return Diagonal_Matrix(out_data, false);
	}

	/**
	 * Returns {@code this} minus {@code m}.
	 *
	 * @param m Matrix to be subtracted.
	 * @return {@code this - m}
	 * @ if {@code m} is not the same
	 * size as {@code this}.
	 */
	Diagonal_Matrix subtract(const Diagonal_Matrix& m)
	{
		//Matrix_Utils::check_subtraction_compatible(this, m);

		const int dim = get_row_dimension();
		auto out_data = std::vector<double>(dim);
		for (int i{}; i < dim; i++)
		{
			out_data[i] = my_data[i] - m.get_data_ref()[i];
		}

		return Diagonal_Matrix(out_data, false);
	}

	/**
	 * Returns the result of postmultiplying {@code this} by {@code m}.
	 *
	 * @param m matrix to postmultiply by
	 * @return {@code this * m}
	 * @ if
	 * {@code column_dimension(this) != row_dimension(m)}
	 */
	Diagonal_Matrix multiply(const Diagonal_Matrix& m)

	{
		Matrix_Utils::check_multiplication_compatible(this, m);

		const int dim = get_row_dimension();
		auto out_data = std::vector<double>(dim);
		for (int i{}; i < dim; i++)
		{
			out_data[i] = my_data[i] * m.get_data_ref()[i];
		}

		return Diagonal_Matrix(out_data, false);
	}

	/** {@inherit_doc} */
	//override
	Real_Matrix multiply(const Real_Matrix& m)
	{
		if (m instanceof Diagonal_Matrix)
		{
			return multiply((Diagonal_Matrix)m);
		}
		Matrix_Utils::check_multiplication_compatible(this, m);
		const Real_Matrix product = m.create_matrix(m.get_row_dimension(), m.get_column_dimension());
		product.walk_in_optimized_order(new Default_Real_Matrix_Changing_Visitor()
			{
				/** {@inherit_doc} */
				//override
				public double visit(const int& row, const int& column, double value)
				{
					return my_data[row] * m.get_entry(row, column);
				}
			});
		return product;
	}

	/**
	 * Returns the result of postmultiplying {@code this} by {@code m^T}.
	 * @param m matrix to first transpose and second postmultiply by
	 * @return {@code this * m}
	 * @ if
	 * {@code column_dimension(this) != column_dimension(m)}
	 * @since 1.3
	 */
	Diagonal_Matrix multiply_transposed(const Diagonal_Matrix& m)
	{
		// transposition is no-op for diagonal matrices
		return multiply(m);
	}

	/** {@inherit_doc} */
	//override
	Real_Matrix multiply_transposed(const Real_Matrix& m)
	{
		if (m instanceof Diagonal_Matrix)
		{
			return multiply_transposed((Diagonal_Matrix)m);
		}

		Matrix_Utils::check_same_column_dimension(this, m);
		const Real_Matrix product = m.create_matrix(m.get_column_dimension(), m.get_row_dimension());
		product.walk_in_optimized_order(new Default_Real_Matrix_Changing_Visitor()
			{
				/** {@inherit_doc} */
				//override
				public double visit(const int& row, const int& column, double value)
				{
					return my_data[row] * m.get_entry(column, row);
				}
			});
		return product;
	}

	/**
	 * Returns the result of postmultiplying {@code this^T} by {@code m}.
	 * @param m matrix to first transpose and second postmultiply by
	 * @return {@code this^T * m}
	 * @ if
	 * {@code column_dimension(this) != column_dimension(m)}
	 * @since 1.3
	 */
	Diagonal_Matrix transpose_multiply(const Diagonal_Matrix& m)
	{
		// transposition is no-op for diagonal matrices
		return multiply(m);
	}

	/** {@inherit_doc} */
	//override
	Real_Matrix transpose_multiply(const Real_Matrix& m)
	{
		if (m instanceof Diagonal_Matrix)
		{
			return transpose_multiply((Diagonal_Matrix)m);
		}

		// transposition is no-op for diagonal matrices
		return multiply(m);
	}

	/** {@inherit_doc} */
	//override
	std::vector<std::vector<double>> get_data()
	{
		const int dim = get_row_dimension();
		auto out = std::vector<std::vector<double>>(dim, std::vector<double>(dim));

		for (int i{}; i < dim; i++)
		{
			out[i][i] = my_data[i];
		}

		return out;
	}

	/**
	 * Gets a reference to the underlying data array.
	 *
	 * @return 1-dimensional array of entries.
	 */
	std::vector<double> get_data_ref() const
	{
		return my_data; // NOPMD - returning an internal array is intentional and documented here
	}

	/** {@inherit_doc} */
	//override
	double get_entry(const int& row, const int column)
	{
		Matrix_Utils::check_matrix_index(this, row, column);
		return row == column
			? my_data[row]
			: 0;
	}

	/** {@inherit_doc}
	 * @ if {@code row != column} and value is non-zero.
	 */
	 //override
	void set_entry(const int& row, const int& column, const double value)
	{
		if (row == column)
		{
			Matrix_Utils::check_row_index(this, row);
			my_data[row] = value;
		}
		else
		{
			ensure_zero(value);
		}
	}

	/** {@inherit_doc}
	 * @ if {@code row != column} and increment is non-zero.
	 */
	 //override
	void add_to_entry(const int& row, const int& column, const double& increment)
	{
		if (row == column)
		{
			Matrix_Utils::check_row_index(this, row);
			my_data[row] += increment;
		}
		else
		{
			ensure_zero(increment);
		}
	}

	/** {@inherit_doc} */
	//override
	void multiply_entry(const int& row, const int& column, const double& factor)
	{
		// we don't care about non-diagonal elements for multiplication
		if (row == column)
		{
			Matrix_Utils::check_row_index(this, row);
			my_data[row] *= factor;
		}
	}

	/** {@inherit_doc} */
	//override
	int get_row_dimension() const
	{
		return my_data.size();
	}

	/** {@inherit_doc} */
	//override
	int get_column_dimension() const
	{
		return my_data.size();
	}

	/** {@inherit_doc} */
	//override
	std::vector<double> operate(const std::vector<double>& v)
	{
		return multiply(new Diagonal_Matrix(v, false)).get_data_ref();
	}

	/** {@inherit_doc} */
	//override
	std::vector<double> pre_multiply(const std::vector<double>& v)
	{
		return operate(v);
	}

	/** {@inherit_doc} */
	//override
	Real_Vector pre_multiply(const Real_Vector v)
	{
		const std::vector<double>& vector_data;
		if (v instanceof Array_Real_Vector)
		{
			vector_data = ((Array_Real_Vector)v).get_data_ref();
		}
		else
		{
			vector_data = v.to_array();
		}
		return Matrix_Utils::create_real__vector(pre_multiply(vector_data));
	}

	/**
	 * Computes the inverse of this diagonal matrix.
	 * <p>
	 * Note: this method will use a singularity threshold of 0, * use {@link #inversestatic_cast<double>(} if a different threshold is needed.
	 *
	 * @return the inverse of {@code m}
	 * @ if the matrix is singular
	 */
	Diagonal_Matrix inverse()
	{
		return inverse(0);
	}

	/**
	 * Computes the inverse of this diagonal matrix.
	 *
	 * @param threshold Singularity threshold.
	 * @return the inverse of {@code m}
	 * @ if the matrix is singular
	 */
	Diagonal_Matrix inverse(double threshold)
	{
		if (is_singular(threshold))
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::SINGULAR_MATRIX);
		}

		auto result = std::vector<double>(my_data.size()];
		for (int i{}; i < my_data.size(); i++)
		{
			result[i] = 1.0 / my_data[i];
		}
		return Diagonal_Matrix(result, false);
	}

	/** Returns whether this diagonal matrix is singular, i.e. any diagonal entry
	 * is equal to {@code 0} within the given threshold.
	 *
	 * @param threshold Singularity threshold.
	 * @return {@code true} if the matrix is singular, {@code false} otherwise
	 */
	bool is_singular(double threshold)
	{
		for (int i{}; i < my_data.size(); i++)
		{
			if (Precision::equals(my_data[i], 0.0, threshold))
			{
				return true;
			}
		}
		return false;
	}
};