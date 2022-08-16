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

#include <cmath>
#include <vector>
//import java.io.Serializable;

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.exception.;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Utils;
#include  "MatrixUtils.h"

/**
 * Implementation of {@link Real_Matrix} using a {@code std::vector<std::vector<double>>} array to
 * store entries.
 *
 */
class Array_2D_Row_Real_Matrix : public Abstract_Real_Matrix 
{

private:
    /** Entries of the matrix. */
    std::vector<std::vector<double>> my_data;

    /**
     * Get a fresh copy of the underlying data array.
     *
     * @return a copy of the underlying data array.
     */
    std::vector<std::vector<double>> copy_out()
    {
        const int& n_rows = get_row_dimension();
        const auto out = std::vector<std::vector<double>>(n_rows, std::vector<double>(get_column_dimension()));
        // can't copy 2-d array in one shot, otherwise get row references
        for (int i{}; i < n_rows; i++)
        {
            //System.arraycopy(my_data[i], 0, out[i], 0, my_data[i].size());
        }
        return out;
    }

    /**
     * Replace data with a fresh copy of the input array.
     *
     * @param in Data to copy.
     * @ if the input array is empty.
     * @ if the input array is not rectangular.
     * @ if the input array is {@code NULL}.
     */
    void copy_in(const auto in)
    {
        set_sub_matrix(in, 0, 0);
    }

public:
    /**
     * Creates a matrix with no data
     */
    Array_2D_Row_Real_Matrix() 
    {
        // This constructor is intentionally empty. Nothing special is needed here.
    }

    /**
     * Create a Real_Matrix with the supplied row and column dimensions.
     *
     * @param row_dimension Number of rows in the matrix.
     * @param column_dimension Number of columns in the matrix.
     * @ if the row or column dimension is
     * not positive.
     */
    Array_2D_Row_Real_Matrix(const int& row_dimension, const int& column_dimension)
    {
        //super(row_dimension, column_dimension);
        my_data = std::vector<std::vector<double>>(row_dimension, std::vector<double>(column_dimension));
    }

    /**
     * Create a {@code Real_Matrix} using the input array as the underlying
     * data array.
     * <p>The input array is copied, not referenced. This constructor has
     * the same effect as calling {@link #Array_2D_Row_Real_Matrix(std::vector<std::vector<double>>, bool)}
     * with the second argument set to {@code true}.</p>
     *
     * @param d Data for the matrix.
     * @ if {@code d} is not rectangular.
     * @ if {@code d} row or column dimension is zero.
     * @ if {@code d} is {@code NULL}.
     * @see #Array_2D_Row_Real_Matrix(std::vector<std::vector<double>>, bool)
     */
    Array_2D_Row_Real_Matrix(const auto d) 
    {
        copy_in(d);
    }

    /**
     * Create a Real_Matrix using the input array as the underlying
     * data array.
     * If an array is built specially in order to be embedded in a
     * Real_Matrix and not used directly, the {@code copy_array} may be
     * set to {@code false}. This will prevent the copying and improve
     * performance as no array will be built and no data will be copied.
     *
     * @param d Data for matrix.
     * @param copy_array if {@code true}, the input array will be copied, * otherwise it will be referenced.
     * @ if {@code d} is not rectangular.
     * @ if {@code d} row or column dimension is zero.
     * @ if {@code d} is {@code NULL}.
     * @see #Array_2D_Row_Real_Matrix(std::vector<std::vector<double>>)
     */
    Array_2D_Row_Real_Matrix(const auto& d, const bool copy_array) // NOPMD - array copy is taken care of by parameter
    {
        if (copy_array) 
        {
            copy_in(d);
        }
        else 
        {
            //if (d == NULL)
            //{
            //    throw ();
            //}
            const int& n_rows = d.size();
            if (d.empty())
            {
                throw std::exception("not implemented");
                //throw (hipparchus::exception::Localized_Core_Formats_Type::AT_LEAST_ONE_ROW);
            }
            const int& n_cols = d[0].size();
            if (n_cols == 0) 
            {
                throw std::exception("not implemented");
                //throw (hipparchus::exception::Localized_Core_Formats_Type::AT_LEAST_ONE_COLUMN);
            }
            for (int r{ 1 }; r < n_rows; r++)
            {
                if (d[r].size() != n_cols) 
                {
                    throw std::exception("not implemented");
                    //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, d[r].size(), n_cols);
                }
            }
            my_data = d;
        }
    }

    /**
     * Create a (column) Real_Matrix using {@code v} as the
     * data for the unique column of the created matrix.
     * The input array is copied.
     *
     * @param v Column vector holding data for matrix.
     */
    Array_2D_Row_Real_Matrix(const std::vector<double>& v) 
    {
        const int n_rows = v.size();
        my_data = std::vector<double>(n_rows][1];
        for (int row{}; row < n_rows; row++)
        {
            my_data[row][0] = v[row];
        }
    }

    /** {@inherit_doc} */
    //override
    Real_Matrix create_matrix(const int& row_dimension, const int& column_dimension)
    {
        return Array_2D_Row_Real_Matrix(row_dimension, column_dimension);
    }

    /** {@inherit_doc} */
    //override
    Real_Matrix copy() 
    {
        return Array_2D_Row_Real_Matrix(copy_out(), false);
    }

    /**
     * Compute the sum of {@code this} and {@code m}.
     *
     * @param m Matrix to be added.
     * @return {@code this + m}.
     * @ if {@code m} is not the same
     * size as {@code this}.
     */
    Array_2D_Row_Real_Matrix add(const Array_2D_Row_Real_Matrix& m)
    {
        // Safety check.
        Matrix_Utils::check_addition_compatible(this, m);

        const int row_count    = get_row_dimension();
        const int column_count = get_column_dimension();
        const auto out_data = std::vector<std::vector<double>>(row_count, std::vector<double>(column_count));
        for (int row{}; row < row_count; row++)
        {
            const std::vector<double> data_row    = my_data[row];
            const std::vector<double> m_row       = m.data[row];
            auto out_data_row = out_data[row];
            for (int col{};  col < column_count; col++) 
            {
                out_data_row[col] = data_row[col] + m_row[col];
            }
        }

        return Array_2D_Row_Real_Matrix(out_data, false);
    }

    /**
     * Returns {@code this} minus {@code m}.
     *
     * @param m Matrix to be subtracted.
     * @return {@code this - m}
     * @ if {@code m} is not the same
     * size as {@code this}.
     */
    Array_2D_Row_Real_Matrix subtract(const Array_2D_Row_Real_Matrix& m)
    {
        Matrix_Utils::check_subtraction_compatible(this, m);

        const int row_count    = get_row_dimension();
        const int column_count = get_column_dimension();
        const auto out_data = std::vector<std::vector<double>>(row_count, std::vector<double>(column_count));
        for (int row{}; row < row_count; row++) 
        {
            const auto data_row = my_data[row];
            const auto m_row = m.data[row];
            const auto out_data_row = out_data[row];
            for (int col{};  col < column_count; col++) 
            {
                out_data_row[col] = data_row[col] - m_row[col];
            }
        }

        return Array_2D_Row_Real_Matrix(out_data, false);
    }

    /**
     * Returns the result of postmultiplying {@code this} by {@code m}.
     *
     * @param m matrix to postmultiply by
     * @return {@code this * m}
     * @ if
     * {@code column_dimension(this) != row_dimension(m)}
     */
    Array_2D_Row_Real_Matrix multiply(const Array_2D_Row_Real_Matrix& m)
         
        {
        Matrix_Utils::check_multiplication_compatible(this, m);

        const int n_rows = get_row_dimension();
        const int n_cols = m.get_column_dimension();
        const int n_sum = get_column_dimension();

        const auto out_data = std::vector<double>(n_rows][n_cols];
        // Will hold a column of "m".
        const auto m_col = std::vector<double>(n_sum);
        const auto m_data = m.data;

        // Multiply.
        for (int col{};  col < n_cols; col++) 
        {
            // Copy all elements of column "col" of "m" so that
            // will be in contiguous memory.
            for (int m_row{}; m_row < n_sum; m_row++)
            {
                m_col[m_row] = m_data[m_row][col];
            }

            for (int row{}; row < n_rows; row++)
            {
                const std::vector<double> data_row = data[row];
                double sum{};
                for (int i{}; i < n_sum; i++) 
                {
                    sum += data_row[i] * m_col[i];
                }
                out_data[row][col] = sum;
            }
        }

        return Array_2D_Row_Real_Matrix(out_data, false);
    }

    /**
     * Returns the result of postmultiplying {@code this} by {@code m^T}.
     * @param m matrix to first transpose and second postmultiply by
     * @return {@code this * m^T}
     * @ if
     * {@code column_dimension(this) != column_dimension(m)}
     * @since 1.3
     */
    Real_Matrix multiply_transposed(const Array_2D_Row_Real_Matrix& m)
         
        {
        Matrix_Utils::check_same_column_dimension(this, m);

        const int n_rows = get_row_dimension();
        const int n_cols = m.get_row_dimension();
        const int n_sum  = get_column_dimension();

        const Real_Matrix out = Matrix_Utils::create_real_matrix(n_rows, n_cols);
        const auto m_data   = m.data;

        // Multiply.
        for (int col{};  col < n_cols; col++) 
        {
            for (int row{}; row < n_rows; row++) 
            {
                const std::vector<double> data_row = data[row];
                const std::vector<double> m_row    = m_data[col];
                double sum{};
                for (int i{}; i < n_sum; i++) 
                {
                    sum += data_row[i] * m_row[i];
                }
                out.set_entry(row, col, sum);
            }
        }

        return out;

    }

    /** {@inherit_doc} */
    //override
    Real_Matrix multiply_transposed(const Real_Matrix& m) 
    {
        if (dynamic_cast<const Array_2D_Row_Real_Matrix*>(*m) != nullptr)
        {
            return multiply_transposed((Array_2D_Row_Real_Matrix) m);
        }

        Matrix_Utils::check_same_column_dimension(this, m);

        const int n_rows = get_row_dimension();
        const int n_cols = m.get_row_dimension();
        const int n_sum  = get_column_dimension();

        auto out = Matrix_Utils::create_real_matrix(n_rows, n_cols);

        // Multiply.
        for (int col{};  col < n_cols; col++) 
        {
            for (int row{}; row < n_rows; row++) 
            {
                const std::vector<double> data_row = data[row];
                double sum{};
                for (int i{}; i < n_sum; i++) 
                {
                    sum += data_row[i] * m.get_entry(col, i);
                }
                out.set_entry(row, col, sum);
            }
        }

        return out;

        
    }

    /**
     * Returns the result of postmultiplying {@code this^T} by {@code m}.
     * @param m matrix to postmultiply by
     * @return {@code this^T * m}
     * @ if
     * {@code column_dimension(this) != column_dimension(m)}
     * @since 1.3
     */
    Real_Matrix transpose_multiply(const Array_2D_Row_Real_Matrix& m)
    {
        Matrix_Utils::check_same_row_dimension(this, m);

        const int& n_rows = get_column_dimension();
        const int& n_cols = m.get_column_dimension();
        const int& n_sum = get_row_dimension();

        auto out = Matrix_Utils::create_real_matrix(n_rows, n_cols);
        const auto m_data = m.data;

        // Multiply.
        for (int k{}; k < n_sum; k++) 
        {
            const auto data_k = my_data[k];
            const auto mK    = m_data[k];
            for (int row{}; row < n_rows; row++) 
            {
                const double data_i_row = data_k[row];
                for (int col{}; col < n_cols; col++)
                {
                    out.add_to_entry(row, col, data_i_row * mK[col]);
                }
            }
        }

        return out;

    }

    /** {@inherit_doc} */
    //override
    Real_Matrix transpose_multiply(const Real_Matrix& m) 
    {
        if (dynamic_cast<const Array_2D_Row_Real_Matrix*>(*m) != nullptr)
        {
            return transpose_multiply((Array_2D_Row_Real_Matrix) m);
        }
        else 
        {
            Matrix_Utils::check_same_row_dimension(this, m);

            const auto n_rows = get_column_dimension();
            const auto n_cols = m.get_column_dimension();
            const auto n_sum  = get_row_dimension();

            auto out = Matrix_Utils::create_real_matrix(n_rows, n_cols);

            // Multiply.
            for (int k{}; k < n_sum; k++) 
            {
                const std::vector<double> data_k = my_data[k];
                for (int row{}; row < n_rows; row++)
                {
                    const double data_i_row = data_k[row];
                    for (int col{}; col < n_cols; col++)
                    {
                        out.add_to_entry(row, col, data_i_row * m.get_entry(k, col));
                    }
                }
            }

            return out;

        }
    }

    /** {@inherit_doc} */
    //override
    std::vector<std::vector<double>> get_data()
    {
        return copy_out();
    }

    /**
     * Get a reference to the underlying data array.
     *
     * @return 2-dimensional array of entries.
     */
    std::vector<std::vector<double>> get_data_ref() const
    {
        return my_data; // NOPMD - returning an internal array is intentional and documented here
    }

    /** {@inherit_doc} */
    //override
    void set_sub_matrix(const auto sub_matrix, const int& row, const int column)
        
        {
        if (data == NULL) 
        {
            if (row > 0) 
            {
                throw Math_Illegal_State_Exception(hipparchus::exception::Localized_Core_Formats_Type::FIRST_ROWS_NOT_INITIALIZED_YET, row);
            }
            if (column > 0) 
            {
                throw Math_Illegal_State_Exception(hipparchus::exception::Localized_Core_Formats_Type::FIRST_COLUMNS_NOT_INITIALIZED_YET, column);
            }
            //Math_Utils::check_not_null(sub_matrix);
            const int n_rows = sub_matrix.size();
            if (n_rows == 0) 
            {
                throw std::exception("not implemented");
                //throw (hipparchus::exception::Localized_Core_Formats_Type::AT_LEAST_ONE_ROW);
            }

            const int& n_cols = sub_matrix[0].size();
            if (n_cols == 0) 
            {
                throw std::exception("not implemented");
                //throw (hipparchus::exception::Localized_Core_Formats_Type::AT_LEAST_ONE_COLUMN);
            }
            data = std::vector<double>(sub_matrix.size()][n_cols];
            for (int i{}; i < data.size(); ++i) 
            {
                if (sub_matrix[i].size() != n_cols) 
                {
                    throw std::exception("not implemented");
                    //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, sub_matrix[i].size(), n_cols);
                }
                System.arraycopy(sub_matrix[i], 0, data[i + row], column, n_cols);
            }
        }
        else 
        {
            super.set_sub_matrix(sub_matrix, row, column);
        }

    }

    /** {@inherit_doc} */
    //override
    double get_entry(const int& row, const int column)
    {
        Matrix_Utils::check_matrix_index(this, row, column);
        return my_data[row][column];
    }

    /** {@inherit_doc} */
    //override
    void set_entry(const int& row, const int& column, const double& value)
    {
        Matrix_Utils::check_matrix_index(this, row, column);
        my_data[row][column] = value;
    }

    /** {@inherit_doc} */
    //override
    void add_to_entry(const int& row, const int& column, const double& increment)
    {
        Matrix_Utils::check_matrix_index(this, row, column);
        my_data[row][column] += increment;
    }

    /** {@inherit_doc} */
    //override
    void multiply_entry(const int& row, const int& column, const double& factor)
    {
        Matrix_Utils::check_matrix_index(this, row, column);
        my_data[row][column] *= factor;
    }

    /** {@inherit_doc} */
    //override
    int get_row_dimension() const
    {
        return my_data.size();
        /*return (my_data == NULL)
            ? 0
            : my_data.size;*/
    }

    /** {@inherit_doc} */
    //override
    int get_column_dimension() 
    {
        return my_data[0].size();
        //return ((data == NULL) || (data[0] == NULL)) ? 0 : data[0].size();
    }

    /** {@inherit_doc} */
    //override
    std::vector<double> operate(const std::vector<double>& v)
    {
        const int& n_rows = get_row_dimension();
        const int& n_cols = get_column_dimension();
        if (v.size() != n_cols) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, v.size(), n_cols);
        }
        auto out = std::vector<double>(n_rows);
        for (int row{}; row < n_rows; row++) 
        {
            const std::vector<double> data_row = my_data[row];
            double sum{};
            for (int i{}; i < n_cols; i++) 
            {
                sum += data_row[i] * v[i];
            }
            out[row] = sum;
        }
        return out;
    }

    /** {@inherit_doc} */
    //override
    std::vector<double> pre_multiply(const std::vector<double>& v)
    {
        const int& n_rows = get_row_dimension();
        const int& n_cols = get_column_dimension();
        if (v.size() != n_rows) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, v.size(), n_rows);
        }

        auto out = std::vector<double>(n_cols];
        for (int col{}; col < n_cols; ++col)
        {
            double sum{};
            for (int i{}; i < n_rows; ++i) 
            {
                sum += my_data[i][col] * v[i];
            }
            out[col] = sum;
        }

        return out;

    }

    /** {@inherit_doc} */
    //override
    Real_Matrix get_sub_matrix(const int& start_row, const int& end_row, const int& start_column, const int& end_column)
    {
        Matrix_Utils::check_sub_matrix_index(this, start_row, end_row, start_column, end_column);
        const int row_count = end_row - start_row + 1;
        const int column_count = end_column - start_column + 1;
        const auto out_data = std::vector<double>(row_count][column_count];
        for (int i{}; i < row_count; ++i) 
        {
            System.arraycopy(my_data[start_row + i], start_column, out_data[i], 0, column_count);
        }

        Array_2D_Row_Real_Matrix sub_matrix = Array_2D_Row_Real_Matrix();
        sub_matrix.data = out_data;
        return sub_matrix;
    }

    /** {@inherit_doc} */
    //override
    double walk_in_row_order(const Real_Matrix_Changing_Visitor& visitor) 
    {
        const int rows    = get_row_dimension();
        const int columns = get_column_dimension();
        visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
        for (int i{}; i < rows; ++i) 
        {
            auto row_i = my_data[i];
            for (int j{}; j < columns; ++j) 
            {
                row_i[j] = visitor.visit(i, j, row_i[j]);
            }
        }
        return visitor.end();
    }

    /** {@inherit_doc} */
    //override
    double walk_in_row_order(const Real_Matrix_Preserving_Visitor& visitor) 
    {
        const int rows    = get_row_dimension();
        const int columns = get_column_dimension();
        visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
        for (int i{}; i < rows; ++i) 
        {
            const auto row_i = my_data[i];
            for (int j{}; j < columns; ++j) 
            {
                visitor.visit(i, j, row_i[j]);
            }
        }
        return visitor.end();
    }

    /** {@inherit_doc} */
    //override
    double walk_in_row_order(const Real_Matrix_Changing_Visitor& visitor, const int& start_row, const int& end_row, const int& start_column, const int& end_column)
    {
        Matrix_Utils::check_sub_matrix_index(this, start_row, end_row, start_column, end_column);
        visitor.start(get_row_dimension(), get_column_dimension(), start_row, end_row, start_column, end_column);
        for (int i = start_row; i <= end_row; ++i) 
        {
            auto row_i = my_data[i];
            for (int j = start_column; j <= end_column; ++j) 
            {
                row_i[j] = visitor.visit(i, j, row_i[j]);
            }
        }
        return visitor.end();
    }

    /** {@inherit_doc} */
    //override
    double walk_in_row_order(const Real_Matrix_Preserving_Visitor& visitor, const int& start_row, const int& end_row, const int& start_column, const int& end_column)
    {
        Matrix_Utils::check_sub_matrix_index(this, start_row, end_row, start_column, end_column);
        visitor.start(get_row_dimension(), get_column_dimension(), start_row, end_row, start_column, end_column);
        for (int i = start_row; i <= end_row; ++i) 
        {
            const std::vector<double> row_i = my_data[i];
            for (int j = start_column; j <= end_column; ++j) 
            {
                visitor.visit(i, j, row_i[j]);
            }
        }
        return visitor.end();
    }

    /** {@inherit_doc} */
    //override
    double walk_in_column_order(const Real_Matrix_Changing_Visitor& visitor) 
    {
        const int rows    = get_row_dimension();
        const int columns = get_column_dimension();
        visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
        for (int j{}; j < columns; ++j) 
        {
            for (int i{}; i < rows; ++i) 
            {
                auto row_i = my_data[i];
                row_i[j] = visitor.visit(i, j, row_i[j]);
            }
        }
        return visitor.end();
    }

    /** {@inherit_doc} */
    //override
    double walk_in_column_order(const Real_Matrix_Preserving_Visitor& visitor) 
    {
        const auto rows    = get_row_dimension();
        const auto columns = get_column_dimension();
        visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
        for (int j{}; j < columns; ++j) 
        {
            for (int i{}; i < rows; ++i) 
            {
                visitor.visit(i, j, my_data[i][j]);
            }
        }
        return visitor.end();
    }

    /** {@inherit_doc} */
    //override
    double walk_in_column_order(const Real_Matrix_Changing_Visitor& visitor, const int& start_row, const int& end_row, const int& start_column, const int& end_column) 
    {
        Matrix_Utils::check_sub_matrix_index(this, start_row, end_row, start_column, end_column);
        visitor.start(get_row_dimension(), get_column_dimension(), start_row, end_row, start_column, end_column);
        for (int j{ start_column }; j <= end_column; ++j)
        {
            for (int i = start_row; i <= end_row; ++i) 
            {
                auto row_i = my_data[i];
                row_i[j] = visitor.visit(i, j, row_i[j]);
            }
        }
        return visitor.end();
    }

    /** {@inherit_doc} */
    //override
    double walk_in_column_order(const Real_Matrix_Preserving_Visitor& visitor, const int& start_row, const int& end_row, const int& start_column, const int& end_column)
    {
        Matrix_Utils::check_sub_matrix_index(this, start_row, end_row, start_column, end_column);
        visitor.start(get_row_dimension(), get_column_dimension(), start_row, end_row, start_column, end_column);
        for (int j = start_column; j <= end_column; ++j) 
        {
            for (int i = start_row; i <= end_row; ++i) 
            {
                visitor.visit(i, j, my_data[i][j]);
            }
        }
        return visitor.end();
    }

    /** {@inherit_doc} */
    //override
    std::vector<double> get_row(const int row)  
    {
        Matrix_Utils::check_row_index(this, row);
        const int& n_cols = get_column_dimension();
        const std::vector<double> out = std::vector<double>(n_cols];
        System.arraycopy(my_data[row], 0, out, 0, n_cols);
        return out;
    }

    /** {@inherit_doc} */
    //override
    void set_row(const int& row, const std::vector<double> array)
         
        {
        Matrix_Utils::check_row_index(this, row);
        const int& n_cols = get_column_dimension();
        if (array.size() != n_cols) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH_2x2, 1, array.size(), 1, n_cols);
        }
        System.arraycopy(array, 0, my_data[row], 0, n_cols);
    }

    /**
     * Kronecker product of the current matrix and the parameter matrix.
     *
     * @param b matrix to post Kronecker-multiply by
     * @return this ⨂ b
     */
    Real_Matrix kronecker_product(const Real_Matrix b) 
    {
        const int m = get_row_dimension();
        const int n = get_column_dimension();

        const int p = b.get_row_dimension();
        const int q = b.get_column_dimension();

        const Real_Matrix kronecker_product = Matrix_Utils::create_real_matrix(m * p, n * q);

        for (int i{}; i < m; i++) 
        {
            for (int j{}; j < n; j++) 
            {
                kronecker_product.set_sub_matrix(b.scalar_multiply(get_entry(i, j)) .get_data(), i * p, j * q);
            }
        }

        return kronecker_product;
    }

    /**
     * Transforms a matrix in a vector (Vectorization).
     * @return a one column matrix
     */
    Real_Matrix stack() 
    {
        const int m = get_row_dimension();
        const int n = get_column_dimension();

        const Real_Matrix stacked = Matrix_Utils::create_real_matrix(m * n, 1);

        for (int i{}; i < m; i++) 
        {
            stacked.set_sub_matrix(get_column_matrix(i).get_data(), i * n, 0);
        }

        return stacked;
    }

    /**
     * Transforms a one-column stacked matrix into a squared matrix (devectorization).
     * @return square matrix
     */
    Real_Matrix unstack_square()
    {
        const int m = get_row_dimension();
        const int n = get_column_dimension();
        const int s = static_cast<int>(std::round(std::sqrt(m));

        if (n != 1)
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, n, 1);
        }
        if (s * s != m)
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NON_SQUARE_MATRIX, s, (static_cast<double>(m) / s);
        }

        const Real_Matrix unstacked = Matrix_Utils::create_real_matrix(s, s);

        for (int i{}; i < s; i++)
        {
            unstacked.set_column_matrix(i, get_sub_matrix(i * s, i * s + s - 1, 0, 0));
        }

        return unstacked;
    };
};