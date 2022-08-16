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
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.exception.Null_Argument_Exception;
//import org.hipparchus.util.Math_Arrays;
//import org.hipparchus.util.Math_Utils;
#include "MatrixUtils.h"
#include "../FieldElement.h"
/**
 * Implementation of Field_Matrix<T> using a {@link Field_Element}[][] array to store entries.
 * <p>
 * As specified in the {@link Field_Matrix} interface, matrix element indexing
 * is 0-based -- e.g., <code>get_entry(0, 0)</code>
 * returns the element in the first row, first column of the matrix.</li></ul>
 * </p>
 *
 * @param <T> the type of the field elements
 */
template<typename T, typename std::enable_if<std::is_base_of<Field_Element<T>, T>::value>::type* = nullptr>
class Array2DRowField_Matrix : public Abstract_Field_Matrix<T>
{
private:
    /** Entries of the matrix */
    std::vector<std::vector<T>> data;

    /**
     * Get a fresh copy of the underlying data array.
     *
     * @return a copy of the underlying data array.
     */
    std::vector<std::vector<T>> copy_out()
    {
        const auto n_rows = this.get_row_dimension();
        const std::vector<std::vector<T>> out = Math_Arrays::build_array(get_field(), n_rows, get_column_dimension());
        // can't copy 2-d array in one shot, otherwise get row references
        for (int i{}; i < n_rows; i++)
        {
            System.arraycopy(data[i], 0, out[i], 0, data[i].size());
        }
        return out;
    }

    /**
     * Replace data with a fresh copy of the input array.
     *
     * @param in Data to copy.
     * @ if the input array is empty.
     * @ if the input array is not rectangular.
     * @Null_Argument_Exception if the input array is {@code NULL}.
     */
    void copy_in(const std::vector<std::vector<T>> in)
        , Null_Argument_Exception
    {
    set_sub_matrix(in, 0, 0);
    }

public:
    /**
     * Creates a matrix with no data
     * @param field field to which the elements belong
     */
    Array2DRowField_Matrix(const Field<T> field) 
    {
        super(field);
    }

    /**
     * Create a {@code Field_Matrix<T>} with the supplied row and column dimensions.
     *
     * @param field Field to which the elements belong.
     * @param row_dimension Number of rows in the matrix.
     * @param column_dimension Number of columns in the matrix.
     * @ if row or column dimension is not positive.
     */
    Array2DRowField_Matrix(const Field<T>& field, const int& row_dimension, const int& column_dimension)
    {
        super(field, row_dimension, column_dimension);
        data = Math_Arrays::build_array(field, row_dimension, column_dimension);
    }

    /**
     * Create a {@code Field_Matrix<T>} using the input array as the underlying
     * data array.
     * <p>The input array is copied, not referenced. This constructor has
     * the same effect as calling {@link #Array2DRowField_Matrix(Field_Element[][], bool)}
     * with the second argument set to {@code true}.</p>
     *
     * @param d Data for the matrix.
     * @ if {@code d} is not rectangular.
     * @Null_Argument_Exception if {@code d} is {@code NULL}.
     * @ if there are not at least one row and one column.
     * @see #Array2DRowField_Matrix(Field_Element[][], bool)
     */
    Array2DRowField_Matrix(const std::vector<std::vector<T>>& d)
    {
        this(extract_field(d), d);
    }

    /**
     * Create a {@code Field_Matrix<T>} using the input array as the underlying
     * data array.
     * <p>The input array is copied, not referenced. This constructor has
     * the same effect as calling {@link #Array2DRowField_Matrix(Field_Element[][], bool)}
     * with the second argument set to {@code true}.</p>
     *
     * @param field Field to which the elements belong.
     * @param d Data for the matrix.
     * @ if {@code d} is not rectangular.
     * @Null_Argument_Exception if {@code d} is {@code NULL}.
     * @ if there are not at least one row and one column.
     * @see #Array2DRowField_Matrix(Field_Element[][], bool)
     */
    Array2DRowField_Matrix(const Field<T>& field, const std::vector<std::vector<T>>& d)
    {
        super(field);
        copy_in(d);
    }

    /**
     * Create a {@code Field_Matrix<T>} using the input array as the underlying
     * data array.
     * <p>If an array is built specially in order to be embedded in a
     * {@code Field_Matrix<T>} and not used directly, the {@code copy_array} may be
     * set to {@code false}. This will prevent the copying and improve
     * performance as no array will be built and no data will be copied.</p>
     *
     * @param d Data for the matrix.
     * @param copy_array Whether to copy or reference the input array.
     * @ if {@code d} is not rectangular.
     * @ if there are not at least one row and one column.
     * @Null_Argument_Exception if {@code d} is {@code NULL}.
     * @see #Array2DRowField_Matrix(Field_Element[][])
     */
    Array2DRowField_Matrix(const std::vector<std::vector<T>>& d, const bool copy_array)
    {
        this(extract_field(d), d, copy_array);
    }

    /**
     * Create a {@code Field_Matrix<T>} using the input array as the underlying
     * data array.
     * <p>If an array is built specially in order to be embedded in a
     * {@code Field_Matrix<T>} and not used directly, the {@code copy_array} may be
     * set to {@code false}. This will prevent the copying and improve
     * performance as no array will be built and no data will be copied.</p>
     *
     * @param field Field to which the elements belong.
     * @param d Data for the matrix.
     * @param copy_array Whether to copy or reference the input array.
     * @ if {@code d} is not rectangular.
     * @ if there are not at least one row and one column.
     * @Null_Argument_Exception if {@code d} is {@code NULL}.
     * @see #Array2DRowField_Matrix(Field_Element[][])
     */
    Array2DRowField_Matrix(const Field<T>& field, const std::vector<std::vector<T>>& d, const bool copy_array) // NOPMD - array copy is taken care of by parameter
    {
        super(field);
        if (copy_array) 
        {
            copy_in(d);
        }
        else 
        {
            //Math_Utils::check_not_null(d);
            const int& n_rows = d.size();
            if (n_rows == 0) 
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
            for (const int& r = 1; r < n_rows; r++) 
            {
                if (d[r].size() != n_cols) 
                {
                    throw std::exception("not implemented");
                    // throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, n_cols, d[r].size());
                }
            }
            data = d;
        }
    }

    /**
     * Create a (column) {@code Field_Matrix<T>} using {@code v} as the
     * data for the unique column of the created matrix.
     * The input array is copied.
     *
     * @param v Column vector holding data for matrix.
     * @ if v is empty
     */
    Array2DRowField_Matrix(const std::vector<T>& v)  
    {
        this(extract_field(v), v);
    }

    /**
     * Create a (column) {@code Field_Matrix<T>} using {@code v} as the
     * data for the unique column of the created matrix.
     * The input array is copied.
     *
     * @param field Field to which the elements belong.
     * @param v Column vector holding data for matrix.
     */
    Array2DRowField_Matrix(const Field<T>& field, const std::vector<T>& v) 
    {
        super(field);
        const int& n_rows = v.size();
        data = Math_Arrays::build_array(get_field(), n_rows, 1);
        for (int row{}; row < n_rows; row++) 
        {
            data[row][0] = v[row];
        }
    }

    /** {@inherit_doc} */
    //override
    Field_Matrix<T> create_matrix(const int& row_dimension, const int& column_dimension)
    {
        return Array2DRowField_Matrix<T>(get_field(), row_dimension, column_dimension);
    }

    /** {@inherit_doc} */
    //override
    Field_Matrix<T> copy() 
    {
        return Array2DRowField_Matrix<T>(get_field(), copy_out(), false);
    }

    /**
     * Add {@code m} to this matrix.
     *
     * @param m Matrix to be added.
     * @return {@code this} + m.
     * @ if {@code m} is not the same
     * size as this matrix.
     */
    Array2DRowField_Matrix<T> add(const Array2DRowField_Matrix<T>& m)
         
        {
        // safety check
        check_addition_compatible(m);

        const int row_count    = get_row_dimension();
        const int column_count = get_column_dimension();
        const std::vector<std::vector<T>> out_data = Math_Arrays::build_array(get_field(), row_count, column_count);
        for (int row{}; row < row_count; row++) 
        {
            const std::vector<T> data_row    = data[row];
            const std::vector<T> m_row       = m.data[row];
            const std::vector<T> out_data_row = out_data[row];
            for (int col{};  col < column_count; col++) 
            {
                out_data_row[col] = data_row[col].add(m_row[col]);
            }
        }

        return Array2DRowField_Matrix<T>(get_field(), out_data, false);
    }

    /**
     * Subtract {@code m} from this matrix.
     *
     * @param m Matrix to be subtracted.
     * @return {@code this} + m.
     * @ if {@code m} is not the same
     * size as this matrix.
     */
    Array2DRowField_Matrix<T> subtract(const Array2DRowField_Matrix<T>& m)
         
        {
        // safety check
        check_subtraction_compatible(m);

        const int row_count    = get_row_dimension();
        const int column_count = get_column_dimension();
        const std::vector<std::vector<T>> out_data = Math_Arrays::build_array(get_field(), row_count, column_count);
        for (int row{}; row < row_count; row++)
        {
            const std::vector<T> data_row    = data[row];
            const std::vector<T> m_row       = m.data[row];
            const std::vector<T> out_data_row = out_data[row];
            for (int col{};  col < column_count; col++) 
            {
                out_data_row[col] = data_row[col].subtract(m_row[col]);
            }
        }

        return Array2DRowField_Matrix<T>(get_field(), out_data, false);

    }

    /**
     * Postmultiplying this matrix by {@code m}.
     *
     * @param m Matrix to postmultiply by.
     * @return {@code this} * m.
     * @ if the number of columns of this
     * matrix is not equal to the number of rows of {@code m}.
     */
    Array2DRowField_Matrix<T> multiply(const Array2DRowField_Matrix<T>& m)
    {
        // safety check
        check_multiplication_compatible(m);

        const int& n_rows = this.get_row_dimension();
        const int& n_cols = m.get_column_dimension();
        const int& n_sum = this.get_column_dimension();
        const std::vector<std::vector<T>> out_data = Math_Arrays::build_array(get_field(), n_rows, n_cols);
        for (int row{}; row < n_rows; row++) 
        {
            const std::vector<T> data_row    = data[row];
            const std::vector<T> out_data_row = out_data[row];
            for (int col{};  col < n_cols; col++) 
            {
                T sum = get_field().get_zero();
                for (int i{}; i < n_sum; i++) 
                {
                    sum = sum.add(data_row[i].multiply(m.data[i][col]));
                }
                out_data_row[col] = sum;
            }
        }

        return Array2DRowField_Matrix<T>(get_field(), out_data, false);

    }

    /**
     * Returns the result of postmultiplying {@code this} by {@code m^T}.
     * @param m matrix to first transpose and second postmultiply by
     * @return {@code this * m^T}
     * @ if
     * {@code column_dimension(this) != column_dimension(m)}
     * @since 1.3
     */
    Field_Matrix<T> multiply_transposed(const Array2DRowField_Matrix<T>& m)
    {
        Matrix_Utils::check_same_column_dimension(this, m);

        const int& n_rows = this.get_row_dimension();
        const int& n_cols = m.get_row_dimension();
        const int& n_sum  = this.get_column_dimension();

        const Field_Matrix<T> out   = Matrix_Utils::create_field_matrix(get_field(), n_rows, n_cols);
        const std::vector<std::vector<T>>          m_data = m.data;

        // Multiply.
        for (int col{};  col < n_cols; col++) 
        {
            for (int row{}; row < n_rows; row++) 
            {
                const std::vector<T> data_row = data[row];
                const std::vector<T> m_row    = m_data[col];
                T sum = get_field().get_zero();
                for (int i{}; i < n_sum; i++) 
                {
                    sum = sum.add(data_row[i].multiply(m_row[i]));
                }
                out.set_entry(row, col, sum);
            }
        }

        return out;

    }

    /** {@inherit_doc} */
    //override
    Field_Matrix<T> multiply_transposed(const Field_Matrix<T>& m) 
    {
        if (dynamic_cast<const Array2DRowField_Matrix*>(*m) != nullptr)
        {
            return multiply_transposed((Array2DRowField_Matrix<T>) m);
        }
        Matrix_Utils::check_same_column_dimension(this, m);

        const int n_rows = this.get_row_dimension();
        const int n_cols = m.get_row_dimension();
        const int n_sum  = this.get_column_dimension();

        auto out = Matrix_Utils::create_field_matrix(get_field(), n_rows, n_cols);

        // Multiply.
        for (int col{};  col < n_cols; col++) 
        {
            for (int row{}; row < n_rows; row++) 
            {
                const std::vector<T> data_row = data[row];
                T sum = get_field().get_zero();
                for (int i{}; i < n_sum; i++) 
                {
                    sum = sum.add(data_row[i].multiply(m.get_entry(col, i)));
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
    Field_Matrix<T> transpose_multiply(const Array2DRowField_Matrix<T>& m)
    {
        Matrix_Utils::check_same_row_dimension(this, m);

        const int& n_rows = this.get_column_dimension();
        const int& n_cols = m.get_column_dimension();
        const int& n_sum  = this.get_row_dimension();

        const Field_Matrix<T> out   = Matrix_Utils::create_field_matrix(get_field(), n_rows, n_cols);
        const std::vector<std::vector<T>>          m_data = m.data;

        // Multiply.
        for (int k{}; k < n_sum; k++) 
        {
            const std::vector<T> data_k = data[k];
            const std::vector<T> mK    = m_data[k];
            for (int row{}; row < n_rows; row++) 
            {
                const T data_i_row = data_k[row];
                for (int col{};  col < n_cols; col++) 
                {
                    out.add_to_entry(row, col, data_i_row.multiply(mK[col]));
                }
            }
        }

        return out;
    }

    /** {@inherit_doc} */
    //override
    Field_Matrix<T> transpose_multiply(const Field_Matrix<T>& m) 
    {
        if (dynamic_cast<const Array2DRowField_Matrix*>(*m) != nullptr)
        {
            return transpose_multiply((Array2DRowField_Matrix<T>) m);
        }
        Matrix_Utils::check_same_row_dimension(this, m);

        const int& n_rows = this.get_column_dimension();
        const int& n_cols = m.get_column_dimension();
        const int& n_sum  = this.get_row_dimension();

        const Field_Matrix<T> out = Matrix_Utils::create_field_matrix(get_field(), n_rows, n_cols);

        // Multiply.
        for (int k{}; k < n_sum; k++) 
        {
            const std::vector<T> data_k = data[k];
            for (int row{}; row < n_rows; row++) 
            {
                const T data_i_row = data_k[row];
                for (int col{};  col < n_cols; col++) 
                {
                    out.add_to_entry(row, col, data_i_row.multiply(m.get_entry(k, col)));
                }
            }
        }

        return out;
    }

    /** {@inherit_doc} */
    //override
    std::vector<std::vector<T>> get_data() 
    {
        return copy_out();
    }

    /**
     * Get a reference to the underlying data array.
     * This methods returns internal data, <strong>not</strong> fresh copy of it.
     *
     * @return the 2-dimensional array of entries.
     */
    std::vector<std::vector<T>> get_data_ref() 
    {
        return data; // NOPMD - returning an internal array is intentional and documented here
    }

    /** {@inherit_doc} */
    //override
    void set_sub_matrix(const std::vector<std::vector<T>>& sub_matrix, const int& row, const int& column)
    {
        if (data == NULL) 
        {
            if (row > 0) 
            {
                throw std::exception("not implemented");
                //throw Math_Illegal_State_Exception(hipparchus::exception::Localized_Core_Formats_Type::FIRST_ROWS_NOT_INITIALIZED_YET, row);
            }
            if (column > 0) 
            {
                throw std::exception("not implemented");
                //throw Math_Illegal_State_Exception(hipparchus::exception::Localized_Core_Formats_Type::FIRST_COLUMNS_NOT_INITIALIZED_YET, column);
            }
            const int& n_rows = sub_matrix.size();
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
            data = Math_Arrays::build_array(get_field(), sub_matrix.size(), n_cols);
            for (int i{}; i < data.size(); ++i) 
            {
                if (sub_matrix[i].size() != n_cols) 
                {
                    throw std::exception("not implemented");
                    //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, n_cols, sub_matrix[i].size());
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
    T get_entry(const int& row, const int& column)
    {
        check_row_index(row);
        check_column_index(column);

        return data[row][column];
    }

    /** {@inherit_doc} */
    //override
    void set_entry(const int& row, const int& column, const T& value)
    {
        check_row_index(row);
        check_column_index(column);

        data[row][column] = value;
    }

    /** {@inherit_doc} */
    //override
    void add_to_entry(const int& row, const int& column, const T& increment)
    {
        check_row_index(row);
        check_column_index(column);

        data[row][column] = data[row][column].add(increment);
    }

    /** {@inherit_doc} */
    //override
    void multiply_entry(const int& row, const int& column, const T& factor)
    {
        check_row_index(row);
        check_column_index(column);

        data[row][column] = data[row][column].multiply(factor);
    }

    /** {@inherit_doc} */
    //override
    int get_row_dimension() 
    {
        return (data == NULL) ? 0 : data.size();
    }

    /** {@inherit_doc} */
    //override
    int get_column_dimension() 
    {
        return ((data == NULL) || (data[0] == NULL)) ? 0 : data[0].size();
    }

    /** {@inherit_doc} */
    //override
    std::vector<T> operate(const std::vector<T>& v)  
    {
        const int n_rows = this.get_row_dimension();
        const int n_cols = this.get_column_dimension();
        if (v.size() != n_cols) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, v.size(), n_cols);
        }
        const std::vector<T> out = Math_Arrays::build_array(get_field(), n_rows);
        for (int row{}; row < n_rows; row++) 
        {
            const std::vector<T> data_row = data[row];
            T sum = get_field().get_zero();
            for (int i{}; i < n_cols; i++) 
            {
                sum = sum.add(data_row[i].multiply(v[i]));
            }
            out[row] = sum;
        }
        return out;
    }

    /** {@inherit_doc} */
    //override
    std::vector<T> pre_multiply(const std::vector<T>& v)  
    {
        const int n_rows = get_row_dimension();
        const int n_cols = get_column_dimension();
        if (v.size() != n_rows) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, v.size(), n_rows);
        }

        const std::vector<T> out = Math_Arrays::build_array(get_field(), n_cols);
        for (int col{};  col < n_cols; ++col) 
        {
            T sum = get_field().get_zero();
            for (int i{}; i < n_rows; ++i) 
            {
                sum = sum.add(data[i][col].multiply(v[i]));
            }
            out[col] = sum;
        }

        return out;
    }

    /** {@inherit_doc} */
    //override
    Field_Matrix<T> get_sub_matrix(const int& start_row, const int& end_row, const int& start_column, const int& end_column)
    {
        Matrix_Utils::check_sub_matrix_index(this, start_row, end_row, start_column, end_column);
        const int row_count = end_row - start_row + 1;
        const int column_count = end_column - start_column + 1;
        const std::vector<std::vector<T>> out_data = Math_Arrays::build_array(get_field(), row_count, column_count);
        for (int i{}; i < row_count; ++i) 
        {
            System.arraycopy(data[start_row + i], start_column, out_data[i], 0, column_count);
        }

        Array2DRowField_Matrix<T> sub_matrix = Array2DRowField_Matrix<>(get_field());
        sub_matrix.data = out_data;
        return sub_matrix;
    }

    /** {@inherit_doc} */
    //override
    T walk_in_row_order(const Field_Matrix_Changing_Visitor<T>& visitor) 
    {
        const int rows = get_row_dimension();
        const int columns = get_column_dimension();
        visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
        for (int i{}; i < rows; ++i) 
        {
            const std::vector<T> row_i = data[i];
            for (int j{}; j < columns; ++j) 
            {
                row_i[j] = visitor.visit(i, j, row_i[j]);
            }
        }
        return visitor.end();
    }

    /** {@inherit_doc} */
    //override
    T walk_in_row_order(const Field_Matrix_Preserving_Visitor<T>& visitor) 
    {
        const int rows    = get_row_dimension();
        const int columns = get_column_dimension();
        visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
        for (int i{}; i < rows; ++i) 
        {
            const std::vector<T> row_i = data[i];
            for (int j{}; j < columns; ++j) 
            {
                visitor.visit(i, j, row_i[j]);
            }
        }
        return visitor.end();
    }

    /** {@inherit_doc} */
    //override
    T walk_in_row_order(const Field_Matrix_Changing_Visitor<T>& visitor, const int& start_row, const int& end_row, const int& start_column, const int& end_column)
    {
        check_sub_matrix_index(start_row, end_row, start_column, end_column);
        visitor.start(get_row_dimension(), get_column_dimension(), start_row, end_row, start_column, end_column);
        for (int i = start_row; i <= end_row; ++i) 
        {
            const std::vector<T> row_i = data[i];
            for (int j = start_column; j <= end_column; ++j) 
            {
                row_i[j] = visitor.visit(i, j, row_i[j]);
            }
        }
        return visitor.end();
    }

    /** {@inherit_doc} */
    //override
    T walk_in_row_order(const Field_Matrix_Preserving_Visitor<T>& visitor, const int& start_row, const int& end_row, const int& start_column, const int& end_column)
    {
        check_sub_matrix_index(start_row, end_row, start_column, end_column);
        visitor.start(get_row_dimension(), get_column_dimension(), start_row, end_row, start_column, end_column);
        for (int i = start_row; i <= end_row; ++i) 
        {
            const std::vector<T> row_i = data[i];
            for (int j = start_column; j <= end_column; ++j) 
            {
                visitor.visit(i, j, row_i[j]);
            }
        }
        return visitor.end();
    }

    /** {@inherit_doc} */
    //override
    T walk_in_column_order(const Field_Matrix_Changing_Visitor<T>& visitor) 
    {
        const int rows    = get_row_dimension();
        const int columns = get_column_dimension();
        visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
        for (int j{}; j < columns; ++j) 
        {
            for (int i{}; i < rows; ++i) 
            {
                const std::vector<T> row_i = data[i];
                row_i[j] = visitor.visit(i, j, row_i[j]);
            }
        }
        return visitor.end();
    }

    /** {@inherit_doc} */
    //override
    T walk_in_column_order(const Field_Matrix_Preserving_Visitor<T>& visitor) 
    {
        const int rows    = get_row_dimension();
        const int columns = get_column_dimension();
        visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
        for (int j{}; j < columns; ++j) 
        {
            for (int i{}; i < rows; ++i) 
            {
                visitor.visit(i, j, data[i][j]);
            }
        }
        return visitor.end();
    }

    /** {@inherit_doc} */
    //override
    T walk_in_column_order(const Field_Matrix_Changing_Visitor<T>& visitor, const int& start_row, const int& end_row, const int& start_column, const int& end_column)
    {
        check_sub_matrix_index(start_row, end_row, start_column, end_column);
        visitor.start(get_row_dimension(), get_column_dimension(), start_row, end_row, start_column, end_column);
        for (int j = start_column; j <= end_column; ++j) 
        {
            for (int i = start_row; i <= end_row; ++i) 
            {
                const std::vector<T> row_i = data[i];
                row_i[j] = visitor.visit(i, j, row_i[j]);
            }
        }
        return visitor.end();
    }

    /** {@inherit_doc} */
    //override
    T walk_in_column_order(const Field_Matrix_Preserving_Visitor<T>& visitor, const int& start_row, const int& end_row, const int& start_column, const int& end_column)
    {
        check_sub_matrix_index(start_row, end_row, start_column, end_column);
        visitor.start(get_row_dimension(), get_column_dimension(), start_row, end_row, start_column, end_column);
        for (int j = start_column; j <= end_column; ++j) 
        {
            for (int i = start_row; i <= end_row; ++i) 
            {
                visitor.visit(i, j, data[i][j]);
            }
        }
        return visitor.end();
    }



    /** {@inherit_doc} */
    //override
    std::vector<T> get_row(const int& row)  
    {
        Matrix_Utils::check_row_index(this, row);
        const int& n_cols = get_column_dimension();
        const std::vector<T> out = Math_Arrays::build_array(get_field(), n_cols);
        System.arraycopy(data[row], 0, out, 0, n_cols);
        return out;
    }

    /** {@inherit_doc} */
    //override
    void set_row(const int& row, const std::vector<T>& arr)
    {
        Matrix_Utils::check_row_index(this, row);
        const int& n_cols = get_column_dimension();
        if (arr.size() != n_cols) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH_2x2, 1, arr.size(), 1, n_cols);
        }
        System.arraycopy(arr, 0, data[row], 0, n_cols);
    }

};