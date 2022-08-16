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

//import java.util.Array_list;

//import org.hipparchus.Field;
//import org.hipparchus.Field_Element;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Null_Argument_Exception;
//import org.hipparchus.util.Math_Arrays;
#include <type_traits>
#include <vector>
#include "MatrixUtils.h"
#include "../Field.h"
#include "../FieldElement.h"

/**
 * Basic implementation of {@link Field_Matrix} methods regardless of the underlying storage.
 * <p>All the methods implemented here use {@link #get_entry(int, int)} to access
 * matrix elements. Derived class can provide faster implementations. </p>
 *
 * @param <T> Type of the field elements.
 *
 */
template<typename T, typename std::enable_if<std::is_base_of<Field_Element<T>, T>::value>::type* = nullptr>
class Abstract_Field_Matrix : public Field_Matrix<T> 
{
private:

    /** Field to which the elements belong. */
    const Field<T> my_field;

protected:
    /**
     * Constructor for use with Serializable
     */
    Abstract_Field_Matrix() 
    {
        field = NULL;
    }

    /**
     * Creates a matrix with no data
     * @param field field to which the elements belong
     */
    Abstract_Field_Matrix(const Field<T> field) : my_field{ field } {};

    /**
     * Create a Field_Matrix<T> with the supplied row and column dimensions.
     *
     * @param field Field to which the elements belong.
     * @param row_dimension Number of rows in the matrix.
     * @param column_dimension Number of columns in the matrix.
     * @ if row or column dimension is not
     * positive.
     */
    Abstract_Field_Matrix(const Field<T>& field, const int& row_dimension, const int& column_dimension)
    {
        if (row_dimension <= 0) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSION, row_dimension);
        }
        if (column_dimension <= 0) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSION, column_dimension);
        }
        my_field = field;
    }

    /**
     * Get the elements type from an array.
     *
     * @param <T> Type of the field elements.
     * @param d Data array.
     * @return the field to which the array elements belong.
     * @Null_Argument_Exception if the array is {@code NULL}.
     * @ if the array is empty.
     */
    static Field<T> extract_field(const std::vector<std::vector<T>> d)
    {
        if (d == NULL) 
        {
            throw std::exception("not implemented");
            //throw Null_Argument_Exception();
        }
        if (d.size() == 0) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::AT_LEAST_ONE_ROW);
        }
        if (d[0].size() == 0) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::AT_LEAST_ONE_COLUMN);
        }
        return d[0][0].get_field();
    }

    /**
     * Get the elements type from an array.
     *
     * @param <T> Type of the field elements.
     * @param d Data array.
     * @return the field to which the array elements belong.
     * @ if array is empty.
     */
    static Field<T> extract_field(const std::vector<T>& d)
    {
        if (d.size() == 0) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::AT_LEAST_ONE_ROW);
        }
        return d[0].get_field();
    }

public:

    /** {@inherit_doc} */
    //override
    public Field<T> get_field() const 
    {
        return my_field;
    }

    /** {@inherit_doc} */
    //override
    virtual Field_Matrix<T> create_matrix(const int& row_dimension, const int& column_dimension);

    /** {@inherit_doc} */
    //override
    virtual Field_Matrix<T> copy();

    /** {@inherit_doc} */
    //override
    Field_Matrix<T> add(Field_Matrix<T> m)     
    {
        // safety check
        check_addition_compatible(m);

        const int row_count    = get_row_dimension();
        const int column_count = get_column_dimension();
        Field_Matrix<T> out = create_matrix(row_count, column_count);
        for (int row{}; row < row_count; ++row)
        {
            for (int col{}; col < column_count; ++col)
            {
                out.set_entry(row, col, get_entry(row, col).add(m.get_entry(row, col)));
            }
        }

        return out;
    }

    /** {@inherit_doc} */
    //override
    Field_Matrix<T> subtract(const Field_Matrix<T>& m)
    {
        // safety check
        check_subtraction_compatible(m);

        const int row_count    = get_row_dimension();
        const int column_count = get_column_dimension();
        Field_Matrix<T> out = create_matrix(row_count, column_count);
        for (int row{}; row < row_count; ++row) 
        {
            for (int col{};  col < column_count; ++col) 
            {
                out.set_entry(row, col, get_entry(row, col).subtract(m.get_entry(row, col)));
            }
        }

        return out;
    }

    /** {@inherit_doc} */
    //override
    Field_Matrix<T> scalar_add(const T d) 
    {
        const int row_count    = get_row_dimension();
        const int column_count = get_column_dimension();
        Field_Matrix<T> out = create_matrix(row_count, column_count);
        for (int row{}; row < row_count; ++row) 
        {
            for (int col{};  col < column_count; ++col) 
            {
                out.set_entry(row, col, get_entry(row, col).add(d));
            }
        }

        return out;
    }

    /** {@inherit_doc} */
    //override
    Field_Matrix<T> scalar_multiply(const T d) 
    {
        const int row_count    = get_row_dimension();
        const int column_count = get_column_dimension();
        const Field_Matrix<T> out = create_matrix(row_count, column_count);
        for (int row{}; row < row_count; ++row) 
        {
            for (int col{};  col < column_count; ++col) 
            {
                out.set_entry(row, col, get_entry(row, col).multiply(d));
            }
        }

        return out;
    }

    /** {@inherit_doc} */
    //override
    Field_Matrix<T> multiply(const Field_Matrix<T>& m)     
    {
        // safety check
        check_multiplication_compatible(m);

        const int n_rows = get_row_dimension();
        const int n_cols = m.get_column_dimension();
        const int& n_sum  = get_column_dimension();
        const Field_Matrix<T> out = create_matrix(n_rows, n_cols);
        for (int row{}; row < n_rows; ++row)
        {
            for (int col{}; col < n_cols; ++col)
            {
                T sum = field.get_zero();
                for (int i{}; i < n_sum; ++i) 
                {
                    sum = sum.add(get_entry(row, i).multiply(m.get_entry(i, col)));
                }
                out.set_entry(row, col, sum);
            }
        }

        return out;
    }

    /** {@inherit_doc} */
    //override
    Field_Matrix<T> pre_multiply(const Field_Matrix<T>& m)
    {
        return m.multiply(this);
    }

    /** {@inherit_doc} */
    //override
    public Field_Matrix<T> power(const int p)  
    {
        if (p < 0) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL, p, 0);
        }

        if (!is_square()) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NON_SQUARE_MATRIX, get_row_dimension(), get_column_dimension());
        }

        if (p == 0) 
        {
            return Matrix_Utils::create_field_identity_matrix(this.get_field(), this.get_row_dimension());
        }

        if (p == 1) 
        {
            return this.copy();
        }

        const int power = p - 1;

        /*
         * Only log_2(p) operations is used by doing as follows:
         * 5^214 = 5^128 * 5^64 * 5^16 * 5^4 * 5^2
         *
         * In general, the same approach is used for A^p.
         */

        const char[] binary_representation = Integer.to_binary_string(power).to_char_array();
        const Array_list<Integer> non_zero_positions = Array_list<>();

        for (int i{}; i < binary_representation.size(); ++i) 
        {
            if (binary_representation[i] == '1') 
            {
                const int pos = binary_representation.size() - i - 1;
                non_zero_positions.add(pos);
            }
        }

        Array_list<Field_Matrix<T>> results = Array_list<>(binary_representation.size());

        results.add(0, this.copy());

        for (int i{ 1 }; i < binary_representation.size(); ++i) 
        {
            const Field_Matrix<T> s = results.get(i - 1);
            const Field_Matrix<T> r = s.multiply(s);
            results.add(i, r);
        }

        Field_Matrix<T> result = this.copy();

        for (Integer i : non_zero_positions) 
        {
            result = result.multiply(results.get(i));
        }

        return result;
    }

    /** {@inherit_doc} */
    //override
    public std::vector<std::vector<T>> get_data() 
    {
        const std::vector<std::vector<T>> data = Math_Arrays::build_array(field, get_row_dimension(), get_column_dimension());

        for (int i{}; i < data.size(); ++i) 
        {
            const std::vector<T> data_i = data[i];
            for (int j{}; j < data_i.size(); ++j) 
            {
                data_i[j] = get_entry(i, j);
            }
        }

        return data;
    }

    /** {@inherit_doc} */
    //override
    public Field_Matrix<T> get_sub_matrix(const int& start_row, const int& end_row, const int& start_column, const int& end_column)
         
        {
        check_sub_matrix_index(start_row, end_row, start_column, end_column);

        const Field_Matrix<T> sub_matrix =
            create_matrix(end_row - start_row + 1, end_column - start_column + 1);
        for (int i = start_row; i <= end_row; ++i) 
        {
            for (int j = start_column; j <= end_column; ++j) 
            {
                sub_matrix.set_entry(i - start_row, j - start_column, get_entry(i, j));
            }
        }

        return sub_matrix;

    }

    /** {@inherit_doc} */
    //override
    public Field_Matrix<T> get_sub_matrix(const std::vector<int> selected_rows, const std::vector<int> selected_columns)
    , Null_Argument_Exception 
    {

        // safety checks
        check_sub_matrix_index(selected_rows, selected_columns);

        // copy entries
        const Field_Matrix<T> sub_matrix =
            create_matrix(selected_rows.size(), selected_columns.size());
        sub_matrix.walk_in_optimized_order(new DefaultField_Matrix_Changing_Visitor<T>(field.get_zero()) 
        {

            /** {@inherit_doc} */
            //override
            public T visit(const int& row, const int& column, const T value) 
            {
                return get_entry(selected_rows[row], selected_columns[column]);
            }

        });

        return sub_matrix;

    }

    /** {@inherit_doc} */
    //override
    public void copy_sub_matrix(const int& start_row, const int& end_row, const int& start_column, const int& end_column, const std::vector<std::vector<T>> destination)
     
    {
        // safety checks
        check_sub_matrix_index(start_row, end_row, start_column, end_column);
        const int rows_count    = end_row + 1 - start_row;
        const int columns_count = end_column + 1 - start_column;
        if ((destination.size() < rows_count) || (destination[0].size() < columns_count)) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH_2x2, destination.size(), destination[0].size(), rows_count, columns_count);
        }

        // copy entries
        walk_in_optimized_order(DefaultField_Matrix_Preserving_Visitor<T>(field.get_zero()) 
        {

            /** Initial row index. */
            private int my_start_row;

            /** Initial column index. */
            private int my_start_column;

            /** {@inherit_doc} */
            //override
            public void start(const int rows, const int columns, const int& start_row, const int& end_row, const int& start_column, const int& end_column) 
            {
                this.my_start_row    = start_row;
                this.my_start_column = start_column;
            }

            /** {@inherit_doc} */
            //override
            public void visit(const int& row, const int& column, const T& value) 
            {
                destination[row - start_row][column - start_column] = value;
            }

        }, start_row, end_row, start_column, end_column);

    }

    /** {@inherit_doc} */
    //override
    public void copy_sub_matrix(std::vector<int> selected_rows, std::vector<int> selected_columns, std::vector<std::vector<T>> destination)
    {
        // safety checks
        check_sub_matrix_index(selected_rows, selected_columns);
        if ((destination.size() < selected_rows.size()) || (destination[0].size() < selected_columns.size())) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH_2x2, destination.size(), destination[0].size(), selected_rows.size(), selected_columns.size());
        }

        // copy entries
        for (int i{}; i < selected_rows.size(); i++) 
        {
            const std::vector<T> destination_i = destination[i];
            for (int j{}; j < selected_columns.size(); j++) 
            {
                destination_i[j] = get_entry(selected_rows[i], selected_columns[j]);
            }
        }

    }

    /** {@inherit_doc} */
    //override
    public void set_sub_matrix(const std::vector<std::vector<T>> sub_matrix, const int& row, const int column)
        , Null_Argument_Exception 
        {
        if (sub_matrix == NULL) 
        {
            throw std::exception("not implemented");
            //throw Null_Argument_Exception();
        }
        const int n_rows = sub_matrix.size();
        if (n_rows == 0) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::AT_LEAST_ONE_ROW);
        }

        const int n_cols = sub_matrix[0].size();
        if (n_cols == 0) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::AT_LEAST_ONE_COLUMN);
        }

        for (const int r{1}; r < n_rows; ++r)
        {
            if (sub_matrix[r].size() != n_cols) 
            {
                throw std::exception("not implemented");
                //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, n_cols, sub_matrix[r].size());
            }
        }

        check_row_index(row);
        check_column_index(column);
        check_row_index(n_rows + row - 1);
        check_column_index(n_cols + column - 1);

        for (int i{}; i < n_rows; ++i) 
        {
            for (int j{}; j < n_cols; ++j) 
            {
                set_entry(row + i, column + j, sub_matrix[i][j]);
            }
        }
    }

    /** {@inherit_doc} */
    //override
    public Field_Matrix<T> get_row_matrix(const int row)  
    {
        check_row_index(row);
        const int n_cols = get_column_dimension();
        const Field_Matrix<T> out = create_matrix(1, n_cols);
        for (int i{}; i < n_cols; ++i) 
        {
            out.set_entry(0, i, get_entry(row, i));
        }

        return out;

    }

    /** {@inherit_doc} */
    //override
    public void set_row_matrix(const int& row, const Field_Matrix<T> matrix)
         
        {
        check_row_index(row);
        const int n_cols = get_column_dimension();
        if ((matrix.get_row_dimension() != 1) || (matrix.get_column_dimension() != n_cols)) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH_2x2, matrix.get_row_dimension(), matrix.get_column_dimension(), 1, n_cols);
        }
        for (int i{}; i < n_cols; ++i) 
        {
            set_entry(row, i, matrix.get_entry(0, i));
        }

    }

    /** {@inherit_doc} */
    //override
    public Field_Matrix<T> get_column_matrix(const int column)
    {
        check_column_index(column);
        const int n_rows = get_row_dimension();
        Field_Matrix<T> out = create_matrix(n_rows, 1);
        for (int i{}; i < n_rows; ++i) 
        {
            out.set_entry(i, 0, get_entry(i, column));
        }

        return out;

    }

    /** {@inherit_doc} */
    //override
    public void set_column_matrix(const int& column, const Field_Matrix<T> matrix)
         
        {
        check_column_index(column);
        const int n_rows = get_row_dimension();
        if ((matrix.get_row_dimension() != n_rows) || (matrix.get_column_dimension() != 1)) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH_2x2, matrix.get_row_dimension(), matrix.get_column_dimension(), n_rows, 1);
        }
        for (int i{}; i < n_rows; ++i) 
        {
            set_entry(i, column, matrix.get_entry(i, 0));
        }

    }

    /** {@inherit_doc} */
    //override
    public Field_Vector<T> get_row_vector(const int& row)
    {
        return ArrayField_Vector<T>(field, get_row(row), false);
    }

    /** {@inherit_doc} */
    //override
    public void set_row_vector(const int& row, const Field_Vector<T>& vector)
         
        {
        check_row_index(row);
        const int n_cols = get_column_dimension();
        if (vector.get_dimension() != n_cols) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH_2x2, 1, vector.get_dimension(), 1, n_cols);
        }
        for (int i{}; i < n_cols; ++i) 
        {
            set_entry(row, i, vector.get_entry(i));
        }

    }

    /** {@inherit_doc} */
    //override
    public Field_Vector<T> get_column_vector(const int column)
         
        {
        return ArrayField_Vector<T>(field, get_column(column), false);
    }

    /** {@inherit_doc} */
    //override
    public void set_column_vector(const int& column, const Field_Vector<T>& vector)
    {
        check_column_index(column);
        const int n_rows = get_row_dimension();
        if (vector.get_dimension() != n_rows) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH_2x2, vector.get_dimension(), 1, n_rows, 1);
        }
        for (int i{}; i < n_rows; ++i) 
        {
            set_entry(i, column, vector.get_entry(i));
        }

    }

    /** {@inherit_doc} */
    //override
    public std::vector<T> get_row(const int& row)  
    {
        check_row_index(row);
        const int n_cols = get_column_dimension();
        std::vector<T> out = Math_Arrays::build_array(field, n_cols);
        for (int i{}; i < n_cols; ++i) 
        {
            out[i] = get_entry(row, i);
        }

        return out;

    }

    /** {@inherit_doc} */
    //override
    public void set_row(const int& row, const std::vector<T>& arr)
    {
        check_row_index(row);
        const int n_cols = get_column_dimension();
        if (array.size() != n_cols) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH_2x2, 1, array.size(), 1, n_cols);
        }
        for (int i{}; i < n_cols; ++i) 
        {
            set_entry(row, i, arr[i]);
        }
    }

    /** {@inherit_doc} */
    //override
    public std::vector<T> get_column(const int& column)  
    {
        check_column_index(column);
        const int n_rows = get_row_dimension();
        std::vector<T> out = Math_Arrays::build_array(field, n_rows);
        for (int i{}; i < n_rows; ++i) 
        {
            out[i] = get_entry(i, column);
        }

        return out;
    }

    /** {@inherit_doc} */
    //override
    public void set_column(const int& column, const std::vector<T>& arr)
         
        {
        check_column_index(column);
        const int n_rows = get_row_dimension();
        if (arr.size() != n_rows) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH_2x2, array.size(), 1, n_rows, 1);
        }
        for (int i{}; i < n_rows; ++i) 
        {
            set_entry(i, column, arr[i]);
        }
    }

    /** {@inherit_doc} */
    //override
    public virtual T get_entry(const int& row, const int& column);

    /** {@inherit_doc} */
    //override
    public virtual void set_entry(const int& row, const int& column, T& value) ;

    /** {@inherit_doc} */
    //override
    public virtual void add_to_entry(const int& row, const int& column, const T& increment) ;

    /** {@inherit_doc} */
    //override
    public virtual void multiply_entry(const int& row, const int& column, const T& factor) ;

    /** {@inherit_doc} */
    //override
    public Field_Matrix<T> transpose() 
    {
        const int n_rows = get_row_dimension();
        const int n_cols = get_column_dimension();
        const Field_Matrix<T> out = create_matrix(n_cols, n_rows);
        walk_in_optimized_order(new DefaultField_Matrix_Preserving_Visitor<T>(field.get_zero()) 
        {
            /** {@inherit_doc} */
            //override
            public void visit(const int& row, const int& column, const T& value) 
            {
                out.set_entry(column, row, value);
            }
        });

        return out;
    }

    /** {@inherit_doc} */
    //override
    public bool is_square() 
    {
        return get_column_dimension() == get_row_dimension();
    }

    /** {@inherit_doc} */
    //override
    public virtual int get_row_dimension();

    /** {@inherit_doc} */
    //override
    public virtual int get_column_dimension();

    /** {@inherit_doc} */
    //override
    public T get_trace()  
    {
        const int n_rows = get_row_dimension();
        const int n_cols = get_column_dimension();
        if (n_rows != n_cols) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NON_SQUARE_MATRIX, n_rows, n_cols);
        }
        T trace = field.get_zero();
        for (int i{}; i < n_rows; ++i) 
        {
            trace = trace.add(get_entry(i, i));
        }
        return trace;
    }

    /** {@inherit_doc} */
    //override
    public std::vector<T> operate(const std::vector<T> v)  
    {

        const int n_rows = get_row_dimension();
        const int n_cols = get_column_dimension();
        if (v.size() != n_cols) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, v.size(), n_cols);
        }

        std::vector<T> out = Math_Arrays::build_array(field, n_rows);
        for (int row{}; row < n_rows; ++row)
        {
            T sum = field.get_zero();
            for (int i{}; i < n_cols; ++i) 
            {
                sum = sum.add(get_entry(row, i).multiply(v[i]));
            }
            out[row] = sum;
        }

        return out;
    }

    /** {@inherit_doc} */
    //override
    public Field_Vector<T> operate(const Field_Vector<T>& v)
    {
        if (dynamic_cast<const ArrayField_Vector*>(*v) != nullptr)
        {
            return ArrayField_Vector<T>(field, operate(((ArrayField_Vector<T>) v).get_data_ref()), false);
        }
        const int n_rows = get_row_dimension();
        const int n_cols = get_column_dimension();
        if (v.get_dimension() != n_cols) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, v.get_dimension(), n_cols);
        }

        std::vector<T> out = Math_Arrays::build_array(field, n_rows);
        for (int row{}; row < n_rows; ++row)
        {
            T sum = field.get_zero();
            for (int i{}; i < n_cols; ++i) 
            {
                sum = sum.add(get_entry(row, i).multiply(v.get_entry(i)));
            }
            out[row] = sum;
        }

        return ArrayField_Vector<T>(field, out, false);
    }

    /** {@inherit_doc} */
    //override
    public std::vector<T> pre_multiply(const std::vector<T>& v)  
    {

        const int n_rows = get_row_dimension();
        const int n_cols = get_column_dimension();
        if (v.size() != n_rows) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, v.size(), n_rows);
        }

        std::vector<T> out = Math_Arrays::build_array(field, n_cols);
        for (int col{}; col < n_cols; ++col)
        {
            T sum = field.get_zero();
            for (int i{}; i < n_rows; ++i) 
            {
                sum = sum.add(get_entry(i, col).multiply(v[i]));
            }
            out[col] = sum;
        }

        return out;
    }

    /** {@inherit_doc} */
    //override
    public Field_Vector<T> pre_multiply(const Field_Vector<T>& v)
    {
        if (dynamic_cast<const ArrayField_Vector*>(*v) != nullptr)
        {
            return ArrayField_Vector<T>(field, pre_multiply(((ArrayField_Vector<T>) v).get_data_ref()), false);
        }

        const int n_rows = get_row_dimension();
        const int n_cols = get_column_dimension();
        if (v.get_dimension() != n_rows) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, v.get_dimension(), n_rows);
        }

        const std::vector<T> out = Math_Arrays::build_array(field, n_cols);
        for (int col{};  col < n_cols; ++col) 
        {
            T sum = field.get_zero();
            for (int i{}; i < n_rows; ++i) 
            {
                sum = sum.add(get_entry(i, col).multiply(v.get_entry(i)));
            }
            out[col] = sum;
        }

        return ArrayField_Vector<T>(field, out, false);
    }

    /** {@inherit_doc} */
    //override
    public T walk_in_row_order(const Field_Matrix_Changing_Visitor<T>& visitor) 
    {
        const int rows = get_row_dimension();
        const int columns = get_column_dimension();
        visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
        for (int row{}; row < rows; ++row) 
        {
            for (const int column{}; column < columns; ++column)
            {
                const T old_value = get_entry(row, column);
                const T new_value = visitor.visit(row, column, old_value);
                set_entry(row, column, new_value);
            }
        }
        return visitor.end();
    }

    /** {@inherit_doc} */
    //override
    public T walk_in_row_order(const Field_Matrix_Preserving_Visitor<T> visitor) 
    {
        const int rows = get_row_dimension();
        const int columns = get_column_dimension();
        visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
        for (int row{}; row < rows; ++row) 
        {
            for (const int& column = 0; column < columns; ++column) 
            {
                visitor.visit(row, column, get_entry(row, column));
            }
        }
        return visitor.end();
    }

    /** {@inherit_doc} */
    //override
    public T walk_in_row_order(const Field_Matrix_Changing_Visitor<T> visitor, const int& start_row, const int& end_row, const int& start_column, const int& end_column)
         
        {
        check_sub_matrix_index(start_row, end_row, start_column, end_column);
        visitor.start(get_row_dimension(), get_column_dimension(), start_row, end_row, start_column, end_column);
        for (int row = start_row; row <= end_row; ++row) 
        {
            for (int column = start_column; column <= end_column; ++column) 
            {
                const T old_value = get_entry(row, column);
                const T new_value = visitor.visit(row, column, old_value);
                set_entry(row, column, new_value);
            }
        }
        return visitor.end();
    }

    /** {@inherit_doc} */
    //override
    public T walk_in_row_order(const Field_Matrix_Preserving_Visitor<T> visitor, const int& start_row, const int& end_row, const int& start_column, const int& end_column)
         
        {
        check_sub_matrix_index(start_row, end_row, start_column, end_column);
        visitor.start(get_row_dimension(), get_column_dimension(), start_row, end_row, start_column, end_column);
        for (int row = start_row; row <= end_row; ++row) 
        {
            for (int column = start_column; column <= end_column; ++column) 
            {
                visitor.visit(row, column, get_entry(row, column));
            }
        }
        return visitor.end();
    }

    /** {@inherit_doc} */
    //override
    public T walk_in_column_order(const Field_Matrix_Changing_Visitor<T> visitor) 
    {
        const int rows    = get_row_dimension();
        const int columns = get_column_dimension();
        visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
        for (const int& column = 0; column < columns; ++column) 
        {
            for (int row{}; row < rows; ++row) 
            {
                const T old_value = get_entry(row, column);
                const T new_value = visitor.visit(row, column, old_value);
                set_entry(row, column, new_value);
            }
        }
        return visitor.end();
    }

    /** {@inherit_doc} */
    //override
    public T walk_in_column_order(const Field_Matrix_Preserving_Visitor<T> visitor) 
    {
        const int rows    = get_row_dimension();
        const int columns = get_column_dimension();
        visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
        for (const int& column = 0; column < columns; ++column) 
        {
            for (int row{}; row < rows; ++row) 
            {
                visitor.visit(row, column, get_entry(row, column));
            }
        }
        return visitor.end();
    }

    /** {@inherit_doc} */
    //override
    public T walk_in_column_order(const Field_Matrix_Changing_Visitor<T> visitor, const int& start_row, const int& end_row, const int& start_column, const int& end_column)
     
    {
        check_sub_matrix_index(start_row, end_row, start_column, end_column);
        visitor.start(get_row_dimension(), get_column_dimension(), start_row, end_row, start_column, end_column);
        for (int column = start_column; column <= end_column; ++column) 
        {
            for (int row = start_row; row <= end_row; ++row) 
            {
                const T old_value = get_entry(row, column);
                const T new_value = visitor.visit(row, column, old_value);
                set_entry(row, column, new_value);
            }
        }
        return visitor.end();
    }

    /** {@inherit_doc} */
    //override
    public T walk_in_column_order(const Field_Matrix_Preserving_Visitor<T> visitor, const int& start_row, const int& end_row, const int& start_column, const int& end_column)
     
    {
        check_sub_matrix_index(start_row, end_row, start_column, end_column);
        visitor.start(get_row_dimension(), get_column_dimension(), start_row, end_row, start_column, end_column);
        for (int column = start_column; column <= end_column; ++column) 
        {
            for (int row = start_row; row <= end_row; ++row) 
            {
                visitor.visit(row, column, get_entry(row, column));
            }
        }
        return visitor.end();
    }

    /** {@inherit_doc} */
    //override
    public T walk_in_optimized_order(const Field_Matrix_Changing_Visitor<T> visitor) 
    {
        return walk_in_row_order(visitor);
    }

    /** {@inherit_doc} */
    //override
    public T walk_in_optimized_order(const Field_Matrix_Preserving_Visitor<T> visitor) 
    {
        return walk_in_row_order(visitor);
    }

    /** {@inherit_doc} */
    //override
    public T walk_in_optimized_order(const Field_Matrix_Changing_Visitor<T> visitor, const int& start_row, const int& end_row, const int& start_column, const int& end_column)
    {
        return walk_in_row_order(visitor, start_row, end_row, start_column, end_column);
    }

    /** {@inherit_doc} */
    //override
    public T walk_in_optimized_order(const Field_Matrix_Preserving_Visitor<T> visitor, const int& start_row, const int& end_row, const int& start_column, const int& end_column)
    {
        return walk_in_row_order(visitor, start_row, end_row, start_column, end_column);
    }

    /**
     * Get a string representation for this matrix.
     * @return a string representation for this matrix
     */
    //override
    public std::string to_string() const 
    {
        const int n_rows = get_row_dimension();
        const int n_cols = get_column_dimension();
        const std::stringstreamres = String_Buffer();
        std::string full_class_name = get_class().get_name();
        std::string short_class_name = full_class_name.substring(full_class_name.last_index_of('.') + 1);
        res.append(short_class_name).append('{');

        for (int i{}; i < n_rows; ++i) 
        {
            if (i > 0) 
            {
                res.append(',');
            }
            res.append('{');
            for (int j{}; j < n_cols; ++j) 
            {
                if (j > 0) 
                {
                    res.append(',');
                }
                res.append(get_entry(i, j));
            }
            res.append('}');
        }

        res.append('}');
        return res.to_string();
    }

    /**
     * Returns true iff <code>object</code> is a
     * <code>Field_Matrix</code> instance with the same dimensions as this
     * and all corresponding matrix entries are equal.
     *
     * @param object the object to test equality against.
     * @return true if object equals this
     */
    //override
    public bool equals(const Object object) 
    {
        if (object == this) 
        {
            return true;
        }
        if (!dynamic_cast<const Field_Matrix*>(*object) != nullptr)
        {
            return false;
        }
        Field_Matrix<?> m = (Field_Matrix<?>) object;
        const int n_rows = get_row_dimension();
        const int n_cols = get_column_dimension();
        if (m.get_column_dimension() != n_cols || m.get_row_dimension() != n_rows) 
        {
            return false;
        }
        for (int row{}; row < n_rows; ++row) 
        {
            for (int col{};  col < n_cols; ++col) 
            {
                if (!get_entry(row, col).equals(m.get_entry(row, col))) 
                {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * Computes a hashcode for the matrix.
     *
     * @return hashcode for matrix
     */
    //override
    public int hash_code() 
    {
        int ret = 322562;
        const int n_rows = get_row_dimension();
        const int n_cols = get_column_dimension();
        ret = ret * 31 + n_rows;
        ret = ret * 31 + n_cols;
        for (int row{}; row < n_rows; ++row) 
        {
            for (int col{};  col < n_cols; ++col) 
            {
               ret = ret * 31 + (11 * (row+1) + 17 * (col+1)) * get_entry(row, col).hash_code();
           }
        }
        return ret;
    }

    /**
     * Check if a row index is valid.
     *
     * @param row Row index to check.
     * @ if {@code index} is not valid.
     */
    protected void check_row_index(const int& row)  
    {
        if (row < 0 || row >= get_row_dimension()) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::ROW_INDEX, row, 0, get_row_dimension() - 1);
        }
    }

    /**
     * Check if a column index is valid.
     *
     * @param column Column index to check.
     * @ if {@code index} is not valid.
     */
    protected void check_column_index(const int column)
         
        {
        if (column < 0 || column >= get_column_dimension()) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::COLUMN_INDEX, column, 0, get_column_dimension() - 1);
        }
    }

    /**
     * Check if submatrix ranges indices are valid.
     * Rows and columns are indicated counting from 0 to n-1.
     *
     * @param start_row Initial row index.
     * @param end_row Final row index.
     * @param start_column Initial column index.
     * @param end_column Final column index.
     * @ if the indices are not valid.
     * @ if {@code end_row < start_row} or
     * {@code end_column < start_column}.
     */
    protected void check_sub_matrix_index(const int& start_row, const int& end_row, const int& start_column, const int& end_column)
         
        {
        check_row_index(start_row);
        check_row_index(end_row);
        if (end_row < start_row) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::INITIAL_ROW_AFTER_FINAL_ROW, end_row, start_row, true);
        }

        check_column_index(start_column);
        check_column_index(end_column);
        if (end_column < start_column) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::INITIAL_COLUMN_AFTER_FINAL_COLUMN, end_column, start_column, true);
        }
    }

    /**
     * Check if submatrix ranges indices are valid.
     * Rows and columns are indicated counting from 0 to n-1.
     *
     * @param selected_rows Array of row indices.
     * @param selected_columns Array of column indices.
     * @Null_Argument_Exception if the arrays are {@code NULL}.
     * @ if the arrays have zero length.
     * @ if row or column selections are not valid.
     */
    protected void check_sub_matrix_index(const std::vector<int> selected_rows, const std::vector<int> selected_columns)
        , Null_Argument_Exception 
        {
        if (selected_rows == NULL || selected_columns == NULL) 
            {
            throw std::exception("not implemented");
            //throw Null_Argument_Exception();
        }
        if (selected_rows.size() == 0 || selected_columns.size() == 0) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NO_DATA);
        }

        for (const auto& row : selected_rows) 
        {
            check_row_index(row);
        }
        for (const auto& column : selected_columns) 
        {
            check_column_index(column);
        }
    }

    /**
     * Check if a matrix is addition compatible with the instance.
     *
     * @param m Matrix to check.
     * @ if the matrix is not
     * addition-compatible with instance.
     */
    protected void check_addition_compatible(const Field_Matrix<T> m)
         
        {
        if ((get_row_dimension() != m.get_row_dimension()) ||
            (get_column_dimension() != m.get_column_dimension())) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH_2x2, m.get_row_dimension(), m.get_column_dimension(), get_row_dimension(), get_column_dimension());
        }
    }

    /**
     * Check if a matrix is subtraction compatible with the instance.
     *
     * @param m Matrix to check.
     * @ if the matrix is not
     * subtraction-compatible with instance.
     */
    protected void check_subtraction_compatible(const Field_Matrix<T> m)
         
        {
        if ((get_row_dimension() != m.get_row_dimension()) ||
            (get_column_dimension() != m.get_column_dimension())) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH_2x2, m.get_row_dimension(), m.get_column_dimension(), get_row_dimension(), get_column_dimension());
        }
    }

    /**
     * Check if a matrix is multiplication compatible with the instance.
     *
     * @param m Matrix to check.
     * @ if the matrix is not
     * multiplication-compatible with instance.
     */
    protected void check_multiplication_compatible(const Field_Matrix<T> m)
    {
        if (get_column_dimension() != m.get_row_dimension()) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, get_column_dimension(), m.get_row_dimension());
        }
    }
};