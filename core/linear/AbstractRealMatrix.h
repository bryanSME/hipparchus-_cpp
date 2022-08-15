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
#include <string>
#include "RealMatrix.h"
#include "RealLinearOperator.h"
#include "RealMatrixFormat.h"
#include "MatrixUtils.h"

/**
 * Basic implementation of Real_Matrix methods regardless of the underlying storage.
 * <p>All the methods implemented here use {@link #get_entry(int, int)} to access
 * matrix elements. Derived class can provide faster implementations.</p>
 *
 */
class Abstract_Real_Matrix : public Real_Matrix, public Real_Linear_Operator 
{

private:
    /** Default format. */
    static const Real_Matrix_Format DEFAULT_FORMAT = Real_Matrix_Format.get_real__matrix_format(Locale.US);
    
    // set the minimum fraction digits to 1 to keep compatibility
    //DEFAULT_FORMAT.get_format().set_minimum_fraction_digits(1);
    

protected:
    /**
     * Creates a matrix with no data
     */
    Abstract_Real_Matrix() 
    {
        // This constructor is intentionally empty. Nothing special is needed here.
    }

    /**
     * Create a Real_Matrix with the supplied row and column dimensions.
     *
     * @param row_dimension  the number of rows in the matrix
     * @param column_dimension  the number of columns in the matrix
     * @ if row or column dimension is not positive
     */
    Abstract_Real_Matrix(const int& row_dimension, const int& column_dimension) 
    {
        if (row_dimension < 1) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::AT_LEAST_ONE_ROW);
        }
        if (column_dimension < 1) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::AT_LEAST_ONE_COLUMN);
        }
    }

public:
    /** {@inherit_doc} */
    // //override
    Real_Matrix add(Real_Matrix m)
    {
        Matrix_Utils::check_addition_compatible(this, m);

        const auto row_count    = get_row_dimension();
        const auto column_count = get_column_dimension();
        auto out = create_matrix(row_count, column_count);
        for (int row{}; row < row_count; ++row) 
        {
            for (int col{}; col < column_count; ++col) 
            {
                out.set_entry(row, col, get_entry(row, col) + m.get_entry(row, col));
            }
        }

        return out;
    }

    /** {@inherit_doc} */
    // //override
    Real_Matrix subtract(const Real_Matrix& m)
    {
        Matrix_Utils::check_subtraction_compatible(this, m);

        const auto row_count = get_row_dimension();
        const auto column_count = get_column_dimension();
        auto out = create_matrix(row_count, column_count);
        for (int row{}; row < row_count; ++row) 
        {
            for (int col{}; col < column_count; ++col) 
            {
                out.set_entry(row, col, get_entry(row, col) - m.get_entry(row, col));
            }
        }

        return out;
    }

    /** {@inherit_doc} */
    // //override
    Real_Matrix scalar_add(const double& d) 
    {
        const auto row_count    = get_row_dimension();
        const auto column_count = get_column_dimension();
        auto out = create_matrix(row_count, column_count);
        for (int row{}; row < row_count; ++row) 
        {
            for (int col{}; col < column_count; ++col) 
            {
                out.set_entry(row, col, get_entry(row, col) + d);
            }
        }
        return out;
    }

    /** {@inherit_doc} */
    // //override
    Real_Matrix scalar_multiply(const double& d) 
    {
        const auto row_count = get_row_dimension();
        const auto column_count = get_column_dimension();
        auto out = create_matrix(row_count, column_count);
        for (int row{}; row < row_count; ++row) 
        {
            for (int col{}; col < column_count; ++col) 
            {
                out.set_entry(row, col, get_entry(row, col) * d);
            }
        }

        return out;
    }

    /** {@inherit_doc} */
    // //override
    Real_Matrix multiply(const Real_Matrix& m)
    {
        Matrix_Utils::check_multiplication_compatible(this, m);

        const auto n_rows = get_row_dimension();
        const auto n_cols = m.get_column_dimension();
        const auto n_sum  = get_column_dimension();
        auto out = create_matrix(n_rows, n_cols);
        for (int row{}; row < n_rows; ++row) 
        {
            for (int col{}; col < n_cols; ++col) 
            {
                double sum{};
                for (int i{}; i < n_sum; ++i) 
                {
                    sum += get_entry(row, i) * m.get_entry(i, col);
                }
                out.set_entry(row, col, sum);
            }
        }

        return out;
    }

    /** {@inherit_doc} */
    // //override
    Real_Matrix pre_multiply(const Real_Matrix& m)
    {
        return m.multiply(this);
    }

    /** {@inherit_doc} */
    // //override
    Real_Matrix power(const int p)
    {
        if (p < 0) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::NOT_POSITIVE_EXPONENT, p);
        }

        if (!is_square()) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::NON_SQUARE_MATRIX, get_row_dimension(), get_column_dimension());
        }

        if (p == 0) 
        {
            return Matrix_Utils::create_real_identity_matrix(this.get_row_dimension());
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
        int max_i = -1;

        for (int i{}; i < binary_representation.size(); ++i) 
        {
            if (binary_representation[i] == '1') 
            {
                const int pos = binary_representation.size() - i - 1;
                non_zero_positions.add(pos);

                // The positions are taken in turn, so max_i is only changed once
                if (max_i == -1) 
                {
                    max_i = pos;
                }
            }
        }

        auto results = Real_Matrix[max_i + 1];
        results[0] = this.copy();

        for (int i{ 1 }; i <= max_i; ++i) 
        {
            results[i] = results[i-1].multiply(results[i-1]);
        }

        Real_Matrix result = this.copy();

        for (auto i : non_zero_positions) 
        {
            result = result.multiply(results[i]);
        }

        return result;
    }

    /** {@inherit_doc} */
    // //override
    std::vector<std::vector<double>> get_data() 
    {
        const auto data = std::vector<std::vector<double>>(get_row_dimension(), std::vector<double>(get_column_dimension());

        for (int i{}; i < data.size(); ++i) 
        {
            auto data_i = data[i];
            for (int j{}; j < data_i.size(); ++j) 
            {
                data_i[j] = get_entry(i, j);
            }
        }

        return data;
    }

    /** {@inherit_doc} */
    // //override
    double get_frobenius_norm() 
    {
        return walk_in_optimized_order(
            Real_Matrix_Preserving_Visitor() 
            {
                /** Sum of squared entries. */
                private double my_sum;

                /** {@inherit_doc} */
                // //override
                public void start(const int& rows, const int& columns, const int& start_row, const int& end_row, const int& start_column, const int& end_column) 
                {
                    my_sum = 0;
                }

                /** {@inherit_doc} */
                // //override
                public void visit(const int& row, const int& column, const double value) 
                {
                    my_sum += value * value;
                }

                /** {@inherit_doc} */
                // //override
                public double end() const
                {
                    return std::sqrt(my_sum);
                }
            }
        );
    }

    /** {@inherit_doc} */
    // //override
    Real_Matrix get_sub_matrix(const int& start_row, const int& end_row, const int& start_column, const int& end_column)
    {
        Matrix_Utils::check_sub_matrix_index(this, start_row, end_row, start_column, end_column);

        const Real_Matrix sub_matrix = create_matrix(end_row - start_row + 1, end_column - start_column + 1);
        for (int i{ start_row }; i <= end_row; ++i)
        {
            for (int j{ start_column }; j <= end_column; ++j)
            {
                sub_matrix.set_entry(i - start_row, j - start_column, get_entry(i, j));
            }
        }

        return sub_matrix;
    }

    /** {@inherit_doc} */
    // //override
    Real_Matrix get_sub_matrix(const std::vector<int> selected_rows, const std::vector<int> selected_columns)
    {
        Matrix_Utils::check_sub_matrix_index(this, selected_rows, selected_columns);

        const auto sub_matrix = create_matrix(selected_rows.size(), selected_columns.size());
        sub_matrix.walk_in_optimized_order(
            new Default_Real_Matrix_Changing_Visitor() 
            {

                /** {@inherit_doc} */
                // //override
                public double visit(const int& row, const int& column, const double value) 
                {
                    return get_entry(selected_rows[row], selected_columns[column]);
            }

            }
        );

        return sub_matrix;
    }

    /** {@inherit_doc} */
    // //override
    void copy_sub_matrix(const int& start_row, const int& end_row, const int& start_column, const int& end_column, const std::vector<std::vector<double>> destination)
         
        {
        Matrix_Utils::check_sub_matrix_index(this, start_row, end_row, start_column, end_column);
        const int rows_count    = end_row + 1 - start_row;
        const int columns_count = end_column + 1 - start_column;
        if ((destination.size() < rows_count) || (destination[0].size() < columns_count)) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH_2x2, destination.size(), destination[0].size(), rows_count, columns_count);
        }

        for (int i{ 1 }; i < rows_count; i++) 
        {
            if (destination[i].size() < columns_count) 
            {
                throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH_2x2, destination.size(), destination[i].size(), rows_count, columns_count);
            }
        }

        walk_in_optimized_order(new DefaultReal_Matrix_Preserving_Visitor() 
        {

            /** Initial row index. */
            private int start_row;

            /** Initial column index. */
            private int start_column;

            /** {@inherit_doc} */
            // //override
            public void start(const int rows, const int columns, const int& start_row, const int& end_row, const int& start_column, const int& end_column) 
            {
                this.start_row    = start_row;
                this.start_column = start_column;
            }

            /** {@inherit_doc} */
            // //override
            public void visit(const int& row, const int& column, const double value) 
            {
                destination[row - start_row][column - start_column] = value;
            }

        }, start_row, end_row, start_column, end_column);
    }

    /** {@inherit_doc} */
    // //override
    void copy_sub_matrix(std::vector<int> selected_rows, std::vector<int> selected_columns, std::vector<std::vector<double>> destination)
    {
        Matrix_Utils::check_sub_matrix_index(this, selected_rows, selected_columns);
        const int& n_cols = selected_columns.size();
        if ((destination.size() < selected_rows.size()) ||
            (destination[0].size() < n_cols)) 
            {
            throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH_2x2, destination.size(), destination[0].size(), selected_rows.size(), selected_columns.size());
        }

        for (int i{}; i < selected_rows.size(); i++) 
        {
            const auto destination_i = destination[i];
            if (destination_i.size() < n_cols) 
            {
                throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH_2x2, destination.size(), destination_i.size(), selected_rows.size(), selected_columns.size());
            }
            for (int j{}; j < selected_columns.size(); j++) 
            {
                destination_i[j] = get_entry(selected_rows[i], selected_columns[j]);
            }
        }
    }

    /** {@inherit_doc} */
    // //override
    void set_sub_matrix(const std::vector<std::vector<double>> sub_matrix, const int& row, const int column) 
    {
        //Math_Utils::check_not_null(sub_matrix);
        const auto n_rows = sub_matrix.size();
        if (n_rows == 0) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::AT_LEAST_ONE_ROW);
        }

        const int n_cols = sub_matrix[0].size();
        if (n_cols == 0) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::AT_LEAST_ONE_COLUMN);
        }

        for (const int r{ 1 }; r < n_rows; ++r)
        {
            if (sub_matrix[r].size() != n_cols) 
            {
                throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, n_cols, sub_matrix[r].size());
            }
        }

        Matrix_Utils::check_row_index(this, row);
        Matrix_Utils::check_column_index(this, column);
        Matrix_Utils::check_row_index(this, n_rows + row - 1);
        Matrix_Utils::check_column_index(this, n_cols + column - 1);

        for (int i{}; i < n_rows; ++i) 
        {
            for (int j{}; j < n_cols; ++j) 
            {
                set_entry(row + i, column + j, sub_matrix[i][j]);
            }
        }
    }

    /** {@inherit_doc} */
    // //override
    Real_Matrix get_row_matrix(const int& row)  
    {
        Matrix_Utils::check_row_index(this, row);
        const auto n_cols = get_column_dimension();
        auto out = create_matrix(1, n_cols);
        for (int i{}; i < n_cols; ++i) 
        {
            out.set_entry(0, i, get_entry(row, i));
        }

        return out;
    }

    /** {@inherit_doc} */
    // //override
    void set_row_matrix(const int& row, const Real_Matrix matrix)
    {
        Matrix_Utils::check_row_index(this, row);
        const auto n_cols = get_column_dimension();
        if ((matrix.get_row_dimension() != 1) || (matrix.get_column_dimension() != n_cols)) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH_2x2, matrix.get_row_dimension(), matrix.get_column_dimension(), 1, n_cols);
        }
        for (int i{}; i < n_cols; ++i) 
        {
            set_entry(row, i, matrix.get_entry(0, i));
        }
    }

    /** {@inherit_doc} */
    // //override
    Real_Matrix get_column_matrix(const int column)
         
        {
        Matrix_Utils::check_column_index(this, column);
        const int n_rows = get_row_dimension();
        auto out = create_matrix(n_rows, 1);
        for (int i{}; i < n_rows; ++i) 
        {
            out.set_entry(i, 0, get_entry(i, column));
        }

        return out;
    }

    /** {@inherit_doc} */
    // //override
    void set_column_matrix(const int& column, const Real_Matrix matrix)
         
        {
        Matrix_Utils::check_column_index(this, column);
        const int n_rows = get_row_dimension();
        if ((matrix.get_row_dimension() != n_rows) ||
            (matrix.get_column_dimension() != 1)) 
            {
            throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH_2x2, matrix.get_row_dimension(), matrix.get_column_dimension(), n_rows, 1);
        }
        for (int i{}; i < n_rows; ++i) 
        {
            set_entry(i, column, matrix.get_entry(i, 0));
        }
    }

    /** {@inherit_doc} */
    // //override
    Real_Vector get_row_vector(const int row)
         
        {
        return Array_Real_Vector(get_row(row), false);
    }

    /** {@inherit_doc} */
    // //override
    void set_row_vector(const int& row, const Real_Vector vector)
         
        {
        Matrix_Utils::check_row_index(this, row);
        const int n_cols = get_column_dimension();
        if (vector.get_dimension() != n_cols) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH_2x2, 1, vector.get_dimension(), 1, n_cols);
        }
        for (int i{}; i < n_cols; ++i) 
        {
            set_entry(row, i, vector.get_entry(i));
        }
    }

    /** {@inherit_doc} */
    // //override
    Real_Vector get_column_vector(const int column)
         
        {
        return Array_Real_Vector(get_column(column), false);
    }

    /** {@inherit_doc} */
    // //override
    void set_column_vector(const int& column, const Real_Vector vector)
         
        {
        Matrix_Utils::check_column_index(this, column);
        const int n_rows = get_row_dimension();
        if (vector.get_dimension() != n_rows) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH_2x2, vector.get_dimension(), 1, n_rows, 1);
        }
        for (int i{}; i < n_rows; ++i) 
        {
            set_entry(i, column, vector.get_entry(i));
        }
    }

    /** {@inherit_doc} */
    // //override
    std::vector<double> get_row(const int row)  
    {
        Matrix_Utils::check_row_index(this, row);
        const int n_cols = get_column_dimension();
        const auto out = std::vector<double>(n_cols);
        for (int i{}; i < n_cols; ++i) 
        {
            out[i] = get_entry(row, i);
        }

        return out;
    }

    /** {@inherit_doc} */
    // //override
    void set_row(const int& row, const std::vector<double> array)
    {
        Matrix_Utils::check_row_index(this, row);
        const int n_cols = get_column_dimension();
        if (array.size() != n_cols) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH_2x2, 1, array.size(), 1, n_cols);
        }
        for (int i{}; i < n_cols; ++i) 
        {
            set_entry(row, i, array[i]);
        }
    }

    /** {@inherit_doc} */
    // //override
    std::vector<double> get_column(const int& column)  
    {
        Matrix_Utils::check_column_index(this, column);
        const int n_rows = get_row_dimension();
        auto out = std::vector<double>(n_rows);
        for (int i{}; i < n_rows; ++i) 
        {
            out[i] = get_entry(i, column);
        }

        return out;
    }

    /** {@inherit_doc} */
    // //override
    void set_column(const int& column, const std::vector<double> array)
    {
        Matrix_Utils::check_column_index(this, column);
        const int n_rows = get_row_dimension();
        if (array.size() != n_rows) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH_2x2, array.size(), 1, n_rows, 1);
        }
        for (int i{}; i < n_rows; ++i) 
        {
            set_entry(i, column, array[i]);
        }
    }

    /** {@inherit_doc} */
    // //override
    void add_to_entry(const int& row, const int& column, double increment)
    {
        Matrix_Utils::check_matrix_index(this, row, column);
        set_entry(row, column, get_entry(row, column) + increment);
    }

    /** {@inherit_doc} */
    // //override
    void multiply_entry(const int& row, const int& column, double factor)
    {
        Matrix_Utils::check_matrix_index(this, row, column);
        set_entry(row, column, get_entry(row, column) * factor);
    }

    /** {@inherit_doc} */
    // //override
    Real_Matrix transpose() 
    {
        const int n_rows = get_row_dimension();
        const int n_cols = get_column_dimension();
        auto out = create_matrix(n_cols, n_rows);
        walk_in_optimized_order(new DefaultReal_Matrix_Preserving_Visitor() 
        {

            /** {@inherit_doc} */
            // //override
            public void visit(const int& row, const int& column, const double value) 
            {
                out.set_entry(column, row, value);
            }

        });

        return out;
    }

    /** {@inherit_doc} */
    // //override
    bool is_square() 
    {
        return get_column_dimension() == get_row_dimension();
    }

    /**
     * Returns the number of rows of this matrix.
     *
     * @return the number of rows.
     */
    // //override
    virtual int get_row_dimension() = 0;

    /**
     * Returns the number of columns of this matrix.
     *
     * @return the number of columns.
     */
    // //override
    virtual int get_column_dimension() = 0;

    /** {@inherit_doc} */
    // //override
    double get_trace()  
    {
        const int n_rows = get_row_dimension();
        const int n_cols = get_column_dimension();
        if (n_rows != n_cols) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::NON_SQUARE_MATRIX, n_rows, n_cols);
       }
        double trace = 0;
        for (int i{}; i < n_rows; ++i) 
        {
            trace += get_entry(i, i);
        }
        return trace;
    }

    /** {@inherit_doc} */
    // //override
    std::vector<double> operate(const std::vector<double>& v)
         
        {
        const int n_rows = get_row_dimension();
        const int n_cols = get_column_dimension();
        if (v.size() != n_cols) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, v.size(), n_cols);
        }

        auto out = std::vector<double>(n_rows);
        for (int row{}; row < n_rows; ++row) 
        {
            double sum{};
            for (int i{}; i < n_cols; ++i) 
            {
                sum += get_entry(row, i) * v[i];
            }
            out[row] = sum;
        }

        return out;
    }

    /** {@inherit_doc} */
    // //override
    Real_Vector operate(const Real_Vector& v)
    {
        if (dynamic_cast<const Array_Real_Vector*>(*v) != nullptr)
        {
            return Array_Real_Vector(operate(((Array_Real_Vector) v).get_data_ref()), false);
        }

        const int n_rows = get_row_dimension();
        const int n_cols = get_column_dimension();
        if (v.get_dimension() != n_cols) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, v.get_dimension(), n_cols);
        }

        auto out = std::vector<double>(n_rows;
        for (int row{}; row < n_rows; ++row) 
        {
            double sum{};
            for (int i{}; i < n_cols; ++i) 
            {
                sum += get_entry(row, i) * v.get_entry(i);
            }
            out[row] = sum;
        }
        return Array_Real_Vector(out, false);
        
    }

    /** {@inherit_doc} */
    // //override
    std::vector<double> pre_multiply(const std::vector<double>& v)  
    {
        const int n_rows = get_row_dimension();
        const int n_cols = get_column_dimension();
        if (v.size() != n_rows) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, v.size(), n_rows);
        }

        auto out = std::vector<double>(n_cols);
        for (int col{}; col < n_cols; ++col)
        {
            double sum{};
            for (int i{}; i < n_rows; ++i) 
            {
                sum += get_entry(i, col) * v[i];
            }
            out[col] = sum;
        }

        return out;
    }

    /** {@inherit_doc} */
    // //override
    Real_Vector pre_multiply(const Real_Vector& v)  
    {
        if (dynamic_cast<const Array_Real_Vector*>(*v) != nullptr)
        {
            return Array_Real_Vector(pre_multiply(((Array_Real_Vector) v).get_data_ref()), false);
        }

        const int n_rows = get_row_dimension();
        const int n_cols = get_column_dimension();
        if (v.get_dimension() != n_rows) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, v.get_dimension(), n_rows);
        }

        auto out = std::vector<double>(n_cols);
        for (int col{}; col < n_cols; ++col)
        {
            double sum{};
            for (int i{}; i < n_rows; ++i) 
            {
                sum += get_entry(i, col) * v.get_entry(i);
            }
            out[col] = sum;
        }

        return Array_Real_Vector(out, false);
        
    }

    /** {@inherit_doc} */
    // //override
    double walk_in_row_order(const Real_Matrix_Changing_Visitor& visitor) 
    {
        const int rows    = get_row_dimension();
        const int columns = get_column_dimension();
        visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
        for (int row{}; row < rows; ++row) 
        {
            for (const int& column = 0; column < columns; ++column) 
            {
                const double old_value = get_entry(row, column);
                const double new_value = visitor.visit(row, column, old_value);
                set_entry(row, column, new_value);
            }
        }
        return visitor.end();
    }

    /** {@inherit_doc} */
    // //override
    double walk_in_row_order(const Real_Matrix_Preserving_Visitor& visitor) 
    {
        const int rows    = get_row_dimension();
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
    // //override
    double walk_in_row_order(const Real_Matrix_Changing_Visitor& visitor, const int& start_row, const int& end_row, const int& start_column, const int& end_column)
         
        {
        Matrix_Utils::check_sub_matrix_index(this, start_row, end_row, start_column, end_column);
        visitor.start(get_row_dimension(), get_column_dimension(), start_row, end_row, start_column, end_column);
        for (int row = start_row; row <= end_row; ++row) 
        {
            for (int column = start_column; column <= end_column; ++column) 
            {
                const double old_value = get_entry(row, column);
                const double new_value = visitor.visit(row, column, old_value);
                set_entry(row, column, new_value);
            }
        }
        return visitor.end();
    }

    /** {@inherit_doc} */
    // //override
    double walk_in_row_order(const Real_Matrix_Preserving_Visitor& visitor, const int& start_row, const int& end_row, const int& start_column, const int& end_column)
    {
        Matrix_Utils::check_sub_matrix_index(this, start_row, end_row, start_column, end_column);
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
    // //override
    double walk_in_column_order(const Real_Matrix_Changing_Visitor& visitor) 
    {
        const int rows    = get_row_dimension();
        const int columns = get_column_dimension();
        visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
        for (const int& column = 0; column < columns; ++column) 
        {
            for (int row{}; row < rows; ++row) 
            {
                const double old_value = get_entry(row, column);
                const double new_value = visitor.visit(row, column, old_value);
                set_entry(row, column, new_value);
            }
        }
        return visitor.end();
    }

    /** {@inherit_doc} */
    // //override
    double walk_in_column_order(const Real_Matrix_Preserving_Visitor& visitor) 
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
    // //override
    double walk_in_column_order(const Real_Matrix_Changing_Visitor& visitor, const int& start_row, const int& end_row, const int& start_column, const int& end_column)
    {
        Matrix_Utils::check_sub_matrix_index(this, start_row, end_row, start_column, end_column);
        visitor.start(get_row_dimension(), get_column_dimension(), start_row, end_row, start_column, end_column);
        for (int column = start_column; column <= end_column; ++column) 
        {
            for (int row = start_row; row <= end_row; ++row) 
            {
                const double old_value = get_entry(row, column);
                const double new_value = visitor.visit(row, column, old_value);
                set_entry(row, column, new_value);
            }
        }
        return visitor.end();
    }

    /** {@inherit_doc} */
    // //override
     double walk_in_column_order(const Real_Matrix_Preserving_Visitor& visitor, const int& start_row, const int& end_row, const int& start_column, const int& end_column)
     {
        Matrix_Utils::check_sub_matrix_index(this, start_row, end_row, start_column, end_column);
        visitor.start(get_row_dimension(), get_column_dimension(), start_row, end_row, start_column, end_column);
        for (int column = start_column; column <= end_column; ++column) 
        {
            for (int row{ start_row }; row <= end_row; ++row)
            {
                visitor.visit(row, column, get_entry(row, column));
            }
        }
        return visitor.end();
    }

    /** {@inherit_doc} */
    // //override
    double walk_in_optimized_order(const Real_Matrix_Changing_Visitor& visitor) 
    {
        return walk_in_row_order(visitor);
    }

    /** {@inherit_doc} */
    // //override
    double walk_in_optimized_order(const Real_Matrix_Preserving_Visitor& visitor) 
    {
        return walk_in_row_order(visitor);
    }

    /** {@inherit_doc} */
    // //override
    double walk_in_optimized_order(const Real_Matrix_Changing_Visitor& visitor, const int& start_row, const int& end_row, const int& start_column, const int& end_column)
    {
        return walk_in_row_order(visitor, start_row, end_row, start_column, end_column);
    }

    /** {@inherit_doc} */
    // //override
    double walk_in_optimized_order(const Real_Matrix_Preserving_Visitor& visitor, const int& start_row, const int& end_row, const int& start_column, const int& end_column)
    {
        return walk_in_row_order(visitor, start_row, end_row, start_column, end_column);
    }

    /**
     * Get a string representation for this matrix.
     * @return a string representation for this matrix
     */
    // //override
    std::string to_string() const 
    {
        std::stringstream res_ss = std::stringstream();
        std::string full_class_name = get_class().get_name();
        std::string short_class_name = full_class_name.substring(full_class_name.last_index_of('.') + 1);
        res_ss << (short_class_name) << DEFAULT_FORMAT.format(this);
        return res_ss.str();
    }

    /**
     * Returns true iff <code>object</code> is a
     * <code>Real_Matrix</code> instance with the same dimensions as this
     * and all corresponding matrix entries are equal.
     *
     * @param object the object to test equality against.
     * @return true if object equals this
     */
    // //override
    bool equals(const Object& object) 
    {
        if (object == this) 
        {
            return true;
        }
        if (!dynamic_cast<const Real_Matrix*>(*object) != nullptr)
        {
            return false;
        }
        auto m = static_cast<Real_Matrix>(object);
        const int n_rows = get_row_dimension();
        const int n_cols = get_column_dimension();
        if (m.get_column_dimension() != n_cols || m.get_row_dimension() != n_rows) 
        {
            return false;
        }
        for (int row{}; row < n_rows; ++row) 
        {
            for (int col{}; col < n_cols; ++col) 
            {
                if (get_entry(row, col) != m.get_entry(row, col)) 
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
    // //override
    int hash_code() 
    {
        int ret = 7;
        const int n_rows = get_row_dimension();
        const int n_cols = get_column_dimension();
        ret = ret * 31 + n_rows;
        ret = ret * 31 + n_cols;
        for (int row{}; row < n_rows; ++row) 
        {
            for (int col{}; col < n_cols; ++col) 
            {
               ret = ret * 31 + (11 * (row+1) + 17 * (col+1)) *
                   Math_Utils::hash(get_entry(row, col));
           }
        }
        return ret;
    }


    /*
     * Empty implementations of these methods are provided in order to allow for
     * the use of the // //override tag with Java 1.5.
     */

    /** {@inherit_doc} */
    // //override
    virtual Real_Matrix create_matrix(const int& row_dimension, const int& column_dimension) = 0;

    /** {@inherit_doc} */
    // //override
    virtual Real_Matrix copy() = 0;

    /** {@inherit_doc} */
    // //override
    virtual double get_entry(const int& row, const int& column) = 0;

    /** {@inherit_doc} */
    // //override
    virtual void set_entry(const int& row, const int& column, double value) = 0;
}


