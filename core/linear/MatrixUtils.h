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

//import java.io.IOException;
//import java.io.Object_Input_Stream;
//import java.io.Object_Output_Stream;
//import java.util.Array_list;
//import java.util.Arrays;
//import java.util.List;

//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.Field;
//import org.hipparchus.Field_Element;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Runtime_Exception;
//import org.hipparchus.exception.Null_Argument_Exception;
//import org.hipparchus.fraction.Big_Fraction;
//import org.hipparchus.fraction.Fraction;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Arrays;
//import org.hipparchus.util.Math_Utils;
//import org.hipparchus.util.Precision;
#include <vector>
#include <cmath>
#include "RealMatrixFormat.h"
#include "DecompositionSolver.h"
#include <type_traits>
#include "../CalculusFieldElement.hpp"

/**
 * A collection of static methods that operate on or return matrices.
 *
 */
class Matrix_Utils 
{

private:
    /** Pade coefficients required for the matrix exponential calculation. */
    static const std::vector<double> PADE_COEFFICIENTS_3 = 
    {
            120.0, 60.0, 12.0, 1.0
    };

    /** Pade coefficients required for the matrix exponential calculation. */
    static const std::vector<double> PADE_COEFFICIENTS_5 = 
    {
            30240.0, 15120.0, 3360.0, 420.0, 30.0, 1
    };

    /** Pade coefficients required for the matrix exponential calculation. */
    static const std::vector<double> PADE_COEFFICIENTS_7 = 
    {
            17297280.0, 8648640.0, 1995840.0, 277200.0, 25200.0, 1512.0, 56.0, 1.0
    };

    /** Pade coefficients required for the matrix exponential calculation. */
    static const std::vector<double> PADE_COEFFICIENTS_9 = 
    {
            17643225600.0, 8821612800.0, 2075673600.0, 302702400.0, 30270240.0, 2162160.0, 110880.0, 3960.0, 90.0, 1.0
    };

    /** Pade coefficients required for the matrix exponential calculation. */
    static const std::vector<double> PADE_COEFFICIENTS_13 = 
    {
            6.476475253248e+16, 3.238237626624e+16, 7.7717703038976e+15, 1.1873537964288e+15, 129060195264000.0, 10559470521600.0, 670442572800.0, 33522128640.0, 1323241920.0, 40840800.0, 960960.0, 16380.0, 182.0, 1.0
    };

    /**
     * Private constructor.
     */
    Matrix_Utils() 
    {
        //super();
    }

    /**
     * Checks whether a matrix is symmetric, within a given relative tolerance.
     *
     * @param matrix Matrix to check.
     * @param relative_tolerance Tolerance of the symmetry check.
     * @param raiseException If {@code true}, an exception will be raised if
     * the matrix is not symmetric.
     * @return {@code true} if {@code matrix} is symmetric.
     * @ if the matrix is not square.
     * @ if the matrix is not symmetric.
     */
    private static bool is_symmetricInternal(Real_Matrix matrix, double relative_tolerance, bool raiseException)
    {
        const int rows = matrix.get_row_dimension();
        if (rows != matrix.get_column_dimension())
        {
            if (raiseException)
            {
                throw std::exception("not implemented");
                // throw (hipparchus::exception::Localized_Core_Formats_Type::NON_SQUARE_MATRIX, rows, matrix.get_column_dimension());
            }
            else
            {
                return false;
            }
        }
        for (int i{}; i < rows; i++)
        {
            for (int j = i + 1; j < rows; j++)
            {
                const double mij = matrix.get_entry(i, j);
                const double mji = matrix.get_entry(j, i);
                if (std::abs(mij - mji) >
                    std::max(std::abs(mij), std::abs(mji)) * relative_tolerance)
                {
                    if (raiseException)
                    {
                        throw std::exception("not implemented");
                        //throw (hipparchus::exception::Localized_Core_Formats_Type::NON_SYMMETRIC_MATRIX, i, j, relative_tolerance);
                    }
                    else
                    {
                        return false;
                    }
                }
            }
        }
        return true;
    }

public:
    /**
     * The default format for {@link Real_Matrix} objects.
     */
    static const Real_Matrix_Format DEFAULT_FORMAT = Real_Matrix_Format.get_real__matrix_format();

    /**
     * A format for {@link Real_Matrix} objects compatible with octave.
     */
    static const Real_Matrix_Format OCTAVE_FORMAT = Real_Matrix_Format("[", "]", "", "", "; ", ", ");


    /**
     * Returns a {@link Real_Matrix} with specified dimensions.
     * <p>The type of matrix returned depends on the dimension. Below
     * 2<sup>12</sup> elements (i.e. 4096 elements or 64&times;64 for a
     * square matrix) which can be stored in a 32kB array, a {@link
     * Array_2D_Row_Real_Matrix} instance is built. Above this threshold a {@link
     * Block_Real_Matrix} instance is built.</p>
     * <p>The matrix elements are all set to 0.0.</p>
     * @param rows number of rows of the matrix
     * @param columns number of columns of the matrix
     * @return  Real_Matrix with specified dimensions
     * @see #create_real_matrix(std::vector<std::vector<double>>)
     */
    static Real_Matrix create_real_matrix(const int rows, const int columns) 
    {
        return (rows * columns <= 4096) ?
                Array_2D_Row_Real_Matrix(rows, columns) : Block_Real_Matrix(rows, columns);
    }

    /**
     * Returns a {@link Field_Matrix} with specified dimensions.
     * <p>The type of matrix returned depends on the dimension. Below
     * 2<sup>12</sup> elements (i.e. 4096 elements or 64&times;64 for a
     * square matrix), a {@link Field_Matrix} instance is built. Above
     * this threshold a {@link BlockField_Matrix} instance is built.</p>
     * <p>The matrix elements are all set to field.get_zero().</p>
     * @param <T> the type of the field elements
     * @param field field to which the matrix elements belong
     * @param rows number of rows of the matrix
     * @param columns number of columns of the matrix
     * @return  Field_Matrix with specified dimensions
     * @see #create_field_matrix(Field_Element[][])
     */
    static <T extends Field_Element<T>> Field_Matrix<T> create_field_matrix(const Field<T> field, const int rows, const int columns) 
    {
        return (rows * columns <= 4096) ?
                Array2DRowField_Matrix<T>(field, rows, columns) : BlockField_Matrix<T>(field, rows, columns);
    }

    /**
     * Returns a {@link Real_Matrix} whose entries are the the values in the
     * the input array.
     * <p>The type of matrix returned depends on the dimension. Below
     * 2<sup>12</sup> elements (i.e. 4096 elements or 64&times;64 for a
     * square matrix) which can be stored in a 32kB array, a {@link
     * Array_2D_Row_Real_Matrix} instance is built. Above this threshold a {@link
     * Block_Real_Matrix} instance is built.</p>
     * <p>The input array is copied, not referenced.</p>
     *
     * @param data input array
     * @return  Real_Matrix containing the values of the array
     * @org.hipparchus.exception.
     * if {@code data} is not rectangular (not all rows have the same length).
     * @ if a row or column is empty.
     * @Null_Argument_Exception if either {@code data} or {@code data[0]}
     * is {@code NULL}.
     * @ if {@code data} is not rectangular.
     * @see #create_real_matrix(int, int)
     */
    static Real_Matrix create_real_matrix(std::vector<std::vector<double>> data)
    {
        if (data == NULL || data[0] == NULL) 
        {
            throw std::exception("not implemented");
            //throw Null_Argument_Exception();
        }
        return (data.size() * data[0].size() <= 4096)
            ? Array_2D_Row_Real_Matrix(data)
            : Block_Real_Matrix(data);
    }

    /**
     * Returns a {@link Field_Matrix} whose entries are the the values in the
     * the input array.
     * <p>The type of matrix returned depends on the dimension. Below
     * 2<sup>12</sup> elements (i.e. 4096 elements or 64&times;64 for a
     * square matrix), a {@link Field_Matrix} instance is built. Above
     * this threshold a {@link BlockField_Matrix} instance is built.</p>
     * <p>The input array is copied, not referenced.</p>
     * @param <T> the type of the field elements
     * @param data input array
     * @return a matrix containing the values of the array.
     * @org.hipparchus.exception.
     * if {@code data} is not rectangular (not all rows have the same length).
     * @ if a row or column is empty.
     * @Null_Argument_Exception if either {@code data} or {@code data[0]}
     * is {@code NULL}.
     * @see #create_field_matrix(Field, int, int)
     */
    static <T extends Field_Element<T>> Field_Matrix<T> create_field_matrix(std::vector<std::vector<T>> data)
    {
        if (data == NULL || data[0] == NULL) 
        {
            throw std::exception("not implemented"); 
            //throw Null_Argument_Exception();
        }
        return (data.size() * data[0].size() <= 4096)
            ? Array2DRowField_Matrix<T>(data)
            : BlockField_Matrix<T>(data);
    }

    /**
     * Returns <code>dimension x dimension</code> identity matrix.
     *
     * @param dimension dimension of identity matrix to generate
     * @return identity matrix
     * @Illegal_Argument_Exception if dimension is not positive
     */
    static Real_Matrix create_real_identity_matrix(const int& dimension) 
    {
        const Real_Matrix m = create_real_matrix(dimension, dimension);
        for (int i{}; i < dimension; ++i) 
        {
            m.set_entry(i, i, 1.0);
        }
        return m;
    }

    /**
     * Returns <code>dimension x dimension</code> identity matrix.
     *
     * @param <T> the type of the field elements
     * @param field field to which the elements belong
     * @param dimension dimension of identity matrix to generate
     * @return identity matrix
     * @Illegal_Argument_Exception if dimension is not positive
     */
    static <T extends Field_Element<T>> Field_Matrix<T>
        create_field_identity_matrix(const Field<T> field, const int& dimension) 
        {
        const T zero = field.get_zero();
        const T one  = field.get_one();
        const std::vector<std::vector<T>> d = Math_Arrays::build_array(field, dimension, dimension);
        for (int row{}; row < dimension; row++) 
        {
            const std::vector<T> d_row = d[row];
            Arrays.fill(d_row, zero);
            d_row[row] = one;
        }
        return Array2DRowField_Matrix<T>(field, d, false);
    }

    /**
     * Returns a diagonal matrix with specified elements.
     *
     * @param diagonal diagonal elements of the matrix (the array elements
     * will be copied)
     * @return diagonal matrix
     */
    static Real_Matrix create_real_diagonal_matrix(const std::vector<double> diagonal) 
    {
        const Real_Matrix m = create_real_matrix(diagonal.size(), diagonal.size());
        for (int i{}; i < diagonal.size(); ++i) 
        {
            m.set_entry(i, i, diagonal[i]);
        }
        return m;
    }

    /**
     * Returns a diagonal matrix with specified elements.
     *
     * @param <T> the type of the field elements
     * @param diagonal diagonal elements of the matrix (the array elements
     * will be copied)
     * @return diagonal matrix
     */
    static <T extends Field_Element<T>> Field_Matrix<T>
        createFieldDiagonal_Matrix(const std::vector<T> diagonal) 
        {
        const Field_Matrix<T> m =
            create_field_matrix(diagonal[0].get_field(), diagonal.size(), diagonal.size());
        for (int i{}; i < diagonal.size(); ++i) 
        {
            m.set_entry(i, i, diagonal[i]);
        }
        return m;
    }

    /**
     * Creates a {@link Real_Vector} using the data from the input array.
     *
     * @param data the input data
     * @return a data.size() Real_Vector
     * @ if {@code data} is empty.
     * @Null_Argument_Exception if {@code data} is {@code NULL}.
     */
    static Real_Vector create_real__vector(const std::vector<double>& data)
    {
        if (data == NULL) 
        {
            throw std::exception("not implemented");
            //throw Null_Argument_Exception();
        }
        return Array_Real_Vector(data, true);
    }

    /**
     * Creates a {@link Real_Vector} with specified dimensions.
     *
     * @param dimension dimension of the vector
     * @return a vector
     * @since 1.3
     */
    static Real_Vector create_real__vector(const int& dimension) 
    {
        return Array_Real_Vector(std::vector<double>(dimension]);
    }

    /**
     * Creates a {@link Field_Vector} using the data from the input array.
     *
     * @param <T> the type of the field elements
     * @param data the input data
     * @return a data.size() Field_Vector
     * @ if {@code data} is empty.
     * @Null_Argument_Exception if {@code data} is {@code NULL}.
     * @ if {@code data} has 0 elements
     */
    static <T extends Field_Element<T>> Field_Vector<T> create_field_vector(const std::vector<T> data)
    {
        if (data == NULL) 
        {
            throw std::exception("not implemented");
            //throw Null_Argument_Exception();
        }
        if (data.size() == 0) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::VECTOR_MUST_HAVE_AT_LEAST_ONE_ELEMENT);
        }
        return ArrayField_Vector<T>(data[0].get_field(), data, true);
    }

    /**
     * Creates a {@link Field_Vector} with specified dimensions.
     *
     * @param <T> the type of the field elements
     * @param field field to which array elements belong
     * @param dimension dimension of the vector
     * @return a vector
     * @since 1.3
     */
    static <T extends Field_Element<T>> Field_Vector<T> create_field_vector(const Field<T> field, const int& dimension) 
    {
        return ArrayField_Vector<>(Math_Arrays::build_array(field, dimension));
    }

    /**
     * Create a row {@link Real_Matrix} using the data from the input
     * array.
     *
     * @param row_data the input row data
     * @return a 1 x row_data.size() Real_Matrix
     * @ if {@code row_data} is empty.
     * @Null_Argument_Exception if {@code row_data} is {@code NULL}.
     */
    static Real_Matrix createRowReal_Matrix(const std::vector<double>& row_data) 
    {
        if (row_data == NULL) 
        {
            throw std::exception("not implemented");
            //throw Null_Argument_Exception();
        }
        const int n_cols = row_data.size();
        const Real_Matrix m = create_real_matrix(1, n_cols);
        for (int i{}; i < n_cols; ++i) 
        {
            m.set_entry(0, i, row_data[i]);
        }
        return m;
    }

    /**
     * Create a row {@link Field_Matrix} using the data from the input
     * array.
     *
     * @param <T> the type of the field elements
     * @param row_data the input row data
     * @return a 1 x row_data.size() Field_Matrix
     * @ if {@code row_data} is empty.
     * @Null_Argument_Exception if {@code row_data} is {@code NULL}.
     */
    static <T extends Field_Element<T>> Field_Matrix<T> createRowField_Matrix(const std::vector<T> row_data)
    {
        if (row_data == NULL) 
        {
            throw std::exception("not implemented");
            //throw Null_Argument_Exception();
        }
        const int& n_cols = row_data.size();
        if (n_cols == 0) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::AT_LEAST_ONE_COLUMN);
        }
        const Field_Matrix<T> m = create_field_matrix(row_data[0].get_field(), 1, n_cols);
        for (int i{}; i < n_cols; ++i) 
        {
            m.set_entry(0, i, row_data[i]);
        }
        return m;
    }

    /**
     * Creates a column {@link Real_Matrix} using the data from the input
     * array.
     *
     * @param column_data  the input column data
     * @return a column_data x 1 Real_Matrix
     * @ if {@code column_data} is empty.
     * @Null_Argument_Exception if {@code column_data} is {@code NULL}.
     */
    static Real_Matrix create_column_real__matrix(const std::vector<double>& column_data)
    {
        if (column_data == NULL) 
        {
            throw std::exception("not implemented");
            //throw Null_Argument_Exception();
        }
        const int& n_rows = column_data.size();
        const Real_Matrix m = create_real_matrix(n_rows, 1);
        for (int i{}; i < n_rows; ++i) 
        {
            m.set_entry(i, 0, column_data[i]);
        }
        return m;
    }

    /**
     * Creates a column {@link Field_Matrix} using the data from the input
     * array.
     *
     * @param <T> the type of the field elements
     * @param column_data  the input column data
     * @return a column_data x 1 Field_Matrix
     * @ if {@code data} is empty.
     * @Null_Argument_Exception if {@code column_data} is {@code NULL}.
     */
    static <T extends Field_Element<T>> Field_Matrix<T> createColumnField_Matrix(const std::vector<T>& column_data)
    {
        const int n_rows = column_data.size();
        if (n_rows == 0) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::AT_LEAST_ONE_ROW);
        }
        const Field_Matrix<T> m = create_field_matrix(column_data[0].get_field(), n_rows, 1);
        for (int i{}; i < n_rows; ++i) 
        {
            m.set_entry(i, 0, column_data[i]);
        }
        return m;
    }

    

    /**
     * Checks whether a matrix is symmetric.
     *
     * @param matrix Matrix to check.
     * @param eps Relative tolerance.
     * @ if the matrix is not square.
     * @ if the matrix is not symmetric.
     */
    static void checkSymmetric(Real_Matrix matrix, double eps) 
    {
        is_symmetricInternal(matrix, eps, true);
    }

    /**
     * Checks whether a matrix is symmetric.
     *
     * @param matrix Matrix to check.
     * @param eps Relative tolerance.
     * @return {@code true} if {@code matrix} is symmetric.
     */
    static bool is_symmetric(Real_Matrix matrix, double eps) 
    {
        return is_symmetricInternal(matrix, eps, false);
    }

    /**
     * Check if matrix indices are valid.
     *
     * @param m Matrix.
     * @param row Row index to check.
     * @param column Column index to check.
     * @ if {@code row} or {@code column} is not
     * a valid index.
     */
    static void check_matrix_index(const Any_Matrix m, const int& row, const int column)
    {
        check_row_index(m, row);
        check_column_index(m, column);
    }

    /**
     * Check if a row index is valid.
     *
     * @param m Matrix.
     * @param row Row index to check.
     * @ if {@code row} is not a valid index.
     */
    public static void check_row_index(const Any_Matrix m, const int row)
    {
        if (row < 0 || row >= m.get_row_dimension()) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::ROW_INDEX, row, 0, m.get_row_dimension() - 1);
        }
    }

    /**
     * Check if a column index is valid.
     *
     * @param m Matrix.
     * @param column Column index to check.
     * @ if {@code column} is not a valid index.
     */
    static void check_column_index(const Any_Matrix m, const int column)
    {
        if (column < 0 || column >= m.get_column_dimension()) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::COLUMN_INDEX, column, 0, m.get_column_dimension() - 1);
        }
    }

    /**
     * Check if submatrix ranges indices are valid.
     * Rows and columns are indicated counting from 0 to {@code n - 1}.
     *
     * @param m Matrix.
     * @param start_row Initial row index.
     * @param end_row Final row index.
     * @param start_column Initial column index.
     * @param end_column Final column index.
     * @ if the indices are invalid.
     * @ if {@code end_row < start_row} or
     * {@code end_column < start_column}.
     */
    static void check_sub_matrix_index(const Any_Matrix m, const int& start_row, const int& end_row, const int& start_column, const int& end_column)
    {
        check_row_index(m, start_row);
        check_row_index(m, end_row);
        if (end_row < start_row) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::INITIAL_ROW_AFTER_FINAL_ROW, end_row, start_row, false);
        }

        check_column_index(m, start_column);
        check_column_index(m, end_column);
        if (end_column < start_column) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::INITIAL_COLUMN_AFTER_FINAL_COLUMN, end_column, start_column, false);
        }
    }

    /**
     * Check if submatrix ranges indices are valid.
     * Rows and columns are indicated counting from 0 to n-1.
     *
     * @param m Matrix.
     * @param selected_rows Array of row indices.
     * @param selected_columns Array of column indices.
     * @Null_Argument_Exception if {@code selected_rows} or
     * {@code selected_columns} are {@code NULL}.
     * @ if the row or column selections are empty (zero
     * length).
     * @ if row or column selections are not valid.
     */
    static void check_sub_matrix_index(const Any_Matrix m, const std::vector<int> selected_rows, const std::vector<int> selected_columns)
    {
        if (selected_rows == NULL) 
        {
            throw std::exception("not implemented");
            //throw Null_Argument_Exception();
        }
        if (selected_columns == NULL) 
        {
            throw std::exception("not implemented");
            //throw Null_Argument_Exception();
        }
        if (selected_rows.size() == 0) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::EMPTY_SELECTED_ROW_INDEX_ARRAY);
        }
        if (selected_columns.size() == 0) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::EMPTY_SELECTED_COLUMN_INDEX_ARRAY);
        }

        for (const int row : selected_rows) 
        {
            check_row_index(m, row);
        }
        for (const int column : selected_columns) 
        {
            check_column_index(m, column);
        }
    }

    /**
     * Check if matrices are addition compatible.
     *
     * @param left Left hand side matrix.
     * @param right Right hand side matrix.
     * @ if the matrices are not addition
     * compatible.
     */
    static void check_addition_compatible(const Any_Matrix left, const Any_Matrix right)
         
        {
        if ((left.get_row_dimension()    != right.get_row_dimension()) ||
            (left.get_column_dimension() != right.get_column_dimension())) 
            {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH_2x2, left.get_row_dimension(), left.get_column_dimension(), right.get_row_dimension(), right.get_column_dimension());
        }
    }

    /**
     * Check if matrices are subtraction compatible
     *
     * @param left Left hand side matrix.
     * @param right Right hand side matrix.
     * @ if the matrices are not addition
     * compatible.
     */
    static void check_subtraction_compatible(const Any_Matrix left, const Any_Matrix right)
    {
        if ((left.get_row_dimension()    != right.get_row_dimension()) ||
            (left.get_column_dimension() != right.get_column_dimension())) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH_2x2, left.get_row_dimension(), left.get_column_dimension(), right.get_row_dimension(), right.get_column_dimension());
        }
    }

    /**
     * Check if matrices are multiplication compatible
     *
     * @param left Left hand side matrix.
     * @param right Right hand side matrix.
     * @ if matrices are not multiplication
     * compatible.
     */
    static void check_multiplication_compatible(const Any_Matrix left, const Any_Matrix right)
    {
        if (left.get_column_dimension() != right.get_row_dimension()) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, left.get_column_dimension(), right.get_row_dimension());
        }
    }

    /**
     * Check if matrices have the same number of columns.
     *
     * @param left Left hand side matrix.
     * @param right Right hand side matrix.
     * @ if matrices don't have the same number of columns.
     * @since 1.3
     */
    static void check_same_column_dimension(const Any_Matrix left, const Any_Matrix right)
    {
        if (left.get_column_dimension() != right.get_column_dimension()) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, left.get_column_dimension(), right.get_column_dimension());
        }
    }

    /**
     * Check if matrices have the same number of rows.
     *
     * @param left Left hand side matrix.
     * @param right Right hand side matrix.
     * @ if matrices don't have the same number of rows.
     * @since 1.3
     */
    static void check_same_row_dimension(const Any_Matrix left, const Any_Matrix right)
    {
        if (left.get_row_dimension() != right.get_row_dimension()) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, left.get_row_dimension(), right.get_row_dimension());
        }
    }

    /**
     * Convert a {@link Field_Matrix}/{@link Fraction} matrix to a {@link Real_Matrix}.
     * @param m Matrix to convert.
     * @return the converted matrix.
     */
    static Array_2D_Row_Real_Matrix fractionMatrixToReal_Matrix(const Field_Matrix<Fraction> m) 
    {
        const FractionMatrix_converter converter = FractionMatrix_converter();
        m.walk_in_optimized_order(converter);
        return converter.get_convertedMatrix();
    }



    //----------------------------------------------
    
    /** Converter for {@link Field_Matrix}/{@link Fraction}. */
    private static class FractionMatrix_converter extends DefaultField_Matrix_Preserving_Visitor<Fraction> 
    {
        /** Converted array. */
        private std::vector<std::vector<double>> data;
        /** Simple constructor. */
        FractionMatrix_converter() 
        {
            super(Fraction.ZERO);
        }

        /** {@inherit_doc} */
        //override
        public void start(const int& rows, int columns, int start_row, int end_row, int start_column, int end_column) 
        {
            data = std::vector<double>(rows][columns];
        }

        /** {@inherit_doc} */
        //override
        public void visit(const int& row, const int& column, Fraction value) 
        {
            data[row][column] = value.double_value();
        }

        /**
         * Get the converted matrix.
         *
         * @return the converted matrix.
         */
        Array_2D_Row_Real_Matrix get_convertedMatrix() 
        {
            return Array_2D_Row_Real_Matrix(data, false);
        }

    }

    /**
     * Convert a {@link Field_Matrix}/{@link Big_Fraction} matrix to a {@link Real_Matrix}.
     *
     * @param m Matrix to convert.
     * @return the converted matrix.
     */
    public static Array_2D_Row_Real_Matrix big_fraction_matrix_to_real__matrix(const Field_Matrix<Big_Fraction> m) 
    {
        const Big_FractionMatrix_converter converter = Big_FractionMatrix_converter();
        m.walk_in_optimized_order(converter);
        return converter.get_convertedMatrix();
    }

    /** Converter for {@link Field_Matrix}/{@link Big_Fraction}. */
    private static class Big_FractionMatrix_converter extends DefaultField_Matrix_Preserving_Visitor<Big_Fraction> 
    {
        /** Converted array. */
        private std::vector<std::vector<double>> data;
        /** Simple constructor. */
        Big_FractionMatrix_converter() 
        {
            super(Big_Fraction.ZERO);
        }

        /** {@inherit_doc} */
        //override
        public void start(const int& rows, int columns, int start_row, int end_row, int start_column, int end_column) 
        {
            data = std::vector<double>(rows][columns];
        }

        /** {@inherit_doc} */
        //override
        public void visit(const int& row, const int& column, Big_Fraction value) 
        {
            data[row][column] = value.double_value();
        }

        /**
         * Get the converted matrix.
         *
         * @return the converted matrix.
         */
        Array_2D_Row_Real_Matrix get_convertedMatrix() 
        {
            return Array_2D_Row_Real_Matrix(data, false);
        }
    }

    /** Serialize a {@link Real_Vector}.
     * <p>
     * This method is intended to be called from within a private
     * <code>write_object</code> method (after a call to
     * <code>oos.default_write_object()</code>) in a class that has a
     * {@link Real_Vector} field, which should be declared <code>transient</code>.
     * This way, the default handling does not serialize the vector (the {@link
     * Real_Vector} interface is not serializable by default) but this method does
     * serialize it specifically.
     * </p>
     * <p>
     * The following example shows how a simple class with a name and a real vector
     * should be written:
     * <pre><code>
     * class NamedVector  
     {
     *
     *     private const std::string name;
     *     private const Real_Vector coefficients;
     *
     *     // omitted constructors, getters ...
     *
     *     private void write_object(Object_Output_Stream oos) IOException 
     {
     *         oos.default_write_object();  // takes care of name field
     *         Matrix_Utils::serialize_real__vector(coefficients, oos);
     *     }
     *
     *     private void read_object(Object_Input_Stream ois) Class_Not_Found_Exception, IOException 
     {
     *         ois.default_read_object();  // takes care of name field
     *         Matrix_Utils::deserialize_real__vector(this, "coefficients", ois);
     *     }
     *
     * }
     * </code></pre>
     * </p>
     *
     * @param vector real vector to serialize
     * @param oos stream where the real vector should be written
     * @exception IOException if object cannot be written to stream
     * @see #deserialize_real__vector(Object, std::string, Object_Input_Stream)
     */
    public static void serialize_real__vector(const Real_Vector vector, const Object_Output_Stream oos)
        IOException 
        {
        const int n = vector.get_dimension();
        oos.writeInt(n);
        for (int i{}; i < n; ++i) 
        {
            oos.writeDouble(vector.get_entry(i));
        }
    }

    /** Deserialize  a {@link Real_Vector} field in a class.
     * <p>
     * This method is intended to be called from within a private
     * <code>read_object</code> method (after a call to
     * <code>ois.default_read_object()</code>) in a class that has a
     * {@link Real_Vector} field, which should be declared <code>transient</code>.
     * This way, the default handling does not deserialize the vector (the {@link
     * Real_Vector} interface is not serializable by default) but this method does
     * deserialize it specifically.
     * </p>
     * @param instance instance in which the field must be set up
     * @param fieldName name of the field within the class (may be private and const)
     * @param ois stream from which the real vector should be read
     * @exception Class_Not_Found_Exception if a class in the stream cannot be found
     * @exception IOException if object cannot be read from the stream
     * @see #serialize_real__vector(Real_Vector, Object_Output_Stream)
     */
    public static void deserialize_real__vector(const Object instance, const std::string fieldName, const Object_Input_Stream ois)
      Class_Not_Found_Exception, IOException 
      {
        try 
        {

            // read the vector data
            const int n = ois.readInt();
            const std::vector<double> data = std::vector<double>(n];
            for (int i{}; i < n; ++i) 
            {
                data[i] = ois.readDouble();
            }

            // create the instance
            const Real_Vector vector = Array_Real_Vector(data, false);

            // set up the field
            const java.lang.reflect.Field f =
                instance.get_class().get_declaredField(fieldName);
            f.set_accessible(true);
            f.set(instance, vector);

        }
        catch (NoSuchFieldException | Illegal_Access_Exception e) 
        {
            IOException ioe = IOException();
            ioe.initCause(e);
            throw ioe;
        }

    }

    /** Serialize a {@link Real_Matrix}.
     * <p>
     * This method is intended to be called from within a private
     * <code>write_object</code> method (after a call to
     * <code>oos.default_write_object()</code>) in a class that has a
     * {@link Real_Matrix} field, which should be declared <code>transient</code>.
     * This way, the default handling does not serialize the matrix (the {@link
     * Real_Matrix} interface is not serializable by default) but this method does
     * serialize it specifically.
     * </p>
     * <p>
     * The following example shows how a simple class with a name and a real matrix
     * should be written:
     * <pre><code>
     * class NamedMatrix  
     {
     *
     *     private const std::string name;
     *     private const Real_Matrix coefficients;
     *
     *     // omitted constructors, getters ...
     *
     *     private void write_object(Object_Output_Stream oos) IOException 
     {
     *         oos.default_write_object();  // takes care of name field
     *         Matrix_Utils::serialize_real__matrix(coefficients, oos);
     *     }
     *
     *     private void read_object(Object_Input_Stream ois) Class_Not_Found_Exception, IOException 
     {
     *         ois.default_read_object();  // takes care of name field
     *         Matrix_Utils::deserialize_real__matrix(this, "coefficients", ois);
     *     }
     *
     * }
     * </code></pre>
     * </p>
     *
     * @param matrix real matrix to serialize
     * @param oos stream where the real matrix should be written
     * @exception IOException if object cannot be written to stream
     * @see #deserialize_real__matrix(Object, std::string, Object_Input_Stream)
     */
    public static void serialize_real__matrix(const Real_Matrix matrix, const Object_Output_Stream oos)
        IOException 
        {
        const int n = matrix.get_row_dimension();
        const int m = matrix.get_column_dimension();
        oos.writeInt(n);
        oos.writeInt(m);
        for (int i{}; i < n; ++i) 
        {
            for (int j{}; j < m; ++j) 
            {
                oos.writeDouble(matrix.get_entry(i, j));
            }
        }
    }

    /** Deserialize  a {@link Real_Matrix} field in a class.
     * <p>
     * This method is intended to be called from within a private
     * <code>read_object</code> method (after a call to
     * <code>ois.default_read_object()</code>) in a class that has a
     * {@link Real_Matrix} field, which should be declared <code>transient</code>.
     * This way, the default handling does not deserialize the matrix (the {@link
     * Real_Matrix} interface is not serializable by default) but this method does
     * deserialize it specifically.
     * </p>
     * @param instance instance in which the field must be set up
     * @param fieldName name of the field within the class (may be private and const)
     * @param ois stream from which the real matrix should be read
     * @exception Class_Not_Found_Exception if a class in the stream cannot be found
     * @exception IOException if object cannot be read from the stream
     * @see #serialize_real__matrix(Real_Matrix, Object_Output_Stream)
     */
    public static void deserialize_real__matrix(const Object instance, const std::string fieldName, const Object_Input_Stream ois)
      Class_Not_Found_Exception, IOException 
      {
        try 
        {

            // read the matrix data
            const int n = ois.readInt();
            const int m = ois.readInt();
            const std::vector<std::vector<double>> data = std::vector<double>(n][m];
            for (int i{}; i < n; ++i) 
            {
                const std::vector<double> data_i = data[i];
                for (int j{}; j < m; ++j) 
                {
                    data_i[j] = ois.readDouble();
                }
            }

            // create the instance
            const Real_Matrix matrix = Array_2D_Row_Real_Matrix(data, false);

            // set up the field
            const java.lang.reflect.Field f =
                instance.get_class().get_declaredField(fieldName);
            f.set_accessible(true);
            f.set(instance, matrix);

        }
        catch (NoSuchFieldException | Illegal_Access_Exception e) 
        {
            IOException ioe = IOException();
            ioe.initCause(e);
            throw ioe;
        }
    }

    /**Solve  a  system of composed of a Lower Triangular Matrix
     * {@link Real_Matrix}.
     * <p>
     * This method is called to solve systems of equations which are
     * of the lower triangular form. The matrix {@link Real_Matrix}
     * is assumed, though not checked, to be in lower triangular form.
     * The vector {@link Real_Vector} is overwritten with the solution.
     * The matrix is checked that it is square and its dimensions match
     * the length of the vector.
     * </p>
     * @param rm Real_Matrix which is lower triangular
     * @param b  Real_Vector this is overwritten
     * @ if the matrix and vector are not
     * conformable
     * @ if the matrix {@code rm} is not square
     * @Math_Runtime_Exception if the absolute value of one of the diagonal
     * coefficient of {@code rm} is lower than {@link Precision#SAFE_MIN}
     */
    public static void solve_uower_triangular_system(Real_Matrix rm, Real_Vector b)
        , Math_Runtime_Exception 
        {
        if ((rm == NULL) || (b == NULL) || ( rm.get_row_dimension() != b.get_dimension())) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, (rm == NULL) ? 0 : rm.get_row_dimension(), (b  == NULL) ? 0 : b.get_dimension());
        }
        if( rm.get_column_dimension() != rm.get_row_dimension() )
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NON_SQUARE_MATRIX, rm.get_row_dimension(), rm.get_column_dimension());
        }
        int rows = rm.get_row_dimension();
        for( int i{}; i < rows ; i++ )
        {
            double diag = rm.get_entry(i, i);
            if( std::abs(diag) < Precision.SAFE_MIN )
            {
                throw std::exception("not implemented");
                //throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::ZERO_DENOMINATOR);
            }
            double bi = b.get_entry(i)/diag;
            b.set_entry(i,  bi );
            for( int j = i+1; j< rows; j++ )
            {
                b.set_entry(j, b.get_entry(j)-bi*rm.get_entry(j,i)  );
            }
        }
    }

    /** Solver a  system composed  of an Upper Triangular Matrix
     * {@link Real_Matrix}.
     * <p>
     * This method is called to solve systems of equations which are
     * of the lower triangular form. The matrix {@link Real_Matrix}
     * is assumed, though not checked, to be in upper triangular form.
     * The vector {@link Real_Vector} is overwritten with the solution.
     * The matrix is checked that it is square and its dimensions match
     * the length of the vector.
     * </p>
     * @param rm Real_Matrix which is upper triangular
     * @param b  Real_Vector this is overwritten
     * @ if the matrix and vector are not
     * conformable
     * @ if the matrix {@code rm} is not
     * square
     * @Math_Runtime_Exception if the absolute value of one of the diagonal
     * coefficient of {@code rm} is lower than {@link Precision#SAFE_MIN}
     */
    public static void solve_upper_triangular_system(Real_Matrix rm, Real_Vector b)
        , Math_Runtime_Exception 
        {
        if ((rm == NULL) || (b == NULL) || ( rm.get_row_dimension() != b.get_dimension())) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, (rm == NULL) ? 0 : rm.get_row_dimension(), (b  == NULL) ? 0 : b.get_dimension());
        }
        if( rm.get_column_dimension() != rm.get_row_dimension() )
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NON_SQUARE_MATRIX, rm.get_row_dimension(), rm.get_column_dimension());
        }
        int rows = rm.get_row_dimension();
        for( int i = rows-1 ; i >-1 ; i-- )
        {
            double diag = rm.get_entry(i, i);
            if( std::abs(diag) < Precision.SAFE_MIN )
            {
                throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::ZERO_DENOMINATOR);
            }
            double bi = b.get_entry(i)/diag;
            b.set_entry(i,  bi );
            for( int j = i-1; j>-1; j-- )
            {
                b.set_entry(j, b.get_entry(j)-bi*rm.get_entry(j,i)  );
            }
        }
    }

    /**
     * Computes the inverse of the given matrix by splitting it into
     * 4 sub-matrices.
     *
     * @param m Matrix whose inverse must be computed.
     * @param splitIndex Index that determines the "split" line and
     * column.
     * The element corresponding to this index will part of the
     * upper-left sub-matrix.
     * @return the inverse of {@code m}.
     * @ if {@code m} is not square.
     */
    public static Real_Matrix blockInverse(Real_Matrix m, int splitIndex) 
    {
        const int n = m.get_row_dimension();
        if (m.get_column_dimension() != n) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NON_SQUARE_MATRIX, m.get_row_dimension(), m.get_column_dimension());
        }

        const int splitIndex1 = splitIndex + 1;

        const Real_Matrix& a = m.get_sub_matrix(0, splitIndex, 0, splitIndex);
        const Real_Matrix b = m.get_sub_matrix(0, splitIndex, splitIndex1, n - 1);
        const Real_Matrix c = m.get_sub_matrix(splitIndex1, n - 1, 0, splitIndex);
        const Real_Matrix d = m.get_sub_matrix(splitIndex1, n - 1, splitIndex1, n - 1);

        const Singular_Value_Decomposition aDec = Singular_Value_Decomposition(a);
        const Decomposition_Solver aSolver = aDec.get_solver();
        if (!aSolver.is_non_singular()) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::SINGULAR_MATRIX);
        }
        const Real_Matrix& aInv = aSolver.get_inverse();

        const Singular_Value_Decomposition dDec = Singular_Value_Decomposition(d);
        const Decomposition_Solver dSolver = dDec.get_solver();
        if (!dSolver.is_non_singular()) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::SINGULAR_MATRIX);
        }
        const Real_Matrix dInv = dSolver.get_inverse();

        const Real_Matrix tmp1 = a.subtract(b.multiply(dInv).multiply(c));
        const Singular_Value_Decomposition tmp1Dec = Singular_Value_Decomposition(tmp1);
        const Decomposition_Solver tmp1Solver = tmp1Dec.get_solver();
        if (!tmp1Solver.is_non_singular()) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::SINGULAR_MATRIX);
        }
        const Real_Matrix result00 = tmp1Solver.get_inverse();

        const Real_Matrix tmp2 = d.subtract(c.multiply(aInv).multiply(b));
        const Singular_Value_Decomposition tmp2_dec = Singular_Value_Decomposition(tmp2);
        const Decomposition_Solver tmp2Solver = tmp2_dec.get_solver();
        if (!tmp2Solver.is_non_singular()) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::SINGULAR_MATRIX);
        }
        const Real_Matrix result11 = tmp2Solver.get_inverse();

        const Real_Matrix result01 = aInv.multiply(b).multiply(result11).scalar_multiply(-1);
        const Real_Matrix result10 = dInv.multiply(c).multiply(result00).scalar_multiply(-1);

        const Real_Matrix result = Array_2D_Row_Real_Matrix(n, n);
        result.set_sub_matrix(result00.get_data(), 0, 0);
        result.set_sub_matrix(result01.get_data(), 0, splitIndex1);
        result.set_sub_matrix(result10.get_data(), splitIndex1, 0);
        result.set_sub_matrix(result11.get_data(), splitIndex1, splitIndex1);

        return result;
    }

    /**
     * Computes the inverse of the given matrix.
     * <p>
     * By default, the inverse of the matrix is computed using the QR-decomposition, * unless a more efficient method can be determined for the input matrix.
     * <p>
     * Note: this method will use a singularity threshold of 0, * use {@link #inverse(Real_Matrix, double)} if a different threshold is needed.
     *
     * @param matrix Matrix whose inverse shall be computed
     * @return the inverse of {@code matrix}
     * @Null_Argument_Exception if {@code matrix} is {@code NULL}
     * @ if m is singular
     * @ if matrix is not square
     */
    public static Real_Matrix inverse(Real_Matrix matrix)
            , Null_Argument_Exception 
            {
        return inverse(matrix, 0);
    }

    /**
     * Computes the inverse of the given matrix.
     * <p>
     * By default, the inverse of the matrix is computed using the QR-decomposition, * unless a more efficient method can be determined for the input matrix.
     *
     * @param matrix Matrix whose inverse shall be computed
     * @param threshold Singularity threshold
     * @return the inverse of {@code m}
     * @Null_Argument_Exception if {@code matrix} is {@code NULL}
     * @ if matrix is singular
     * @ if matrix is not square
     */
    public static Real_Matrix inverse(Real_Matrix matrix, double threshold)
            , Null_Argument_Exception 
            {

        //Math_Utils::check_not_null(matrix);

        if (!matrix.is_square()) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NON_SQUARE_MATRIX, matrix.get_row_dimension(), matrix.get_column_dimension());
        }

        if (matrix instanceof Diagonal_Matrix) 
        {
            throw std::exception("not implemented");
            //return ((Diagonal_Matrix) matrix).inverse(threshold);
        }
else 
        {
            QR_Decomposition decomposition = QR_Decomposition(matrix, threshold);
            return decomposition.get_solver().get_inverse();
        }
    }

    /**
     * Computes the <a href="https://mathworld.wolfram.com/MatrixExponential.html">
     * matrix exponential</a> of the given matrix.
     *
     * The algorithm implementation follows the Pade approximant method of
     * <p>Higham, Nicholas J. \xe2\x80\x9cThe Scaling and Squaring Method for the Matrix Exponential
     * Revisited.\xe2\x80\x9d SIAM Journal on Matrix Analysis and Applications 26, no. 4 (January 2005): 1179\xe2\x80\x9393.</p>
     *
     * @param rm Real_Matrix whose inverse shall be computed
     * @return The inverse of {@code rm}
     * @ if matrix is not square
     */
    public static Real_Matrix matrixExponential(const Real_Matrix rm) 
    {

        // Check that the input matrix is square
        if (!rm.is_square()) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NON_SQUARE_MATRIX, rm.get_row_dimension(), rm.get_column_dimension());
        }

        // Preprocessing to reduce the norm
        int dim = rm.get_row_dimension();
        const Real_Matrix identity = Matrix_Utils::create_real_identity_matrix(dim);
        const double preprocessScale = rm.get_trace() / dim;
        Real_Matrix scaledMatrix = rm.copy();
        scaledMatrix = scaledMatrix.subtract(identity.scalar_multiply(preprocessScale));

        // Select pade degree required
        const double l1Norm = rm.get_norm1();
        std::vector<double> padeCoefficients;
        int squaringCount = 0;

        if (l1Norm < 1.495585217958292e-2) 
        {
            padeCoefficients = PADE_COEFFICIENTS_3;
        }
else if (l1Norm < 2.539398330063230e-1) 
        {
            padeCoefficients = PADE_COEFFICIENTS_5;
        }
else if (l1Norm < 9.504178996162932e-1) 
        {
            padeCoefficients = PADE_COEFFICIENTS_7;
        }
else if (l1Norm < 2.097847961257068) 
        {
            padeCoefficients = PADE_COEFFICIENTS_9;
        }
else 
        {
            padeCoefficients = PADE_COEFFICIENTS_13;

            // Calculate scaling factor
            const double normScale = 5.371920351148152;
            squaringCount = Math.max(0, Math.get_exponent(l1Norm / normScale));

            // Scale matrix by power of 2
            const int const_squaring_count = squaringCount;
            scaledMatrix.walk_in_optimized_order(new Default_Real_Matrix_Changing_Visitor() 
            {
                //override
                public double visit(const int& row, const int& column, double value) 
                {
                    return Math.scalb(value, -const_squaring_count);
                }
            });
        }

        // Calculate U and V using Horner
        // See Golub, Gene H., and Charles F. van Loan. Matrix Computations. 4th ed.
        // John Hopkins University Press, 2013.  pages 530/531
        const Real_Matrix scaledMatrix2 = scaledMatrix.multiply(scaledMatrix);
        const int coeffLength = padeCoefficients.size();

        // Calculate V
        Real_Matrix padeV = Matrix_Utils::create_real_matrix(dim, dim);
        for (int i = coeffLength - 1; i > 1; i -= 2) 
        {
            padeV = scaledMatrix2.multiply(padeV.add(identity.scalar_multiply(padeCoefficients[i])));
        }
        padeV = scaledMatrix.multiply(padeV.add(identity.scalar_multiply(padeCoefficients[1])));

        // Calculate U
        Real_Matrix padeU = Matrix_Utils::create_real_matrix(dim, dim);
        for (int i = coeffLength - 2; i > 1; i -= 2) 
        {
            padeU = scaledMatrix2.multiply(padeU.add(identity.scalar_multiply(padeCoefficients[i])));
        }
        padeU = padeU.add(identity.scalar_multiply(padeCoefficients[0]));

        // Calculate pade approximate by solving (U-V) F = (U+V) for F
        Real_Matrix padeNumer = padeU.add(padeV);
        Real_Matrix padeDenom = padeU.subtract(padeV);

        // Calculate the matrix ratio
        QR_Decomposition decomposition = QR_Decomposition(padeDenom);
        Real_Matrix result = decomposition.get_solver().solve(padeNumer);

        // Repeated squaring if matrix was scaled
        for (int i{}; i < squaringCount; i++) 
        {
            result = result.multiply(result);
        }

        // Undo preprocessing
        result = result.scalar_multiply(Math.exp(preprocessScale));

        return result;
    }

    /** Orthonormalize a list of vectors.
     * <p>
     * Orthonormalization is performed by using the Modified Gram-Schmidt process.
     * </p>
     * @param independent list of independent vectors
     * @param threshold projected vectors with a norm less than or equal to this threshold
     * are considered to have zero norm, hence the vectors they come from are not independent from
     * previous vectors
     * @param handler handler for dependent vectors
     * @return orthonormal basis having the same span as {@code independent}
     * @since 2.1
     */
    public static List<Real_Vector> orthonormalize(const List<Real_Vector> independent, const double threshold, const Dependent_Vectors_Handler handler) 
    {

        // create separate list
        const List<Real_Vector> basis = Array_list<>(independent);

        // loop over basis vectors
        int index = 0;
        while (index < basis.size()) 
        {

            // check dependency
            const Real_Vector vi = basis.get(index);
            const double norm = vi.get_norm();
            if (norm <= threshold) 
            {
                // the current vector is dependent from the previous ones
                index = handler.manage_dependent(index, basis);
            }
else 
            {

                // normalize basis vector in place
                vi.map_divide_to_self(vi.get_norm());

                // project remaining vectors in place
                for (int j = index + 1; j < basis.size(); ++j) 
                {
                    const Real_Vector vj  = basis.get(j);
                    const double     dot = vi.dot_product(vj);
                    for (int k{}; k < vj.get_dimension(); ++k) 
                    {
                        vj.set_entry(k, vj.get_entry(k) - dot * vi.get_entry(k));
                    }
                }

                ++index;

            }

        }

        return basis;

    }

    /** Orthonormalize a list of vectors.
     * <p>
     * Orthonormalization is performed by using the Modified Gram-Schmidt process.
     * </p>
     * @param <T> type of the field elements
     * @param independent list of independent vectors
     * @param threshold projected vectors with a norm less than or equal to this threshold
     * are considered to have zero norm, hence the vectors they come from are not independent from
     * previous vectors
     * @param field type of the files elements
     * @param handler handler for dependent vectors
     * @return orthonormal basis having the same span as {@code independent}
     * @since 2.1
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    public static List<Field_Vector<T>> orthonormalize(const Field<T> field, const List<Field_Vector<T>> independent, const T threshold, const Dependent_Vectors_Handler handler) 
    {
        // create separate list
        const List<Field_Vector<T>> basis = Array_list<>(independent);

        // loop over basis vectors
        int index = 0;
        while (index < basis.size()) 
        {

            // check dependency
            const Field_Vector<T> vi = basis.get(index);
            const T norm = vi.dot_product(vi).sqrt();
            if (norm.subtract(threshold).get_real() <= 0) 
            {
                // the current vector is dependent from the previous ones
                index = handler.manage_dependent(field, index, basis);
            }
else 
            {

                // normalize basis vector in place
                vi.map_divide_to_self(norm);

                // project remaining vectors in place
                for (int j = index + 1; j < basis.size(); ++j) 
                {
                    const Field_Vector<T> vj  = basis.get(j);
                    const T              dot = vi.dot_product(vj);
                    for (int k{}; k < vj.get_dimension(); ++k) 
                    {
                        vj.set_entry(k, vj.get_entry(k).subtract(dot.multiply(vi.get_entry(k))));
                    }
                }

                ++index;

            }
        }

        return basis;

    }

}


