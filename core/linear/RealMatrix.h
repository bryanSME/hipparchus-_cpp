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

//import org.hipparchus.analysis.Univariate_Function;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Null_Argument_Exception;
//import org.hipparchus.util.FastMath;

#include "AnyMatrix.h"
#include "RealVector.h"
#include "RealMatrixPreservingVisitor.h"
#include "RealMatrixChangingVisitor.h"
#include "../analysis/UnivariateFunction.h"

#include <vector>

/**
 * Interface defining a real-valued matrix with basic algebraic operations.
 * <p>
 * Matrix element indexing is 0-based -- e.g., <code>get_entry(0, 0)</code>
 * returns the element in the first row, first column of the matrix.</p>
 *
 */
class Real_Matrix : Any_Matrix 
{

    /**
     * Create a Real_Matrix of the same type as the instance with the
     * supplied
     * row and column dimensions.
     *
     * @param row_dimension the number of rows in the matrix
     * @param column_dimension the number of columns in the matrix
     * @return a matrix of the same type as the instance
     * @ if row or column dimension is not
     * positive.
     */
    Real_Matrix create_matrix(const int& row_dimension, const int& column_dimension);

    /**
     * Returns a (deep) copy of this.
     *
     * @return matrix copy
     */
    virtual Real_Matrix copy() = 0;

    /**
     * Returns the sum of {@code this} and {@code m}.
     *
     * @param m matrix to be added
     * @return {@code this + m}
     * @ if {@code m} is not the same
     * size as {@code this}.
     */
    virtual Real_Matrix add(Real_Matrix m) = 0;

    /**
     * Returns {@code this} minus {@code m}.
     *
     * @param m matrix to be subtracted
     * @return {@code this - m}
     * @ if {@code m} is not the same
     * size as {@code this}.
     */
    virtual Real_Matrix subtract(Real_Matrix m) = 0;

    /**
     * Returns the result of adding {@code d} to each entry of {@code this}.
     *
     * @param d value to be added to each entry
     * @return {@code d + this}
     */
    virtual Real_Matrix scalar_add(double d) = 0;

    /**
     * Returns the result of multiplying each entry of {@code this} by
     * {@code d}.
     *
     * @param d value to multiply all entries by
     * @return {@code d * this}
     */
    virtual Real_Matrix scalar_multiply(double d) = 0;

    /**
     * Returns the result of postmultiplying {@code this} by {@code m}.
     *
     * @param m matrix to postmultiply by
     * @return {@code this * m}
     * @ if
     * {@code column_dimension(this) != row_dimension(m)}
     */
    virtual Real_Matrix multiply(Real_Matrix m) = 0;

    /**
     * Returns the result of postmultiplying {@code this} by {@code m^T}.
     * <p>
     * This is equivalent to call {@link #multiply(Real_Matrix) multiply}(m.{@link #transpose()}), * but some implementations may avoid building the intermediate transposed matrix.
     * </p>
     * @param m matrix to first transpose and second postmultiply by
     * @return {@code this * m^T}
     * @ if
     * {@code column_dimension(this) != column_dimension(m)}
     * @since 1.3
     */
    Real_Matrix multiply_transposed(const Real_Matrix& m)
    {
        return multiply(m.transpose());
    }

    /**
     * Returns the result of postmultiplying {@code this^T} by {@code m}.
     * <p>
     * This is equivalent to call {@link #transpose()}.{@link #multiply(Real_Matrix) multiply(m)}, * but some implementations may avoid building the intermediate transposed matrix.
     * </p>
     * @param m matrix to postmultiply by
     * @return {@code this^T * m}
     * @ if
     * {@code column_dimension(this) != column_dimension(m)}
     * @since 1.3
     */
    Real_Matrix transpose_multiply(const Real_Matrix& m)
    {
        return transpose().multiply(m);
    }

    /**
     * Returns the result of premultiplying {@code this} by {@code m}.
     *
     * @param m matrix to premultiply by
     * @return {@code m * this}
     * @ if
     * {@code row_dimension(this) != column_dimension(m)}
     */
    virtual Real_Matrix pre_multiply(Real_Matrix m) = 0;

    /**
     * Returns the result of multiplying {@code this} with itself {@code p}
     * times. Depending on the underlying storage, instability for high powers
     * might occur.
     *
     * @param p raise {@code this} to power {@code p}
     * @return {@code this^p}
     * @ if {@code p < 0}
     * @ if the matrix is not square
     */
    virtual Real_Matrix power(const int& p) = 0;

    /**
     * Returns matrix entries as a two-dimensional array.
     *
     * @return 2-dimensional array of entries
     */
    virtual std::vector<std::vector<double>> get_data() = 0;

    /**
     * Returns the <a href="http://mathworld.wolfram.com/MaximumAbsoluteColumn_sumNorm.html">
     * maximum absolute column sum norm</a> (L<sub>1</sub>) of the matrix.
     *
     * @return norm
     */
    double get_norm1() 
    {
        return walk_in_column_order(Real_Matrix_Preserving_Visitor() 
        {
            /** Last row index. */
            private int end_row;

            /** Sum of absolute values on one column. */
            private double column_sum;

            /** Maximal sum across all columns. */
            private double max_col_sum;

            /** {@inherit_doc} */
            //override
            public void start(const int rows, const int columns, const int& start_row, const int& end_row, const int& start_column, const int& end_column) 
            {
                this.end_row = end_row;
                column_sum   = 0;
                max_col_sum   = 0;
            }

            /** {@inherit_doc} */
            //override
            public void visit(const int& row, const int& column, const double value) 
            {
                column_sum += std::abs(value);
                if (row == end_row) 
                {
                    max_col_sum = std::max(max_col_sum, column_sum);
                    column_sum = 0;
                }
            }

            /** {@inherit_doc} */
            //override
            public double end() 
            {
                return max_col_sum;
            }
        });
    }

    /**
     * Returns the <a href="http://mathworld.wolfram.com/MaximumAbsoluteRowSumNorm.html">
     * maximum absolute row sum norm</a> (L<sub>&infin;</sub>) of the matrix.
     *
     * @return norm
     */
    double get_norm_infty() 
    {
        return walk_in_row_order(
            Real_Matrix_Preserving_Visitor() 
            {
            private:
                /** Last column index. */
                int my_end_column;

                /** Sum of absolute values on one row. */
                double my_row_sum;

                /** Maximal sum across all rows. */
                double my_max_row_sum;

            public:
                /** {@inherit_doc} */
                //override
                void start(const int& rows, const int& columns, const int& start_row, const int& end_row, const int& start_column, const int& end_column) 
                {
                    my_end_column = end_column;
                    row_sum = 0;
                    max_row_sum = 0;
                }

                /** {@inherit_doc} */
                //override
                void visit(const int& row, const int& column, const double& value) 
                {
                    my_row_sum += std::abs(value);
                    if (column == my_end_column) 
                    {
                        max_row_sum = std::max(max_row_sum, my_row_sum);
                        my_row_sum = 0;
                    }
                }

                /** {@inherit_doc} */
                //override
                double end() const
                {
                    return max_row_sum;
                }
            }
        );
    }

    /**
     * Returns the <a href="http://mathworld.wolfram.com/FrobeniusNorm.html">
     * Frobenius norm</a> of the matrix.
     *
     * @return norm
     */
    virtual  double get_frobenius_norm() = 0;

    /**
     * Gets a submatrix. Rows and columns are indicated
     * counting from 0 to n-1.
     *
     * @param start_row Initial row index
     * @param end_row Final row index (inclusive)
     * @param start_column Initial column index
     * @param end_column Final column index (inclusive)
     * @return The sub_matrix containing the data of the
     * specified rows and columns.
     * @ if the indices are not valid.
     * @ if {@code end_row < start_row} or
     * {@code end_column < start_column}.
     */
    virtual Real_Matrix get_sub_matrix(const int& start_row, const int& end_row, const int& start_column, const int& end_column) = 0;

    /**
     * Gets a submatrix. Rows and columns are indicated counting from 0 to n-1.
     *
     * @param selected_rows Array of row indices.
     * @param selected_columns Array of column indices.
     * @return The sub_matrix containing the data in the specified rows and
     * columns
     * @Null_Argument_Exception if the row or column selections are
     * {@code NULL}
     * @ if the row or column selections are empty (zero
     * length).
     * @ if the indices are not valid.
     */
    virtual Real_Matrix get_sub_matrix(const std::vector<int>& selected_rows, const std::vector<int>& selected_columns) = 0;

    /**
     * Copy a submatrix. Rows and columns are indicated counting from 0 to n-1.
     *
     * @param start_row Initial row index
     * @param end_row Final row index (inclusive)
     * @param start_column Initial column index
     * @param end_column Final column index (inclusive)
     * @param destination The arrays where the submatrix data should be copied
     * (if larger than rows/columns counts, only the upper-left part will be
     * used)
     * @ if the indices are not valid.
     * @ if {@code end_row < start_row} or
     * {@code end_column < start_column}.
     * @ if the destination array is too
     * small.
     */
    virtual void copy_sub_matrix(const int& start_row, const int& end_row, const int& start_column, const int& end_column, std::vector<std::vector<double>>& destination) = 0;

    /**
     * Copy a submatrix. Rows and columns are indicated counting from 0 to n-1.
     *
     * @param selected_rows Array of row indices.
     * @param selected_columns Array of column indices.
     * @param destination The arrays where the submatrix data should be copied
     * (if larger than rows/columns counts, only the upper-left part will be
     * used)
     * @Null_Argument_Exception if the row or column selections are
     * {@code NULL}
     * @ if the row or column selections are empty (zero
     * length).
     * @ if the indices are not valid.
     * @ if the destination array is too
     * small.
     */
    virtual void copy_sub_matrix(std::vector<int> selected_rows, std::vector<int> selected_columns, std::vector<std::vector<double>> destination) = 0;

   /**
    * Replace the submatrix starting at {@code row, column} using data in the
    * input {@code sub_matrix} array. Indexes are 0-based.
    * <p>
    * Example:<br>
    * Starting with <pre>
    * 1  2  3  4
    * 5  6  7  8
    * 9  0  1  2
    * </pre>
    * and <code>sub_matrix = {{3, 4} {5,6}}</code>, invoking
    * {@code set_sub_matrix(sub_matrix,1,1))} will result in <pre>
    * 1  2  3  4
    * 5  3  4  8
    * 9  5  6  2
    * </pre></p>
    *
    * @param sub_matrix  array containing the submatrix replacement data
    * @param row  row coordinate of the top, left element to be replaced
    * @param column  column coordinate of the top, left element to be replaced
    * @ if {@code sub_matrix} is empty.
    * @ if {@code sub_matrix} does not fit into
    * this matrix from element in {@code (row, column)}.
    * @ if {@code sub_matrix} is not rectangular
    * (not all rows have the same length) or empty.
    * @Null_Argument_Exception if {@code sub_matrix} is {@code NULL}.
    */
    virtual void set_sub_matrix(std::vector<std::vector<double>> sub_matrix, int row, int column) = 0;

   /**
    * Get the entries at the given row index as a row matrix.  Row indices start
    * at 0.
    *
    * @param row Row to be fetched.
    * @return row Matrix.
    * @ if the specified row index is invalid.
    */
    virtual Real_Matrix get_row_matrix(const int& row) = 0;

    /**
     * Sets the specified {@code row} of {@code this} matrix to the entries of
     * the specified row {@code matrix}. Row indices start at 0.
     *
     * @param row Row to be set.
     * @param matrix Row matrix to be copied (must have one row and the same
     * number of columns as the instance).
     * @ if the specified row index is invalid.
     * @ if the row dimension of the
     * {@code matrix} is not {@code 1}, or the column dimensions of {@code this}
     * and {@code matrix} do not match.
     */
    virtual void set_row_matrix(const int& row, Real_Matrix matrix) = 0;

    /**
     * Get the entries at the given column index as a column matrix. Column
     * indices start at 0.
     *
     * @param column Column to be fetched.
     * @return column Matrix.
     * @ if the specified column index is invalid.
     */
    virtual Real_Matrix get_column_matrix(const int& column) = 0;

    /**
     * Sets the specified {@code column} of {@code this} matrix to the entries
     * of the specified column {@code matrix}. Column indices start at 0.
     *
     * @param column Column to be set.
     * @param matrix Column matrix to be copied (must have one column and the
     * same number of rows as the instance).
     * @ if the specified column index is invalid.
     * @ if the column dimension of the
     * {@code matrix} is not {@code 1}, or the row dimensions of {@code this}
     * and {@code matrix} do not match.
     */
    virtual void set_column_matrix(const int& column, Real_Matrix matrix) = 0;

    /**
     * Returns the entries in row number {@code row} as a vector. Row indices
     * start at 0.
     *
     * @param row Row to be fetched.
     * @return a row vector.
     * @ if the specified row index is invalid.
     */
    virtual Real_Vector get_row_vector(const int& row) = 0;

    /**
     * Sets the specified {@code row} of {@code this} matrix to the entries of
     * the specified {@code vector}. Row indices start at 0.
     *
     * @param row Row to be set.
     * @param vector row vector to be copied (must have the same number of
     * column as the instance).
     * @ if the specified row index is invalid.
     * @ if the {@code vector} dimension
     * does not match the column dimension of {@code this} matrix.
     */
    virtual void set_row_vector(const int& row, Real_Vector& vector) = 0;

    /**
     * Get the entries at the given column index as a vector. Column indices
     * start at 0.
     *
     * @param column Column to be fetched.
     * @return a column vector.
     * @ if the specified column index is invalid
     */
    virtual Real_Vector get_column_vector(const int& column) = 0;

    /**
     * Sets the specified {@code column} of {@code this} matrix to the entries
     * of the specified {@code vector}. Column indices start at 0.
     *
     * @param column Column to be set.
     * @param vector column vector to be copied (must have the same number of
     * rows as the instance).
     * @ if the specified column index is invalid.
     * @ if the {@code vector} dimension
     * does not match the row dimension of {@code this} matrix.
     */
    virtual void set_column_vector(const int& column, Real_Vector& vector) = 0;

    /**
     * Get the entries at the given row index. Row indices start at 0.
     *
     * @param row Row to be fetched.
     * @return the array of entries in the row.
     * @ if the specified row index is not valid.
     */
    virtual std::vector<double> get_row(const int& row) = 0;

    /**
     * Sets the specified {@code row} of {@code this} matrix to the entries
     * of the specified {@code array}. Row indices start at 0.
     *
     * @param row Row to be set.
     * @param array Row matrix to be copied (must have the same number of
     * columns as the instance)
     * @ if the specified row index is invalid.
     * @ if the {@code array} length does
     * not match the column dimension of {@code this} matrix.
     */
    virtual void set_row(const int& row, std::vector<double> array) = 0;

    /**
     * Get the entries at the given column index as an array. Column indices
     * start at 0.
     *
     * @param column Column to be fetched.
     * @return the array of entries in the column.
     * @ if the specified column index is not valid.
     */
    virtual std::vector<double> get_column(const int& column) = 0;

    /**
     * Sets the specified {@code column} of {@code this} matrix to the entries
     * of the specified {@code array}. Column indices start at 0.
     *
     * @param column Column to be set.
     * @param array Column array to be copied (must have the same number of
     * rows as the instance).
     * @ if the specified column index is invalid.
     * @ if the {@code array} length does
     * not match the row dimension of {@code this} matrix.
     */
    virtual void set_column(const int& column, std::vector<double> array) = 0;

    /**
     * Get the entry in the specified row and column. Row and column indices
     * start at 0.
     *
     * @param row Row index of entry to be fetched.
     * @param column Column index of entry to be fetched.
     * @return the matrix entry at {@code (row, column)}.
     * @ if the row or column index is not valid.
     */
    virtual double get_entry(const int& row, const int& column) = 0;

    /**
     * Set the entry in the specified row and column. Row and column indices
     * start at 0.
     *
     * @param row Row index of entry to be set.
     * @param column Column index of entry to be set.
     * @param value the value of the entry.
     * @ if the row or column index is not valid
     */
    virtual void set_entry(const int& row, const int& column, double value) = 0;

    /**
     * Adds (in place) the specified value to the specified entry of
     * {@code this} matrix. Row and column indices start at 0.
     *
     * @param row Row index of the entry to be modified.
     * @param column Column index of the entry to be modified.
     * @param increment value to add to the matrix entry.
     * @ if the row or column index is not valid.
     */
    virtual void add_to_entry(const int& row, const int& column, double increment) = 0;

    /**
     * Multiplies (in place) the specified entry of {@code this} matrix by the
     * specified value. Row and column indices start at 0.
     *
     * @param row Row index of the entry to be modified.
     * @param column Column index of the entry to be modified.
     * @param factor Multiplication factor for the matrix entry.
     * @ if the row or column index is not valid.
     */
    virtual void multiply_entry(const int& row, const int& column, double factor) = 0;

    /**
     * Returns the transpose of this matrix.
     *
     * @return transpose matrix
     */
    virtual Real_Matrix transpose() = 0;

    /**
     * Returns the <a href="http://mathworld.wolfram.com/MatrixTrace.html">
     * trace</a> of the matrix (the sum of the elements on the main diagonal).
     *
     * @return the trace.
     * @ if the matrix is not square.
     */
    virtual double get_trace() = 0;

    /**
     * Returns the result of multiplying this by the vector {@code v}.
     *
     * @param v the vector to operate on
     * @return {@code this * v}
     * @ if the length of {@code v} does not
     * match the column dimension of {@code this}.
     */
    virtual std::vector<double> operate(std::vector<double> v) = 0;

    /**
     * Returns the result of multiplying this by the vector {@code v}.
     *
     * @param v the vector to operate on
     * @return {@code this * v}
     * @ if the dimension of {@code v} does not
     * match the column dimension of {@code this}.
     */
    virtual Real_Vector operate(Real_Vector v) = 0;

    /**
     * Returns the (row) vector result of premultiplying this by the vector {@code v}.
     *
     * @param v the row vector to premultiply by
     * @return {@code v * this}
     * @ if the length of {@code v} does not
     * match the row dimension of {@code this}.
     */
    virtual std::vector<double> pre_multiply(std::vector<double> v) = 0;

    /**
     * Returns the (row) vector result of premultiplying this by the vector {@code v}.
     *
     * @param v the row vector to premultiply by
     * @return {@code v * this}
     * @ if the dimension of {@code v} does not
     * match the row dimension of {@code this}.
     */
    virtual Real_Vector pre_multiply(Real_Vector v) = 0;

    /**
     * Visit (and possibly change) all matrix entries in row order.
     * <p>Row order starts at upper left and iterating through all elements
     * of a row from left to right before going to the leftmost element
     * of the next row.</p>
     * @param visitor visitor used to process all matrix entries
     * @see #walk_in_row_order(Real_Matrix_Preserving_Visitor)
     * @see #walk_in_row_order(Real_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_row_order(Real_Matrix_Preserving_Visitor, int, int, int, int)
     * @see #walk_in_column_order(Real_Matrix_Changing_Visitor)
     * @see #walk_in_column_order(Real_Matrix_Preserving_Visitor)
     * @see #walk_in_column_order(Real_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_column_order(Real_Matrix_Preserving_Visitor, int, int, int, int)
     * @see #walk_in_optimized_order(Real_Matrix_Changing_Visitor)
     * @see #walk_in_optimized_order(Real_Matrix_Preserving_Visitor)
     * @see #walk_in_optimized_order(Real_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_optimized_order(Real_Matrix_Preserving_Visitor, int, int, int, int)
     * @return the value returned by {@link Real_Matrix_Changing_Visitor#end()} at the end
     * of the walk
     */
    virtual double walk_in_row_order(Real_Matrix_Changing_Visitor visitor) = 0;

    /**
     * Visit (but don't change) all matrix entries in row order.
     * <p>Row order starts at upper left and iterating through all elements
     * of a row from left to right before going to the leftmost element
     * of the next row.</p>
     * @param visitor visitor used to process all matrix entries
     * @see #walk_in_row_order(Real_Matrix_Changing_Visitor)
     * @see #walk_in_row_order(Real_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_row_order(Real_Matrix_Preserving_Visitor, int, int, int, int)
     * @see #walk_in_column_order(Real_Matrix_Changing_Visitor)
     * @see #walk_in_column_order(Real_Matrix_Preserving_Visitor)
     * @see #walk_in_column_order(Real_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_column_order(Real_Matrix_Preserving_Visitor, int, int, int, int)
     * @see #walk_in_optimized_order(Real_Matrix_Changing_Visitor)
     * @see #walk_in_optimized_order(Real_Matrix_Preserving_Visitor)
     * @see #walk_in_optimized_order(Real_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_optimized_order(Real_Matrix_Preserving_Visitor, int, int, int, int)
     * @return the value returned by {@link Real_Matrix_Preserving_Visitor#end()} at the end
     * of the walk
     */
    virtual double walk_in_row_order(Real_Matrix_Preserving_Visitor visitor) = 0;

    /**
     * Visit (and possibly change) some matrix entries in row order.
     * <p>Row order starts at upper left and iterating through all elements
     * of a row from left to right before going to the leftmost element
     * of the next row.</p>
     * @param visitor visitor used to process all matrix entries
     * @param start_row Initial row index
     * @param end_row Final row index (inclusive)
     * @param start_column Initial column index
     * @param end_column Final column index
     * @ if the indices are not valid.
     * @ if {@code end_row < start_row} or
     * {@code end_column < start_column}.
     * @see #walk_in_row_order(Real_Matrix_Changing_Visitor)
     * @see #walk_in_row_order(Real_Matrix_Preserving_Visitor)
     * @see #walk_in_row_order(Real_Matrix_Preserving_Visitor, int, int, int, int)
     * @see #walk_in_column_order(Real_Matrix_Changing_Visitor)
     * @see #walk_in_column_order(Real_Matrix_Preserving_Visitor)
     * @see #walk_in_column_order(Real_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_column_order(Real_Matrix_Preserving_Visitor, int, int, int, int)
     * @see #walk_in_optimized_order(Real_Matrix_Changing_Visitor)
     * @see #walk_in_optimized_order(Real_Matrix_Preserving_Visitor)
     * @see #walk_in_optimized_order(Real_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_optimized_order(Real_Matrix_Preserving_Visitor, int, int, int, int)
     * @return the value returned by {@link Real_Matrix_Changing_Visitor#end()} at the end
     * of the walk
     */
    virtual double walk_in_row_order(Real_Matrix_Changing_Visitor visitor, int start_row, int end_row, int start_column, int end_column) = 0;

    /**
     * Visit (but don't change) some matrix entries in row order.
     * <p>Row order starts at upper left and iterating through all elements
     * of a row from left to right before going to the leftmost element
     * of the next row.</p>
     * @param visitor visitor used to process all matrix entries
     * @param start_row Initial row index
     * @param end_row Final row index (inclusive)
     * @param start_column Initial column index
     * @param end_column Final column index
     * @ if the indices are not valid.
     * @ if {@code end_row < start_row} or
     * {@code end_column < start_column}.
     * @see #walk_in_row_order(Real_Matrix_Changing_Visitor)
     * @see #walk_in_row_order(Real_Matrix_Preserving_Visitor)
     * @see #walk_in_row_order(Real_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_column_order(Real_Matrix_Changing_Visitor)
     * @see #walk_in_column_order(Real_Matrix_Preserving_Visitor)
     * @see #walk_in_column_order(Real_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_column_order(Real_Matrix_Preserving_Visitor, int, int, int, int)
     * @see #walk_in_optimized_order(Real_Matrix_Changing_Visitor)
     * @see #walk_in_optimized_order(Real_Matrix_Preserving_Visitor)
     * @see #walk_in_optimized_order(Real_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_optimized_order(Real_Matrix_Preserving_Visitor, int, int, int, int)
     * @return the value returned by {@link Real_Matrix_Preserving_Visitor#end()} at the end
     * of the walk
     */
    virtual double walk_in_row_order(Real_Matrix_Preserving_Visitor visitor, int start_row, int end_row, int start_column, int end_column) = 0;

    /**
     * Visit (and possibly change) all matrix entries in column order.
     * <p>Column order starts at upper left and iterating through all elements
     * of a column from top to bottom before going to the topmost element
     * of the next column.</p>
     * @param visitor visitor used to process all matrix entries
     * @see #walk_in_row_order(Real_Matrix_Changing_Visitor)
     * @see #walk_in_row_order(Real_Matrix_Preserving_Visitor)
     * @see #walk_in_row_order(Real_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_row_order(Real_Matrix_Preserving_Visitor, int, int, int, int)
     * @see #walk_in_column_order(Real_Matrix_Preserving_Visitor)
     * @see #walk_in_column_order(Real_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_column_order(Real_Matrix_Preserving_Visitor, int, int, int, int)
     * @see #walk_in_optimized_order(Real_Matrix_Changing_Visitor)
     * @see #walk_in_optimized_order(Real_Matrix_Preserving_Visitor)
     * @see #walk_in_optimized_order(Real_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_optimized_order(Real_Matrix_Preserving_Visitor, int, int, int, int)
     * @return the value returned by {@link Real_Matrix_Changing_Visitor#end()} at the end
     * of the walk
     */
    virtual double walk_in_column_order(Real_Matrix_Changing_Visitor visitor) = 0;

    /**
     * Visit (but don't change) all matrix entries in column order.
     * <p>Column order starts at upper left and iterating through all elements
     * of a column from top to bottom before going to the topmost element
     * of the next column.</p>
     * @param visitor visitor used to process all matrix entries
     * @see #walk_in_row_order(Real_Matrix_Changing_Visitor)
     * @see #walk_in_row_order(Real_Matrix_Preserving_Visitor)
     * @see #walk_in_row_order(Real_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_row_order(Real_Matrix_Preserving_Visitor, int, int, int, int)
     * @see #walk_in_column_order(Real_Matrix_Changing_Visitor)
     * @see #walk_in_column_order(Real_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_column_order(Real_Matrix_Preserving_Visitor, int, int, int, int)
     * @see #walk_in_optimized_order(Real_Matrix_Changing_Visitor)
     * @see #walk_in_optimized_order(Real_Matrix_Preserving_Visitor)
     * @see #walk_in_optimized_order(Real_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_optimized_order(Real_Matrix_Preserving_Visitor, int, int, int, int)
     * @return the value returned by {@link Real_Matrix_Preserving_Visitor#end()} at the end
     * of the walk
     */
    virtual double walk_in_column_order(Real_Matrix_Preserving_Visitor visitor) = 0;

    /**
     * Visit (and possibly change) some matrix entries in column order.
     * <p>Column order starts at upper left and iterating through all elements
     * of a column from top to bottom before going to the topmost element
     * of the next column.</p>
     * @param visitor visitor used to process all matrix entries
     * @param start_row Initial row index
     * @param end_row Final row index (inclusive)
     * @param start_column Initial column index
     * @param end_column Final column index
     * @ if the indices are not valid.
     * @ if {@code end_row < start_row} or
     * {@code end_column < start_column}.
     * @see #walk_in_row_order(Real_Matrix_Changing_Visitor)
     * @see #walk_in_row_order(Real_Matrix_Preserving_Visitor)
     * @see #walk_in_row_order(Real_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_row_order(Real_Matrix_Preserving_Visitor, int, int, int, int)
     * @see #walk_in_column_order(Real_Matrix_Changing_Visitor)
     * @see #walk_in_column_order(Real_Matrix_Preserving_Visitor)
     * @see #walk_in_column_order(Real_Matrix_Preserving_Visitor, int, int, int, int)
     * @see #walk_in_optimized_order(Real_Matrix_Changing_Visitor)
     * @see #walk_in_optimized_order(Real_Matrix_Preserving_Visitor)
     * @see #walk_in_optimized_order(Real_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_optimized_order(Real_Matrix_Preserving_Visitor, int, int, int, int)
     * @return the value returned by {@link Real_Matrix_Changing_Visitor#end()} at the end
     * of the walk
     */
    virtual double walk_in_column_order(Real_Matrix_Changing_Visitor visitor, int start_row, int end_row, int start_column, int end_column) = 0;

    /**
     * Visit (but don't change) some matrix entries in column order.
     * <p>Column order starts at upper left and iterating through all elements
     * of a column from top to bottom before going to the topmost element
     * of the next column.</p>
     * @param visitor visitor used to process all matrix entries
     * @param start_row Initial row index
     * @param end_row Final row index (inclusive)
     * @param start_column Initial column index
     * @param end_column Final column index
     * @ if the indices are not valid.
     * @ if {@code end_row < start_row} or
     * {@code end_column < start_column}.
     * @see #walk_in_row_order(Real_Matrix_Changing_Visitor)
     * @see #walk_in_row_order(Real_Matrix_Preserving_Visitor)
     * @see #walk_in_row_order(Real_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_row_order(Real_Matrix_Preserving_Visitor, int, int, int, int)
     * @see #walk_in_column_order(Real_Matrix_Changing_Visitor)
     * @see #walk_in_column_order(Real_Matrix_Preserving_Visitor)
     * @see #walk_in_column_order(Real_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_optimized_order(Real_Matrix_Changing_Visitor)
     * @see #walk_in_optimized_order(Real_Matrix_Preserving_Visitor)
     * @see #walk_in_optimized_order(Real_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_optimized_order(Real_Matrix_Preserving_Visitor, int, int, int, int)
     * @return the value returned by {@link Real_Matrix_Preserving_Visitor#end()} at the end
     * of the walk
     */
    virtual double walk_in_column_order(Real_Matrix_Preserving_Visitor visitor, int start_row, int end_row, int start_column, int end_column) = 0;

    /**
     * Visit (and possibly change) all matrix entries using the fastest possible order.
     * <p>The fastest walking order depends on the exact matrix class. It may be
     * different from traditional row or column orders.</p>
     * @param visitor visitor used to process all matrix entries
     * @see #walk_in_row_order(Real_Matrix_Changing_Visitor)
     * @see #walk_in_row_order(Real_Matrix_Preserving_Visitor)
     * @see #walk_in_row_order(Real_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_row_order(Real_Matrix_Preserving_Visitor, int, int, int, int)
     * @see #walk_in_column_order(Real_Matrix_Changing_Visitor)
     * @see #walk_in_column_order(Real_Matrix_Preserving_Visitor)
     * @see #walk_in_column_order(Real_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_column_order(Real_Matrix_Preserving_Visitor, int, int, int, int)
     * @see #walk_in_optimized_order(Real_Matrix_Preserving_Visitor)
     * @see #walk_in_optimized_order(Real_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_optimized_order(Real_Matrix_Preserving_Visitor, int, int, int, int)
     * @return the value returned by {@link Real_Matrix_Changing_Visitor#end()} at the end
     * of the walk
     */
    virtual double walk_in_optimized_order(Real_Matrix_Changing_Visitor visitor) = 0;

    /**
     * Visit (but don't change) all matrix entries using the fastest possible order.
     * <p>The fastest walking order depends on the exact matrix class. It may be
     * different from traditional row or column orders.</p>
     * @param visitor visitor used to process all matrix entries
     * @see #walk_in_row_order(Real_Matrix_Changing_Visitor)
     * @see #walk_in_row_order(Real_Matrix_Preserving_Visitor)
     * @see #walk_in_row_order(Real_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_row_order(Real_Matrix_Preserving_Visitor, int, int, int, int)
     * @see #walk_in_column_order(Real_Matrix_Changing_Visitor)
     * @see #walk_in_column_order(Real_Matrix_Preserving_Visitor)
     * @see #walk_in_column_order(Real_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_column_order(Real_Matrix_Preserving_Visitor, int, int, int, int)
     * @see #walk_in_optimized_order(Real_Matrix_Changing_Visitor)
     * @see #walk_in_optimized_order(Real_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_optimized_order(Real_Matrix_Preserving_Visitor, int, int, int, int)
     * @return the value returned by {@link Real_Matrix_Preserving_Visitor#end()} at the end
     * of the walk
     */
    virtual double walk_in_optimized_order(Real_Matrix_Preserving_Visitor visitor) = 0;

    /**
     * Visit (and possibly change) some matrix entries using the fastest possible order.
     * <p>The fastest walking order depends on the exact matrix class. It may be
     * different from traditional row or column orders.</p>
     * @param visitor visitor used to process all matrix entries
     * @param start_row Initial row index
     * @param end_row Final row index (inclusive)
     * @param start_column Initial column index
     * @param end_column Final column index (inclusive)
     * @ if the indices are not valid.
     * @ if {@code end_row < start_row} or
     * {@code end_column < start_column}.
     * @see #walk_in_row_order(Real_Matrix_Changing_Visitor)
     * @see #walk_in_row_order(Real_Matrix_Preserving_Visitor)
     * @see #walk_in_row_order(Real_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_row_order(Real_Matrix_Preserving_Visitor, int, int, int, int)
     * @see #walk_in_column_order(Real_Matrix_Changing_Visitor)
     * @see #walk_in_column_order(Real_Matrix_Preserving_Visitor)
     * @see #walk_in_column_order(Real_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_column_order(Real_Matrix_Preserving_Visitor, int, int, int, int)
     * @see #walk_in_optimized_order(Real_Matrix_Changing_Visitor)
     * @see #walk_in_optimized_order(Real_Matrix_Preserving_Visitor)
     * @see #walk_in_optimized_order(Real_Matrix_Preserving_Visitor, int, int, int, int)
     * @return the value returned by {@link Real_Matrix_Changing_Visitor#end()} at the end
     * of the walk
     */
    virtual double walk_in_optimized_order(Real_Matrix_Changing_Visitor visitor, int start_row, int end_row, int start_column, int end_column) = 0;

    /**
     * Visit (but don't change) some matrix entries using the fastest possible order.
     * <p>The fastest walking order depends on the exact matrix class. It may be
     * different from traditional row or column orders.</p>
     * @param visitor visitor used to process all matrix entries
     * @param start_row Initial row index
     * @param end_row Final row index (inclusive)
     * @param start_column Initial column index
     * @param end_column Final column index (inclusive)
     * @ if the indices are not valid.
     * @ if {@code end_row < start_row} or
     * {@code end_column < start_column}.
     * @see #walk_in_row_order(Real_Matrix_Changing_Visitor)
     * @see #walk_in_row_order(Real_Matrix_Preserving_Visitor)
     * @see #walk_in_row_order(Real_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_row_order(Real_Matrix_Preserving_Visitor, int, int, int, int)
     * @see #walk_in_column_order(Real_Matrix_Changing_Visitor)
     * @see #walk_in_column_order(Real_Matrix_Preserving_Visitor)
     * @see #walk_in_column_order(Real_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_column_order(Real_Matrix_Preserving_Visitor, int, int, int, int)
     * @see #walk_in_optimized_order(Real_Matrix_Changing_Visitor)
     * @see #walk_in_optimized_order(Real_Matrix_Preserving_Visitor)
     * @see #walk_in_optimized_order(Real_Matrix_Changing_Visitor, int, int, int, int)
     * @return the value returned by {@link Real_Matrix_Preserving_Visitor#end()} at the end
     * of the walk
     */
    virtual double walk_in_optimized_order(Real_Matrix_Preserving_Visitor visitor, int start_row, int end_row, int start_column, int end_column) = 0;

    /**
     * Acts as if implemented as:
     * <pre>
     *  return copy().map_to_self(function);
     * </pre>
     * Returns a matrix. Does not change instance data.
     *
     * @param function Function to apply to each entry.
     * @return a matrix.
     * @since 1.7
     */
    Real_Matrix map(Univariate_Function function) 
    {
        return copy().map_to_self(function);
    }

    /**
     * Replace each entry by the result of applying the function to it.
     *
     * @param function Function to apply to each entry.
     * @return a reference to this matrix.
     * @since 1.7
     */
    Real_Matrix map_to_self(const Univariate_Function function) 
    {
        walk_in_optimized_order(
            Real_Matrix_Changing_Visitor() 
            {
            public:
                /** {@inherit_doc} */
                //override
                double visit(const int& row, const int& column, double value)
                {
                    // apply the function to the current entry
                    return function.value(value);
                }

                /** {@inherit_doc} */
                //override
                void start(const int& rows, int columns, int start_row, int end_row, int start_column, int end_column) const
                {
                }

                /** {@inherit_doc} */
                //override
                double end() const
                {
                    return 0;
                }

            }
        );

        return this;

    }

};