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


//import java.util.function.Function;

//import org.hipparchus.Field;
//import org.hipparchus.Field_Element;
//import org.hipparchus.exception.;
#include <type_traits>
#include <vector>
#include "../FieldElement.h"
#include "../Field.h"
//import org.hipparchus.exception.;

/**
 * Interface defining field-valued matrix with basic algebraic operations.
 * <p>
 * Matrix element indexing is 0-based -- e.g., <code>get_entry(0, 0)</code>
 * returns the element in the first row, first column of the matrix.</p>
 *
 * @param <T> the type of the field elements
 */
template<typename T, typename std::enable_if<std::is_base_of<Field_Element<T>, T>::value>::type* = nullptr>
class Field_Matrix : public Any_Matrix 
{
    /**
     * Get the type of field elements of the matrix.
     *
     * @return the type of field elements of the matrix.
     */
    virtual Field<T> get_field();

    /**
     * Create a Field_Matrix<T> of the same type as the instance with
     * the supplied row and column dimensions.
     *
     * @param row_dimension  the number of rows in the matrix
     * @param column_dimension  the number of columns in the matrix
     * @return a matrix of the same type as the instance
     * @ if row or column dimension is not
     * positive.
     */
    virtual Field_Matrix<T> create_matrix(const int& row_dimension, const int& column_dimension);

    /**
     * Make a (deep) copy of this.
     *
     * @return a copy of this matrix.
     */
    virtual Field_Matrix<T> copy();

    /**
     * Compute the sum of this and m.
     *
     * @param m Matrix to be added.
     * @return {@code this} + {@code m}.
     * @ if {@code m} is not the same
     * size as {@code this} matrix.
     */
    virtual Field_Matrix<T> add(const Field_Matrix<T>& m);

    /**
     * Subtract {@code m} from this matrix.
     *
     * @param m Matrix to be subtracted.
     * @return {@code this} - {@code m}.
     * @ if {@code m} is not the same
     * size as {@code this} matrix.
     */
    virtual Field_Matrix<T> subtract(const Field_Matrix<T>& m);

     /**
     * Increment each entry of this matrix.
     *
     * @param d Value to be added to each entry.
     * @return {@code d} + {@code this}.
     */
    virtual Field_Matrix<T> scalar_add(const T& d);

    /**
     * Multiply each entry by {@code d}.
     *
     * @param d Value to multiply all entries by.
     * @return {@code d} * {@code this}.
     */
    virtual Field_Matrix<T> scalar_multiply(const T& d);

    /**
     * Postmultiply this matrix by {@code m}.
     *
     * @param m  Matrix to postmultiply by.
     * @return {@code this} * {@code m}.
     * @ if the number of columns of
     * {@code this} matrix is not equal to the number of rows of matrix
     * {@code m}.
     */
    virtual Field_Matrix<T> multiply(const Field_Matrix<T>& m);

    /**
     * Returns the result of postmultiplying {@code this} by {@code m^T}.
     * <p>
     * This is equivalent to call {@link #multiply(Field_Matrix) multiply}(m.{@link #transpose()}), * but some implementations may avoid building the intermediate transposed matrix.
     * </p>
     * @param m matrix to first transpose and second postmultiply by
     * @return {@code this * m^T}
     * @ if
     * {@code column_dimension(this) != column_dimension(m)}
     * @since 1.3
     */
    default Field_Matrix<T> multiply_transposed(const Field_Matrix<T>& m)
    {
        return multiply(m.transpose());
    }

    /**
     * Returns the result of postmultiplying {@code this^T} by {@code m}.
     * <p>
     * This is equivalent to call {@link #transpose()}.{@link #multiply(Field_Matrix) multiply(m)}, * but some implementations may avoid building the intermediate transposed matrix.
     * </p>
     * @param m matrix to postmultiply by
     * @return {@code this^T * m}
     * @ if
     * {@code column_dimension(this) != column_dimension(m)}
     * @since 1.3
     */
    default Field_Matrix<T> transpose_multiply(const Field_Matrix<T>& m)
    {
        return transpose().multiply(m);
    }

    /**
     * Premultiply this matrix by {@code m}.
     *
     * @param m Matrix to premultiply by.
     * @return {@code m} * {@code this}.
     * @ if the number of columns of {@code m}
     * differs from the number of rows of {@code this} matrix.
     */
    virtual Field_Matrix<T> pre_multiply(const Field_Matrix<T>& m);

    /**
     * Returns the result multiplying this with itself <code>p</code> times.
     * Depending on the type of the field elements, T, instability for high
     * powers might occur.
     *
     * @param p raise this to power p
     * @return this^p
     * @ if {@code p < 0}
     * @ if {@code this matrix} is not square
     */
    virtual Field_Matrix<T> power(const int& p);

    /**
     * Returns matrix entries as a two-dimensional array.
     *
     * @return a 2-dimensional array of entries.
     */
    virtual std::vector<std::vector<T>> get_data();

    /**
     * Get a submatrix. Rows and columns are indicated
     * counting from 0 to n - 1.
     *
     * @param start_row Initial row index
     * @param end_row Final row index (inclusive)
     * @param start_column Initial column index
     * @param end_column Final column index (inclusive)
     * @return the matrix containing the data of the specified rows and columns.
     * @ is {@code end_row < start_row} of
     * {@code end_column < start_column}.
     * @ if the indices are not valid.
     */
    virtual Field_Matrix<T> get_sub_matrix(const int& start_row, const int& end_row, const int& start_column, const int& end_column);

   /**
    * Get a submatrix. Rows and columns are indicated
    * counting from 0 to n - 1.
    *
    * @param selected_rows Array of row indices.
    * @param selected_columns Array of column indices.
    * @return the matrix containing the data in the
    * specified rows and columns.
    * @ if {@code selected_rows} or
    * {@code selected_columns} is empty
    * @ if {@code selected_rows} or
    * {@code selected_columns} is {@code NULL}.
    * @ if row or column selections are not valid.
    */
   virtual Field_Matrix<T> get_sub_matrix(const std::vector<int>& selected_rows, const std::vector<int>& selected_columns);

   /**
    * Copy a submatrix. Rows and columns are 0-based. The designated submatrix
    * is copied into the top left portion of the destination array.
    *
    * @param start_row Initial row index.
    * @param end_row Final row index (inclusive).
    * @param start_column Initial column index.
    * @param end_column Final column index (inclusive).
    * @param destination The array where the submatrix data should be copied
    * (if larger than rows/columns counts, only the upper-left part will be modified).
    * @ if the dimensions of
    * {@code destination} are not large enough to hold the submatrix.
    * @ if {@code end_row < start_row} or
    * {@code end_column < start_column}.
    * @ if the indices are not valid.
    */
   virtual void copy_sub_matrix(const int& start_row, const int& end_row, const int& start_column, const int& end_column, const std::vector<std::vector<T>>& destination);

  /**
   * Copy a submatrix. Rows and columns are indicated
   * counting from 0 to n - 1.
   *
   * @param selected_rows Array of row indices.
   * @param selected_columns Array of column indices.
   * @param destination Arrays where the submatrix data should be copied
   * (if larger than rows/columns counts, only the upper-left part will be used)
   * @ if the dimensions of
   * {@code destination} do not match those of {@code this}.
   * @ if {@code selected_rows} or
   * {@code selected_columns} is empty
   * @ if {@code selected_rows} or
   * {@code selected_columns} is {@code NULL}.
   * @ if the indices are not valid.
   */
   virtual void copy_sub_matrix(const std::vector<int>& selected_rows, const std::vector<int>& selected_columns, const std::vector<std::vector<T>>& destination);

    /**
     * Replace the submatrix starting at {@code (row, column)} using data in the
     * input {@code sub_matrix} array. Indexes are 0-based.
     * <p>
     * Example:<br>
     * Starting with
     *
     * <pre>
     * 1  2  3  4
     * 5  6  7  8
     * 9  0  1  2
     * </pre>
     *
     * and <code>sub_matrix = {{3, 4} {5,6}}</code>, invoking
     * <code>set_sub_matrix(sub_matrix,1,1))</code> will result in
     *
     * <pre>
     * 1  2  3  4
     * 5  3  4  8
     * 9  5  6  2
     * </pre>
     *
     * </p>
     *
     * @param sub_matrix Array containing the submatrix replacement data.
     * @param row Row coordinate of the top-left element to be replaced.
     * @param column Column coordinate of the top-left element to be replaced.
     * @ if {@code sub_matrix} does not fit into this
     * matrix from element in {@code (row, column)}.
     * @ if a row or column of {@code sub_matrix} is empty.
     * @ if {@code sub_matrix} is not
     * rectangular (not all rows have the same length).
     * @ if {@code sub_matrix} is {@code NULL}.
     */
   virtual void set_sub_matrix(const std::vector<std::vector<T>>& sub_matrix, const int& row, const int& column);

   /**
    * Get the entries in row number {@code row}
    * as a row matrix.
    *
    * @param row Row to be fetched.
    * @return a row matrix.
    * @ if the specified row index is invalid.
    */
   virtual Field_Matrix<T> get_row_matrix(const int& row);

   /**
    * Set the entries in row number {@code row}
    * as a row matrix.
    *
    * @param row Row to be set.
    * @param matrix Row matrix (must have one row and the same number
    * of columns as the instance).
    * @ if the specified row index is invalid.
    * @
    * if the matrix dimensions do not match one instance row.
    */
   virtual void set_row_matrix(const int& row, const Field_Matrix<T>& matrix);

   /**
    * Get the entries in column number {@code column}
    * as a column matrix.
    *
    * @param column Column to be fetched.
    * @return a column matrix.
    * @ if the specified column index is invalid.
    */
   virtual Field_Matrix<T> get_column_matrix(const int& column);

   /**
    * Set the entries in column number {@code column}
    * as a column matrix.
    *
    * @param column Column to be set.
    * @param matrix column matrix (must have one column and the same
    * number of rows as the instance).
    * @ if the specified column index is invalid.
    * @ if the matrix dimensions do
    * not match one instance column.
    */
   virtual void set_column_matrix(const int& column, const Field_Matrix<T>& matrix);

   /**
    * Get the entries in row number {@code row}
    * as a vector.
    *
    * @param row Row to be fetched
    * @return a row vector.
    * @ if the specified row index is invalid.
    */
   virtual Field_Vector<T> get_row_vector(const int& row);

   /**
    * Set the entries in row number {@code row}
    * as a vector.
    *
    * @param row Row to be set.
    * @param vector row vector (must have the same number of columns
    * as the instance).
    * @ if the specified row index is invalid.
    * @ if the vector dimension does not
    * match one instance row.
    */
   virtual void set_row_vector(const int& row, const Field_Vector<T>& vector);

   /**
    * Returns the entries in column number {@code column}
    * as a vector.
    *
    * @param column Column to be fetched.
    * @return a column vector.
    * @ if the specified column index is invalid.
    */
   virtual Field_Vector<T> get_column_vector(const int& column);

   /**
    * Set the entries in column number {@code column}
    * as a vector.
    *
    * @param column Column to be set.
    * @param vector Column vector (must have the same number of rows
    * as the instance).
    * @ if the specified column index is invalid.
    * @ if the vector dimension does not
    * match one instance column.
    */
   virtual void set_column_vector(const int& column, const Field_Vector<T>& vector);

    /**
     * Get the entries in row number {@code row} as an array.
     *
     * @param row Row to be fetched.
     * @return array of entries in the row.
     * @ if the specified row index is not valid.
     */
   virtual  std::vector<T> get_row(const int& row);

    /**
     * Set the entries in row number {@code row}
     * as a row matrix.
     *
     * @param row Row to be set.
     * @param array Row matrix (must have the same number of columns as
     * the instance).
     * @ if the specified row index is invalid.
     * @ if the array size does not match
     * one instance row.
     */
   virtual void set_row(const int& row, const std::vector<T>& arr);

    /**
     * Get the entries in column number {@code col} as an array.
     *
     * @param column the column to be fetched
     * @return array of entries in the column
     * @ if the specified column index is not valid.
     */
   virtual std::vector<T> get_column(const int& column);

    /**
     * Set the entries in column number {@code column}
     * as a column matrix.
     *
     * @param column the column to be set
     * @param array column array (must have the same number of rows as the instance)
     * @ if the specified column index is invalid.
     * @ if the array size does not match
     * one instance column.
     */
   virtual void set_column(const int& column, const std::vector<T>& arr);

    /**
     * Returns the entry in the specified row and column.
     *
     * @param row  row location of entry to be fetched
     * @param column  column location of entry to be fetched
     * @return matrix entry in row,column
     * @ if the row or column index is not valid.
     */
   virtual  T get_entry(const int& row, const int& column);

    /**
     * Set the entry in the specified row and column.
     *
     * @param row  row location of entry to be set
     * @param column  column location of entry to be set
     * @param value matrix entry to be set in row,column
     * @ if the row or column index is not valid.
     */
   virtual void set_entry(const int& row, const int& column, const T& value);

    /**
     * Change an entry in the specified row and column.
     *
     * @param row Row location of entry to be set.
     * @param column Column location of entry to be set.
     * @param increment Value to add to the current matrix entry in
     * {@code (row, column)}.
     * @ if the row or column index is not valid.
     */
   virtual  void add_to_entry(const int& row, const int& column, const T& increment);

    /**
     * Change an entry in the specified row and column.
     *
     * @param row Row location of entry to be set.
     * @param column Column location of entry to be set.
     * @param factor Multiplication factor for the current matrix entry
     * in {@code (row,column)}
     * @ if the row or column index is not valid.
     */
   virtual void multiply_entry(const int& row, const int& column, const T& factor);

    /**
     * Returns the transpose of this matrix.
     *
     * @return transpose matrix
     */
   virtual  Field_Matrix<T> transpose();

    /**
     * Returns the <a href="http://mathworld.wolfram.com/MatrixTrace.html">
     * trace</a> of the matrix (the sum of the elements on the main diagonal).
     *
     * @return trace
     * @ if the matrix is not square.
     */
   virtual T get_trace();

    /**
     * Returns the result of multiplying this by the vector {@code v}.
     *
     * @param v the vector to operate on
     * @return {@code this * v}
     * @ if the number of columns of
     * {@code this} matrix is not equal to the size of the vector {@code v}.
     */
   virtual std::vector<T> operate(const std::vector<T>& v);

    /**
     * Returns the result of multiplying this by the vector {@code v}.
     *
     * @param v the vector to operate on
     * @return {@code this * v}
     * @ if the number of columns of
     * {@code this} matrix is not equal to the size of the vector {@code v}.
     */
   virtual Field_Vector<T> operate(const Field_Vector<T>& v);

    /**
     * Returns the (row) vector result of premultiplying this by the vector
     * {@code v}.
     *
     * @param v the row vector to premultiply by
     * @return {@code v * this}
     * @ if the number of rows of {@code this}
     * matrix is not equal to the size of the vector {@code v}
     */
   virtual std::vector<T> pre_multiply(const std::vector<T>& v);

    /**
     * Returns the (row) vector result of premultiplying this by the vector
     * {@code v}.
     *
     * @param v the row vector to premultiply by
     * @return {@code v * this}
     * @ if the number of rows of {@code this}
     * matrix is not equal to the size of the vector {@code v}
     */
   virtual Field_Vector<T> pre_multiply(const Field_Vector<T>& v);

    /**
     * Visit (and possibly change) all matrix entries in row order.
     * <p>Row order starts at upper left and iterating through all elements
     * of a row from left to right before going to the leftmost element
     * of the next row.</p>
     * @param visitor visitor used to process all matrix entries
     * @see #walk_in_row_order(Field_Matrix_Preserving_Visitor)
     * @see #walk_in_row_order(Field_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_row_order(Field_Matrix_Preserving_Visitor, int, int, int, int)
     * @see #walk_in_column_order(Field_Matrix_Changing_Visitor)
     * @see #walk_in_column_order(Field_Matrix_Preserving_Visitor)
     * @see #walk_in_column_order(Field_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_column_order(Field_Matrix_Preserving_Visitor, int, int, int, int)
     * @see #walk_in_optimized_order(Field_Matrix_Changing_Visitor)
     * @see #walk_in_optimized_order(Field_Matrix_Preserving_Visitor)
     * @see #walk_in_optimized_order(Field_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_optimized_order(Field_Matrix_Preserving_Visitor, int, int, int, int)
     * @return the value returned by {@link Field_Matrix_Changing_Visitor#end()} at the end
     * of the walk
     */
   virtual T walk_in_row_order(const Field_Matrix_Changing_Visitor<T>& visitor);

    /**
     * Visit (but don't change) all matrix entries in row order.
     * <p>Row order starts at upper left and iterating through all elements
     * of a row from left to right before going to the leftmost element
     * of the next row.</p>
     * @param visitor visitor used to process all matrix entries
     * @see #walk_in_row_order(Field_Matrix_Changing_Visitor)
     * @see #walk_in_row_order(Field_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_row_order(Field_Matrix_Preserving_Visitor, int, int, int, int)
     * @see #walk_in_column_order(Field_Matrix_Changing_Visitor)
     * @see #walk_in_column_order(Field_Matrix_Preserving_Visitor)
     * @see #walk_in_column_order(Field_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_column_order(Field_Matrix_Preserving_Visitor, int, int, int, int)
     * @see #walk_in_optimized_order(Field_Matrix_Changing_Visitor)
     * @see #walk_in_optimized_order(Field_Matrix_Preserving_Visitor)
     * @see #walk_in_optimized_order(Field_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_optimized_order(Field_Matrix_Preserving_Visitor, int, int, int, int)
     * @return the value returned by {@link Field_Matrix_Preserving_Visitor#end()} at the end
     * of the walk
     */
   virtual T walk_in_row_order(const Field_Matrix_Preserving_Visitor<T>& visitor);

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
     * @see #walk_in_row_order(Field_Matrix_Changing_Visitor)
     * @see #walk_in_row_order(Field_Matrix_Preserving_Visitor)
     * @see #walk_in_row_order(Field_Matrix_Preserving_Visitor, int, int, int, int)
     * @see #walk_in_column_order(Field_Matrix_Changing_Visitor)
     * @see #walk_in_column_order(Field_Matrix_Preserving_Visitor)
     * @see #walk_in_column_order(Field_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_column_order(Field_Matrix_Preserving_Visitor, int, int, int, int)
     * @see #walk_in_optimized_order(Field_Matrix_Changing_Visitor)
     * @see #walk_in_optimized_order(Field_Matrix_Preserving_Visitor)
     * @see #walk_in_optimized_order(Field_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_optimized_order(Field_Matrix_Preserving_Visitor, int, int, int, int)
     * @return the value returned by {@link Field_Matrix_Changing_Visitor#end()} at the end
     * of the walk
     */
   virtual T walk_in_row_order(const Field_Matrix_Changing_Visitor<T>& visitor, const int& start_row, const int& end_row, const int& start_column, const int& end_column);

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
     * @see #walk_in_row_order(Field_Matrix_Changing_Visitor)
     * @see #walk_in_row_order(Field_Matrix_Preserving_Visitor)
     * @see #walk_in_row_order(Field_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_column_order(Field_Matrix_Changing_Visitor)
     * @see #walk_in_column_order(Field_Matrix_Preserving_Visitor)
     * @see #walk_in_column_order(Field_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_column_order(Field_Matrix_Preserving_Visitor, int, int, int, int)
     * @see #walk_in_optimized_order(Field_Matrix_Changing_Visitor)
     * @see #walk_in_optimized_order(Field_Matrix_Preserving_Visitor)
     * @see #walk_in_optimized_order(Field_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_optimized_order(Field_Matrix_Preserving_Visitor, int, int, int, int)
     * @return the value returned by {@link Field_Matrix_Preserving_Visitor#end()} at the end
     * of the walk
     */
   virtual T walk_in_row_order(const Field_Matrix_Preserving_Visitor<T>& visitor, const int& start_row, const int end_row, const int& start_column, const int& end_column);

    /**
     * Visit (and possibly change) all matrix entries in column order.
     * <p>Column order starts at upper left and iterating through all elements
     * of a column from top to bottom before going to the topmost element
     * of the next column.</p>
     * @param visitor visitor used to process all matrix entries
     * @see #walk_in_row_order(Field_Matrix_Changing_Visitor)
     * @see #walk_in_row_order(Field_Matrix_Preserving_Visitor)
     * @see #walk_in_row_order(Field_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_row_order(Field_Matrix_Preserving_Visitor, int, int, int, int)
     * @see #walk_in_column_order(Field_Matrix_Preserving_Visitor)
     * @see #walk_in_column_order(Field_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_column_order(Field_Matrix_Preserving_Visitor, int, int, int, int)
     * @see #walk_in_optimized_order(Field_Matrix_Changing_Visitor)
     * @see #walk_in_optimized_order(Field_Matrix_Preserving_Visitor)
     * @see #walk_in_optimized_order(Field_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_optimized_order(Field_Matrix_Preserving_Visitor, int, int, int, int)
     * @return the value returned by {@link Field_Matrix_Changing_Visitor#end()} at the end
     * of the walk
     */
   virtual T walk_in_column_order(const Field_Matrix_Changing_Visitor<T>& visitor);

    /**
     * Visit (but don't change) all matrix entries in column order.
     * <p>Column order starts at upper left and iterating through all elements
     * of a column from top to bottom before going to the topmost element
     * of the next column.</p>
     * @param visitor visitor used to process all matrix entries
     * @see #walk_in_row_order(Field_Matrix_Changing_Visitor)
     * @see #walk_in_row_order(Field_Matrix_Preserving_Visitor)
     * @see #walk_in_row_order(Field_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_row_order(Field_Matrix_Preserving_Visitor, int, int, int, int)
     * @see #walk_in_column_order(Field_Matrix_Changing_Visitor)
     * @see #walk_in_column_order(Field_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_column_order(Field_Matrix_Preserving_Visitor, int, int, int, int)
     * @see #walk_in_optimized_order(Field_Matrix_Changing_Visitor)
     * @see #walk_in_optimized_order(Field_Matrix_Preserving_Visitor)
     * @see #walk_in_optimized_order(Field_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_optimized_order(Field_Matrix_Preserving_Visitor, int, int, int, int)
     * @return the value returned by {@link Field_Matrix_Preserving_Visitor#end()} at the end
     * of the walk
     */
   virtual T walk_in_column_order(const Field_Matrix_Preserving_Visitor<T>& visitor);

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
     * @ if {@code end_row < start_row} or
     * {@code end_column < start_column}.
     * @ if the indices are not valid.
     * @see #walk_in_row_order(Field_Matrix_Changing_Visitor)
     * @see #walk_in_row_order(Field_Matrix_Preserving_Visitor)
     * @see #walk_in_row_order(Field_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_row_order(Field_Matrix_Preserving_Visitor, int, int, int, int)
     * @see #walk_in_column_order(Field_Matrix_Changing_Visitor)
     * @see #walk_in_column_order(Field_Matrix_Preserving_Visitor)
     * @see #walk_in_column_order(Field_Matrix_Preserving_Visitor, int, int, int, int)
     * @see #walk_in_optimized_order(Field_Matrix_Changing_Visitor)
     * @see #walk_in_optimized_order(Field_Matrix_Preserving_Visitor)
     * @see #walk_in_optimized_order(Field_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_optimized_order(Field_Matrix_Preserving_Visitor, int, int, int, int)
     * @return the value returned by {@link Field_Matrix_Changing_Visitor#end()} at the end
     * of the walk
     */
   virtual T walk_in_column_order(const Field_Matrix_Changing_Visitor<T>& visitor, const int& start_row, const int& end_row, const int& start_column, const int& end_column);

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
     * @ if {@code end_row < start_row} or
     * {@code end_column < start_column}.
     * @ if the indices are not valid.
     * @see #walk_in_row_order(Field_Matrix_Changing_Visitor)
     * @see #walk_in_row_order(Field_Matrix_Preserving_Visitor)
     * @see #walk_in_row_order(Field_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_row_order(Field_Matrix_Preserving_Visitor, int, int, int, int)
     * @see #walk_in_column_order(Field_Matrix_Changing_Visitor)
     * @see #walk_in_column_order(Field_Matrix_Preserving_Visitor)
     * @see #walk_in_column_order(Field_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_optimized_order(Field_Matrix_Changing_Visitor)
     * @see #walk_in_optimized_order(Field_Matrix_Preserving_Visitor)
     * @see #walk_in_optimized_order(Field_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_optimized_order(Field_Matrix_Preserving_Visitor, int, int, int, int)
     * @return the value returned by {@link Field_Matrix_Preserving_Visitor#end()} at the end
     * of the walk
     */
   virtual T walk_in_column_order(const Field_Matrix_Preserving_Visitor<T>& visitor, const int& start_row, const int& end_row, const int& start_column, const int& end_column);

    /**
     * Visit (and possibly change) all matrix entries using the fastest possible order.
     * <p>The fastest walking order depends on the exact matrix class. It may be
     * different from traditional row or column orders.</p>
     * @param visitor visitor used to process all matrix entries
     * @see #walk_in_row_order(Field_Matrix_Changing_Visitor)
     * @see #walk_in_row_order(Field_Matrix_Preserving_Visitor)
     * @see #walk_in_row_order(Field_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_row_order(Field_Matrix_Preserving_Visitor, int, int, int, int)
     * @see #walk_in_column_order(Field_Matrix_Changing_Visitor)
     * @see #walk_in_column_order(Field_Matrix_Preserving_Visitor)
     * @see #walk_in_column_order(Field_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_column_order(Field_Matrix_Preserving_Visitor, int, int, int, int)
     * @see #walk_in_optimized_order(Field_Matrix_Preserving_Visitor)
     * @see #walk_in_optimized_order(Field_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_optimized_order(Field_Matrix_Preserving_Visitor, int, int, int, int)
     * @return the value returned by {@link Field_Matrix_Changing_Visitor#end()} at the end
     * of the walk
     */
   virtual T walk_in_optimized_order(const Field_Matrix_Changing_Visitor<T>& visitor);

    /**
     * Visit (but don't change) all matrix entries using the fastest possible order.
     * <p>The fastest walking order depends on the exact matrix class. It may be
     * different from traditional row or column orders.</p>
     * @param visitor visitor used to process all matrix entries
     * @see #walk_in_row_order(Field_Matrix_Changing_Visitor)
     * @see #walk_in_row_order(Field_Matrix_Preserving_Visitor)
     * @see #walk_in_row_order(Field_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_row_order(Field_Matrix_Preserving_Visitor, int, int, int, int)
     * @see #walk_in_column_order(Field_Matrix_Changing_Visitor)
     * @see #walk_in_column_order(Field_Matrix_Preserving_Visitor)
     * @see #walk_in_column_order(Field_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_column_order(Field_Matrix_Preserving_Visitor, int, int, int, int)
     * @see #walk_in_optimized_order(Field_Matrix_Changing_Visitor)
     * @see #walk_in_optimized_order(Field_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_optimized_order(Field_Matrix_Preserving_Visitor, int, int, int, int)
     * @return the value returned by {@link Field_Matrix_Preserving_Visitor#end()} at the end
     * of the walk
     */
   virtual T walk_in_optimized_order(const Field_Matrix_Preserving_Visitor<T>& visitor);

    /**
     * Visit (and possibly change) some matrix entries using the fastest possible order.
     * <p>The fastest walking order depends on the exact matrix class. It may be
     * different from traditional row or column orders.</p>
     * @param visitor visitor used to process all matrix entries
     * @param start_row Initial row index
     * @param end_row Final row index (inclusive)
     * @param start_column Initial column index
     * @param end_column Final column index (inclusive)
     * @ if {@code end_row < start_row} or
     * {@code end_column < start_column}.
     * @ if the indices are not valid.
     * @see #walk_in_row_order(Field_Matrix_Changing_Visitor)
     * @see #walk_in_row_order(Field_Matrix_Preserving_Visitor)
     * @see #walk_in_row_order(Field_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_row_order(Field_Matrix_Preserving_Visitor, int, int, int, int)
     * @see #walk_in_column_order(Field_Matrix_Changing_Visitor)
     * @see #walk_in_column_order(Field_Matrix_Preserving_Visitor)
     * @see #walk_in_column_order(Field_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_column_order(Field_Matrix_Preserving_Visitor, int, int, int, int)
     * @see #walk_in_optimized_order(Field_Matrix_Changing_Visitor)
     * @see #walk_in_optimized_order(Field_Matrix_Preserving_Visitor)
     * @see #walk_in_optimized_order(Field_Matrix_Preserving_Visitor, int, int, int, int)
     * @return the value returned by {@link Field_Matrix_Changing_Visitor#end()} at the end
     * of the walk
     */
   virtual T walk_in_optimized_order(const Field_Matrix_Changing_Visitor<T>& visitor, const int& start_row, const int& end_row, const int& start_column, const int& end_column);

    /**
     * Visit (but don't change) some matrix entries using the fastest possible order.
     * <p>The fastest walking order depends on the exact matrix class. It may be
     * different from traditional row or column orders.</p>
     * @param visitor visitor used to process all matrix entries
     * @param start_row Initial row index
     * @param end_row Final row index (inclusive)
     * @param start_column Initial column index
     * @param end_column Final column index (inclusive)
     * @ if {@code end_row < start_row} or
     * {@code end_column < start_column}.
     * @ if the indices are not valid.
     * @see #walk_in_row_order(Field_Matrix_Changing_Visitor)
     * @see #walk_in_row_order(Field_Matrix_Preserving_Visitor)
     * @see #walk_in_row_order(Field_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_row_order(Field_Matrix_Preserving_Visitor, int, int, int, int)
     * @see #walk_in_column_order(Field_Matrix_Changing_Visitor)
     * @see #walk_in_column_order(Field_Matrix_Preserving_Visitor)
     * @see #walk_in_column_order(Field_Matrix_Changing_Visitor, int, int, int, int)
     * @see #walk_in_column_order(Field_Matrix_Preserving_Visitor, int, int, int, int)
     * @see #walk_in_optimized_order(Field_Matrix_Changing_Visitor)
     * @see #walk_in_optimized_order(Field_Matrix_Preserving_Visitor)
     * @see #walk_in_optimized_order(Field_Matrix_Changing_Visitor, int, int, int, int)
     * @return the value returned by {@link Field_Matrix_Preserving_Visitor#end()} at the end
     * of the walk
     */
   virtual T walk_in_optimized_order(const Field_Matrix_Preserving_Visitor<T>& visitor, const int& start_row, const int& end_row, const int& start_column, const int& end_column);

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
    default Field_Matrix<T> map(Function<T, T> function) 
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
    default Field_Matrix<T> map_to_self(const Function<T, T>&function)
    {
        throw std::exception("not implemented");
        //walk_in_optimized_order(Field_Matrix_Changing_Visitor<T>()
        //{
        //public:
        //    /** {@inherit_doc} */
        //    //override
        //    T visit(const int& row, const int& column, const T& value)
        //    {
        //        // apply the function to the current entry
        //        return function.apply(value);
        //    }

        //    /** {@inherit_doc} */
        //    //override
        //    void start(const int& rows, const int& columns, const int& start_row, const int& end_row, const int& start_column, const int& end_column)
        //    {
        //    }

        //    /** {@inherit_doc} */
        //    //override
        //    T end()
        //    {
        //        return get_field().get_zero();
        //    }

        //});

        //return *this;

    };
};