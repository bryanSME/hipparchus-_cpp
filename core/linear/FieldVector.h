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

//import org.hipparchus.Field;
//import org.hipparchus.Field_Element;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Runtime_Exception;
//import org.hipparchus.exception.;

/**
 * Interface defining a field-valued vector with basic algebraic operations.
 * <p>
 * vector element indexing is 0-based -- e.g., <code>get_entry(0)</code>
 * returns the first element of the vector.
 * </p>
 * <p>
 * The various <code>map_xxx</code> and <code>map_xxx_to_self</code> methods operate
 * on vectors element-wise, i.e. they perform the same operation (adding a scalar, * applying a function ...) on each element in turn. The <code>map_xxx</code>
 * versions create a vector to hold the result and do not change the instance.
 * The <code>map_xxx_to_self</code> versions use the instance itself to store the
 * results, so the instance is changed by these methods. In both cases, the result
 * vector is returned by the methods, this allows to use the <i>fluent API</i>
 * style, like this:
 * </p>
 * <pre>
 *   Real_Vector result = v.map_add_to_self(3.0).map_tan_to_self().map_square_to_self();
 * </pre>
 * <p>
 * Note that as almost all operations on {@link Field_Element} throw {@link
 * } when operating on a NULL element, it is the responsibility
 * of <code>Field_Vector</code> implementations to make sure no NULL elements
 * are inserted into the vector. This must be done in all constructors and
 * all setters.
 * <p>
 *
 * @param <T> the type of the field elements
 */
class Field_Vector<T extends Field_Element<T>>  
{

    /**
     * Get the type of field elements of the vector.
     * @return type of field elements of the vector
     */
    Field<T> get_field();

    /**
     * Returns a (deep) copy of this.
     * @return vector copy
     */
    Field_Vector<T> copy();

    /**
     * Compute the sum of {@code this} and {@code v}.
     * @param v vector to be added
     * @return {@code this + v}
     * @ if {@code v} is not the same size as {@code this}
     */
    Field_Vector<T> add(Field_Vector<T> v) ;

    /**
     * Compute {@code this} minus {@code v}.
     * @param v vector to be subtracted
     * @return {@code this - v}
     * @ if {@code v} is not the same size as {@code this}
     */
    Field_Vector<T> subtract(Field_Vector<T> v) ;

    /**
     * Map an addition operation to each entry.
     * @param d value to be added to each entry
     * @return {@code this + d}
     * @ if {@code d} is {@code NULL}.
     */
    Field_Vector<T> map_add(T d) ;

    /**
     * Map an addition operation to each entry.
     * <p>The instance <strong>is</strong> changed by this method.</p>
     * @param d value to be added to each entry
     * @return for convenience, return {@code this}
     * @ if {@code d} is {@code NULL}.
     */
    Field_Vector<T> map_add_to_self(T d) ;

    /**
     * Map a subtraction operation to each entry.
     * @param d value to be subtracted to each entry
     * @return {@code this - d}
     * @ if {@code d} is {@code NULL}
     */
    Field_Vector<T> map_subtract(T d) ;

    /**
     * Map a subtraction operation to each entry.
     * <p>The instance <strong>is</strong> changed by this method.</p>
     * @param d value to be subtracted to each entry
     * @return for convenience, return {@code this}
     * @ if {@code d} is {@code NULL}
     */
    Field_Vector<T> map_subtract_to_self(T d) ;

    /**
     * Map a multiplication operation to each entry.
     * @param d value to multiply all entries by
     * @return {@code this * d}
     * @ if {@code d} is {@code NULL}.
     */
    Field_Vector<T> map_multiply(T d) ;

    /**
     * Map a multiplication operation to each entry.
     * <p>The instance <strong>is</strong> changed by this method.</p>
     * @param d value to multiply all entries by
     * @return for convenience, return {@code this}
     * @ if {@code d} is {@code NULL}.
     */
    Field_Vector<T> map_multiply_to_self(T d) ;

    /**
     * Map a division operation to each entry.
     * @param d value to divide all entries by
     * @return {@code this / d}
     * @ if {@code d} is {@code NULL}.
     * @Math_Runtime_Exception if {@code d} is zero.
     */
    Field_Vector<T> map_divide(T d)
        , Math_Runtime_Exception;

    /**
     * Map a division operation to each entry.
     * <p>The instance <strong>is</strong> changed by this method.</p>
     * @param d value to divide all entries by
     * @return for convenience, return {@code this}
     * @ if {@code d} is {@code NULL}.
     * @Math_Runtime_Exception if {@code d} is zero.
     */
    Field_Vector<T> map_divide_to_self(T d)
        , Math_Runtime_Exception;

    /**
     * Map the 1/x function to each entry.
     * @return a vector containing the result of applying the function to each entry.
     * @Math_Runtime_Exception if one of the entries is zero.
     */
    Field_Vector<T> map_inv() Math_Runtime_Exception;

    /**
     * Map the 1/x function to each entry.
     * <p>The instance <strong>is</strong> changed by this method.</p>
     * @return for convenience, return {@code this}
     * @Math_Runtime_Exception if one of the entries is zero.
     */
    Field_Vector<T> map_inv_to_self() Math_Runtime_Exception;

    /**
     * Element-by-element multiplication.
     * @param v vector by which instance elements must be multiplied
     * @return a vector containing {@code this[i] * v[i]} for all {@code i}
     * @ if {@code v} is not the same size as {@code this}
     */
    Field_Vector<T> ebe_multiply(Field_Vector<T> v)
        ;

    /**
     * Element-by-element division.
     * @param v vector by which instance elements must be divided
     * @return a vector containing {@code this[i] / v[i]} for all {@code i}
     * @ if {@code v} is not the same size as {@code this}
     * @Math_Runtime_Exception if one entry of {@code v} is zero.
     */
    Field_Vector<T> ebe_divide(Field_Vector<T> v)
        , Math_Runtime_Exception;

    /**
     * Compute the dot product.
     * @param v vector with which dot product should be computed
     * @return the scalar dot product of {@code this} and {@code v}
     * @ if {@code v} is not the same size as {@code this}
     */
    T dot_product(Field_Vector<T> v) ;

    /**
     * Find the orthogonal projection of this vector onto another vector.
     * @param v vector onto which {@code this} must be projected
     * @return projection of {@code this} onto {@code v}
     * @ if {@code v} is not the same size as {@code this}
     * @Math_Runtime_Exception if {@code v} is the NULL vector.
     */
    Field_Vector<T> projection(Field_Vector<T> v)
        , Math_Runtime_Exception;

    /**
     * Compute the outer product.
     * @param v vector with which outer product should be computed
     * @return the matrix outer product between instance and v
     */
    Field_Matrix<T> outer_product(Field_Vector<T> v);

    /**
     * Returns the entry in the specified index.
     *
     * @param index Index location of entry to be fetched.
     * @return the vector entry at {@code index}.
     * @ if the index is not valid.
     * @see #set_entry(int, Field_Element)
     */
    T get_entry(const int& index) ;

    /**
     * Set a single element.
     * @param index element index.
     * @param value value for the element.
     * @ if the index is not valid.
     * @see #get_entrystatic_cast<int>(
     */
    void set_entry(const int& index, T value) ;

    /**
     * Returns the size of the vector.
     * @return size
     */
    int get_dimension();

    /**
     * Construct a vector by appending a vector to this vector.
     * @param v vector to append to this one.
     * @return a vector
     */
    Field_Vector<T> append(Field_Vector<T> v);

    /**
     * Construct a vector by appending a T to this vector.
     * @param d T to append.
     * @return a vector
     */
    Field_Vector<T> append(T d);

    /**
     * Get a subvector from consecutive elements.
     * @param index index of first element.
     * @param n number of elements to be retrieved.
     * @return a vector containing n elements.
     * @ if the index is not valid.
     * @ if the number of elements if not positive.
     */
    Field_Vector<T> get_sub_vector(const int& index, int n)
        ;

    /**
     * Set a set of consecutive elements.
     * @param index index of first element to be set.
     * @param v vector containing the values to set.
     * @ if the index is not valid.
     */
    void set_sub_vector(const int& index, Field_Vector<T> v) ;

    /**
     * Set all elements to a single value.
     * @param value single value to set for all elements
     */
    void set(T value);

    /**
     * Convert the vector to a T array.
     * <p>The array is independent from vector data, it's elements
     * are copied.</p>
     * @return array containing a copy of vector elements
     */
    std::vector<T> to_array();

}


