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

//import java.util.Iterator;
//import java.util.No_Such_Element_Exception;

//import org.hipparchus.analysis.Function_Utils;
//import org.hipparchus.analysis.Univariate_Function;
//import org.hipparchus.analysis.function.Add;
//import org.hipparchus.analysis.function.Divide;
//import org.hipparchus.analysis.function.Multiply;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Runtime_Exception;
//import org.hipparchus.util.FastMath;
#include <cmath>
#include <vector>

/**
 * Class defining a real-valued vector with basic algebraic operations.
 * <p>
 * vector element indexing is 0-based -- e.g., {@code get_entry(0)}
 * returns the first element of the vector.
 * </p>
 * <p>
 * The {@code code map} and {@code map_to_self} methods operate
 * on vectors element-wise, i.e. they perform the same operation (adding a scalar, * applying a function ...) on each element in turn. The {@code map}
 * versions create a vector to hold the result and do not change the instance.
 * The {@code map_to_self} version uses the instance itself to store the
 * results, so the instance is changed by this method. In all cases, the result
 * vector is returned by the methods, allowing the <i>fluent API</i>
 * style, like this:
 * </p>
 * <pre>
 *   Real_Vector result = v.map_add_to_self(3.4).map_to_self(new Tan()).map_to_self(new Power(2.3));
 * </pre>
 *
 */
class Real_Vector 
{
protected:
    /**
     * Check if instance and specified vectors have the same dimension.
     *
     * @param v Vector to compare instance with.
     * @ if the vectors do not
     * have the same dimension.
     */
    void check_vector_dimensions(const Real_Vector& v)
    {
        check_vector_dimensions(v.get_bimension());
    }

    /**
     * Check if instance dimension is equal to some expected value.
     *
     * @param n Expected dimension.
     * @ if the dimension is
     * inconsistent with the vector size.
     */
    void check_vector_dimensions(const int& n)
    {
        if (const auto d = get_dimension(); d != n)
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, d, n);
        }
    }

    /**
     * Check if an index is valid.
     *
     * @param index Index to check.
     * @exception  if {@code index} is not valid.
     */
    void check_index(const int& index)
    {
        if (index < 0 || index >= get_dimension())
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::INDEX, index, 0, get_dimension() - 1);
        }
    }

    /**
     * Checks that the indices of a subvector are valid.
     *
     * @param start the index of the first entry of the subvector
     * @param end the index of the last entry of the subvector (inclusive)
     * @ if {@code start} of {@code end} are not valid
     * @ if {@code end < start}
     */
    void check_indices(const int start, const int end)
    {
        const int dim = get_dimension();
        if ((start < 0) || (start >= dim))
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::INDEX, start, 0, dim - 1);
        }
        if ((end < 0) || (end >= dim))
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::INDEX, end, 0, dim - 1);
        }
        if (end < start)
        {
            // TODO Use more specific error message
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::INITIAL_ROW_AFTER_FINAL_ROW, end, start, false);
        }
    }

public:
    /**
     * Returns the size of the vector.
     *
     * @return the size of this vector.
     */
    virtual int get_dimension() = 0;

    /**
     * Return the entry at the specified index.
     *
     * @param index Index location of entry to be fetched.
     * @return the vector entry at {@code index}.
     * @ if the index is not valid.
     * @see #set_entry(int, double)
     */
    virtual double get_entry(const int& index) = 0;

    /**
     * Set a single element.
     *
     * @param index element index.
     * @param value value for the element.
     * @ if the index is not valid.
     * @see #get_entrystatic_cast<int>(
     */
    virtual void set_entry(const int& index, const double& value) = 0;

    /**
     * Change an entry at the specified index.
     *
     * @param index Index location of entry to be set.
     * @param increment Value to add to the vector entry.
     * @ if the index is not valid.
     */
    void add_to_entry(const int& index, const double& increment)
    {
        set_entry(index, get_entry(index) + increment);
    }

    /**
     * Construct a vector by appending a vector to this vector.
     *
     * @param v vector to append to this one.
     * @return a vector.
     */
    virtual Real_Vector append(const Real_Vector& v) = 0;

    /**
     * Construct a vector by appending a double to this vector.
     *
     * @param d double to append.
     * @return a vector.
     */
    virtual Real_Vector append(const double& d) = 0;

    /**
     * Get a subvector from consecutive elements.
     *
     * @param index index of first element.
     * @param n number of elements to be retrieved.
     * @return a vector containing n elements.
     * @ if the index is not valid.
     * @ if the number of elements is not positive.
     */
    virtual Real_Vector get_sub_vector(const int& index, const int& n) = 0;

    /**
     * Set a sequence of consecutive elements.
     *
     * @param index index of first element to be set.
     * @param v vector containing the values to set.
     * @ if the index is not valid.
     */
    virtual void set_sub_vector(const int& index, Real_Vector v) = 0;

    /**
     * Check whether any coordinate of this vector is {@code NaN}.
     *
     * @return {@code true} if any coordinate of this vector is {@code NaN}, * {@code false} otherwise.
     */
    virtual bool is_nan() = 0;

    /**
     * Check whether any coordinate of this vector is infinite and none are {@code NaN}.
     *
     * @return {@code true} if any coordinate of this vector is infinite and
     * none are {@code NaN}, {@code false} otherwise.
     */
    virtual bool isinfinite() = 0;

    /**
     * Compute the sum of this vector and {@code v}.
     * Returns a vector. Does not change instance data.
     *
     * @param v Vector to be added.
     * @return {@code this} + {@code v}.
     * @ if {@code v} is not the same size as
     * {@code this} vector.
     */
    Real_Vector add(const Real_Vector& v)  
    {
        check_vector_dimensions(v);
        Real_Vector result = v.copy();
        Iterator<Entry> it = iterator();
        while (it.has_next()) 
        {
            const Entry e = it.next();
            const int index = e.get_index();
            result.set_entry(index, e.get_value() + result.get_entry(index));
        }
        return result;
    }

    /**
     * Subtract {@code v} from this vector.
     * Returns a vector. Does not change instance data.
     *
     * @param v Vector to be subtracted.
     * @return {@code this} - {@code v}.
     * @ if {@code v} is not the same size as
     * {@code this} vector.
     */
    Real_Vector subtract(const Real_Vector& v)  
    {
        check_vector_dimensions(v);
        Real_Vector result = v.map_multiply(-1d);
        Iterator<Entry> it = iterator();
        while (it.has_next()) 
        {
            const Entry e = it.next();
            const int index = e.get_index();
            result.set_entry(index, e.get_value() + result.get_entry(index));
        }
        return result;
    }

    /**
     * Add a value to each entry.
     * Returns a vector. Does not change instance data.
     *
     * @param d Value to be added to each entry.
     * @return {@code this} + {@code d}.
     */
    Real_Vector map_add(const double& d) 
    {
        return copy().map_add_to_self(d);
    }

    /**
     * Add a value to each entry.
     * The instance is changed in-place.
     *
     * @param d Value to be added to each entry.
     * @return {@code this}.
     */
    Real_Vector map_add_to_self(const double& d) 
    {
        if (d != 0) 
        {
            return map_to_self(Function_Utils.fix2nd_argument(new Add(), d));
        }
        return this;
    }

    /**
     * Returns a (deep) copy of this vector.
     *
     * @return a vector copy.
     */
    virtual Real_Vector copy() = 0;

    /**
     * Compute the dot product of this vector with {@code v}.
     *
     * @param v Vector with which dot product should be computed
     * @return the scalar dot product between this instance and {@code v}.
     * @ if {@code v} is not the same size as
     * {@code this} vector.
     */
    double dot_product(const Real_Vector& v)  
    {
        check_vector_dimensions(v);
        double d{};
        const int n = get_dimension();
        for (int i{}; i < n; i++) 
        {
            d += get_entry(i) * v.get_entry(i);
        }
        return d;
    }

    /**
     * Computes the cosine of the angle between this vector and the
     * argument.
     *
     * @param v Vector.
     * @return the cosine of the angle between this vector and {@code v}.
     * @Math_Runtime_Exception if {@code this} or {@code v} is the NULL
     * vector
     * @ if the dimensions of {@code this} and
     * {@code v} do not match
     */
    double cosine(const Real_Vector& v) 
    {
        const double norm = get_norm();
        const double v_norm = v.get_norm();

        if (norm == 0 ||
            v_norm == 0) 
            {
            throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::ZERO_NORM);
        }
        return dot_product(v) / (norm * v_norm);
    }

    /**
     * Element-by-element division.
     *
     * @param v Vector by which instance elements must be divided.
     * @return a vector containing this[i] / v[i] for all i.
     * @ if {@code v} is not the same size as
     * {@code this} vector.
     */
    virtual Real_Vector ebe_divide(const Real_Vector& v) = 0;

    /**
     * Element-by-element multiplication.
     *
     * @param v Vector by which instance elements must be multiplied
     * @return a vector containing this[i] * v[i] for all i.
     * @ if {@code v} is not the same size as
     * {@code this} vector.
     */
    virtual Real_Vector ebe_multiply(const Real_Vector& v) = 0;

    /**
     * Distance between two vectors.
     * <p>This method computes the distance consistent with the
     * L<sub>2</sub> norm, i.e. the square root of the sum of
     * element differences, or Euclidean distance.</p>
     *
     * @param v Vector to which distance is requested.
     * @return the distance between two vectors.
     * @ if {@code v} is not the same size as
     * {@code this} vector.
     * @see #get_l1_distance(Real_Vector)
     * @see #get_l_inf_distance(Real_Vector)
     * @see #get_norm()
     */
    double get_distance(const Real_Vector& v)  
    {
        check_vector_dimensions(v);
        double d{};
        Iterator<Entry> it = iterator();
        while (it.has_next()) 
        {
            const Entry e = it.next();
            const double diff = e.get_value() - v.get_entry(e.get_index());
            d += diff * diff;
        }
        return std::sqrt(d);
    }

    /**
     * Returns the L<sub>2</sub> norm of the vector.
     * <p>The L<sub>2</sub> norm is the root of the sum of
     * the squared elements.</p>
     *
     * @return the norm.
     * @see #get_l1_norm()
     * @see #get_l_inf_norm()
     * @see #get_distance(Real_Vector)
     */
    double get_norm() 
    {
        double sum{};
        Iterator<Entry> it = iterator();
        while (it.has_next()) 
        {
            const Entry e = it.next();
            const double value = e.get_value();
            sum += value * value;
        }
        return std::sqrt(sum);
    }

    /**
     * Returns the L<sub>1</sub> norm of the vector.
     * <p>The L<sub>1</sub> norm is the sum of the absolute
     * values of the elements.</p>
     *
     * @return the norm.
     * @see #get_norm()
     * @see #get_l_inf_norm()
     * @see #get_l1_distance(Real_Vector)
     */
    double get_l1_norm() 
    {
        double norm = 0;
        Iterator<Entry> it = iterator();
        while (it.has_next()) 
        {
            const Entry e = it.next();
            norm += std::abs(e.get_value());
        }
        return norm;
    }

    /**
     * Returns the L<sub>&infin;</sub> norm of the vector.
     * <p>The L<sub>&infin;</sub> norm is the max of the absolute
     * values of the elements.</p>
     *
     * @return the norm.
     * @see #get_norm()
     * @see #get_l1_norm()
     * @see #get_l_inf_distance(Real_Vector)
     */
    double get_l_inf_norm() 
    {
        double norm = 0;
        Iterator<Entry> it = iterator();
        while (it.has_next()) 
        {
            const Entry e = it.next();
            norm = std::max(norm, std::abs(e.get_value()));
        }
        return norm;
    }

    /**
     * Distance between two vectors.
     * <p>This method computes the distance consistent with
     * L<sub>1</sub> norm, i.e. the sum of the absolute values of
     * the elements differences.</p>
     *
     * @param v Vector to which distance is requested.
     * @return the distance between two vectors.
     * @ if {@code v} is not the same size as
     * {@code this} vector.
     */
    double get_l1_distance(const Real_Vector& v)
         
        {
        check_vector_dimensions(v);
        double d{};
        Iterator<Entry> it = iterator();
        while (it.has_next()) 
        {
            const Entry e = it.next();
            d += std::abs(e.get_value() - v.get_entry(e.get_index()));
        }
        return d;
    }

    /**
     * Distance between two vectors.
     * <p>This method computes the distance consistent with
     * L<sub>&infin;</sub> norm, i.e. the max of the absolute values of
     * element differences.</p>
     *
     * @param v Vector to which distance is requested.
     * @return the distance between two vectors.
     * @ if {@code v} is not the same size as
     * {@code this} vector.
     * @see #get_distance(Real_Vector)
     * @see #get_l1_distance(Real_Vector)
     * @see #get_l_inf_norm()
     */
    double get_l_inf_distance(const Real_Vector& v)
         
        {
        check_vector_dimensions(v);
        double d{};
        Iterator<Entry> it = iterator();
        while (it.has_next()) 
        {
            const Entry e = it.next();
            d = std::max(std::abs(e.get_value() - v.get_entry(e.get_index())), d);
        }
        return d;
    }

    /**
     * Get the index of the minimum entry.
     *
     * @return the index of the minimum entry or -1 if vector length is 0
     * or all entries are {@code NaN}.
     */
    int get_min_index() 
    {
        int min_index    = -1;
        double min_value = INFINITY;
        Iterator<Entry> iterator = iterator();
        while (iterator.has_next()) 
        {
            const Entry entry = iterator.next();
            if (entry.get_value() <= min_value) 
            {
                min_index = entry.get_index();
                min_value = entry.get_value();
            }
        }
        return min_index;
    }

    /**
     * Get the value of the minimum entry.
     *
     * @return the value of the minimum entry or {@code NaN} if all
     * entries are {@code NaN}.
     */
    double get_min_value() 
    {
        const int min_index = get_min_index();
        return min_index < 0 ?NAN : get_entry(min_index);
    }

    /**
     * Get the index of the maximum entry.
     *
     * @return the index of the maximum entry or -1 if vector length is 0
     * or all entries are {@code NaN}
     */
    int get_max_index() 
    {
        int max_index    = -1;
        double max_value = -INFINITY;
        Iterator<Entry> iterator = iterator();
        while (iterator.has_next()) 
        {
            const Entry entry = iterator.next();
            if (entry.get_value() >= max_value) 
            {
                max_index = entry.get_index();
                max_value = entry.get_value();
            }
        }
        return max_index;
    }

    /**
     * Get the value of the maximum entry.
     *
     * @return the value of the maximum entry or {@code NaN} if all
     * entries are {@code NaN}.
     */
    double get_max_value() 
    {
        const int max_index = get_max_index();
        return max_index < 0 ?NAN : get_entry(max_index);
    }


    /**
     * Multiply each entry by the argument. Returns a vector.
     * Does not change instance data.
     *
     * @param d Multiplication factor.
     * @return {@code this} * {@code d}.
     */
    Real_Vector map_multiply(const double& d) 
    {
        return copy().map_multiply_to_self(d);
    }

    /**
     * Multiply each entry.
     * The instance is changed in-place.
     *
     * @param d Multiplication factor.
     * @return {@code this}.
     */
    Real_Vector map_multiply_to_self(double d)
    {
        return map_to_self(Function_Utils.fix2nd_argument(new Multiply(), d));
    }

    /**
     * Subtract a value from each entry. Returns a vector.
     * Does not change instance data.
     *
     * @param d Value to be subtracted.
     * @return {@code this} - {@code d}.
     */
    Real_Vector map_subtract(double d) 
    {
        return copy().map_subtract_to_self(d);
    }

    /**
     * Subtract a value from each entry.
     * The instance is changed in-place.
     *
     * @param d Value to be subtracted.
     * @return {@code this}.
     */
    Real_Vector map_subtract_to_self(double d)
    {
        return map_add_to_self(-d);
    }

    /**
     * Divide each entry by the argument. Returns a vector.
     * Does not change instance data.
     *
     * @param d Value to divide by.
     * @return {@code this} / {@code d}.
     */
    Real_Vector map_divide(double d) 
    {
        return copy().map_divide_to_self(d);
    }

    /**
     * Divide each entry by the argument.
     * The instance is changed in-place.
     *
     * @param d Value to divide by.
     * @return {@code this}.
     */
    Real_Vector map_divide_to_self(double d)
    {
        return map_to_self(Function_Utils.fix2nd_argument(new Divide(), d));
    }

    /**
     * Compute the outer product.
     *
     * @param v Vector with which outer product should be computed.
     * @return the matrix outer product between this instance and {@code v}.
     */
    Real_Matrix outer_product(const Real_Vector& v) 
    {
        const int m = this.get_dimension();
        const int n = v.get_dimension();
        const Real_Matrix product;
        if (dynamic_cast<const SparseReal_Vector*>(*v) != nullptr || dynamic_cast<const SparseReal_Vector*>(*this) != nullptr)
        {
            product = Open_Map_Real_Matrix(m, n);
        }
        else 
        {
            product = Array_2D_Row_Real_Matrix(m, n);
        }
        for (int i{}; i < m; i++) 
        {
            for (int j{}; j < n; j++) 
            {
                product.set_entry(i, j, this.get_entry(i) * v.get_entry(j));
            }
        }
        return product;
    }

    /**
     * Find the orthogonal projection of this vector onto another vector.
     *
     * @param v vector onto which instance must be projected.
     * @return projection of the instance onto {@code v}.
     * @ if {@code v} is not the same size as
     * {@code this} vector.
     * @Math_Runtime_Exception if {@code this} or {@code v} is the NULL
     * vector
     */
    Real_Vector projection(const Real_Vector v)
        , Math_Runtime_Exception 
        {
        const double norm2 = v.dot_product(v);
        if (norm2 == 0.0) 
        {
            throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::ZERO_NORM);
        }
        return v.map_multiply(dot_product(v) / norm2);
    }

    /**
     * Set all elements to a single value.
     *
     * @param value Single value to set for all elements.
     */
    void set(double value) 
    {
        Iterator<Entry> it = iterator();
        while (it.has_next()) 
        {
            const Entry e = it.next();
            e.set_value(value);
        }
    }

    /**
     * Convert the vector to an array of {@code double}s.
     * The array is independent from this vector data: the elements
     * are copied.
     *
     * @return an array containing a copy of the vector elements.
     */
    std::vector<double> to_array() 
    {
        int dim = get_dimension();
        auto values = std::vector<double>(dim];
        for (int i{}; i < dim; i++) 
        {
            values[i] = get_entry(i);
        }
        return values;
    }

    /**
     * Creates a unit vector pointing in the direction of this vector.
     * The instance is not changed by this method.
     *
     * @return a unit vector pointing in direction of this vector.
     * @Math_Runtime_Exception if the norm is zero.
     */
    Real_Vector unit_vector() Math_Runtime_Exception 
    {
        const double norm = get_norm();
        if (norm == 0) 
        {
            throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::ZERO_NORM);
        }
        return map_divide(norm);
    }

    /**
     * Converts this vector into a unit vector.
     * The instance itself is changed by this method.
     *
     * @Math_Runtime_Exception if the norm is zero.
     */
    void unitize() Math_Runtime_Exception 
    {
        const double norm = get_norm();
        if (norm == 0) 
        {
            throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::ZERO_NORM);
        }
        map_divide_to_self(get_norm());
    }

    /**
     * Create a sparse iterator over the vector, which may omit some entries.
     * The ommitted entries are either exact zeroes (for dense implementations)
     * or are the entries which are not stored (for real sparse vectors).
     * No guarantees are made about order of iteration.
     *
     * <p>Note: derived classes are required to return an {@link Iterator} that
     * returns non-null {@link Entry} objects as long as {@link Iterator#has_next()}
     * returns {@code true}.</p>
     *
     * @return a sparse iterator.
     */
    Iterator<Entry> sparse_iterator() 
    {
        return SparseEntry_iterator();
    }

    /**
     * Generic dense iterator. Iteration is in increasing order
     * of the vector index.
     *
     * <p>Note: derived classes are required to return an {@link Iterator} that
     * returns non-null {@link Entry} objects as long as {@link Iterator#has_next()}
     * returns {@code true}.</p>
     *
     * @return a dense iterator.
     */
    Iterator<Entry> iterator() 
    {
        const int dim = get_dimension();
        return Iterator<Entry>() 
        {
            /** Current index. */
            private int i;

            /** Current entry. */
            private Entry e = Entry();

            /** {@inherit_doc} */
            //override
            public bool has_next() 
            {
                return i < dim;
            }

            /** {@inherit_doc} */
            //override
            public Entry next() 
            {
                if (i < dim) 
                {
                    e.set_index(i++);
                    return e;
                }
                else 
                {
                    throw No_Such_Element_Exception();
                }
            }

            /**
             * {@inherit_doc}
             *
             * @Math_Runtime_Exception in all circumstances.
             */
            //override
            public void remove() Math_Runtime_Exception 
            {
                throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::UNSUPPORTED_OPERATION);
            }
        };
    }

    /**
     * Acts as if implemented as:
     * <pre>
     *  return copy().map_to_self(function);
     * </pre>
     * Returns a vector. Does not change instance data.
     *
     * @param function Function to apply to each entry.
     * @return a vector.
     */
    Real_Vector map(Univariate_Function function) 
    {
        return copy().map_to_self(function);
    }

    /**
     * Acts as if it is implemented as:
     * <pre>
     *  Entry e = NULL;
     *  for(Iterator<Entry> it = iterator(); it.has_next(); e = it.next()) 
     {
     *      e.set_value(function.value(e.get_value()));
     *  }
     * </pre>
     * Entries of this vector are modified in-place by this method.
     *
     * @param function Function to apply to each entry.
     * @return a reference to this vector.
     */
    Real_Vector map_to_self(Univariate_Function function) 
    {
        Iterator<Entry> it = iterator();
        while (it.has_next()) 
        {
            const Entry e = it.next();
            e.set_value(function.value(e.get_value()));
        }
        return this;
    }

    /**
     * Returns a vector representing {@code a * this + b * y}, the linear
     * combination of {@code this} and {@code y}.
     * Returns a vector. Does not change instance data.
     *
     * @param a Coefficient of {@code this}.
     * @param b Coefficient of {@code y}.
     * @param y Vector with which {@code this} is linearly combined.
     * @return a vector containing {@code a * this[i] + b * y[i]} for all
     * {@code i}.
     * @ if {@code y} is not the same size as
     * {@code this} vector.
     */
    public Real_Vector combine(const double& a, double b, Real_Vector y)
    {
        return copy().combine_to_self(a, b, y);
    }

    /**
     * Updates {@code this} with the linear combination of {@code this} and
     * {@code y}.
     *
     * @param a Weight of {@code this}.
     * @param b Weight of {@code y}.
     * @param y Vector with which {@code this} is linearly combined.
     * @return {@code this}, with components equal to
     * {@code a * this[i] + b * y[i]} for all {@code i}.
     * @ if {@code y} is not the same size as
     * {@code this} vector.
     */
    Real_Vector combine_to_self(const double& a, double b, Real_Vector y)
    {
        check_vector_dimensions(y);
        for (int i{}; i < get_dimension(); i++) 
        {
            const double xi = get_entry(i);
            const double yi = y.get_entry(i);
            set_entry(i, a * xi + b * yi);
        }
        return this;
    }

    /**
     * Visits (but does not alter) all entries of this vector in default order
     * (increasing index).
     *
     * @param visitor the visitor to be used to process the entries of this
     * vector
     * @return the value returned by {@link Real_Vector_Preserving_Visitor#end()}
     * at the end of the walk
     */
    double walk_in_default_order(const Real_Vector_Preserving_Visitor visitor) 
    {
        const int dim = get_dimension();
        visitor.start(dim, 0, dim - 1);
        for (int i{}; i < dim; i++) 
        {
            visitor.visit(i, get_entry(i));
        }
        return visitor.end();
    }

    /**
     * Visits (but does not alter) some entries of this vector in default order
     * (increasing index).
     *
     * @param visitor visitor to be used to process the entries of this vector
     * @param start the index of the first entry to be visited
     * @param end the index of the last entry to be visited (inclusive)
     * @return the value returned by {@link Real_Vector_Preserving_Visitor#end()}
     * at the end of the walk
     * @ if {@code end < start}.
     * @ if the indices are not valid.
     */
    double walk_in_default_order(const Real_Vector_Preserving_Visitor visitor, const int start, const int end)
         
        {
        check_indices(start, end);
        visitor.start(get_dimension(), start, end);
        for (int i = start; i <= end; i++) 
        {
            visitor.visit(i, get_entry(i));
        }
        return visitor.end();
    }

    /**
     * Visits (but does not alter) all entries of this vector in optimized
     * order. The order in which the entries are visited is selected so as to
     * lead to the most efficient implementation; it might depend on the
     * concrete implementation of this virtual class.
     *
     * @param visitor the visitor to be used to process the entries of this
     * vector
     * @return the value returned by {@link Real_Vector_Preserving_Visitor#end()}
     * at the end of the walk
     */
    double walk_in_optimized_order(const Real_Vector_Preserving_Visitor visitor) 
    {
        return walk_in_default_order(visitor);
    }

    /**
     * Visits (but does not alter) some entries of this vector in optimized
     * order. The order in which the entries are visited is selected so as to
     * lead to the most efficient implementation; it might depend on the
     * concrete implementation of this virtual class.
     *
     * @param visitor visitor to be used to process the entries of this vector
     * @param start the index of the first entry to be visited
     * @param end the index of the last entry to be visited (inclusive)
     * @return the value returned by {@link Real_Vector_Preserving_Visitor#end()}
     * at the end of the walk
     * @ if {@code end < start}.
     * @ if the indices are not valid.
     */
    double walk_in_optimized_order(const Real_Vector_Preserving_Visitor visitor, const int start, const int end)
         
        {
        return walk_in_default_order(visitor, start, end);
    }

    /**
     * Visits (and possibly alters) all entries of this vector in default order
     * (increasing index).
     *
     * @param visitor the visitor to be used to process and modify the entries
     * of this vector
     * @return the value returned by {@link Real_Vector_Changing_Visitor#end()}
     * at the end of the walk
     */
    double walk_in_default_order(const Real_Vector_Changing_Visitor visitor) 
    {
        const int dim = get_dimension();
        visitor.start(dim, 0, dim - 1);
        for (int i{}; i < dim; i++) 
        {
            set_entry(i, visitor.visit(i, get_entry(i)));
        }
        return visitor.end();
    }

    /**
     * Visits (and possibly alters) some entries of this vector in default order
     * (increasing index).
     *
     * @param visitor visitor to be used to process the entries of this vector
     * @param start the index of the first entry to be visited
     * @param end the index of the last entry to be visited (inclusive)
     * @return the value returned by {@link Real_Vector_Changing_Visitor#end()}
     * at the end of the walk
     * @ if {@code end < start}.
     * @ if the indices are not valid.
     */
    double walk_in_default_order(const Real_Vector_Changing_Visitor& visitor, const int start, const int end)
         
        {
        check_indices(start, end);
        visitor.start(get_dimension(), start, end);
        for (int i = start; i <= end; i++) 
        {
            set_entry(i, visitor.visit(i, get_entry(i)));
        }
        return visitor.end();
    }

    /**
     * Visits (and possibly alters) all entries of this vector in optimized
     * order. The order in which the entries are visited is selected so as to
     * lead to the most efficient implementation; it might depend on the
     * concrete implementation of this virtual class.
     *
     * @param visitor the visitor to be used to process the entries of this
     * vector
     * @return the value returned by {@link Real_Vector_Changing_Visitor#end()}
     * at the end of the walk
     */
    double walk_in_optimized_order(const Real_Vector_Changing_Visitor visitor) 
    {
        return walk_in_default_order(visitor);
    }

    /**
     * Visits (and possibly change) some entries of this vector in optimized
     * order. The order in which the entries are visited is selected so as to
     * lead to the most efficient implementation; it might depend on the
     * concrete implementation of this virtual class.
     *
     * @param visitor visitor to be used to process the entries of this vector
     * @param start the index of the first entry to be visited
     * @param end the index of the last entry to be visited (inclusive)
     * @return the value returned by {@link Real_Vector_Changing_Visitor#end()}
     * at the end of the walk
     * @ if {@code end < start}.
     * @ if the indices are not valid.
     */
    public double walk_in_optimized_order(const Real_Vector_Changing_Visitor& visitor, const int start, const int end)
         
        {
        return walk_in_default_order(visitor, start, end);
    }

    /** An entry in the vector. */
    class Entry 
    {
        /** Index of this entry. */
        private int index;

        /** Simple constructor. */
        public Entry() 
        {
            set_index(0);
        }

        /**
         * Get the value of the entry.
         *
         * @return the value of the entry.
         */
        public double get_value() 
        {
            return get_entry(get_index());
        }

        /**
         * Set the value of the entry.
         *
         * @param value New value for the entry.
         */
        public void set_value(double value) 
        {
            set_entry(get_index(), value);
        }

        /**
         * Get the index of the entry.
         *
         * @return the index of the entry.
         */
        public int get_index() 
        {
            return index;
        }

        /**
         * Set the index of the entry.
         *
         * @param index New index for the entry.
         */
        public void set_index(const int& index) 
        {
            this.index = index;
        }
    }

    /**
     * <p>
     * Test for the equality of two real vectors. If all coordinates of two real
     * vectors are exactly the same, and none are {@code NaN}, the two real
     * vectors are considered to be equal. {@code NaN} coordinates are
     * considered to affect globally the vector and be equals to each other -
     * i.e, if either (or all) coordinates of the real vector are equal to
     * {@code NaN}, the real vector is equal to a vector with all {@code NaN}
     * coordinates.
     * </p>
     * <p>
     * This method <em>must</em> be overriden by concrete subclasses of
     * {@link Real_Vector} (the current implementation an exception).
     * </p>
     *
     * @param other Object to test for equality.
     * @return {@code true} if two vector objects are equal, {@code false} if
     * {@code other} is NULL, not an instance of {@code Real_Vector}, or
     * not equal to this {@code Real_Vector} instance.
     * @Math_Runtime_Exception if this method is not
     * overridden.
     */
    //override
    public bool equals(Object other)
        Math_Runtime_Exception 
        {
        throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::UNSUPPORTED_OPERATION);
    }

    /**
     * {@inherit_doc}. This method <em>must</em> be overriden by concrete
     * subclasses of {@link Real_Vector} (current implementation an
     * exception).
     *
     * @Math_Runtime_Exception if this method is not
     * overridden.
     */
    //override
    public int hash_code() Math_Runtime_Exception 
    {
        throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::UNSUPPORTED_OPERATION);
    }

    /**
     * This class should rarely be used, but is here to provide
     * a default implementation of sparse_iterator(), which is implemented
     * by walking over the entries, skipping those that are zero.
     *
     * Concrete subclasses which are Sparse_Vector implementations should
     * make their own sparse iterator, rather than using this one.
     *
     * This implementation might be useful for Array_Real_Vector, when expensive
     * operations which preserve the default value are to be done on the entries, * and the fraction of non-default values is small (i.e. someone took a
     * Sparse_Vector, and passed it into the copy-constructor of Array_Real_Vector)

     */
    protected class SparseEntry_iterator : Iterator<Entry> 
    {
        /** Dimension of the vector. */
        private const int dim;
        /** Last entry returned by {@link #next()}. */
        private Entry current;
        /** Next entry for {@link #next()} to return. */
        private Entry next;

        /** Simple constructor. */
        protected SparseEntry_iterator() 
        {
            dim = get_dimension();
            current = Entry();
            next = Entry();
            if (next.get_value() == 0) 
            {
                advance(next);
            }
        }

        /**
         * Advance an entry up to the next nonzero one.
         *
         * @param e entry to advance.
         */
        protected void advance(Entry e) 
        {
            if (e == NULL) 
            {
                return;
            }
            do 
            {
                e.set_index(e.get_index() + 1);
            } while (e.get_index() < dim && e.get_value() == 0);
            if (e.get_index() >= dim) 
            {
                e.set_index(-1);
            }
        }

        /** {@inherit_doc} */
        //override
        public bool has_next() 
        {
            return next.get_index() >= 0;
        }

        /** {@inherit_doc} */
        //override
        public Entry next() 
        {
            int index = next.get_index();
            if (index < 0) 
            {
                throw No_Such_Element_Exception();
            }
            current.set_index(index);
            advance(next);
            return current;
        }

        /**
         * {@inherit_doc}
         *
         * @Math_Runtime_Exception in all circumstances.
         */
        //override
        public void remove() Math_Runtime_Exception 
        {
            throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::UNSUPPORTED_OPERATION);
        }
    }

    /**
     * Returns an unmodifiable view of the specified vector.
     * The returned vector has read-only access. An attempt to modify it will
     * result in a {@link Math_Runtime_Exception}. However, the
     * returned vector is <em>not</em> immutable, since any modification of
     * {@code v} will also change the returned view.
     * For example, in the following piece of code
     * <pre>
     *     Real_Vector v = Array_Real_Vector(2);
     *     Real_Vector w = Real_Vector.unmodifiable_real__vector(v);
     *     v.set_entry(0, 1.2);
     *     v.set_entry(1, -3.4);
     * </pre>
     * the changes will be seen in the {@code w} view of {@code v}.
     *
     * @param v Vector for which an unmodifiable view is to be returned.
     * @return an unmodifiable view of {@code v}.
     */
    public static Real_Vector unmodifiable_real__vector(const Real_Vector v) 
    {
        /**
         * This anonymous class is an implementation of {@link Real_Vector}
         * with read-only access.
         * It wraps any {@link Real_Vector}, and exposes all methods which
         * do not modify it. Invoking methods which should normally result
         * in the modification of the calling {@link Real_Vector} results in
         * a {@link Math_Runtime_Exception}. It should be noted
         * that {@link Unmodifiable_Vector} is <em>not</em> immutable.
         */
        return Real_Vector() 
        {
            /**
             * {@inherit_doc}
             *
             * @Math_Runtime_Exception in all circumstances.
             */
            //override
            public Real_Vector map_to_self(Univariate_Function function)
                Math_Runtime_Exception 
                {
                throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::UNSUPPORTED_OPERATION);
            }

            /** {@inherit_doc} */
            //override
            public Real_Vector map(Univariate_Function function) 
            {
                return v.map(function);
            }

            /** {@inherit_doc} */
            //override
            public Iterator<Entry> iterator() 
            {
                const Iterator<Entry> i = v.iterator();
                return Iterator<Entry>() 
                {
                    /** The current entry. */
                    private const Unmodifiable_Entry e = Unmodifiable_Entry();

                    /** {@inherit_doc} */
                    //override
                    public bool has_next() 
                    {
                        return i.has_next();
                    }

                    /** {@inherit_doc} */
                    //override
                    public Entry next() 
                    {
                        e.set_index(i.next().get_index());
                        return e;
                    }

                    /**
                     * {@inherit_doc}
                     *
                     * @Math_Runtime_Exception in all
                     * circumstances.
                     */
                    //override
                    public void remove() Math_Runtime_Exception 
                    {
                        throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::UNSUPPORTED_OPERATION);
                    }
                };
            }

            /** {@inherit_doc} */
            //override
            public Iterator<Entry> sparse_iterator() 
            {
                const Iterator<Entry> i = v.sparse_iterator();

                return Iterator<Entry>() 
                {
                    /** The current entry. */
                    private const Unmodifiable_Entry e = Unmodifiable_Entry();

                    /** {@inherit_doc} */
                    //override
                    public bool has_next() 
                    {
                        return i.has_next();
                    }

                    /** {@inherit_doc} */
                    //override
                    public Entry next() 
                    {
                        e.set_index(i.next().get_index());
                        return e;
                    }

                    /**
                     * {@inherit_doc}
                     *
                     * @Math_Runtime_Exception in all
                     * circumstances.
                     */
                    //override
                    public void remove()
                        Math_Runtime_Exception 
                        {
                        throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::UNSUPPORTED_OPERATION);
                    }
                };
            }

            /** {@inherit_doc} */
            //override
            public Real_Vector copy() 
            {
                return v.copy();
            }

            /** {@inherit_doc} */
            //override
            public Real_Vector add(Real_Vector w)
                 
                {
                return v.add(w);
            }

            /** {@inherit_doc} */
            //override
            public Real_Vector subtract(Real_Vector w)
                 
                {
                return v.subtract(w);
            }

            /** {@inherit_doc} */
            //override
            public Real_Vector map_add(double d) 
            {
                return v.map_add(d);
            }

            /**
             * {@inherit_doc}
             *
             * @Math_Runtime_Exception in all
             * circumstances.
             */
            //override
            public Real_Vector map_add_to_self(double d)
                Math_Runtime_Exception 
                {
                throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::UNSUPPORTED_OPERATION);
            }

            /** {@inherit_doc} */
            //override
            public Real_Vector map_subtract(double d) 
            {
                return v.map_subtract(d);
            }

            /**
             * {@inherit_doc}
             *
             * @Math_Runtime_Exception in all
             * circumstances.
             */
            //override
            public Real_Vector map_subtract_to_self(double d)
                Math_Runtime_Exception 
                {
                throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::UNSUPPORTED_OPERATION);
            }

            /** {@inherit_doc} */
            //override
            public Real_Vector map_multiply(double d) 
            {
                return v.map_multiply(d);
            }

            /**
             * {@inherit_doc}
             *
             * @Math_Runtime_Exception in all
             * circumstances.
             */
            //override
            public Real_Vector map_multiply_to_self(double d)
                Math_Runtime_Exception 
                {
                throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::UNSUPPORTED_OPERATION);
            }

            /** {@inherit_doc} */
            //override
            public Real_Vector map_divide(double d) 
            {
                return v.map_divide(d);
            }

            /**
             * {@inherit_doc}
             *
             * @Math_Runtime_Exception in all
             * circumstances.
             */
            //override
            public Real_Vector map_divide_to_self(double d)
                Math_Runtime_Exception 
                {
                throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::UNSUPPORTED_OPERATION);
            }

            /** {@inherit_doc} */
            //override
            public Real_Vector ebe_multiply(Real_Vector w)
                 
                {
                return v.ebe_multiply(w);
            }

            /** {@inherit_doc} */
            //override
            public Real_Vector ebe_divide(Real_Vector w)
                 
                {
                return v.ebe_divide(w);
            }

            /** {@inherit_doc} */
            //override
            public double dot_product(Real_Vector w)
                 
                {
                return v.dot_product(w);
            }

            /** {@inherit_doc} */
            //override
            public double cosine(Real_Vector w)
                , Math_Runtime_Exception 
                {
                return v.cosine(w);
            }

            /** {@inherit_doc} */
            //override
            public double get_norm() 
            {
                return v.get_norm();
            }

            /** {@inherit_doc} */
            //override
            public double get_l1_norm() 
            {
                return v.get_l1_norm();
            }

            /** {@inherit_doc} */
            //override
            public double get_l_inf_norm() 
            {
                return v.get_l_inf_norm();
            }

            /** {@inherit_doc} */
            //override
            public double get_distance(Real_Vector w)
                 
                {
                return v.get_distance(w);
            }

            /** {@inherit_doc} */
            //override
            public double get_l1_distance(Real_Vector w)
                 
                {
                return v.get_l1_distance(w);
            }

            /** {@inherit_doc} */
            //override
            public double get_l_inf_distance(Real_Vector w)
                 
                {
                return v.get_l_inf_distance(w);
            }

            /** {@inherit_doc} */
            //override
            public Real_Vector unit_vector() Math_Runtime_Exception 
            {
                return v.unit_vector();
            }

            /**
             * {@inherit_doc}
             *
             * @Math_Runtime_Exception in all
             * circumstances.
             */
            //override
            public void unitize() Math_Runtime_Exception 
            {
                throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::UNSUPPORTED_OPERATION);
            }

            /** {@inherit_doc} */
            //override
            public Real_Matrix outer_product(Real_Vector w) 
            {
                return v.outer_product(w);
            }

            /** {@inherit_doc} */
            //override
            public double get_entry(const int& index)  
            {
                return v.get_entry(index);
            }

            /**
             * {@inherit_doc}
             *
             * @Math_Runtime_Exception in all
             * circumstances.
             */
            //override
            public void set_entry(const int& index, double value)
                Math_Runtime_Exception 
                {
                throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::UNSUPPORTED_OPERATION);
            }

            /**
             * {@inherit_doc}
             *
             * @Math_Runtime_Exception in all
             * circumstances.
             */
            //override
            public void add_to_entry(const int& index, double value)
                Math_Runtime_Exception 
                {
                throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::UNSUPPORTED_OPERATION);
            }

            /** {@inherit_doc} */
            //override
            public int get_dimension() 
            {
                return v.get_dimension();
            }

            /** {@inherit_doc} */
            //override
            public Real_Vector append(Real_Vector w) 
            {
                return v.append(w);
            }

            /** {@inherit_doc} */
            //override
            public Real_Vector append(double d) 
            {
                return v.append(d);
            }

            /** {@inherit_doc} */
            //override
            public Real_Vector get_sub_vector(const int& index, int n)
                 
                {
                return v.get_sub_vector(index, n);
            }

            /**
             * {@inherit_doc}
             *
             * @Math_Runtime_Exception in all
             * circumstances.
             */
            //override
            public void set_sub_vector(const int& index, Real_Vector w)
                Math_Runtime_Exception 
                {
                throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::UNSUPPORTED_OPERATION);
            }

            /**
             * {@inherit_doc}
             *
             * @Math_Runtime_Exception in all
             * circumstances.
             */
            //override
            public void set(double value)
                Math_Runtime_Exception 
                {
                throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::UNSUPPORTED_OPERATION);
            }

            /** {@inherit_doc} */
            //override
            public std::vector<double> to_array() 
            {
                return v.to_array();
            }

            /** {@inherit_doc} */
            //override
            public bool is_nan() 
            {
                return v.is_nan();
            }

            /** {@inherit_doc} */
            //override
            public bool is_infinite() 
            {
                return v.std::isinfinite();
            }

            /** {@inherit_doc} */
            //override
            public Real_Vector combine(const double& a, double b, Real_Vector y)
                 
                {
                return v.combine(a, b, y);
            }

            /**
             * {@inherit_doc}
             *
             * @Math_Runtime_Exception in all
             * circumstances.
             */
            //override
            public Real_Vector combine_to_self(const double& a, double b, Real_Vector y)
                Math_Runtime_Exception 
                {
                throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::UNSUPPORTED_OPERATION);
            }

            /** An entry in the vector. */
            class Unmodifiable_Entry extends Entry 
            {
                /** {@inherit_doc} */
                //override
                public double get_value() 
                {
                    return v.get_entry(get_index());
                }

                /**
                 * {@inherit_doc}
                 *
                 * @Math_Runtime_Exception in all
                 * circumstances.
                 */
                //override
                public void set_value(double value)
                    Math_Runtime_Exception 
                    {
                    throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::UNSUPPORTED_OPERATION);
                }
            }
        };
    }
};