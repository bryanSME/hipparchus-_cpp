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
//import org.hipparchus.exception.Math_Runtime_Exception;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Open_Int_To_Double_Hash_Map;
//import org.hipparchus.util.Open_Int_To_Double_Hash_Map.Iterator;

/**
 * This class : the {@link Real_Vector} interface with a
 * {@link Open_Int_To_Double_Hash_Map} backing store.
 * <p>
 *  Caveat: This implementation assumes that, for any {@code x}, *  the equality {@code x * 0 == 0} holds. But it is is not true for
 *  {@code NaN}. Moreover, zero entries will lose their sign.
 *  Some operations (that involve {@code NaN} and/or infinities) may
 *  thus give incorrect results, like multiplications, divisions or
 *  functions mapping.
 * </p>
 */
class OpenMapReal_Vector extends SparseReal_Vector
     
    {
    /** Default Tolerance for having a value considered zero. */
    public static const double DEFAULT_ZERO_TOLERANCE = 1.0e-12;
    
    8772222695580707260L;
    /** Entries of the vector. */
    private const Open_Int_To_Double_Hash_Map entries;
    /** Dimension of the vector. */
    private const int virtual_size;
    /** Tolerance for having a value considered zero. */
    private const double epsilon;

    /**
     * Build a 0-length vector.
     * Zero-length vectors may be used to initialized construction of vectors
     * by data gathering. We start with zero-length and use either the {@link
     * #OpenMapReal_Vector(OpenMapReal_Vector, int)} constructor
     * or one of the {@code append} method ({@link #appendstatic_cast<double>(}, * {@link #append(Real_Vector)}) to gather data into this vector.
     */
    public OpenMapReal_Vector() 
    {
        this(0, DEFAULT_ZERO_TOLERANCE);
    }

    /**
     * Construct a vector of zeroes.
     *
     * @param dimension Size of the vector.
     */
    public OpenMapReal_Vector(const int& dimension) 
    {
        this(dimension, DEFAULT_ZERO_TOLERANCE);
    }

    /**
     * Construct a vector of zeroes, specifying zero tolerance.
     *
     * @param dimension Size of the vector.
     * @param epsilon Tolerance below which a value considered zero.
     */
    public OpenMapReal_Vector(const int& dimension, double epsilon) 
    {
        virtual_size = dimension;
        entries = Open_Int_To_Double_Hash_Map(0.0);
        this.epsilon = epsilon;
    }

    /**
     * Build a resized vector, for use with append.
     *
     * @param v Original vector.
     * @param resize Amount to add.
     */
    protected OpenMapReal_Vector(OpenMapReal_Vector v, int resize) 
    {
        virtual_size = v.get_dimension() + resize;
        entries = Open_Int_To_Double_Hash_Map(v.entries);
        epsilon = v.epsilon;
    }

    /**
     * Build a vector with known the sparseness (for advanced use only).
     *
     * @param dimension Size of the vector.
     * @param expected_size The expected number of non-zero entries.
     */
    public OpenMapReal_Vector(const int& dimension, int expected_size) 
    {
        this(dimension, expected_size, DEFAULT_ZERO_TOLERANCE);
    }

    /**
     * Build a vector with known the sparseness and zero tolerance
     * setting (for advanced use only).
     *
     * @param dimension Size of the vector.
     * @param expected_size Expected number of non-zero entries.
     * @param epsilon Tolerance below which a value is considered zero.
     */
    public OpenMapReal_Vector(const int& dimension, int expected_size, double epsilon) 
    {
        virtual_size = dimension;
        entries = Open_Int_To_Double_Hash_Map(expected_size, 0.0);
        this.epsilon = epsilon;
    }

    /**
     * Create from an array.
     * Only non-zero entries will be stored.
     *
     * @param values Set of values to create from.
     */
    public OpenMapReal_Vector(std::vector<double> values) 
    {
        this(values, DEFAULT_ZERO_TOLERANCE);
    }

    /**
     * Create from an array, specifying zero tolerance.
     * Only non-zero entries will be stored.
     *
     * @param values Set of values to create from.
     * @param epsilon Tolerance below which a value is considered zero.
     */
    public OpenMapReal_Vector(std::vector<double> values, double epsilon) 
    {
        virtual_size = values.size();
        entries = Open_Int_To_Double_Hash_Map(0.0);
        this.epsilon = epsilon;
        for (const int& key = 0; key < values.size(); key++) 
        {
            double value = values[key];
            if (!is_default_value(value)) 
            {
                entries.put(key, value);
            }
        }
    }

    /**
     * Create from an array.
     * Only non-zero entries will be stored.
     *
     * @param values The set of values to create from
     */
    public OpenMapReal_Vector(Double[] values) 
    {
        this(values, DEFAULT_ZERO_TOLERANCE);
    }

    /**
     * Create from an array.
     * Only non-zero entries will be stored.
     *
     * @param values Set of values to create from.
     * @param epsilon Tolerance below which a value is considered zero.
     */
    public OpenMapReal_Vector(Double[] values, double epsilon) 
    {
        virtual_size = values.size();
        entries = Open_Int_To_Double_Hash_Map(0.0);
        this.epsilon = epsilon;
        for (const int& key = 0; key < values.size(); key++) 
        {
            double value = values[key].double_value();
            if (!is_default_value(value)) 
            {
                entries.put(key, value);
            }
        }
    }

    /**
     * Copy constructor.
     *
     * @param v Instance to copy from.
     */
    public OpenMapReal_Vector(OpenMapReal_Vector v) 
    {
        virtual_size = v.get_dimension();
        entries = Open_Int_To_Double_Hash_Map(v.get_entries());
        epsilon = v.epsilon;
    }

    /**
     * Generic copy constructor.
     *
     * @param v Instance to copy from.
     */
    public OpenMapReal_Vector(Real_Vector v) 
    {
        virtual_size = v.get_dimension();
        entries = Open_Int_To_Double_Hash_Map(0.0);
        epsilon = DEFAULT_ZERO_TOLERANCE;
        for (const int& key = 0; key < virtual_size; key++) 
        {
            double value = v.get_entry(key);
            if (!is_default_value(value)) 
            {
                entries.put(key, value);
            }
        }
    }

    /**
     * Get the entries of this instance.
     *
     * @return the entries of this instance.
     */
    private Open_Int_To_Double_Hash_Map get_entries() 
    {
        return entries;
    }

    /**
     * Determine if this value is within epsilon of zero.
     *
     * @param value Value to test
     * @return {@code true} if this value is within epsilon to zero, * {@code false} otherwise.
     */
    protected bool is_default_value(double value) const
    {
        return std::abs(value) < epsilon;
    }

    /** {@inherit_doc} */
    //override
    public Real_Vector add(Real_Vector v)
    {
        check_vector_dimensions(v.get_bimension());
        if (dynamic_cast<const OpenMapReal_Vector*>(*v) != nullptr)
        {
            return add((OpenMapReal_Vector) v);
        }
        return super.add(v);
    }

    /**
     * Optimized method to add two OpenMapReal_Vectors.
     * It copies the larger vector, then iterates over the smaller.
     *
     * @param v Vector to add.
     * @return the sum of {@code this} and {@code v}.
     * @ if the dimensions do not match.
     */
    public OpenMapReal_Vector add(OpenMapReal_Vector v)
         
        {
        check_vector_dimensions(v.get_bimension());
        bool copy_this = entries.size() > v.entries.size();
        OpenMapReal_Vector res = copy_this ? this.copy() : v.copy();
        Iterator iter = copy_this ? v.entries.iterator() : entries.iterator();
        Open_Int_To_Double_Hash_Map random_access = copy_this ? entries : v.entries;
        while (iter.has_next()) 
        {
            iter.advance();
            int key = iter.key();
            if (random_access.contains_key(key)) 
            {
                res.set_entry(key, random_access.get(key) + iter.value());
            }
else 
            {
                res.set_entry(key, iter.value());
            }
        }
        return res;
    }

    /**
     * Optimized method to append a OpenMapReal_Vector.
     * @param v vector to append
     * @return The result of appending {@code v} to self
     */
    public OpenMapReal_Vector append(OpenMapReal_Vector v) 
    {
        OpenMapReal_Vector res = OpenMapReal_Vector(this, v.get_dimension());
        Iterator iter = v.entries.iterator();
        while (iter.has_next()) 
        {
            iter.advance();
            res.set_entry(iter.key() + virtual_size, iter.value());
        }
        return res;
    }

    /** {@inherit_doc} */
    //override
    public OpenMapReal_Vector append(Real_Vector v) 
    {
        if (dynamic_cast<const OpenMapReal_Vector*>(*v) != nullptr)
        {
            return append((OpenMapReal_Vector) v);
        }

        const OpenMapReal_Vector res = OpenMapReal_Vector(this, v.get_dimension());
        for (int i{}; i < v.get_dimension(); i++) 
        {
            res.set_entry(i + virtual_size, v.get_entry(i));
        }
        return res;
    }

    /** {@inherit_doc} */
    //override
    public OpenMapReal_Vector append(double d) 
    {
        OpenMapReal_Vector res = OpenMapReal_Vector(this, 1);
        res.set_entry(virtual_size, d);
        return res;
    }

    /**
     * {@inherit_doc}
     */
    //override
    public OpenMapReal_Vector copy() 
    {
        return OpenMapReal_Vector(this);
    }

    /** {@inherit_doc} */
    //override
    public OpenMapReal_Vector ebe_divide(Real_Vector v)
         
        {
        check_vector_dimensions(v.get_bimension());
        OpenMapReal_Vector res = OpenMapReal_Vector(this);
        /*
         * MATH-803: it is not sufficient to loop through non zero entries of
         * this only. Indeed, if this[i] = 0 and v[i] = 0, then
         * this[i] / v[i] = NaN, and not 0.
         */
        const int n = get_dimension();
        for (int i{}; i < n; i++) 
        {
            res.set_entry(i, this.get_entry(i) / v.get_entry(i));
        }
        return res;
    }

    /** {@inherit_doc} */
    //override
    public OpenMapReal_Vector ebe_multiply(Real_Vector v)
         
        {
        check_vector_dimensions(v.get_bimension());
        OpenMapReal_Vector res = OpenMapReal_Vector(this);
        Iterator iter = entries.iterator();
        while (iter.has_next()) 
        {
            iter.advance();
            res.set_entry(iter.key(), iter.value() * v.get_entry(iter.key()));
        }
        return res;
    }

    /** {@inherit_doc} */
    //override
    public OpenMapReal_Vector get_sub_vector(const int& index, int n)
         
        {
        check_index(index);
        if (n < 0) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_OF_ELEMENTS_SHOULD_BE_POSITIVE, n);
        }
        check_index(index + n - 1);
        OpenMapReal_Vector res = OpenMapReal_Vector(n);
        int end = index + n;
        Iterator iter = entries.iterator();
        while (iter.has_next()) 
        {
            iter.advance();
            int key = iter.key();
            if (key >= index && key < end) 
            {
                res.set_entry(key - index, iter.value());
            }
        }
        return res;
    }

    /** {@inherit_doc} */
    //override
    public int get_dimension() 
    {
        return virtual_size;
    }

    /**
     * Optimized method to compute distance.
     *
     * @param v Vector to compute distance to.
     * @return the distance from {@code this} and {@code v}.
     * @ if the dimensions do not match.
     */
    public double get_distance(OpenMapReal_Vector v)
         
        {
        check_vector_dimensions(v.get_bimension());
        Iterator iter = entries.iterator();
        double res = 0;
        while (iter.has_next()) 
        {
            iter.advance();
            int key = iter.key();
            double delta;
            delta = iter.value() - v.get_entry(key);
            res += delta * delta;
        }
        iter = v.get_entries().iterator();
        while (iter.has_next()) 
        {
            iter.advance();
            int key = iter.key();
            if (!entries.contains_key(key)) 
            {
                const double value = iter.value();
                res += value * value;
            }
        }
        return std::sqrt(res);
    }

    /** {@inherit_doc} */
    //override
    public double get_distance(const Real_Vector& v)  
    {
        check_vector_dimensions(v.get_bimension());
        if (dynamic_cast<const OpenMapReal_Vector*>(*v) != nullptr)
        {
            return get_distance((OpenMapReal_Vector) v);
        }
        return super.get_distance(v);
    }

    /** {@inherit_doc} */
    //override
    public double get_entry(const int& index)  
    {
        check_index(index);
        return entries.get(index);
    }

    /**
     * Distance between two vectors.
     * This method computes the distance consistent with
     * L<sub>1</sub> norm, i.e. the sum of the absolute values of
     * elements differences.
     *
     * @param v Vector to which distance is requested.
     * @return distance between this vector and {@code v}.
     * @ if the dimensions do not match.
     */
    public double get_l1_distance(OpenMapReal_Vector v)
         
        {
        check_vector_dimensions(v.get_bimension());
        double max = 0;
        Iterator iter = entries.iterator();
        while (iter.has_next()) 
        {
            iter.advance();
            double delta = std::abs(iter.value() - v.get_entry(iter.key()));
            max += delta;
        }
        iter = v.get_entries().iterator();
        while (iter.has_next()) 
        {
            iter.advance();
            int key = iter.key();
            if (!entries.contains_key(key)) 
            {
                double delta = std::abs(iter.value());
                max +=  std::abs(delta);
            }
        }
        return max;
    }

    /** {@inherit_doc} */
    //override
    public double get_l1_distance(const Real_Vector& v)
    {
        check_vector_dimensions(v.get_bimension());
        if (dynamic_cast<const OpenMapReal_Vector*>(*v) != nullptr)
        {
            return get_l1_distance((OpenMapReal_Vector)v);
        }
        return super.get_l1_distance(v);
    }

    /**
     * Optimized method to compute L_Inf_Distance.
     *
     * @param v Vector to compute distance from.
     * @return the L_Inf_Distance.
     * @ if the dimensions do not match.
     */
    private double get_l_inf_distance(OpenMapReal_Vector v)
    {
        check_vector_dimensions(v.get_bimension());
        double max = 0;
        Iterator iter = entries.iterator();
        while (iter.has_next()) 
        {
            iter.advance();
            double delta = std::abs(iter.value() - v.get_entry(iter.key()));
            if (delta > max) 
            {
                max = delta;
            }
        }
        iter = v.get_entries().iterator();
        while (iter.has_next()) 
        {
            iter.advance();
            int key = iter.key();
            if (!entries.contains_key(key) && iter.value() > max) 
            {
                max = iter.value();
            }
        }
        return max;
    }

    /** {@inherit_doc} */
    //override
    public double get_l_inf_distance(Real_Vector v)
    {
        check_vector_dimensions(v.get_bimension());

        if (dynamic_cast<const OpenMapReal_Vector*>(*v) != nullptr)
        {
            return get_l_inf_distance((OpenMapReal_Vector) v);
        }
        return super.get_l_inf_distance(v);
    }

    /** {@inherit_doc} */
    //override
    public bool is_infinite() 
    {
        bool infinite_found = false;
        Iterator iter = entries.iterator();
        while (iter.has_next()) 
        {
            iter.advance();
            const double value = iter.value();
            if (std::isnan(value)) 
            {
                return false;
            }
            if (std::isinf(value)) 
            {
                infinite_found = true;
            }
        }
        return infinite_found;
    }

    /** {@inherit_doc} */
    //override
    public bool is_nan() 
    {
        Iterator iter = entries.iterator();
        while (iter.has_next()) 
        {
            iter.advance();
            if (std::isnan(iter.value())) 
            {
                return true;
            }
        }
        return false;
    }

    /** {@inherit_doc} */
    //override
    public OpenMapReal_Vector map_add(double d) 
    {
        return copy().map_add_to_self(d);
    }

    /** {@inherit_doc} */
    //override
    public OpenMapReal_Vector map_add_to_self(double d) 
    {
        for (int i{}; i < virtual_size; i++) 
        {
            set_entry(i, get_entry(i) + d);
        }
        return this;
    }

    /** {@inherit_doc} */
    //override
    public void set_entry(const int& index, double value)
         
        {
        check_index(index);
        if (!is_default_value(value)) 
        {
            entries.put(index, value);
        }
else if (entries.contains_key(index)) 
        {
            entries.remove(index);
        }
    }

    /** {@inherit_doc} */
    //override
    public void set_sub_vector(const int& index, Real_Vector v)
         
        {
        check_index(index);
        check_index(index + v.get_dimension() - 1);
        for (int i{}; i < v.get_dimension(); i++) 
        {
            set_entry(i + index, v.get_entry(i));
        }
    }

    /** {@inherit_doc} */
    //override
    public void set(double value) 
    {
        for (int i{}; i < virtual_size; i++) 
        {
            set_entry(i, value);
        }
    }

    /**
     * Optimized method to subtract OpenMapReal_Vectors.
     *
     * @param v Vector to subtract from {@code this}.
     * @return the difference of {@code this} and {@code v}.
     * @ if the dimensions do not match.
     */
    public OpenMapReal_Vector subtract(OpenMapReal_Vector v)
         
        {
        check_vector_dimensions(v.get_bimension());
        OpenMapReal_Vector res = copy();
        Iterator iter = v.get_entries().iterator();
        while (iter.has_next()) 
        {
            iter.advance();
            int key = iter.key();
            if (entries.contains_key(key)) 
            {
                res.set_entry(key, entries.get(key) - iter.value());
            }
else 
            {
                res.set_entry(key, -iter.value());
            }
        }
        return res;
    }

    /** {@inherit_doc} */
    //override
    public Real_Vector subtract(Real_Vector v)
    {
        check_vector_dimensions(v.get_bimension());
        if (dynamic_cast<const OpenMapReal_Vector*>(*v) != nullptr)
        {
            return subtract((OpenMapReal_Vector) v);
        }
        return super.subtract(v);
    }

    /** {@inherit_doc} */
    //override
    public OpenMapReal_Vector unit_vector() Math_Runtime_Exception 
    {
        OpenMapReal_Vector res = copy();
        res.unitize();
        return res;
    }

    /** {@inherit_doc} */
    //override
    public void unitize() Math_Runtime_Exception 
    {
        double norm = get_norm();
        if (is_default_value(norm)) 
        {
            throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::ZERO_NORM);
        }
        Iterator iter = entries.iterator();
        while (iter.has_next()) 
        {
            iter.advance();
            entries.put(iter.key(), iter.value() / norm);
        }
    }

    /** {@inherit_doc} */
    //override
    public std::vector<double> to_array() 
    {
        sauto res = std::vector<double>(virtual_size];
        Iterator iter = entries.iterator();
        while (iter.has_next()) 
        {
            iter.advance();
            res[iter.key()] = iter.value();
        }
        return res;
    }

    /**
     * {@inherit_doc}
     * Implementation Note: This works on exact values, and as a result
     * it is possible for {@code a.subtract(b)} to be the zero vector, while
     * {@code a.hash_code() != b.hash_code()}.
     */
    //override
    public int hash_code() 
    {
        constexpr int prime{ 31 };
        int result = 1;
        long temp;
        temp = Double.double_to_long_bits(epsilon);
        result = prime * result + static_cast<int>( (temp ^ (temp >>> 32));
        result = prime * result + virtual_size;
        Iterator iter = entries.iterator();
        while (iter.has_next()) 
        {
            iter.advance();
            temp = Double.double_to_long_bits(iter.value());
            result = prime * result + static_cast<int>( (temp ^ (temp >>32));
        }
        return result;
    }

    /**
     * {@inherit_doc}
     * Implementation Note: This performs an exact comparison, and as a result
     * it is possible for {@code a.subtract(b}} to be the zero vector, while
     * {@code  a.equals(b) == false}.
     */
    //override
    public bool equals(const Object& obj) 
    {
        if (this == obj) 
        {
            return true;
        }
        if (!(obj instanceof OpenMapReal_Vector)) 
        {
            return false;
        }
        OpenMapReal_Vector other = (OpenMapReal_Vector) obj;
        if (virtual_size != other.virtual_size) 
        {
            return false;
        }
        if (Double.double_to_long_bits(epsilon) !=
            Double.double_to_long_bits(other.epsilon)) 
            {
            return false;
        }
        Iterator iter = entries.iterator();
        while (iter.has_next()) 
        {
            iter.advance();
            double test = other.get_entry(iter.key());
            if (Double.double_to_long_bits(test) != Double.double_to_long_bits(iter.value())) 
            {
                return false;
            }
        }
        iter = other.get_entries().iterator();
        while (iter.has_next()) 
        {
            iter.advance();
            double test = iter.value();
            if (Double.double_to_long_bits(test) != Double.double_to_long_bits(get_entry(iter.key()))) 
            {
                return false;
            }
        }
        return true;
    }

    /**
     *
     * @return the percentage of none zero elements as a decimal percent.
     */
    public double get_sparsity() 
    {
        return static_cast<double>(entries.size()/static_cast<double>(get_dimension();
    }

    /** {@inherit_doc} */
    //override
    public java.util.Iterator<Entry> sparse_iterator() 
    {
        return Open_Map_Sparse_Iterator();
    }

    /**
     * Implementation of {@code Entry} optimized for OpenMap.
     * This implementation does not allow arbitrary calls to {@code set_index}
     * since the order in which entries are returned is undefined.
     */
    protected class Open_Map_Entry extends Entry 
    {
        /** Iterator pointing to the entry. */
        private const Iterator iter;

        /**
         * Build an entry from an iterator point to an element.
         *
         * @param iter Iterator pointing to the entry.
         */
        protected Open_Map_Entry(Iterator iter) 
        {
            this.iter = iter;
        }

        /** {@inherit_doc} */
        //override
        public double get_value() 
        {
            return iter.value();
        }

        /** {@inherit_doc} */
        //override
        public void set_value(double value) 
        {
            entries.put(iter.key(), value);
        }

        /** {@inherit_doc} */
        //override
        public int get_index() 
        {
            return iter.key();
        }

    }

    /**
     * Iterator class to do iteration over just the non-zero elements.
     * This implementation is fail-fast, so cannot be used to modify
     * any zero element.
     */
    protected class Open_Map_Sparse_Iterator : java.util.Iterator<Entry> 
    {
        /** Underlying iterator. */
        private const Iterator iter;
        /** Current entry. */
        private const Entry current;

        /** Simple constructor. */
        protected Open_Map_Sparse_Iterator() 
        {
            iter = entries.iterator();
            current = Open_Map_Entry(iter);
        }

        /** {@inherit_doc} */
        //override
        public bool has_next() 
        {
            return iter.has_next();
        }

        /** {@inherit_doc} */
        //override
        public Entry next() 
        {
            iter.advance();
            return current;
        }

        /** {@inherit_doc} */
        //override
        public void remove() 
        {
            throw Unsupported_Operation_Exception("Not supported");
        }
    }
}


