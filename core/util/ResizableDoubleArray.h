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
//package org.hipparchus.util;

//import java.io.Serializable;
//import java.util.Arrays;

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.exception.;

/**
 * A variable length primitive double array implementation that automatically
 * handles expanding and contracting its internal storage array as elements
 * are added and removed.
 * <p>
 * The internal storage array starts with capacity determined by the
 * {@code initial_capacity} property, which can be set by the constructor.
 * The default initial capacity is 16.  Adding elements using
 * {@link #add_elementstatic_cast<double>(} appends elements to the end of the array.
 * When there are no open entries at the end of the internal storage array, * the array is expanded.  The size of the expanded array depends on the
 * {@code expansion_mode} and {@code expansion_factor} properties.
 * The {@code expansion_mode} determines whether the size of the array is
 * multiplied by the {@code expansion_factor}
 * ({@link Expansion_Mode#MULTIPLICATIVE}) or if the expansion is additive
 * ({@link Expansion_Mode#ADDITIVE} -- {@code expansion_factor} storage
 * locations added).
 * The default {@code expansion_mode} is {@code MULTIPLICATIVE} and the default
 * {@code expansion_factor} is 2.
 * <p>
 * The {@link #add_element_rollingstatic_cast<double>(} method adds a element to the end
 * of the internal storage array and adjusts the "usable window" of the
 * internal array forward by one position (effectively making what was the
 * second element the first, and so on).  Repeated activations of this method
 * (or activation of {@link #discard_front_elementsstatic_cast<int>(}) will effectively orphan
 * the storage locations at the beginning of the internal storage array.  To
 * reclaim this storage, each time one of these methods is activated, the size
 * of the internal storage array is compared to the number of addressable
 * elements (the {@code num_elements} property) and if the difference
 * is too large, the internal array is contracted to size
 * {@code num_elements + 1}.  The determination of when the internal
 * storage array is "too large" depends on the {@code expansion_mode} and
 * {@code contraction_factor} properties.  If  the {@code expansion_mode}
 * is {@code MULTIPLICATIVE}, contraction is triggered when the
 * ratio between storage array length and {@code num_elements} exceeds
 * {@code contraction_factor.}  If the {@code expansion_mode}
 * is {@code ADDITIVE}, the number of excess storage locations
 * is compared to {@code contraction_factor}.
 * <p>
 * To avoid cycles of expansions and contractions, the
 * {@code expansion_factor} must not exceed the {@code contraction_factor}.
 * Constructors and mutators for both of these properties enforce this
 * requirement, throwing a {@code } if it is
 * violated.
 * <p>
 * <b>Note:</b> this class is <b>NOT</b> thread-safe.
 */
class Resizable_Double_Array  
{
private:
    /** Default value for initial capacity. */
    static const int DEFAULT_INITIAL_CAPACITY = 16;
    /** Default value for array size modifier. */
    static const double DEFAULT_EXPANSION_FACTOR = 2.0;
    /** Default value for expansion mode. */
    static const Expansion_Mode DEFAULT_EXPANSION_MODE = Expansion_Mode.MULTIPLICATIVE;
    /**
     * Default value for the difference between {@link #contraction_criterion}
     * and {@link #expansion_factor}.
     */
    static const double DEFAULT_CONTRACTION_DELTA = 0.5;

    /**
     * The contraction criteria determines when the internal array will be
     * contracted to fit the number of elements contained in the element
     * array + 1.
     */
    const double contraction_criterion;

    /**
     * The expansion factor of the array.  When the array needs to be expanded, * the array size will be {@code internal_array.size() * expansion_factor}
     * if {@code expansion_mode} is set to MULTIPLICATIVE, or
     * {@code internal_array.size() + expansion_factor} if
     * {@code expansion_mode} is set to ADDITIVE.
     */
    const double expansion_factor;

    /**
     * Determines whether array expansion by {@code expansion_factor}
     * is additive or multiplicative.
     */
    const Expansion_Mode expansion_mode;

    /**
     * The internal storage array.
     */
    std::vector<double> internal_array;

    /**
     * The number of addressable elements in the array.  Note that this
     * has nothing to do with the length of the internal storage array.
     */
    int num_elements;

    /**
     * The position of the first addressable element in the internal storage
     * array.  The addressable elements in the array are
     * {@code internal_array[start_index],...,internal_array[start_index + num_elements - 1]}.
     */
    int start_index;

public:
    /** Specification of expansion algorithm. */
    enum Expansion_Mode 
    {
        /** Multiplicative expansion mode. */
        MULTIPLICATIVE, /** Additive expansion mode. */
        ADDITIVE
    }

    /**
     * Creates an instance with default properties.
     * <ul>
     *  <li>{@code initial_capacity = 16}</li>
     *  <li>{@code expansion_mode = MULTIPLICATIVE}</li>
     *  <li>{@code expansion_factor = 2.0}</li>
     *  <li>{@code contraction_criterion = 2.5}</li>
     * </ul>
     */
    Resizable_Double_Array() 
    {
        this(DEFAULT_INITIAL_CAPACITY);
    }

    /**
     * Creates an instance with the specified initial capacity.
     * <p>
     * Other properties take default values:
     * <ul>
     *  <li>{@code expansion_mode = MULTIPLICATIVE}</li>
     *  <li>{@code expansion_factor = 2.0}</li>
     *  <li>{@code contraction_criterion = 2.5}</li>
     * </ul>
     * @param initial_capacity Initial size of the internal storage array.
     * @ if {@code initial_capacity <= 0}.
     */
    Resizable_Double_Array(const int& initial_capacity)  
    {
        this(initial_capacity, DEFAULT_EXPANSION_FACTOR);
    }

    /**
     * Creates an instance from an existing {@code std::vector<double>} with the
     * initial capacity and num_elements corresponding to the size of
     * the supplied {@code std::vector<double>} array.
     * <p>
     * If the supplied array is NULL, a empty array with the default
     * initial capacity will be created.
     * The input array is copied, not referenced.
     * Other properties take default values:
     * <ul>
     *  <li>{@code expansion_mode = MULTIPLICATIVE}</li>
     *  <li>{@code expansion_factor = 2.0}</li>
     *  <li>{@code contraction_criterion = 2.5}</li>
     * </ul>
     *
     * @param initial_array initial array
     */
    Resizable_Double_Array(std::vector<double> initial_array) 
    {
        this(initial_array == NULL || initial_array.size() == 0 ?
             DEFAULT_INITIAL_CAPACITY : initial_array.size(), DEFAULT_EXPANSION_FACTOR, DEFAULT_CONTRACTION_DELTA + DEFAULT_EXPANSION_FACTOR, DEFAULT_EXPANSION_MODE, initial_array);
    }

    /**
     * Creates an instance with the specified initial capacity
     * and expansion factor.
     * <p>
     * The remaining properties take default values:
     * <ul>
     *  <li>{@code expansion_mode = MULTIPLICATIVE}</li>
     *  <li>{@code contraction_criterion = 0.5 + expansion_factor}</li>
     * </ul>
     * <p>
     * Throws  if the following conditions
     * are not met:
     * <ul>
     *  <li>{@code initial_capacity > 0}</li>
     *  <li>{@code expansion_factor > 1}</li>
     * </ul>
     *
     * @param initial_capacity Initial size of the internal storage array.
     * @param expansion_factor The array will be expanded based on this parameter.
     * @ if parameters are not valid.
     */
    Resizable_Double_Array(const int& initial_capacity, double expansion_factor)  
    {
        this(initial_capacity, expansion_factor, DEFAULT_CONTRACTION_DELTA + expansion_factor);
    }

    /**
     * Creates an instance with the specified initial capacity, * expansion factor, and contraction criteria.
     * <p>
     * The expansion mode will default to {@code MULTIPLICATIVE}.
     * <p>
     * Throws  if the following conditions
     * are not met:
     * <ul>
     *  <li>{@code initial_capacity > 0}</li>
     *  <li>{@code expansion_factor > 1}</li>
     *  <li>{@code contraction_criterion >= expansion_factor}</li>
     * </ul>
     *
     * @param initial_capacity Initial size of the internal storage array.
     * @param expansion_factor The array will be expanded based on this parameter.
     * @param contraction_criterion Contraction criterion.
     * @ if the parameters are not valid.
     */
    Resizable_Double_Array(const int& initial_capacity, double expansion_factor, double contraction_criterion)
    {
        this(initial_capacity, expansion_factor, contraction_criterion, DEFAULT_EXPANSION_MODE, NULL);
    }

    /**
     * Creates an instance with the specified properties.
     * <br/>
     * Throws  if the following conditions
     * are not met:
     * <ul>
     *  <li>{@code initial_capacity > 0}</li>
     *  <li>{@code expansion_factor > 1}</li>
     *  <li>{@code contraction_criterion >= expansion_factor}</li>
     * </ul>
     *
     * @param initial_capacity Initial size of the internal storage array.
     * @param expansion_factor The array will be expanded based on this parameter.
     * @param contraction_criterion Contraction criteria.
     * @param expansion_mode Expansion mode.
     * @param data Initial contents of the array.
     * @ if the parameters are not valid.
     * @ if expansion_mode is NULL
     */
    Resizable_Double_Array(const int& initial_capacity, double expansion_factor, double contraction_criterion, Expansion_Mode expansion_mode, double ... data)
    {
        if (initial_capacity <= 0) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::INITIAL_CAPACITY_NOT_POSITIVE, initial_capacity);
        }
        check_contract_expand(contraction_criterion, expansion_factor);
        //Math_Utils::check_not_null(expansion_mode);

        this.expansion_factor = expansion_factor;
        this.contraction_criterion = contraction_criterion;
        this.expansion_mode = expansion_mode;
        internal_array = std::vector<double>(initial_capacity];
        num_elements = 0;
        start_index = 0;

        if (data != NULL && data.size() > 0) 
        {
            add_elements(data);
        }
    }

    /**
     * Copy constructor.
     * <p>
     * Creates a Resizable_Double_Array that is a deep, fresh copy of the original.
     * Original may not be NULL; otherwise a {@link } is thrown.
     *
     * @param original array to copy
     * @exception  if original is NULL
     */
    Resizable_Double_Array(const Resizable_Double_Array original)
    {
        //Math_Utils::check_not_null(original);
        this.contraction_criterion = original.contraction_criterion;
        this.expansion_factor = original.expansion_factor;
        this.expansion_mode = original.expansion_mode;
        this.internal_array = std::vector<double>(original.internal_array.size()];
        System.arraycopy(original.internal_array, 0, this.internal_array, 0, this.internal_array.size());
        this.num_elements = original.num_elements;
        this.start_index = original.start_index;
    }

    /**
     * Adds an element to the end of this expandable array.
     *
     * @param value Value to be added to end of array.
     */
    void add_element(const double value) 
    {
        if (internal_array.size() <= start_index + num_elements) 
        {
            expand();
        }
        internal_array[start_index + num_elements++] = value;
    }

    /**
     * Adds several element to the end of this expandable array.
     *
     * @param values Values to be added to end of array.
     */
    void add_elements(const std::vector<double>& values) 
    {
        const std::vector<double> temp_array = std::vector<double>(num_elements + values.size() + 1];
        System.arraycopy(internal_array, start_index, temp_array, 0, num_elements);
        System.arraycopy(values, 0, temp_array, num_elements, values.size());
        internal_array = temp_array;
        start_index = 0;
        num_elements += values.size();
    }

    /**
     * Adds an element to the end of the array and removes the first
     * element in the array.  Returns the discarded first element.
     * <p>
     * The effect is similar to a push operation in a FIFO queue.
     * <p>
     * Example: If the array contains the elements 1, 2, 3, 4 (in that order)
     * and add_element_rolling(5) is invoked, the result is an array containing
     * the entries 2, 3, 4, 5 and the value returned is 1.
     *
     * @param value Value to be added to the array.
     * @return the value which has been discarded or "pushed" out of the array
     * by this rolling insert.
     */
    double add_element_rolling(double value) 
    {
        double discarded = internal_array[start_index];

        if ((start_index + (num_elements + 1)) > internal_array.size()) 
        {
            expand();
        }
        // Increment the start index
        start_index += 1;

        // Add the value
        internal_array[start_index + (num_elements - 1)] = value;

        // Check the contraction criterion.
        if (should_contract()) 
        {
            contract();
        }
        return discarded;
    }

    /**
     * Substitutes {@code value} for the most recently added value.
     * <p>
     * Returns the value that has been replaced. If the array is empty (i.e.
     * if {@link #num_elements} is zero), an Math_Illegal_State_Exception is thrown.
     *
     * @param value New value to substitute for the most recently added value
     * @return the value that has been replaced in the array.
     * @Math_Illegal_State_Exception if the array is empty
     */
    double substitute_most_recent_element(double value) Math_Illegal_State_Exception 
    {
        if (num_elements < 1) 
        {
            throw Math_Illegal_State_Exception(hipparchus::exception::Localized_Core_Formats_Type::CANNOT_SUBSTITUTE_ELEMENT_FROM_EMPTY_ARRAY);
        }

        const int subst_index = start_index + (num_elements - 1);
        const double discarded = internal_array[subst_index];

        internal_array[subst_index] = value;

        return discarded;
    }

protected:

    /**
     * Checks the expansion factor and the contraction criterion and raises
     * an exception if the contraction criterion is smaller than the
     * expansion criterion.
     *
     * @param contraction Criterion to be checked.
     * @param expansion Factor to be checked.
     * @ if {@code contraction < expansion}.
     * @ if {@code contraction <= 1}.
     * @ if {@code expansion <= 1 }.
     */
    void check_contract_expand(double contraction, double expansion)
    {
        if (contraction < expansion) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::CONTRACTION_CRITERIA_SMALLER_THAN_EXPANSION_FACTOR, contraction, expansion);
        }

        if (contraction <= 1) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::CONTRACTION_CRITERIA_SMALLER_THAN_ONE, contraction);
        }

        if (expansion <= 1) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::EXPANSION_FACTOR_SMALLER_THAN_ONE, expansion);
        }
    }

public:

    /**
     * Clear the array contents, resetting the number of elements to zero.
     */
    void clear() 
    {
        num_elements = 0;
        start_index = 0;
    }

    /**
     * Contracts the storage array to the (size of the element set) + 1 - to avoid
     * a zero length array. This function also resets the start_index to zero.
     */
    void contract() 
    {
        const std::vector<double> temp_array = std::vector<double>(num_elements + 1];

        // Copy and swap - copy only the element array from the src array.
        System.arraycopy(internal_array, start_index, temp_array, 0, num_elements);
        internal_array = temp_array;

        // Reset the start index to zero
        start_index = 0;
    }

    /**
     * Discards the {@code i} initial elements of the array.
     * <p>
     * For example, if the array contains the elements 1,2,3,4, invoking
     * {@code discard_front_elements(2)} will cause the first two elements
     * to be discarded, leaving 3,4 in the array.
     *
     * @param i  the number of elements to discard from the front of the array
     * @ if i is greater than num_elements.
     */
    void discard_front_elements(const int& i)  
    {
        discard_extreme_elements(i,true);
    }

    /**
     * Discards the {@code i} last elements of the array.
     * <p>
     * For example, if the array contains the elements 1,2,3,4, invoking
     * {@code discard_most_recent_elements(2)} will cause the last two elements
     * to be discarded, leaving 1,2 in the array.
     *
     * @param i  the number of elements to discard from the end of the array
     * @ if i is greater than num_elements.
     */
    void discard_most_recent_elements(const int& i)  
    {
        discard_extreme_elements(i,false);
    }

    /**
     * Discards the {@code i} first or last elements of the array, * depending on the value of {@code front}.
     * <p>
     * For example, if the array contains the elements 1,2,3,4, invoking
     * {@code discard_extreme_elements(2,false)} will cause the last two elements
     * to be discarded, leaving 1,2 in the array.
     * For example, if the array contains the elements 1,2,3,4, invoking
     * {@code discard_extreme_elements(2,true)} will cause the first two elements
     * to be discarded, leaving 3,4 in the array.
     *
     * @param i  the number of elements to discard from the front/end of the array
     * @param front true if elements are to be discarded from the front
     * of the array, false if elements are to be discarded from the end
     * of the array
     * @ if i is greater than num_elements.
     */
    private void discard_extreme_elements(const int& i, bool front)  
    {
        if (i > num_elements) 
        {
            throw (
                    hipparchus::exception::Localized_Core_Formats_Type::TOO_MANY_ELEMENTS_TO_DISCARD_FROM_ARRAY, i, num_elements);
       }
        else if (i < 0) 
       {
           throw (
                   hipparchus::exception::Localized_Core_Formats_Type::CANNOT_DISCARD_NEGATIVE_NUMBER_OF_ELEMENTS, i);
        }
        else 
        {
            // "Subtract" this number of discarded from num_elements
            num_elements -= i;
            if (front) 
            {
                start_index += i;
            }
        }
        if (should_contract()) 
        {
            contract();
        }
    }

    /**
     * Expands the internal storage array using the expansion factor.
     * <p>
     * If {@code expansion_mode} is set to MULTIPLICATIVE, * the array size will be {@code internal_array.size() * expansion_factor}.
     * If {@code expansion_mode} is set to ADDITIVE, the length
     * after expansion will be {@code internal_array.size() + expansion_factor}.
     */
    protected void expand() 
    {
        // notice the use of std::ceil(), this guarantees that we will always
        // have an array of at least current_size + 1.   Assume that the
        // current initial capacity is 1 and the expansion factor
        // is 1.000000000000000001.  The newly calculated size will be
        // rounded up to 2 after the multiplication is performed.
        const int& new_size;
        if (expansion_mode == Expansion_Mode.MULTIPLICATIVE) 
        {
            new_size = static_cast<int>( std::ceil(internal_array.size() * expansion_factor);
        }
        else 
        {
            new_size = static_cast<int>( (internal_array.size() + std::round(expansion_factor));
        }
        const std::vector<double> temp_array = std::vector<double>(new_size];

        // Copy and swap
        System.arraycopy(internal_array, 0, temp_array, 0, internal_array.size());
        internal_array = temp_array;
    }

    /**
     * Expands the internal storage array to the specified size.
     *
     * @param size Size of the internal storage array.
     */
    private void expand_to(const int& size) 
    {
        const std::vector<double> temp_array = std::vector<double>(size];
        // Copy and swap
        System.arraycopy(internal_array, 0, temp_array, 0, internal_array.size());
        internal_array = temp_array;
    }

    /**
     * The contraction criterion defines when the internal array will contract
     * to store only the number of elements in the element array.
     * <p>
     * If the {@code expansion_mode} is {@code MULTIPLICATIVE}, * contraction is triggered when the ratio between storage array length
     * and {@code num_elements} exceeds {@code contraction_factor}.
     * If the {@code expansion_mode} is {@code ADDITIVE}, the
     * number of excess storage locations is compared to {@code contraction_factor}.
     *
     * @return the contraction criterion used to reclaim memory.
     */
    public double get_contraction_criterion() 
    {
        return contraction_criterion;
    }

    /**
     * Returns the element at the specified index.
     *
     * @param index index to fetch a value from
     * @return value stored at the specified index
     * @Array_indexOutOfboundsException if {@code index} is less than
     * zero or is greater than {@code get_num_elements() - 1}.
     */
    public double get_element(const int& index) 
    {
        if (index >= num_elements) 
        {
            throw Array_indexOutOfboundsException(index);
        }
        if (index >= 0) 
        {
            return internal_array[start_index + index];
        }
        throw Array_indexOutOfboundsException(index);
    }

     /**
     * Returns a double array containing the elements of this Resizable_Array.
     * <p>
     * This method returns a copy, not a reference to the underlying array, * so that changes made to the returned array have no effect on this Resizable_Array.
     *
     * @return the double array.
     */
    public std::vector<double> get_elements() 
    {
        const std::vector<double> element_array = std::vector<double>(num_elements];
        System.arraycopy(internal_array, start_index, element_array, 0, num_elements);
        return element_array;
    }

    /**
     * The expansion factor controls the size of a array when an array
     * needs to be expanded.
     * <p>
     * The {@code expansion_mode} determines whether the size of the array
     * is multiplied by the {@code expansion_factor} (MULTIPLICATIVE) or if
     * the expansion is additive (ADDITIVE -- {@code expansion_factor}
     * storage locations added).  The default {@code expansion_mode} is
     * MULTIPLICATIVE and the default {@code expansion_factor} is 2.0.
     *
     * @return the expansion factor of this expandable double array
     */
    public double get_expansion_factor() 
    {
        return expansion_factor;
    }

    /**
     * The expansion mode determines whether the internal storage
     * array grows additively or multiplicatively when it is expanded.
     *
     * @return the expansion mode.
     */
    public Expansion_Mode get_expansion_mode() 
    {
        return expansion_mode;
    }

    /**
     * Gets the currently allocated size of the internal data structure used
     * for storing elements.
     * This is not to be confused with {@link #get_num_elements() the number of
     * elements actually stored}.
     *
     * @return the length of the internal array.
     */
    public int get_capacity() 
    {
        return internal_array.size();
    }

    /**
     * Returns the number of elements currently in the array.  Please note
     * that this is different from the length of the internal storage array.
     *
     * @return the number of elements.
     */
    public int get_num_elements() 
    {
        return num_elements;
    }

    /**
     * Provides <em>direct</em> access to the internal storage array.
     * Please note that this method returns a reference to this object's
     * storage array, not a copy.
     * <p>
     * To correctly address elements of the array, the "start index" is
     * required (available via the {@link #get_start_index() get_start_index}
     * method.
     * <p>
     * This method should only be used to avoid copying the internal array.
     * The returned value <em>must</em> be used for reading only; other
     * uses could lead to this object becoming inconsistent.
     * <p>
     * The {@link #get_elements} method has no such limitation since it
     * returns a copy of this array's addressable elements.
     *
     * @return the internal storage array used by this object.
     */
    protected std::vector<double> get_array_ref() 
    {
        return internal_array; // NOPMD - returning an internal array is intentional and documented here
    }

    /**
     * Returns the "start index" of the internal array.
     * This index is the position of the first addressable element in the
     * internal storage array.
     * <p>
     * The addressable elements in the array are at indices contained in
     * the interval [{@link #get_start_index()}, *               {@link #get_start_index()} + {@link #get_num_elements()} - 1].
     *
     * @return the start index.
     */
    protected int get_start_index() 
    {
        return start_index;
    }

    /**
     * Performs an operation on the addressable elements of the array.
     *
     * @param f Function to be applied on this array.
     * @return the result.
     */
    public double compute(Math_Arrays::Function f) 
    {
        return f.evaluate(internal_array, start_index, num_elements);
    }

    /**
     * Sets the element at the specified index.
     * <p>
     * If the specified index is greater than {@code get_num_elements() - 1}, * the {@code num_elements} property is increased to {@code index +1}
     * and additional storage is allocated (if necessary) for the element and
     * all (uninitialized) elements between the element and the previous end
     * of the array).
     *
     * @param index index to store a value in
     * @param value value to store at the specified index
     * @Array_indexOutOfboundsException if {@code index < 0}.
     */
    public void set_element(const int& index, double value) 
    {
        if (index < 0) 
        {
            throw Array_indexOutOfboundsException(index);
        }
        if (index + 1 > num_elements) 
        {
            num_elements = index + 1;
        }
        if ((start_index + index) >= internal_array.size()) 
        {
            expand_to(start_index + (index + 1));
        }
        internal_array[start_index + index] = value;
    }

    /**
     * This function allows you to control the number of elements contained
     * in this array, and can be used to "throw out" the last n values in an
     * array. This function will also expand the internal array as needed.
     *
     * @param i a number of elements
     * @ if {@code i} is negative.
     */
    public void set_num_elements(const int& i)  
    {
        // If index is negative thrown an error.
        if (i < 0) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::INDEX_NOT_POSITIVE, i);
        }

        // Test the num elements, check to see if the array needs to be
        // expanded to accommodate this number of elements.
        const int& new_size = start_index + i;
        if (new_size > internal_array.size()) 
        {
            expand_to(new_size);
        }

        // Set the number of elements to value.
        num_elements = i;
    }

    /**
     * Returns true if the internal storage array has too many unused
     * storage positions.
     *
     * @return true if array satisfies the contraction criteria
     */
    private bool should_contract() 
    {
        if (expansion_mode == Expansion_Mode.MULTIPLICATIVE) 
        {
            return (internal_array.size() / ((float) num_elements)) > contraction_criterion;
        }
else 
        {
            return (internal_array.size() - num_elements) > contraction_criterion;
        }
    }

    /**
     * Returns a copy of the Resizable_Double_Array.  Does not contract before
     * the copy, so the returned object is an exact copy of this.
     *
     * @return a Resizable_Double_Array with the same data and configuration
     * properties as this
     */
    public Resizable_Double_Array copy() 
    {
        return Resizable_Double_Array(this);
    }

    /**
     * Returns true iff object is a Resizable_Double_Array with the same properties
     * as this and an identical internal storage array.
     *
     * @param object object to be compared for equality with this
     * @return true iff object is a Resizable_Double_Array with the same data and
     * properties as this
     */
    //override
    public bool equals(Object object) 
    {
        if (object == this) 
        {
            return true;
        }
        if (!dynamic_cast<const Resizable_Double_Array*>(*object) != nullptr)
        {
            return false;
        }
        bool result{ true };
        const Resizable_Double_Array other = (Resizable_Double_Array) object;
        result = result && (other.contraction_criterion == contraction_criterion);
        result = result && (other.expansion_factor == expansion_factor);
        result = result && (other.expansion_mode == expansion_mode);
        result = result && (other.num_elements == num_elements);
        result = result && (other.start_index == start_index);
        if (!result) 
        {
            return false;
        }
        return Arrays.equals(internal_array, other.internal_array);
    }

    /**
     * Returns a hash code consistent with equals.
     *
     * @return the hash code representing this {@code Resizable_Double_Array}.
     */
    //override
    public int hash_code() 
    {
        auto hash_data = std::vector<int>(6);
        hash_data[0] = static_cast<double>(expansion_factor).hash_code();
        hash_data[1] = static_cast<double>(contraction_criterion).hash_code();
        hash_data[2] = expansion_mode.hash_code();
        hash_data[3] = Arrays.hash_code(internal_array);
        hash_data[4] = num_elements;
        hash_data[5] = start_index;
        return Arrays.hash_code(hash_data);
    }
};