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

//import java.io.IOException;
//import java.io.Object_Input_Stream;
//import java.io.Serializable;
//import java.util.Concurrent_Modification_Exception;
//import java.util.No_Such_Element_Exception;
#include <vector>
#include <cmath>

/**
 * Open addressed map from int to double.
 * <p>This class provides a dedicated map from integers to doubles with a
 * much smaller memory overhead than standard <code>java.util.Map</code>.</p>
 * <p>This class is not synchronized. The specialized iterators returned by
 * {@link #iterator()} are fail-fast: they throw a
 * <code>Concurrent_Modification_Exception</code> when they detect the map has been
 * modified during iteration.</p>
 */
class Open_Int_To_Double_Hash_Map  
{
protected:
    /** Status indicator for free table entries. */
    static constexpr std::byte FREE{};

    /** Status indicator for full table entries. */
    static constexpr std::byte FULL{ 1 };

    /** Status indicator for removed table entries. */
    static constexpr std::byte REMOVED{ 2 };;

private:

    /** Load factor for the map. */
    static constexpr float LOAD_FACTOR{ 0.5 };

    /** Default starting size.
     * <p>This must be a power of two for bit mask to work properly. </p>
     */
    static constexpr int DEFAULT_EXPECTED_SIZE{ 16 };

    /** Multiplier for size growth when map fills up.
     * <p>This must be a power of two for bit mask to work properly. </p>
     */
    static constexpr int RESIZE_MULTIPLIER{ 2 };

    /** Number of bits to perturb the index when probing for collision resolution. */
    static constexpr int PERTURB_SHIFT{ 5 };

    /** Keys table. */
    std::vector<int> my_keys;

    /** Values table. */
    std::vector<double> my_values;

    /** States table. */
    std::vector<std::byte> my_states;

    /** Return value for missing entries. */
    const double my_missing_entries;

    /** Current size of the map. */
    int my_size;

    /** Bit mask for hash values. */
    int my_mask;

    /** Modifications count. */
    int my_count;

    /**
     * Check if the tables contain an element associated with specified key
     * at specified index.
     * @param key key to check
     * @param index index to check
     * @return true if an element is associated with key at index
     */
    bool contains_key(const int& key, const int& index) const
    {
        return (key != 0 || my_states[index] == FULL) && my_keys[index] == key;
    }

    /**
     * Remove an element at specified index.
     * @param index index of the element to remove
     * @return removed value
     */
    double do_remove(const int& index)
    {
        my_keys[index] = 0;
        my_states[index] = REMOVED;
        const auto previous = my_values[index];
        my_values[index] = my_missing_entries;
        --my_size;
        ++my_count;
        return previous;
    }

    /**
     * Grow the tables.
     */
    void grow_table()
    {
        const auto old_length = my_states.size();
        const auto old_keys = my_keys;
        const auto old_values = my_values;
        const auto old_states = my_states;

        const int new_length = RESIZE_MULTIPLIER * old_length;
        auto new_keys = std::vector<int>(new_length);
        auto new_values = std::vector<double>(new_length);
        auto new_states = std::vector<std::byte>(new_length);
        const int new_mask = new_length - 1;
        for (int i{}; i < old_length; ++i)
        {
            if (old_states[i] == FULL)
            {
                const int key = old_keys[i];
                const int index = find_insertion_index(new_keys, new_states, key, new_mask);
                new_keys[index] = key;
                new_values[index] = old_values[i];
                new_states[index] = FULL;
            }
        }

        my_mask = new_mask;
        my_keys = new_keys;
        my_values = new_values;
        my_states = new_states;
    }

    /**
     * Check if tables should grow due to increased size.
     * @return true if  tables should grow
     */
    bool should_grow_table() const
    {
        return my_size > (my_mask + 1) * LOAD_FACTOR;
    }

    /**
     * Compute the hash value of a key
     * @param key key to hash
     * @return hash value of the key
     */
    static int hash_of(const int& key)
    {
        const int h = key ^ ((key >> > 20) ^ (key >> > 12));
        return h ^ (h >> > 7) ^ (h >> > 4);
    }

    /**
     * Compute the capacity needed for a given size.
     * @param expected_size expected size of the map
     * @return capacity to use for the specified size
     */
    static int compute_capacity(const int& expected_size)
    {
        if (expected_size == 0)
        {
            return 1;
        }
        const auto capacity = static_cast<int>(std::ceil(expected_size / LOAD_FACTOR));
        const int power_of_two = Integer.highest_one_bit(capacity);
        if (power_of_two == capacity)
        {
            return capacity;
        }
        return next_power_of_two(capacity);
    }

    /**
     * Find the smallest power of two greater than the input value
     * @param i input value
     * @return smallest power of two greater than the input value
     */
    static int next_power_of_two(const int& i)
    {
        throw std::exception("NOT IMPLEMENTED");
        //return Integer.highest_one_bit(i) << 1;
    }

    /**
     * Perturb the hash for starting probing.
     * @param hash initial hash
     * @return perturbed hash
     */
    static int calc_perturb(const int& hash)
    {
        return hash & 0x7fffffff;
    }

    /**
     * Find the index at which a key should be inserted
     * @param key key to lookup
     * @return index at which key should be inserted
     */
    int find_insertion_index(const int& key) const
    {
        return find_insertion_index(my_keys, my_states, key, my_mask);
    }

    /**
     * Find the index at which a key should be inserted
     * @param keys keys table
     * @param states states table
     * @param key key to lookup
     * @param mask bit mask for hash values
     * @return index at which key should be inserted
     */
    static int find_insertion_index(const std::vector<int>& keys, const std::vector<std::byte>& states, const int& key, const int& mask)
    {
        const int hash = hash_of(key);
        int index = hash & mask;
        if (states[index] == FREE)
        {
            return index;
        }
        if (states[index] == FULL && keys[index] == key)
        {
            return change_index_sign(index);
        }

        auto perturb = calc_perturb(hash);
        int j = index;
        if (states[index] == FULL)
        {
            while (true)
            {
                j = probe(perturb, j);
                index = j & mask;
                perturb >>= PERTURB_SHIFT;

                if (states[index] != FULL || keys[index] == key)
                {
                    break;
                }
            }
        }

        if (states[index] == FREE)
        {
            return index;
        }
        else if (states[index] == FULL)
        {
            // due to the loop exit condition, // if (states[index] == FULL) then keys[index] == key
            return change_index_sign(index);
        }

        const int first_removed = index;
        while (true)
        {
            j = probe(perturb, j);
            index = j & mask;

            if (states[index] == FREE)
            {
                return first_removed;
            }
            else if (states[index] == FULL && keys[index] == key)
            {
                return change_index_sign(index);
            }

            perturb >>= PERTURB_SHIFT;

        }

    }

    /**
     * Compute next probe for collision resolution
     * @param perturb perturbed hash
     * @param j previous probe
     * @return next probe
     */
    static int probe(const int& perturb, const int& j)
    {
        return (j << 2) + j + perturb + 1;
    }

    /**
     * Change the index sign
     * @param index initial index
     * @return changed index
     */
    static int change_index_sign(const int& index)
    {
        return -index - 1;
    }

    /**
     * Read a serialized object.
     * @param stream input stream
     * @IOException if object cannot be read
     * @Class_Not_Found_Exception if the class corresponding
     * to the serialized object cannot be found
     */
    void read_object(const Object_Input_Stream stream)
    {
        stream.default_read_object();
        my_count = 0;
    }

public:
    /**
     * Build an empty map with default size and using NaN for missing entries.
     */
    Open_Int_To_Double_Hash_Map() 
    {
        Open_Int_To_Double_Hash_Map(DEFAULT_EXPECTED_SIZE, NAN);
    }

    /**
     * Build an empty map with default size
     * @param missing_entries value to return when a missing entry is fetched
     */
    Open_Int_To_Double_Hash_Map(const double& missing_entries) 
    {
        Open_Int_To_Double_Hash_Map(DEFAULT_EXPECTED_SIZE, missing_entries);
    }

    /**
     * Build an empty map with specified size and using NaN for missing entries.
     * @param expected_size expected number of elements in the map
     */
    Open_Int_To_Double_Hash_Map(const int& expected_size) 
    {
        Open_Int_To_Double_Hash_Map(expected_size, NAN);
    }

    /**
     * Build an empty map with specified size.
     * @param expected_size expected number of elements in the map
     * @param missing_entries value to return when a missing entry is fetched
     */
    Open_Int_To_Double_Hash_Map(const int& expected_size, const double& missing_entries) 
    {
        const auto capacity = compute_capacity(expected_size);
        my_keys   = std::vector<int>(capacity);
        my_values = std::vector<double>(capacity);
        my_states = std::vector<std::byte>(capacity);
        my_missing_entries = missing_entries;
        my_mask   = capacity - 1;
    }

    /**
     * Copy constructor.
     * @param source map to copy
     */
    Open_Int_To_Double_Hash_Map(const Open_Int_To_Double_Hash_Map& source) 
    {
        const auto length = source.keys.size();
        my_keys = std::vector<int>(length);
        System.arraycopy(source.keys, 0, my_keys, 0, length);
        my_values = std::vector<double>(length);
        System.arraycopy(source.values, 0, my_values, 0, length);
        states = byte[length];
        System.arraycopy(source.states, 0, my_states, 0, length);
        my_missing_entries = source.missing_entries;
        my_size  = source.size;
        my_mask  = source.mask;
        my_count = source.count;
    }

    /**
     * Get the stored value associated with the given key
     * @param key key associated with the data
     * @return data associated with the key
     */
    double get(const int& key) 
    {

        const int hash  = hash_of(key);
        int index = hash & my_mask;
        if (contains_key(key, index)) 
        {
            return my_values[index];
        }

        if (my_states[index] == FREE)
        {
            return missing_entries;
        }

        int j = index;
        for (const int& perturb = calc_perturb(hash); my_states[index] != FREE; perturb >>= PERTURB_SHIFT)
        {
            j = probe(perturb, j);
            index = j & my_mask;
            if (contains_key(key, index)) 
            {
                return my_values[index];
            }
        }

        return my_missing_entries;

    }

    /**
     * Check if a value is associated with a key.
     * @param key key to check
     * @return true if a value is associated with key
     */
    bool contains_key(const int& key) 
    {

        const int hash  = hash_of(key);
        int index = hash & my_mask;
        if (contains_key(key, index)) 
        {
            return true;
        }

        if (my_states[index] == FREE)
        {
            return false;
        }

        int j = index;
        for (const int& perturb = calc_perturb(hash); my_states[index] != FREE; perturb >>= PERTURB_SHIFT)
        {
            j = probe(perturb, j);
            index = j & my_mask;
            if (contains_key(key, index)) 
            {
                return true;
            }
        }

        return false;

    }

    /**
     * Get an iterator over map elements.
     * <p>The specialized iterators returned are fail-fast: they throw a
     * <code>Concurrent_Modification_Exception</code> when they detect the map
     * has been modified during iteration.</p>
     * @return iterator over the map elements
     */
    Iterator iterator() 
    {
        return Iterator();
    }

    /**
     * Get the number of elements stored in the map.
     * @return number of elements stored in the map
     */
    int size() 
    {
        return my_size;
    }


    /**
     * Remove the value associated with a key.
     * @param key key to which the value is associated
     * @return removed value
     */
    double remove(const int& key) 
    {

        const int hash  = hash_of(key);
        int index = hash & my_mask;
        if (contains_key(key, index)) 
        {
            return do_remove(index);
        }

        if (my_states[index] == FREE)
        {
            return missing_entries;
        }

        int j = index;
        for (const int& perturb = calc_perturb(hash); my_states[index] != FREE; perturb >>= PERTURB_SHIFT)
        {
            j = probe(perturb, j);
            index = j & my_mask;
            if (contains_key(key, index)) 
            {
                return do_remove(index);
            }
        }

        return my_missing_entries;

    }



    /**
     * Put a value associated with a key in the map.
     * @param key key to which value is associated
     * @param value value to put in the map
     * @return previous value associated with the key
     */
    double put(const int& key, const double value) 
    {
        int index = find_insertion_index(key);
        double previous = my_missing_entries;
        bool new_mapping = true;
        if (index < 0) 
        {
            index = change_index_sign(index);
            previous = my_values[index];
            new_mapping = false;
        }
        my_keys[index]   = key;
        my_states[index] = FULL;
        my_values[index] = value;
        if (new_mapping) 
        {
            ++my_size;
            if (should_grow_table()) 
            {
                grow_table();
            }
            ++my_count;
        }
        return previous;

    }




    /** Iterator class for the map. */
    class Iterator 
    {
    private:
        /** Reference modification count. */
        int my_reference_count;

        /** Index of current element. */
        int my_current;

        /** Index of next element. */
        int my_next;

        /**
         * Simple constructor.
         */
        Iterator() 
        {

            // preserve the modification count of the map to detect concurrent modifications later
            my_reference_count = count;

            // initialize current index
            my_next = -1;
            try 
            {
                advance();
            }
            catch (No_Such_Element_Exception nsee) { // NOPMD
                // ignored
            }

        }
    public:
        /**
         * Check if there is a next element in the map.
         * @return true if there is a next element
         */
        bool has_next()
        {
            return my_next >= 0;
        };

        /**
         * Get the key of current entry.
         * @return key of current entry
         * @exception Concurrent_Modification_Exception if the map is modified during iteration
         * @exception No_Such_Element_Exception if there is no element left in the map
         */
        int key()
        {
            if (my_reference_count != count)
            {
                throw Concurrent_Modification_Exception();
            }
            if (my_current < 0)
            {
                throw No_Such_Element_Exception();
            }
            return my_keys[my_current];
        };

        /**
            * Get the value of current entry.
            * @return value of current entry
            * @exception Concurrent_Modification_Exception if the map is modified during iteration
            * @exception No_Such_Element_Exception if there is no element left in the map
            */
        double value()
        {
            if (my_reference_count != count)
            {
                throw Concurrent_Modification_Exception();
            }
            if (current < 0)
            {
                throw No_Such_Element_Exception();
            }
            return my_values[my_current];
        };

        /**
         * Advance iterator one step further.
         * @exception Concurrent_Modification_Exception if the map is modified during iteration
         * @exception No_Such_Element_Exception if there is no element left in the map
         */
        void advance()
            Concurrent_Modification_Exception, No_Such_Element_Exception 
            {

            if (reference_count != count) 
            {
                throw Concurrent_Modification_Exception();
            }

            // advance on step
            current = next;

            // prepare next step
            try 
            {
                while (states[++next] != FULL) { // NOPMD
                    // nothing to do
                }
            }
            catch (Array_indexOutOfboundsException e) 
            {
                next = -2;
                if (current < 0) 
                {
                    throw No_Such_Element_Exception(); // NOPMD
                }
            }

        };

    };
};

