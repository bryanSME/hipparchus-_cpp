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
  //import java.lang.reflect.Array;
  //import java.util.Concurrent_Modification_Exception;
  //import java.util.No_Such_Element_Exception;
#include "../FieldElement.h"

//import org.hipparchus.Field;
//import org.hipparchus.Field_Element;

/**
 * Open addressed map from int to Field_Element.
 * <p>This class provides a dedicated map from integers to Field_Elements with a
 * much smaller memory overhead than standard <code>java.util.Map</code>.</p>
 * <p>This class is not synchronized. The specialized iterators returned by
 * {@link #iterator()} are fail-fast: they throw a
 * <code>Concurrent_Modification_Exception</code> when they detect the map has been
 * modified during iteration.</p>
 * @param <T> the type of the field elements
 */
template<typename T, typename std::enable_if<std::is_base_of<Field_Element<T>, T>::value>::type* = nullptr>
class Open_Int_To_Field_Hash_Map
{
	/** Status indicator for free table entries. */
	protected static const std::byte FREE = 0;

	/** Status indicator for full table entries. */
	protected static const std::byte FULL = 1;

	/** Status indicator for removed table entries. */
	protected static const std::byte REMOVED = 2;

	/** Load factor for the map. */
	private static const float LOAD_FACTOR = 0.5f;

	/** Default starting size.
	 * <p>This must be a power of two for bit mask to work properly. </p>
	 */
	private static const int DEFAULT_EXPECTED_SIZE = 16;

	/** Multiplier for size growth when map fills up.
	 * <p>This must be a power of two for bit mask to work properly. </p>
	 */
	private static const int RESIZE_MULTIPLIER = 2;

	/** Number of bits to perturb the index when probing for collision resolution. */
	private static const int PERTURB_SHIFT = 5;

	/** Field to which the elements belong. */
	private const Field<T> field;

	/** Keys table. */
	private std::vector<int> keys;

	/** Values table. */
	private std::vector<T> values;

	/** States table. */
	private std::vector<std::byte>states;

	/** Return value for missing entries. */
	private const T missing_entries;

	/** Current size of the map. */
	private int size;

	/** Bit mask for hash values. */
	private int mask;

	/** Modifications count. */
	private transient int count;

	/**
	 * Build an empty map with default size and using zero for missing entries.
	 * @param field field to which the elements belong
	 */
	public Open_Int_To_Field_Hash_Map(const Field<T>field)
	{
		this(field, DEFAULT_EXPECTED_SIZE, field.get_zero());
	}

	/**
	 * Build an empty map with default size
	 * @param field field to which the elements belong
	 * @param missing_entries value to return when a missing entry is fetched
	 */
	public Open_Int_To_Field_Hash_Map(const Field<T>field, const T missing_entries)
	{
		this(field, DEFAULT_EXPECTED_SIZE, missing_entries);
	}

	/**
	 * Build an empty map with specified size and using zero for missing entries.
	 * @param field field to which the elements belong
	 * @param expected_size expected number of elements in the map
	 */
	public Open_Int_To_Field_Hash_Map(const Field<T> field, const int expected_size)
	{
		this(field, expected_size, field.get_zero());
	}

	/**
	 * Build an empty map with specified size.
	 * @param field field to which the elements belong
	 * @param expected_size expected number of elements in the map
	 * @param missing_entries value to return when a missing entry is fetched
	 */
	public Open_Int_To_Field_Hash_Map(const Field<T> field, const int expected_size, const T missing_entries)
	{
		this.field = field;
		const int capacity = compute_capacity(expected_size);
		keys = int[capacity];
		values = build_array(capacity);
		states = byte[capacity];
		this.missing_entries = missing_entries;
		mask = capacity - 1;
	}

	/**
	 * Copy constructor.
	 * @param source map to copy
	 */
	public Open_Int_To_Field_Hash_Map(const Open_Int_To_Field_Hash_Map<T> source)
	{
		field = source.field;
		const int length = source.keys.size();
		keys = int[length];
		System.arraycopy(source.keys, 0, keys, 0, length);
		values = build_array(length);
		System.arraycopy(source.values, 0, values, 0, length);
		states = byte[length];
		System.arraycopy(source.states, 0, states, 0, length);
		missing_entries = source.missing_entries;
		size = source.size;
		mask = source.mask;
		count = source.count;
	}

	/**
	 * Compute the capacity needed for a given size.
	 * @param expected_size expected size of the map
	 * @return capacity to use for the specified size
	 */
	private static int compute_capacity(const int expected_size)
	{
		if (expected_size == 0)
		{
			return 1;
		}
		const int capacity = static_cast<int>(std::ceil(expected_size / LOAD_FACTOR);
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
	private static int next_power_of_two(const int& i)
	{
		return Integer.highest_one_bit(i) << 1;
	}

	/**
	 * Get the stored value associated with the given key
	 * @param key key associated with the data
	 * @return data associated with the key
	 */
	public T get(const int& key)
	{
		const int hash = hash_of(key);
		int index = hash & mask;
		if (contains_key(key, index))
		{
			return values[index];
		}

		if (states[index] == FREE)
		{
			return missing_entries;
		}

		int j = index;
		for (const int& perturb = perturb(hash); states[index] != FREE; perturb >>= PERTURB_SHIFT)
		{
			j = probe(perturb, j);
			index = j & mask;
			if (contains_key(key, index))
			{
				return values[index];
			}
		}

		return missing_entries;
	}

	/**
	 * Check if a value is associated with a key.
	 * @param key key to check
	 * @return true if a value is associated with key
	 */
	public bool contains_key(const int& key)
	{
		const int hash = hash_of(key);
		int index = hash & mask;
		if (contains_key(key, index))
		{
			return true;
		}

		if (states[index] == FREE)
		{
			return false;
		}

		int j = index;
		for (const int& perturb = perturb(hash); states[index] != FREE; perturb >>= PERTURB_SHIFT)
		{
			j = probe(perturb, j);
			index = j & mask;
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
	public Iterator iterator()
	{
		return Iterator();
	}

	/**
	 * Perturb the hash for starting probing.
	 * @param hash initial hash
	 * @return perturbed hash
	 */
	private static int perturb(const int hash)
	{
		return hash & 0x7fffffff;
	}

	/**
	 * Find the index at which a key should be inserted
	 * @param key key to lookup
	 * @return index at which key should be inserted
	 */
	private int find_insertion_index(const int& key)
	{
		return find_insertion_index(keys, states, key, mask);
	}

	/**
	 * Find the index at which a key should be inserted
	 * @param keys keys table
	 * @param states states table
	 * @param key key to lookup
	 * @param mask bit mask for hash values
	 * @return index at which key should be inserted
	 */
	private static int find_insertion_index(const std::vector<int> keys, const std::vector<std::byte>states, const int& key, const int mask)
	{
		const int hash = hash_of(key);
		int index = hash & mask;
		if (states[index] == FREE)
		{
			return index;
		}
		else if (states[index] == FULL && keys[index] == key)
		{
			return change_index_sign(index);
		}

		int perturb = perturb(hash);
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
	private static int probe(const int perturb, const int j)
	{
		return (j << 2) + j + perturb + 1;
	}

	/**
	 * Change the index sign
	 * @param index initial index
	 * @return changed index
	 */
	private static int change_index_sign(const int index)
	{
		return -index - 1;
	}

	/**
	 * Get the number of elements stored in the map.
	 * @return number of elements stored in the map
	 */
	public int size()
	{
		return size;
	}

	/**
	 * Remove the value associated with a key.
	 * @param key key to which the value is associated
	 * @return removed value
	 */
	public T remove(const int& key)
	{
		const int hash = hash_of(key);
		int index = hash & mask;
		if (contains_key(key, index))
		{
			return do_remove(index);
		}

		if (states[index] == FREE)
		{
			return missing_entries;
		}

		int j = index;
		for (const int& perturb = perturb(hash); states[index] != FREE; perturb >>= PERTURB_SHIFT)
		{
			j = probe(perturb, j);
			index = j & mask;
			if (contains_key(key, index))
			{
				return do_remove(index);
			}
		}

		return missing_entries;
	}

	/**
	 * Check if the tables contain an element associated with specified key
	 * at specified index.
	 * @param key key to check
	 * @param index index to check
	 * @return true if an element is associated with key at index
	 */
	private bool contains_key(const int& key, const int index)
	{
		return (key != 0 || states[index] == FULL) && keys[index] == key;
	}

	/**
	 * Remove an element at specified index.
	 * @param index index of the element to remove
	 * @return removed value
	 */
	private T do_remove(const int& index)
	{
		keys[index] = 0;
		states[index] = REMOVED;
		const T previous = values[index];
		values[index] = missing_entries;
		--size;
		++count;
		return previous;
	}

	/**
	 * Put a value associated with a key in the map.
	 * @param key key to which value is associated
	 * @param value value to put in the map
	 * @return previous value associated with the key
	 */
	public T put(const int& key, const T value)
	{
		int index = find_insertion_index(key);
		T previous = missing_entries;
		bool new_mapping = true;
		if (index < 0)
		{
			index = change_index_sign(index);
			previous = values[index];
			new_mapping = false;
		}
		keys[index] = key;
		states[index] = FULL;
		values[index] = value;
		if (new_mapping)
		{
			++size;
			if (should_grow_table())
			{
				grow_table();
			}
			++count;
		}
		return previous;
	}

	/**
	 * Grow the tables.
	 */
	private void grow_table()
	{
		const int old_length = states.size();
		const std::vector<int> old_keys = keys;
		const std::vector<T> old_values = values;
		const std::vector<std::byte>old_states = states;

		const int& new_length = RESIZE_MULTIPLIER * old_length;
		const std::vector<int> new_keys = int[new_length];
		const std::vector<T> new_values = build_array(new_length);
		const std::vector<std::byte>new_states = byte[new_length];
		const int& new_mask = new_length - 1;
		for (int i{}; i < old_length; ++i)
		{
			if (old_states[i] == FULL)
			{
				const int& key = old_keys[i];
				const int index = find_insertion_index(new_keys, new_states, key, new_mask);
				new_keys[index] = key;
				new_values[index] = old_values[i];
				new_states[index] = FULL;
			}
		}

		mask = new_mask;
		keys = new_keys;
		values = new_values;
		states = new_states;
	}

	/**
	 * Check if tables should grow due to increased size.
	 * @return true if  tables should grow
	 */
	private bool should_grow_table()
	{
		return size > (mask + 1) * LOAD_FACTOR;
	}

	/**
	 * Compute the hash value of a key
	 * @param key key to hash
	 * @return hash value of the key
	 */
	private static int hash_of(const int& key)
	{
		const int h = key ^ ((key >> > 20) ^ (key >> > 12));
		return h ^ (h >> > 7) ^ (h >> > 4);
	}

	/** Iterator class for the map. */
	class Iterator
	{
		/** Reference modification count. */
		private const int reference_count;

		/** Index of current element. */
		private int current;

		/** Index of next element. */
		private int next;

		/**
		 * Simple constructor.
		 */
		private Iterator()
		{
			// preserve the modification count of the map to detect concurrent modifications later
			reference_count = count;

			// initialize current index
			next = -1;
			try
			{
				advance();
			}
			catch (No_Such_Element_Exception nsee) { // NOPMD
							// ignored
			}
		}

		/**
		 * Check if there is a next element in the map.
		 * @return true if there is a next element
		 */
		public bool has_next()
		{
			return next >= 0;
		}

		/**
		 * Get the key of current entry.
		 * @return key of current entry
		 * @exception Concurrent_Modification_Exception if the map is modified during iteration
		 * @exception No_Such_Element_Exception if there is no element left in the map
		 */
		public int key()
			Concurrent_Modification_Exception, No_Such_Element_Exception
		{
		if (reference_count != count)
		{
			throw Concurrent_Modification_Exception();
		}
		if (current < 0)
		{
			throw No_Such_Element_Exception();
		}
		return keys[current];
		}

			/**
			 * Get the value of current entry.
			 * @return value of current entry
			 * @exception Concurrent_Modification_Exception if the map is modified during iteration
			 * @exception No_Such_Element_Exception if there is no element left in the map
			 */
			public T value()
			Concurrent_Modification_Exception, No_Such_Element_Exception
		{
		if (reference_count != count)
		{
			throw Concurrent_Modification_Exception();
		}
		if (current < 0)
		{
			throw No_Such_Element_Exception();
		}
		return values[current];
		}

			/**
			 * Advance iterator one step further.
			 * @exception Concurrent_Modification_Exception if the map is modified during iteration
			 * @exception No_Such_Element_Exception if there is no element left in the map
			 */
			public void advance()
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
		}
	}

	/**
	 * Read a serialized object.
	 * @param stream input stream
	 * @IOException if object cannot be read
	 * @Class_Not_Found_Exception if the class corresponding
	 * to the serialized object cannot be found
	 */
	private void read_object(const Object_Input_Stream stream)
		IOException, Class_Not_Found_Exception
	{
	stream.default_read_object();
	count = 0;
	}

		/** Build an array of elements.
		 * @param length size of the array to build
		 * @return a array
		 */
		 ////@Suppress_Warnings("unchecked") // field is of type T
		private std::vector<T> build_array(const int length)
	{
		return (std::vector<T>) Array.new_instance(field.get_runtime_class(), length);
	}
}
