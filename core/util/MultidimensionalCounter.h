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

  //import java.util.No_Such_Element_Exception;

  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
#include <vector>
#include <string>
#include <sstream>

/**
 * Converter between unidimensional storage structure and multidimensional
 * conceptual structure.
 * This utility will convert from indices in a multidimensional structure
 * to the corresponding index in a one-dimensional array. For example, * assuming that the ranges (in 3 dimensions) of indices are 2, 4 and 3, * the following correspondences, between 3-tuples indices and unidimensional
 * indices, will hold:
 * <ul>
 *  <li>(0, 0, 0) corresponds to 0</li>
 *  <li>(0, 0, 1) corresponds to 1</li>
 *  <li>(0, 0, 2) corresponds to 2</li>
 *  <li>(0, 1, 0) corresponds to 3</li>
 *  <li>...</li>
 *  <li>(1, 0, 0) corresponds to 12</li>
 *  <li>...</li>
 *  <li>(1, 3, 2) corresponds to 23</li>
 * </ul>
 */
class Multidimensional_Counter //: public Iterable<int>
{
private:
	/**
	 * Number of dimensions.
	 */
	const int my_dimension;
	/**
	 * Offset for each dimension.
	 */
	const std::vector<int> uni_counter_offset;
	/**
	 * Counter sizes.
	 */
	const std::vector<int> size;
	/**
	 * Total number of (one-dimensional) slots.
	 */
	const int total_size;
	/**
	 * Index of last dimension.
	 */
	const int last;

public:
	/**
	 * Perform iteration over the multidimensional counter.
	 */
	class Iterator //: java.util.auto
	{
	private:
		/**
		 * Multidimensional counter.
		 */
		const auto my_counter = std::vector<int>(dimension);
		/**
		 * Unidimensional counter.
		 */
		int my_count{ -1 };
		/**
		 * Maximum value for {@link #count}.
		 */
		const int my_max_count{ total_size - 1 };
	public:
		/**
		 * Create an iterator
		 * @see #iterator()
		 */
		Iterator()
		{
			counter[last] = -1;
		}

		/**
		 * {@inherit_doc}
		 */
		 //override
		bool has_next()
		{
			return count < max_count;
		}

		/**
		 * @return the unidimensional count after the counter has been
		 * incremented by {@code 1}.
		 * @No_Such_Element_Exception if {@link #has_next()} would have
		 * returned {@code false}.
		 */
		 //override
		Integer next()
		{
			if (!has_next())
			{
				throw std::exception("not implemented");
				//throw No_Such_Element_Exception();
			}

			for (int i{ last }; i >= 0; i--)
			{
				if (counter[i] == size[i] - 1)
				{
					counter[i] = 0;
				}
				else
				{
					++counter[i];
					break;
				}
			}

			return ++count;
		}

		/**
		 * Get the current unidimensional counter slot.
		 *
		 * @return the index within the unidimensionl counter.
		 */
		int get_count() const
		{
			return my_count;
		}
		/**
		 * Get the current multidimensional counter slots.
		 *
		 * @return the indices within the multidimensional counter.
		 */
		std::vector<int> get_counts() const
		{
			return my_counter;
		}

		/**
		 * Get the current count in the selected dimension.
		 *
		 * @param dim Dimension index.
		 * @return the count at the corresponding index for the current state
		 * of the iterator.
		 * @Index_Out_Of_Bounds_Exception if {@code index} is not in the
		 * correct interval (as defined by the length of the argument in the
		 * {@link Multidimensional_Counter#Multidimensional_Counter(std::vector<int>)
		 * constructor of the enclosing class}).
		 */
		int get_count(const int& dim) const
		{
			return my_counter[dim];
		}

		/**
		 * @Unsupported_Operation_Exception
		 */
		 //override
		void remove()
		{
			throw std::exception("not implmented");
			//throw Unsupported_Operation_Exception();
		}
	}

	/**
	 * Create a counter.
	 *
	 * @param size Counter sizes (number of slots in each dimension).
	 * @ if one of the sizes is
	 * negative or zero.
	 */
	Multidimensional_Counter(int... size)
	{
		my_dimension = size.size();
		my_size = size.clone();

		uni_counter_offset = std::vector<int>(dimension);

		last = dimension - 1;
		int tS = size[last];
		for (int i{}; i < last; i++)
		{
			int count{ 1 };
			for (int j{ i + 1 }; j < dimension; j++)
			{
				count *= size[j];
			}
			uni_counter_offset[i] = count;
			tS *= size[i];
		}
		uni_counter_offset[last] = 0;

		if (tS <= 0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL_BOUND_EXCLUDED, tS, 0);
		}

		total_size = tS;
	}

	/**
	 * Create an iterator over this counter.
	 *
	 * @return the iterator.
	 */
	 //override
	Iterator iterator()
	{
		return Iterator();
	}

	/**
	 * Get the number of dimensions of the multidimensional counter.
	 *
	 * @return the number of dimensions.
	 */
	int get_dimension() const
	{
		return my_dimension;
	}

	/**
	 * Convert to multidimensional counter.
	 *
	 * @param index Index in unidimensional counter.
	 * @return the multidimensional counts.
	 * @ if {@code index} is not between
	 * {@code 0} and the value returned by {@link #get_size()} (excluded).
	 */
	std::vector<int> get_counts(const int& index)
	{
		Math_Utils::check_range_inclusive(index, 0, total_size - 1);

		auto indices = std::vector<int>(my_dimension);

		int count{};
		for (int i{}; i < last; i++)
		{
			int idx{};
			const auto offset = uni_counter_offset[i];
			while (count <= index)
			{
				count += offset;
				++idx;
			}
			--idx;
			count -= offset;
			indices[i] = idx;
		}

		indices[last] = index - count;

		return indices;
	}

	/**
	 * Convert to unidimensional counter.
	 *
	 * @param c Indices in multidimensional counter.
	 * @return the index within the unidimensionl counter.
	 * @ if the size of {@code c}
	 * does not match the size of the array given in the constructor.
	 * @ if a value of {@code c} is not in
	 * the range of the corresponding dimension, as defined in the
	 * {@link Multidimensional_Counter#Multidimensional_Counter(int...) constructor}.
	 */
	int get_count(const int& ... c)
	{
		Math_Utils::check_dimension(c.size(), dimension);
		int count{};
		for (int i{}; i < dimension; i++)
		{
			const int index = c[i];
			Math_Utils::check_range_inclusive(index, 0, size[i] - 1);
			count += uni_counter_offset[i] * c[i];
		}
		return count + c[last];
	}

	/**
	 * Get the total number of elements.
	 *
	 * @return the total size of the unidimensional counter.
	 */
	int get_size() const
	{
		return my_total_size;
	}
	/**
	 * Get the number of multidimensional counter slots in each dimension.
	 *
	 * @return the sizes of the multidimensional counter in each dimension.
	 */
	std::vector<int> get_sizes() const
	{
		return my_size;
	}

	/**
	 * {@inherit_doc}
	 */
	 //override
	std::string to_string() const
	{
		const std::stringstream sb = std::stringstream();
		for (int i{}; i < dimension; i++)
		{
			sb.append('[').append(get_count(i)).append(']');
		}
		return sb.to_string();
	}
};