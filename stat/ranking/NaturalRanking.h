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

  //package org.hipparchus.stat.ranking;

  //import java.util.Array_list;
  //import java.util.Arrays;
  //import java.util.Iterator;
  //import java.util.List;

  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.exception.Math_Runtime_Exception;
  //import org.hipparchus.random.Random_Data_Generator;
  //import org.hipparchus.random.Random_Generator;
  //import org.hipparchus.util.FastMath;
#include "RankingAlgorithm.h"
#include "NaNStrategy.h"
#include "TiesStrategy.h"
#include <vector>
#include <cmath>
#include "../../core/random/RandomDataGenerator.h"

/**
 * <p> Ranking based on the natural ordering on doubles.</p>
 * <p>NaNs are treated according to the configured {@link NaN_Strategy} and ties
 * are handled using the selected {@link Ties_Strategy}.
 * Configuration settings are supplied in optional constructor arguments.
 * Defaults are {@link NaN_Strategy#FAILED} and {@link Ties_Strategy#AVERAGE}, * respectively. When using {@link Ties_Strategy#RANDOM}, a
 * {@link Random_Generator} may be supplied as a constructor argument.</p>
 * <p>Examples:
 * <table border="1" cellpadding="3">
 * <tr><th colspan="3">
 * Input data: (20, 17, 30, 42.3, 17, 50,NAN, -INFINITY, 17)
 * </th></tr>
 * <tr><th>NaN_Strategy</th><th>Ties_Strategy</th>
 * <th><code>rank(data)</code></th>
 * <tr>
 * <td>default (NaNs maximal)</td>
 * <td>default (ties averaged)</td>
 * <td>(5, 3, 6, 7, 3, 8, 9, 1, 3)</td></tr>
 * <tr>
 * <td>default (NaNs maximal)</td>
 * <td>MINIMUM</td>
 * <td>(5, 2, 6, 7, 2, 8, 9, 1, 2)</td></tr>
 * <tr>
 * <td>MINIMAL</td>
 * <td>default (ties averaged)</td>
 * <td>(6, 4, 7, 8, 4, 9, 1.5, 1.5, 4)</td></tr>
 * <tr>
 * <td>REMOVED</td>
 * <td>SEQUENTIAL</td>
 * <td>(5, 2, 6, 7, 3, 8, 1, 4)</td></tr>
 * <tr>
 * <td>MINIMAL</td>
 * <td>MAXIMUM</td>
 * <td>(6, 5, 7, 8, 5, 9, 2, 2, 5)</td></tr></table></p>
 *
 */
class Natural_Ranking : public Ranking_Algorithm
{
private:
	/** NaN strategy - defaults to NaNs maximal */
	const NaN_Strategy my_nan_strategy;

	/** Ties strategy - defaults to ties averaged */
	const Ties_Strategy my_ties_strategy;

	/** Source of random data - used only when ties strategy is RANDOM */
	const Random_Data_Generator my_random_data;

public:

	/** default NaN strategy */
	static const NaN_Strategy DEFAULT_NAN_STRATEGY = NaN_Strategy::FAILED;

	/** default ties strategy */
	static const Ties_Strategy DEFAULT_TIES_STRATEGY = Ties_Strategy::AVERAGE;

	/**
	 * Create a Natural_Ranking with default strategies for handling ties and NaNs.
	 */
	Natural_Ranking()
		:
		my_ties_strategy{ DEFAULT_TIES_STRATEGY },
		my_nan_strategy{ DEFAULT_NAN_STRATEGY },
		my_random_data{ NULL }
	{
		super();
	}

	/**
	 * Create a Natural_Ranking with the given Ties_Strategy.
	 *
	 * @param ties_strategy the Ties_Strategy to use
	 */
	Natural_Ranking(const Ties_Strategy& ties_strategy)
		:
		my_ties_strategy{ ties_strategy },
		my_nan_strategy{ DEFAULT_NAN_STRATEGY },
		my_random_data{ Random_Data_Generator() }
	{
		super();
	}

	/**
	 * Create a Natural_Ranking with the given NaN_Strategy.
	 *
	 * @param nan_strategy the NaN_Strategy to use
	 */
	Natural_Ranking(NaN_Strategy nan_strategy)
		:
		my_nan_strategy{ nan_strategy },
		my_ties_strategy{ DEFAULT_TIES_STRATEGY },
		my_random_data{ NULL }
	{
		super();
	}

	/**
	 * Create a Natural_Ranking with the given NaN_Strategy and Ties_Strategy.
	 *
	 * @param nan_strategy NaN_Strategy to use
	 * @param ties_strategy Ties_Strategy to use
	 */
	Natural_Ranking(const NaN_Strategy& nan_strategy, const Ties_Strategy& ties_strategy)
		:
		my_nan_strategy{ nan_strategy },
		my_ties_strategy{ ties_strategy },
		my_random_data{ Random_Data_Generator() }
	{
		super();
	}

	/**
	 * Create a Natural_Ranking with Ties_Strategy.RANDOM and the given
	 * Random_Generator as the source of random data.
	 *
	 * @param random_generator source of random data
	 */
	Natural_Ranking(Random_Generator random_generator)
		:
		my_nan_strategy{ DEFAULT_NAN_STRATEGY },
		my_ties_strategy{ Ties_Strategy::RANDOM },
		my_random_data{ Random_Data_Generator::of(random_generator) }
	{
		super();
	}

	/**
	 * Create a Natural_Ranking with the given NaN_Strategy, Ties_Strategy.RANDOM
	 * and the given source of random data.
	 *
	 * @param nan_strategy NaN_Strategy to use
	 * @param random_generator source of random data
	 */
	Natural_Ranking(NaN_Strategy nan_strategy, Random_Generator random_generator)
		:
		my_nan_strategy{ nan_strategy },
		my_ties_strategy{ Ties_Strategy::RANDOM },
		my_random_data{ Random_Data_Generator::of(random_generator) }
	{
		super();
	}

	/**
	 * Return the NaN_Strategy
	 *
	 * @return returns the NaN_Strategy
	 */
	NaN_Strategy get_nan_strategy() const
	{
		return my_nan_strategy;
	}

	/**
	 * Return the Ties_Strategy
	 *
	 * @return the Ties_Strategy
	 */
	Ties_Strategy get_ties_strategy() const
	{
		return my_ties_strategy;
	}

	/**
	 * Rank <code>data</code> using the natural ordering on Doubles, with
	 * NaN values handled according to <code>nan_strategy</code> and ties
	 * resolved using <code>ties_strategy.</code>
	 *
	 * @param data array to be ranked
	 * @return array of ranks
	 * @ if the selected {@link NaN_Strategy} is {@code FAILED}
	 * and a {@link Double#NaN} is encountered in the input data
	 */
	 //override
	std::vector<double> rank(std::vector<double> data)
	{
		// Array recording initial positions of data to be ranked
		auto ranks = std::vector<Int_Double_Pair>(data.size());
		for (int i{}; i < data.size(); i++)
		{
			ranks[i] = Int_Double_Pair(data[i], i);
		}

		// Recode, remove or record positions of NaNs
		std::vector<int> nan_positions = NULL;
		switch (my_nan_strategy)
		{
		case MAXIMAL: // Replace NaNs with +INFs
			recode_nans(ranks, INFINITY);
			break;
		case MINIMAL: // Replace NaNs with -INFs
			recode_nans(ranks, -INFINITY);
			break;
		case REMOVED: // Drop NaNs from data
			ranks = remove_na_ns(ranks);
			break;
		case FIXED:   // Record positions of NaNs
			nan_positions = get_nan_positions(ranks);
			break;
		case FAILED:
			nan_positions = get_nan_positions(ranks);
			if (!nan_positions.is_empty())
			{
				throw std::exception("not implemented");
				//throw (hipparchus::exception::Localized_Core_Formats_Type::NAN_NOT_ALLOWED);
			}
			break;
		default: // this should not happen unless NaN_Strategy enum is changed
			throw std::exception("not implemented");
			//throw Math_Runtime_Exception.create_internal_error();
		}

		// Sort the Int_Double_Pairs
		Arrays.sort(ranks, (p1, p2)->Double.compare(p1.value, p2.value));

		// Walk the sorted array, filling output array using sorted positions, // resolving ties as we go
		std::vector<double> out = std::vector<double>(ranks.size()];
		int pos{ 1 };  // position in sorted array
		out[ranks[0].get_position()] = pos;
		auto ties_trace = std::vector<int>();
		ties_trace.push_back(ranks[0].get_position());
		for (int i{ 1 }; i < ranks.size(); i++)
		{
			if (Double.compare(ranks[i].get_value(), ranks[i - 1].get_value()) > 0)
			{
				// tie sequence has ended (or had length 1)
				pos = i + 1;
				if (ties_trace.size() > 1) {  // if seq is nontrivial, resolve
					resolve_tie(out, ties_trace);
				}
				ties_trace = Array_list<>();
				ties_trace.add(ranks[i].get_position());
			}
			else
			{
				// tie sequence continues
				ties_trace.add(ranks[i].get_position());
			}
			out[ranks[i].get_position()] = pos;
		}
		if (ties_trace.size() > 1) {  // handle tie sequence at end
			resolve_tie(out, ties_trace);
		}
		if (nan_strategy == NaN_Strategy::FIXED)
		{
			restore_na_ns(out, nan_positions);
		}
		return out;
	}

private:
	/**
	 * Returns an array that is a copy of the input array with Int_Double_Pairs
	 * having NaN values removed.
	 *
	 * @param ranks input array
	 * @return array with NaN-valued entries removed
	 */
	std::vector<Int_Double_Pair> remove_na_ns(std::vector<Int_Double_Pair>& ranks)
	{
		if (!contains_nans(ranks))
		{
			return ranks;
		}
		auto out_ranks = std::vector<Int_Double_Pair>(ranks.size());
		int j{};
		for (int i{}; i < ranks.size(); i++)
		{
			if (std::isnan(ranks[i].get_value()))
			{
				// drop, but adjust original ranks of later elements
				for (int k{ i + 1 }; k < ranks.size(); k++)
				{
					ranks[k] = Int_Double_Pair(ranks[k].get_value(), ranks[k].get_position() - 1);
				}
			}
			else
			{
				out_ranks[j] = Int_Double_Pair(
					ranks[i].get_value(), ranks[i].get_position());
				j++;
			}
		}
		auto return_ranks = std::vector<Int_Double_Pair>(j);
		System.arraycopy(out_ranks, 0, return_ranks, 0, j);
		return return_ranks;
	}

	/**
	 * Recodes NaN values to the given value.
	 *
	 * @param ranks array to recode
	 * @param value the value to replace NaNs with
	 */
	void recode_nans(const std::vector<Int_Double_Pair>& ranks, const double& value)
	{
		for (int i{}; i < ranks.size(); i++)
		{
			if (std::isnan(ranks[i].get_value()))
			{
				ranks[i] = Int_Double_Pair(value, ranks[i].get_position());
			}
		}
	}

	/**
	 * Checks for presence of NaNs in <code>ranks.</code>
	 *
	 * @param ranks array to be searched for NaNs
	 * @return true iff ranks contains one or more NaNs
	 */
	bool contains_nans(const std::vector<Int_Double_Pair>& ranks)
	{
		for (int i{}; i < ranks.size(); i++)
		{
			if (std::isnan(ranks[i].get_value()))
			{
				return true;
			}
		}
		return false;
	}

	/**
	 * Resolve a sequence of ties, using the configured {@link Ties_Strategy}.
	 * The input <code>ranks</code> array is expected to take the same value
	 * for all indices in <code>ties_trace</code>.  The common value is recoded
	 * according to the ties_strategy. For example, if ranks = <5,8,2,6,2,7,1,2>, * ties_trace = <2,4,7> and ties_strategy is MINIMUM, ranks will be unchanged.
	 * The same array and trace with ties_strategy AVERAGE will come out
	 * <5,8,3,6,3,7,1,3>.
	 *
	 * @param ranks array of ranks
	 * @param ties_trace list of indices where <code>ranks</code> is constant
	 * -- that is, for any i and j in Ties_Trace, <code> ranks[i] == ranks[j]
	 * </code>
	 */
	void resolve_tie(const std::vector<double>& ranks, std::vector<int>& ties_trace)
	{
		// constant value of ranks over ties_trace
		const double c = ranks[ties_trace.at(0)];

		// length of sequence of tied ranks
		const int length = ties_trace.size();

		switch (ties_strategy)
		{
		case  AVERAGE:  // Replace ranks with average
			fill(ranks, ties_trace, (2 * c + length - 1) / 2.0);
			break;
		case MAXIMUM:   // Replace ranks with maximum values
			fill(ranks, ties_trace, c + length - 1);
			break;
		case MINIMUM:   // Replace ties with minimum
			fill(ranks, ties_trace, c);
			break;
		case RANDOM:    // Fill with random integral values in [c, c + length - 1]
			auto iterator = ties_trace.iterator();
			long f = std::round(c);
			while (iterator.has_next())
			{
				// No advertised exception because args are guaranteed valid
				ranks[iterator.next()] =
					random_data.next_long(f, f + length - 1);
			}
			break;
		case SEQUENTIAL:  // Fill sequentially from c to c + length - 1
			// walk and fill
			iterator = ties_trace.iterator();
			f = std::round(c);
			int i = 0;
			while (iterator.has_next())
			{
				ranks[iterator.next()] = f + i++;
			}
			break;
		default: // this should not happen unless Ties_Strategy enum is changed
			throw Math_Runtime_Exception.create_internal_error();
		}
	}

	/**
	 * Sets<code>data[i] = value</code> for each i in <code>ties_trace.</code>
	 *
	 * @param data array to modify
	 * @param ties_trace list of index values to set
	 * @param value value to set
	 */
	void fill(const std::vector<double>& data, std::vector<int>& ties_trace, const double& value)
	{
		auto iterator = ties_trace.iterator();
		while (iterator.has_next())
		{
			data[iterator.next()] = value;
		}
	}

	/**
	 * Set <code>ranks[i] =NAN</code> for each i in <code>nan_positions.</code>
	 *
	 * @param ranks array to modify
	 * @param nan_positions list of index values to set to <code>Double.NaN</code>
	 */
	void restore_na_ns(std::vector<double>& ranks, std::vector<int>& nan_positions)
	{
		if (nan_positions.is_empty())
		{
			return;
		}
		auto iterator = nan_positions.iterator();
		while (iterator.has_next())
		{
			ranks[iterator.next().int_value()] = std::numeric_limits<double>::quiet_NaN();
		}
	}

	/**
	 * Returns a list of indexes where <code>ranks</code> is <code>NaN.</code>
	 *
	 * @param ranks array to search for <code>NaNs</code>
	 * @return list of indexes i such that <code>ranks[i] = NaN</code>
	 */
	std::vector<int> get_nan_positions(const std::vector<Int_Double_Pair>& ranks)
	{
		auto out = std::vector<int>();
		const auto ranks_size{ ranks.size() };
		for (int i{}; i < ranks_size; i++)
		{
			if (std::isnan(ranks[i].get_value()))
			{
				out.add(static_cast<int>(i));
			}
		}
		return out;
	}
};

/**
 * Represents the position of a double value in an ordering.
 */
class Int_Double_Pair
{
private:
	/** Value of the pair */
	const double my_value;

	/** Original position of the pair */
	const int my_position;

public:
	/**
	 * Construct an Int_Double_Pair with the given value and position.
	 * @param value the value of the pair
	 * @param position the original position
	 */
	Int_Double_Pair(const double& value, const int& position)
		:
		my_value{ value },
		my_position{ position }
	{};

	/**
	 * Returns the value of the pair.
	 * @return value
	 */
	double get_value() const
	{
		return my_value;
	}

	/**
	 * Returns the original position of the pair.
	 * @return position
	 */
	int get_position() const
	{
		return my_position;
	}
};