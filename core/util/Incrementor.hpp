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

#include <limits>
#include "../exception/LocalizedCoreFormats.h"
#include "../../migration/exception/MaxCountExceededException.hpp"

  /**
   * Utility that increments a counter until a maximum is reached, at
   * which point, the instance will by default throw a
   * {@link Math_Illegal_State_Exception}.
   * However, the user is able to //override this behaviour by defining a
   * custom {@link Max_Count_Exceeded_Callback callback}, in order to e.g.
   * select which exception must be thrown.
   */
template<
	typename T, //real type
	typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type
>
class Incrementor
{
private:
	/** Default callback. */
	static const Max_Count_Exceeded_Callback<T> DEFAULT_CALLBACK = (const int& max)
	{
		throw std::exception("not implemented");
		//throw Math_Illegal_State_Exception(hipparchus::exception::Localized_Core_Formats_Type::MAX_COUNT_EXCEEDED, max);
	};

	/** Upper limit for the counter. */
	const int my_maximal_count;
	/** Function called at counter exhaustion. */
	const Max_Count_Exceeded_Callback<T> my_max_count_callback;
	/** Current count. */
	int my_count;

public:
	/**
	 * Defines a method to be called at counter exhaustion.
	 * The {@link #triggerstatic_cast<int>( trigger} method should usually throw an exception.
	 */
	class Max_Count_Exceeded_Callback
	{
		/**
		 * Function called when the maximal count has been reached.
		 *
		 * @param maximal_count Maximal count.
		 * @Math_Illegal_State_Exception at counter exhaustion
		 */
		void trigger(const int& maximal_count);
	};

	/**
	 * Creates an Incrementor.
	 * <p>
	 * The maximal value will be set to {@code std::numeric_limits<int>::max()}.
	 */
	Incrementor<T>()
	{
		Incrementor<T>(std::numeric_limits<T>::max());
	}

	/**
	 * Creates an Incrementor.
	 *
	 * @param max Maximal count.
	 * @ if {@code max} is negative.
	 */
	Incrementor<T>(const T& max)
	{
		Incrementor<T>(max, DEFAULT_CALLBACK);
	}

	/**
	 * Creates an Incrementor.
	 *
	 * @param max Maximal count.
	 * @param cb Function to be called when the maximal count has been reached.
	 * @Null_Argument_Exception if {@code cb} is {@code NULL}.
	 * @ if {@code max} is negative.
	 */
	Incrementor<T>(const int& max, const Max_Count_Exceeded_Callback<T>& cb)
	{
		Incrementor<T>(0, max, cb);
	}

	/**
	 * Creates an Incrementor.
	 *
	 * @param count Initial counter value.
	 * @param max Maximal count.
	 * @param cb Function to be called when the maximal count has been reached.
	 * @Null_Argument_Exception if {@code cb} is {@code NULL}.
	 * @ if {@code max} is negative.
	 */
	Incrementor<T>(const int& count, const int& max, const Max_Count_Exceeded_Callback<T>& cb)
		:
		my_maximal_count{max},
		my_max_count_callback{cb},
		my_count{count}
	{
		throw std::exception("not implemented");
		//if (cb == NULL)
		//{
		//	throw Null_Argument_Exception();
		//}
		//if (max < 0)
		//{
		//	throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL, max, 0);
		//}
	}

	/**
	 * Creates a instance and set the counter to the given value.
	 *
	 * @param value Value of the counter.
	 * @return a instance.
	 */
	Incrementor<T> with_count(const int& value)
	{
		return Incrementor<T>(value, my_maximal_count, my_max_count_callback);
	}

	/**
	 * Creates a instance with a given maximal count.
	 * The counter is reset to 0.
	 *
	 * @param max Maximal count.
	 * @return a instance.
	 * @ if {@code max} is negative.
	 */
	Incrementor<T> with_maximal_count(const int& max)
	{
		return Incrementor<T>(0, max, my_max_count_callback);
	}

	/**
	 * Creates a instance with a given callback.
	 * The counter is reset to 0.
	 *
	 * @param cb Callback to be called at counter exhaustion.
	 * @return a instance.
	 */
	Incrementor<T> with_callback(const Max_Count_Exceeded_Callback<T>& cb)
	{
		return Incrementor<T>(0, my_maximal_count, cb);
	}

	/**
	 * Gets the upper limit of the counter.
	 *
	 * @return the counter upper limit.
	 */
	int get_maximal_count() const
	{
		return my_maximal_count;
	}

	/**
	 * Gets the current count.
	 *
	 * @return the current count.
	 */
	int get_count() const
	{
		return my_count;
	}

	/**
	 * Checks whether incrementing the counter {@code n_times} is allowed.
	 *
	 * @return {@code false} if calling {@link #increment()}
	 * will trigger a {@code Math_Illegal_State_Exception}, * {@code true} otherwise.
	 */
	bool can_increment()
	{
		return can_increment(1);
	}

	/**
	 * Checks whether incrementing the counter several times is allowed.
	 *
	 * @param n_times Number of increments.
	 * @return {@code false} if calling {@link #incrementstatic_cast<int>(
	 * increment(n_times)} would call the {@link Max_Count_Exceeded_Callback callback}
	 * {@code true} otherwise.
	 * @ if {@code n_times} is negative.
	 */
	bool can_increment(const int& n_times)
	{
		if (n_times < 0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL, n_times, 0);
		}
		return my_count <= my_maximal_count - n_times;
	}

	/**
	 * Performs multiple increments.
	 *
	 * @param n_times Number of increments.
	 * @ if {@code n_times} is negative.
	 *
	 * @see #increment()
	 */
	void increment(const int& n_times)
	{
		if (n_times < 0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL, n_times, 0);
		}

		for (int i{}; i < n_times; i++)
		{
			increment();
		}
	}

	/**
	 * Adds the increment value to the current iteration count.
	 * At counter exhaustion, this method will call the
	 * {@link Max_Count_Exceeded_Callback#triggerstatic_cast<int>( trigger} method of the
	 * callback object passed to the
	 * {@link #with_callback(Max_Count_Exceeded_Callback)} method.
	 *
	 * @see #incrementstatic_cast<int>(
	 */
	void increment()
	{
		if (my_count > my_maximal_count - 1)
		{
			my_max_count_callback.trigger(my_maximal_count);
		}
		++my_count;
	}

	/**
	 * Resets the counter to 0.
	 */
	void reset()
	{
		my_count = 0;
	}
};