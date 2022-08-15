#pragma once
/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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
//package org.hipparchus.stat.descriptive.rank;

//import java.io.Serializable;
//import java.util.Array_list;
//import java.util.Arrays;
//import java.util.Collection;
//import java.util.Hash_Map;
//import java.util.Iterator;
//import java.util.List;
//import java.util.Map;
//import java.util.No_Such_Element_Exception;
//import java.util.UUID;


//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.exception.Null_Argument_Exception;
//import org.hipparchus.random.Random_Generator;
//import org.hipparchus.random.Well19937c;
//import org.hipparchus.stat.Stat_Utils;
//import org.hipparchus.stat.descriptive.Abstract_Storeless_Univariate_Statistic;
//import org.hipparchus.stat.descriptive.Aggregatable_Statistic;
//import org.hipparchus.stat.descriptive.Storeless_Univariate_Statistic;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Arrays;
#include <cmath>
#include <algorithm>
#include <vector>
#include "Percentile.h"
#include "../../../core/random/RandomGenerator.h"
#include "../AggregatableStatistic.hpp"
#include "../StorelessUnivariateStatistic.h"
#include "../AbstractStorelessUnivariateStatistic.h"
#include "../../StatUtils.h"

/**
 * A {@link Storeless_Univariate_Statistic} estimating percentiles using the
 * <a href=http://dimacs.rutgers.edu/~graham/pubs/papers/nquantiles.pdf>RANDOM</a>
 * Algorithm.
 * <p>
 * Storage requirements for the RANDOM algorithm depend on the desired
 * accuracy of quantile estimates. Quantile estimate accuracy is defined as follows.
 * <p>
 * Let \(X\) be the set of all data values consumed from the stream and let \(q\)
 * be a quantile (measured between 0 and 1) to be estimated. If <ul>
 * <li>\(\epsilon\) is the configured accuracy</li>
 * <li> \(\hat{q}\) is a Random_Percentile estimate for \(q\) (what is returned
 *      by {@link #get_result()} or {@link #get_resultstatic_cast<double>(}) with \(100q\) as
 *      actual parameter)</li>
 * <li> \(rank(\hat{q}) = |\{x \in X : x &lt; \hat{q}\}|\) is the actual rank of
 *      \(\hat{q}\) in the full data stream</li>
 * <li>\(n = |X|\) is the number of observations</li></ul>
 * then we can expect \((q - \epsilon)n &lt; rank(\hat{q}) &lt; (q + \epsilon)n\).
 * <p>
 * The algorithm maintains \(\left\lceil {log_{2}(1/\epsilon)}\right\rceil + 1\) buffers
 * of size \(\left\lceil {1/\epsilon \sqrt{log_2(1/\epsilon)}}\right\rceil\).  When
 * {@code epsilon} is set to the default value of \(10^{-4}\), this makes 15 buffers
 * of size 36,453.
 * <p>
 * The algorithm uses the buffers to maintain samples of data from the stream.  Until
 * all buffers are full, the entire sample is stored in the buffers.
 * If one of the {@code get_result} methods is called when all data are available in memory
 * and there is room to make a copy of the data (meaning the combined set of buffers is
 * less than half full), the {@code get_result} method delegates to a {@link Percentile}
 * instance to compute and return the exact value for the desired quantile.
 * For default {@code epsilon}, this means exact values will be returned whenever fewer than
 * \(\left\lceil {15 	imes 36453 / 2} \right\rceil = 273,398\) values have been consumed
 * from the data stream.
 * <p>
 * When buffers become full, the algorithm merges buffers so that they effectively represent
 * a larger set of values than they can hold. Subsequently, data values are sampled from the
 * stream to fill buffers freed by merge operations.  Both the merging and the sampling
 * require random selection, which is done using a {@code Random_Generator}.  To get
 * repeatable results for large data streams, users should provide {@code Random_Generator}
 * instances with fixed seeds. {@code Random_Percentile} itself does not reseed or otherwise
 * initialize the {@code Random_Generator} provided to it.  By default, it uses a
 * {@link Well19937c} generator with the default seed.
 * <p>
 * Note: This implementation is not thread-safe.
 */
class Random_Percentile
    :
    Abstract_Storeless_Univariate_Statistic,
    Storeless_Univariate_Statistic,
    Aggregatable_Statistic<Random_Percentile>
{

public:
    /** Default quantile estimation error setting */
    static constexpr double DEFAULT_EPSILON{ 1e-4 };
    /** Serialization version id */
 
private:
    /** Storage size of each buffer */
    int my_s;
    /** Maximum number of buffers minus 1 */
    int my_h;
    /** Data structure used to manage buffers */
    Buffer_Map my_buffer_map;
    /** Bound on the quantile estimation error */
    double epsilon;
    /** Source of random data */
    Random_Generator my_random_generator;
    /** Number of elements consumed from the input data stream */
    long my_n;
    /** Buffer currently being filled */
    Buffer my_current_buffer;

    /**
     * Computes base 2 log of the argument.
     *
     * @param x input value
     * @return the value y such that 2^y = x
     */
    static double log2(double x)
    {
        return std::log(x) / std::log(2);
    }

public:
    /**
     * Constructs a {@code Random_Percentile} with quantile estimation error
     * {@code epsilon} using {@code random_generator} as its source of random data.
     *
     * @param epsilon bound on quantile estimation error (see class javadoc)
     * @param random_generator PRNG used in sampling and merge operations
     * @ if percentile is not in the range [0, 100]
     */
    Random_Percentile(const double& epsilon, const Random_Generator& random_generator) 
    {
        if (epsilon <= 0) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL, epsilon, 0);
        }
        my_h = static_cast<int>( std::ceil(log2(1/epsilon));
        my_s = static_cast<int>( std::ceil(std::sqrt(log2(1/epsilon)) / epsilon);
        my_random_generator = random_generator;
        my_buffer_map = Buffer_Map(my_h + 1, my_s, random_generator);
        my_current_buffer = my_buffer_map.create(0);
        my_epsilon = epsilon;
    }

    /**
     * Constructs a {@code Random_Percentile} with default estimation error
     * using {@code random_generator} as its source of random data.
     *
     * @param random_generator PRNG used in sampling and merge operations
     * @ if percentile is not in the range [0, 100]
     */
    Random_Percentile(const Random_Generator& random_generator) 
    {
        Random_Percentile(DEFAULT_EPSILON, random_generator);
    }

    /**
     * Constructs a {@code Random_Percentile} with quantile estimation error
     * {@code epsilon} using the default PRNG as source of random data.
     *
     * @param epsilon bound on quantile estimation error (see class javadoc)
     * @ if percentile is not in the range [0, 100]
     */
    Random_Percentile(const double& epsilon) 
    {
        Random_Percentile(epsilon, Well19937c());
    }

    /**
     * Constructs a {@code Random_Percentile} with quantile estimation error
     * set to the default ({@link #DEFAULT_EPSILON}), using the default PRNG
     * as source of random data.
     */
    Random_Percentile() 
    {
        Random_Percentile(DEFAULT_EPSILON, Well19937c());
    }

    /**
     * Copy constructor, creates a {@code Random_Percentile} identical
     * to the {@code original}.  Note: the Random_Generator used by the new
     * instance is referenced, not copied - i.e., the instance shares
     * a generator with the original.
     *
     * @param original the {@code P_Square_Percentile} instance to copy
     */
    Random_Percentile(Random_Percentile original) 
    {
        super();
        my_h = original.h;
        my_n = original.n;
        my_s = original.s;
        my_epsilon = original.epsilon;
        my_buffer_map = Buffer_Map(original.buffer_map);
        my_random_generator = original.random_generator;
        Iterator<Buffer> iterator = my_buffer_map.iterator();
        Buffer current = NULL;
        Buffer curr = NULL;
        // See if there is a partially filled buffer - that will be my_current_buffer
        while (current == NULL && iterator.has_next()) 
        {
            curr = iterator.next();
            if (curr.has_capacity()) 
            {
                current = curr;
            }
        }
        // If there is no partially filled buffer, just assign the last one.
        // Next increment() will find no capacity and create a one or trigger
        // a merge.
        my_current_buffer = current == NULL
            ? curr
            : current;
    }

    //override
    long get_n() 
    {
        return n;
    }

    /**
     * Returns an estimate of the given percentile, computed using the designated
     * array segment as input data.
     *
     * @param values source of input data
     * @param begin position of the first element of the values array to include
     * @param length number of array elements to include
     * @param percentile desired percentile (scaled 0 - 100)
     *
     * @return estimated percentile
     * @ if percentile is out of the range [0, 100]
     */
    double evaluate(const double& percentile, const std::vector<double>& values, const int& begin, const int& length)
    {
        if (Math_Arrays::verify_values(values, begin, length)) 
        {
            Random_Percentile random_percentile = Random_Percentile(my_epsilon, my_random_generator);
            random_percentile.increment_all(values, begin, length);
            return random_percentile.get_result(percentile);
        }
        return std::numeric_limits<double>::quiet_NaN();
    }

    /**
     * Returns an estimate of the median, computed using the designated
     * array segment as input data.
     *
     * @param values source of input data
     * @param begin position of the first element of the values array to include
     * @param length number of array elements to include
     *
     * @return estimated percentile
     * @ if percentile is out of the range [0, 100]
     */
    //override
    double evaluate(const std::vector<double>& values, const int& begin, const int& length) 
    {
        return evaluate(50, values, begin, length);
    }

    /**
     * Returns an estimate of percentile over the given array.
     *
     * @param values source of input data
     * @param percentile desired percentile (scaled 0 - 100)
     *
     * @return estimated percentile
     * @ if percentile is out of the range [0, 100]
     */
    double evaluate(const double& percentile, const std::vector<double>& values) 
    {
        return evaluate(percentile, values, 0, values.size());
    }

    //override
    Random_Percentile copy() 
    {
        return Random_Percentile(*this);
    }

    //override
    void clear() 
    {
        my_n = 0;
        my_buffer_map.clear();
        my_current_buffer = my_buffer_map.create(0);
    }

    /**
     * Returns an estimate of the median.
     */
    //override
    double get_result() 
    {
        return get_result(50);

    }

    /**
     * Returns an estimate of the given percentile.
     *
     * @param percentile desired percentile (scaled 0 - 100)
     * @return estimated percentile
     * @ if percentile is out of the range [0, 100]
     */
    double get_result(const double& percentile) 
    {
        if (percentile > 100 || percentile < 0) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::OUT_OF_RANGE, percentile, 0, 100);
        }
        // Convert to internal quantile scale
        const double q = percentile / 100;
        // First get global min and max to bound search.
        double min = INFINITY;
        double max = -INFINITY;
        double b_min;
        double b_max;
        Iterator<Buffer> buffer_iterator = my_buffer_map.iterator();
        while (buffer_iterator.has_next()) 
        {
            Buffer buffer = buffer_iterator.next();
            b_min = buffer.min();
            if (b_min < min) 
            {
                min = b_min;
            }
            b_max = buffer.max();
            if (b_max > max) 
            {
                max = b_max;
            }
        }

        // Handle degenerate cases
        if (Double.compare(q, 0.0) == 0 || my_n == 1) 
        {
            return min;
        }
        if (Double.compare(q, 1) == 0) 
        {
            return max;
        }
        if (my_n == 0) 
        {
            return std::numeric_limits<double>::quiet_NaN();
        }

        // See if we have all data in memory and enough free memory to copy.
        // If so, use Percentile to perform exact computation.
        if (my_buffer_map.half_empty()) 
        {
            return Percentile(percentile).evaluate(my_buffer_map.level_zero_data());
        }

        // Compute target rank
        const double target_rank = q * my_n;

        // Start with initial guess min + quantile * (max - min).
        double estimate = min + q * (max - min);
        double estimate_rank = get_rank(estimate);
        double lower;
        double upper;
        if (estimate_rank > target_rank) 
        {
            upper = estimate;
            lower = min;
        }
        else if (estimate_rank < target_rank) 
        {
            lower = estimate;
            upper = max;
        }
        else 
        {
            return estimate;
        }
        const double eps = epsilon / 2;
        const double rank_tolerance = eps * my_n;
        const double min_width = eps / my_n;
        double interval_width = std::abs(upper - lower);
        while (std::abs(estimate_rank - target_rank) > rank_tolerance && interval_width > min_width) 
        {
            if (estimate_rank > target_rank) 
            {
                upper = estimate;
            }
            else 
            {
                lower = estimate;
            }
            interval_width = upper - lower;
            estimate = lower + interval_width / 2;
            estimate_rank = get_rank(estimate);
        }
        return estimate;
    }

    /**
     * Gets the estimated rank of {@code value}, i.e.  \(|\{x \in X : x &lt; value\}|\)
     * where \(X\) is the set of values that have been consumed from the stream.
     *
     * @param value value whose overall rank is sought
     * @return estimated number of sample values that are strictly less than {@code value}
     */
    double get_rank(double value) 
    {
        double rank_sum = 0;
        Iterator<Buffer> buffer_iterator = my_buffer_map.iterator();
        while (buffer_iterator.has_next()) 
        {
            Buffer buffer = buffer_iterator.next();
            rank_sum += buffer.rank_of(value) * std::pow(2, buffer.level);
        }
        return rank_sum;
    }

    /**
     * Returns the estimated quantile position of value in the dataset.
     * Specifically, what is returned is an estimate of \(|\{x \in X : x &lt; value\}| / |X|\)
     * where \(X\) is the set of values that have been consumed from the stream.
     *
     * @param value value whose quantile rank is sought.
     * @return estimated proportion of sample values that are strictly less than {@code value}
     */
    double get_quantile_rank(const double& value) 
    {
        return get_rank(value) / get_n();
    }

    //override
    void increment(const double& d) 
    {
        my_n++;
        if (!my_current_buffer.has_capacity())
        { // Need to get a buffer to fill
            // First see if we have not yet created all the buffers
            if (my_buffer_map.can_create()) 
            {
                const int level = static_cast<int>( std::ceil(std::max(0, log2(my_n/(my_s * std::pow(2, my_h - 1))))));
                my_current_buffer = my_buffer_map.create(level);
            }
        else
        { // All buffers have been created - need to merge to free one
                my_current_buffer = my_buffer_map.merge();
            }
        }
        my_current_buffer.consume(d);
    }

    /**
     * Maintains a buffer of values sampled from the input data stream.
     * <p>
     * The {@link #level} of a buffer determines its sampling frequency.
     * The {@link #consumestatic_cast<double>(} method retains 1 out of every 2^level values
     * read in from the stream.
     * <p>
     * The {@link #size} of the buffer is the number of values that it can store
     * The buffer is considered full when it has consumed 2^level * size values.
     * <p>
     * The {@link #block_size} of a buffer is 2^level.
     * The consume method starts each block by generating a random integer in
     * [0, block_size - 1].  It then skips over all but the element with that offset
     * in the block, retaining only the selected value. So after 2^level * size
     * elements have been consumed, it will have retained size elements - one
     * from each 2^level block.
     * <p>
     * The {@link #merge_with(Buffer)} method merges this buffer with another one, * The merge operation merges the data from the other buffer into this and clears
     * the other buffer (so it can consume data). Both buffers have their level
     * incremented by the merge. This operation is only used on full buffers.
     */
    static class Buffer  
    {
    private:
        /** Number of values actually stored in the buffer */
        const int my_size;
        /** Data sampled from the stream */
        std::vector<double> my_data;
        /** PRNG used for merges and stream sampling */
        const Random_Generator my_random_generator;
        /** Level of the buffer */
        int my_level;
        /** Block size  = 2^level */
        long my_block_size;
        /** Next location in backing array for stored (taken) value */
        int my_next;
        /** Number of values consumed in current 2^level block of values from the stream */
        long my_consumed;
        /** Index of next value to take in current 2^level block */
        long my_next_to_take;
        /** ID */
        const UUID my_id;

        /**
         * Sets block_size and next_to_take based on level.
         */
        void compute_block_size()
        {
            if (my_level == 0)
            {
                my_block_size = 1;
            }
            else
            {
                long product = 1;
                for (int i{}; i < my_level; i++)
                {
                    product *= 2;
                }
                my_block_size = product;
            }
            if (my_block_size > 1)
            {
                my_next_to_take = my_random_generator.next_long(my_block_size);
            }
        }

    public:
        /**
         * Creates a buffer capable of retaining size values with the given level.
         *
         * @param size number of values the buffer can retain
         * @param level the base 2 log of the sampling frequency
         *        (one out of every 2^level values is retained)
         * @param random_generator PRNG used for sampling and merge operations
         */
        Buffer(const int& size, const int& level, const Random_Generator& random_generator) 
            :
            my_size{ size },
            my_data{ std::vector<double>(size) },
            my_level{ level },
            my_random_generator{ random_generator },
            my_id{ UUID.random_uuid() }
        {
            compute_block_size();
        }

        /**
         * Consumes a value from the input stream.
         * <p>
         * For each 2^level values consumed, one is added to the buffer.
         * The buffer is not considered full until 2^level * size values
         * have been consumed.
         * <p>
         * Sorts the data array if the consumption renders the buffer full.
         * <p>
         * There is no capacity check in this method.  Clients are expected
         * to use {@link #has_capacity()} before invoking this method.  If
         * it is invoked on a full buffer, an Array_indexOutOfbounds exception
         * will result.
         *
         * @param value value to consume from the stream
         */
        void consume(const double& value) 
        {
            if (my_consumed == my_next_to_take) 
            {
                my_data[my_next] = value;
                my_next++;
            }
            consumed++;
            if (consumed == my_block_size) 
            {
                if (my_next == my_size) 
                {   // Buffer is full
                    Arrays.sort(my_data);
                }
                else
                {              // Reset in-block counter and next_to_take
                    consumed = 0;
                    if (my_block_size > 1) 
                    {
                        my_next_to_take = my_random_generator.next_long(my_block_size);
                    }
                }
            }
        }

        /**
         * Merges this with other.
         * <p>
         * After the merge, this will be the merged buffer and other will be free.
         * Both will have level+1.
         * Post-merge, other can be used to accept data.
         * <p>
         * The contents of the merged buffer (this after the merge) are determined
         * by randomly choosing one of the two retained elements in each of the
         * [0...size - 1] positions in the two buffers.
         * <p>
         * This and other must have the same level and both must be full.
         *
         * @param other initially full other buffer at the same level as this.
         * @ if either buffer is not full or they
         * have different levels
         */
        void merge_with(const Buffer& other) 
        {
            // Make sure both this and other are full and have the same level
            if (has_capacity() || other.has_capacity() || other.level != my_level) 
            {
                throw std::exception("not implemented");
                //throw (hipparchus::exception::Localized_Core_Formats_Type::INTERNAL_ERROR);
            }
            // Randomly select one of the two entries for each slot
            for (int i{}; i < my_size; i++) 
            {
                if (my_random_generator.next_boolean()) 
                {
                    my_data[i] = other.data[i];
                }
            }
            // Re-sort data
            Arrays.sort(my_data);
            // Bump level of both buffers
            other.set_level(my_level + 1);
            set_level(my_level + 1);
            // Clear the free one (and compute blocksize)
            other.clear();
        }

        /**
         * Merge this into a higher-level buffer.
         * <p>
         * Does not alter this; but after the merge, higher may have some of its
         * data replaced by data from this.  Levels are not changed for either buffer.
         * <p>
         * Probability of selection into the newly constituted higher buffer is weighted
         * according to level. So for example, if this has level 0 and higher has level
         * 2, the ith element of higher is 4 times more likely than the corresponding
         * element of this to retain its spot.
         * <p>
         * This method is only used when aggregating Random_Percentile instances.
         * <p>
         * Preconditions:
         * <ol><li> this.level < higher.level </li>
         *     <li> this.size = higher.size </li>
         *     <li> Both buffers are full </li>
         * </ol>
         *
         * @param higher higher-level buffer to merge this into
         * @ if the buffers have different sizes, * either buffer is not full or this has level greater than or equal to higher
         */
        void merge_into(const Buffer& higher) 
        {
            // Check preconditions
            if (my_size != higher.my_size || has_capacity() || higher.has_capacity() || my_level >= higher.my_level) 
            {
                throw std::exception("not implemented");
                //throw (hipparchus::exception::Localized_Core_Formats_Type::INTERNAL_ERROR);
            }
            const int level_difference = higher.my_level - my_level;
            int m{ 1 };
            for (int i{}; i < level_difference; i++) 
            {
                m *= 2;
            }
            // Randomly select one of the two entries for each slot in higher, giving
            // m-times higher weight to the entries of higher.
            for (int i{}; i < my_size; i++) 
            {
                // data[i] <-> {0}, higher.data[i] <-> {1, ..., m}
                if (my_random_generator.next_int(m + 1) == 0) 
                {
                    higher.data[i] = my_data[i];
                }
            }
            // Resort higher's data
            Arrays.sort(higher.my_data);
        }

        /**
         * @return true if the buffer has capacity; false if it is full
         */
        bool has_capacity() 
        {
            // Buffer has capacity if it has not yet set all of its data
            // values or if it has but still has not finished its last block
            return next < my_size || consumed < my_block_size;
        }

        /**
         * Sets the level of the buffer.
         *
         * @param level level value
         */
        void set_level(const int& level) 
        {
            my_level = level;
        }

        /**
         * Clears data, recomputes block_size and resets consumed and next_to_take.
         */
        void clear() 
        {
            my_consumed = 0;
            my_next = 0;
            compute_block_size();
        }

        /**
         * Returns a copy of the data that has been added to the buffer
         *
         * @return possibly unsorted copy of the portion of the buffer that has been filled
         */
        std::vector<double> get_data() 
        {
            auto out = std::vector<double>(my_next);
            System.arraycopy(my_data, 0, out, 0, my_next);
            return out;
        }

        /**
         * Returns the ordinal rank of value among the sampled values in this buffer.
         *
         * @param value value whose rank is sought
         * @return |{v in data : v < value}|
         */
        int rank_of(const double& value) 
        {
            int ret{};
            if (!has_capacity()) { // Full sorted buffer, can do binary search
                ret = Arrays.binary_search(my_data, value);
                return ret < 0
                    ? -ret - 1
                    : ret;
            }
            else { // have to count - not sorted yet and can't sort yet
                for (int i{}; i < my_next; i++) 
                {
                    if (my_data[i] < value) 
                    {
                        ret++;
                    }
                }
                return ret;
            }
        }

        /**
         * @return the smallest value held in this buffer
         */
        double min() 
        {
            return !has_capacity() 
                ? my_data[0]
                : Stat_Utils.min(get_data());
        }

        /**
         * @return the largest value held in this buffer
         */
        double max() 
        {
            return !has_capacity()
                ? my_data[my_data.size() - 1]
                : Stat_Utils.max(get_data());
        }

        /**
         * @return the level of this buffer
         */
        int get_level() const
        {
            return my_level;
        }

        /**
         * @return the id
         */
        UUID get_id() const 
        {
            return my_id;
        }
    };

    /**
     * A map structure to hold the buffers.
     * Keys are levels and values are lists of buffers at the given level.
     * Overall capacity is limited by the total number of buffers.
     */
    static class Buffer_Map //: Iterable<Buffer>
    {
    private:
        /** Total number of buffers that can be created - cap for count */
        const int my_capacity;
        /** PRNG used in merges */
        const Random_Generator my_random_generator;
        /** Total count of all buffers */
        int my_count;
        /** Uniform buffer size */
        const int my_buffer_size;
        /** Backing store for the buffer map. Keys are levels, values are lists of registered buffers. */
        const std::unordered_map<int, std::vector<Buffer>> my_registry;
        /** Maximum buffer level */
        int my_max_level;
    public:
        /**
         * Creates a Buffer_Map that can manage up to capacity buffers.
         * Buffers created by the pool with have size = buffersize.
         *
         * @param capacity cap on the number of buffers
         * @param buffer_size size of each buffer
         * @param random_generator Random_Generator to use in merges
         */
        Buffer_Map(const int& capacity, const int& buffer_size, const Random_Generator& random_generator) 
            :
            my_buffer_size{ buffer_size },
            my_capacity{ capacity },
            my_random_generator{ random_generator },
            my_registry{ std::unordered_map<int, std::vector<Buffer>> }
        {
        }

        /**
         * Copy constructor.
         *
         * @param original Buffer_Map to copy
         */
        Buffer_Map(const Buffer_Map& original)
            : my_buffer_size{  original }
        {
            super();
            this.buffer_size = original.buffer_size;
            this.capacity = original.capacity;
            this.count = 0;
            this.random_generator = original.random_generator;
            this.registry = Hash_Map<>();
            Iterator<Buffer> iterator = original.iterator();
            while (iterator.has_next()) 
            {
                const Buffer current = iterator.next();
                // Create and register a buffer at the same level
                const Buffer new_copy = create(current.get_level());
                // Consume the data
                const std::vector<double> data = current.get_data();
                for (double value : data) 
                {
                    new_copy.consume(value);
                }
            }
        }

        /**
         * Tries to create a buffer with the given level.
         * <p>
         * If there is capacity to create a buffer (i.e., fewer than
         * count have been created), a buffer is created with the given
         * level, registered and returned.  If capacity has been reached, * NULL is returned.
         *
         * @param level level of the buffer
         * @return an empty buffer or NULL if a buffer can't be provided
         */
        Buffer create(const int& level) 
        {
            if (!can_create()) 
            {
                return NULL;
            }
            my_count++;
            auto buffer = Buffer(buffer_size, level, random_generator);
            std::vector<Buffer> buffer_list = my_registry.at(level);
            if (buffer_list.empty())
            {
                buffer_list = Array_list<>();
                registry.put(level, buffer_list);
            }
            buffer_list.push_back(buffer);
            if (level > my_max_level) 
            {
                my_max_level = level;
            }
            return buffer;
        }

        /**
         * Returns true if there is capacity to create a buffer.
         *
         * @return true if fewer than capacity buffers have been created.
         */
        bool can_create() const
        {
            return count < my_capacity;
        }

        /**
         * Returns true if we have used less than half of the allocated storage.
         * <p>
         * Includes a check to make sure all buffers have level 0;
         * but this should always be the case.
         * <p>
         * When this method returns true, we have all consumed data in storage
         * and enough space to make a copy of the combined dataset.
         *
         * @return true if all buffers have level 0 and less than half of the
         * available storage has been used
         */
        bool half_empty() 
        {
            return my_count * 2 < my_capacity && my_registry.size() == 1 && my_registry.contains_key(0);
        }

        /**
         * Returns a fresh copy of all data from level 0 buffers.
         *
         * @return combined data stored in all level 0 buffers
         */
        std::vector<double> level_zero_data() 
        {
            auto level_zero_buffers = my_registry.get(0);
            // First determine the combined size of the data
            int length = 0;
            for (const auto& buffer : level_zero_buffers) 
            {
                if (!buffer.has_capacity())
                { // full buffer
                    length += buffer.size;
                }
                else 
                {
                    length += buffer.next;  // filled amount
                }
            }
            // Copy the data
            int pos{};
            int curr_len;
            auto out = std::vector<double>(length);
            for (const auto& buffer : level_zero_buffers) 
            {
                curr_len = !buffer.has_capacity()
                    ? buffer.size
                    : buffer.next;
                System.arraycopy(buffer.data, 0, out, pos, curr_len);
                pos += curr_len;
            }
            return out;
        }

        /**
         * Finds the lowest level l where there exist at least two buffers, * merges them to create a buffer with level l+1 and returns
         * a free buffer with level l+1.
         *
         * @return free buffer that can accept data
         */
        Buffer merge() 
        {
            int l{};
            std::vector<Buffer> merge_candidates;
            // Find the lowest level containing at least two buffers
            while (merge_candidates.empty() && l <= my_max_level)
            {
                auto buffer_list = my_registry.at(l);
                if (buffer_list.empty() && buffer_list.size() > 1)
                {
                    merge_candidates = buffer_list;
                }
                else 
                {
                    l++;
                }
            }
            if (merge_candidates == NULL) 
            {
                throw std::exception("not implemented");
                // Should never happen
                //throw Math_Illegal_State_Exception(hipparchus::exception::Localized_Core_Formats_Type::INTERNAL_ERROR);
            }
            Buffer buffer1 = merge_candidates.get(0);
            Buffer buffer2 = merge_candidates.get(1);
            // Remove buffers to be merged
            merge_candidates.remove(0);
            merge_candidates.remove(0);
            // If these are the last level-l buffers, remove the empty list
            if (registry.get(l).size() == 0) 
            {
                registry.remove(l);
            }
            // Merge the buffers
            buffer1.merge_with(buffer2);
            // Now both buffers have level l+1; buffer1 is full and buffer2 is free.
            // Register both buffers
            register(buffer1);
            register(buffer2);

            // Return the free one
            return buffer2;
        }

        /**
         * Clears the buffer map.
         */
        void clear() 
        {
            for (std::vector<Buffer> buffer_list : my_registry.values()) 
            {
                buffer_list.clear();
            }
            my_registry.clear();
            my_count = 0;
        }

        /**
         * Registers a buffer.
         *
         * @param buffer Buffer to be registered.
         */
        void register(const Buffer& buffer) 
        {
            const int level = buffer.get_level();
            std::vector<Buffer> list = my_registry.get(level);
            if (list.empty()) 
            {
                list = Array_list<>();
                registry.put(level, list);
                if (level > max_level) 
                {
                    max_level = level;
                }
            }
            list.add(buffer);
        }

        /**
         * De-register a buffer.
         *
         * @param buffer Buffer to be de-registered
         * @Illegal_State_Exception if the buffer is not registered
         */
        void de_register(const Buffer& buffer) 
        {
            std::vector<Buffer> buffer_list = registry.get(buffer.get_level());
            const UUID target_id = buffer.get_id();
            int i{};
            bool found = false;
            while (i < buffer_list.size() && !found) 
            {
                if (buffer_list.get(i).get_id().equals(target_id)) 
                {
                    buffer_list.remove(i);
                    found = true;
                    buffer.clear();
                }
            }
            if (!found) 
            {
                throw std::exception("not implemented");
                //throw Math_Illegal_State_Exception(hipparchus::exception::Localized_Core_Formats_Type::INTERNAL_ERROR);
            }
        }

        /**
         * Returns an iterator over all of the buffers. Iteration goes by level, with
         * level 0 first.  Assumes there are no empty buffer lists.
         */
        //override
        Iterator<Buffer> iterator()
        {
            return Iterator<Buffer>()
            {
            private:
                /** Outer loop iterator, from level to level. */
                const auto level_iterator = registry.key_set().iterator();

                /** List of buffers at current level. */
                List<Buffer> current_list = registry.get(level_iterator.next());

                /** Inner loop iterator, from buffer to buffer. */
                Iterator<Buffer> buffer_iterator = current_list == NULL
                    ? NULL
                    : current_list.iterator();

            public:
                //override
                bool has_next()
                {
                    if (buffer_iterator == NULL)
                    {
                        return false;
                    }
                    if (buffer_iterator.has_next())
                    {
                        return true;
                    }
                    // The current level iterator has just finished, try to bump level
                    if (level_iterator.has_next())
                    {
                        List<Buffer> current_list = registry.get(level_iterator.next());
                        buffer_iterator = current_list.iterator();
                        return true;
                    }
                    else
                    {
                        // Nothing left, signal this by NULLing buffer_iterator
                        buffer_iterator = NULL;
                        return false;
                    }
                }

                //override
                Buffer next()
                {
                    if (has_next())
                    {
                        return buffer_iterator.next();
                    }
                    throw std::exception("not implemented");
                    //throw No_Such_Element_Exception();
                }

                //override
                void remove()
                {
                    throw std::exception("not implemented");
                    //throw Unsupported_Operation_Exception();
                }
            };
        };

        /**
         * Absorbs the data in other into this, merging buffers as necessary to trim
         * the aggregate down to capacity. This method is only used when aggregating
         * Random_Percentile instances.
         *
         * @param other other Buffer_Map to merge in
         */
        void absorb(Buffer_Map other) 
        {
            // Add all of other's buffers to the map - possibly exceeding cap
            int full_count{};
            Buffer not_full = NULL;
            Iterator<Buffer> other_iterator = other.iterator();
            while (other_iterator.has_next()) 
            {
                Buffer buffer = other_iterator.next();
                if (buffer.has_capacity()) 
                {
                    not_full = buffer;
                }
                else 
                {
                    full_count++;
                }
                register(buffer);
                my_count++;
            }
            // Determine how many extra buffers we now have: + old - capacity
            const int excess = full_count + (not_full == NULL ? 0 : 1) + my_count - capacity;
            // Now eliminate the excess by merging
            for (int i{}; i < excess - 1; i++) 
            {
                merge_up();
                my_count--;
            }
        }

        /**
         * Find two buffers, first and second, of minimal level. Then merge
         * first into second and discard first.
         * <p>
         * If the buffers have different levels, make second the higher level
         * buffer and make probability of selection in the merge proportional
         * to level weight ratio.
         * <p>
         * This method is only used when aggregating Random_Percentile instances.
         */
        void merge_up() 
        {
            // Find two minimum-level buffers to merge
            // Loop depends on two invariants:
            //   0) iterator goes in level order
            //   1) there are no empty lists in the registry
            Iterator<Buffer> buffer_iterator = iterator();
            Buffer first = NULL;
            Buffer second = NULL;
            while ((first == NULL || second == NULL) && buffer_iterator.has_next()) 
            {
                Buffer buffer = buffer_iterator.next();
                if (!buffer.has_capacity()) { // Skip not full buffers
                    if (first == NULL) 
                    {
                        first = buffer;
                    }
                    else 
                    {
                        second = buffer;
                    }
                }
            }
            if (first == NULL || second == NULL || first.level > second.level) 
            {
                throw std::exception("not implemented");
                //throw Math_Illegal_State_Exception(hipparchus::exception::Localized_Core_Formats_Type::INTERNAL_ERROR);
            }
            // Merge first into second and deregister first.
            // Assumes that first has level <= second (checked above).
            if (first.get_level() == second.get_level()) 
            {
                second.merge_with(first);
            }
            else 
            {
                first.merge_into(second);
            }
            de_register(first);
        }
    }

    /**
     * Computes the given percentile by combining the data from the collection
     * of aggregates. The result describes the combined sample of all data added
     * to any of the aggregates.
     *
     * @param percentile desired percentile (scaled 0-100)
     * @param aggregates Random_Percentile instances to combine data from
     * @return estimate of the given percentile using combined data from the aggregates
     * @ if percentile is out of the range [0, 100]
     */
    double reduce(double percentile, Collection<Random_Percentile> aggregates) 
    {
        if (percentile > 100 || percentile < 0) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::OUT_OF_RANGE, percentile, 0, 100);
        }

        // First see if we can copy all data and just compute exactly.
        // The following could be improved to verify that all have only level 0 buffers
        // and the sum of the data sizes is less than 1/2 total capacity.  Here we
        // just check that each of the aggregates is less than half full.
        Iterator<Random_Percentile> iterator = aggregates.iterator();
        bool small = true;
        while (small && iterator.has_next()) 
        {
            small = iterator.next().buffer_map.half_empty();
        }
        if (small) 
        {
            iterator = aggregates.iterator();
            std::vector<double> combined = {};
            while (iterator.has_next()) 
            {
               combined = Math_Arrays::concatenate(combined, iterator.next().buffer_map.level_zero_data());
            }
            const Percentile exact_p = Percentile(percentile);
            return exact_p.evaluate(combined);
        }

        // Below largely duplicates code in get_result(percentile).
        // Common binary search code with function parameter should be factored out.

        // Get global max and min to bound binary search and total N
        double min = INFINITY;
        double max = -INFINITY;
        double combined_n = 0;
        iterator = aggregates.iterator();
        while (iterator.has_next()) 
        {
            const Random_Percentile curr = iterator.next();
            const double cur_min = curr.get_result(0);
            const double cur_max = curr.get_result(100);
            if (cur_min < min) 
            {
                min = cur_min;
            }
            if (cur_max > max) 
            {
                max = cur_max;
            }
            combined_n += curr.get_n();
        }

        const double q = percentile / 100;
        // Handle degenerate cases
        if (Double.compare(q, 0.0) == 0) 
        {
            return min;
        }
        if (Double.compare(q, 1) == 0) 
        {
            return max;
        }

        // Compute target rank
        const double target_rank = q * combined_n;

        // Perform binary search using aggregated rank computation
        // Start with initial guess min + quantile * (max - min).
        double estimate = min + q * (max - min);
        double estimate_rank = get_aggregate_rank(estimate, aggregates);
        double lower;
        double upper;
        if (estimate_rank > target_rank) 
        {
            upper = estimate;
            lower = min;
        }
        else if (estimate_rank < target_rank) 
        {
            lower = estimate;
            upper = max;
        }
        else 
        {
            return estimate;
        }
        const double eps = epsilon / 2;
        double interval_width = std::abs(upper - lower);
        while (std::abs(estimate_rank / combined_n - q) > eps && interval_width > eps / combined_n) 
        {
            if (estimate_rank == target_rank) 
            {
                return estimate;
            }
            if (estimate_rank > target_rank) 
            {
                upper = estimate;
            }
            else 
            {
                lower = estimate;
            }
            interval_width = std::abs(upper - lower);
            estimate = lower + interval_width / 2;
            estimate_rank = get_aggregate_rank(estimate, aggregates);
        }
        return estimate;
    }

    /**
     * Computes the estimated rank of value in the combined dataset of the aggregates.
     * Sums the values from {@link #get_rankstatic_cast<double>(}.
     *
     * @param value value whose rank is sought
     * @param aggregates collection to aggregate rank over
     * @return estimated number of elements in the combined dataset that are less than value
     */
    double get_aggregate_rank(const double& value, Collection<Random_Percentile> aggregates) 
    {
        double result{};
        const Iterator<Random_Percentile> iterator = aggregates.iterator();
        while (iterator.has_next()) 
        {
            result += iterator.next().get_rank(value);
        }
        return result;
    }

    /**
     * Returns the estimated quantile position of value in the combined dataset of the aggregates.
     * Specifically, what is returned is an estimate of \(|\{x \in X : x &lt; value\}| / |X|\)
     * where \(X\) is the set of values that have been consumed from all of the datastreams
     * feeding the aggregates.
     *
     * @param value value whose quantile rank is sought.
     * @param aggregates collection of Random_Percentile instances being combined
     * @return estimated proportion of combined sample values that are strictly less than {@code value}
     */
    double get_aggregate_quantile_rank(const double& value, Collection<Random_Percentile> aggregates) 
    {
        return get_aggregate_rank(value, aggregates) / get_aggregate_n(aggregates);
    }

    /**
     * Returns the total number of values that have been consumed by the aggregates.
     *
     * @param aggregates collection of Random_Percentile instances whose combined sample size is sought
     * @return total number of values that have been consumed by the aggregates
     */
    double get_aggregate_n(Collection<Random_Percentile> aggregates) 
    {
        double result{};
        const Iterator<Random_Percentile> iterator = aggregates.iterator();
        while (iterator.has_next()) 
        {
            result += iterator.next().get_n();
        }
        return result;
    }

    /**
     * Aggregates the provided instance into this instance.
     * <p>
     * Other must have the same buffer size as this. If the combined data size
     * exceeds the maximum storage configured for this instance, buffers are
     * merged to create capacity. If all that is needed is computation of
     * aggregate results, {@link #reduce(double, Collection)} is faster, * may be more accurate and does not require the buffer sizes to be the same.
     *
     * @param other the instance to aggregate into this instance
     * @Null_Argument_Exception if the input is NULL
     * @Illegal_Argument_Exception if other has different buffer size than this
     */
    //override
    void aggregate(Random_Percentile other)
    {
        if (other == NULL) 
        {
            throw std::exception("not implemented");
            //throw Null_Argument_Exception();
        }
        if (other.s != my_s) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::INTERNAL_ERROR);
        }
        my_buffer_map.absorb(other.buffer_map);
        my_n += other.n;
    }

    /**
     * Returns the maximum number of {@code double} values that a {@code Random_Percentile}
     * instance created with the given {@code epsilon} value will retain in memory.
     * <p>
     * If the number of values that have been consumed from the stream is less than 1/2
     * of this value, reported statistics are exact.
     *
     * @param epsilon bound on the relative quantile error (see class javadoc)
     * @return upper bound on the total number of primitive double values retained in memory
     * @ if epsilon is not in the interval (0,1)
     */
    static long max_values_retained(const double& epsilon) 
    {
        if (epsilon >= 1) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_LARGE_BOUND_EXCLUDED, epsilon, 1);
        }
        if (epsilon <= 0) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL_BOUND_EXCLUDED, epsilon, 0);
        }
        const long h = static_cast<long>( std::ceil(log2(1/epsilon));
        const long s = static_cast<long>( std::ceil(std::sqrt(log2(1/epsilon)) / epsilon);
        return (h+1) * s;
    }
};