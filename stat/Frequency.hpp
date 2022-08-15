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

#include <unordered_map>
#include "../core/util/Comparable.h"

//import java.text.Number_Format;
//import java.util.Collection;
//import java.util.Comparator;
//import java.util.Iterator;
//import java.util.List;
//import java.util.Map;
//import java.util.Navigable_Map;
//import java.util.Objects;
//import java.util.Tree_Map;
//import java.util.stream.Collectors;
#include <map>

//import org.hipparchus.exception.Null_Argument_Exception;
//import org.hipparchus.util.Math_Utils;

/**
 * Maintains a frequency distribution of Comparable values.
 * <p>
 * The values are ordered using the default (natural order), unless a
 * {@code Comparator} is supplied in the constructor.
 *
 * @see Long_Frequency
 * @param <T> the element type
 */
template<typename T, typename std::enable_if<std::is_base_of<Comparable, T>::value>::type* = nullptr>
class Frequency 
{
private:

    /** underlying collection */
    const std::map<T, long> my_freq_table;

public:
    /**
     * Default constructor.
     */
    Frequency() : my_freq_table{}
    {
        //my_freq_table = Tree_Map<>();
    }

    /**
     * Constructor allowing values Comparator to be specified.
     *
     * @param comparator Comparator used to order values
     */
    Frequency(Comparator<T> comparator) 
    {
        my_freq_table = Tree_Map<>(comparator);
    }

    std::map<T, long> get_freq_table() const
    {
        return my_freq_table;
    }

    bool operator==(const Comparable& other) override const
    {
        my_freq_table == other.get_freq_table();
    }

    /**
     * Adds 1 to the frequency count for v.
     *
     * @param v the value to add.
     */
    void add_value(T v) 
    {
        increment_value(v, 1);
    }

    /**
     * Increments the frequency count for v.
     *
     * @param v the value to add.
     * @param increment the amount by which the value should be incremented
     */
    void increment_value(const T& v, long increment) 
    {
        long count = my_freq_table.get_or_default(v, 0L);
        my_freq_table.put(v, static_cast<long>(count.long_value() + increment));
    }

    /** Clears the frequency table */
    void clear() 
    {
        my_freq_table.clear();
    }

    /**
     * Returns an Iterator over the set of values that have been added.
     *
     * @return values Iterator
     */
    Iterator<T> values_iterator() 
    {
        return my_freq_table.key_set().iterator();
    }

    /**
     * Return an Iterator over the set of keys and values that have been added.
     * Using the entry set to iterate is more efficient in the case where you
     * need to access respective counts as well as values, since it doesn't
     * require a "get" for every key...the value is provided in the Map.Entry.
     *
     * @return entry set Iterator
     */

    Iterator<Map.Entry<T, long>> entry_set_iterator() 
    {
        return my_freq_table.entry_set().iterator();
    }

    //-------------------------------------------------------------------------

    /**
     * Returns the sum of all frequencies.
     *
     * @return the total frequency count.
     */
    long get_sum_freq() 
    {
        return my_freq_table.values()
                        .stream()
                        .map_to_long(long::long_value)
                        .sum();
    }

    /**
     * Returns the number of values equal to v.
     * Returns 0 if the value is not comparable.
     *
     * @param v the value to lookup.
     * @return the frequency of v.
     */
    long get_count(T v) 
    {
        return my_freq_table.get_or_default(v, 0L);
    }

    /**
     * Returns the number of values in the frequency table.
     *
     * @return the number of unique values that have been added to the frequency table.
     * @see #values_iterator()
     */
    int get_unique_count()
    {
        return my_freq_table.key_set().size();
    }

    /**
     * Returns the percentage of values that are equal to v
     * (as a proportion between 0 and 1).
     * <p>
     * Returns {@codeNAN} if no values have been added.
     *
     * @param v the value to lookup
     * @return the proportion of values equal to v
     */
    double get_pct(T v) 
    {
        const long sum_freq = get_sum_freq();
        if (sum_freq == 0) 
        {
            return std::numeric_limits<double>::quiet_NaN();
        }
        return static_cast<double>( get_count(v) / static_cast<double>( sum_freq;
    }

    //-----------------------------------------------------------------------------------------

    /**
     * Returns the cumulative frequency of values less than or equal to v.
     *
     * @param v the value to lookup.
     * @return the proportion of values equal to v
     */
    long get_cum_freq(T v) 
    {
        if (get_sum_freq() == 0) 
        {
            return 0;
        }

        std::map<T, long> head_map = my_freq_table.head_map(v, true);

        if (head_map.is_empty()) 
        {
            // v is less than first value
            return 0;
        }
        if (head_map.size() == my_freq_table.size())
        {
            // v is greater than or equal to last value
            return get_sum_freq();
        }

        return head_map.values()
                      .stream()
                      .map_to_long(long::long_value)
                      .sum();
    }

    //----------------------------------------------------------------------------------------------

    /**
     * Returns the cumulative percentage of values less than or equal to v
     * (as a proportion between 0 and 1).
     * <p>
     * Returns {@codeNAN} if no values have been added.
     *
     * @param v the value to lookup
     * @return the proportion of values less than or equal to v
     */
    double get_cum_pct(T v) 
    {
        const long sum_freq = get_sum_freq();
        if (sum_freq == 0) 
        {
            return std::numeric_limits<double>::quiet_NaN();
        }
        return static_cast<double>( get_cum_freq(v) / static_cast<double>( sum_freq;
    }

    /**
     * Returns the mode value(s) in comparator order.
     *
     * @return a list containing the value(s) which appear most often.
     */
    std::vector<T> get_mode() 
    {
        // Get the max count first
        const long most_popular =
            my_freq_table.values()
                         .stream()
                         .map_to_long(long::long_value)
                         .max()
                         .or_else(0L);

        return my_freq_table.entry_set()
                        .stream()
                        .filter(entry -> entry.get_value() == most_popular)
                        .map(entry -> entry.get_key())
                        .collect(Collectors.to_list());
    }

    //----------------------------------------------------------------------------------------------

    /**
     * Merge another Frequency object's counts into this instance.
     * This Frequency's counts will be incremented (or set when not already set)
     * by the counts represented by other.
     *
     * @param other the other {@link Frequency} object to be merged
     * @Null_Argument_Exception if {@code other} is NULL
     */
    void merge(const Frequency<? extends T> other) 
    {
        //Math_Utils::check_not_null(other);

        Iterator<? extends Map.Entry<? extends T, long>> iter = other.entry_set_iterator();
        while (iter.has_next()) 
        {
            const Map.Entry<? extends T, long> entry = iter.next();
            increment_value(entry.get_key(), entry.get_value().long_value());
        }
    }

    /**
     * Merge a {@link Collection} of {@link Frequency} objects into this instance.
     * This Frequency's counts will be incremented (or set when not already set)
     * by the counts represented by each of the others.
     *
     * @param others the other {@link Frequency} objects to be merged
     * @Null_Argument_Exception if the collection is NULL
     */
    //template<typename T, typename std::enable_if<std::is_base_of<Frequency, T>::value>::type* = nullptr>
    //void merge(const std::vector<T extends Frequency< extends T>> others)
    void merge(const std::vector<int>& others)
    {
        throw std::exception("Not implemented");
        //Math_Utils::check_not_null(others);

        //for (const Frequency<T> freq : others) 
        //{
        //    merge(freq);
        //}
    }

    /**
     * Return a string representation of this frequency distribution.
     *
     * @return a string representation.
     */
    //override
    std::string to_string() const 
    {
        Number_Format nf = Number_Format.get_percent_instance();
        std::stringBuilder out_buffer = std::stringBuilder(200); // this size is just a wild guess
        out_buffer.append("Value 	Freq. 	Pct. 	Cum Pct. \n");
        Iterator<T> iter = my_freq_table.key_set().iterator();
        while (iter.has_next()) 
        {
            T value = iter.next();
            out_buffer.append(value).
                      append('	').
                      append(get_count(value)).
                      append('	').
                      append(nf.format(get_pct(value))).
                      append('	').
                      append(nf.format(get_cum_pct(value))).
                      append('\n');
        }
        return out_buffer.to_string();
    }

    /** {@inherit_doc} */
    //override
    int hash_code() 
    {
        constexpr int prime{ 31 };
        int result{ 1 };
        result = prime * result + ((my_freq_table == NULL)
            ? 0
            : my_freq_table.hash_code());
        return result;
    }

    /** {@inherit_doc} */
    //override
    bool equals(const T& obj) 
    {
        return *this == obj;
        ///*if (this == obj) 
        //{
        //    return true;
        //}
        //if (!(obj instanceof Frequency)) 
        //{
        //    return false;
        //}
        //Frequency<?> other = (Frequency<?>) obj;
        //return Objects.equals(my_freq_table, other.get_freq_table(*/));
    }
};