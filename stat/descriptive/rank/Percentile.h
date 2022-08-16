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
//package org.hipparchus.stat.descriptive.rank;

//import java.io.Serializable;
//import java.util.Arrays;
//import java.util.Bit_Set;

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.;
//import org.hipparchus.stat.Localized_Stat_Formats;
//import org.hipparchus.stat.descriptive.Abstract_Univariate_Statistic;
//import org.hipparchus.stat.ranking.NaN_Strategy;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Kth_Selector;
//import org.hipparchus.util.Math_Arrays;
//import org.hipparchus.util.Math_Utils;
//import org.hipparchus.util.Pivoting_Strategy;
//import org.hipparchus.util.Precision;
#include <vector>
#include "../../../core/util/MathArrays.h"
#include "../../../core/util/MathUtils.h"
#include "../AbstractUnivariateStatistic.h"
#include "../../../core/util/KthSelector.h"

/**
 * Provides percentile computation.
 * <p>
 * There are several commonly used methods for estimating percentiles (a.k.a.
 * quantiles) based on sample data.  For large samples, the different methods
 * agree closely, but when sample sizes are small, different methods will give
 * significantly different results.  The algorithm implemented here works as follows:
 * <ol>
 * <li>Let <code>n</code> be the length of the (sorted) array and
 * <code>0 &lt; p &lt;= 100</code> be the desired percentile.</li>
 * <li>If <code> n = 1 </code> return the unique array element (regardless of
 * the value of <code>p</code>); otherwise </li>
 * <li>Compute the estimated percentile position
 * <code> pos = p * (n + 1) / 100</code> and the difference, <code>d</code>
 * between <code>pos</code> and <code>floor(pos)</code> (i.e. the fractional
 * part of <code>pos</code>).</li>
 * <li> If <code>pos &lt; 1</code> return the smallest element in the array.</li>
 * <li> Else if <code>pos &gt;= n</code> return the largest element in the array.</li>
 * <li> Else let <code>lower</code> be the element in position
 * <code>floor(pos)</code> in the array and let <code>upper</code> be the
 * next element in the array.  Return <code>lower + d * (upper - lower)</code>
 * </li>
 * </ol>
 * <p>
 * To compute percentiles, the data must be at least partially ordered.  Input
 * arrays are copied and recursively partitioned using an ordering definition.
 * The ordering used by <code>Arrays.sort(std::vector<double>)</code> is the one determined
 * by {@link java.lang.Double#compare_tostatic_cast<double>(}.  This ordering makes
 * <code>Double.NaN</code> larger than any other value (including
 * <code>INFINITY</code>).  Therefore, for example, the median
 * (50th percentile) of
 * <code>{0, 1, 2, 3, 4,NAN}</code> evaluates to <code>2.5.</code>
 * <p>
 * sin_ce percentile estimation usually involves interpolation between array
 * elements, arrays containing  <code>NaN</code> or infinite values will often
 * result in <code>NaN</code> or infinite values returned.
 * <p>
 * Further, to include different estimation types such as R1, R2 as mentioned in
 * <a href="http://en.wikipedia.org/wiki/Quantile">Quantile page(wikipedia)</a>, * a type specific NaN handling strategy is used to closely match with the
 * typically observed results from popular tools like R(R1-R9), Excel(R7).
 * <p>
 * Percentile uses only selection instead of complete sorting and caches selection
 * algorithm state between calls to the various {@code evaluate} methods. This
 * greatly improves efficiency, both for a single percentile and multiple percentile
 * computations. To maximize performance when multiple percentiles are computed
 * based on the same data, users should set the data array once using either one
 * of the {@link #evaluate(std::vector<double>, double)} or {@link #set_data(std::vector<double>)} methods
 * and thereafter {@link #evaluatestatic_cast<double>(} with just the percentile provided.
 * <p>
 * <strong>Note that this implementation is not synchronized.</strong> If
 * multiple threads access an instance of this class concurrently, and at least
 * one of the threads invokes the <code>increment()</code> or
 * <code>clear()</code> method, it must be synchronized externally.
 */
class Percentile : public Abstract_Univariate_Statistic  
{
    /** Maximum number of partitioning pivots cached (each level double the number of pivots). */
    static constexpr int MAX_CACHED_LEVELS{ 10 };

    /** Maximum number of cached pivots in the pivots cached array */
    static constexpr int PIVOTS_HEAP_LENGTH{ 0x1 << MAX_CACHED_LEVELS - 1 };

    /** Default Kth_Selector used with default pivoting strategy */
    const Kth_Selector my_kth_selector;

    /** Any of the {@link Estimation_Type}s such as {@link Estimation_Type#LEGACY CM} can be used. */
    const Estimation_Type my_estimation_type;

    /** NaN Handling of the input as defined by {@link NaN_Strategy} */
    const NaN_Strategy my_nan_strategy;

    /**
     * Determines what percentile is computed when evaluate() is activated
     * with no quantile argument.
     */
    double quantile;

    /** Cached pivots. */
    std::vector<int> cached_pivots;

protected:
    /**
     * Constructs a Percentile with the specific quantile value, * {@link Estimation_Type}, {@link NaN_Strategy} and {@link Kth_Selector}.
     *
     * @param quantile the quantile to be computed
     * @param estimation_type one of the percentile {@link Estimation_Type  estimation types}
     * @param nan_strategy one of {@link NaN_Strategy} to handle with NaNs
     * @param kth_selector a {@link Kth_Selector} to use for pivoting during search
     * @ if p is not within (0,100]
     * @ if type or NaN_Strategy passed is NULL
     */
    Percentile(const double& quantile, const Estimation_Type& estimation_type, const NaN_Strategy& nan_strategy, const Kth_Selector& kth_selector)
    {
        set_quantile(quantile);
        cached_pivots = NULL;
        //Math_Utils::check_not_null(estimation_type);
        //Math_Utils::check_not_null(nan_strategy);
        //Math_Utils::check_not_null(kth_selector);
        my_estimation_type = estimation_type;
        my_nan_strategy = nan_strategy;
        my_kth_selector = kth_selector;
    }

    /**
     * Get the work array to operate. Makes use of prior {@code stored_data} if
     * it exists or else do a check on NaNs and copy a subset of the array
     * defined by begin and length parameters. The set {@link #nan_strategy} will
     * be used to either retain/remove/replace any NaNs present before returning
     * the resultant array.
     *
     * @param values the array of numbers
     * @param begin index to start reading the array
     * @param length the length of array to be read from the begin index
     * @return work array sliced from values in the range [begin,begin+length)
     * @ if values or indices are invalid
     */
    std::vector<double> get_work_array(const std::vector<double>& values, const int& begin, const int& length)
    {
        const std::vector<double> work;
        if (values == get_data_ref())
        {
            work = get_data_ref();
        }
        else
        {
            switch (nan_strategy)
            {
            case MAXIMAL:// Replace NaNs with +INFs
                work = replace_and_slice(values, begin, length, NAN, INFINITY);
                break;
            case MINIMAL:// Replace NaNs with -INFs
                work = replace_and_slice(values, begin, length,NAN, -INFINITY);
                break;
            case REMOVED:// Drop NaNs from data
                work = remove_and_slice(values, begin, length,NAN);
                break;
            case FAILED:// just throw exception as NaN is un-acceptable
                work = copy_of(values, begin, length);
                Math_Arrays::check_not_na_n(work);
                break;
            default: //FIXED
                work = copy_of(values, begin, length);
                break;
            }
        }
        return work;
    }

public:
    /**
     * Constructs a Percentile with the following defaults.
     * <ul>
     *   <li>default quantile: 50.0, can be reset with {@link #set_quantilestatic_cast<double>(}</li>
     *   <li>default estimation type: {@link Estimation_Type#LEGACY}, *   can be reset with {@link #with_estimation_type(Estimation_Type)}</li>
     *   <li>default NaN strategy: {@link NaN_Strategy#REMOVED}, *   can be reset with {@link #with_na_n_strategy(NaN_Strategy)}</li>
     *   <li>a Kth_Selector that makes use of {@link Pivoting_Strategy#MEDIAN_OF_3}, *   can be reset with {@link #with_kth_selector(Kth_Selector)}</li>
     * </ul>
     */
    Percentile() 
    {
        // No try-catch or advertised exception here - arg is valid
        Percentile(50.0);
    }

    /**
     * Constructs a Percentile with the specific quantile value and the following
     * <ul>
     *   <li>default method type: {@link Estimation_Type#LEGACY}</li>
     *   <li>default NaN strategy: {@link NaN_Strategy#REMOVED}</li>
     *   <li>a Kth Selector : {@link Kth_Selector}</li>
     * </ul>
     * @param quantile the quantile
     * @  if p is not greater than 0 and less
     * than or equal to 100
     */
    Percentile(const double& quantile)  
    {
        Percentile(quantile, Estimation_Type.LEGACY, NaN_Strategy.REMOVED, Kth_Selector(Pivoting_Strategy.MEDIAN_OF_3));
    }

    /**
     * Copy constructor, creates a {@code Percentile} identical
     * to the {@code original}
     *
     * @param original the {@code Percentile} instance to copy
     * @ if original is NULL
     */
    Percentile(const Percentile& original)
    {
        super(original);
        my_estimation_type   = original.get_estimation_type();
        my_nan_strategy      = original.get_nan_strategy();
        my_kth_selector      = original.get_kth_selector();

        set_data(original.get_data_ref());
 
        cached_pivots = original.cached_pivots;
        set_quantile(original.quantile);
    }



    /** {@inherit_doc} */
    //override
    void set_data(const std::vector<double>& values) 
    {
        cached_pivots = std::vector<int>(PIVOTS_HEAP_LENGTH, -1);
        super.set_data(values);
    }

    /** {@inherit_doc} */
    //override
    void set_data(const std::vector<double>& values, const int& begin, const int& length)
    {
        //Math_Utils::check_not_null(values, hipparchus::exception::Localized_Core_Formats_Type::INPUT_ARRAY);
        cached_pivots = std::vector<int>(PIVOTS_HEAP_LENGTH, -1);
        super.set_data(values, begin, length);
    }

    /**
     * Returns the result of evaluating the statistic over the stored data.
     * <p>
     * The stored array is the one which was set by previous calls to
     * {@link #set_data(std::vector<double>)}
     *
     * @param p the percentile value to compute
     * @return the value of the statistic applied to the stored data
     * @ if p is not a valid quantile value
     * (p must be greater than 0 and less than or equal to 100)
     */
    double evaluate(const double& p)  
    {
        return evaluate(get_data_ref(), p);
    }

    /**
     * Returns an estimate of the <code>quantile</code>th percentile of the
     * designated values in the <code>values</code> array.  The quantile
     * estimated is determined by the <code>quantile</code> property.
     * <p>
     * <ul>
     * <li>Returns <code>Double.NaN</code> if <code>length = 0</code></li>
     * <li>Returns (for any value of <code>quantile</code>)
     * <code>values[begin]</code> if <code>length = 1 </code></li>
     * <li>Throws <code></code> if <code>values</code>
     * is NULL, or <code>start</code> or <code>length</code> is invalid</li>
     * </ul>
     * <p>
     * See {@link Percentile} for a description of the percentile estimation
     * algorithm used.
     *
     * @param values the input array
     * @param start index of the first array element to include
     * @param length the number of elements to include
     * @return the percentile value
     * @ if the parameters are not valid
     *
     */
    //override
    double evaluate(const std::vector<double>& values, const int& start, const int& length)
    {
        return evaluate(values, start, length, quantile);
    }

    /**
     * Returns an estimate of the <code>p</code>th percentile of the values
     * in the <code>values</code> array.
     * <p>
     * <ul>
     * <li>Returns <code>Double.NaN</code> if <code>values</code> has length
     * <code>0</code></li>
     * <li>Returns (for any value of <code>p</code>) <code>values[0]</code>
     *  if <code>values</code> has length <code>1</code></li>
     * <li>Throws <code></code> if <code>values</code>
     * is NULL or p is not a valid quantile value (p must be greater than 0
     * and less than or equal to 100) </li>
     * </ul>
     * <p>
     * The default implementation delegates to
     * <code>evaluate(std::vector<double>, int, int, double)</code> in the natural way.
     *
     * @param values input array of values
     * @param p the percentile value to compute
     * @return the percentile value orNAN if the array is empty
     * @ if <code>values</code> is NULL or p is invalid
     */
    double evaluate(const std::vector<double>& values, const double& p)
    {
        //Math_Utils::check_not_null(values, hipparchus::exception::Localized_Core_Formats_Type::INPUT_ARRAY);
        return evaluate(values, 0, values.size(), p);
    }

    /**
     * Returns an estimate of the <code>p</code>th percentile of the values
     * in the <code>values</code> array, starting with the element in (0-based)
     * position <code>begin</code> in the array and including <code>length</code>
     * values.
     * <p>
     * Calls to this method do not modify the internal <code>quantile</code>
     * state of this statistic.
     * <p>
     * <ul>
     * <li>Returns <code>Double.NaN</code> if <code>length = 0</code></li>
     * <li>Returns (for any value of <code>p</code>) <code>values[begin]</code>
     *  if <code>length = 1 </code></li>
     * <li>Throws <code></code> if <code>values</code>
     *  is NULL , <code>begin</code> or <code>length</code> is invalid, or
     * <code>p</code> is not a valid quantile value (p must be greater than 0
     * and less than or equal to 100)</li>
     * </ul>
     * <p>
     * See {@link Percentile} for a description of the percentile estimation
     * algorithm used.
     *
     * @param values array of input values
     * @param p  the percentile to compute
     * @param begin  the first (0-based) element to include in the computation
     * @param length  the number of array elements to include
     * @return  the percentile value
     * @ if the parameters are not valid or the
     * input array is NULL
     */
    double evaluate(const std::vector<double>& values, const int& begin, const int& length, const double& p)     
    {

        Math_Arrays::verify_values(values, begin, length);
        if (p > 100 || p <= 0) 
        {
            throw std::exception("error not implemented, p > 100 || p <= 0");
            //throw (Localized_Stat_Formats.OUT_OF_BOUNDS_QUANTILE_VALUE, p, 0, 100);
        }
        if (length == 0) 
        {
            return std::numeric_limits<double>::quiet_NaN();
        }
        if (length == 1) 
        {
            return values[begin]; // always return single value for n = 1
        }

        const auto work = get_work_array(values, begin, length);
        const auto pivots_heap = get_pivots(values);
        return work.size() == 0
            ? NAN
            : estimation_type.evaluate(work, pivots_heap, p, kth_selector);
    }

    /**
     * Returns the value of the quantile field (determines what percentile is
     * computed when evaluate() is called with no quantile argument).
     *
     * @return quantile set while construction or {@link #set_quantilestatic_cast<double>(}
     */
    double get_quantile() 
    {
        return quantile;
    }

    /**
     * Sets the value of the quantile field (determines what percentile is
     * computed when evaluate() is called with no quantile argument).
     *
     * @param p a value between 0 &lt; p &lt;= 100
     * @  if p is not greater than 0 and less
     * than or equal to 100
     */
    void set_quantile(const double& p)  
    {
        if (p <= 0 || p > 100) 
        {
            throw std::exception("error not implemented, (p <= 0 || p > 100) ");
            //throw (Localized_Stat_Formats.OUT_OF_BOUNDS_QUANTILE_VALUE, p, 0, 100);
        }
        quantile = p;
    }

    /** {@inherit_doc} */
    //override
    Percentile copy() 
    {
        return Percentile(*this);
    }

    

    /**
     * Make a copy of the array for the slice defined by array part from
     * [begin, begin+length)
     * @param values the input array
     * @param begin start index of the array to include
     * @param length number of elements to include from begin
     * @return copy of a slice of the original array
     */
    private static std::vector<double> copy_of(const std::vector<double>& values, const int& begin, const int& length) 
    {
        Math_Arrays::verify_values(values, begin, length);
        return Arrays.copy_of_range(values, begin, begin + length);
    }

    /**
     * Replace every occurrence of a given value with a replacement value in a
     * copied slice of array defined by array part from [begin, begin+length).
     * @param values the input array
     * @param begin start index of the array to include
     * @param length number of elements to include from begin
     * @param original the value to be replaced with
     * @param replacement the value to be used for replacement
     * @return the copy of sliced array with replaced values
     */
    private static std::vector<double> replace_and_slice(const std::vector<double>& values, const int& begin, const int& length, const double& original, const double& replacement) 
    {
        auto temp = copy_of(values, begin, length);
        for (const int i{}; i < length; i++)
        {
            temp[i] = Precision::equals_including_nan(original, temp[i]) ?
                      replacement : temp[i];
        }
        return temp;
    }

    /**
     * Remove the occurrence of a given value in a copied slice of array
     * defined by the array part from [begin, begin+length).
     * @param values the input array
     * @param begin start index of the array to include
     * @param length number of elements to include from begin
     * @param removed_value the value to be removed from the sliced array
     * @return the copy of the sliced array after removing the removed_value
     */
    private static std::vector<double> remove_and_slice(const std::vector<double>& values, const int& begin, const int& length, const double removed_value) 
    {
        Math_Arrays::verify_values(values, begin, length);
        const std::vector<double> temp;
        //Bit_Set(length) to indicate where the removed_value is located
        const Bit_Set bits = Bit_Set(length);
        for (int i{ begin }; i < begin+length; i++) 
        {
            if (Precision::equals_including_nan(removed_value, values[i])) 
            {
                bits.set(i - begin);
            }
        }
        //Check if empty then create a copy
        if (bits.is_empty()) 
        {
            temp = copy_of(values, begin, length); // Nothing removed, just copy
        }
        else if (bits.cardinality() == length) 
        {
            temp = std::vector<double>(0);                 // All removed, just empty
        }
        else {                                   // Some removable, so new
            temp = std::vector<double>(length - bits.cardinality());
            int start = begin;  //start index from source array (i.e values)
            int dest = 0;       //dest index in destination array(i.e temp)
            int bit_set_ptr = 0;  //bit_set_ptr is start index pointer of bitset
            for (const int& next_one = bits.next_set_bit(bit_set_ptr); next_one != -1; next_one = bits.next_set_bit(bit_set_ptr)) 
            {
                const int length_to_copy = next_one - bit_set_ptr;
                System.arraycopy(values, start, temp, dest, length_to_copy);
                dest += length_to_copy;
                start = begin + (bit_set_ptr = bits.next_clear_bit(next_one));
            }
            //Copy any residue past start index till begin+length
            if (start < begin + length) 
            {
                System.arraycopy(values,start,temp,dest,begin + length - start);
            }
        }
        return temp;
    }

    /**
     * Get pivots which is either cached or a newly created one
     *
     * @param values array containing the input numbers
     * @return cached pivots or a newly created one
     */
    private std::vector<int> get_pivots(const std::vector<double>& values) 
    {
        const std::vector<int> pivots_heap;
        if (values == get_data_ref()) 
        {
            pivots_heap = cached_pivots;
        }
        else 
        {
            pivots_heap = int[PIVOTS_HEAP_LENGTH];
            Arrays.fill(pivots_heap, -1);
        }
        return pivots_heap;
    }

    /**
     * Get the estimation {@link Estimation_Type type} used for computation.
     *
     * @return the {@code estimation_type} set
     */
    public Estimation_Type get_estimation_type() 
    {
        return estimation_type;
    }

    /**
     * Build a instance similar to the current one except for the
     * {@link Estimation_Type estimation type}.
     * <p>
     * This method is intended to be used as part of a fluent-type builder
     * pattern. Building finely tune instances should be done as follows:
     * <pre>
     *   Percentile customized = Percentile(quantile).
     *                           with_estimation_type(estimation_type).
     *                           with_na_n_strategy(nan_strategy).
     *                           with_kth_selector(kth_selector);
     * </pre>
     * <p>
     * If any of the {@code with_xxx} method is omitted, the default value for
     * the corresponding customization parameter will be used.
     *
     * @param new_estimation_type estimation type for the instance
     * @return a instance, with changed estimation type
     * @ when new_estimation_type is NULL
     */
    public Percentile with_estimation_type(const Estimation_Type new_estimation_type) 
    {
        return Percentile(quantile, new_estimation_type, nan_strategy, kth_selector);
    }

    /**
     * Get the {@link NaN_Strategy NaN Handling} strategy used for computation.
     * @return {@code NaN Handling} strategy set during construction
     */
    public NaN_Strategy get_nan_strategy() 
    {
        return nan_strategy;
    }

    /**
     * Build a instance similar to the current one except for the
     * {@link NaN_Strategy NaN handling} strategy.
     * <p>
     * This method is intended to be used as part of a fluent-type builder
     * pattern. Building finely tune instances should be done as follows:
     * <pre>
     *   Percentile customized = Percentile(quantile).
     *                           with_estimation_type(estimation_type).
     *                           with_na_n_strategy(nan_strategy).
     *                           with_kth_selector(kth_selector);
     * </pre>
     * <p>
     * If any of the {@code with_xxx} method is omitted, the default value for
     * the corresponding customization parameter will be used.
     *
     * @param new_nan_strategy NaN strategy for the instance
     * @return a instance, with changed NaN handling strategy
     * @ when new_nan_strategy is NULL
     */
    public Percentile with_na_n_strategy(const NaN_Strategy new_nan_strategy) 
    {
        return Percentile(quantile, estimation_type, new_nan_strategy, kth_selector);
    }

    /**
     * Get the {@link Kth_Selector kth_selector} used for computation.
     * @return the {@code kth_selector} set
     */
    public Kth_Selector get_kth_selector() 
    {
        return kth_selector;
    }

    /**
     * Get the {@link Pivoting_Strategy} used in Kth_Selector for computation.
     * @return the pivoting strategy set
     */
    public Pivoting_Strategy get_pivoting_strategy() 
    {
        return kth_selector.get_pivoting_strategy();
    }

    /**
     * Build a instance similar to the current one except for the
     * {@link Kth_Selector kth_selector} instance specifically set.
     * <p>
     * This method is intended to be used as part of a fluent-type builder
     * pattern. Building finely tune instances should be done as follows:
     * <pre>
     *   Percentile customized = Percentile(quantile).
     *                           with_estimation_type(estimation_type).
     *                           with_na_n_strategy(nan_strategy).
     *                           with_kth_selector(new_kth_selector);
     * </pre>
     * <p>
     * If any of the {@code with_xxx} method is omitted, the default value for
     * the corresponding customization parameter will be used.
     *
     * @param new_kth_selector Kth_Selector for the instance
     * @return a instance, with changed Kth_Selector
     * @ when new_kth_selector is NULL
     */
    public Percentile with_kth_selector(const Kth_Selector new_kth_selector) 
    {
        return Percentile(quantile, estimation_type, nan_strategy, new_kth_selector);
    }

    /**
     * An enum for various estimation strategies of a percentile referred in
     * <a href="http://en.wikipedia.org/wiki/Quantile">wikipedia on quantile</a>
     * with the names of enum matching those of types mentioned in
     * wikipedia.
     * <p>
     * Each enum corresponding to the specific type of estimation in wikipedia
     * :  the respective formulae that specializes in the below aspects
     * <ul>
     * <li>An <b>index method</b> to calculate approximate index of the
     * estimate</li>
     * <li>An <b>estimate method</b> to estimate a value found at the earlier
     * computed index</li>
     * <li>A <b> min_limit</b> on the quantile for which first element of sorted
     * input is returned as an estimate </li>
     * <li>A <b> max_limit</b> on the quantile for which last element of sorted
     * input is returned as an estimate </li>
     * </ul>
     * <p>
     * Users can now create {@link Percentile} by explicitly passing this enum;
     * such as by invoking {@link Percentile#with_estimation_type(Estimation_Type)}
     * <p>
     * References:
     * <ol>
     * <li>
     * <a href="http://en.wikipedia.org/wiki/Quantile">Wikipedia on quantile</a>
     * </li>
     * <li>
     * <a href="https://www.amherst.edu/media/view/129116/.../Sample+Quantiles.pdf">
     * Hyndman, R. J. and Fan, Y. (1996) Sample quantiles in statistical
     * //packages, American Statistician 50, 361–365</a> </li>
     * <li>
     * <a href="http://stat.ethz.ch/R-manual/R-devel/library/stats/html/quantile.html">
     * R-Manual </a></li>
     * </ol>
     */
    enum Estimation_Type 
    {
        /**
         * This is the default type used in the {@link Percentile}.This method
         * has the following formulae for index and estimates<br>
         * \( \begin{align}
         * &amp;index    = (N+1)p\ \\
         * &amp;estimate = x_{\lceil h\,-\,1/2 \rceil} \\
         * &amp;min_limit = 0 \\
         * &amp;max_limit = 1 \\
         * \end{align}\)
         */
        LEGACY("Legacy Hipparchus") 
        {
            /**
             * {@inherit_doc}.This method in particular makes use of existing
             * Hipparchus style of picking up the index.
             */
            //override
            protected double index(const double p, const int length) 
            {
                const double min_limit = 0;
                const double max_limit = 1d;
                return Double.compare(p, min_limit) == 0 ? 0 :
                       Double.compare(p, max_limit) == 0 ?
                               length : p * (length + 1);
            }
        }, /**
         * The method R_1 has the following formulae for index and estimates<br>
         * \( \begin{align}
         * &amp;index= Np + 1/2\,  \\
         * &amp;estimate= x_{\lceil h\,-\,1/2 \rceil} \\
         * &amp;min_limit = 0 \\
         * \end{align}\)
         */
        R_1("R-1") 
        {

            //override
            protected double index(const double p, const int length) 
            {
                const double min_limit = 0;
                return Double.compare(p, min_limit) == 0 ? 0 : length * p + 0.5;
            }

            /**
             * {@inherit_doc}This method in particular for R_1 uses ceil(pos-0.5)
             */
            //override
            protected double estimate(const std::vector<double>& values, const std::vector<int> pivots_heap, const double pos, const int length, const Kth_Selector selector) 
            {
                return super.estimate(values, pivots_heap, std::ceil(pos - 0.5), length, selector);
            }

        }, /**
         * The method R_2 has the following formulae for index and estimates<br>
         * \( \begin{align}
         * &amp;index= Np + 1/2\, \\
         * &amp;estimate=\frac{x_{\lceil h\,-\,1/2 \rceil} +
         * x_{\lfloor h\,+\,1/2 \rfloor}}{2} \\
         * &amp;min_limit = 0 \\
         * &amp;max_limit = 1 \\
         * \end{align}\)
         */
        R_2("R-2") 
        {

            //override
            protected double index(const double p, const int length) 
            {
                const double min_limit = 0;
                const double max_limit = 1d;
                return Double.compare(p, max_limit) == 0 ? length :
                       Double.compare(p, min_limit) == 0 ? 0 : length * p + 0.5;
            }

            /**
             * {@inherit_doc}This method in particular for R_2 averages the
             * values at ceil(p+0.5) and floor(p-0.5).
             */
            //override
            protected double estimate(const std::vector<double>& values, const std::vector<int> pivots_heap, const double pos, const int length, const Kth_Selector selector) 
            {
                const double low =
                        super.estimate(values, pivots_heap, std::ceil(pos - 0.5), length, selector);
                const double high =
                        super.estimate(values, pivots_heap,std::floor(pos + 0.5), length, selector);
                return (low + high) / 2;
            }

        }, /**
         * The method R_3 has the following formulae for index and estimates<br>
         * \( \begin{align}
         * &amp;index= Np \\
         * &amp;estimate= x_{\lfloor h \rceil}\, \\
         * &amp;min_limit = 0.5/N \\
         * \end{align}\)
         */
        R_3("R-3") 
        {
            //override
            protected double index(const double p, const int length) 
            {
                const double min_limit = 1d/2 / length;
                return Double.compare(p, min_limit) <= 0 ?
                        0 : std::rint(length * p);
            }

        }, /**
         * The method R_4 has the following formulae for index and estimates<br>
         * \( \begin{align}
         * &amp;index= Np\, \\
         * &amp;estimate= x_{\lfloor h \rfloor} + (h -
         * \lfloor h \rfloor) (x_{\lfloor h \rfloor + 1} - x_{\lfloor h
         * \rfloor}) \\
         * &amp;min_limit = 1/N \\
         * &amp;max_limit = 1 \\
         * \end{align}\)
         */
        R_4("R-4") 
        {
            //override
            protected double index(const double p, const int length) 
            {
                const double min_limit = 1.0/ length;
                const double max_limit = 1d;
                return Double.compare(p, min_limit) < 0 ? 0 :
                       Double.compare(p, max_limit) == 0 ? length : length * p;
            }

        }, /**
         * The method R_5 has the following formulae for index and estimates<br>
         * \( \begin{align}
         * &amp;index= Np + 1/2\\
         * &amp;estimate= x_{\lfloor h \rfloor} + (h -
         * \lfloor h \rfloor) (x_{\lfloor h \rfloor + 1} - x_{\lfloor h
         * \rfloor}) \\
         * &amp;min_limit = 0.5/N \\
         * &amp;max_limit = (N-0.5)/N
         * \end{align}\)
         */
        R_5("R-5") 
        {

            //override
            protected double index(const double p, const int length) 
            {
                const double min_limit = 1d/2 / length;
                const double max_limit = (length - 0.5) / length;
                return Double.compare(p, min_limit) < 0 ? 0 :
                       Double.compare(p, max_limit) >= 0 ?
                               length : length * p + 0.5;
            }
        }, /**
         * The method R_6 has the following formulae for index and estimates<br>
         * \( \begin{align}
         * &amp;index= (N + 1)p \\
         * &amp;estimate= x_{\lfloor h \rfloor} + (h -
         * \lfloor h \rfloor) (x_{\lfloor h \rfloor + 1} - x_{\lfloor h
         * \rfloor}) \\
         * &amp;min_limit = 1/(N+1) \\
         * &amp;max_limit = N/(N+1) \\
         * \end{align}\)
         * <p>
         * <b>Note:</b> This method computes the index in a manner very close to
         * the default Hipparchus Percentile existing implementation. However
         * the difference to be noted is in picking up the limits with which
         * first element (p&lt;1(N+1)) and last elements (p&gt;N/(N+1))are done.
         * While in default case; these are done with p=0 and p=1 respectively.
         */
        R_6("R-6") 
        {

            //override
            protected double index(const double p, const int length) 
            {
                const double min_limit = 1.0/ (length + 1);
                const double max_limit = 1.0* length / (length + 1);
                return Double.compare(p, min_limit) < 0 ? 0 :
                       Double.compare(p, max_limit) >= 0 ?
                               length : (length + 1) * p;
            }
        }, 
        /**
         * The method R_7 : Microsoft Excel style computation has the
         * following formulae for index and estimates.<br>
         * \( \begin{align}
         * &amp;index = (N-1)p + 1 \\
         * &amp;estimate = x_{\lfloor h \rfloor} + (h -
         * \lfloor h \rfloor) (x_{\lfloor h \rfloor + 1} - x_{\lfloor h
         * \rfloor}) \\
         * &amp;min_limit = 0 \\
         * &amp;max_limit = 1 \\
         * \end{align}\)
         */
        R_7("R-7") 
        {
            //override
            protected double index(const double p, const int length) 
            {
                const double min_limit = 0;
                const double max_limit = 1d;
                return Double.compare(p, min_limit) == 0 ? 0 :
                       Double.compare(p, max_limit) == 0 ?
                               length : 1 + (length - 1) * p;
            }

        }, 
        /**
         * The method R_8 has the following formulae for index and estimates<br>
         * \( \begin{align}
         * &amp;index = (N + 1/3)p + 1/3  \\
         * &amp;estimate = x_{\lfloor h \rfloor} + (h -
           \lfloor h \rfloor) (x_{\lfloor h \rfloor + 1} - x_{\lfloor h
         * \rfloor}) \\
         * &amp;min_limit = (2/3)/(N+1/3) \\
         * &amp;max_limit = (N-1/3)/(N+1/3) \\
         * \end{align}\)
         * <p>
         * As per Ref [2,3] this approach is most recommended as it provides
         * an approximate median-unbiased estimate regardless of distribution.
         */
        R_8("R-8") 
        {
            //override
            protected double index(const double p, const int length) 
            {
                const double min_limit = 2 * (1.0/ 3) / (length + 1.0/ 3);
                const double max_limit =
                        (length - 1.0/ 3) / (length + 1.0/ 3);
                return Double.compare(p, min_limit) < 0 ? 0 :
                       Double.compare(p, max_limit) >= 0 ? length :
                           (length + 1.0/ 3) * p + 1.0/ 3;
            }
        }, 
        /**
         * The method R_9 has the following formulae for index and estimates<br>
         * \( \begin{align}
         * &amp;index = (N + 1/4)p + 3/8\\
         * &amp;estimate = x_{\lfloor h \rfloor} + (h -
           \lfloor h \rfloor) (x_{\lfloor h \rfloor + 1} - x_{\lfloor h
         * \rfloor}) \\
         * &amp;min_limit = (5/8)/(N+1/4) \\
         * &amp;max_limit = (N-3/8)/(N+1/4) \\
         * \end{align}\)
         */
        R_9("R-9") 
        {
            //override
            protected double index(const double p, const int length) 
            {
                const double min_limit = 5d/8 / (length + 0.25);
                const double max_limit = (length - 3d/8) / (length + 0.25);
                return Double.compare(p, min_limit) < 0 ? 0 :
                       Double.compare(p, max_limit) >= 0 ? length :
                               (length + 0.25) * p + 3d/8;
            }

        }, ;

        /** Simple name such as R-1, R-2 corresponding to those in wikipedia. */
        private const std::string name;

        /**
         * Constructor
         *
         * @param type name of estimation type as per wikipedia
         */
        Estimation_Type(const std::string type) 
        {
            this.name = type;
        }

        /**
         * Finds the index of array that can be used as starting index to
         * {@link #estimate(std::vector<double>, std::vector<int>, double, int, Kth_Selector) estimate}
         * percentile. The calculation of index calculation is specific to each
         * {@link Estimation_Type}.
         *
         * @param p the p<sup>th</sup> quantile
         * @param length the total number of array elements in the work array
         * @return a computed real valued index as explained in the wikipedia
         */
        protected virtual double index(double p, int length);

        /**
         * Estimation based on K<sup>th</sup> selection. This may be overridden
         * in specific enums to compute slightly different estimations.
         *
         * @param work array of numbers to be used for finding the percentile
         * @param pos indicated positional index prior computed from calling
         *            {@link #index(double, int)}
         * @param pivots_heap an earlier populated cache if exists; will be used
         * @param length size of array considered
         * @param selector a {@link Kth_Selector} used for pivoting during search
         * @return estimated percentile
         */
        protected double estimate(const std::vector<double> work, const std::vector<int> pivots_heap, const double pos, const int length, const Kth_Selector selector) 
        {

            const double fpos = std::floor(pos);
            const int int_pos = static_cast<int>( fpos;
            const double dif = pos - fpos;

            if (pos < 1) 
            {
                return selector.select(work, pivots_heap, 0);
            }
            if (pos >= length) 
            {
                return selector.select(work, pivots_heap, length - 1);
            }

            const double lower = selector.select(work, pivots_heap, int_pos - 1);
            const double upper = selector.select(work, pivots_heap, int_pos);
            return lower + dif * (upper - lower);
        }

        /**
         * Evaluate method to compute the percentile for a given bounded array
         * using earlier computed pivots heap.<br>
         * This basically calls the {@link #index(double, int) index} and then
         * {@link #estimate(std::vector<double>, std::vector<int>, double, int, Kth_Selector) estimate}
         * functions to return the estimated percentile value.
         *
         * @param work array of numbers to be used for finding the percentile
         * @param pivots_heap a prior cached heap which can speed up estimation
         * @param p the p<sup>th</sup> quantile to be computed
         * @param selector a {@link Kth_Selector} used for pivoting during search
         * @return estimated percentile
         * @ if p is out of range
         * @ if work array is NULL
         */
        protected double evaluate(const std::vector<double> work, const std::vector<int> pivots_heap, const double p, const Kth_Selector selector) 
        {
            //Math_Utils::check_not_null(work);
            if (p > 100 || p <= 0) 
            {
                throw (Localized_Stat_Formats.OUT_OF_BOUNDS_QUANTILE_VALUE, p, 0, 100);
            }
            return estimate(work, pivots_heap, index(p/100d, work.size()), work.size(), selector);
        }

        /**
         * Evaluate method to compute the percentile for a given bounded array.
         * This basically calls the {@link #index(double, int) index} and then
         * {@link #estimate(std::vector<double>, std::vector<int>, double, int, Kth_Selector) estimate}
         * functions to return the estimated percentile value. Please
         * note that this method does not make use of cached pivots.
         *
         * @param work array of numbers to be used for finding the percentile
         * @param p the p<sup>th</sup> quantile to be computed
         * @return estimated percentile
         * @param selector a {@link Kth_Selector} used for pivoting during search
         * @ if length or p is out of range
         * @ if work array is NULL
         */
        public double evaluate(const std::vector<double> work, const double p, const Kth_Selector selector) 
        {
            return this.evaluate(work, NULL, p, selector);
        }

        /**
         * Gets the name of the enum
         *
         * @return the name
         */
        std::string get_name() 
        {
            return name;
        }
    }
}


