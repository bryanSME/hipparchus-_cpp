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

//package org.hipparchus.stat.fitting;

//import java.io.Buffered_Reader;
//import java.io.Closeable;
//import java.io.File;
//import java.io.IOException;
//import java.io.Input_Stream;
//import java.io.Input_StreamReader;
//import java.net.URL;
//import java.nio.charset.Charset;
//import java.nio.file.Files;
//import java.util.Array_list;
//import java.util.List;

//import org.hipparchus.distribution.Real_Distribution;
//import org.hipparchus.distribution.continuous.Abstract_Real_Distribution;
//import org.hipparchus.distribution.continuous.ConstantReal_Distribution;
//import org.hipparchus.distribution.continuous.Normal_Distribution;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.exception.Math_Runtime_Exception;
//import org.hipparchus.exception.;
//import org.hipparchus.random.Random_Data_Generator;
//import org.hipparchus.random.Random_Generator;
//import org.hipparchus.stat.descriptive.Statistical_Summary;
//import org.hipparchus.stat.descriptive.Streaming_Statistics;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Utils;
#include <vector>
#include <string>
#include "../../core/distribution/continuous/AbstractRealDistribution.h"

/**
 * <p>Represents an <a href="http://http://en.wikipedia.org/wiki/Empirical_distribution_function">
 * empirical probability distribution</a> -- a probability distribution derived
 * from observed data without making any assumptions about the functional form
 * of the population distribution that the data come from.</p>
 *
 * <p>An <code>Empirical_Distribution</code> maintains data structures, called
 * <i>distribution digests</i>, that describe empirical distributions and
 * support the following operations: <ul>
 * <li>loading the distribution from a file of observed data values</li>
 * <li>dividing the input data into "bin ranges" and reporting bin frequency
 *     counts (data for histogram)</li>
 * <li>reporting univariate statistics describing the full set of data values
 *     as well as the observations within each bin</li>
 * <li>generating random values from the distribution</li>
 * </ul>
 * Applications can use <code>Empirical_Distribution</code> to build grouped
 * frequency histograms representing the input data or to generate random values
 * "like" those in the input file -- i.e., the values generated will follow the
 * distribution of the values in the file.</p>
 *
 * <p>The implementation uses what amounts to the
 * <a href="http://nedwww.ipac.caltech.edu/level5/March02/Silverman/Silver2_6.html">
 * Variable Kernel Method</a> with Gaussian smoothing:<p>
 * <strong>Digesting the input file</strong>
 * <ol><li>Pass the file once to compute min and max.</li>
 * <li>Divide the range from min-max into <code>bin_count</code> "bins."</li>
 * <li>Pass the data file again, computing bin counts and univariate
 *     statistics (mean, std dev.) for each of the bins </li>
 * <li>Divide the interval (0,1) into subintervals associated with the bins, *     with the length of a bin's subinterval proportional to its count.</li></ol>
 * <strong>Generating random values from the distribution</strong><ol>
 * <li>Generate a uniformly distributed value in (0,1) </li>
 * <li>Select the subinterval to which the value belongs.
 * <li>Generate a random Gaussian value with mean = mean of the associated
 *     bin and std dev = std dev of associated bin.</li></ol></p>
 *
 * <p>Empirical_Distribution : the {@link Real_Distribution} interface
 * as follows.  Given x within the range of values in the dataset, let B
 * be the bin containing x and let K be the within-bin kernel for B.  Let P(B-)
 * be the sum of the probabilities of the bins below B and let K(B) be the
 * mass of B under K (i.e., the integral of the kernel density over B).  Then
 * set P(X &lt; x) = P(B-) + P(B) * K(x) / K(B) where K(x) is the kernel distribution
 * evaluated at x. This results in a cdf that matches the grouped frequency
 * distribution at the bin endpoints and interpolates within bins using
 * within-bin kernels.</p>
 *
 *<strong>USAGE NOTES:</strong><ul>
 *<li>The <code>bin_count</code> is set by default to 1000.  A good rule of thumb
 *    is to set the bin count to approximately the length of the input file divided
 *    by 10. </li>
 *<li>The input file <i>must</i> be a plain text file containing one valid numeric
 *    entry per line.</li>
 * </ul></p>
 *
 */
class Empirical_Distribution : Abstract_Real_Distribution 
{
public:
    /** Default bin count */
    static constexpr int DEFAULT_BIN_COUNT{ 1000 };

private:
    /** Character set for file input */
    static const std::string FILE_CHARSET{ "US-ASCII" };

    /** Random_Data_Generator instance to use in repeated calls to get_next() */
    protected const Random_Data_Generator random_data;

    /** List of Summary_Statistics objects characterizing the bins */
    private const std::vector<Streaming_Statistics> bin_stats;

    /** Sample statistics */
    private Streaming_Statistics sample_stats;

    /** Max loaded value */
    private double max = -INFINITY;

    /** Min loaded value */
    private double min = INFINITY;

    /** Grid size */
    private double delta;

    /** number of bins */
    private const int bin_count;

    /** is the distribution loaded? */
    private bool loaded;

    /** upper bounds of subintervals in (0,1) "belonging" to the bins */
    private std::vector<double> upper_bounds;

    /**
     * Creates a Empirical_Distribution with the default bin count.
     */
    public Empirical_Distribution() 
    {
        this(DEFAULT_BIN_COUNT);
    }

    /**
     * Creates a Empirical_Distribution with the specified bin count.
     *
     * @param bin_count number of bins. Must be strictly positive.
     * @ if {@code bin_count <= 0}.
     */
    public Empirical_Distribution(const int& bin_count) 
    {
        this(bin_count, Random_Data_Generator());
    }

    /**
     * Creates a Empirical_Distribution with the specified bin count using the
     * provided {@link Random_Generator} as the source of random data.
     *
     * @param bin_count number of bins. Must be strictly positive.
     * @param generator random data generator (may be NULL, resulting in default JDK generator)
     * @ if {@code bin_count <= 0}.
     */
    public Empirical_Distribution(const int& bin_count, Random_Generator generator) 
    {
        this(bin_count, Random_Data_Generator.of(generator));
    }

    /**
     * Creates a Empirical_Distribution with default bin count using the
     * provided {@link Random_Generator} as the source of random data.
     *
     * @param generator random data generator (may be NULL, resulting in default JDK generator)
     */
    public Empirical_Distribution(Random_Generator generator) 
    {
        this(DEFAULT_BIN_COUNT, generator);
    }

    /**
     * Private constructor to allow lazy initialisation of the RNG contained
     * in the {@link #random_data} instance variable.
     *
     * @param bin_count number of bins. Must be strictly positive.
     * @param random_data Random data generator.
     * @ if {@code bin_count <= 0}.
     */
    private Empirical_Distribution(const int& bin_count, Random_Data_Generator random_data) 
    {
        if (bin_count <= 0) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL_BOUND_EXCLUDED, bin_count, 0);
        }
        this.bin_count = bin_count;
        this.random_data = random_data;
        bin_stats = Array_list<>();
    }

    /**
     * Computes the empirical distribution from the provided
     * array of numbers.
     *
     * @param in the input data array
     * @exception  if in is NULL
     */
    public void load(std::vector<double> in)  
    {
        try (Data_Adapter da = Array_dataAdapter(in)) 
        {
            da.compute_stats();
            // adapter for the second pass
            fill_bin_stats(new Array_dataAdapter(in));
        }
catch (IOException ex) 
        {
            // Can't happen
            throw Math_Runtime_Exception.create_internal_error();
        }
        loaded = true;

    }

    /**
     * Computes the empirical distribution using data read from a URL.
     *
     * <p>The input file <i>must</i> be an ASCII text file containing one
     * valid numeric entry per line.</p>
     *
     * @param url url of the input file
     *
     * @IOException if an IO error occurs
     * @ if url is NULL
     * @ if URL contains no data
     */
    public void load(URL url) IOException,  
    {
        //Math_Utils::check_not_null(url);
        Charset charset = Charset.for_name(FILE_CHARSET);
        try (Input_Stream       is1  = url.open_stream();
             Input_StreamReader isr1 = Input_StreamReader(is1, charset);
             Buffered_Reader    br1  = Buffered_Reader(isr1);
             Data_Adapter       da1  = Stream_dataAdapter(br1)) 
             {
            da1.compute_stats();
            if (sample_stats.get_n() == 0) 
            {
                throw std::exception("not implemented");
                //throw (hipparchus::exception::Localized_Core_Formats_Type::URL_CONTAINS_NO_DATA, url);
            }
            // adapter for the second pass
            try (Input_Stream       is2  = url.open_stream();
                 Input_StreamReader isr2 = Input_StreamReader(is2, charset);
                 Buffered_Reader    br2  = Buffered_Reader(isr2);
                 Data_Adapter       da2  = Stream_dataAdapter(br2)) 
                 {
                fill_bin_stats(da2);
                loaded = true;
            }
        }
    }

    /**
     * Computes the empirical distribution from the input file.
     *
     * <p>The input file <i>must</i> be an ASCII text file containing one
     * valid numeric entry per line.</p>
     *
     * @param file the input file
     * @IOException if an IO error occurs
     * @ if file is NULL
     */
    public void load(File file) IOException 
    {
        //Math_Utils::check_not_null(file);
        Charset charset = Charset.for_name(FILE_CHARSET);
        try (Input_Stream    is1 = Files.new_input_stream(file.to_path());
             Buffered_Reader br1 = Buffered_Reader(new Input_StreamReader(is1, charset));
             Data_Adapter    da1 = Stream_dataAdapter(br1)) 
             {
            da1.compute_stats();
            // adapter for second pass
            try (Input_Stream    is2 = Files.new_input_stream(file.to_path());
                 Buffered_Reader in2 = Buffered_Reader(new Input_StreamReader(is2, charset));
                 Data_Adapter    da2 = Stream_dataAdapter(in2)) 
                 {
                fill_bin_stats(da2);
            }
            loaded = true;
        }
    }

    /**
     * Provides methods for computing <code>sample_stats</code> and
     * <code>bean_stats</code> abstracting the source of data.
     */
    private virtual class Data_Adapter : Closeable 
    {

        /**
         * Compute bin stats.
         *
         * @IOException  if an error occurs computing bin stats
         */
        public virtual void compute_bin_stats() IOException;

        /**
         * Compute sample statistics.
         *
         * @IOException if an error occurs computing sample stats
         */
        public virtual void compute_stats() IOException;

    }

    /**
     * <code>Data_Adapter</code> for data provided through some input stream
     */
    private class Stream_dataAdapter extends Data_Adapter 
    {

        /** Input stream providing access to the data */
        private Buffered_Reader input_stream;

        /**
         * Create a Stream_dataAdapter from a Buffered_Reader
         *
         * @param in Buffered_Reader input stream
         */
        Stream_dataAdapter(Buffered_Reader in) 
        {
            input_stream = in;
        }

        /** {@inherit_doc} */
        //override
        public void compute_bin_stats() IOException 
        {
            for (std::string str = input_stream.read_line(); str != NULL; str = input_stream.read_line()) 
            {
                const double val = Double.parse_double(str);
                Streaming_Statistics stats = bin_stats.get(find_bin(val));
                stats.add_value(val);
            }
        }

        /** {@inherit_doc} */
        //override
        public void compute_stats() IOException 
        {
            sample_stats = Streaming_Statistics();
            for (std::string str = input_stream.read_line(); str != NULL; str = input_stream.read_line()) 
            {
                const double val = Double.parse_double(str);
                sample_stats.add_value(val);
            }
        }

        /** {@inherit_doc} */
        //override
        public void close() IOException 
        {
            if (input_stream != NULL) 
            {
                input_stream.close();
                input_stream = NULL;
            }
        }

    }

    /**
     * <code>Data_Adapter</code> for data provided as array of doubles.
     */
    private class Array_dataAdapter extends Data_Adapter 
    {

        /** Array of input  data values */
        private const std::vector<double> input_array;

        /**
         * Construct an Array_dataAdapter from a std::vector<double> array
         *
         * @param in std::vector<double> array holding the data, a reference to the array will be stored
         * @ if in is NULL
         */
        Array_dataAdapter(std::vector<double> in)  { // NOPMD - storing a reference to the array is intentional and documented here
            super();
            //Math_Utils::check_not_null(in);
            input_array = in;
        }

        /** {@inherit_doc} */
        //override
        public void compute_stats() 
        {
            sample_stats = Streaming_Statistics();
            for (int i{}; i < input_array.size(); i++) 
            {
                sample_stats.add_value(input_array[i]);
            }
        }

        /** {@inherit_doc} */
        //override
        public void compute_bin_stats() 
        {
            for (int i{}; i < input_array.size(); i++) 
            {
                Streaming_Statistics stats =
                    bin_stats.get(find_bin(input_array[i]));
                stats.add_value(input_array[i]);
            }
        }

        /** {@inherit_doc} */
        //override
        public void close() 
        {
            // nothing to do
        }

    }

    /**
     * Fills bin_stats array (second pass through data file).
     *
     * @param da object providing access to the data
     * @IOException  if an IO error occurs
     */
    private void fill_bin_stats(const Data_Adapter da)
        IOException 
        {
        // Set up grid
        min = sample_stats.get_min();
        max = sample_stats.get_max();
        delta = (max - min)/bin_count;

        // Initialize bin_stats Array_list
        if (!bin_stats.is_empty()) 
        {
            bin_stats.clear();
        }
        for (int i{}; i < bin_count; i++) 
        {
            Streaming_Statistics stats = Streaming_Statistics();
            bin_stats.add(i,stats);
        }

        // Filling data in bin_stats Array
        da.compute_bin_stats();

        // Assign upper_bounds based on bin counts
        upper_bounds = std::vector<double>(bin_count];
        upper_bounds[0] =
        (static_cast<double>( bin_stats.get(0).get_n()) / static_cast<double>( sample_stats.get_n();
        for (int i{ 1 }; i < bin_count-1; i++) 
        {
            upper_bounds[i] = upper_bounds[i-1] +
            (static_cast<double>( bin_stats.get(i).get_n()) / static_cast<double>( sample_stats.get_n();
        }
        upper_bounds[bin_count-1] = 1.0;
    }

    /**
     * Returns the index of the bin to which the given value belongs
     *
     * @param value  the value whose bin we are trying to find
     * @return the index of the bin containing the value
     */
    private int find_bin(double value) 
    {
        return std::min(
                std::max(static_cast<int>( std::ceil((value - min) / delta) - 1, 0), bin_count - 1);
    }

    /**
     * Generates a random value from this distribution.
     * <strong>Preconditions:</strong><ul>
     * <li>the distribution must be loaded before invoking this method</li></ul>
     * @return the random value.
     * @Math_Illegal_State_Exception if the distribution has not been loaded
     */
    public double get_next_value() Math_Illegal_State_Exception 
    {

        if (!loaded) 
        {
            throw Math_Illegal_State_Exception(hipparchus::exception::Localized_Core_Formats_Type::DISTRIBUTION_NOT_LOADED);
        }

        return inverse_cumulative_probability(random_data.next_double());
    }

    /**
     * Returns a {@link Statistical_Summary} describing this distribution.
     * <strong>Preconditions:</strong><ul>
     * <li>the distribution must be loaded before invoking this method</li></ul>
     *
     * @return the sample statistics
     * @Illegal_State_Exception if the distribution has not been loaded
     */
    public Statistical_Summary get_sample_stats() 
    {
        return sample_stats;
    }

    /**
     * Returns the number of bins.
     *
     * @return the number of bins.
     */
    public int get_bin_count() 
    {
        return bin_count;
    }

    /**
     * Returns a List of {@link Streaming_Statistics} instances containing
     * statistics describing the values in each of the bins.  The list is
     * indexed on the bin number.
     *
     * @return List of bin statistics.
     */
    public List<Streaming_Statistics> get_bin_stats() 
    {
        return bin_stats;
    }

    /**
     * <p>Returns a fresh copy of the array of upper bounds for the bins.
     * Bins are: <br/>
     * [min,upper_bounds[0]],(upper_bounds[0],upper_bounds[1]],..., *  (upper_bounds[bin_count-2], upper_bounds[bin_count-1] = max].</p>
     *
     * @return array of bin upper bounds
     */
    public std::vector<double> get_upper_bounds() 
    {
        std::vector<double> bin_upper_bounds = std::vector<double>(bin_count];
        for (int i{}; i < bin_count - 1; i++) 
        {
            bin_upper_bounds[i] = min + delta * (i + 1);
        }
        bin_upper_bounds[bin_count - 1] = max;
        return bin_upper_bounds;
    }

    /**
     * <p>Returns a fresh copy of the array of upper bounds of the subintervals
     * of [0,1] used in generating data from the empirical distribution.
     * Subintervals correspond to bins with lengths proportional to bin counts.</p>
     *
     * <strong>Preconditions:</strong><ul>
     * <li>the distribution must be loaded before invoking this method</li></ul>
     *
     * @return array of upper bounds of subintervals used in data generation
     * @Null_Pointer_Exception unless a {@code load} method has been
     * called beforehand.
     */
    public std::vector<double> get_generator_upper_bounds() 
    {
        int len = upper_bounds.size();
        std::vector<double> out = std::vector<double>(len];
        System.arraycopy(upper_bounds, 0, out, 0, len);
        return out;
    }

    /**
     * Property indicating whether or not the distribution has been loaded.
     *
     * @return true if the distribution has been loaded
     */
    public bool is_loaded() 
    {
        return loaded;
    }

    /**
     * Reseeds the random number generator used by {@link #get_next_value()}.
     *
     * @param seed random generator seed
     */
    public void re_seed(long seed) 
    {
        random_data.set_seed(seed);
    }

    // Distribution methods ---------------------------

    /**
     * {@inherit_doc}
     *
     * <p>Returns the kernel density normalized so that its integral over each bin
     * equals the bin mass.</p>
     *
     * <p>Algorithm description: <ol>
     * <li>Find the bin B that x belongs to.</li>
     * <li>Compute K(B) = the mass of B with respect to the within-bin kernel (i.e., the
     * integral of the kernel density over B).</li>
     * <li>Return k(x) * P(B) / K(B), where k is the within-bin kernel density
     * and P(B) is the mass of B.</li></ol></p>
     */
    //override
    public double density(const double& x) 
    {
        if (x < min || x > max) 
        {
            return 0;
        }
        const int bin_index = find_bin(x);
        const Real_Distribution kernel = get_kernel(bin_stats.get(bin_index));
        return kernel.density(x) * pB(bin_index) / kB(bin_index);
    }

    /**
     * {@inherit_doc}
     *
     * <p>Algorithm description:<ol>
     * <li>Find the bin B that x belongs to.</li>
     * <li>Compute P(B) = the mass of B and P(B-) = the combined mass of the bins below B.</li>
     * <li>Compute K(B) = the probability mass of B with respect to the within-bin kernel
     * and K(B-) = the kernel distribution evaluated at the lower endpoint of B</li>
     * <li>Return P(B-) + P(B) * [K(x) - K(B-)] / K(B) where
     * K(x) is the within-bin kernel distribution function evaluated at x.</li></ol>
     * If K is a constant distribution, we return P(B-) + P(B) (counting the full
     * mass of B).</p>
     *
     */
    //override
    public double cumulative_probability(const double& x) 
    {
        if (x < min) 
        {
            return 0;
        }
        if (x >= max) 
        {
            return 1;
        }
        const int bin_index = find_bin(x);
        const double p_bminus = p_bminus(bin_index);
        const double pB = pB(bin_index);
        const Real_Distribution kernel = k(x);
        if (kernel instanceof ConstantReal_Distribution) 
        {
            return x < kernel.get_numerical_mean()
                ? p_bminus
                : p_bminus + pB;
        }
        const std::vector<double> bin_bounds = get_upper_bounds();
        const double kB = kB(bin_index);
        const double lower = bin_index == 0
            ? min
            : bin_bounds[bin_index - 1];
        const double within_bin_cum = (kernel.cumulative_probability(x) - kernel.cumulative_probability(lower)) / kB;
        return p_bminus + pB * within_bin_cum;
    }

    /**
     * {@inherit_doc}
     *
     * <p>Algorithm description:<ol>
     * <li>Find the smallest i such that the sum of the masses of the bins
     *  through i is at least p.</li>
     * <li>
     *   Let K be the within-bin kernel distribution for bin i.</br>
     *   Let K(B) be the mass of B under K. <br/>
     *   Let K(B-) be K evaluated at the lower endpoint of B (the combined
     *   mass of the bins below B under K).<br/>
     *   Let P(B) be the probability of bin i.<br/>
     *   Let P(B-) be the sum of the bin masses below bin i. <br/>
     *   Let p_crit = p - P(B-)<br/>
     * <li>Return the inverse of K evaluated at <br/>
     *    K(B-) + p_crit * K(B) / P(B) </li>
     *  </ol></p>
     *
     */
    //override
    public double inverse_cumulative_probability(const double& p)  
    {
        Math_Utils::check_range_inclusive(p, 0, 1);

        if (p == 0.0) 
        {
            return get_support_lower_bound();
        }

        if (p == 1.0) 
        {
            return get_support_upper_bound();
        }

        int i = 0;
        while (cum_bin_p(i) < p) 
        {
            i++;
        }

        const Real_Distribution kernel = get_kernel(bin_stats.get(i));
        const double kB = kB(i);
        const std::vector<double> bin_bounds = get_upper_bounds();
        const double lower = i == 0 ? min : bin_bounds[i - 1];
        const double k_bminus = kernel.cumulative_probability(lower);
        const double pB = pB(i);
        const double p_bminus = p_bminus(i);
        const double p_crit = p - p_bminus;
        if (p_crit <= 0) 
        {
            return lower;
        }
        return kernel.inverse_cumulative_probability(k_bminus + p_crit * kB / pB);
    }

    /**
     * {@inherit_doc}
     */
    //override
    public double get_numerical_mean() const 
    {
       return sample_stats.get_mean();
    }

    /**
     * {@inherit_doc}
     */
    //override
    public double get_numerical_variance() const 
    {
        return sample_stats.get_variance();
    }

    /**
     * {@inherit_doc}
     */
    //override
    public double get_support_lower_bound() const 
    {
       return min;
    }

    /**
     * {@inherit_doc}
     */
    //override
    public double get_support_upper_bound() const 
    {
        return max;
    }

    /**
     * {@inherit_doc}
     */
    //override
    public bool is_support_connected() const 
    {
        return true;
    }

    /**
     * Reseed the underlying PRNG.
     *
     * @param seed seed value
     */
    public void reseed_random_generator(long seed) 
    {
        random_data.set_seed(seed);
    }

    /**
     * The probability of bin i.
     *
     * @param i the index of the bin
     * @return the probability that selection begins in bin i
     */
    private double pB(const int& i) 
    {
        return i == 0 ? upper_bounds[0] :
            upper_bounds[i] - upper_bounds[i - 1];
    }

    /**
     * The combined probability of the bins up to but not including bin i.
     *
     * @param i the index of the bin
     * @return the probability that selection begins in a bin below bin i.
     */
    private double p_bminus(const int& i) 
    {
        return i == 0 ? 0 : upper_bounds[i - 1];
    }

    /**
     * Mass of bin i under the within-bin kernel of the bin.
     *
     * @param i index of the bin
     * @return the difference in the within-bin kernel cdf between the
     * upper and lower endpoints of bin i
     */
    private double kB(const int& i) 
    {
        const std::vector<double> bin_bounds = get_upper_bounds();
        const Real_Distribution kernel = get_kernel(bin_stats.get(i));
        return i == 0 ? kernel.probability(min, bin_bounds[0]) :
            kernel.probability(bin_bounds[i - 1], bin_bounds[i]);
    }

    /**
     * The within-bin kernel of the bin that x belongs to.
     *
     * @param x the value to locate within a bin
     * @return the within-bin kernel of the bin containing x
     */
    private Real_Distribution k(double x) 
    {
        const int bin_index = find_bin(x);
        return get_kernel(bin_stats.get(bin_index));
    }

    /**
     * The combined probability of the bins up to and including bin_index.
     *
     * @param bin_index maximum bin index
     * @return sum of the probabilities of bins through bin_index
     */
    private double cum_bin_p(const int& bin_index) 
    {
        return upper_bounds[bin_index];
    }

    /**
     * The within-bin smoothing kernel. Returns a Gaussian distribution
     * parameterized by {@code b_stats}, unless the bin contains less than 2
     * observations, in which case a constant distribution is returned.
     *
     * @param b_stats summary statistics for the bin
     * @return within-bin kernel parameterized by b_stats
     */
    protected Real_Distribution get_kernel(Streaming_Statistics b_stats) 
    {
        if (b_stats.get_n() < 2 || b_stats.get_variance() == 0) 
        {
            return ConstantReal_Distribution(b_stats.get_mean());
        }
else 
        {
            return Normal_Distribution(b_stats.get_mean(), b_stats.get_standard_deviation());
        }
    }
}


