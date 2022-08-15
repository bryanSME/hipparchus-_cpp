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

//import java.io.IOException;
//import java.io.Object_Input_Stream;
//import java.io.Serializable;
//import java.text.Decimal_Format;
//import java.util.Array_list;
//import java.util.Arrays;
//import java.util.Collection;
//import java.util.Collections;
//import java.util.List;

//import org.hipparchus.analysis.Univariate_Function;
//import org.hipparchus.analysis.interpolation.Linear_Interpolator;
//import org.hipparchus.analysis.interpolation.Neville_Interpolator;
//import org.hipparchus.analysis.interpolation.Univariate_Interpolator;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.stat.descriptive.Abstract_Storeless_Univariate_Statistic;
//import org.hipparchus.stat.descriptive.Storeless_Univariate_Statistic;
//import org.hipparchus.util.Math_Arrays;
//import org.hipparchus.util.Math_Utils;
//import org.hipparchus.util.Precision;
#include <vector>
#include <string>

/**
 * A {@link Storeless_Univariate_Statistic} estimating percentiles using the
 * <a href=http://www.cs.wustl.edu/~jain/papers/ftp/psqr.pdf>P<SUP>2</SUP></a>
 * Algorithm as explained by <a href=http://www.cse.wustl.edu/~jain/>Raj
 * Jain</a> and Imrich Chlamtac in
 * <a href=http://www.cse.wustl.edu/~jain/papers/psqr.htm>P<SUP>2</SUP> Algorithm
 * for Dynamic Calculation of Quantiles and Histogram Without Storing
 * Observations</a>.
 * <p>
 * Note: This implementation is not synchronized and produces an approximate
 * result. For small samples, where data can be stored and processed in memory, * {@link Percentile} should be used.
 */
class P_Square_Percentile : public Abstract_Storeless_Univariate_Statistic, public Storeless_Univariate_Statistic
{
private:
    /** The maximum array size used for psquare algorithm */
    static const int PSQUARE_CONSTANT = 5;

    /**
     * A Default quantile needed in case if user prefers to use default no
     * argument constructor.
     */
    static constexpr double DEFAULT_QUANTILE_DESIRED{ 50 };

    /** A decimal formatter for print convenience */
    static const Decimal_Format DECIMAL_FORMAT = Decimal_Format("00.00");

    /**
     * Initial list of 5 numbers corresponding to 5 markers. <b>NOTE:</b>watch
     * out for the add methods that are overloaded
     */
    const std::vector<double> my_initial_five{ std::vector<double>(PSQUARE_CONSTANT) };

    /**
     * The quantile needed should be in range of 0-1. The constructor
     * {@link #P_Square_Percentilestatic_cast<double>(} ensures that passed in percentile is
     * divided by 100.
     */
    double my_quantile;

    /**
     * last_observation is the last observation value/input sample. No need to
     * serialize.
     */
    double my_last_observation;

    /**
     * Markers is the marker collection object which comes to effect
     * only after 5 values are inserted
     */
    P_Square_Markers my_markers;

    /**
     * Computed p value (i,e percentile value of data set hither to received)
     */
    double my_p_value = std::numeric_limits<double>::quiet_NaN();

    /**
     * Counter to count the values/observations accepted into this data set
     */
    long my_count_of_observations;

public:
    /**
     * Constructs a P_Square_Percentile with the specific percentile value.
     * @param p the percentile
     * @  if p is not greater than 0 and less
     * than or equal to 100
     */
    P_Square_Percentile(const double& p) 
    {
        if (p > 100 || p < 0) 
        {
            throw std::exception("hipparchus::exception::Localized_Core_Formats_Type::OUT_OF_RANGE, p, 0, 100)");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::OUT_OF_RANGE, p, 0, 100);
        }
        my_quantile = p / 100.0;// always set it within (0,1]
    }

    /**
     * Default constructor that assumes a {@link #DEFAULT_QUANTILE_DESIRED
     * default quantile} needed.
     */
    P_Square_Percentile() 
    {
        P_Square_Percentile(DEFAULT_QUANTILE_DESIRED);
    }

    /**
     * Copy constructor, creates a {@code P_Square_Percentile} identical
     * to the {@code original}.
     *
     * @param original the {@code P_Square_Percentile} instance to copy
     * @org.hipparchus.exception.Null_Argument_Exception if original is NULL
     */
    P_Square_Percentile(P_Square_Percentile original) 
    {
        //super();

        my_quantile = original.quantile;

        if (original.markers != NULL) 
        {
            my_markers = original.markers.copy_self();
        }

        my_count_of_observations = original.count_of_observations;
        my_p_value = original.p_value;
        my_initial_five.add_all(original.initial_five);
    }

    /** {@inherit_doc} */
    //override
    int hash_code() 
    {
        double result = get_result();
        result = std::isnan(result) ? 37 : result;
        const double markers_hash = my_markers == NULL ? 0 : my_markers.hash_code();
        const std::vector<double> to_hash = {result, my_quantile, markers_hash, my_count_of_observations};
        return Arrays.hash_code(to_hash);
    }

    /**
     * Returns true iff {@code o} is a {@code P_Square_Percentile} returning the
     * same values as this for {@code get_result()} and {@code get_n()} and also
     * having equal markers
     *
     * @param o object to compare
     * @return true if {@code o} is a {@code P_Square_Percentile} with
     * equivalent internal state
     */
    //override
    bool equals(const Object& other) 
    {
        if (this == other)
        {
            return true;
        }
        bool result{};
        if (other instanceof P_Square_Percentile)
        {
            auto that = static_cast<P_Square_Percentile>(other);
            // markers as in the case of first
            // five observations
            result = my_markers.equals(that.markers) && get_n() == that.get_n();
        }
        return result;
    }

    /**
     * {@inherit_doc}The internal state updated due to the value in this
     * context is basically of the marker positions and computation of the
     * approximate quantile.
     *
     * @param observation the observation currently being added.
     */
    //override
    void increment(const double& observation) 
    {
        // Increment counter
        count_of_observations++;

        // Store last observation
        my_last_observation = observation;

        // 0. Use Brute force for <5
        if (my_markers == NULL) 
        {
            if (my_initial_five.add(observation)) 
            {
                Collections.sort(my_initial_five);
                my_p_value = my_initial_five.at(static_cast<int>( (my_quantile * (my_initial_five.size() - 1))));
                return;
            }
            // 1. Initialize once after 5th observation
            my_markers = new_markers(my_initial_five, my_quantile);
        }
        // 2. process a Data Point and return p_value
        my_p_value = my_markers.process_data_point(observation);
    }

    /**
     * Returns a string containing the last observation, the current estimate
     * of the quantile and all markers.
     *
     * @return string representation of state data
     */
    //override
    std::string to_string() const 
    {
        synchronized (this) 
        {
            synchronized (DECIMAL_FORMAT) 
            {
                if (markers == NULL) 
                {
                    return std::string.format("obs=%s p_value=%s", DECIMAL_FORMAT.format(last_observation), DECIMAL_FORMAT.format(p_value));
                }
                else 
                {
                    return std::string.format("obs=%s markers=%s", DECIMAL_FORMAT.format(last_observation), markers.to_string());
                }
            }
        }
   }

    /** {@inherit_doc} */
    //override
    long get_n() const
    {
        return my_count_of_observations;
    }

    /** {@inherit_doc} */
    //override
    P_Square_Percentile copy() 
    {
        return P_Square_Percentile(this);
    }

    /**
     * Returns the quantile estimated by this statistic in the range [0.0-1.0]
     *
     * @return quantile estimated by {@link #get_result()}
     */
    double quantile() const
    {
        return my_quantile;
    }

    /**
     * {@inherit_doc}. This basically clears all the markers, the
     * initial_five list and sets count_of_observations to 0.
     */
    //override
    void clear() 
    {
        my_markers = NULL;
        my_initial_five.clear();
        my_count_of_observations = 0L;
        my_p_value = std::numeric_limits<double>::quiet_NaN();
    }

    /**
     * {@inherit_doc}
     */
    //override
    double get_result() 
    {
        if (Double.compare(quantile, 1.0) == 0) 
        {
            my_p_value = maximum();
        }
        else if (Double.compare(quantile, 0.0) == 0) 
        {
            my_p_value = minimum();
        }
        return my_p_value;
    }

    /**
     * @return the quantile estimated by this statistic
     */
    double get_quantile() const
    {
        return my_quantile;
    }

    /**
     * @return maximum in the data set added to this statistic
     */
    private double maximum() 
    {
        double val = std::numeric_limits<double>::quiet_NaN();
        if (my_markers != NULL) 
        {
            val = my_markers.height(PSQUARE_CONSTANT);
        }
        else if (!my_initial_five.empty()) 
        {
            val = my_initial_five.at(my_initial_five.size() - 1);
        }
        return val;
    }

    /**
     * @return minimum in the data set added to this statistic
     */
    private double minimum() 
    {
        double val = std::numeric_limits<double>::quiet_NaN();
        if (my_markers != NULL)
        {
            val = my_markers.height(1);
        }
        else if (!my_initial_five.empty()) 
        {
            val = my_initial_five.at(0);
        }
        return val;
    }

    /**
     * Markers is an encapsulation of the five markers/buckets as indicated in
     * the original works.
     */
    private static class Markers : P_Square_Markers
    {
        /** Low marker index */
        private static constexpr int LOW{ 2 };

        /** High marker index */
        private static constexpr int HIGH{ 4 };

        /**
         * Array of 5+1 Markers (The first marker is dummy just so we
         * can match the rest of indexes [1-5] indicated in the original works
         * which follows unit based index)
         */
        private const std::vector<Marker> my_marker_array;

        /**
         * Kth cell belonging to [1-5] of the marker_array. No need for
         * this to be serialized
         */
        private int k{ -1 };

        /**
         * Constructor
         *
         * @param the_marker_array marker array to be used, a reference to the array will be stored
         */
        private Markers(const std::vector<Marker> the_marker_array) { // NOPMD - storing a reference to the array is intentional and documented here
            //Math_Utils::check_not_null(the_marker_array);
            marker_array = the_marker_array;
            for (int i{ 1 }; i < PSQUARE_CONSTANT; i++) 
            {
                marker_array[i].previous(marker_array[i - 1])
                        .next(marker_array[i + 1]).index(i);
            }
            marker_array[0].previous(marker_array[0])
                          .next(marker_array[1])
                          .index(0);
            marker_array[5].previous(marker_array[4])
                          .next(marker_array[5])
                          .index(5);
        }

        /**
         * Constructor
         *
         * @param initial_five elements required to build Marker
         * @param p quantile required to be computed
         */
        private Markers(const std::vector<double>& initial_five, const double& p) 
        {
            this(create_marker_array(initial_five, p));
        }

        /**
         * Creates a marker array using initial five elements and a quantile
         *
         * @param initial_five list of initial five elements
         * @param p the pth quantile
         * @return Marker array
         */
        private static std::vector<Marker> create_marker_array(const std::vector<double> initial_five, const double p) 
        {
            const int count_observed = initial_five == NULL ? -1 : initial_five.size();
            if (count_observed < PSQUARE_CONSTANT) 
            {
                throw (hipparchus::exception::Localized_Core_Formats_Type::INSUFFICIENT_OBSERVED_POINTS_IN_SAMPLE, count_observed, PSQUARE_CONSTANT);
            }
            Collections.sort(initial_five);
            return std::vector<Marker>{
                    Marker(),// Null Marker
                    Marker(initial_five.get(0), 1, 0, 1),
                    Marker(initial_five.get(1), 1 + 2 * p, p / 2, 2),
                    Marker(initial_five.get(2), 1 + 4 * p, p, 3),
                    Marker(initial_five.get(3), 3 + 2 * p, (1 + p) / 2, 4), 
                    Marker(initial_five.get(4), 5, 1, 5)
            };
        }

        /**
         * {@inherit_doc}
         */
        //override
        public int hash_code() 
        {
            return Arrays.deep_hash_code(marker_array);
        }

        /**
         * {@inherit_doc}.This equals method basically checks for marker array to
         * be deep equals.
         *
         * @param o is the other object
         * @return true if the object compares with this object are equivalent
         */
        //override
        public bool equals(Object o) 
        {
            bool result = false;
            if (this == o) 
            {
                result = true;
            }
            else if (o instanceof Markers) 
            {
                Markers that = (Markers) o;
                result = Arrays.deep_equals(marker_array, that.marker_array);
            }
            return result;
        }

        /**
         * Process a data point
         *
         * @param input_data_point is the data point passed
         * @return computed percentile
         */
        //override
        public double process_data_point(const double input_data_point) 
        {

            // 1. Find cell and update minima and maxima
            const int& kth_cell = find_cell_and_update_min_max(input_data_point);

            // 2. Increment positions
            increment_positions(1, kth_cell + 1, 5);

            // 2a. Update desired position with increments
            update_desired_positions();

            // 3. Adjust heights of m[2-4] if necessary
            adjust_heights_of_markers();

            // 4. Return percentile
            return get_percentile_value();
        }

        /**
         * Returns the percentile computed thus far.
         *
         * @return height of mid point marker
         */
        //override
        public double get_percentile_value() 
        {
            return height(3);
        }

        /**
         * Finds the cell where the input observation / value fits.
         *
         * @param observation the input value to be checked for
         * @return kth cell (of the markers ranging from 1-5) where observed
         *         sample fits
         */
        private int find_cell_and_update_min_max(const double observation) 
        {
            k = -1;
            if (observation < height(1)) 
            {
                marker_array[1].marker_height = observation;
                k = 1;
            }
            else if (observation < height(2)) 
            {
                k = 1;
            }
            else if (observation < height(3)) 
            {
                k = 2;
            }
            else if (observation < height(4)) 
            {
                k = 3;
            }
            else if (observation <= height(5)) 
            {
                k = 4;
            }
            else 
            {
                marker_array[5].marker_height = observation;
                k = 4;
            }
            return k;
        }

        /**
         * Adjust marker heights by setting quantile estimates to middle markers.
         */
        private void adjust_heights_of_markers() 
        {
            for (int i = LOW; i <= HIGH; i++) 
            {
                estimate(i);
            }
        }

        /**
         * {@inherit_doc}
         */
        //override
        public double estimate(const int index) 
        {
            Math_Utils::check_range_inclusive(index, LOW, HIGH);
            return marker_array[index].estimate();
        }

        /**
         * Increment positions by d. Refer to algorithm paper for the
         * definition of d.
         *
         * @param d The increment value for the position
         * @param start_index start index of the marker array
         * @param end_index end index of the marker array
         */
        private void increment_positions(const int d, const int start_index, const int end_index) 
        {
            for (int i = start_index; i <= end_index; i++) 
            {
                marker_array[i].increment_position(d);
            }
        }

        /**
         * Desired positions incremented by bucket width. The bucket width is
         * basically the desired increments.
         */
        private void update_desired_positions() 
        {
            for (int i{ 1 }; i < marker_array.size(); i++) 
            {
                marker_array[i].update_desired_position();
            }
        }

        /**
         * Sets previous and next markers after default read is done.
         *
         * @param an_input_stream the input stream to be deserialized
         * @Class_Not_Found_Exception thrown when a desired class not found
         * @IOException thrown due to any io errors
         */
        private void read_object(Object_Input_Stream an_input_stream)
                Class_Not_Found_Exception, IOException 
                {
            // always perform the default de-serialization first
            an_input_stream.default_read_object();
            // Build links
            for (int i{ 1 }; i < PSQUARE_CONSTANT; i++) 
            {
                marker_array[i].previous(marker_array[i - 1]).next(marker_array[i + 1]).index(i);
            }
            marker_array[0].previous(marker_array[0]).next(marker_array[1]).index(0);
            marker_array[5].previous(marker_array[4]).next(marker_array[5]).index(5);
        }

        /**
         * Return marker height given index
         *
         * @param marker_index index of marker within (1,6)
         * @return marker height
         */
        //override
        public double height(const int marker_index) 
        {
            Math_Utils::check_range_inclusive(marker_index, 1, marker_array.size() - 1);
            return marker_array[marker_index].marker_height;
        }

        /** {@inherit_doc} */
        //override
        public Markers copy_self() 
        {
            return Markers(new Marker[] 
            {
                Marker(), marker_array[1].copy_self(), marker_array[2].copy_self(), marker_array[3].copy_self(), marker_array[4].copy_self(), marker_array[5].copy_self()
            });

        }

        /**
         * Returns string representation of the Marker array.
         *
         * @return Markers as a string
         */
        //override
        public std::string to_string() const 
        {
            return std::string.format("m1=[%s],m2=[%s],m3=[%s],m4=[%s],m5=[%s]", marker_array[1].to_string(), marker_array[2].to_string(), marker_array[3].to_string(), marker_array[4].to_string(), marker_array[5].to_string());
        }

    }

    /**
     * The class modeling the attributes of the marker of the P-square algorithm
     */
    private static class Marker  
    {

        /**
         * Serial Version ID
         */
        -3575879478288538431L;

        /**
         * The marker index which is just a serial number for the marker in the
         * marker array of 5+1.
         */
        private int index;

        /**
         * The integral marker position. Refer to the variable n in the original
         * works.
         */
        private double int_marker_position;

        /**
         * Desired marker position. Refer to the variable n' in the original
         * works.
         */
        private double desired_marker_position;

        /**
         * Marker height or the quantile. Refer to the variable q in the
         * original works.
         */
        private double marker_height;

        /**
         * Desired marker increment. Refer to the variable dn' in the original
         * works.
         */
        private double desired_marker_increment;

        /**
         * Next and previous markers for easy linked navigation in loops. this
         * is not serialized as they can be rebuilt during deserialization.
         */
        private transient Marker next;

        /**
         * The previous marker links
         */
        private transient Marker previous;

        /**
         * Nonlinear interpolator
         */
        private const Univariate_Interpolator non_linear = Neville_Interpolator();

        /**
         * Linear interpolator which is not serializable
         */
        private transient Univariate_Interpolator linear = Linear_Interpolator();

        /**
         * Default constructor
         */
        private Marker() 
        {
            this.next = this.previous = this;
        }

        /**
         * Constructor of the marker with parameters
         *
         * @param height_of_marker represent the quantile value
         * @param maker_position_desired represent the desired marker position
         * @param marker_position_increment represent increments for position
         * @param marker_position_number represent the position number of marker
         */
        private Marker(double height_of_marker, double maker_position_desired, double marker_position_increment, double marker_position_number) 
        {
            this();
            this.marker_height = height_of_marker;
            this.desired_marker_position = maker_position_desired;
            this.desired_marker_increment = marker_position_increment;
            this.int_marker_position = marker_position_number;
        }

        /**
         * Sets the previous marker.
         *
         * @param previous_marker the previous marker to the current marker in
         *            the array of markers
         * @return this instance
         */
        private Marker previous(const Marker previous_marker) 
        {
            //Math_Utils::check_not_null(previous_marker);
            this.previous = previous_marker;
            return this;
        }

        /**
         * Sets the next marker.
         *
         * @param next_marker the next marker to the current marker in the array
         *            of markers
         * @return this instance
         */
        private Marker next(const Marker next_marker) 
        {
            //Math_Utils::check_not_null(next_marker);
            this.next = next_marker;
            return this;
        }

        /**
         * Sets the index of the marker.
         *
         * @param index_of_marker the array index of the marker in marker array
         * @return this instance
         */
        private Marker index(const int index_of_marker) 
        {
            this.index = index_of_marker;
            return this;
        }

        /**
         * Update desired Position with increment.
         */
        private void update_desired_position() 
        {
            desired_marker_position += desired_marker_increment;
        }

        /**
         * Increment Position by d.
         *
         * @param d a delta value to increment
         */
        private void increment_position(const int d) 
        {
            int_marker_position += d;
        }

        /**
         * Difference between desired and actual position
         *
         * @return difference between desired and actual position
         */
        private double difference() 
        {
            return desired_marker_position - int_marker_position;
        }

        /**
         * Estimate the quantile for the current marker.
         *
         * @return estimated quantile
         */
        private double estimate() 
        {
            const double di = difference();
            const bool is_next_higher =
                    next.int_marker_position - int_marker_position > 1;
            const bool is_previous_lower =
                    previous.int_marker_position - int_marker_position < -1;

            if (di >= 1 && is_next_higher || di <= -1 && is_previous_lower) 
            {
                const int d = di >= 0 ? 1 : -1;
                const std::vector<double> xval = { previous.int_marker_position, int_marker_position, next.int_marker_position };
                const std::vector<double> yval = { previous.marker_height, marker_height, next.marker_height };
                const double xD = int_marker_position + d;

                Univariate_Function univariate_function =
                        non_linear.interpolate(xval, yval);
                marker_height = univariate_function.value(xD);

                // If parabolic estimate is bad then turn linear
                if (is_estimate_bad(yval, marker_height)) 
                {
                    int delta = xD - xval[1] > 0 ? 1 : -1;
                    const std::vector<double> x_bad = { xval[1], xval[1 + delta] };
                    const std::vector<double> y_bad = { yval[1], yval[1 + delta] };
                    Math_Arrays::sort_in_place(x_bad, y_bad);// since d can be +/- 1
                    univariate_function = linear.interpolate(x_bad, y_bad);
                    marker_height = univariate_function.value(xD);
                }
                increment_position(d);
            }
            return marker_height;
        }

        /**
         * Check if parabolic/nonlinear estimate is bad by checking if the
         * ordinate found is beyond the y[0] and y[2].
         *
         * @param y the array to get the bounds
         * @param yD the estimate
         * @return true if yD is a bad estimate
         */
        private bool is_estimate_bad(const std::vector<double> y, const double yD) 
        {
            return yD <= y[0] || yD >= y[2];
        }

        /**
         * {@inherit_doc}<i>This equals method checks for marker attributes and
         * as well checks if navigation pointers (next and previous) are the same
         * between this and passed in object</i>
         *
         * @param o Other object
         * @return true if this equals passed in other object o
         */
        //override
        public bool equals(Object o) 
        {
            bool result = false;
            if (this == o) 
            {
                result = true;
            }
            else if (o instanceof Marker) 
            {
                Marker that = (Marker) o;

                result = Double.compare(marker_height, that.marker_height) == 0;
                result =
                        result &&
                                Double.compare(int_marker_position, that.int_marker_position) == 0;
                result =
                        result &&
                                Double.compare(desired_marker_position, that.desired_marker_position) == 0;
                result =
                        result &&
                                Double.compare(desired_marker_increment, that.desired_marker_increment) == 0;

                result = result && next.index == that.next.index;
                result = result && previous.index == that.previous.index;
            }
            return result;
        }

        /** {@inherit_doc} */
        //override
        public int hash_code() 
        {
            return Arrays.hash_code(std::vector<double> {marker_height, int_marker_position, desired_marker_increment, desired_marker_position, previous.index, next.index});
        }

        /**
         * Read Object to deserialize.
         *
         * @param an_instream Stream Object data
         * @IOException thrown for IO Errors
         * @Class_Not_Found_Exception thrown for class not being found
         */
        private void read_object(Object_Input_Stream an_instream)
                Class_Not_Found_Exception, IOException 
                {
            an_instream.default_read_object();
            previous=next=this;
            linear = Linear_Interpolator();
        }

        /** Copy this instance.
         * @return copy of the instance
         */
        public Marker copy_self() 
        {
            return Marker(marker_height, desired_marker_position, desired_marker_increment, int_marker_position);
        }

        /**
         * {@inherit_doc}
         */
        //override
        public std::string to_string() const 
        {
            return std::string.format(
                    "index=%.0f,n=%.0f,np=%.2f,q=%.2f,dn=%.2f,prev=%d,next=%d", static_cast<double>( index, Precision.round(int_marker_position, 0), Precision.round(desired_marker_position, 2), Precision.round(marker_height, 2), Precision.round(desired_marker_increment, 2), previous.index, next.index);
        }
    }

    /**
     * A simple fixed capacity list that has an upper bound to growth.
     * Once its capacity is reached, {@code add} is a no-op, returning
     * {@code false}.
     *
     * @param <E>
     */
    private static class FixedCapacity_list<E> extends Array_list<E>  
    {

        /**
         * Serialization Version Id
         */
        2283952083075725479L;
        /**
         * Capacity of the list
         */
        private const int capacity;

        /**
         * This constructor constructs the list with given capacity and as well
         * as stores the capacity
         *
         * @param fixed_capacity the capacity to be fixed for this list
         */
        FixedCapacity_list(const int fixed_capacity) 
        {
            super(fixed_capacity);
            this.capacity = fixed_capacity;
        }

        /**
         * {@inherit_doc} In addition it checks if the {@link #size()} returns a
         * size that is within capacity and if true it adds; otherwise the list
         * contents are unchanged and {@code false} is returned.
         *
         * @return true if addition is successful and false otherwise
         */
        //override
        public bool add(const E e) 
        {
            return size() < capacity && super.add(e);
        }

        /**
         * {@inherit_doc} In addition it checks if the sum of Collection size and
         * this instance's {@link #size()} returns a value that is within
         * capacity and if true it adds the collection; otherwise the list
         * contents are unchanged and {@code false} is returned.
         *
         * @return true if addition is successful and false otherwise
         */
        //override
        public bool add_all(Collection<? extends E> collection) 
        {
            bool is_collection_less =
                    collection != NULL &&
                            collection.size() + size() <= capacity;
            return is_collection_less && super.add_all(collection);
        }

        /** {@inherit_doc} */
        //override
        public bool equals(const Object& other) 
        {
            return super.equals(other) && capacity == ((FixedCapacity_list<?>) other).capacity;
        }

        /** {@inherit_doc} */
        //override
        public int hash_code() 
        {
            return super.hash_code() + capacity;
        }

    }

    /**
     * A creation method to build Markers
     *
     * @param initial_five list of initial five elements
     * @param p the quantile desired
     * @return an instance of P_Square_Markers
     */
    public static P_Square_Markers new_markers(const List<Double> initial_five, const double p) 
    {
        return Markers(initial_five, p);
    }

    /**
     * An interface that encapsulates abstractions of the
     * P-square algorithm markers as is explained in the original works. This
     * interface is exposed with protected access to help in testability.
     */
    protected interface P_Square_Markers 
    {
        /**
         * Returns Percentile value computed thus far.
         *
         * @return percentile
         */
        double get_percentile_value();

        /**
         * A deep copy function to clone the current instance.
         *
         * @return deep copy of this instance
         */
        P_Square_Markers copy_self();

        /**
         * Returns the marker height (or percentile) of a given marker index.
         *
         * @param marker_index is the index of marker in the marker array
         * @return percentile value of the marker index passed
         * @ in case the index is not within [1-5]
         */
        double height(const int& marker_index);

        /**
         * Process a data point by moving the marker heights based on estimator.
         *
         * @param input_data_point is the data point passed
         * @return computed percentile
         */
        double process_data_point(double input_data_point);

        /**
         * An Estimate of the percentile value of a given Marker
         *
         * @param index the marker's index in the array of markers
         * @return percentile estimate
         * @ in case if index is not within [1-5]
         */
        double estimate(const int& index);
    }
};