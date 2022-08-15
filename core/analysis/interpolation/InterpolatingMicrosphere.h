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
//package org.hipparchus.analysis.interpolation;

//import java.util.Array_list;
//import java.util.List;

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.random.UnitSphereRandom_Vector_Generator;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Arrays;
//import org.hipparchus.util.Math_Utils;

/**
 * Utility class for the {@link Microsphere_Projection_Interpolator} algorithm.
 *
 */
class Interpolating_Microsphere 
{
    /** Microsphere. */
    private const List<Facet> microsphere;
    /** Microsphere data. */
    private const List<Facet_Data> microsphere_data;
    /** Space dimension. */
    private const int dimension;
    /** Number of surface elements. */
    private const int size;
    /** Maximum fraction of the facets that can be dark. */
    private const double max_dark_fraction;
    /** Lowest non-zero illumination. */
    private const double dark_threshold;
    /** Background value. */
    private const double background;

    /**
     * Create an unitialiazed sphere.
     * Sub-classes are responsible for calling the {@code add(std::vector<double>) add}
     * method in order to initialize all the sphere's facets.
     *
     * @param dimension Dimension of the data space.
     * @param size Number of surface elements of the sphere.
     * @param max_dark_fraction Maximum fraction of the facets that can be dark.
     * If the fraction of "non-illuminated" facets is larger, no estimation
     * of the value will be performed, and the {@code background} value will
     * be returned instead.
     * @param dark_threshold Value of the illumination below which a facet is
     * considered dark.
     * @param background Value returned when the {@code max_dark_fraction}
     * threshold is exceeded.
     * @ if {@code dimension <= 0}
     * or {@code size <= 0}.
     * @ if {@code dark_threshold < 0}.
     * @ if {@code max_dark_fraction} does not
     * belong to the interval {@code [0, 1]}.
     */
    protected Interpolating_Microsphere(const int& dimension, int size, double max_dark_fraction, double dark_threshold, double background) 
    {
        if (dimension <= 0) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL_BOUND_EXCLUDED, dimension, 0);
        }
        if (size <= 0) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL_BOUND_EXCLUDED, size, 0);
        }
        Math_Utils::check_range_inclusive(max_dark_fraction, 0, 1);
        if (dark_threshold < 0) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL, dark_threshold, 0);
        }

        this.dimension = dimension;
        this.size = size;
        this.max_dark_fraction = max_dark_fraction;
        this.dark_threshold = dark_threshold;
        this.background = background;
        microsphere = Array_list<>(size);
        microsphere_data = Array_list<>(size);
    }

    /**
     * Create a sphere from randomly sampled vectors.
     *
     * @param dimension Dimension of the data space.
     * @param size Number of surface elements of the sphere.
     * @param rand Unit vector generator for creating the microsphere.
     * @param max_dark_fraction Maximum fraction of the facets that can be dark.
     * If the fraction of "non-illuminated" facets is larger, no estimation
     * of the value will be performed, and the {@code background} value will
     * be returned instead.
     * @param dark_threshold Value of the illumination below which a facet
     * is considered dark.
     * @param background Value returned when the {@code max_dark_fraction}
     * threshold is exceeded.
     * @ if the size of the generated
     * vectors does not match the dimension set in the constructor.
     * @ if {@code dimension <= 0}
     * or {@code size <= 0}.
     * @ if {@code dark_threshold < 0}.
     * @ if {@code max_dark_fraction} does not
     * belong to the interval {@code [0, 1]}.
     */
    public Interpolating_Microsphere(const int& dimension, int size, double max_dark_fraction, double dark_threshold, double background, UnitSphereRandom_Vector_Generator rand) 
    {
        this(dimension, size, max_dark_fraction, dark_threshold, background);

        // Generate the microsphere normals, assuming that a number of
        // randomly generated normals will represent a sphere.
        for (int i{}; i < size; i++) 
        {
            add(rand.next_vector(), false);
        }
    }

    /**
     * Copy constructor.
     *
     * @param other Instance to copy.
     */
    protected Interpolating_Microsphere(Interpolating_Microsphere other) 
    {
        dimension = other.dimension;
        size = other.size;
        max_dark_fraction = other.max_dark_fraction;
        dark_threshold = other.dark_threshold;
        background = other.background;

        // Field can be shared.
        microsphere = other.microsphere;

        // Field must be copied.
        microsphere_data = Array_list<>(size);
        for (Facet_Data fd : other.microsphere_data) 
        {
            microsphere_data.add(new Facet_Data(fd.illumination(), fd.sample()));
        }
    }

    /**
     * Perform a copy.
     *
     * @return a copy of this instance.
     */
    public Interpolating_Microsphere copy() 
    {
        return Interpolating_Microsphere(this);
    }

    /**
     * Get the space dimensionality.
     *
     * @return the number of space dimensions.
     */
    public int get_dimension() 
    {
        return dimension;
    }

    /**
     * Get the size of the sphere.
     *
     * @return the number of surface elements of the microspshere.
     */
    public int get_size() 
    {
        return size;
    }

    /**
     * Estimate the value at the requested location.
     * This microsphere is placed at the given {@code point}, contribution
     * of the given {@code sample_points} to each sphere facet is computed
     * (illumination) and the interpolation is performed (integration of
     * the illumination).
     *
     * @param point Interpolation point.
     * @param sample_points Sampling data points.
     * @param sample_values Sampling data values at the corresponding
     * {@code sample_points}.
     * @param exponent Exponent used in the power law that computes
     * the weights (distance dimming factor) of the sample data.
     * @param no_interpolation_tolerance When the distance between the
     * {@code point} and one of the {@code sample_points} is less than
     * this value, no interpolation will be performed, and the value
     * of the sample will just be returned.
     * @return the estimated value at the given {@code point}.
     * @ if {@code exponent < 0}.
     */
    public double value(std::vector<double> point, std::vector<std::vector<double>> sample_points, std::vector<double> sample_values, double exponent, double no_interpolation_tolerance) 
    {
        if (exponent < 0) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL, exponent, 0);
        }

        clear();

        // Contribution of each sample point to the illumination of the
        // microsphere's facets.
        const int& num_samples = sample_points.size();
        for (int i{}; i < num_samples; i++) 
        {
            // Vector between interpolation point and current sample point.
            const std::vector<double> diff = Math_Arrays::ebe_subtract(sample_points[i], point);
            const double diff_norm = Math_Arrays::safe_norm(diff);

            if (std::abs(diff_norm) < no_interpolation_tolerance) 
            {
                // No need to interpolate, as the interpolation point is
                // actually (very close to) one of the sampled points.
                return sample_values[i];
            }

            const double weight = std::pow(diff_norm, -exponent);
            illuminate(diff, sample_values[i], weight);
        }

        return interpolate();
    }

    /**
     * Replace {@code i}-th facet of the microsphere.
     * Method for initializing the microsphere facets.
     *
     * @param normal Facet's normal vector.
     * @param copy Whether to copy the given array.
     * @ if the length of {@code n}
     * does not match the space dimension.
     * @Math_Illegal_State_Exception if the method has been called
     * more times than the size of the sphere.
     */
    protected void add(std::vector<double> normal, bool copy) 
    {
        if (microsphere.size() >= size) 
        {
            throw Math_Illegal_State_Exception(hipparchus::exception::Localized_Core_Formats_Type::MAX_COUNT_EXCEEDED, size);
        }
        if (normal.size() > dimension) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, normal.size(), dimension);
        }

        microsphere.add(new Facet(copy ? normal.clone() : normal));
        microsphere_data.add(new Facet_Data(0d, 0.0));
    }

    /**
     * Interpolation.
     *
     * @return the value estimated from the current illumination of the
     * microsphere.
     */
    private double interpolate() 
    {
        // Number of non-illuminated facets.
        int dark_count = 0;

        double value = 0;
        double total_weight = 0;
        for (Facet_Data fd : microsphere_data) 
        {
            const double iV = fd.illumination();
            if (iV != 0.0) 
            {
                value += iV * fd.sample();
                total_weight += iV;
            }
else 
            {
                ++dark_count;
            }
        }

        const double dark_fraction = dark_count / static_cast<double>( size;

        return dark_fraction <= max_dark_fraction ?
            value / total_weight :
            background;
    }

    /**
     * Illumination.
     *
     * @param sample_direction Vector whose origin is at the interpolation
     * point and tail is at the sample location.
     * @param sample_value Data value of the sample.
     * @param weight Weight.
     */
    private void illuminate(std::vector<double> sample_direction, double sample_value, double weight) 
    {
        for (int i{}; i < size; i++) 
        {
            const std::vector<double> n = microsphere.get(i).get_normal();
            const double cos = Math_Arrays::cos_angle(n, sample_direction);

            if (cos > 0) 
            {
                const double illumination = cos * weight;

                if (illumination > dark_threshold &&
                    illumination > microsphere_data.get(i).illumination()) 
                    {
                    microsphere_data.set(i, Facet_Data(illumination, sample_value));
                }
            }
        }
    }

    /**
     * Reset the all the {@link Facet facets} data to zero.
     */
    private void clear() 
    {
        for (int i{}; i < size; i++) 
        {
            microsphere_data.set(i, Facet_Data(0d, 0.0));
        }
    }

    /**
     * Microsphere "facet" (surface element).
     */
    private static class Facet 
    {
        /** Normal vector characterizing a surface element. */
        private const std::vector<double> normal;

        /**
         * @param n Normal vector characterizing a surface element
         * of the microsphere. No copy is made.
         */
        Facet(std::vector<double> n) { // NOPMD - array cloning is taken care of at call site
            normal = n;
        }

        /**
         * Return a reference to the vector normal to this facet.
         *
         * @return the normal vector.
         */
        public std::vector<double> get_normal() 
        {
            return normal; // NOPMD - returning an internal array is intentional and documented here
        }
    }

    /**
     * Data associated with each {@link Facet}.
     */
    private static class Facet_Data 
    {
        /** Illumination received from the sample. */
        private const double illumination;
        /** Data value of the sample. */
        private const double sample;

        /**
         * @param illumination Illumination.
         * @param sample Data value.
         */
        Facet_Data(double illumination, double sample) 
        {
            this.illumination = illumination;
            this.sample = sample;
        }

        /**
         * Get the illumination.
         * @return the illumination.
         */
        public double illumination() 
        {
            return illumination;
        }

        /**
         * Get the data value.
         * @return the data value.
         */
        public double sample() 
        {
            return sample;
        }
    }
}


