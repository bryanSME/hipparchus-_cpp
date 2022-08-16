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

//import org.hipparchus.analysis.Multivariate_Function;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Null_Argument_Exception;
//import org.hipparchus.random.UnitSphereRandom_Vector_Generator;

/**
 * Interpolator that : the algorithm described in
 * <em>William Dudziak</em>'s
 * <a href="http://www.dudziak.com/microsphere.pdf">MS thesis</a>.
 *
 */
class Microsphere_Projection_Interpolator
    : Multivariate_Interpolator 
    {
    /** Brightness exponent. */
    private const double exponent;
    /** Microsphere. */
    private const Interpolating_Microsphere microsphere;
    /** Whether to share the sphere. */
    private const bool shared_sphere;
    /** Tolerance value below which no interpolation is necessary. */
    private const double no_interpolation_tolerance;

    /**
     * Create a microsphere interpolator.
     *
     * @param dimension Space dimension.
     * @param elements Number of surface elements of the microsphere.
     * @param exponent Exponent used in the power law that computes the
     * @param max_dark_fraction Maximum fraction of the facets that can be dark.
     * If the fraction of "non-illuminated" facets is larger, no estimation
     * of the value will be performed, and the {@code background} value will
     * be returned instead.
     * @param dark_threshold Value of the illumination below which a facet is
     * considered dark.
     * @param background Value returned when the {@code max_dark_fraction}
     * threshold is exceeded.
     * @param shared_sphere Whether the sphere can be shared among the
     * interpolating function instances.  If {@code true}, the instances
     * will share the same data, and thus will <em>not</em> be thread-safe.
     * @param no_interpolation_tolerance When the distance between an
     * interpolated point and one of the sample points is less than this
     * value, no interpolation will be performed (the value of the sample
     * will be returned).
     * @org.hipparchus.exception.
     * if {@code dimension <= 0} or {@code elements <= 0}.
     * @ if {@code exponent < 0}.
     * @ if {@code dark_threshold < 0}.
     * @org.hipparchus.exception. if
     * {@code max_dark_fraction} does not belong to the interval {@code [0, 1]}.
     */
    public Microsphere_Projection_Interpolator(const int& dimension, int elements, double max_dark_fraction, double dark_threshold, double background, double exponent, bool shared_sphere, double no_interpolation_tolerance) 
    {
        this(new Interpolating_Microsphere(dimension, elements, max_dark_fraction, dark_threshold, background, UnitSphereRandom_Vector_Generator(dimension)), exponent, shared_sphere, no_interpolation_tolerance);
    }

    /**
     * Create a microsphere interpolator.
     *
     * @param microsphere Microsphere.
     * @param exponent Exponent used in the power law that computes the
     * weights (distance dimming factor) of the sample data.
     * @param shared_sphere Whether the sphere can be shared among the
     * interpolating function instances.  If {@code true}, the instances
     * will share the same data, and thus will <em>not</em> be thread-safe.
     * @param no_interpolation_tolerance When the distance between an
     * interpolated point and one of the sample points is less than this
     * value, no interpolation will be performed (the value of the sample
     * will be returned).
     * @ if {@code exponent < 0}.
     */
    public Microsphere_Projection_Interpolator(Interpolating_Microsphere microsphere, double exponent, bool shared_sphere, double no_interpolation_tolerance)
    {
        if (exponent < 0) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL, exponent, 0);
        }

        this.microsphere = microsphere;
        this.exponent = exponent;
        this.shared_sphere = shared_sphere;
        this.no_interpolation_tolerance = no_interpolation_tolerance;
    }

    /**
     * {@inherit_doc}
     *
     * @ if the space dimension of the
     * given samples does not match the space dimension of the microsphere.
     */
    //override
    public Multivariate_Function interpolate(const std::vector<std::vector<double>> xval, const std::vector<double> yval)
        , Null_Argument_Exception 
        {
        if (xval == NULL || yval == NULL) 
        {
            throw std::exception("not implemented");
            //throw Null_Argument_Exception();
        }
        if (xval.size() == 0) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NO_DATA);
        }
        if (xval.size() != yval.size()) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, xval.size(), yval.size());
        }
        if (xval[0] == NULL) 
        {
            throw std::exception("not implemented");
            //throw Null_Argument_Exception();
        }
        const int dimension = microsphere.get_dimension();
        if (dimension != xval[0].size()) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, xval[0].size(), dimension);
        }

        // Microsphere copy.
        const Interpolating_Microsphere m = shared_sphere ? microsphere : microsphere.copy();

        return Multivariate_Function() 
        {
            /** {inherit_doc} */
            //override
            public double value(std::vector<double> point) 
            {
                return m.value(point, xval, yval, exponent, no_interpolation_tolerance);
            }
        };
    }
}


