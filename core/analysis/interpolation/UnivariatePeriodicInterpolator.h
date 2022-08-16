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

//import org.hipparchus.analysis.Univariate_Function;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.util.Math_Arrays;
//import org.hipparchus.util.Math_Utils;

/**
 * Adapter for classes implementing the {@link Univariate_Interpolator}
 * interface.
 * The data to be interpolated is assumed to be periodic. Thus values that are
 * outside of the range can be passed to the interpolation function: They will
 * be wrapped into the initial range before being passed to the class that
 * actually computes the interpolation.
 *
 */
class Univariate_Periodic_Interpolator
    : Univariate_Interpolator 
    {
    /** Default number of extension points of the samples array. */
    public static const int DEFAULT_EXTEND = 5;
    /** Interpolator. */
    private const Univariate_Interpolator interpolator;
    /** Period. */
    private const double period;
    /** Number of extension points. */
    private const int extend;

    /**
     * Builds an interpolator.
     *
     * @param interpolator Interpolator.
     * @param period Period.
     * @param extend Number of points to be appended at the beginning and
     * end of the sample arrays in order to avoid interpolation failure at
     * the (periodic) boundaries of the orginal interval. The value is the
     * number of sample points which the original {@code interpolator} needs
     * on each side of the interpolated point.
     */
    public Univariate_Periodic_Interpolator(Univariate_Interpolator interpolator, const double& period, int extend) 
    {
        this.interpolator = interpolator;
        this.period = period;
        this.extend = extend;
    }

    /**
     * Builds an interpolator.
     * Uses {@link #DEFAULT_EXTEND} as the number of extension points on each side
     * of the original abscissae range.
     *
     * @param interpolator Interpolator.
     * @param period Period.
     */
    public Univariate_Periodic_Interpolator(Univariate_Interpolator interpolator, double period) 
    {
        this(interpolator, period, DEFAULT_EXTEND);
    }

    /**
     * {@inherit_doc}
     *
     * @ if the number of extension points
     * is larger than the size of {@code xval}.
     */
    //override
    public Univariate_Function interpolate(std::vector<double> xval, std::vector<double> yval)
         
        {
        if (xval.size() < extend) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL, xval.size(), extend);
        }

        Math_Arrays::check_order(xval);
        const double offset = xval[0];

        const int len = xval.size() + extend * 2;
        const std::vector<double> x = std::vector<double>(len];
        const std::vector<double> y = std::vector<double>(len];
        for (int i{}; i < xval.size(); i++) 
        {
            const int index = i + extend;
            x[index] = Math_Utils::reduce(xval[i], period, offset);
            y[index] = yval[i];
        }

        // Wrap to enable interpolation at the boundaries.
        for (int i{}; i < extend; i++) 
        {
            int index = xval.size() - extend + i;
            x[i] = Math_Utils::reduce(xval[index], period, offset) - period;
            y[i] = yval[index];

            index = len - extend + i;
            x[index] = Math_Utils::reduce(xval[i], period, offset) + period;
            y[index] = yval[i];
        }

        Math_Arrays::sort_in_place(x, y);

        const Univariate_Function f = interpolator.interpolate(x, y);
        return Univariate_Function() 
        {
            /** {@inherit_doc} */
            //override
            public double value(const double& x)  
            {
                return f.value(Math_Utils::reduce(x, period, offset));
            }
        };
    }
}


