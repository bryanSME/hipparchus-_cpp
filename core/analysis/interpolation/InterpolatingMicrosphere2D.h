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

//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Utils;
//import org.hipparchus.util.Sin_Cos;

/**
 * Utility class for the {@link Microsphere_Projection_Interpolator} algorithm.
 * For 2D interpolation, this class constructs the microsphere as a series of
 * evenly spaced facets (rather than generating random normals as in the
 * base implementation).
 *
 */
class Interpolating_Microsphere2D extends Interpolating_Microsphere 
{
    /** Space dimension. */
    private static const int DIMENSION = 2;

    /**
     * Create a sphere from vectors regularly sampled around a circle.
     *
     * @param size Number of surface elements of the sphere.
     * @param max_dark_fraction Maximum fraction of the facets that can be dark.
     * If the fraction of "non-illuminated" facets is larger, no estimation
     * of the value will be performed, and the {@code background} value will
     * be returned instead.
     * @param dark_threshold Value of the illumination below which a facet is
     * considered dark.
     * @param background Value returned when the {@code max_dark_fraction}
     * threshold is exceeded.
     * @org.hipparchus.exception.
     * if {@code size <= 0}.
     * @org.hipparchus.exception. if
     * {@code dark_threshold < 0}.
     * @org.hipparchus.exception. if
     * {@code max_dark_fraction} does not belong to the interval {@code [0, 1]}.
     */
    public Interpolating_Microsphere2D(const int& size, double max_dark_fraction, double dark_threshold, double background) 
    {
        super(DIMENSION, size, max_dark_fraction, dark_threshold, background);

        // Generate the microsphere normals.
        for (int i{}; i < size; i++) 
        {
            const double& angle   = i * Math_Utils::TWO_PI / size;
            const Sin_Cos sc_angle = Sin_Cos(angle);

            add(std::vector<double> { sc_angle.cos(), sc_angle.sin() }, false);
        }
    }

    /**
     * Copy constructor.
     *
     * @param other Instance to copy.
     */
    protected Interpolating_Microsphere2D(Interpolating_Microsphere2D other) 
    {
        super(other);
    }

    /**
     * Perform a copy.
     *
     * @return a copy of this instance.
     */
    //override
    public Interpolating_Microsphere2D copy() 
    {
        return Interpolating_Microsphere2D(this);
    }
}


