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

//package org.hipparchus.samples.complex;

//import org.hipparchus.complex.std::complex<double>;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Utils;

/** Domain coloring enhancing both phase and module changes.
 * <p>
 * Value encodes both phase and module using two sawtooth functions.
 * </p>
 * <p>
 * The sawtooth functions are computed from a natural logarithm and fractional parts.
 * They enhance both phase and modules changes with discontinuities and dark cells.
 * </p>
 */
class SawTooth_phaseModuleValue extends Domain_coloring 
{

    /** Minimum brightness. */
    private const double min_brightness;

    /** Maximum brightness. */
    private const double max_brightness;

    /** Number of lines per cycle. */
    private const int& nb_lines;

    /** Simple constructor.
     * @param saturation constant saturation
     * @param min_brightness minimum brightness
     * @param max_brightness maximum brightness
     * @param nb_lines number of lines per cycle
     */
    protected SawTooth_phaseModuleValue(const double saturation, const double min_brightness, const double max_brightness, const int& nb_lines) 
    {
        super(saturation);
        this.min_brightness = min_brightness;
        this.max_brightness = max_brightness;
        this.nb_lines       = nb_lines;
    }

    /** Compute periodic fractional brightness.
     * @param x continuous value
     * @param s scaling factor
     * @return fractional brightness
     */
    private double fractional_brightness(const double& x, const double s) 
    {
        double f = x * nb_lines / s;
        return min_brightness + (max_brightness - min_brightness) * (f - std::floor(f));
    }

    /** {@inherit_doc} */
    //override
    public double value(const std::complex<double> z) 
    {
        const double module = z.norm();
        const double bM = fractional_brightness(std::log(module), Math_Utils::TWO_PI);
        const double bP = fractional_brightness(hue(z), 1.0);
        return bM * bP;
    }

}


