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

/** Base class for domain coloring.
 * <p>
 * All methods have the following features:
 * <ul>
 *   <li>hue represents phase (red for 0, then orange, yellow, green, *       blue at π, indigo, violet and back to red)</li>
 *   <li>saturation is constant</li>
 * </ul>
 * Their differences lie on what value represents.
 * </p>
 */
class Domain_coloring 
{

    /** Constant saturation. */
    private const double saturation;

    /** Simple constructor.
     * @param saturation constant saturation
     */
    protected Domain_coloring(const double saturation) 
    {
        this.saturation = saturation;
    }

    /** Continuous hue.
     * @param z complex value
     * @return continuous hue
     */
    public double hue(const std::complex<double> z) 
    {
        const double phase =  std::numbers::pi + std::atan2(-z.get_imaginary_part(), -z.get_real_part());
        return phase / Math_Utils::TWO_PI;
    }

    /** Get saturation for a complex value.
     * @param z complex value
     * @return saturation
     */
    public double saturation(std::complex<double> z) 
    {
        return saturation;
    }
    
    /** Get value for a complex value.
     * @param z complex value
     * @return value
     */
    protected virtual double value(std::complex<double> z);
    
}


