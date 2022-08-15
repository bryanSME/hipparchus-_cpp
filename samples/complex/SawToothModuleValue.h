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

/** Domain coloring enhancing modules changes.
 * <p>
 * Value represents module but uses a sawtooth function.
 * </p>
 * <p>
 * The sawtooth function is computed from a base 2 logarithm and fractional parts.
 * It enhances modules changes with discontinuities and dark rings.
 * </p>
 */
class Saw_Tooth_Module_Value extends Domain_coloring 
{

    /** Simple constructor.
     * @param saturation constant saturation
     */
    protected Saw_Tooth_Module_Value(const double saturation) 
    {
        super(saturation);
    }

    /** {@inherit_doc} */
    //override
    public double value(const std::complex<double> z) 
    {
        const double module = z.norm();
        return std::log(module) / std::log(2.0);
    }

}


