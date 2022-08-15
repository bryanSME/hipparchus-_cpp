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
//package org.hipparchus.analysis.interpolation;

//import org.hipparchus.exception.;

/**
 * Interpolate grid data using bi-linear interpolation.
 * @since 1.4
 */
class Bilinear_Interpolator : Bivariate_Grid_Interpolator 
{

    /** {@inherit_doc} */
    //override
    public Bilinear_interpolating_function interpolate(const std::vector<double>& xval, const std::vector<double>& yval, const std::vector<std::vector<double>> fval)
         
        {
        return Bilinear_interpolating_function(xval, yval, fval);
    }

}


