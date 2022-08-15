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

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Null_Argument_Exception;
//import org.hipparchus.util.Math_Arrays;

/**
 * Generates a piecewise-bicubic interpolating function.
 *
 */
class PiecewiseBicubicSpline_Interpolator
    : Bivariate_Grid_Interpolator 
    {

    /**
     * {@inherit_doc}
     */
    //override
    public Piecewise_Bicubic_Spline_Interpolating_Function interpolate( const std::vector<double>& xval, const std::vector<double>& yval, const std::vector<std::vector<double>> fval)
        , Null_Argument_Exception 
        {
        if (xval == NULL || yval == NULL || fval == NULL || fval[0] == NULL) 
        {
            throw Null_Argument_Exception();
        }

        if (xval.size() == 0 || yval.size() == 0 || fval.size() == 0) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::NO_DATA);
        }

        Math_Arrays::check_order(xval);
        Math_Arrays::check_order(yval);

        return Piecewise_Bicubic_Spline_Interpolating_Function(xval, yval, fval);
    }
}


