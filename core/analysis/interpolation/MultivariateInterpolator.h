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
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Null_Argument_Exception;

/**
 * Interface representing a univariate real interpolating function.
 *
 */
class Multivariate_Interpolator 
{

    /**
     * Computes an interpolating function for the data set.
     *
     * @param xval the arguments for the interpolation points.
     * {@code xval[i][0]} is the first component of interpolation point
     * {@code i}, {@code xval[i][1]} is the second component, and so on
     * until {@code xval[i][d-1]}, the last component of that interpolation
     * point (where {@code d} is thus the dimension of the space).
     * @param yval the values for the interpolation points
     * @return a function which interpolates the data set
     * @ if the arguments violate assumptions
     * made by the interpolation algorithm.
     * @ when the array dimensions are not consistent.
     * @ if an array has zero-length.
     * @Null_Argument_Exception if the arguments are {@code NULL}.
     */
    Multivariate_Function interpolate(std::vector<std::vector<double>> xval, std::vector<double> yval)
        , Null_Argument_Exception;
}


