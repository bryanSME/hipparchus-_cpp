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

#include  <vector>
#include "../BivariateFunction.h"

  /**
   * Interface representing a bivariate real interpolating function where the
   * sample points must be specified on a regular grid.
   *
   */
class Bivariate_Grid_Interpolator
{
	/**
	 * Compute an interpolating function for the dataset.
	 *
	 * @param xval All the x-coordinates of the interpolation points, sorted
	 * in increasing order.
	 * @param yval All the y-coordinates of the interpolation points, sorted
	 * in increasing order.
	 * @param fval The values of the interpolation points on all the grid knots:
	 * {@code fval[i][j] = f(xval[i], yval[j])}.
	 * @return a function which interpolates the dataset.
	 * @ if any of the arrays has zero length.
	 * @ if the array lengths are inconsistent.
	 * @ if the array is not sorted.
	 * @ if the number of points is too small for
	 * the order of the interpolation
	 */
	virtual Bivariate_Function interpolate(const std::vector<double>& xval, const std::vector<double>& yval, const std::vector<std::vector<double>>& fval);
};