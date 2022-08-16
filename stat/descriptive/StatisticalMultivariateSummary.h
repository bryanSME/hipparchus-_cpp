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
  //package org.hipparchus.stat.descriptive;
#include <vector>
#include "../../core/linear/RealMatrix.h"
//import org.hipparchus.linear.Real_Matrix;

/**
 * Reporting interface for basic multivariate statistics.
 */
class Statistical_Multivariate_Summary
{
	/**
	 * Returns the dimension of the data
	 * @return The dimension of the data
	 */
	virtual int get_dimension() = 0;

	/**
	 * Returns an array whose i<sup>th</sup> entry is the
	 * mean of the i<sup>th</sup> entries of the arrays
	 * that correspond to each multivariate sample
	 *
	 * @return the array of component means
	 */
	virtual std::vector<double> get_mean() = 0;

	/**
	 * Returns the covariance of the available values.
	 * @return The covariance, NULL if no multivariate sample
	 * have been added or a zeroed matrix for a single value set.
	 */
	virtual Real_Matrix get_covariance() = 0;

	/**
	 * Returns an array whose i<sup>th</sup> entry is the
	 * standard deviation of the i<sup>th</sup> entries of the arrays
	 * that correspond to each multivariate sample
	 *
	 * @return the array of component standard deviations
	 */
	virtual std::vector<double> get_standard_deviation() = 0;

	/**
	 * Returns an array whose i<sup>th</sup> entry is the
	 * maximum of the i<sup>th</sup> entries of the arrays
	 * that correspond to each multivariate sample
	 *
	 * @return the array of component maxima
	 */
	virtual std::vector<double> get_max() = 0;

	/**
	 * Returns an array whose i<sup>th</sup> entry is the
	 * minimum of the i<sup>th</sup> entries of the arrays
	 * that correspond to each multivariate sample
	 *
	 * @return the array of component minima
	 */
	virtual std::vector<double> get_min() = 0;

	/**
	 * Returns the number of available values
	 * @return The number of available values
	 */
	virtual long get_n() = 0;

	/**
	 * Returns an array whose i<sup>th</sup> entry is the
	 * geometric mean of the i<sup>th</sup> entries of the arrays
	 * that correspond to each multivariate sample
	 *
	 * @return the array of component geometric means
	 */
	virtual std::vector<double> get_geometric_mean() = 0;

	/**
	 * Returns an array whose i<sup>th</sup> entry is the
	 * sum of the i<sup>th</sup> entries of the arrays
	 * that correspond to each multivariate sample
	 *
	 * @return the array of component sums
	 */
	virtual std::vector<double> get_sum() = 0;

	/**
	 * Returns an array whose i<sup>th</sup> entry is the
	 * sum of squares of the i<sup>th</sup> entries of the arrays
	 * that correspond to each multivariate sample
	 *
	 * @return the array of component sums of squares
	 */
	virtual std::vector<double> get_sum_sq() = 0;

	/**
	 * Returns an array whose i<sup>th</sup> entry is the
	 * sum of logs of the i<sup>th</sup> entries of the arrays
	 * that correspond to each multivariate sample
	 *
	 * @return the array of component log sums
	 */
	virtual std::vector<double> get_sum_log() = 0;
};