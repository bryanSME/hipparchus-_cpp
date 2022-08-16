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

#include <vector>
#include "../../../util/MathArrays.h"
#include "../../UnivariateFunction.h"

  /**
   * Class that : the Gaussian rule for
   * {@link #integrate(Univariate_Function) integrating} a weighted
   * function.
   *
   */
class Gauss_Integrator
{
private:
	/** Nodes. */
	const std::vector<double> my_points;
	/** Nodes weights. */
	const std::vector<double> my_weights;

public:
	/**
	 * Creates an integrator from the given {@code points} and {@code weights}.
	 * The integration interval is defined by the first and last value of
	 * {@code points} which must be sorted in increasing order.
	 *
	 * @param points Integration points.
	 * @param weights Weights of the corresponding integration nodes.
	 * @ if the {@code points} are not
	 * sorted in increasing order.
	 * @ if points and weights don't have the same length
	 */
	Gauss_Integrator(const std::vector<double>& points, const std::vector<double>& weights)
		: my_points{ points }, my_weights{ weights }
	{
		Math_Arrays::check_equal_length(points, weights);
		//Math_Arrays::check_order(points, Math_Arrays::Order_Direction::INCREASING, true, true);
	}

	/**
	 * Creates an integrator from the given pair of points (first element of
	 * the pair) and weights (second element of the pair.
	 *
	 * @param points_and_weights Integration points and corresponding weights.
	 * @ if the {@code points} are not
	 * sorted in increasing order.
	 *
	 * @see #Gauss_Integrator(std::vector<double>, std::vector<double>)
	 */
	Gauss_Integrator(const std::pair<std::vector<double>, std::vector<double>>& points_and_weights)
	{
		Gauss_Integrator(points_and_weights.first, points_and_weights.second);
	}

	/**
	 * Returns an estimate of the integral of {@code f(x) * w(x)}, * where {@code w} is a weight function that depends on the actual
	 * flavor of the Gauss integration scheme.
	 * The algorithm uses the points and associated weights, as passed
	 * to the {@link #Gauss_Integrator(std::vector<double>,std::vector<double>) constructor}.
	 *
	 * @param f Function to integrate.
	 * @return the integral of the weighted function.
	 */
	double integrate(const Univariate_Function& f)
	{
		double s{};
		double c{};
		for (int i{}; i < my_points.size(); i++)
		{
			const double x = my_points[i];
			const double w = my_weights[i];
			const double y = w * f.value(x) - c;
			const double t = s + y;
			c = (t - s) - y;
			s = t;
		}
		return s;
	}

	/**
	 * @return the order of the integration rule (the number of integration
	 * points).
	 */
	int get_number_of_points() const
	{
		return my_points.size();
	}

	/**
	 * Gets the integration point at the given index.
	 * The index must be in the valid range but no check is performed.
	 * @param index index of the integration point
	 * @return the integration point.
	 */
	double get_point(const int& index) const
	{
		return my_points[index];
	}

	/**
	 * Gets the weight of the integration point at the given index.
	 * The index must be in the valid range but no check is performed.
	 * @param index index of the integration point
	 * @return the weight.
	 */
	double get_weight(const int& index) const
	{
		return my_weights[index];
	}
};