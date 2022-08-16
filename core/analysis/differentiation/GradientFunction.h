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
  //package org.hipparchus.analysis.differentiation;

  //import org.hipparchus.analysis.Multivariate_Vector_function;
#include "../MultivariateVectorFunction.h"
#include "MultivariateDifferentiableFunction.h"
#include "DSFactory.h"

/** Class representing the gradient of a multivariate function.
 * <p>
 * The vectorial components of the function represent the derivatives
 * with respect to each function parameters.
 * </p>
 */
class Gradient_Function : public Multivariate_Vector_function
{
private:
	/** Underlying real-valued function. */
	const Multivariate_Differentiable_Function my_f;

public:
	/** Simple constructor.
	 * @param f underlying real-valued function
	 */
	Gradient_Function(const Multivariate_Differentiable_Function& f)
		: my_f{ f }
	{};

	/** {@inherit_doc} */
	//override
	std::vector<double> value(const std::vector<double>& point)
	{
		// set up parameters
		DS_Factory factory = DS_Factory(point.size(), 1);
		auto ds_x = std::vector<Derivative_Structure>(point.size());
		for (int i{}; i < point.size(); ++i)
		{
			ds_x[i] = factory.variable(i, point[i]);
		}

		// compute the derivatives
		const Derivative_Structure ds_y = my_f.value(ds_x);

		// extract the gradient
		auto y = std::vector<double>(point.size());
		auto orders = std::vector<int>(point.size());
		for (int i{}; i < point.size(); ++i)
		{
			orders[i] = 1;
			y[i] = ds_y.get_partial_derivative(orders);
			orders[i] = 0;
		}

		return y;
	}
};