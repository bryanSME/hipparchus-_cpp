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
  //package org.hipparchus.fitting;

  //import java.util.Collection;
#include <vector>
#include "../core/distribution/multivariate/MultivariateNormalDistribution.h"
#include "WeightedObservedPoint.h"
#include "../core/analysis/MultivariateVectorFunction.h"
#include "../optim/nonlinear/vector/leastsquares/LeastSquaresOptimizer.h"
#include "../optim/nonlinear/vector/leastsquares/LevenbergMarquardtOptimizer.h"
#include "../optim/nonlinear/vector/leastsquares/LeastSquaresProblem.h"
#include "../core/analysis/ParametricUnivariateFunction.h"
#include "../core/analysis/MultivariateVectorFunction.h"
#include "../core/analysis/MultivariateMatrixFunction.h"

  //import org.hipparchus.analysis.Multivariate_Matrix_Function;
  //import org.hipparchus.analysis.Multivariate_Vector_function;
  //import org.hipparchus.analysis.Parametric_Univariate_Function ;
  //import org.hipparchus.optim.nonlinear.vector.leastsquares.Least_Squares_Optimizer;
  //import org.hipparchus.optim.nonlinear.vector.leastsquares.Least_Squares_Problem;
  //import org.hipparchus.optim.nonlinear.vector.leastsquares.Levenberg_Marquardt_Optimizer;

  /**
   * Base class that contains common code for fitting parametric univariate
   * real functions <code>y = f(p<sub>i</sub>;x)</code>, where {@code x} is
   * the independent variable and the <code>p<sub>i</sub></code> are the
   * <em>parameters</em>.
   * <br/>
   * A fitter will find the optimal values of the parameters by
   * <em>fitting</em> the curve so it remains very close to a set of
   * {@code N} observed points <code>(x<sub>k</sub>, y<sub>k</sub>)</code>, * {@code 0 <= k < N}.
   * <br/>
   * An algorithm usually performs the fit by finding the parameter
   * values that minimizes the objective function
   * <pre><code>
   *  &sum;y<sub>k</sub> - f(x<sub>k</sub>)<sup>2</sup>, * </code></pre>
   * which is actually a least-squares problem.
   * This class contains boilerplate code for calling the
   * {@link #fit(Collection)} method for obtaining the parameters.
   * The problem setup, such as the choice of optimization algorithm
   * for fitting a specific function is delegated to subclasses.
   *
   */
class Abstract_Curve_Fitter
{
public:
	/**
	 * Fits a curve.
	 * This method computes the coefficients of the curve that best
	 * fit the sample of observed points.
	 *
	 * @param points Observations.
	 * @return the fitted parameters.
	 */
	std::vector<double> fit(const std::vector<Weighted_Observed_Point>& points)
	{
		// Perform the fit.
		return get_optimizer().optimize(get_problem(points)).get_point().to_array();
	}
protected:
	/**
	 * Creates an optimizer set up to fit the appropriate curve.
	 * <p>
	 * The default implementation uses a {@link Levenberg_Marquardt_Optimizer
	 * Levenberg-_Marquardt} optimizer.
	 * </p>
	 * @return the optimizer to use for fitting the curve to the
	 * given {@code points}.
	 */
	Least_Squares_Optimizer get_optimizer()
	{
		return Levenberg_Marquardt_Optimizer();
	}

	/**
	 * Creates a least squares problem corresponding to the appropriate curve.
	 *
	 * @param points Sample points.
	 * @return the least squares problem to use for fitting the curve to the
	 * given {@code points}.
	 */
	virtual Least_Squares_Problem get_problem(const std::vector<Weighted_Observed_Point>& points);

	/**
	 * Vector function for computing function theoretical values.
	 */
	static class Theoretical_Values_Function
	{
	private:
		/** Function to fit. */
		const Parametric_Univariate_Function my_f;
		/** Observations. */
		std::vector<double> my_points;

	public:
		/**
		 * @param f function to fit.
		 * @param observations Observations.
		 */
		Theoretical_Values_Function(const Parametric_Univariate_Function& f, const std::vector<Weighted_Observed_Point>& observations)
		{
			my_f = f;

			const int len = observations.size();
			my_points = std::vector<double>(len);
			int i{};
			for (const auto& obs : observations)
			{
				my_points[i++] = obs.get_x();
			}
		}

		/**
		 * @return the model function values.
		 */
		Multivariate_Vector_function get_model_function()
		{
			return Multivariate_Vector_function()
			{
				/** {@inherit_doc} */
			 //override
				public std::vector<double> value(std::vector<double> p)
				{
					const int len = points.size();
					const std::vector<double>& values = std::vector<double>(len];
					for (int i{}; i < len; i++)
					{
						values[i] = f.value(points[i], p);
					}

					return values;
				}
			};
		}

		/**
		 * @return the model function Jacobian.
		 */
		Multivariate_Matrix_Function get_model_function_jacobian()
		{
			return Multivariate_Matrix_Function()
			{
				/** {@inherit_doc} */
				//override
				public std::vector<std::vector<double>> value(const std::vector<double>& p)
				{
					const int len = points.size();
					auto jacobian = std::vector<std::vector<double>>(len);
					for (int i{}; i < len; i++)
					{
						jacobian[i] = f.gradient(points[i], p);
					}
					return jacobian;
				}
			};
		}
	};
};