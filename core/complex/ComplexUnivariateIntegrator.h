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

 //import java.util.function.Double_Function;

#include "../analysis/integration/UnivariateIntegrator.h"
#include "../analysis/CalculusFieldUnivariateFunction.h"
#include <complex>
#include <vector>

/**
 * Wrapper to perform univariate complex integration using an underlying real integration algorithms.
 * @since 2.0
 */
class std::complex<double>_Univariate_Integrator
{
private:
	/** Underlying real integrator. */
	Univariate_Integrator my_integrator;

public:
	/** Crate a complex integrator from a real integrator.
	 * @param integrator underlying real integrator to use
	 */
	std::complex<double>_Univariate_Integrator(const Univariate_Integrator integrator) : my_integrator{ integrator } {};

	/**
	 * Integrate a function along a straight path between points.
	 *
	 * @param max_eval maximum number of evaluations (real and imaginary
	 * parts are evaluated separately, so up to twice this number may be used)
	 * @param f the integrand function
	 * @param start start point of the integration path
	 * @param end end point of the integration path
	 * @return the value of integral along the straight path
	 */
	std::complex<double> integrate(const int& max_eval, const Calculus_Field_Univariate_Function<std::complex<double>>& f, const std::complex<double>& start, const std::complex<double>& end)
	{
		throw std::exception("std::complex<double>_Univariate_Integrator -- not fully implemented");
		//// linear mapping from real interval [0; 1] to function value along complex straight path from start to end
		//const std::complex<double> rate = end - start;
		//const Double_Function<std::complex> mapped = t->f.value(start + rate * t);

		//// integrate real and imaginary parts separately
		//const double real = integrator.integrate(max_eval, t->mapped.apply(t).get_real_part(), 0.0, 1.0);
		//const double imaginary = integrator.integrate(max_eval, t->mapped.apply(t).get_imaginary_part(), 0.0, 1.0);

		//// combine integrals
		//return std::complex<double>(real, imaginary) * rate;
	}

	/**
	 * Integrate a function along a polyline path between any number of points.
	 *
	 * @param max_eval maximum number of evaluations (real and imaginary
	 * parts are evaluated separately and each path segments are also evaluated
	 * separately, so up to 2n times this number may be used for n segments)
	 * @param f the integrand function
	 * @param start start point of the integration path
	 * @param path successive points defining the path vertices
	 * @return the value of integral along the polyline path
	 */
	std::complex<double> integrate(const int& max_eval, const Calculus_Field_Univariate_Function<std::complex<double>>& f, const std::complex<double>& start, const std::vector<std::complex<double>>& path)
	{
		std::complex<double> sum = std::complex<double>(0); // .ZERO;
		std::complex<double> previous = start;
		for (const auto& current : path)
		{
			sum = sum + integrate(max_eval, f, previous, current);
			previous = current;
		}
		return sum;
	}
};