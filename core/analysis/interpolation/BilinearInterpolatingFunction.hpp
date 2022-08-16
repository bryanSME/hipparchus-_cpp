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

 //import java.io.Serializable;

 //import org.hipparchus.Calculus_Field_Element;
 //import org.hipparchus.analysis.Bivariate_Function;
 //import org.hipparchus.analysis.Field_Bivariate_Function;
 //import org.hipparchus.exception.;
#include <type_traits>
#include <vector>
#include "../../CalculusFieldElement.hpp"
#include "GridAxis.h"
#include "../BivariateFunction.h"
#include "../FieldBivariateFunction.hpp"

/**
 * Interpolate grid data using bi-linear interpolation.
 * <p>
 * This interpolator is thread-safe.
 * </p>
 * @since 1.4
 */
class Bilinear_interpolating_function : public Bivariate_Function, public Field_Bivariate_Function
{
private:

	/** Grid along the x axis. */
	const Grid_Axis my_x_grid;

	/** Grid along the y axis. */
	const Grid_Axis my_y_grid;

	/** Grid size along the y axis. */
	const int my_y_size;

	/** Values of the interpolation points on all the grid knots (in a flatten array). */
	std::vector<double> my_f_val;

public:
	/** Simple constructor.
	 * @param x_val All the x-coordinates of the interpolation points, sorted
	 * in increasing order.
	 * @param y_val All the y-coordinates of the interpolation points, sorted
	 * in increasing order.
	 * @param f_val The values of the interpolation points on all the grid knots:
	 * {@code f_val[i][j] = f(x_val[i], y_val[j])}.
	 * @exception  if grid size is smaller than 2
	 * or if the grid is not sorted in strict increasing order
	 */
	Bilinear_interpolating_function(const std::vector<double>& x_val, const std::vector<double>& y_val, const std::vector<std::vector<double>>& f_val)
		:
		my_x_grid{ Grid_Axis(x_val, 2) },
		my_y_grid{ Grid_Axis(y_val, 2) },
		my_y_size{ y_val.size() }
		my_f_val{ std::vector<double>(x_val.size() * y_val.size()) }
	{
		int k{};
		for (int i{}; i < x_val.size(); ++i)
		{
			const std::vector<double> fi = f_val[i];
			for (int j{}; j < y_size; ++j)
			{
				my_f_val[k++] = fi[j];
			}
		}
	}

	/** Get the lowest grid x coordinate.
	 * @return lowest grid x coordinate
	 */
	double get_x_inf()
	{
		return x_grid.node(0);
	}

	/** Get the highest grid x coordinate.
	 * @return highest grid x coordinate
	 */
	double get_x_sup()
	{
		return x_grid.node(x_grid.size() - 1);
	}

	/** Get the lowest grid y coordinate.
	 * @return lowest grid y coordinate
	 */
	double get_y_inf()
	{
		return y_grid.node(0);
	}

	/** Get the highest grid y coordinate.
	 * @return highest grid y coordinate
	 */
	double get_y_sup()
	{
		return y_grid.node(y_grid.size() - 1);
	}

	/** {@inherit_doc} */
	//override
	double value(const double& x, const double& y)
	{
		// get the interpolation nodes
		const int    i = x_grid.interpolation_index(x);
		const int    j = y_grid.interpolation_index(y);
		const double x0 = x_grid.node(i);
		const double x1 = x_grid.node(i + 1);
		const double y0 = y_grid.node(j);
		const double y1 = y_grid.node(j + 1);

		// get the function values at interpolation nodes
		const int    k0 = i * y_size + j;
		const int    k1 = k0 + y_size;
		const double z00 = f_val[k0];
		const double z01 = f_val[k0 + 1];
		const double z10 = f_val[k1];
		const double z11 = f_val[k1 + 1];

		// interpolate
		const double dx0 = x - x0;
		const double dx1 = x1 - x;
		const double dx10 = x1 - x0;
		const double dy0 = y - y0;
		const double dy1 = y1 - y;
		const double dy10 = y1 - y0;
		return (dx0 * (dy0 * z11 + dy1 * z10) + dx1 * (dy0 * z01 + dy1 * z00)) /
			(dx10 * dy10);
	}

	/** {@inherit_doc}
	 * @since 1.5
	 */
	 //override
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	T value(const T& x, const T& y)
	{
		// get the interpolation nodes
		const int    i = x_grid.interpolation_index(x.get_real());
		const int    j = y_grid.interpolation_index(y.get_real());
		const double x0 = x_grid.node(i);
		const double x1 = x_grid.node(i + 1);
		const double y0 = y_grid.node(j);
		const double y1 = y_grid.node(j + 1);

		// get the function values at interpolation nodes
		const int    k0 = i * y_size + j;
		const int    k1 = k0 + y_size;
		const double z00 = f_val[k0];
		const double z01 = f_val[k0 + 1];
		const double z10 = f_val[k1];
		const double z11 = f_val[k1 + 1];

		// interpolate
		const T      dx0 = x.subtract(x0);
		const T      mdx1 = x.subtract(x1);
		const double dx10 = x1 - x0;
		const T      dy0 = y.subtract(y0);
		const T      mdy1 = y.subtract(y1);
		const double dy10 = y1 - y0;
		return dy0.multiply(z11).subtract(mdy1.multiply(z10)).multiply(dx0).
			subtract(dy0.multiply(z01).subtract(mdy1.multiply(z00)).multiply(mdx1)).
			divide(dx10 * dy10);
	}
};