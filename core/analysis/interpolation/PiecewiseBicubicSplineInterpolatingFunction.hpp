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

//import java.util.Arrays;

//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.analysis.Bivariate_Function;
//import org.hipparchus.analysis.Field_Bivariate_Function;
//import org.hipparchus.analysis.polynomials.Field_Polynomial_Spline_Function;
//import org.hipparchus.analysis.polynomials.Polynomial_Spline_Function;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Null_Argument_Exception;
//import org.hipparchus.util.Math_Arrays;
#include <type_traits>
#include <vector>
#include "../../analysis/BivariateFunction.h"
#include "../../analysis/FieldBivariateFunction.h"
#include "../../exception/LocalizedCoreFormats.h"
#include "../../util/MathArrays.h"
#include "AkimaSplineInterpolator.h"
#include "../../CalculusFieldElement.hpp"
#include "../polynomials/PolynomialSplineFunction.h"

/**
 * Function that : the
 * <a href="http://www.paulinternet.nl/?page=bicubic">bicubic spline</a>
 * interpolation.
 * This implementation currently uses {@link Akima_Spline_Interpolator} as the
 * underlying one-dimensional interpolator, which requires 5 sample points;
 * insufficient data will raise an exception when the
 * {@link #value(double,double) value} method is called.
 *
 */
class Piecewise_Bicubic_Spline_Interpolating_Function : public Bivariate_Function, public Field_Bivariate_Function
{
private:
    /** The minimum number of points that are needed to compute the function. */
    static constexpr int MIN_NUM_POINTS{ 5 };
    /** Samples x-coordinates */
    std::vector<double> my_xval;
    /** Samples y-coordinates */
    std::vector<double> my_yval;
    /** Set of cubic splines patching the whole data grid */
    std::vector<std::vector<double>> my_fval;

    /**
 * @param c Coordinate.
 * @param val Coordinate samples.
 * @param offset how far back from found value to offset for querying
 * @param count total number of elements forward from beginning that will be
 *        queried
 * @return the index in {@code val} corresponding to the interval containing
 *         {@code c}.
 * @ if {@code c} is out of the range defined by
 *         the boundary values of {@code val}.
 */
    private int search_index(const double& c, const std::vector<double>& val, const int& offset, const int& count)
    {
        int r = Arrays.binary_search(val, c);

        if (r == -1 || r == -val.size() - 1)
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::OUT_OF_RANGE_SIMPLE, c, val[0], val[val.size() - 1]);
        }

        if (r < 0)
        {
            // "c" in within an interpolation sub-interval, which returns
            // negative
            // need to remove the negative sign for consistency
            r = -r - offset - 1;
        }
        else
        {
            r -= offset;
        }

        if (r < 0)
        {
            r = 0;
        }

        if ((r + count) >= val.size())
        {
            // "c" is the last sample of the range: Return the index
            // of the sample at the lower end of the last sub-interval.
            r = val.size() - count;
        }

        return r;
    }

public:
    /**
     * @param x Sample values of the x-coordinate, in increasing order.
     * @param y Sample values of the y-coordinate, in increasing order.
     * @param f Values of the function on every grid point. the expected number
     *        of elements.
     * @ if {@code x} or {@code y} are not
     *         strictly increasing.
     * @Null_Argument_Exception if any of the arguments are NULL
     * @ if any of the arrays has zero length.
     * @ if the length of x and y don't match the row, column
     *         height of f
     */
    Piecewise_Bicubic_Spline_Interpolating_Function(const std::vector<double>& x, const std::vector<double>& y, const std::vector<std::vector<double>>& f)
    {
        if (x == NULL ||
            y == NULL ||
            f == NULL ||
            f[0] == NULL)
        {
            throw std::exception("not implemented");
            //throw Null_Argument_Exception();
        }

        const int x_len = x.size();
        const int y_len = y.size();

        if (x_len == 0 ||
            y_len == 0 ||
            f.size() == 0 ||
            f[0].size() == 0)
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::NO_DATA);
        }

        if (x_len < MIN_NUM_POINTS ||
            y_len < MIN_NUM_POINTS ||
            f.size() < MIN_NUM_POINTS ||
            f[0].size() < MIN_NUM_POINTS)
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::INSUFFICIENT_DATA);
        }

        if (x_len != f.size())
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, x_len, f.size());
        }

        if (y_len != f[0].size())
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, y_len, f[0].size());
        }

        Math_Arrays::check_order(x);
        Math_Arrays::check_order(y);

        my_xval = x;
        my_yval = y;
        my_fval = f;
    }

    /**
     * {@inherit_doc}
     */
     //override
    double value(const double& x, const double& y)
    {
        const auto interpolator = Akima_Spline_Interpolator();
        const int offset{ 2 };
        const int count{ offset + 3 };
        const int i = search_index(x, my_xval, offset, count);
        const int j = search_index(y, my_yval, offset, count);

        auto x_array = std::vector<double>(count);
        auto y_array = std::vector<double>(count);
        auto z_array = std::vector<double>(count);
        auto interp_array = std::vector<double>(count);

        for (int index = 0; index < count; index++)
        {
            x_array[index] = my_xval[i + index];
            y_array[index] = my_yval[j + index];
        }

        for (int z_index{}; z_index < count; z_index++)
        {
            for (int index = 0; index < count; index++)
            {
                z_array[index] = fval[i + index][j + z_index];
            }
            const auto spline = interpolator.interpolate(x_array, z_array);
            interp_array[z_index] = spline.value(x);
        }

        const auto spline = interpolator.interpolate(y_array, interp_array);

        return spline.value(y);

    }

    /**
     * {@inherit_doc}
     * @since 1.5
     */
     //override
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    T value(const T& x, const T& y)

    {
        const auto interpolator = Akima_Spline_Interpolator();
        const int offset{ 2 };
        const int count{ offset + 3 };
        const int i = search_index(x.get_real(), my_xval, offset, count);
        const int j = search_index(y.get_real(), my_yval, offset, count);

        const auto x_array = std::vector<double>(count);
        const std::vector<T> y_array = Math_Arrays::build_array(x.get_field(), count);
        const auto z_array = std::vector<double>(count);
        const std::vector<T> interp_array = Math_Arrays::build_array(x.get_field(), count);

        const T zero = x.get_field().get_zero();
        for (int index{}; index < count; index++)
        {
            x_array[index] = xval[i + index];
            y_array[index] = zero.add(yval[j + index]);
        }

        for (int z_index{}; z_index < count; z_index++)
        {
            for (int index{}; index < count; index++)
            {
                z_array[index] = fval[i + index][j + z_index];
            }
            const Polynomial_Spline_Function spline = interpolator.interpolate(x_array, z_array);
            interp_array[z_index] = spline.value(x);
        }

        const Field_Polynomial_Spline_Function<T> spline = interpolator.interpolate(y_array, interp_array);

        return spline.value(y);

    }

    /**
     * Indicates whether a point is within the interpolation range.
     *
     * @param x First coordinate.
     * @param y Second coordinate.
     * @return {@code true} if (x, y) is a valid point.
     */
    bool is_valid_point(const double& x, const double& y)
    {
        return !(
            x < my_xval[0] ||
            x > my_xval[my_xval.size() - 1] ||
            y < my_yval[0] ||
            y > my_yval[my_yval.size() - 1]
        );
    }
};