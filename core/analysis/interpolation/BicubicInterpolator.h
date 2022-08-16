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
//import org.hipparchus.util.Math_Arrays;
//import org.hipparchus.util.Math_Utils;
#include <vector>
#include "BivariateGridInterpolator.h"
#include "../../util/MathArrays.h"
#include "../../util/MathUtils.h"

/**
 * Generates a {@link Bicubic_Interpolating_Function bicubic interpolating
 * function}.
 * <p>
 *  Caveat: Because the interpolation scheme requires that derivatives be
 *  specified at the sample points, those are approximated with finite
 *  differences (using the 2-points symmetric formulae).
 *  sin_ce their values are undefined at the borders of the provided
 *  interpolation ranges, the interpolated values will be wrong at the
 *  edges of the patch.
 *  The {@code interpolate} method will return a function that overrides
 *  {@link Bicubic_Interpolating_Function#is_valid_point(double,double)} to
 *  indicate points where the interpolation will be inaccurate.
 * </p>
 *
 */
class Bicubic_Interpolator : public Bivariate_Grid_Interpolator
{
public:
    /**
     * {@inherit_doc}
     */
     //override
    Bicubic_Interpolating_Function interpolate(const std::vector<double>& xval, const std::vector<double>& yval, const std::vector<std::vector<double>> fval)
    {
        if (xval.size() == 0 || yval.size() == 0 || fval.size() == 0)
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NO_DATA);
        }
        Math_Utils::check_dimension(xval.size(), fval.size());
        Math_Arrays::check_order(xval);
        Math_Arrays::check_order(yval);

        const int x_len = xval.size();
        const int y_len = yval.size();

        // Approximation to the partial derivatives using finite differences.
        auto dFdX = std::vector<std::vector<double>>(x_len, std::vector<double>(y_len));
        auto d_fd_y = std::vector<std::vector<double>>(x_len, std::vector<double>(y_len));
        auto d2FdXdY = std::vector<std::vector<double>>(x_len, std::vector<double>(y_len));
        for (int i{ 1 }; i < x_len - 1; i++)
        {
            const int nI = i + 1;
            const int pI = i - 1;

            const double nX = xval[nI];
            const double pX = xval[pI];

            const double delta_x = nX - pX;

            for (int j{ 1 }; j < y_len - 1; j++)
            {
                const int& nJ = j + 1;
                const int pJ = j - 1;

                const double nY = yval[nJ];
                const double pY = yval[pJ];

                const double delta_y = nY - pY;

                dFdX[i][j] = (fval[nI][j] - fval[pI][j]) / delta_x;
                d_fd_y[i][j] = (fval[i][nJ] - fval[i][pJ]) / delta_y;

                const double delta_x_y = delta_x * delta_y;

                d2FdXdY[i][j] = (fval[nI][nJ] - fval[nI][pJ] - fval[pI][nJ] + fval[pI][pJ]) / delta_x_y;
            }
        }

        // Create the interpolating function.
        return Bicubic_Interpolating_Function(xval, yval, fval, dFdX, d_fd_y, d2FdXdY)
        {
        public:
            /** {@inherit_doc} */
            //override
            bool is_valid_point(const double& x, const double& y)
            {
                if (x < xval[1] ||
                    x > xval[xval.size() - 2] ||
                    y < yval[1] ||
                    y > yval[yval.size() - 2])
                {
                    return false;
                }
                return true;
            }
        };
    }
};