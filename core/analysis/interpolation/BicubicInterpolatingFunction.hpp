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

#include "../BivariateFunction.h"
#include <vector>
#include "../../util/MathArrays.h"
#include "../../util/MathUtils.h"

/**
 * Function that : the
 * <a href="http://en.wikipedia.org/wiki/Bicubic_interpolation">
 * bicubic spline interpolation</a>.
 *
 */
class Bicubic_Interpolating_Function
    : Bivariate_Function 
{
private:
    /** Number of coefficients. */
    static constexpr int NUM_COEFF{ 16 };
    /**
     * Matrix to compute the spline coefficients from the function values
     * and function derivatives values
     */
    static const std::vector<std::vector<double>> AINV = {
        { 1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
        { 0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0 }, 
        { -3,3,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0 }, 
        { 2,-2,0,0,1,1,0,0,0,0,0,0,0,0,0,0 },
        { 0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0 }, 
        { 0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0 }, 
        { 0,0,0,0,0,0,0,0,-3,3,0,0,-2,-1,0,0 },
        { 0,0,0,0,0,0,0,0,2,-2,0,0,1,1,0,0 }, 
        { -3,0,3,0,0,0,0,0,-2,0,-1,0,0,0,0,0 }, 
        { 0,0,0,0,-3,0,3,0,0,0,0,0,-2,0,-1,0 }, 
        { 9,-9,-9,9,6,3,-6,-3,6,-6,3,-3,4,2,2,1 },
        { -6,6,6,-6,-3,-3,3,3,-4,4,-2,2,-2,-2,-1,-1 },
        { 2,0,-2,0,0,0,0,0,1,0,1,0,0,0,0,0 },
        { 0,0,0,0,2,0,-2,0,0,0,0,0,1,0,1,0 },
        { -6,6,6,-6,-4,-2,4,2,-3,3,-3,3,-2,-1,-2,-1 },
        { 4,-4,-4,4,2,2,-2,-2,2,-2,2,-2,1,1,1,1 }
    };

    /** Samples x-coordinates */
    const std::vector<double> my_xval;
    /** Samples y-coordinates */
    const std::vector<double> my_yval;
    /** Set of cubic splines patching the whole data grid */
    const std::vector<std::vector<Bicubic_Function>> my_splines;

public:
    /**
     * @param x Sample values of the x-coordinate, in increasing order.
     * @param y Sample values of the y-coordinate, in increasing order.
     * @param f Values of the function on every grid point.
     * @param dFdX Values of the partial derivative of function with respect
     * to x on every grid point.
     * @param d_fd_y Values of the partial derivative of function with respect
     * to y on every grid point.
     * @param d2FdXdY Values of the cross partial derivative of function on
     * every grid point.
     * @ if the various arrays do not contain
     * the expected number of elements.
     * @ if {@code x} or {@code y} are
     * not strictly increasing.
     * @ if any of the arrays has zero length.
     */
    Bicubic_Interpolating_Function(const std::vector<double>& x, const std::vector<double>& y, const std::vector<std::vector<double>>& f, const std::vector<std::vector<double>>& dFdX, const std::vector<std::vector<double>>& d_fd_y, const std::vector<std::vector<double>>& d2FdXdY)
    {
        const int x_len = x.size();
        const int y_len = y.size();

        if (x_len == 0 || y_len == 0 || f.size() == 0 || f[0].size() == 0) 
        {
            throw std::exception("not implmented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NO_DATA);
        }
        Math_Utils::check_dimension(x_len, f.size());
        Math_Utils::check_dimension(x_len, dFdX.size());
        Math_Utils::check_dimension(x_len, d_fd_y.size());
        Math_Utils::check_dimension(x_len, d2FdXdY.size());
        Math_Arrays::check_order(x);
        Math_Arrays::check_order(y);

        xval = x;
        yval = y;

        const int last_i = x_len - 1;
        const int last_j = y_len - 1;
        splines = Bicubic_Function[last_i][last_j];

        for (int i{}; i < last_i; i++) 
        {
            Math_Utils::check_dimension(f[i].size(), y_len);
            Math_Utils::check_dimension(dFdX[i].size(), y_len);
            Math_Utils::check_dimension(d_fd_y[i].size(), y_len);
            Math_Utils::check_dimension(d2FdXdY[i].size(), y_len);

            const int ip1 = i + 1;
            const double x_r = xval[ip1] - xval[i];
            for (int j{}; j < last_j; j++) 
            {
                const int jp1 = j + 1;
                const double y_r = yval[jp1] - yval[j];
                const double x_ry_r = x_r * y_r;
                const auto beta = std::vector<double>
                {
                    f[i][j], f[ip1][j], f[i][jp1], f[ip1][jp1], dFdX[i][j] * x_r, dFdX[ip1][j] * x_r, dFdX[i][jp1] * x_r, dFdX[ip1][jp1] * x_r, d_fd_y[i][j] * y_r, d_fd_y[ip1][j] * y_r, d_fd_y[i][jp1] * y_r, d_fd_y[ip1][jp1] * y_r, d2FdXdY[i][j] * x_ry_r, d2FdXdY[ip1][j] * x_ry_r, d2FdXdY[i][jp1] * x_ry_r, d2FdXdY[ip1][jp1] * x_ry_r
                };

                splines[i][j] = Bicubic_Function(compute_spline_coefficients(beta));
            }
        }
    }

    /**
     * {@inherit_doc}
     */
    //override
    public double value(const double& x, double y)
         
        {
        const int i = search_index(x, xval);
        const int j = search_index(y, yval);

        const double xN = (x - xval[i]) / (xval[i + 1] - xval[i]);
        const double yN = (y - yval[j]) / (yval[j + 1] - yval[j]);

        return splines[i][j].value(xN, yN);
    }

    /**
     * Indicates whether a point is within the interpolation range.
     *
     * @param x First coordinate.
     * @param y Second coordinate.
     * @return {@code true} if (x, y) is a valid point.
     */
    public bool is_valid_point(const double& x, double y) 
    {
        if (x < xval[0] ||
            x > xval[xval.size() - 1] ||
            y < yval[0] ||
            y > yval[yval.size() - 1]) 
            {
            return false;
        }
else 
        {
            return true;
        }
    }

    /**
     * @param c Coordinate.
     * @param val Coordinate samples.
     * @return the index in {@code val} corresponding to the interval
     * containing {@code c}.
     * @ if {@code c} is out of the
     * range defined by the boundary values of {@code val}.
     */
    private int search_index(cosnt double& c, const std::vector<double>& val) 
    {
        const int r = Arrays.binary_search(val, c);

        if (r == -1 ||
            r == -val.size() - 1) 
            {
            throw std::exception("not implmented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::OUT_OF_RANGE_SIMPLE, c, val[0], val[val.size() - 1]);
        }

        if (r < 0) 
        {
            // "c" in within an interpolation sub-interval: Return the
            // index of the sample at the lower end of the sub-interval.
            return -r - 2;
        }
        const int last = val.size() - 1;
        if (r == last) 
        {
            // "c" is the last sample of the range: Return the index
            // of the sample at the lower end of the last sub-interval.
            return last - 1;
        }

        // "c" is another sample point.
        return r;
    }

    /**
     * Compute the spline coefficients from the list of function values and
     * function partial derivatives values at the four corners of a grid
     * element. They must be specified in the following order:
     * <ul>
     *  <li>f(0,0)</li>
     *  <li>f(1,0)</li>
     *  <li>f(0,1)</li>
     *  <li>f(1,1)</li>
     *  <li>f<sub>x</sub>(0,0)</li>
     *  <li>f<sub>x</sub>(1,0)</li>
     *  <li>f<sub>x</sub>(0,1)</li>
     *  <li>f<sub>x</sub>(1,1)</li>
     *  <li>f<sub>y</sub>(0,0)</li>
     *  <li>f<sub>y</sub>(1,0)</li>
     *  <li>f<sub>y</sub>(0,1)</li>
     *  <li>f<sub>y</sub>(1,1)</li>
     *  <li>f<sub>xy</sub>(0,0)</li>
     *  <li>f<sub>xy</sub>(1,0)</li>
     *  <li>f<sub>xy</sub>(0,1)</li>
     *  <li>f<sub>xy</sub>(1,1)</li>
     * </ul>
     * where the subscripts indicate the partial derivative with respect to
     * the corresponding variable(s).
     *
     * @param beta List of function values and function partial derivatives
     * values.
     * @return the spline coefficients.
     */
    private std::vector<double> compute_spline_coefficients(const std::vector<double>& beta) 
    {
        const std::vector<double> a = std::vector<double>(NUM_COEFF];

        for (int i{}; i < NUM_COEFF; i++) 
        {
            double result{};
            const std::vector<double> row = AINV[i];
            for (int j{}; j < NUM_COEFF; j++) 
            {
                result += row[j] * beta[j];
            }
            a[i] = result;
        }

        return a;
    }
}

/**
 * Bicubic function.
 */
class Bicubic_Function : Bivariate_Function 
{
    /** Number of points. */
    private static const short N = 4;
    /** Coefficients */
    private const std::vector<std::vector<double>> a;

    /**
     * Simple constructor.
     *
     * @param coeff Spline coefficients.
     */
    Bicubic_Function(std::vector<double> coeff) 
    {
        a = std::vector<double>(N][N];
        for (int j{}; j < N; j++) 
        {
            const std::vector<double> aJ = a[j];
            for (int i{}; i < N; i++) 
            {
                aJ[i] = coeff[i * N + j];
            }
        }
    }

    /**
     * {@inherit_doc}
     */
    //override
    public double value(const double& x, double y) 
    {
        Math_Utils::check_range_inclusive(x, 0, 1);
        Math_Utils::check_range_inclusive(y, 0, 1);

        const double x2 = x * x;
        const double x3 = x2 * x;
        const std::vector<double> pX = {1, x, x2, x3};

        const double y2 = y * y;
        const double y3 = y2 * y;
        const std::vector<double> pY = {1, y, y2, y3};

        return apply(pX, pY, a);
    }

    /**
     * Compute the value of the bicubic polynomial.
     *
     * @param pX Powers of the x-coordinate.
     * @param pY Powers of the y-coordinate.
     * @param coeff Spline coefficients.
     * @return the interpolated value.
     */
    private double apply(std::vector<double> pX, std::vector<double> pY, std::vector<std::vector<double>> coeff) 
    {
        double result{};
        for (int i{}; i < N; i++) 
        {
            const double r = Math_Arrays::linear_combination(coeff[i], pY);
            result += r * pX[i];
        }

        return result;
    }
}


