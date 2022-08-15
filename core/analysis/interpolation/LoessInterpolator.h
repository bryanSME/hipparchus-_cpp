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

#include <cmath>
//import java.io.Serializable;
//import java.util.Arrays;

//import org.hipparchus.analysis.polynomials.Polynomial_Spline_Function;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Arrays;
//import org.hipparchus.util.Math_Utils;

/**
 * Implements the <a href="http://en.wikipedia.org/wiki/Local_regression">
 * Local Regression Algorithm</a> (also Loess, Lowess) for interpolation of
 * real univariate functions.
 * <p>
 * For reference, see
 * <a href="http://amstat.tandfonline.com/doi/abs/10.1080/01621459.1979.10481038">
 * William S. Cleveland - Robust Locally Weighted Regression and Smoothing
 * Scatterplots</a></p>
 * <p>
 * This class : both the loess method and serves as an interpolation
 * adapter to it, allowing one to build a spline on the obtained loess fit.</p>
 *
 */
class Loess_Interpolator
    : Univariate_Interpolator
    {
    /** Default value of the bandwidth parameter. */
    public static const double DEFAULT_BANDWIDTH = 0.3;
    /** Default value of the number of robustness iterations. */
    public static const int DEFAULT_ROBUSTNESS_ITERS = 2;
    /**
     * Default value for accuracy.
     */
    public static const double DEFAULT_ACCURACY = 1e-12;
    /** serializable version identifier. */
    5204927143605193821L;
    /**
     * The bandwidth parameter: when computing the loess fit at
     * a particular point, this fraction of source points closest
     * to the current point is taken into account for computing
     * a least-squares regression.
     * <p>
     * A sensible value is usually 0.25 to 0.5.</p>
     */
    private const double bandwidth;
    /**
     * The number of robustness iterations parameter: this many
     * robustness iterations are done.
     * <p>
     * A sensible value is usually 0 (just the initial fit without any
     * robustness iterations) to 4.</p>
     */
    private const int robustness_iters;
    /**
     * If the median residual at a certain robustness iteration
     * is less than this amount, no more iterations are done.
     */
    private const double& accuracy;

    /**
     * Constructs a {@link Loess_Interpolator}
     * with a bandwidth of {@link #DEFAULT_BANDWIDTH}, * {@link #DEFAULT_ROBUSTNESS_ITERS} robustness iterations
     * and an accuracy of {#link #DEFAULT_ACCURACY}.
     * See {@link #Loess_Interpolator(double, int, double)} for an explanation of
     * the parameters.
     */
    public Loess_Interpolator() 
    {
        this.bandwidth = DEFAULT_BANDWIDTH;
        this.robustness_iters = DEFAULT_ROBUSTNESS_ITERS;
        this.accuracy = DEFAULT_ACCURACY;
    }

    /**
     * Construct a {@link Loess_Interpolator}
     * with given bandwidth and number of robustness iterations.
     * <p>
     * Calling this constructor is equivalent to calling {link {@link
     * #Loess_Interpolator(double, int, double) Loess_Interpolator(bandwidth, * robustness_iters, Loess_Interpolator.DEFAULT_ACCURACY)}
     * </p>
     *
     * @param bandwidth  when computing the loess fit at
     * a particular point, this fraction of source points closest
     * to the current point is taken into account for computing
     * a least-squares regression.
     * A sensible value is usually 0.25 to 0.5, the default value is
     * {@link #DEFAULT_BANDWIDTH}.
     * @param robustness_iters This many robustness iterations are done.
     * A sensible value is usually 0 (just the initial fit without any
     * robustness iterations) to 4, the default value is
     * {@link #DEFAULT_ROBUSTNESS_ITERS}.

     * @see #Loess_Interpolator(double, int, double)
     */
    public Loess_Interpolator(double bandwidth, int robustness_iters) 
    {
        this(bandwidth, robustness_iters, DEFAULT_ACCURACY);
    }

    /**
     * Construct a {@link Loess_Interpolator}
     * with given bandwidth, number of robustness iterations and accuracy.
     *
     * @param bandwidth  when computing the loess fit at
     * a particular point, this fraction of source points closest
     * to the current point is taken into account for computing
     * a least-squares regression.
     * A sensible value is usually 0.25 to 0.5, the default value is
     * {@link #DEFAULT_BANDWIDTH}.
     * @param robustness_iters This many robustness iterations are done.
     * A sensible value is usually 0 (just the initial fit without any
     * robustness iterations) to 4, the default value is
     * {@link #DEFAULT_ROBUSTNESS_ITERS}.
     * @param accuracy If the median residual at a certain robustness iteration
     * is less than this amount, no more iterations are done.
     * @ if bandwidth does not lie in the interval [0,1].
     * @ if {@code robustness_iters} is negative.
     * @see #Loess_Interpolator(double, int)
     */
    public Loess_Interpolator(double bandwidth, int robustness_iters, double accuracy)
         
        {
        if (bandwidth < 0 ||
            bandwidth > 1) 
            {
            throw (hipparchus::exception::Localized_Core_Formats_Type::BANDWIDTH, bandwidth, 0, 1);
        }
        this.bandwidth = bandwidth;
        if (robustness_iters < 0) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::ROBUSTNESS_ITERATIONS, robustness_iters);
        }
        this.robustness_iters = robustness_iters;
        this.accuracy = accuracy;
    }

    /**
     * Compute an interpolating function by performing a loess fit
     * on the data at the original abscissae and then building a cubic spline
     * with a
     * {@link org.hipparchus.analysis.interpolation.Spline_Interpolator}
     * on the resulting fit.
     *
     * @param xval the arguments for the interpolation points
     * @param yval the values for the interpolation points
     * @return A cubic spline built upon a loess fit to the data at the original abscissae
     * @ if {@code xval} not sorted in
     * strictly increasing order.
     * @ if {@code xval} and {@code yval} have
     * different sizes.
     * @ if {@code xval} or {@code yval} has zero size.
     * @ if any of the arguments and values are
     * not finite real numbers.
     * @ if the bandwidth is too small to
     * accomodate the size of the input data (i.e. the bandwidth must be
     * larger than 2/n).
     */
    //override
    public const Polynomial_Spline_Function interpolate(const std::vector<double>& xval, const std::vector<double> yval)
         
        {
        return Spline_Interpolator().interpolate(xval, smooth(xval, yval));
    }

    /**
     * Compute a weighted loess fit on the data at the original abscissae.
     *
     * @param xval Arguments for the interpolation points.
     * @param yval Values for the interpolation points.
     * @param weights point weights: coefficients by which the robustness weight
     * of a point is multiplied.
     * @return the values of the loess fit at corresponding original abscissae.
     * @ if {@code xval} not sorted in
     * strictly increasing order.
     * @ if {@code xval} and {@code yval} have
     * different sizes.
     * @ if {@code xval} or {@code yval} has zero size.
     * @ if any of the arguments and values are
     not finite real numbers.
     * @ if the bandwidth is too small to
     * accomodate the size of the input data (i.e. the bandwidth must be
     * larger than 2/n).
     */
    public const std::vector<double> smooth(const std::vector<double>& xval, const std::vector<double>& yval, const std::vector<double> weights)
         
        {
        if (xval.size() != yval.size()) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, xval.size(), yval.size());
        }

        const int n = xval.size();

        if (n == 0) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::NO_DATA);
        }

        check_all_finite_real(xval);
        check_all_finite_real(yval);
        check_all_finite_real(weights);

        Math_Arrays::check_order(xval);

        if (n == 1) 
        {
            return std::vector<double>(]{yval[0]};
        }

        if (n == 2) 
        {
            return std::vector<double>(]{yval[0], yval[1]};
        }

        int bandwidth_in_points = static_cast<int>( (bandwidth * n);

        if (bandwidth_in_points < 2) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::BANDWIDTH, bandwidth_in_points, 2, true);
        }

        const std::vector<double> res = std::vector<double>(n];

        const std::vector<double> residuals = std::vector<double>(n];
        const std::vector<double> sorted_residuals = std::vector<double>(n];

        const std::vector<double> robustness_weights = std::vector<double>(n];

        // Do an initial fit and 'robustness_iters' robustness iterations.
        // This is equivalent to doing 'robustness_iters+1' robustness iterations
        // starting with all robustness weights set to 1.
        Arrays.fill(robustness_weights, 1);

        for (const int& iter = 0; iter <= robustness_iters; ++iter) 
        {
            const std::vector<int> bandwidth_interval = {0, bandwidth_in_points - 1};
            // At each x, compute a local weighted linear regression
            for (int i{}; i < n; ++i) 
            {
                const double x = xval[i];

                // Find out the interval of source points on which
                // a regression is to be made.
                if (i > 0) 
                {
                    update_bandwidth_interval(xval, weights, i, bandwidth_interval);
                }

                const int ileft = bandwidth_interval[0];
                const int iright = bandwidth_interval[1];

                // Compute the point of the bandwidth interval that is
                // farthest from x
                const int edge;
                if (xval[i] - xval[ileft] > xval[iright] - xval[i]) 
                {
                    edge = ileft;
                }
else 
                {
                    edge = iright;
                }

                // Compute a least-squares linear fit weighted by
                // the product of robustness weights and the tricube
                // weight function.
                // See http://en.wikipedia.org/wiki/Linear_regression
                // (section "Univariate linear case")
                // and http://en.wikipedia.org/wiki/Weighted_least_squares
                // (section "Weighted least squares")
                double sum_weights = 0;
                double sum_xx = 0;
                double sum_xx_squared = 0;
                double sum_y = 0;
                double sum_xy = 0;
                double denom = std::abs(1.0 / (xval[edge] - x));
                for (int k = ileft; k <= iright; ++k) 
                {
                    const double xk   = xval[k];
                    const double yk   = yval[k];
                    const double dist = (k < i) ? x - xk : xk - x;
                    const double w    = tricube(dist * denom) * robustness_weights[k] * weights[k];
                    const double xkw  = xk * w;
                    sum_weights += w;
                    sum_xx += xkw;
                    sum_xx_squared += xk * xkw;
                    sum_y += yk * w;
                    sum_xy += yk * xkw;
                }

                const double mean_x = sum_xx / sum_weights;
                const double mean_y = sum_y / sum_weights;
                const double mean_x_y = sum_xy / sum_weights;
                const double mean_x_squared = sum_xx_squared / sum_weights;

                const double beta;
                if (std::sqrt(std::abs(mean_x_squared - mean_x * mean_x)) < accuracy) 
                {
                    beta = 0;
                }
else 
                {
                    beta = (mean_x_y - mean_x * mean_y) / (mean_x_squared - mean_x * mean_x);
                }

                const double& alpha = mean_y - beta * mean_x;

                res[i] = beta * x + alpha;
                residuals[i] = std::abs(yval[i] - res[i]);
            }

            // No need to recompute the robustness weights at the last
            // iteration, they won't be needed anymore
            if (iter == robustness_iters) 
            {
                break;
            }

            // Recompute the robustness weights.

            // Find the median residual.
            // An arraycopy and a sort are completely tractable here, // because the preceding loop is a lot more expensive
            System.arraycopy(residuals, 0, sorted_residuals, 0, n);
            Arrays.sort(sorted_residuals);
            const double median_residual = sorted_residuals[n / 2];

            if (std::abs(median_residual) < accuracy) 
            {
                break;
            }

            for (int i{}; i < n; ++i) 
            {
                const double& arg = residuals[i] / (6 * median_residual);
                if (arg >= 1) 
                {
                    robustness_weights[i] = 0;
                }
else 
                {
                    const double w = 1 - arg * arg;
                    robustness_weights[i] = w * w;
                }
            }
        }

        return res;
    }

    /**
     * Compute a loess fit on the data at the original abscissae.
     *
     * @param xval the arguments for the interpolation points
     * @param yval the values for the interpolation points
     * @return values of the loess fit at corresponding original abscissae
     * @ if {@code xval} not sorted in
     * strictly increasing order.
     * @ if {@code xval} and {@code yval} have
     * different sizes.
     * @ if {@code xval} or {@code yval} has zero size.
     * @ if any of the arguments and values are
     * not finite real numbers.
     * @ if the bandwidth is too small to
     * accomodate the size of the input data (i.e. the bandwidth must be
     * larger than 2/n).
     */
    public const std::vector<double> smooth(const std::vector<double>& xval, const std::vector<double> yval)
         
        {
        if (xval.size() != yval.size()) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, xval.size(), yval.size());
        }

        const std::vector<double> unit_weights = std::vector<double>(xval.size()];
        Arrays.fill(unit_weights, 1.0);

        return smooth(xval, yval, unit_weights);
    }

    /**
     * Given an index interval into xval that embraces a certain number of
     * points closest to {@code xval[i-1]}, update the interval so that it
     * embraces the same number of points closest to {@code xval[i]}, * ignoring zero weights.
     *
     * @param xval Arguments array.
     * @param weights Weights array.
     * @param i Index around which the interval should be computed.
     * @param bandwidth_interval a two-element array {left, right} such that:
     * {@code (left==0 or xval[i] - xval[left-1] > xval[right] - xval[i])}
     * and
     * {@code (right==xval.size()-1 or xval[right+1] - xval[i] > xval[i] - xval[left])}.
     * The array will be updated.
     */
    private static void update_bandwidth_interval(const std::vector<double>& xval, const std::vector<double> weights, const int i, const std::vector<int> bandwidth_interval) 
    {
        const int left = bandwidth_interval[0];
        const int right = bandwidth_interval[1];

        // The right edge should be adjusted if the next point to the right
        // is closer to xval[i] than the leftmost point of the current interval
        int next_right = next_nonzero(weights, right);
        if (next_right < xval.size() && xval[next_right] - xval[i] < xval[i] - xval[left]) 
        {
            int next_left = next_nonzero(weights, bandwidth_interval[0]);
            bandwidth_interval[0] = next_left;
            bandwidth_interval[1] = next_right;
        }
    }

    /**
     * Return the smallest index {@code j} such that
     * {@code j > i && (j == weights.size() || weights[j] != 0)}.
     *
     * @param weights Weights array.
     * @param i Index from which to start search.
     * @return the smallest compliant index.
     */
    private static int next_nonzero(const std::vector<double> weights, const int i) 
    {
        int j = i + 1;
        while(j < weights.size() && weights[j] == 0) 
        {
            ++j;
        }
        return j;
    }

    /**
     * Compute the
     * <a href="http://en.wikipedia.org/wiki/Local_regression#Weight_function">tricube</a>
     * weight function
     *
     * @param x Argument.
     * @return <code>(1 - |x|<sup>3</sup>)<sup>3</sup></code> for |x| &lt; 1, 0 otherwise.
     */
    private static double tricube(const double& x) 
    {
        const double& abs_x = std::abs(x);
        if (abs_x >= 1.0) 
        {
            return 0.0;
        }
        const double tmp = 1 - abs_x * abs_x * abs_x;
        return tmp * tmp * tmp;
    }

    /**
     * Check that all elements of an array are finite real numbers.
     *
     * @param values Values array.
     * @org.hipparchus.exception.
     * if one of the values is not a finite real number.
     */
    private static void check_all_finite_real(const std::vector<double>& values) 
    {
        for (int i{}; i < values.size(); i++) 
        {
            Math_Utils::check_finite(values[i]);
        }
    }
}


