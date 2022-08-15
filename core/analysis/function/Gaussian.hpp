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

#include <exception>
#include <numbers>
#include <cmath>
#include "../differentiation/UnivariateDifferentiableFunction.h"
#include "../ParametricUnivariateFunction.h"
#include "../../util/MathUtils.h"

//import java.util.Arrays;

//import org.hipparchus.analysis.Parametric_Univariate_Function ;
//import org.hipparchus.analysis.differentiation.Derivative;
//import org.hipparchus.analysis.differentiation.Univariate_Differentiable_Function;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Null_Argument_Exception;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Utils;
//import org.hipparchus.util.Precision;

/**
 * <a href="http://en.wikipedia.org/wiki/Gaussian_function">
 *  Gaussian</a> function.
 *
 */
class Gaussian : public Univariate_Differentiable_Function 
{
private:
    /** Mean. */
    const double my_mean;
    /** Inverse of the standard deviation. */
    double my_is;
    /** Inverse of twice the square of the standard deviation. */
    double my_i2s2;
    /** Normalization factor. */
    const double my_norm;

    /**
     * @param x_minus_mean {@code x - mean}.
     * @param norm Normalization factor.
     * @param i2s2 Inverse of twice the square of the standard deviation.
     * @return the value of the Gaussian at {@code x}.
     */
    static double value(const double& x_minus_mean, const double& norm, const double& i2s2)
    {
        return norm * std::exp(-x_minus_mean * x_minus_mean * i2s2);
    }

public:
    /**
     * Gaussian with given normalization factor, mean and standard deviation.
     *
     * @param norm Normalization factor.
     * @param mean Mean.
     * @param sigma Standard deviation.
     * @ if {@code sigma <= 0}.
     */
    Gaussian(const double& norm, const double& mean, const double& sigma)
        :
        my_norm{ norm }, my_mean{ mean }
    {
        if (sigma <= 0) 
        {
            throw std::exception("(hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL_BOUND_EXCLUDED, sigma, 0)");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL_BOUND_EXCLUDED, sigma, 0);
        }

        my_is   = 1 / sigma;
        my_i2s2 = 0.5 * my_is * my_is;
    }

    /**
     * Normalized gaussian with given mean and standard deviation.
     *
     * @param mean Mean.
     * @param sigma Standard deviation.
     * @ if {@code sigma <= 0}.
     */
    Gaussian(const double& mean, const double& sigma)
    {
        this(1 / (sigma * std::sqrt(2 * std::numbers::pi)), mean, sigma);
    }

    /**
     * Normalized gaussian with zero mean and unit standard deviation.
     */
    Gaussian() 
    {
        this(0, 1);
    }

    /** {@inherit_doc} */
    //override
    double value(double x) 
    {
        return value(x - my_mean, my_norm, my_i2s2);
    }

    /** {@inherit_doc}
     */
     //override
    template<typename T, typename std::enable_if<std::is_base_of<Derivative<T>, T>::value>::type* = nullptr>
    T value(T t)
    {
        const double u = is * (t.get_value() - mean);
        auto f = std::vector<double>(t.get_order() + 1];

        // the nth order derivative of the Gaussian has the form:
        // dn(g(x)/dxn = (norm / s^n) P_n(u) exp(-u^2/2) with u=(x-m)/s
        // where P_n(u) is a degree n polynomial with same parity as n
        // P_0(u) = 1, P_1(u) = -u, P_2(u) = u^2 - 1, P_3(u) = -u^3 + 3 u...
        // the general recurrence relation for P_n is:
        // P_n(u) = P_(n-1)'(u) - u P_(n-1)(u)
        // as per polynomial parity, we can store coefficients of both P_(n-1) and P_n in the same array
        const std::vector<double> p = std::vector<double>(f.size()];
        p[0] = 1;
        const double u2 = u * u;
        double coeff = norm * std::exp(-0.5 * u2);
        if (coeff <= Precision.SAFE_MIN)
        {
            Arrays.fill(f, 0.0);
        }
        else
        {
            f[0] = coeff;
            for (int n{ 1 }; n < f.size(); ++n)
            {
                // update and evaluate polynomial P_n(x)
                double v{};
                p[n] = -p[n - 1];
                for (int k{n}; k >= 0; k -= 2)
                {
                    v = v * u2 + p[k];
                    if (k > 2)
                    {
                        p[k - 2] = (k - 1) * p[k - 1] - p[k - 3];
                    }
                    else if (k == 2)
                    {
                        p[0] = p[1];
                    }
                }
                if ((n & 0x1) == 1)
                {
                    v *= u;
                }

                coeff *= is;
                f[n] = coeff * v;

            }
        }

        return t.compose(f);

    }

    /**
     * Parametric function where the input array contains the parameters of
     * the Gaussian, ordered as follows:
     * <ul>
     *  <li>Norm</li>
     *  <li>Mean</li>
     *  <li>Standard deviation</li>
     * </ul>
     */
    static class Parametric : Parametric_Univariate_Function  
    {
    public:
        /**
         * Computes the value of the Gaussian at {@code x}.
         *
         * @param x Value for which the function must be computed.
         * @param param Values of norm, mean and standard deviation.
         * @return the value of the function.
         * @Null_Argument_Exception if {@code param} is {@code NULL}.
         * @ if the size of {@code param} is
         * not 3.
         * @ if {@code param[2]} is negative.
         */
        //override
        double value(const double& x, double ... param)
        {
            validate_parameters(param);

            const double diff = x - param[1];
            const double i2s2 = 1 / (2 * param[2] * param[2]);
            return Gaussian::value(diff, param[0], i2s2);
        }

        /**
         * Computes the value of the gradient at {@code x}.
         * The components of the gradient vector are the partial
         * derivatives of the function with respect to each of the
         * <em>parameters</em> (norm, mean and standard deviation).
         *
         * @param x Value at which the gradient must be computed.
         * @param param Values of norm, mean and standard deviation.
         * @return the gradient vector at {@code x}.
         * @Null_Argument_Exception if {@code param} is {@code NULL}.
         * @ if the size of {@code param} is
         * not 3.
         * @ if {@code param[2]} is negative.
         */
        //override
        std::vector<double> gradient(const double& x, double ... param)
        {
            validate_parameters(param);

            const double norm = param[0];
            const double diff = x - param[1];
            const double sigma = param[2];
            const double i2s2 = 1 / (2 * sigma * sigma);

            const double n = Gaussian::value(diff, 1, i2s2);
            const double m = norm * n * 2 * i2s2 * diff;
            const double s = m * diff / sigma;

            return std::vector<double> { n, m, s };
        }
    private:
        /**
         * Validates parameters to ensure they are appropriate for the evaluation of
         * the {@link #value(double,std::vector<double>)} and {@link #gradient(double,std::vector<double>)}
         * methods.
         *
         * @param param Values of norm, mean and standard deviation.
         * @Null_Argument_Exception if {@code param} is {@code NULL}.
         * @ if the size of {@code param} is
         * not 3.
         * @ if {@code param[2]} is negative.
         */
        void validate_parameters(std::vector<double> param)
        {
            if (param == NULL) 
            {
                throw std::exception("not implemented");
                //throw Null_Argument_Exception();
            }
            Math_Utils::check_dimension(param.size(), 3);
            if (param[2] <= 0) 
            {
                throw std::exception("not implemented");
                //throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL_BOUND_EXCLUDED, param[2], 0);
            }
        }
    };
};