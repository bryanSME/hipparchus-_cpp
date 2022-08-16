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
#include "../differentiation/UnivariateDifferentiableFunction.h"
#include "../ParametricUnivariateFunction.h"
#include "../../util/MathUtils.h"

/**
 * <a href="http://en.wikipedia.org/wiki/Harmonic_oscillator">
 *  simple harmonic oscillator</a> function.
 *
 */
class Harmonic_Oscillator : public Univariate_Differentiable_Function 
{
private:
    /** Amplitude. */
    const double my_amplitude;
    /** Angular frequency. */
    const double my_omega;
    /** Phase. */
    const double my_phase;

    /**
     * @param x_times_omega_plus_phase {@code omega * x + phase}.
     * @param amplitude Amplitude.
     * @return the value of the harmonic oscillator function at {@code x}.
     */
    static double value(const double& x_times_omega_plus_phase, const double& amplitude)
    {
        return amplitude * std::cos(x_times_omega_plus_phase);
    }

public:
    /**
     * Harmonic oscillator function.
     *
     * @param amplitude Amplitude.
     * @param omega Angular frequency.
     * @param phase Phase.
     */
    Harmonic_Oscillator(double amplitude, double omega, double phase) : my_amplitude{ amplitude }, my_omega{ omega }, my_phase{ phase } {};

    /** {@inherit_doc} */
    //override
    double value(double x) const
    {
        return value(my_omega * x + my_phase, my_amplitude);
    }

    /**
     * Parametric function where the input array contains the parameters of
     * the harmonic oscillator function, ordered as follows:
     * <ul>
     *  <li>Amplitude</li>
     *  <li>Angular frequency</li>
     *  <li>Phase</li>
     * </ul>
     */
    class Parametric : Parametric_Univariate_Function
    {
    public:
        /**
         * Computes the value of the harmonic oscillator at {@code x}.
         *
         * @param x Value for which the function must be computed.
         * @param param Values of norm, mean and standard deviation.
         * @return the value of the function.
         * @ if {@code param} is {@code NULL}.
         * @ if the size of {@code param} is
         * not 3.
         */
         //override
        double value(const double& x, double ... param)
        {
            validate_parameters(param);
            return Harmonic_Oscillator::value(x * param[1] + param[2], param[0]);
        }

        /**
         * Computes the value of the gradient at {@code x}.
         * The components of the gradient vector are the partial
         * derivatives of the function with respect to each of the
         * <em>parameters</em> (amplitude, angular frequency and phase).
         *
         * @param x Value at which the gradient must be computed.
         * @param param Values of amplitude, angular frequency and phase.
         * @return the gradient vector at {@code x}.
         * @ if {@code param} is {@code NULL}.
         * @ if the size of {@code param} is
         * not 3.
         */
         //override
        std::vector<double> gradient(const double& x, double ... param)
        {
            validate_parameters(param);

            const auto amplitude = param[0];
            const auto omega = param[1];
            const auto phase = param[2];

            const double x_times_omega_plus_phase = omega * x + phase;
            const double a = Harmonic_Oscillator::value(x_times_omega_plus_phase, 1);
            const double p = -amplitude * std::sin(x_times_omega_plus_phase);
            const double w = p * x;

            return std::vector<double> { a, w, p };
        }

    private:
        /**
         * Validates parameters to ensure they are appropriate for the evaluation of
         * the {@link #value(double,std::vector<double>)} and {@link #gradient(double,std::vector<double>)}
         * methods.
         *
         * @param param Values of norm, mean and standard deviation.
         * @ if {@code param} is {@code NULL}.
         * @ if the size of {@code param} is
         * not 3.
         */
        void validate_parameters(const std::vector<double>& param)
        {
            //Math_Utils::check_not_null(param);
            Math_Utils::check_dimension(param.size(), 3);
        }
    };

    /** {@inherit_doc}
     */
    //override
    template<typename T, typename std::enable_if<std::is_base_of<Derivative<T>, T>::value>::type* = nullptr>
    T value(T t)
    {
        const double x = t.get_value();
        auto f = std::vector<double>(t.get_order() + 1];

        const double alpha   = omega * x + phase;
        const auto sc_alpha = Sin_Cos(alpha);
        f[0] = amplitude * sc_alpha.cos();
        if (f.size() > 1) 
        {
            f[1] = -amplitude * omega * sc_alpha.sin();
            const double mo2 = - omega * omega;
            for (int i{ 2 }; i < f.size(); ++i) 
            {
                f[i] = mo2 * f[i - 2];
            }
        }

        return t.compose(f);

    }

};