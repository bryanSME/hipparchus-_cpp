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

#include "../../special/Gamma.h"
#include "../../util/Precision.h"

/**
 * Utility class used by various distributions to accurately compute their
 * respective probability mass functions. The implementation for this class is
 * based on the Catherine Loader's <a target="_blank"
 * href="http://www.herine.net/stat/software/dbinom.html">dbinom</a> routines.
 * <p>
 * This class is not intended to be called directly.
 * <p>
 * References:
 * <ol>
 * <li>Catherine Loader (2000). "Fast and Accurate Computation of Binomial
 * Probabilities.". <a target="_blank"
 * href="http://www.herine.net/stat/papers/dbinom.pdf">
 * http://www.herine.net/stat/papers/dbinom.pdf</a></li>
 * </ol>
 */
const class Saddle_Point_Expansion 
{
private:
    /** 1/2 * log(2 &#960;). */
    static const double HALF_LOG_2_PI = 0.5 * std::log(Math_Utils::TWO_PI);

    /** exact Stirling expansion error for certain values. */
    static const std::vector<double> EXACT_STIRLING_ERRORS = 
    {
        0.0,                           /* 0.0 */
        0.1534264097200273452913848,   /* 0.5 */
        0.0810614667953272582196702,   /* 1.0 */
        0.0548141210519176538961390,   /* 1.5 */
        0.0413406959554092940938221,   /* 2.0 */
        0.03316287351993628748511048,  /* 2.5 */
        0.02767792568499833914878929,  /* 3.0 */
        0.02374616365629749597132920,  /* 3.5 */
        0.02079067210376509311152277,  /* 4.0 */
        0.01848845053267318523077934,  /* 4.5 */
        0.01664469118982119216319487,  /* 5.0 */
        0.01513497322191737887351255,  /* 5.5 */
        0.01387612882307074799874573,  /* 6.0 */
        0.01281046524292022692424986,  /* 6.5 */
        0.01189670994589177009505572,  /* 7.0 */
        0.01110455975820691732662991,  /* 7.5 */
        0.010411265261972096497478567, /* 8.0 */
        0.009799416126158803298389475, /* 8.5 */
        0.009255462182712732917728637, /* 9.0 */
        0.008768700134139385462952823, /* 9.5 */
        0.008330563433362871256469318, /* 10.0 */
        0.007934114564314020547248100, /* 10.5 */
        0.007573675487951840794972024, /* 11.0 */
        0.007244554301320383179543912, /* 11.5 */
        0.006942840107209529865664152, /* 12.0 */
        0.006665247032707682442354394, /* 12.5 */
        0.006408994188004207068439631, /* 13.0 */
        0.006171712263039457647532867, /* 13.5 */
        0.005951370112758847735624416, /* 14.0 */
        0.005746216513010115682023589, /* 14.5 */
        0.005554733551962801371038690  /* 15.0 */
    };

    /**
     * Default constructor.
     */
    Saddle_Point_Expansion() {}

public:
    /**
     * Compute the error of Stirling's series at the given value.
     * <p>
     * References:
     * <ol>
     * <li>Eric W. Weisstein. "Stirling's Series." From MathWorld--A Wolfram Web
     * Resource. <a target="_blank"
     * href="http://mathworld.wolfram.com/StirlingsSeries.html">
     * http://mathworld.wolfram.com/StirlingsSeries.html</a></li>
     * </ol>
     *
     * @param z the value.
     * @return the Striling's series error.
     */
    static double get_stirling_error(const double& z) 
    {
        if (z < 15.0) 
        {
            double z2 = 2.0 * z;
            if (Precision::is_mathematical_integer(z2)) 
            {
                return EXACT_STIRLING_ERRORS[static_cast<int>(z2)];
            }
            return Gamma::log_gamma(z + 1.0) - (z + 0.5) * std::log(z) +
                   z - HALF_LOG_2_PI;
        }
        const double z2 = z * z;
        return (0.083333333333333333333 -
                           (0.00277777777777777777778 -
                             (0.00079365079365079365079365 -
                               (0.000595238095238095238095238 -
                                 0.0008417508417508417508417508 /
                                   z2) / z2) / z2) / z2) / z;
    }

    /**
     * A part of the deviance portion of the saddle point approximation.
     * <p>
     * References:
     * <ol>
     * <li>Catherine Loader (2000). "Fast and Accurate Computation of Binomial
     * Probabilities.". <a target="_blank"
     * href="http://www.herine.net/stat/papers/dbinom.pdf">
     * http://www.herine.net/stat/papers/dbinom.pdf</a></li>
     * </ol>
     *
     * @param x the x value.
     * @param mu the average.
     * @return a part of the deviance.
     */
    static double get_deviance_part(const double& x, const double& mu) 
    {
        if (std::abs(x - mu) < 0.1 * (x + mu)) 
        {
            double d = x - mu;
            double v = d / (x + mu);
            double s1 = v * d;
            double s = std::numeric_limits<double>::quiet_NaN();
            double ej = 2.0 * x * v;
            v *= v;
            int j = 1;
            while (s1 != s) 
            {
                s = s1;
                ej *= v;
                s1 = s + ej / ((j * 2) + 1);
                ++j;
            }
            return s1;
        }
        return x * std::log(x / mu) + mu - x;
    }

    /**
     * Compute the logarithm of the PMF for a binomial distribution
     * using the saddle point expansion.
     *
     * @param x the value at which the probability is evaluated.
     * @param n the number of trials.
     * @param p the probability of success.
     * @param q the probability of failure (1 - p).
     * @return log(p(x)).
     */
    static double log_binomial_probability(const int& x, const int& n, const double& p, const double& q) 
    {
        if (n == 0) 
        {
            return x == 0
                ? 0
                : -INFINITY;
        }
        if (x == 0) 
        {
            return p < 0.1 
                ? -get_deviance_part(n, n * q) - n * p
                : n * std::log(q);
        }
        if (x == n) 
        {
            return q < 0.1
                ? -get_deviance_part(n, n * p) - n * q
                : n * std::log(p);
        }

        double ret = get_stirling_error(n) - get_stirling_error(x) -
                        get_stirling_error(n - x) - get_deviance_part(x, n * p) -
                        get_deviance_part(n - x, n * q);
        double f = (Math_Utils::TWO_PI * x * (n - x)) / n;
        ret = -0.5 * std::log(f) + ret;
        return ret;
        
    }
};