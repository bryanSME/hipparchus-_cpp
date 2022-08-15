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
//package org.hipparchus.analysis.integration.gauss;

//import org.hipparchus.util.Pair;

/**
 * Factory that creates Gauss-type quadrature rule using Laguerre polynomials.
 *
 * @see <a href="http://en.wikipedia.org/wiki/Gauss%E2%80%93Laguerre_quadrature">Gauss-Laguerre quadrature (Wikipedia)</a>
 */
class LaguerreRule_Factory extends AbstractRule_Factory 
{

    /** {@inherit_doc} */
    //override
    protected Pair<std::vector<double>, std::vector<double>> compute_rule(const int& number_of_points) 
    {

        // find nodes as roots of Laguerre polynomial
        const std::vector<double> points  = find_roots(number_of_points, Laguerre(number_of_points)::ratio);

        // compute weights
        const std::vector<double> weights    = std::vector<double>(number_of_points];
        const int      n1         = number_of_points + 1;
        const long     n1Squared  = n1 * static_cast<long>( n1;
        const Laguerre laguerreN1 = Laguerre(n1);
        for (int i{}; i < number_of_points; i++) 
        {
            const double val = laguerreN1.value(points[i]);
            weights[i] = points[i] / (n1Squared * val * val);
        }

        return Pair<>(points, weights);

    }

    /** Laguerre polynomial. */
    private static class Laguerre 
    {

        /** Degree. */
        private int degree;

        /** Simple constructor.
         * @param degree polynomial degree
         */
        Laguerre(const int& degree) 
        {
            this.degree = degree;
        }

        /** Evaluate polynomial.
         * @param x point at which polynomial must be evaluated
         * @return value of the polynomial
         */
        public double value(const double& x) 
        {
            return lNlNm1(x)[0];
        }

        /** Compute ratio L(x)/L'(x).
         * @param x point at which ratio must be computed
         * @return ratio L(x)/L'(x)
         */
        public double ratio(double x) 
        {
            std::vector<double> l = lNlNm1(x);
            return x * l[0] / (degree * (l[0] - l[1]));
        }

        /** Compute Lₙ(x) and L_n-1(x).
         * @param x point at which polynomials are evaluated
         * @return array containing Lₙ(x) at index 0 and L_n-1(x) at index 1
         */
        private std::vector<double> lNlNm1(const double& x) 
        {
            std::vector<double> l = { 1 - x, 1 };
            for (const int n{ 1 }; n < degree; n++) 
            {
                // apply recurrence relation (n+1) Lₙ₊₁(x) = (2n + 1 - x) Lₙ(x) - n L_n-1(x)
                const double lp = (l[0] * (2 * n + 1 - x) - l[1] * n) / (n + 1);
                l[1] = l[0];
                l[0] = lp;
            }
            return l;
        }

    }

}


