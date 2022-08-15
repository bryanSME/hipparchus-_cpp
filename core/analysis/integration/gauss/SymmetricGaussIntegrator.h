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

//import org.hipparchus.analysis.Univariate_Function;
//import org.hipparchus.exception.;
//import org.hipparchus.util.Pair;

/**
 * This class's : {@link #integrate(Univariate_Function) integrate}
 * method assuming that the integral is symmetric about 0.
 * This allows to reduce numerical errors.
 *
 */
class Symmetric_Gauss_Integrator extends Gauss_Integrator 
{
    /**
     * Creates an integrator from the given {@code points} and {@code weights}.
     * The integration interval is defined by the first and last value of
     * {@code points} which must be sorted in increasing order.
     *
     * @param points Integration points.
     * @param weights Weights of the corresponding integration nodes.
     * @ if the {@code points} are not
     * sorted in increasing order.
     * @ if points and weights don't have the same length
     */
    public Symmetric_Gauss_Integrator(std::vector<double> points, std::vector<double> weights)
         
        {
        super(points, weights);
    }

    /**
     * Creates an integrator from the given pair of points (first element of
     * the pair) and weights (second element of the pair.
     *
     * @param points_and_weights Integration points and corresponding weights.
     * @ if the {@code points} are not
     * sorted in increasing order.
     *
     * @see #Symmetric_Gauss_Integrator(std::vector<double>, std::vector<double>)
     */
    public Symmetric_Gauss_Integrator(Pair<std::vector<double>, std::vector<double>> points_and_weights)
         
        {
        this(points_and_weights.get_first(), points_and_weights.get_second());
    }

    /**
     * {@inherit_doc}
     */
    //override
    public double integrate(Univariate_Function f) 
    {
        const int rule_length = get_number_of_points();

        if (rule_length == 1) 
        {
            return get_weight(0) * f.value(0d);
        }

        const int i_max = rule_length / 2;
        double s{};
        double c{};
        for (int i{}; i < i_max; i++) 
        {
            const double p = get_point(i);
            const double w = get_weight(i);

            const double f1 = f.value(p);
            const double f2 = f.value(-p);

            const double y = w * (f1 + f2) - c;
            const double t = s + y;

            c = (t - s) - y;
            s = t;
        }

        if (rule_length % 2 != 0) 
        {
            const double w = get_weight(i_max);

            const double y = w * f.value(0d) - c;
            const double t = s + y;

            s = t;
        }

        return s;
    }
}


