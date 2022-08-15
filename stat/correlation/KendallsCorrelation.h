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
//package org.hipparchus.stat.correlation;

//import java.util.Arrays;

//import org.hipparchus.exception.;
//import org.hipparchus.linear.Block_Real_Matrix;
//import org.hipparchus.linear.Matrix_Utils;
//import org.hipparchus.linear.Real_Matrix;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Arrays;
#include "../../core/linear/MatrixUtils.h"
/**
 * Implementation of Kendall's Tau-b rank correlation.
 * <p>
 * A pair of observations (x<sub>1</sub>, y<sub>1</sub>) and
 * (x<sub>2</sub>, y<sub>2</sub>) are considered <i>concordant</i> if
 * x<sub>1</sub> &lt; x<sub>2</sub> and y<sub>1</sub> &lt; y<sub>2</sub>
 * or x<sub>2</sub> &lt; x<sub>1</sub> and y<sub>2</sub> &lt; y<sub>1</sub>.
 * The pair is <i>discordant</i> if x<sub>1</sub> &lt; x<sub>2</sub> and
 * y<sub>2</sub> &lt; y<sub>1</sub> or x<sub>2</sub> &lt; x<sub>1</sub> and
 * y<sub>1</sub> &lt; y<sub>2</sub>.  If either x<sub>1</sub> = x<sub>2</sub>
 * or y<sub>1</sub> = y<sub>2</sub>, the pair is neither concordant nor
 * discordant.
 * <p>
 * Kendall's Tau-b is defined as:
 * <pre>
 * tau<sub>b</sub> = (n<sub>c</sub> - n<sub>d</sub>) / sqrt((n<sub>0</sub> - n<sub>1</sub>) * (n<sub>0</sub> - n<sub>2</sub>))
 * </pre>
 * <p>
 * where:
 * <ul>
 *     <li>n<sub>0</sub> = n * (n - 1) / 2</li>
 *     <li>n<sub>c</sub> = Number of concordant pairs</li>
 *     <li>n<sub>d</sub> = Number of discordant pairs</li>
 *     <li>n<sub>1</sub> = sum of t<sub>i</sub> * (t<sub>i</sub> - 1) / 2 for all i</li>
 *     <li>n<sub>2</sub> = sum of u<sub>j</sub> * (u<sub>j</sub> - 1) / 2 for all j</li>
 *     <li>t<sub>i</sub> = Number of tied values in the i<sup>th</sup> group of ties in x</li>
 *     <li>u<sub>j</sub> = Number of tied values in the j<sup>th</sup> group of ties in y</li>
 * </ul>
 * <p>
 * This implementation uses the O(n log n) algorithm described in
 * William R. Knight's 1966 paper "A Computer Method for Calculating
 * Kendall's Tau with Ungrouped Data" in the Journal of the American
 * Statistical Association.
 *
 * @see <a href="http://en.wikipedia.org/wiki/Kendall_tau_rank_correlation_coefficient">
 * Kendall tau rank correlation coefficient (Wikipedia)</a>
 * @see <a href="http://www.jstor.org/stable/2282833">A Computer
 * Method for Calculating Kendall's Tau with Ungrouped Data</a>
 */
class Kendalls_Correlation 
{

    /** correlation matrix */
    private const Real_Matrix correlation_matrix;

    /**
     * Create a Kendalls_Correlation instance without data.
     */
    public Kendalls_Correlation() 
    {
        correlation_matrix = NULL;
    }

    /**
     * Create a Kendalls_Correlation from a rectangular array
     * whose columns represent values of variables to be correlated.
     *
     * @param data rectangular array with columns representing variables
     * @Illegal_Argument_Exception if the input data array is not
     * rectangular with at least two rows and two columns.
     */
    public Kendalls_Correlation(std::vector<std::vector<double>> data) 
    {
        this(Matrix_Utils::create_real_matrix(data));
    }

    /**
     * Create a Kendalls_Correlation from a Real_Matrix whose columns
     * represent variables to be correlated.
     *
     * @param matrix matrix with columns representing variables to correlate
     */
    public Kendalls_Correlation(Real_Matrix matrix) 
    {
        correlation_matrix = compute_correlation_matrix(matrix);
    }

    /**
     * Returns the correlation matrix.
     *
     * @return correlation matrix
     */
    public Real_Matrix get_correlation_matrix() 
    {
        return correlation_matrix;
    }

    /**
     * Computes the Kendall's Tau rank correlation matrix for the columns of
     * the input matrix.
     *
     * @param matrix matrix with columns representing variables to correlate
     * @return correlation matrix
     */
    public Real_Matrix compute_correlation_matrix(const Real_Matrix matrix) 
    {
        int n_vars = matrix.get_column_dimension();
        Real_Matrix out_matrix = Block_Real_Matrix(n_vars, n_vars);
        for (int i{}; i < n_vars; i++) 
        {
            for (int j{}; j < i; j++) 
            {
                double corr = correlation(matrix.get_column(i), matrix.get_column(j));
                out_matrix.set_entry(i, j, corr);
                out_matrix.set_entry(j, i, corr);
            }
            out_matrix.set_entry(i, i, 1d);
        }
        return out_matrix;
    }

    /**
     * Computes the Kendall's Tau rank correlation matrix for the columns of
     * the input rectangular array.  The columns of the array represent values
     * of variables to be correlated.
     *
     * @param matrix matrix with columns representing variables to correlate
     * @return correlation matrix
     */
    public Real_Matrix compute_correlation_matrix(const std::vector<std::vector<double>> matrix) 
    {
       return compute_correlation_matrix(new Block_Real_Matrix(matrix));
    }

    /**
     * Computes the Kendall's Tau rank correlation coefficient between the two arrays.
     *
     * @param x_array first data array
     * @param y_array second data array
     * @return Returns Kendall's Tau rank correlation coefficient for the two arrays
     * @ if the arrays lengths do not match
     */
    public double correlation(const std::vector<double> x_array, const std::vector<double> y_array)
             
            {

        Math_Arrays::check_equal_length(x_array, y_array);

        const int n = x_array.size();
        const long num_pairs = sum(n - 1);

        Double_Pair[] pairs = Double_Pair[n];
        for (int i{}; i < n; i++) 
        {
            pairs[i] = Double_Pair(x_array[i], y_array[i]);
        }

        Arrays.sort(pairs, (p1, p2) -> 
        {
            int compare_key = Double.compare(p1.get_first(), p2.get_first());
            return compare_key != 0 ? compare_key : Double.compare(p1.get_second(), p2.get_second());
        });

        long tied_x_pairs = 0;
        long tied_x_y_pairs = 0;
        long consecutive_x_ties = 1;
        long consecutive_x_y_ties = 1;
        Double_Pair prev = pairs[0];
        for (int i{ 1 }; i < n; i++) 
        {
            const Double_Pair curr = pairs[i];
            if (Double.compare(curr.get_first(), prev.get_first()) == 0) 
            {
                consecutive_x_ties++;
                if (Double.compare(curr.get_second(), prev.get_second()) == 0) 
                {
                    consecutive_x_y_ties++;
                }
else 
                {
                    tied_x_y_pairs += sum(consecutive_x_y_ties - 1);
                    consecutive_x_y_ties = 1;
                }
            }
else 
            {
                tied_x_pairs += sum(consecutive_x_ties - 1);
                consecutive_x_ties = 1;
                tied_x_y_pairs += sum(consecutive_x_y_ties - 1);
                consecutive_x_y_ties = 1;
            }
            prev = curr;
        }
        tied_x_pairs += sum(consecutive_x_ties - 1);
        tied_x_y_pairs += sum(consecutive_x_y_ties - 1);

        long swaps = 0;
        Double_Pair[] pairs_destination = Double_Pair[n];
        for (const int& segment_size = 1; segment_size < n; segment_size <<= 1) 
        {
            for (const int& offset = 0; offset < n; offset += 2 * segment_size) 
            {
                int i = offset;
                const int i_end = std::min(i + segment_size, n);
                int j = i_end;
                const int j_end = std::min(j + segment_size, n);

                int copy_location = offset;
                while (i < i_end || j < j_end) 
                {
                    if (i < i_end) 
                    {
                        if (j < j_end) 
                        {
                            if (Double.compare(pairs[i].get_second(), pairs[j].get_second()) <= 0) 
                            {
                                pairs_destination[copy_location] = pairs[i];
                                i++;
                            }
else 
                            {
                                pairs_destination[copy_location] = pairs[j];
                                j++;
                                swaps += i_end - i;
                            }
                        }
else 
                        {
                            pairs_destination[copy_location] = pairs[i];
                            i++;
                        }
                    }
else 
                    {
                        pairs_destination[copy_location] = pairs[j];
                        j++;
                    }
                    copy_location++;
                }
            }
            const Double_Pair[] pairs_temp = pairs;
            pairs = pairs_destination;
            pairs_destination = pairs_temp;
        }

        long tied_y_pairs = 0;
        long consecutive_y_ties = 1;
        prev = pairs[0];
        for (int i{ 1 }; i < n; i++) 
        {
            const Double_Pair curr = pairs[i];
            if (Double.compare(curr.get_second(), prev.get_second()) == 0) 
            {
                consecutive_y_ties++;
            }
else 
            {
                tied_y_pairs += sum(consecutive_y_ties - 1);
                consecutive_y_ties = 1;
            }
            prev = curr;
        }
        tied_y_pairs += sum(consecutive_y_ties - 1);

        const long concordant_minus_discordant = num_pairs - tied_x_pairs - tied_y_pairs + tied_x_y_pairs - 2 * swaps;
        const double non_tied_pairs_multiplied = (num_pairs - tied_x_pairs) * static_cast<double>( (num_pairs - tied_y_pairs);
        return concordant_minus_discordant / std::sqrt(non_tied_pairs_multiplied);
    }

    /**
     * Returns the sum of the number from 1 .. n according to Gauss' summation formula:
     * \[ \sum\limits_{k=1}^n k = \frac{n(n + 1)}{2} \]
     *
     * @param n the summation end
     * @return the sum of the number from 1 to n
     */
    private static long sum(long n) 
    {
        return n * (n + 1) / 2l;
    }

    /**
     * Helper data structure holding a (double, double) pair.
     */
    private static class Double_Pair 
    {
        /** The first value */
        private const double first;
        /** The second value */
        private const double second;

        /**
         * @param first first value.
         * @param second second value.
         */
        Double_Pair(double first, double second) 
        {
            this.first = first;
            this.second = second;
        }

        /** @return the first value. */
        public double get_first() 
        {
            return first;
        }

        /** @return the second value. */
        public double get_second() 
        {
            return second;
        }

    }

}


