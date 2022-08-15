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
//package org.hipparchus.stat.descriptive.vector;

//import java.io.Serializable;
//import java.util.Arrays;
#include <vector>

//import org.hipparchus.exception.;
//import org.hipparchus.linear.Matrix_Utils;
//import org.hipparchus.linear.Real_Matrix;
//import org.hipparchus.util.Math_Arrays;
#include "../../../core/linear/MatrixUtils.h"
/**
 * Returns the covariance matrix of the available vectors.
 */
class Vectorial_Covariance  
{
private:

    /** Sums for each component. */
    std::vector<double> my_sums;

    /** Sums of products for each component. */
    const std::vector<double> my_products_sums;

    /** Indicator for bias correction. */
    const bool is_bias_corrected;

    /** Number of vectors in the sample. */
    long my_n;

public:
    /** Constructs a Vectorial_Covariance.
     * @param dimension vectors dimension
     * @param is_bias_corrected if true, computed the unbiased sample covariance, * otherwise computes the biased population covariance
     */
    Vectorial_Covariance(const int& dimension, bool is_bias_corrected) 
    {
        my_sums = std::vector<double>(dimension];
        my_products_sums = std::vector<double>(dimension * (dimension + 1) / 2];
        my_n = 0;
        this.is_bias_corrected = is_bias_corrected;
    }

    /**
     * Add a vector to the sample.
     * @param v vector to add
     * @ if the vector does not have the right dimension
     */
    void increment(std::vector<double> v)  
    {
        Math_Arrays::check_equal_length(v, my_sums);
        int k = 0;
        for (int i{}; i < v.size(); ++i) 
        {
            my_sums[i] += v[i];
            for (int j{}; j <= i; ++j) 
            {
                my_products_sums[k++] += v[i] * v[j];
            }
        }
        my_n++;
    }

    /**
     * Get the covariance matrix.
     * @return covariance matrix
     */
    Real_Matrix get_result() 
    {

        int dimension = my_sums.size();
        Real_Matrix result = Matrix_Utils::create_real_matrix(dimension, dimension);

        if (n > 1) 
        {
            double c = 1.0 / (n * (is_bias_corrected ? (n - 1) : my_n));
            int k = 0;
            for (int i{}; i < dimension; ++i) 
            {
                for (int j{}; j <= i; ++j) 
                {
                    double e = c * (n * my_products_sums[k++] - my_sums[i] * my_sums[j]);
                    result.set_entry(i, j, e);
                    result.set_entry(j, i, e);
                }
            }
        }

        return result;

    }

    /**
     * Get the number of vectors in the sample.
     * @return number of vectors in the sample
     */
    long get_n() const
    {
        return my_n;
    }

    std::vector<double> get_sums() const
    {
        return my_sums;
    }

    /**
     * Clears the internal state of the Statistic
     */
    void clear() 
    {
        my_n = 0;
        Arrays.fill(my_sums, 0.0);
        Arrays.fill(my_products_sums, 0.0);
    }

    /** {@inherit_doc} */
    //override
    int hash_code() 
    {
        constexpr int prime{ 31 };
        int result{ 1 };
        result = prime * result + (
            is_bias_corrected
                ? 1231
                : 1237
        );
        result = prime * result + static_cast<int>(n ^ (n >>> 32));
        result = prime * result + Arrays.hash_code(my_products_sums);
        result = prime * result + Arrays.hash_code(my_sums);
        return result;
    }

    /** {@inherit_doc} */
    //override
    bool equals(const Object& obj) 
    {
        if (this == obj) 
        {
            return true;
        }
        if (!(obj instanceof Vectorial_Covariance)) 
        {
            return false;
        }
        Vectorial_Covariance other = (Vectorial_Covariance) obj;
        if (is_bias_corrected != other.is_bias_corrected) 
        {
            return false;
        }
        if (n != other.n) 
        {
            return false;
        }
        if (!Arrays.equals(my_products_sums, other.products_sums))
        {
            return false;
        }
        if (!Arrays.equals(my_sums, other.get_sums()))
        {
            return false;
        }
        return true;
    }

}


