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

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.stat.descriptive.Storeless_Multivariate_Statistic;
//import org.hipparchus.stat.descriptive.Storeless_Univariate_Statistic;
//import org.hipparchus.util.Math_Utils;

#include "../StorelessMultivariateStatistic.h"
#include "../StorelessUnivariateStatistic.h"
#include <vector>

/**
 * Uses an independent {@link Storeless_Univariate_Statistic} instance
 * for each component of a vector.
 */
class Vectorial_Storeless_Statistic : public Storeless_Multivariate_Statistic
{
private:
    /** Statistic for each component */
     std::vector<Storeless_Univariate_Statistic> my_stats;

    /**
     * Create a Vectorial_Storeless_Statistic with the given dimension
     * and statistic implementation. A copy of the provided statistic
     * will be created for each component of the vector.
     *
     * @param dimension the vector dimension
     * @param univariate_statistic the prototype statistic
     * @ if dimension &lt; 1
     */
    public Vectorial_Storeless_Statistic(const int& dimension, Storeless_Univariate_Statistic univariate_statistic) 
    {
        if (dimension < 1) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL, dimension, 1);
        }
        my_stats = Storeless_Univariate_Statistic[dimension];
        for (int i{}; i < dimension; i++) 
        {
            my_stats[i] = univariate_statistic.copy();
        }
    }

    /** {@inherit_doc} */
    //override
    public void increment(std::vector<double> d) 
    {
        Math_Utils::check_dimension(d.size(), my_stats.size());
        for (int i{}; i < my_stats.size(); i++)
        {
            my_stats[i].increment(d[i]);
        }
    }

    /** {@inherit_doc} */
    //override
    public std::vector<double> get_result() 
    {
        std::vector<double> result = std::vector<double>(my_stats.size()];
        for (int i{}; i < result.size(); ++i) 
        {
            result[i] = my_stats[i].get_result();
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    public long get_n() 
    {
        return my_stats[0].get_n();
    }

    /** {@inherit_doc} */
    //override
    public void clear() 
    {
        for (auto& stat : my_stats)
        {
            stat.clear();
        }
    }

    /** {@inherit_doc} */
    //override
    public int get_dimension() 
    {
        return my_stats.size();
    }

    /** {@inherit_doc} */
    //override
    public int hash_code() 
    {
        constexpr int prime{ 31 };
        int result{ 1 };
        result = prime * result + Arrays.hash_code(my_stats);
        return result;
    }

    /** {@inherit_doc} */
    //override
    public bool equals(const Object& obj) 
    {
        if (this == obj) 
        {
            return true;
        }
        if (!(obj instanceof Vectorial_Storeless_Statistic)) 
        {
            return false;
        }
        Vectorial_Storeless_Statistic other = (Vectorial_Storeless_Statistic) obj;
        if (!Arrays.equals(my_stats, other.stats))
        {
            return false;
        }
        return true;
    }

};