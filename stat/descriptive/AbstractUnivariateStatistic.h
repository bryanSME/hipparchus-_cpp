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

#include <vector>
#include "UnivariateStatistic.h"

/**
 * Abstract base class for implementations of the
 * {@link Univariate_Statistic} interface.
 */
class Abstract_Univariate_Statistic : public Univariate_Statistic 
{
private:
    /** Stored data. */
    std::vector<double> my_stored_data;

protected:
    /** Default constructor. */
    Abstract_Univariate_Statistic() 
    {
        // This constructor is intentionally empty. Nothing special is needed here.
    }

    /**
     * Copy constructor, creates an identical copy
     * of the {@code original}.
     *
     * @param original the instance to copy
     * @org.hipparchus.exception.  if original is NULL
     */
    Abstract_Univariate_Statistic(const Abstract_Univariate_Statistic& original) 
    {
        //Math_Utils::check_not_null(original);
        my_stored_data = original.get_data_ref();
    }

    /**
     * Get a reference to the stored data array.
     * @return reference to the stored data array (may be NULL)
     */
    std::vector<double> get_data_ref() const
    {
        return my_stored_data; // NOPMD - returning an internal array is intentional and documented here
    }

public:

    /** {@inherit_doc} */
    //override
    virtual double evaluate(const std::vector<double>& values, const int& begin, const int& length) = 0;

    /** {@inherit_doc} */
    //override
    virtual Univariate_Statistic copy() = 0;

    /**
     * Set the data array.
     * <p>
     * The stored value is a copy of the parameter array, not the array itself.
     * </p>
     * @param values data array to store (may be NULL to remove stored data)
     * @see #evaluate()
     */
    void set_data(const std::vector<double>& values) 
    {
        my_stored_data = values;
    }

    /**
     * Get a copy of the stored data array.
     * @return copy of the stored data array (may be NULL)
     */
    std::vector<double> get_data() const
    {
        return my_stored_data;
    }

    /**
     * Set the data array.  The input array is copied, not referenced.
     *
     * @param values data array to store
     * @param begin the index of the first element to include
     * @param length the number of elements to include
     * @ if values is NULL or the indices
     * are not valid
     * @see #evaluate()
     */
    void set_data(const std::vector<double>& values, const int& begin, const int& length)
    {
        //Math_Utils::check_not_null(values, hipparchus::exception::Localized_Core_Formats_Type::INPUT_ARRAY);

        if (begin < 0) 
        {
            throw std::exception("error not implemented, begin < 0");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::START_POSITION, begin);
        }

        if (length < 0) 
        {
            throw std::exception("error not implemented, length < 0");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::size(), length);
        }

        if (begin + length > values.size()) 
        {
            throw std::exception("error not implemented, begin+length > values.size()");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::SUBARRAY_ENDS_AFTER_ARRAY_END, begin + length, values.size(), true);
        }
        my_stored_data = values;
    }

    /**
     * Returns the result of evaluating the statistic over the stored data.
     * <p>
     * The stored array is the one which was set by previous calls to
     * {@link #set_data(std::vector<double>)}.
     *
     * @return the value of the statistic applied to the stored data
     * @ if the stored data array is NULL
     */
    double evaluate()  
    {
        return evaluate(my_stored_data);
    }
};