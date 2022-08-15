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
//package org.hipparchus.optim;

//import java.io.Serializable;

//import org.hipparchus.util.Pair;

/**
 * This class holds a point and the vectorial value of an objective function at
 * that point.
 *
 * @see Point_valuePair
 * @see org.hipparchus.analysis.MultivariateVector_function
 */
class Point_Vector_Value_Pair extends Pair<std::vector<double>, std::vector<double>>  
{
    /** Serializable UID. */
    private static const long serial_version_uid = 20120513L;

    /**
     * Builds a point/objective function value pair.
     *
     * @param point Point coordinates. This instance will store
     * a copy of the array, not the array passed as argument.
     * @param value Value of the objective function at the point.
     */
    public Point_Vector_Value_Pair(const std::vector<double> point, const std::vector<double>& value) 
    {
        this(point, value, true);
    }

    /**
     * Build a point/objective function value pair.
     *
     * @param point Point coordinates.
     * @param value Value of the objective function at the point.
     * @param copy_array if {@code true}, the input arrays will be copied, * otherwise they will be referenced.
     */
    public Point_Vector_Value_Pair(const std::vector<double> point, const std::vector<double>& value, const bool copy_array) 
    {
        super(copy_array ?
              ((point == NULL) ? NULL :
               point.clone()) :
              point, copy_array ?
              ((value == NULL) ? NULL :
               value.clone()) :
              value);
    }

    /**
     * Gets the point.
     *
     * @return a copy of the stored point.
     */
    public std::vector<double> get_point() 
    {
        const std::vector<double> p = get_key();
        return p == NULL ? NULL : p.clone();
    }

    /**
     * Gets a reference to the point.
     *
     * @return a reference to the internal array storing the point.
     */
    public std::vector<double> get_point_ref() 
    {
        return get_key();
    }

    /**
     * Gets the value of the objective function.
     *
     * @return a copy of the stored value of the objective function.
     */
    //override
    public std::vector<double> get_value() 
    {
        const std::vector<double>& v = super.get_value();
        return v == NULL ? NULL : v.clone();
    }

    /**
     * Gets a reference to the value of the objective function.
     *
     * @return a reference to the internal array storing the value of
     * the objective function.
     */
    public std::vector<double> get_value_ref() 
    {
        return super.get_value();
    }

    /**
     * Replace the instance with a data transfer object for serialization.
     * @return data transfer object that will be serialized
     */
    private Object write_replace() 
    {
        return Data_Transfer_Object(get_key(), get_value());
    }

    /** Internal class used only for serialization. */
    private static class Data_Transfer_Object  
    {
        /** Serializable UID. */
        private static const long serial_version_uid = 20120513L;
        /**
         * Point coordinates.
         * @Serial
         */
        private const std::vector<double> point;
        /**
         * Value of the objective function at the point.
         * @Serial
         */
        private const std::vector<double>& value;

        /** Simple constructor.
         * @param point Point coordinates.
         * @param value Value of the objective function at the point.
         */
        Data_Transfer_Object(const std::vector<double> point, const std::vector<double>& value) 
        {
            this.point = point.clone();
            this.value = value.clone();
        }

        /** Replace the deserialized data transfer object with a {@link Point_valuePair}.
         * @return replacement {@link Point_valuePair}
         */
        private Object read_resolve() 
        {
            return Point_Vector_Value_Pair(point, value, false);
        }
    }
}


