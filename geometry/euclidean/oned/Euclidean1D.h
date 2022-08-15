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

//package org.hipparchus.geometry.euclidean.oned;

//import java.io.Serializable;

//import org.hipparchus.exception.Math_Runtime_Exception;
//import org.hipparchus.geometry.Localized_Geometry_Formats;
//import org.hipparchus.geometry.Space;

/**
 * This class : a one-dimensional space.
 */
class Euclidean_1D : public Space 
{
private:
    /** Private constructor for the singleton.
     */
    Euclidean_1D() = default;

    /** Handle deserialization of the singleton.
    * @return the singleton instance
    */
    Object read_resolve()
    {
        // return the singleton instance
        return Lazy_Holder.INSTANCE;
    }

public:
    /** Get the unique instance.
     * @return the unique instance
     */
    static Euclidean_1D get_instance() 
    {
        return Lazy_Holder.INSTANCE;
    }

    /** {@inherit_doc} */
    //override
    int get_dimension() const
    {
        return 1;
    }

    /** {@inherit_doc}
     * <p>
     * As the 1-dimension Euclidean space does not have proper sub-spaces, * this method always a {@link No_Sub_Space_Exception}
     * </p>
     * @return nothing
     * @No_Sub_Space_Exception in all cases
     */
    //override
    Space get_sub_space() 
    {
        std::exception("not implemented");
        //throw No_Sub_Space_Exception();
    }

    // CHECKSTYLE: stop Hide_Utility_Class_Constructor
    /** Holder for the instance.
     * <p>We use here the Initialization On Demand Holder Idiom.</p>
     */
    private static class Lazy_Holder
    {
    private:
        /** Cached field instance. */
        static const Euclidean_1D INSTANCE = Euclidean_1D();
    };
    // CHECKSTYLE: resume Hide_Utility_Class_Constructor


    /** Specialized exception for inexistent sub-space.
     * <p>
     * This exception is thrown when attempting to get the sub-space of a one-dimensional space
     * </p>
     */
    public static class No_Sub_Space_Exception extends Math_Runtime_Exception 
    {

        
        20140225L;

        /** Simple constructor.
         */
        public No_Sub_Space_Exception() 
        {
            super(Localized_Geometry_Formats.NOT_SUPPORTED_IN_DIMENSION_N, 1);
        }

    }

}


