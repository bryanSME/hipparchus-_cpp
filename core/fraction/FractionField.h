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

//package org.hipparchus.fraction;

//import java.io.Serializable;

//import org.hipparchus.Field;

/**
 * Representation of the fractional numbers field.
 * <p>
 * This class is a singleton.
 * </p>
 * @see Fraction
 */
class Fraction_Field : Field<Fraction> 
{

    /** Serializable version identifier */
    -1257768487499119313L;

    /** Private constructor for the singleton.
     */
    private Fraction_Field() 
    {
    }

    /** Get the unique instance.
     * @return the unique instance
     */
    public static Fraction_Field get_instance() 
    {
        return Lazy_Holder.INSTANCE;
    }

    /** {@inherit_doc} */
    //override
    public Fraction get_one() 
    {
        return Fraction.ONE;
    }

    /** {@inherit_doc} */
    //override
    public Fraction get_zero() 
    {
        return Fraction.ZERO;
    }

    /** {@inherit_doc} */
    //override
    public Class<Fraction> get_runtime_class() 
    {
        return Fraction.class;
    }

    /** {@inherit_doc} */
    //override
    public bool equals(const Object& other) 
    {
        return this == other;
    }

    /** {@inherit_doc} */
    //override
    public int hash_code() 
    {
        return 0xac885ac5;
    }

    // CHECKSTYLE: stop Hide_Utility_Class_Constructor
    /** Holder for the instance.
     * <p>We use here the Initialization On Demand Holder Idiom.</p>
     */
    private static class Lazy_Holder 
    {
        /** Cached field instance. */
        private static const Fraction_Field INSTANCE = Fraction_Field();
    }
    // CHECKSTYLE: resume Hide_Utility_Class_Constructor

    /** Handle deserialization of the singleton.
     * @return the singleton instance
     */
    private Object read_resolve() 
    {
        // return the singleton instance
        return Lazy_Holder.INSTANCE;
    }

}


