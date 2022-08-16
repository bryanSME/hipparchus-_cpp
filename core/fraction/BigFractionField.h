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
 * Representation of the fractional numbers  without any overflow field.
 * <p>
 * This class is a singleton.
 * </p>
 * @see Fraction
 */
class Big_Fraction_Field : public Field<Big_Fraction>
{
private:
    /** Private constructor for the singleton.
     */
    Big_Fraction_Field()
    {
    }

    // CHECKSTYLE: stop Hide_Utility_Class_Constructor
    /** Holder for the instance.
     * <p>We use here the Initialization On Demand Holder Idiom.</p>
     */
    static class Lazy_Holder
    {
    private:
        /** Cached field instance. */
        static const Big_Fraction_Field INSTANCE = Big_Fraction_Field();
    }
    // CHECKSTYLE: resume Hide_Utility_Class_Constructor

    /** Handle deserialization of the singleton.
     * @return the singleton instance
     */
    Object read_resolve() const
    {
        // return the singleton instance
        return Lazy_Holder::INSTANCE;
    }

public:
    /** Get the unique instance.
     * @return the unique instance
     */
    static Big_Fraction_Field get_instance()
    {
        return Lazy_Holder::INSTANCE;
    }

    /** {@inherit_doc} */
    //override
    Big_Fraction get_one() const
    {
        return Big_Fraction::ONE;
    }

    /** {@inherit_doc} */
    //override
    Big_Fraction get_zero() const
    {
        return Big_Fraction::ZERO;
    }

    /** {@inherit_doc} */
    //override
    Class<Big_Fraction> get_runtime_class()
    {
        return Big_Fraction.class;
    }

    /** {@inherit_doc} */
    //override
    bool equals(const Object& other) const
    {
        return *this == other;
    }

    /** {@inherit_doc} */
    //override
    int hash_code() const
    {
        return 0x7666e832;
    }
};