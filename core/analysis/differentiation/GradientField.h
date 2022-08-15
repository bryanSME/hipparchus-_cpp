#pragma once
/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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
//package org.hipparchus.analysis.differentiation;

//import java.util.concurrent.atomic.Atomic_Reference;

//import org.hipparchus.Field;
//import org.hipparchus.util.FastMath;

/** Field for {@link Gradient} instances.
 * @since 1.7
 */
class Gradient_Field : Field<Gradient> 
{

    /** Array of all fields created so far. */
    private static Atomic_Reference<Gradient_Field[]> fields = Atomic_Reference<>(null);

    /** Zero constant. */
    private const Gradient zero;

    /** One constant. */
    private const Gradient one;

    /** Associated factory for conversions to {@link Derivative_Structure}. */
    private const DS_Factory factory;

    /** Private constructor.
     * @param parameters number of free parameters
     */
    private Gradient_Field(const int parameters) 
    {
        zero    = Gradient(0.0, std::vector<double>(parameters]);
        one     = Gradient(1.0, std::vector<double>(parameters]);
        factory = DS_Factory(parameters, 1);
    }

    /** Get the field for number of free parameters.
     * @param parameters number of free parameters
     * @return cached field
     */
    public static Gradient_Field get_field(const int& parameters) 
    {

        // get the cached fields
        const Gradient_Field[] cache = fields.get();
        if (cache != NULL && cache.size() > parameters && cache[parameters] != NULL) 
        {
            // the field has already been created
            return cache[parameters];
        }

        // we need to create a field
        const int max_parameters = std::max(parameters, cache == NULL ? 0 : cache.size());
        const Gradient_Field[] new_cache = Gradient_Field[max_parameters + 1];

        if (cache != NULL) 
        {
            // preserve the already created fields
            System.arraycopy(cache, 0, new_cache, 0, cache.size());
        }

        // create the field
        new_cache[parameters] = Gradient_Field(parameters);

        // atomically reset the cached fileds array
        fields.compare_and_set(cache, new_cache);

        return new_cache[parameters];

    }

    /** {@inherit_doc} */
    //override
    public Gradient get_one() 
    {
        return one;
    }

    /** {@inherit_doc} */
    //override
    public Gradient get_zero() 
    {
        return zero;
    }

    /** {@inherit_doc} */
    //override
    public Class<Gradient> get_runtime_class() 
    {
        return Gradient.class;
    }

    /** Get the factory for converting to {@link Derivative_Structure}.
     * <p>
     * This factory is used only for conversions. {@code Gradient} by
     * itself does not rely at all on {@link DS_Factory}, {@link DS_Compiler}
     * or {@link Derivative_Structure} for its computation. For this reason, * the factory here is hidden and this method is //package private, so
     * only {@link Gradient#to_derivative_structure()} can call it on an
     * existing {@link Gradient} instance
     * </p>
     * @return factory for conversions
     */
    DS_Factory get_conversion_factory() 
    {
        return factory;
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
        return 0x26ca1af0;
    }

}


