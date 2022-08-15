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

//import java.lang.reflect.Array;
//import java.util.Hash_Map;
//import java.util.Map;
#include <type_traits>
#include <unordered_map>
#include "../../CalculusFieldElement.hpp"
#include "../../Field.h"


//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.Field;
//import org.hipparchus.util.Math_Arrays;

/** Field for {@link Gradient} instances.
 * @param <T> the type of the function parameters and value
 * @since 1.7
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
class Field_Gradient_Field : Field<Field_Gradient<T>> 
{
private:
    /** Cached fields. */
    static const std::unordered_map<Field<>, std::vector<Field_Gradient_Field<>> MYCACHE = Hash_Map<>();

    /** Zero constant. */
    const Field_Gradient<T> my_zero;

    /** One constant. */
    const Field_Gradient<T> my_one;

    /** Associated factory for conversions to {@link Derivative_Structure}. */
    const FDS_Factory<T> my_factory;

    /** Private constructor.
     * @param value_field field for the function parameters and value
     * @param parameters number of free parameters
     */
    Field_Gradient_Field(const Field<T>& value_field, const int& parameters)
        :
        my_zero{ Field_Gradient<>(value_field.get_zero(), Math_Arrays::build_array(value_field, parameters)) },
        my_one{ Field_Gradient<>(value_field.get_one(), Math_Arrays::build_array(value_field, parameters)) },
        my_factory{ FDS_Factory<>(value_field, parameters, 1) }
    {
    };

public:
    /** Get the field for number of free parameters.
     * @param value_field field for the function parameters and value
     * @param parameters number of free parameters
     * @param <T> the type of the function parameters and value
     * @return cached field
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static Field_Gradient_Field<T> get_field(const Field<T>& value_field, const int& parameters) 
    {
        std::vector<Field_Gradient_Field<>> cached_fields;
        synchronized (CACHE) 
        {
            cached_fields = CACHE.get(value_field);
            if (cached_fields == NULL || cached_fields.size() <= parameters) 
            {
                std::vector<Field_Gradient_Field<>> new_cached_fields =
                                (std::vector<Field_Gradient_Field<>>) Array.new_instance(Field_Gradient_Field.class, parameters + 1);
                if (cached_fields != NULL) 
                {
                    // preserve the already created fields
                    System.arraycopy(cached_fields, 0, new_cached_fields, 0, cached_fields.size());
                }
                cached_fields = new_cached_fields;
                CACHE.put(value_field, cached_fields);
            }
        }

        if (cached_fields[parameters] == NULL) 
        {
            // we need to create a field
            cached_fields[parameters] = Field_Gradient_Field<>(value_field, parameters);
        }

        //@Suppress_Warnings("unchecked")
        const Field_Gradient_Field<T> t_cached = (Field_Gradient_Field<T>) cached_fields[parameters];
        return t_cached;
    }

    /** {@inherit_doc} */
    //override
    Field_Gradient<T>  get_one() const
    {
        return one;
    }

    /** {@inherit_doc} */
    //override
    Field_Gradient<T>  get_zero() const
    {
        return zero;
    }

    /** {@inherit_doc} */
    //@Suppress_Warnings("unchecked")
    //override
    Class<Field_Gradient<T>> get_runtime_class() 
    {
        return (Class<Field_Gradient<T>>) get_zero().get_class();
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
    FDS_Factory<T> get_conversion_factory() const
    {
        return my_factory;
    }

    /** {@inherit_doc} */
    //override
    bool equals(const Object& other) const
    {
        return this == other;
    }

    /** {@inherit_doc} */
    //override
    int hash_code() const
    {
        return 0xcd3e92ee;
    }

};