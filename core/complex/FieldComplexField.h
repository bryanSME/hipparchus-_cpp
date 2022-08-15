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
//package org.hipparchus.complex;

//import java.util.Hash_Map;
//import java.util.Map;

//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.Field;
#include <type_traits>
#include "../CalculusFieldElement.hpp"
#include <unordered_map>

/**
 * Representation of the complex numbers field.
 * <p>
 *
 * @param <T> the type of the field elements
 * @see Field_Complex<double>
 * @since 2.0
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
class Field_Complex_Field : public Field<Field_Complex<T>>
{
private:
    /** Cached fields. */
    static const Map<Field< ? >, Field_Complex_Field< ? >> CACHE = std::unordered_map<>();

    /** Constant 0. */
    const Field_Complex<T> my_zero;

    /** Constant 1. */
    const Field_Complex<T> my_one;

    /** Simple constructor.
     * @param field type of the field element
     */
    Field_Complex_Field(const Field<T> field) 
    {
        zero = Field_Complex<double>.get_zero(field);
        one = Field_Complex<double>.get_one(field);
    }

public:
    /** Get the field for complex numbers.
     * @param parts_field field for the real and imaginary parts
     * @param <T> the type of the field elements
     * @return cached field
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  Field_Complex_Field<T> get_field(const Field<T> parts_field)
    {
        Field_Complex_Field< ? > cached_field;
        synchronized(CACHE)
        {
            cached_field = CACHE.get(parts_field);
            if (cached_field == NULL)
            {
                cached_field = Field_Complex_Field<>(parts_field);
                CACHE.put(parts_field, cached_field);
            }
        }

        //@Suppress_Warnings("unchecked")
        const Field_Complex_Field<T> t_cached = (Field_Complex_Field<T>) cached_field;
        return t_cached;

    }

    /** {@inherit_doc} */
    //override
    Field_Complex<T> get_one() const
    {
        return my_one;
    }

    /** {@inherit_doc} */
    //override
    Field_Complex<T> get_zero() const
    {
        return my_zero;
    }

    /** {@inherit_doc} */
    //@Suppress_Warnings("unchecked")
    //override
    Class<Field_Complex<T>> get_runtime_class()
    {
        return (Class<Field_Complex<T>>) get_zero().get_class();
    }

    /** {@inherit_doc} */
    //override
    bool equals(const Object& other)
    {
        return this == other;
    }

    /** {@inherit_doc} */
    //override
    int hash_code() const
    {
        return 0xd368f208;
    }

};