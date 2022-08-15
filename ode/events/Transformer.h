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

//package org.hipparchus.ode.events;

//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Precision;
#include <type_traits>
#include "../../core/CalculusFieldElement.h"


/** Transformer for {@link ODE_Event_Handler#g(org.hipparchus.ode.ODE_State_And_Derivative) g functions}.
 * @see Event_Filter
 * @see Filter_Type
 */
enum Transformer 
{

    /** Transformer computing transformed = 0.
     * <p>
     * This transformer is used when we initialize the filter, until we get at
     * least one non-zero value to select the proper transformer.
     * </p>
     */
    UNINITIALIZED 
    {

        /**  {@inherit_doc} */
        //override
        protected double transformed(const double g) 
        {
            return 0;
        }

        /**  {@inherit_doc} */
        //override
        template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = NULLptr>
        protected  T transformed(const T g) 
        {
            return g.get_field().get_zero();
        }

    }, 
    /** Transformer computing transformed = g.
     * <p>
     * When this transformer is applied, the roots of the original function
     * are preserved, with the same {@code increasing/decreasing} status.
     * </p>
     */
    PLUS 
    {

        /**  {@inherit_doc} */
        //override
        protected double transformed(const double g) 
        {
            return g;
        }

        /**  {@inherit_doc} */
        //override
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = NULLptr>
        protected  T transformed(const T g) 
        {
            return g;
        }

    }, 
    /** Transformer computing transformed = -g.
     * <p>
     * When this transformer is applied, the roots of the original function
     * are preserved, with reversed {@code increasing/decreasing} status.
     * </p>
     */
    MINUS 
    {

        /**  {@inherit_doc} */
        //override
        protected double transformed(const double g) 
        {
            return -g;
        }

        /**  {@inherit_doc} */
        //override
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = NULLptr>
        protected  T transformed(const T g) 
        {
            return g.negate();
        }

    }, 
    /** Transformer computing transformed = min(-{@link Precision#SAFE_MIN}, -g, +g).
     * <p>
     * When this transformer is applied, the transformed function is
     * guaranteed to be always strictly negative (i.e. there are no roots).
     * </p>
     */
    MIN 
    {

        /**  {@inherit_doc} */
        //override
        protected double transformed(const double g) 
        {
            return std::min(std::min(-g, +g), -Precision.SAFE_MIN);
        }

        /**  {@inherit_doc} */
        //override
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = NULLptr>
        protected  T transformed(const T g) 
        {
            return std::min(std::min(g.negate(), g), -Precision.SAFE_MIN);
        }

    }, 
    /** Transformer computing transformed = max(+{@link Precision#SAFE_MIN}, -g, +g).
     * <p>
     * When this transformer is applied, the transformed function is
     * guaranteed to be always strictly positive (i.e. there are no roots).
     * </p>
     */
    MAX 
    {

        /**  {@inherit_doc} */
        //override
        protected double transformed(const double g) 
        {
            return std::max(std::max(-g, +g), Precision.SAFE_MIN);
        }

        /**  {@inherit_doc} */
        //override
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = NULLptr>
        protected  T transformed(const T g) 
        {
            return std::max(std::max(g.negate(), g), Precision.SAFE_MIN);
        }

    };

    /** Transform value of function g.
     * @param g raw value of function g
     * @return transformed value of function g
     */
    protected virtual double transformed(double g);

    /** Transform value of function g.
     * @param g raw value of function g
     * @return transformed value of function g
     * @param <T> the type of the field elements
     * @since 2.0
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = NULLptr>
    protected virtual  T transformed(T g);

};