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

#include "../../CalculusFieldElement.hpp"
#include "Derivative.h"
#include "DerivativeStructure.h"


/** Abstract class representing both the value and the differentials of a function.
 * @param <T> the type of the function derivative
 * @since 1.7
 */
//<T extends Univariate_Derivative<T>>
template<typename T, typename std::enable_if<std::is_base_of<Univariate_Derivative, T>::value>::type* = nullptr>
class Univariate_Derivative : public Derivative<T>, public Calculus_Field_Element<T>
{
public:
    /** {@inherit_doc} */
    //override
    int get_free_parameters()
    {
        return 1;
    }

    /** {@inherit_doc} */
    //override
    double get_partial_derivative(const int ... orders)
    {
        if (orders.size() != 1)
        {
            throw std::exception("not implmented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, orders.size(), 1);
        }
        return get_derivative(orders[0]);
    }

    /** Get a derivative from the univariate derivative.
     * @param n derivation order (must be between 0 and {@link #get_order()}, both inclusive)
     * @return n<sup>th</sup> derivative
     * @exception  if n is
     * either negative or strictly larger than {@link #get_order()}
     */
    virtual double get_derivative(const int& n);

    /** Convert the instance to a {@link Derivative_Structure}.
     * @return derivative structure with same value and derivative as the instance
     */
    virtual Derivative_Structure to_derivative_structure();

};