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
//package org.hipparchus.util;

//import java.util.Arrays;
#include <vector>
#include <cmath>
#include "SinCos.h"
#include "../CalculusFieldElement.hpp"
#include "../Field.h"
//#include <algorithm>
//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.Field;
//import org.hipparchus.exception.;

/**
 * This class allows to perform the same computation of all components of a Tuple at once.
 * @since 1.2
 */
class Tuple : Calculus_Field_Element<Tuple> 
{
private:
    /** Components of the tuple. */
    const std::vector<double> my_values;

    /** Field the instance belongs to. */
    const Tuple_Field my_field;

    /** Creates a instance from its components.
     * @param field field the instance belongs to
     * @param x components of the tuple (beware, it is <em>not</em> copied, it is shared with caller)
     */
    Tuple(const Tuple_Field field, const std::vector<double> x) : my_values{ x }, my_field{field}
    {
        // NOPMD - storing user-supplied array is intentional and documented here
    }

public:
    /** Creates a instance from its components.
     * @param x components of the tuple
     */
    Tuple(const double... x) 
    {
        Tuple(Tuple_Field(x.size()), x.clone());
    }

    /** {@inherit_doc} */
    //override
    Tuple new_instance(const double value) 
    {
        const auto t = Tuple(my_field, std::vector<double>(my_values.size()));
        Arrays.fill(t.get_values(), value);
        return t;
    }

    /** Get the dimension of the tuple.
     * @return dimension of the tuple
     */
    int get_dimension() 
    {
        return my_values.size();
    }

    /** Get one component of the tuple.
     * @param index index of the component, between 0 and {@link #get_dimension() get_dimension()} - 1
     * @return value of the component
     */
    double get_component(const int& index) const
    {
        return my_values[index];
    }

    /** Get all components of the tuple.
     * @return all components
     */
    std::vector<double> get_components() const 
    {
        return my_values;
    }

    /** {@inherit_doc} */
    //override
    Field<Tuple> get_field() const
    {
        return my_field;
    }

    /** {@inherit_doc} */
    //override
    Tuple add(const Tuple& a) 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = my_values[i] + a.get_values()[i];
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple subtract(const Tuple& a) 
    {
        Tuple result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = my_values[i] - a.get_values()[i];
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple negate() 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i)
        {
            result.get_values()[i] = -my_values[i];
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple multiply(const Tuple a) 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i)
        {
            result.get_values()[i] = my_values[i] * a.get_values()[i];
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple multiply(const int& n) 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()]);
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = my_values[i] * n;
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple divide(const Tuple a) 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = my_values[i] / a.get_values()[i];
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple reciprocal() 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = 1.0 / my_values[i];
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    bool equals(const Object obj) 
    {
        if (dynamic_cast<const Tuple*>(*other) != nullptr)
        {
            const Tuple that = (Tuple) obj;
            if (get_dimension() == that.get_dimension()) 
            {
                bool equals = true;
                for (int i{}; i < my_values.size(); ++i) 
                {
                    equals &= Double.double_to_raw_long_bits(my_values[i]) == Double.double_to_raw_long_bits(that.values[i]);
                }
                return equals;
            }
        }
        return false;
    }

    /** {@inherit_doc} */
    //override
    int hash_code() 
    {
        return  0x34b1a444 + Arrays.hash_code(my_values);
    }

    std::vector<double> get_values() const
    {
        return my_values;
    }

    /** {@inherit_doc} */
    //override
    double get_real() const
    {
        return my_values[0];
    }

    /** {@inherit_doc} */
    //override
    Tuple add(const double& a) 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = my_values[i] + a;
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple subtract(const double& a) 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = my_values[i] - a;
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple multiply(const double& a) 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = my_values[i] * a;
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple divide(const double& a) 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = my_values[i] / a;
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple remainder(const double& a) 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = std::remainder(get_values()[i], a);
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple remainder(const Tuple a) 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = std::remainder(get_values()[i], a.get_values()[i]);
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple abs() 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = std::abs(my_values[i]);
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple ceil() 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = std::ceil(my_values[i]);
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple floor() 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = std::floor(my_values[i]);
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple rint() 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = std::rint(my_values[i]);
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple sign() 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = FastMath.signum(my_values[i]);
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple copy_sign(const Tuple sign) 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = std::copysign(my_values[i], sign.get_values()[i]);
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple copy_sign(const double sign) 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = std::copysign(my_values[i], sign);
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple scalb(const int& n) 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = std::scalbn(my_values[i], n);
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple ulp() 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = FastMath.ulp(my_values[i]);
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple hypot(const Tuple y) 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = std::hypot(my_values[i], y.get_values()[i]);
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple sqrt() 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = std::sqrt(my_values[i]);
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple cbrt() 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = std::cbrt(my_values[i]);
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple root_n(const int& n) 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            if (my_values[i] < 0) 
            {
                result.get_values()[i] = -std::pow(-my_values[i], 1.0 / n);
            }
            else 
            {
                result.get_values()[i] = std::pow(my_values[i], 1.0 / n);
            }
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple pow(const double& p) 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = std::pow(my_values[i], p);
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple pow(const int& n) 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = std::pow(my_values[i], n);
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple pow(const Tuple e) 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = std::pow(my_values[i], e.get_values()[i]);
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple exp() 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = std::exp(my_values[i]);
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple expm1() 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = std::expm1(my_values[i]);
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple log() 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = std::log(my_values[i]);
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple log1p() 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = std::log1p(my_values[i]);
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple log10() 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = std::log10(my_values[i]);
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple cos() 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = std::cos(my_values[i]);
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple sin() 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = std::sin(my_values[i]);
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Field_Sin_Cos<Tuple> sin_cos() 
    {
        const Tuple sin = Tuple(my_field, std::vector<double>(my_values.size()));
        const Tuple cos = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            const Sin_Cos sc = Sin_Cos(my_values[i]);
            sin.get_values()[i] = sc.sin();
            cos.get_values()[i] = sc.cos();
        }
        return Field_Sin_Cos<>(sin, cos);
    }

    /** {@inherit_doc} */
    //override
    Tuple tan() 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = std::tan(my_values[i]);
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple acos() 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = std::acos(my_values[i]);
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple asin() 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = std::asin(my_values[i]);
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple atan() 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = std::atan(my_values[i]);
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple atan2(const Tuple x) 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = std::atan2(my_values[i], x.get_values()[i]);
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple cosh() 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = std::cosh(my_values[i]);
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple sinh() 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = std::sinh(my_values[i]);
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Field_Sinh_Cosh<Tuple> sinh_cosh() 
    {
        const Tuple sinh = Tuple(my_field, std::vector<double>(my_values.size()));
        const Tuple cosh = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            const Sinh_Cosh sch = std::sinh_cosh(my_values[i]);
            sinh.get_values()[i] = sch.sinh();
            cosh.get_values()[i] = sch.cosh();
        }
        return Field_Sinh_Cosh<>(sinh, cosh);
    }

    /** {@inherit_doc} */
    //override
    Tuple tanh() 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = std::tanh(my_values[i]);
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple acosh() 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = std::acosh(my_values[i]);
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple asinh() 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = std::asinh(my_values[i]);
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple atanh() 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = std::atanh(my_values[i]);
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple to_degrees() 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = FastMath.to_degrees(my_values[i]);
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple to_radians() 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = FastMath.to_radians(my_values[i]);
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple linear_combination(const Tuple[] a, const Tuple[] b)
         
        {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        Math_Utils::check_dimension(a.size(), b.size());
        const std::vector<double> a_double = std::vector<double>(a.size()];
        const std::vector<double> b_double = std::vector<double>(b.size()];
        for (int i{}; i < my_values.size(); ++i) 
        {
            for (int j{}; j < a.size(); ++j) 
            {
                a_double[j] = a[j].values[i];
                b_double[j] = b[j].values[i];
            }
            result.get_values()[i] = Math_Arrays::linear_combination(a_double, b_double);
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple linear_combination(const std::vector<double> a, const Tuple[] b)
         
        {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        Math_Utils::check_dimension(a.size(), b.size());
        const std::vector<double> b_double = std::vector<double>(b.size()];
        for (int i{}; i < my_values.size(); ++i) 
        {
            for (int j{}; j < a.size(); ++j) 
            {
                b_double[j] = b[j].values[i];
            }
            result.get_values()[i] = Math_Arrays::linear_combination(a, b_double);
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple linear_combination(const Tuple a1, const Tuple b1, const Tuple a2, const Tuple b2) 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = Math_Arrays::linear_combination(a1.values[i], b1.values[i], a2.values[i], b2.values[i]);
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple linear_combination(const double& a1, const Tuple b1, const double& a2, const Tuple b2) 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = Math_Arrays::linear_combination(a1, b1.get_values()[i], a2, b2.values[i]);
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple linear_combination(const Tuple a1, const Tuple b1, const Tuple a2, const Tuple b2, const Tuple a3, const Tuple b3) 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = Math_Arrays::linear_combination(a1.get_values()[i], b1.get_values()[i], a2.get_values()[i], b2.get_values()[i], a3.get_values()[i], b3.get_values()[i]);
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple linear_combination(const double& a1, const Tuple b1, const double& a2, const Tuple b2, const double& a3, const Tuple b3) 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = Math_Arrays::linear_combination(a1, b1.get_values()[i], a2, b2.get_values()[i], a3, b3.get_values()[i]);
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple linear_combination(const Tuple a1, const Tuple b1, const Tuple a2, const Tuple b2, const Tuple a3, const Tuple b3, const Tuple a4, const Tuple b4) 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = Math_Arrays::linear_combination(a1.get_values()[i], b1.get_values()[i], a2.get_values()[i], b2.get_values()[i], a3.get_values()[i], b3.get_values()[i], a4.get_values()[i], b4.get_values()[i]);
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple linear_combination(const double& a1, const Tuple b1, const double& a2, const Tuple b2, const double& a3, const Tuple b3, const double& a4, const Tuple b4) 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        for (int i{}; i < my_values.size(); ++i) 
        {
            result.get_values()[i] = Math_Arrays::linear_combination(a1, b1.get_values()[i], a2, b2.values[i], a3, b3.values[i], a4, b4.values[i]);
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    Tuple get_pi() 
    {
        auto result = Tuple(my_field, std::vector<double>(my_values.size()));
        Arrays.fill(result.get_values(), std::numbers::pi);
        return result;
    }

    /** Field for {link Tuple} instances.
     */
    static class Tuple_Field : Field<Tuple> 
    {
    private:
        /** Constant function evaluating to 0.0. */
        const Tuple my_zero;

        /** Constant function evaluating to 1.0. */
        const Tuple my_one;

    public:
        /** Simple constructor.
         * @param dimension dimension of the tuple
         */
        Tuple_Field(const int& dimension) 
        {
            const auto zero_data = std::vector<double>(dimension);
            const auto one_data  = std::vector<double>(dimension);
            Arrays.fill(one_data, 1.0);
            my_zero = Tuple(this, zero_data);
            my_one  = Tuple(this, one_data);
        }

        /** {@inherit_doc} */
        //override
        Tuple get_zero() const 
        {
            return my_zero;
        }

        /** {@inherit_doc} */
        //override
        Tuple get_one() const
        {
            return my_one;
        }

        /** {@inherit_doc} */
        //override
        Class<Tuple> get_runtime_class() 
        {
            return Tuple.class;
        }

        /** {@inherit_doc} */
        //override
        bool equals(const Object& other) 
        {
            if (dynamic_cast<const Tuple_Field*>(*other) != nullptr)
            {
                return my_zero.get_dimension() == ((Tuple_Field) other).zero.get_dimension();
            }
            return false;
        }

        /** {@inherit_doc} */
        //override
        int hash_code() 
        {
            return 0x6672493d ^ my_zero.get_dimension();
        }

    }

}


