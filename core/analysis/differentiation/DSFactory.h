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

//import java.io.Serializable;

//import org.hipparchus.Field;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.util.FastMath;
#include <numbers>
#include "DSCompiler.h"

/** Factory for {@link Derivative_Structure}.
 * <p>This class is a factory for {@link Derivative_Structure} instances.</p>
 * <p>Instances of this class are guaranteed to be immutable.</p>
 * @see Derivative_Structure
 * @since 1.1
 */
class DS_Factory  
{
private:
    /** Compiler for the current dimensions. */
    const DS_Compiler my_compiler;

    /** Field the {@link Derivative_Structure} instances belong to. */
    const DS_Field my_derivative_field;

public:
    /** Simple constructor.
     * @param parameters number of free parameters
     * @param order derivation order
     */
    DS_Factory(const int& parameters, const int& order)
        :
        my_compiler{ DS_Compiler.get_compiler(parameters, order) },
        my_derivative_field{ DS_Field(constant(0.0), constant(1.0), constant(std::numbers::pi)) }
    {
    }

    /** Get the {@link Field} the {@link Derivative_Structure} instances belong to.
     * @return {@link Field} the {@link Derivative_Structure} instances belong to
     */
    DS_Field get_derivative_field() const 
    {
        return my_derivative_field;
    }

    /** Build a {@link Derivative_Structure} representing a constant value.
     * @param value value of the constant
     * @return a {@link Derivative_Structure} representing a constant value
     */
    Derivative_Structure constant(const double& value) 
    {
        auto ds = Derivative_Structure(*this);
        ds.set_derivative_component(0, value);
        return ds;
    }

    /** Build a {@link Derivative_Structure} representing a variable.
     * <p>Instances built using this method are considered
     * to be the free variables with respect to which differentials
     * are computed. As such, their differential with respect to
     * themselves is +1.</p>
     * @param index index of the variable (from 0 to
     * {@link #get_compiler()}.{@link DS_Compiler#get_free_parameters() get_free_parameters()} - 1)
     * @param value value of the variable
     * @exception  if index if greater or
     * equal to {@link #get_compiler()}.{@link DS_Compiler#get_free_parameters() get_free_parameters()}.
     * @return a {@link Derivative_Structure} representing a variable
     */
    Derivative_Structure variable(const int& index, const double& value)
    {
        if (index >= get_compiler().get_free_parameters()) 
        {
            throw std::exception("not implmented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_LARGE_BOUND_EXCLUDED, index, get_compiler().get_free_parameters());
        }

        const auto ds = Derivative_Structure(*this);
        ds.set_derivative_component(0, value);

        if (get_compiler().get_order() > 0) 
        {
            // the derivative of the variable with respect to itself is 1.
            ds.set_derivative_component(DS_Compiler.get_compiler(index, get_compiler().get_order()).get_size(), 1.0);
        }

        return ds;

    }

    /** Build a {@link Derivative_Structure} from all its derivatives.
     * @param derivatives derivatives sorted according to
     * {@link DS_Compiler#get_partial_derivative_index(int...)}
     * @return a {@link Derivative_Structure} with specified derivatives
     * @exception  if derivatives array does not match the
     * {@link DS_Compiler#get_size() size} expected by the compiler
     * @exception  if order is too large
     * @see Derivative_Structure#get_all_derivatives()
     */
   //@Safe_Varargs
    const Derivative_Structure build(const double ... derivatives)
    {
        if (derivatives.size() != compiler.get_size()) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, derivatives.size(), compiler.get_size());
        }

        return Derivative_Structure(this, derivatives);

    }

    /** Build a {@link Derivative_Structure} with an uninitialized array.
     * <p>This method is intended only for Derivative_Structure internal use.</p>
     * @return a {@link Derivative_Structure} with an uninitialized array
     */
    Derivative_Structure build() 
    {
        return Derivative_Structure(this);
    }

    /** Get the compiler for the current dimensions.
     * @return compiler for the current dimensions
     */
    DS_Compiler get_compiler() 
    {
        return compiler;
    }

    /** Check rules set compatibility.
     * @param factory other factory field to check against instance
     * @exception  if number of free parameters or orders are inconsistent
     */
    void check_compatibility(const DS_Factory& factory)  
    {
        compiler.check_compatibility(factory.compiler);
    }

    /**
     * Replace the instance with a data transfer object for serialization.
     * @return data transfer object that will be serialized
     */
    private Object write_replace() 
    {
        return Data_Transfer_Object(compiler.get_free_parameters(), compiler.get_order());
    }

    /** Internal class used only for serialization. */
    private static class Data_Transfer_Object
    {

        
        20161222L;

        /** Number of variables.
         * @serial
         */
        private const int variables;

        /** Derivation order.
         * @serial
         */
        private const int order;

        /** Simple constructor.
         * @param variables number of variables
         * @param order derivation order
         */
        Data_Transfer_Object(const int variables, const int order) : my_variables{ variables }, my_order{ order } {};

        /** Replace the deserialized data transfer object with a {@link DS_Factory}.
         * @return replacement {@link DS_Factory}
         */
        private Object read_resolve() 
        {
            return DS_Factory(my_variables, my_order);
        }

    }

    /** Field for {link Derivative_Structure} instances.
     */
    public static class DS_Field : Field<Derivative_Structure> 
    {

        /** Constant function evaluating to 0.0. */
        private const Derivative_Structure zero;

        /** Constant function evaluating to 1.0. */
        private const Derivative_Structure one;

        /** Constant function evaluating to π. */
        private const Derivative_Structure pi;

        /** Simple constructor.
         * @param zero constant function evaluating to 0.0
         * @param one constant function evaluating to 1.0
         * @param pi constant function evaluating to π
         */
        DS_Field(const Derivative_Structure zero, const Derivative_Structure one, const Derivative_Structure pi) 
        {
            this.zero = zero;
            this.one  = one;
            this.pi   = pi;
        }

        /** {@inherit_doc} */
        //override
        public Derivative_Structure get_zero() 
        {
            return zero;
        }

        /** {@inherit_doc} */
        //override
        public Derivative_Structure get_one() 
        {
            return one;
        }

        /** Get the Archimedes constant π.
         * <p>
         * Archimedes constant is the ratio of a circle's circumference to its diameter.
         * </p>
         * @return Archimedes constant π
         * @since 2.0
         */
        public Derivative_Structure get_pi() 
        {
            return pi;
        }

        /** {@inherit_doc} */
        //override
        public Class<Derivative_Structure> get_runtime_class() 
        {
            return Derivative_Structure.class;
        }

        /** {@inherit_doc} */
        //override
        public bool equals(const Object& other) 
        {
            if (this == other) 
            {
                return true;
            }
            if (dynamic_cast<const DS_Field*>(*other) != nullptr)
            {
                DS_Factory lhs_factory = zero.get_factory();
                DS_Factory rhs_factory = ((DS_Field) other).zero.get_factory();
                return lhs_factory.compiler == rhs_factory.compiler;
            }
            return false;
        }

        /** {@inherit_doc} */
        //override
        public int hash_code() 
        {
            const DS_Compiler compiler = zero.get_factory().get_compiler();
            return 0x9943b886 ^ (compiler.get_free_parameters() << 16 & compiler.get_order());
        }

    }

};