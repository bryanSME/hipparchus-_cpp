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
//package org.hipparchus.optim.linear;

//import java.io.IOException;
//import java.io.Object_Input_Stream;
//import java.io.Object_Output_Stream;
//import java.io.Serializable;

//import org.hipparchus.linear.Array_Real_Vector;
//import org.hipparchus.linear.Matrix_Utils;
//import org.hipparchus.linear.Real_Vector;
#include "../../core/linear/MatrixUtils.h"
/**
 * A linear constraint for a linear optimization problem.
 * <p>
 * A linear constraint has one of the forms:
 * <ul>
 *   <li>c<sub>1</sub>x<sub>1</sub> + ... c<sub>n</sub>x<sub>n</sub> = v</li>
 *   <li>c<sub>1</sub>x<sub>1</sub> + ... c<sub>n</sub>x<sub>n</sub> &lt;= v</li>
 *   <li>c<sub>1</sub>x<sub>1</sub> + ... c<sub>n</sub>x<sub>n</sub> &gt;= v</li>
 *   <li>l<sub>1</sub>x<sub>1</sub> + ... l<sub>n</sub>x<sub>n</sub> + l<sub>cst</sub> =
 *       r<sub>1</sub>x<sub>1</sub> + ... r<sub>n</sub>x<sub>n</sub> + r<sub>cst</sub></li>
 *   <li>l<sub>1</sub>x<sub>1</sub> + ... l<sub>n</sub>x<sub>n</sub> + l<sub>cst</sub> &lt;=
 *       r<sub>1</sub>x<sub>1</sub> + ... r<sub>n</sub>x<sub>n</sub> + r<sub>cst</sub></li>
 *   <li>l<sub>1</sub>x<sub>1</sub> + ... l<sub>n</sub>x<sub>n</sub> + l<sub>cst</sub> &gt;=
 *       r<sub>1</sub>x<sub>1</sub> + ... r<sub>n</sub>x<sub>n</sub> + r<sub>cst</sub></li>
 * </ul>
 * The c<sub>i</sub>, l<sub>i</sub> or r<sub>i</sub> are the coefficients of the constraints, the x<sub>i</sub>
 * are the coordinates of the current point and v is the value of the constraint.
 * </p>
 *
 */
class Linear_Constraint  
{
private:
    /** Coefficients of the constraint (left hand side). */
    const Real_Vector my_coefficients;
    /** Relationship between left and right hand sides (=, &lt;=, &gt;=). */
    const Relationship my_relationship;
    /** Value of the constraint (right hand side). */
    const double my_value;

    /**
     * Serialize the instance.
     * @param oos stream where object should be written
     * @IOException if object cannot be written to stream
     */
    void write_object(Object_Output_Stream oos)
    {
        oos.default_write_object();
        Matrix_Utils.serialize_real__vector(coefficients, oos);
    }

    /**
     * Deserialize the instance.
     * @param ois stream from which the object should be read
     * @Class_Not_Found_Exception if a class in the stream cannot be found
     * @IOException if object cannot be read from the stream
     */
    void read_object(Object_Input_Stream ois)
    {
      ois.default_read_object();
      Matrix_Utils.deserialize_real__vector(this, "coefficients", ois);
    }

public:
    /**
     * Build a constraint involving a single linear equation.
     * <p>
     * A linear constraint with a single linear equation has one of the forms:
     * <ul>
     *   <li>c<sub>1</sub>x<sub>1</sub> + ... c<sub>n</sub>x<sub>n</sub> = v</li>
     *   <li>c<sub>1</sub>x<sub>1</sub> + ... c<sub>n</sub>x<sub>n</sub> &lt;= v</li>
     *   <li>c<sub>1</sub>x<sub>1</sub> + ... c<sub>n</sub>x<sub>n</sub> &gt;= v</li>
     * </ul>
     * </p>
     * @param coefficients The coefficients of the constraint (left hand side)
     * @param relationship The type of (in)equality used in the constraint
     * @param value The value of the constraint (right hand side)
     */
    Linear_Constraint(const std::vector<double> coefficients, const Relationship relationship, const double& value) 
    {
        this(new Array_Real_Vector(coefficients), relationship, value);
    }

    /**
     * Build a constraint involving a single linear equation.
     * <p>
     * A linear constraint with a single linear equation has one of the forms:
     * <ul>
     *   <li>c<sub>1</sub>x<sub>1</sub> + ... c<sub>n</sub>x<sub>n</sub> = v</li>
     *   <li>c<sub>1</sub>x<sub>1</sub> + ... c<sub>n</sub>x<sub>n</sub> &lt;= v</li>
     *   <li>c<sub>1</sub>x<sub>1</sub> + ... c<sub>n</sub>x<sub>n</sub> &gt;= v</li>
     * </ul>
     * </p>
     * @param coefficients The coefficients of the constraint (left hand side)
     * @param relationship The type of (in)equality used in the constraint
     * @param value The value of the constraint (right hand side)
     */
    Linear_Constraint(const Real_Vector& coefficients, const Relationship& relationship, const double& value) 
        : my_coefficients{ coefficients }, my_relationship{ relationship }, my_value{ value }
    {}

    /**
     * Build a constraint involving two linear equations.
     * <p>
     * A linear constraint with two linear equation has one of the forms:
     * <ul>
     *   <li>l<sub>1</sub>x<sub>1</sub> + ... l<sub>n</sub>x<sub>n</sub> + l<sub>cst</sub> =
     *       r<sub>1</sub>x<sub>1</sub> + ... r<sub>n</sub>x<sub>n</sub> + r<sub>cst</sub></li>
     *   <li>l<sub>1</sub>x<sub>1</sub> + ... l<sub>n</sub>x<sub>n</sub> + l<sub>cst</sub> &lt;=
     *       r<sub>1</sub>x<sub>1</sub> + ... r<sub>n</sub>x<sub>n</sub> + r<sub>cst</sub></li>
     *   <li>l<sub>1</sub>x<sub>1</sub> + ... l<sub>n</sub>x<sub>n</sub> + l<sub>cst</sub> &gt;=
     *       r<sub>1</sub>x<sub>1</sub> + ... r<sub>n</sub>x<sub>n</sub> + r<sub>cst</sub></li>
     * </ul>
     * </p>
     * @param lhs_coefficients The coefficients of the linear expression on the left hand side of the constraint
     * @param lhs_constant The constant term of the linear expression on the left hand side of the constraint
     * @param relationship The type of (in)equality used in the constraint
     * @param rhs_coefficients The coefficients of the linear expression on the right hand side of the constraint
     * @param rhs_constant The constant term of the linear expression on the right hand side of the constraint
     */
    Linear_Constraint(const std::vector<double>& lhs_coefficients, const double& lhs_constant, const Relationship& relationship, const std::vector<double> rhs_coefficients, const double rhs_constant) 
    {
        auto sub = std::vector<double>(lhs_coefficients.size()];
        for (int i{}; i < sub.size(); ++i) 
        {
            sub[i] = lhs_coefficients[i] - rhs_coefficients[i];
        }
        my_coefficients = Array_Real_Vector(sub, false);
        my_relationship = relationship;
        my_value        = rhs_constant - lhs_constant;
    }

    /**
     * Build a constraint involving two linear equations.
     * <p>
     * A linear constraint with two linear equation has one of the forms:
     * <ul>
     *   <li>l<sub>1</sub>x<sub>1</sub> + ... l<sub>n</sub>x<sub>n</sub> + l<sub>cst</sub> =
     *       r<sub>1</sub>x<sub>1</sub> + ... r<sub>n</sub>x<sub>n</sub> + r<sub>cst</sub></li>
     *   <li>l<sub>1</sub>x<sub>1</sub> + ... l<sub>n</sub>x<sub>n</sub> + l<sub>cst</sub> &lt;=
     *       r<sub>1</sub>x<sub>1</sub> + ... r<sub>n</sub>x<sub>n</sub> + r<sub>cst</sub></li>
     *   <li>l<sub>1</sub>x<sub>1</sub> + ... l<sub>n</sub>x<sub>n</sub> + l<sub>cst</sub> &gt;=
     *       r<sub>1</sub>x<sub>1</sub> + ... r<sub>n</sub>x<sub>n</sub> + r<sub>cst</sub></li>
     * </ul>
     * </p>
     * @param lhs_coefficients The coefficients of the linear expression on the left hand side of the constraint
     * @param lhs_constant The constant term of the linear expression on the left hand side of the constraint
     * @param relationship The type of (in)equality used in the constraint
     * @param rhs_coefficients The coefficients of the linear expression on the right hand side of the constraint
     * @param rhs_constant The constant term of the linear expression on the right hand side of the constraint
     */
    Linear_Constraint(const Real_Vector& lhs_coefficients, const double& lhs_constant, const Relationship& relationship, const Real_Vector& rhs_coefficients, const double& rhs_constant) 
        : my_coefficients{ lhs_coefficients.subtract(rhs_coefficients) }, my_relationship{ relationship }, my_value{ rhs_constant - lhs_constant }
    {}

    /**
     * Gets the coefficients of the constraint (left hand side).
     *
     * @return the coefficients of the constraint (left hand side).
     */
    Real_Vector get_coefficients() const
    {
        return my_coefficients;
    }

    /**
     * Gets the relationship between left and right hand sides.
     *
     * @return the relationship between left and right hand sides.
     */
    Relationship get_relationship() const 
    {
        return my_relationship;
    }

    /**
     * Gets the value of the constraint (right hand side).
     *
     * @return the value of the constraint (right hand side).
     */
    double get_value() const
    {
        return my_value;
    }

    /** {@inherit_doc} */
    //override
    bool equals(Object other) 
    {
        if (this == other) 
        {
            return true;
        }
        if (other instanceof Linear_Constraint) 
        {
            Linear_Constraint rhs = (Linear_Constraint) other;
            return relationship == rhs.relationship &&
                value == rhs.value &&
                coefficients.equals(rhs.coefficients);
        }
        return false;
    }

    /** {@inherit_doc} */
    //override
    int hash_code() 
    {
        return relationship.hash_code() ^
            static_cast<double>(value).hash_code() ^
            coefficients.hash_code();
    }
};