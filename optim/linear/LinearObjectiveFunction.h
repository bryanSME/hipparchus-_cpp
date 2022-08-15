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
#include <vector>
//import java.io.IOException;
//import java.io.Object_Input_Stream;
//import java.io.Object_Output_Stream;
//import java.io.Serializable;

//import org.hipparchus.analysis.Multivariate_Function;
//import org.hipparchus.linear.Array_Real_Vector;
//import org.hipparchus.linear.Matrix_Utils;
//import org.hipparchus.linear.Real_Vector;
//import org.hipparchus.optim.Optimization_data;
#include "../../core/linear/MatrixUtils.h"
/**
 * An objective function for a linear optimization problem.
 * <p>
 * A linear objective function has one the form:
 * <pre>
 * c<sub>1</sub>x<sub>1</sub> + ... c<sub>n</sub>x<sub>n</sub> + d
 * </pre>
 * The c<sub>i</sub> and d are the coefficients of the equation, * the x<sub>i</sub> are the coordinates of the current point.
 * </p>
 *
 */
class Linear_Objective_Function
    : Multivariate_Function, Optimization_data
{

private:
    /** Coefficients of the linear equation (c<sub>i</sub>). */
    const Real_Vector my_coefficients;
    /** Constant term of the linear equation. */
    const double constant_term;

    /**
     * Serialize the instance.
     * @param oos stream where object should be written
     * @IOException if object cannot be written to stream
     */
    void write_object(Object_Output_Stream oos)
        IOException
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
     * @param coefficients Coefficients for the linear equation being optimized.
     * @param constant_term Constant term of the linear equation.
     */
    Linear_Objective_Function(std::vector<double> coefficients, double constant_term)
    {
        this(new Array_Real_Vector(coefficients), constant_term);
    }

    /**
     * @param coefficients Coefficients for the linear equation being optimized.
     * @param constant_term Constant term of the linear equation.
     */
    Linear_Objective_Function(const Real_Vector& coefficients, const double& constant_term) : my_coefficients{ coefficients }, my_constant_term{ constant_term }
    {}

    /**
     * Gets the coefficients of the linear equation being optimized.
     *
     * @return coefficients of the linear equation being optimized.
     */
    Real_Vector get_coefficients() const
    {
        return my_coefficients;
    }

    /**
     * Gets the constant of the linear equation being optimized.
     *
     * @return constant of the linear equation being optimized.
     */
    double get_constant_term() const
    {
        return constant_term;
    }

    /**
     * Computes the value of the linear equation at the current point.
     *
     * @param point Point at which linear equation must be evaluated.
     * @return the value of the linear equation at the current point.
     */
     //override
    double value(const std::vector<double>& point)
    {
        return value(new Array_Real_Vector(point, false));
    }

    /**
     * Computes the value of the linear equation at the current point.
     *
     * @param point Point at which linear equation must be evaluated.
     * @return the value of the linear equation at the current point.
     */
    double value(const Real_Vector& point)
    {
        return my_coefficients.dot_product(point) + constant_term;
    }

    /** {@inherit_doc} */
    //override
    bool equals(Object other)
    {
        if (this == other)
        {
            return true;
        }
        if (other instanceof Linear_Objective_Function)
        {
            Linear_Objective_Function rhs = (Linear_Objective_Function)other;
            return (constant_term == rhs.constant_term) && my_coefficients.equals(rhs.coefficients);
        }

        return false;
    }

    /** {@inherit_doc} */
    //override
    int hash_code()
    {
        return static_cast<double>(my_constant_term).hash_code() ^ my_coefficients.hash_code();
    }
};