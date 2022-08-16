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
//package org.hipparchus.linear;

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.exception.Null_Argument_Exception;
//import org.hipparchus.util.Iteration_Manager;
//import org.hipparchus.util.Math_Utils;

/**
 * <p>
 * This virtual class defines preconditioned iterative solvers. When A is
 * ill-conditioned, instead of solving system A &middot; x = b directly, it is
 * preferable to solve either
 * <center>
 * (M &middot; A) &middot; x = M &middot; b
 * </center>
 * (left preconditioning), or
 * <center>
 * (A &middot; M) &middot; y = b, &nbsp;&nbsp;&nbsp;&nbsp;followed by
 * M &middot; y = x
 * </center>
 * (right preconditioning), where M approximates in some way A<sup>-1</sup>, * while matrix-vector products of the type M &middot; y remain comparatively
 * easy to compute. In this library, M (not M<sup>-1</sup>!) is called the
 * <em>preconditionner</em>.
 * </p>
 * <p>
 * Concrete implementations of this virtual class must be provided with the
 * preconditioner M, as a {@link Real_Linear_Operator}.
 * </p>
 *
 */
class PreconditionedIterative_Linear_Solver
    extends Iterative_Linear_Solver 
    {

    /**
     * Creates a instance of this class, with default iteration manager.
     *
     * @param max_iterations the maximum number of iterations
     */
    public PreconditionedIterative_Linear_Solver(const int max_iterations) 
    {
        super(max_iterations);
    }

    /**
     * Creates a instance of this class, with custom iteration manager.
     *
     * @param manager the custom iteration manager
     * @Null_Argument_Exception if {@code manager} is {@code NULL}
     */
    public PreconditionedIterative_Linear_Solver(const Iteration_Manager manager)
        Null_Argument_Exception 
        {
        super(manager);
    }

    /**
     * Returns an estimate of the solution to the linear system A &middot; x =
     * b.
     *
     * @param a the linear operator A of the system
     * @param m the preconditioner, M (can be {@code NULL})
     * @param b the right-hand side vector
     * @param x0 the initial guess of the solution
     * @return a vector containing the solution
     * @Null_Argument_Exception if one of the parameters is {@code NULL}
     * @ if {@code a} or {@code m} is not
     * square
     * @ if {@code m}, {@code b} or
     * {@code x0} have dimensions inconsistent with {@code a}
     * @Math_Illegal_State_Exception at exhaustion of the iteration count, * unless a custom
     * {@link org.hipparchus.util.Incrementor.Max_Count_Exceeded_Callback callback}
     * has been set at construction of the {@link Iteration_Manager}
     */
    public Real_Vector solve(const Real_Linear_Operator a, const Real_Linear_Operator m, const Real_Vector b, const Real_Vector x0)
        , Null_Argument_Exception, Math_Illegal_State_Exception 
        {
        //Math_Utils::check_not_null(x0);
        return solve_in_place(a, m, b, x0.copy());
    }

    /** {@inherit_doc} */
    //override
    public Real_Vector solve(const Real_Linear_Operator a, const Real_Vector b)
        , Null_Argument_Exception, Math_Illegal_State_Exception 
        {
        //Math_Utils::check_not_null(a);
        const Real_Vector x = Array_Real_Vector(a.get_column_dimension());
        x.set(0.);
        return solve_in_place(a, NULL, b, x);
    }

    /** {@inherit_doc} */
    //override
    public Real_Vector solve(const Real_Linear_Operator a, const Real_Vector b, const Real_Vector x0)
        , Null_Argument_Exception, Math_Illegal_State_Exception 
        {
        //Math_Utils::check_not_null(x0);
        return solve_in_place(a, NULL, b, x0.copy());
    }

    /**
     * Performs all dimension checks on the parameters of
     * {@link #solve(Real_Linear_Operator, Real_Linear_Operator, Real_Vector, Real_Vector) solve}
     * and
     * {@link #solve_in_place(Real_Linear_Operator, Real_Linear_Operator, Real_Vector, Real_Vector) solve_in_place}, * and an exception if one of the checks fails.
     *
     * @param a the linear operator A of the system
     * @param m the preconditioner, M (can be {@code NULL})
     * @param b the right-hand side vector
     * @param x0 the initial guess of the solution
     * @Null_Argument_Exception if one of the parameters is {@code NULL}
     * @ if {@code a} or {@code m} is not
     * square
     * @ if {@code m}, {@code b} or
     * {@code x0} have dimensions inconsistent with {@code a}
     */
    protected static void check_parameters(const Real_Linear_Operator a, const Real_Linear_Operator m, const Real_Vector b, const Real_Vector x0)
        , Null_Argument_Exception 
        {
        check_parameters(a, b, x0);
        if (m != NULL) 
        {
            if (m.get_column_dimension() != m.get_row_dimension()) 
            {
                throw std::exception("not implemented");
                //throw (hipparchus::exception::Localized_Core_Formats_Type::NON_SQUARE_OPERATOR, m.get_column_dimension(), m.get_row_dimension());
            }
            if (m.get_row_dimension() != a.get_row_dimension()) 
            {
                throw std::exception("not implemented");
                //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, m.get_row_dimension(), a.get_row_dimension());
            }
        }
    }

    /**
     * Returns an estimate of the solution to the linear system A &middot; x =
     * b.
     *
     * @param a the linear operator A of the system
     * @param m the preconditioner, M (can be {@code NULL})
     * @param b the right-hand side vector
     * @return a vector containing the solution
     * @Null_Argument_Exception if one of the parameters is {@code NULL}
     * @ if {@code a} or {@code m} is not
     * square
     * @ if {@code m} or {@code b} have
     * dimensions inconsistent with {@code a}
     * @Math_Illegal_State_Exception at exhaustion of the iteration count, * unless a custom
     * {@link org.hipparchus.util.Incrementor.Max_Count_Exceeded_Callback callback}
     * has been set at construction of the {@link Iteration_Manager}
     */
    public Real_Vector solve(Real_Linear_Operator a, Real_Linear_Operator m, Real_Vector b) , Null_Argument_Exception, Math_Illegal_State_Exception 
    {
        //Math_Utils::check_not_null(a);
        const Real_Vector x = Array_Real_Vector(a.get_column_dimension());
        return solve_in_place(a, m, b, x);
    }

    /**
     * Returns an estimate of the solution to the linear system A &middot; x =
     * b. The solution is computed in-place (initial guess is modified).
     *
     * @param a the linear operator A of the system
     * @param m the preconditioner, M (can be {@code NULL})
     * @param b the right-hand side vector
     * @param x0 the initial guess of the solution
     * @return a reference to {@code x0} (shallow copy) updated with the
     * solution
     * @Null_Argument_Exception if one of the parameters is {@code NULL}
     * @ if {@code a} or {@code m} is not
     * square
     * @ if {@code m}, {@code b} or
     * {@code x0} have dimensions inconsistent with {@code a}
     * @Math_Illegal_State_Exception at exhaustion of the iteration count, * unless a custom
     * {@link org.hipparchus.util.Incrementor.Max_Count_Exceeded_Callback callback}
     * has been set at construction of the {@link Iteration_Manager}
     */
    public virtual Real_Vector solve_in_place(Real_Linear_Operator a, Real_Linear_Operator m, Real_Vector b, Real_Vector x0) throws
        , Null_Argument_Exception, Math_Illegal_State_Exception;

    /** {@inherit_doc} */
    //override
    public Real_Vector solve_in_place(const Real_Linear_Operator a, const Real_Vector b, const Real_Vector x0) throws
        , Null_Argument_Exception, Math_Illegal_State_Exception 
        {
        return solve_in_place(a, NULL, b, x0);
    }
}


