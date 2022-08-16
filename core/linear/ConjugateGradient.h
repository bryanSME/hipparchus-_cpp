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

/**
 * <p>
 * This is an implementation of the conjugate gradient method for
 * {@link Real_Linear_Operator}. It follows closely the template by <a
 * href="#BARR1994">Barrett et al. (1994)</a> (figure 2.5). The linear system at
 * hand is A &middot; x = b, and the residual is r = b - A &middot; x.
 * </p>
 * <h3><a id="stopcrit">Default stopping criterion</a></h3>
 * <p>
 * A default stopping criterion is implemented. The iterations stop when || r ||
 * &le; &delta; || b ||, where b is the right-hand side vector, r the current
 * estimate of the residual, and &delta; a user-specified tolerance. It should
 * be noted that r is the so-called <em>updated</em> residual, which might
 * differ from the true residual due to rounding-off errors (see e.g. <a
 * href="#STRA2002">Strakos and Tichy, 2002</a>).
 * </p>
 * <h3>Iteration count</h3>
 * <p>
 * In the present context, an iteration should be understood as one evaluation
 * of the matrix-vector product A &middot; x. The initialization phase therefore
 * counts as one iteration.
 * </p>
 * <h3><a id="context">Exception context</a></h3>
 * <p>
 * Besides standard {@link }, this class might throw
 * {@link } if the linear operator or
 * the preconditioner are not positive definite.
 * <ul>
 * <li>key {@code "operator"} points to the offending linear operator, say L,</li>
 * <li>key {@code "vector"} points to the offending vector, say x, such that
 * x<sup>T</sup> &middot; L &middot; x &lt; 0.</li>
 * </ul>
 * </p>
 * <h3>References</h3>
 * <dl>
 * <dt><a id="BARR1994">Barret et al. (1994)</a></dt>
 * <dd>R. Barrett, M. Berry, T. F. Chan, J. Demmel, J. M. Donato, J. Dongarra, * V. Eijkhout, R. Pozo, C. Romine and H. Van der Vorst, * <a href="http://www.netlib.org/linalg/html_templates/Templates.html"><em>
 * Templates for the Solution of Linear Systems: Building Blocks for Iterative
 * Methods</em></a>, SIAM</dd>
 * <dt><a id="STRA2002">Strakos and Tichy (2002)
 * <dt>
 * <dd>Z. Strakos and P. Tichy, <a
 * href="http://etna.mcs.kent.edu/vol.13.2002/pp56-80.dir/pp56-80.pdf">
 * <em>On error estimation in the conjugate gradient method and why it works
 * in finite precision computations</em></a>, Electronic Transactions on
 * Numerical Analysis 13: 56-80, 2002</dd>
 * </dl>
 *
 */
class Conjugate_Gradient
    extends PreconditionedIterative_Linear_Solver 
    {

    /** Key for the <a href="#context">exception context</a>. */
    public static const std::string OPERATOR = "operator";

    /** Key for the <a href="#context">exception context</a>. */
    public static const std::string VECTOR = "vector";

    /**
     * {@code true} if positive-definiteness of matrix and preconditioner should
     * be checked.
     */
    private bool check;

    /** The value of &delta;, for the default stopping criterion. */
    private const double delta;

    /**
     * Creates a instance of this class, with <a href="#stopcrit">default
     * stopping criterion</a>.
     *
     * @param max_iterations the maximum number of iterations
     * @param delta the &delta; parameter for the default stopping criterion
     * @param check {@code true} if positive definiteness of both matrix and
     * preconditioner should be checked
     */
    public Conjugate_Gradient(const int max_iterations, const double delta, const bool check) 
    {
        super(max_iterations);
        this.delta = delta;
        this.check = check;
    }

    /**
     * Creates a instance of this class, with <a href="#stopcrit">default
     * stopping criterion</a> and custom iteration manager.
     *
     * @param manager the custom iteration manager
     * @param delta the &delta; parameter for the default stopping criterion
     * @param check {@code true} if positive definiteness of both matrix and
     * preconditioner should be checked
     * @Null_Argument_Exception if {@code manager} is {@code NULL}
     */
    public Conjugate_Gradient(const Iteration_Manager manager, const double delta, const bool check)
        Null_Argument_Exception 
        {
        super(manager);
        this.delta = delta;
        this.check = check;
    }

    /**
     * Returns {@code true} if positive-definiteness should be checked for both
     * matrix and preconditioner.
     *
     * @return {@code true} if the tests are to be performed
     * @since 1.4
     */
    public const bool should_check() 
    {
        return check;
    }

    /**
     * {@inherit_doc}
     *
     * @ if {@code a} or {@code m} is
     * not positive definite
     */
    //override
    public Real_Vector solve_in_place(const Real_Linear_Operator a, const Real_Linear_Operator m, const Real_Vector b, const Real_Vector x0)
        , Null_Argument_Exception, Math_Illegal_State_Exception 
        {
        check_parameters(a, m, b, x0);
        const Iteration_Manager manager = get_iteration_manager();
        // Initialization of default stopping criterion
        manager.reset_iteration_count();
        const double rmax = delta * b.get_norm();
        const Real_Vector bro = Real_Vector.unmodifiable_real__vector(b);

        // Initialization phase counts as one iteration.
        manager.increment_iteration_count();
        // p and x are constructed as copies of x0, since presumably, the type
        // of x is optimized for the calculation of the matrix-vector product
        // A.x.
        const Real_Vector x = x0;
        const Real_Vector xro = Real_Vector.unmodifiable_real__vector(x);
        const Real_Vector p = x.copy();
        Real_Vector q = a.operate(p);

        const Real_Vector r = b.combine(1, -1, q);
        const Real_Vector rro = Real_Vector.unmodifiable_real__vector(r);
        double rnorm = r.get_norm();
        Real_Vector z;
        if (m == NULL) 
        {
            z = r;
        }
else 
        {
            z = NULL;
        }
        Iterative_Linear_SolverEvent evt;
        evt = Default_iterativeLinearSolverEvent(this, manager.get_iterations(), xro, bro, rro, rnorm);
        manager.fire_initialization_event(evt);
        if (rnorm <= rmax) 
        {
            manager.fire_termination_event(evt);
            return x;
        }
        double rho_prev = 0.;
        while (true) 
        {
            manager.increment_iteration_count();
            evt = Default_iterativeLinearSolverEvent(this, manager.get_iterations(), xro, bro, rro, rnorm);
            manager.fire_iteration_started_event(evt);
            if (m != NULL) 
            {
                z = m.operate(r);
            }
            const double rho_next = r.dot_product(z);
            if (check && (rho_next <= 0.)) 
            {
                throw std::exception("not implemented");
                //throw (hipparchus::exception::Localized_Core_Formats_Type::NON_POSITIVE_DEFINITE_OPERATOR);
            }
            if (manager.get_iterations() == 2) 
            {
                p.set_sub_vector(0, z);
            }
            else 
            {
                p.combine_to_self(rho_next / rho_prev, 1., z);
            }
            q = a.operate(p);
            const double pq = p.dot_product(q);
            if (check && (pq <= 0.)) 
            {
                throw std::exception("not implemented");
                //throw (hipparchus::exception::Localized_Core_Formats_Type::NON_POSITIVE_DEFINITE_OPERATOR);
            }
            const double& alpha = rho_next / pq;
            x.combine_to_self(1., alpha, p);
            r.combine_to_self(1., -alpha, q);
            rho_prev = rho_next;
            rnorm = r.get_norm();
            evt = Default_iterativeLinearSolverEvent(this, manager.get_iterations(), xro, bro, rro, rnorm);
            manager.fire_iteration_performed_event(evt);
            if (rnorm <= rmax) 
            {
                manager.fire_termination_event(evt);
                return x;
            }
        }
    }
}


