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
//import org.hipparchus.exception.;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Iteration_Manager;
//import org.hipparchus.util.Math_Utils;

/**
 * <p>
 * Implementation of the SYMMLQ iterative linear solver proposed by <a
 * href="#PAIG1975">Paige and Saunders (1975)</a>. This implementation is
 * largely based on the FORTRAN code by Pr. Michael A. Saunders, available <a
 * href="http://www.stanford.edu/group/SOL/software/symmlq/f77/">here</a>.
 * </p>
 * <p>
 * SYMMLQ is designed to solve the system of linear equations A &middot; x = b
 * where A is an n &times; n self-adjoint linear operator (defined as a
 * {@link Real_Linear_Operator}), and b is a given vector. The operator A is not
 * required to be positive definite. If A is known to be definite, the method of
 * conjugate gradients might be preferred, since it will require about the same
 * number of iterations as SYMMLQ but slightly less work per iteration.
 * </p>
 * <p>
 * SYMMLQ is designed to solve the system (A - shift &middot; I) &middot; x = b, * where shift is a specified scalar value. If shift and b are suitably chosen, * the computed vector x may approximate an (unnormalized) eigenvector of A, as
 * in the methods of inverse iteration and/or Rayleigh-quotient iteration.
 * Again, the linear operator (A - shift &middot; I) need not be positive
 * definite (but <em>must</em> be self-adjoint). The work per iteration is very
 * slightly less if shift = 0.
 * </p>
 * <h3>Preconditioning</h3>
 * <p>
 * Preconditioning may reduce the number of iterations required. The solver may
 * be provided with a positive definite preconditioner
 * M = P<sup>T</sup> &middot; P
 * that is known to approximate
 * (A - shift &middot; I)<sup>-1</sup> in some sense, where matrix-vector
 * products of the form M &middot; y = x can be computed efficiently. Then
 * SYMMLQ will implicitly solve the system of equations
 * P &middot; (A - shift &middot; I) &middot; P<sup>T</sup> &middot;
 * x<sub>hat</sub> = P &middot; b, i.e.
 * A<sub>hat</sub> &middot; x<sub>hat</sub> = b<sub>hat</sub>, * where
 * A<sub>hat</sub> = P &middot; (A - shift &middot; I) &middot; P<sup>T</sup>, * b<sub>hat</sub> = P &middot; b, * and return the solution
 * x = P<sup>T</sup> &middot; x<sub>hat</sub>.
 * The associated residual is
 * r<sub>hat</sub> = b<sub>hat</sub> - A<sub>hat</sub> &middot; x<sub>hat</sub>
 *                 = P &middot; [b - (A - shift &middot; I) &middot; x]
 *                 = P &middot; r.
 * </p>
 * <p>
 * In the case of preconditioning, the {@link Iterative_Linear_SolverEvent}s that
 * this solver fires are such that
 * {@link Iterative_Linear_SolverEvent#get_norm_of_residual()} returns the norm of
 * the <em>preconditioned</em>, updated residual, ||P &middot; r||, not the norm
 * of the <em>true</em> residual ||r||.
 * </p>
 * <h3><a id="stopcrit">Default stopping criterion</a></h3>
 * <p>
 * A default stopping criterion is implemented. The iterations stop when || rhat
 * || &le; &delta; || Ahat || || xhat ||, where xhat is the current estimate of
 * the solution of the transformed system, rhat the current estimate of the
 * corresponding residual, and &delta; a user-specified tolerance.
 * </p>
 * <h3>Iteration count</h3>
 * <p>
 * In the present context, an iteration should be understood as one evaluation
 * of the matrix-vector product A &middot; x. The initialization phase therefore
 * counts as one iteration. If the user requires checks on the symmetry of A, * this entails one further matrix-vector product in the initial phase. This
 * further product is <em>not</em> accounted for in the iteration count. In
 * other words, the number of iterations required to reach convergence will be
 * identical, whether checks have been required or not.
 * </p>
 * <p>
 * The present definition of the iteration count differs from that adopted in
 * the original FOTRAN code, where the initialization phase was <em>not</em>
 * taken into account.
 * </p>
 * <h3><a id="initguess">Initial guess of the solution</a></h3>
 * <p>
 * The {@code x} parameter in
 * <ul>
 * <li>{@link #solve(Real_Linear_Operator, Real_Vector, Real_Vector)},</li>
 * <li>{@link #solve(Real_Linear_Operator, Real_Linear_Operator, Real_Vector, Real_Vector)}},</li>
 * <li>{@link #solve_in_place(Real_Linear_Operator, Real_Vector, Real_Vector)},</li>
 * <li>{@link #solve_in_place(Real_Linear_Operator, Real_Linear_Operator, Real_Vector, Real_Vector)},</li>
 * <li>{@link #solve_in_place(Real_Linear_Operator, Real_Linear_Operator, Real_Vector, Real_Vector, bool, double)},</li>
 * </ul>
 * should not be considered as an initial guess, as it is set to zero in the
 * initial phase. If x<sub>0</sub> is known to be a good approximation to x, one
 * should compute r<sub>0</sub> = b - A &middot; x, solve A &middot; dx = r0, * and set x = x<sub>0</sub> + dx.
 * </p>
 * <h3><a id="context">Exception context</a></h3>
 * <p>
 * Besides standard {@link }, this class might throw
 * {@link } if the linear operator or the
 * preconditioner are not symmetric.
 * <ul>
 * <li>key {@code "operator"} points to the offending linear operator, say L,</li>
 * <li>key {@code "vector1"} points to the first offending vector, say x, * <li>key {@code "vector2"} points to the second offending vector, say y, such
 * that x<sup>T</sup> &middot; L &middot; y &ne; y<sup>T</sup> &middot; L
 * &middot; x (within a certain accuracy).</li>
 * </ul>
 * </p>
 * <p>
 * {@link } might also be thrown in case the
 * preconditioner is not positive definite.
 * </p>
 * <h3>References</h3>
 * <dl>
 * <dt><a id="PAIG1975">Paige and Saunders (1975)</a></dt>
 * <dd>C. C. Paige and M. A. Saunders, <a
 * href="http://www.stanford.edu/group/SOL/software/symmlq/PS75.pdf"><em>
 * Solution of Sparse Indefinite Systems of Linear Equations</em></a>, SIAM
 * Journal on Numerical Analysis 12(4): 617-629, 1975</dd>
 * </dl>
 *
 */
class Symm_LQ
    extends PreconditionedIterative_Linear_Solver 
    {

    /** {@code true} if symmetry of matrix and conditioner must be checked. */
    private const bool check;

    /**
     * The value of the custom tolerance &delta; for the default stopping
     * criterion.
     */
    private const double delta;

    /*
     * IMPLEMENTATION NOTES
     * --------------------
     * The implementation follows as closely as possible the notations of Paige
     * and Saunders (1975). Attention must be paid to the fact that some
     * quantities which are relevant to iteration k can only be computed in
     * iteration (k+1). Therefore, minute attention must be paid to the index of
     * each state variable of this algorithm.
     *
     * 1. Preconditioning
     *    ---------------
     * The Lanczos iterations associated with Ahat and bhat read
     *   beta[1] = ||P * b||
     *   v[1] = P * b / beta[1]
     *   beta[k+1] * v[k+1] = Ahat * v[k] - alpha[k] * v[k] - beta[k] * v[k-1]
     *                      = P * (A - shift * I) * P' * v[k] - alpha[k] * v[k]
     *                        - beta[k] * v[k-1]
     * Multiplying both sides by P', we get
     *   beta[k+1] * (P' * v)[k+1] = M * (A - shift * I) * (P' * v)[k]
     *                               - alpha[k] * (P' * v)[k]
     *                               - beta[k] * (P' * v[k-1]), * and
     *   alpha[k+1] = v[k+1]' * Ahat * v[k+1]
     *              = v[k+1]' * P * (A - shift * I) * P' * v[k+1]
     *              = (P' * v)[k+1]' * (A - shift * I) * (P' * v)[k+1].
     *
     * In other words, the Lanczos iterations are unchanged, except for the fact
     * that we really compute (P' * v) instead of v. It can easily be checked
     * that all other formulas are unchanged. It must be noted that P is never
     * explicitly used, only matrix-vector products involving are invoked.
     *
     * 2. Accounting for the shift parameter
     *    ----------------------------------
     * Is trivial: each time A.operate(x) is invoked, one must subtract shift * x
     * to the result.
     *
     * 3. Accounting for the goodb flag
     *    -----------------------------
     * When goodb is set to true, the component of x_l along b is computed
     * separately. From Paige and Saunders (1975), equation (5.9), we have
     *   wbar[k+1] = s[k] * wbar[k] - c[k] * v[k+1], *   wbar[1] = v[1].
     * Introducing wbar2[k] = wbar[k] - s[1] * ... * s[k-1] * v[1], it can
     * easily be verified by induction that wbar2 follows the same recursive
     * relation
     *   wbar2[k+1] = s[k] * wbar2[k] - c[k] * v[k+1], *   wbar2[1] = 0, * and we then have
     *   w[k] = c[k] * wbar2[k] + s[k] * v[k+1]
     *          + s[1] * ... * s[k-1] * c[k] * v[1].
     * Introducing w2[k] = w[k] - s[1] * ... * s[k-1] * c[k] * v[1], we find, * from (5.10)
     *   x_l[k] = zeta[1] * w[1] + ... + zeta[k] * w[k]
     *         = zeta[1] * w2[1] + ... + zeta[k] * w2[k]
     *           + (s[1] * c[2] * zeta[2] + ...
     *           + s[1] * ... * s[k-1] * c[k] * zeta[k]) * v[1]
     *         = x_l2[k] + bstep[k] * v[1], * where x_l2[k] is defined by
     *   x_l2[0] = 0, *   x_l2[k+1] = x_l2[k] + zeta[k+1] * w2[k+1], * and bstep is defined by
     *   bstep[1] = 0, *   bstep[k] = bstep[k-1] + s[1] * ... * s[k-1] * c[k] * zeta[k].
     * We also have, from (5.11)
     *   x_c[k] = x_l[k-1] + zbar[k] * wbar[k]
     *         = x_l2[k-1] + zbar[k] * wbar2[k]
     *           + (bstep[k-1] + s[1] * ... * s[k-1] * zbar[k]) * v[1].
     */

    /**
     * <p>
     * A simple container holding the non-const variables used in the
     * iterations. Making the current state of the solver visible from the
     * outside is necessary, because during the iterations, {@code x} does not
     * <em>exactly</em> hold the current estimate of the solution. Indeed, * {@code x} needs in general to be moved from the LQ point to the CG point.
     * Besides, additional upudates must be carried out in case {@code goodb} is
     * set to {@code true}.
     * </p>
     * <p>
     * In all subsequent comments, the description of the state variables refer
     * to their value after a call to {@link #update()}. In these comments, k is
     * the current number of evaluations of matrix-vector products.
     * </p>
     */
    private static class State 
    {
        /** The cubic root of {@link #MACH_PREC}. */
        static const double CBRT_MACH_PREC;

        /** The machine precision. */
        static const double MACH_PREC;

        /** Reference to the linear operator. */
        private const Real_Linear_Operator a;

        /** Reference to the right-hand side vector. */
        private const Real_Vector b;

        /** {@code true} if symmetry of matrix and conditioner must be checked. */
        private const bool check;

        /**
         * The value of the custom tolerance &delta; for the default stopping
         * criterion.
         */
        private const double delta;

        /** The value of beta[k+1]. */
        private double beta;

        /** The value of beta[1]. */
        private double beta1;

        /** The value of bstep[k-1]. */
        private double bstep;

        /** The estimate of the norm of P * rC[k]. */
        private double cgnorm;

        /** The value of dbar[k+1] = -beta[k+1] * c[k-1]. */
        private double dbar;

        /**
         * The value of gamma[k] * zeta[k]. Was called {@code rhs1} in the
         * initial code.
         */
        private double gamma_zeta;

        /** The value of gbar[k]. */
        private double gbar;

        /** The value of max(|alpha[1]|, gamma[1], ..., gamma[k-1]). */
        private double gmax;

        /** The value of min(|alpha[1]|, gamma[1], ..., gamma[k-1]). */
        private double gmin;

        /** Copy of the {@code goodb} parameter. */
        private const bool goodb;

        /** {@code true} if the default convergence criterion is verified. */
        private bool has_converged;

        /** The estimate of the norm of P * rL[k-1]. */
        private double lqnorm;

        /** Reference to the preconditioner, M. */
        private const Real_Linear_Operator m;

        /**
         * The value of (-eps[k+1] * zeta[k-1]). Was called {@code rhs2} in the
         * initial code.
         */
        private double minus_eps_zeta;

        /** The value of M * b. */
        private const Real_Vector mb;

        /** The value of beta[k]. */
        private double oldb;

        /** The value of beta[k] * M^(-1) * P' * v[k]. */
        private Real_Vector r1;

        /** The value of beta[k+1] * M^(-1) * P' * v[k+1]. */
        private Real_Vector r2;

        /**
         * The value of the updated, preconditioned residual P * r. This value is
         * given by {@code min(}{@link #cgnorm}{@code , }{@link #lqnorm}{@code )}.
         */
        private double rnorm;

        /** Copy of the {@code shift} parameter. */
        private const double shift;

        /** The value of s[1] * ... * s[k-1]. */
        private double snprod;

        /**
         * An estimate of the square of the norm of A * V[k], based on Paige and
         * Saunders (1975), equation (3.3).
         */
        private double tnorm;

        /**
         * The value of P' * wbar[k] or P' * (wbar[k] - s[1] * ... * s[k-1] *
         * v[1]) if {@code goodb} is {@code true}. Was called {@code w} in the
         * initial code.
         */
        private Real_Vector wbar;

        /**
         * A reference to the vector to be updated with the solution. Contains
         * the value of x_l[k-1] if {@code goodb} is {@code false}, (x_l[k-1] -
         * bstep[k-1] * v[1]) otherwise.
         */
        private const Real_Vector x_l;

        /** The value of beta[k+1] * P' * v[k+1]. */
        private Real_Vector y;

        /** The value of zeta[1]^2 + ... + zeta[k-1]^2. */
        private double ynorm2;

        /** The value of {@code b == 0} (exact floating-point equality). */
        private bool b_is_null;

        static 
        {
            MACH_PREC = FastMath.ulp(1.);
            CBRT_MACH_PREC = std::cbrt(MACH_PREC);
        }

        /**
         * Creates and inits to k = 1 a instance of this class.
         *
         * @param a the linear operator A of the system
         * @param m the preconditioner, M (can be {@code NULL})
         * @param b the right-hand side vector
         * @param goodb usually {@code false}, except if {@code x} is expected
         * to contain a large multiple of {@code b}
         * @param shift the amount to be subtracted to all diagonal elements of
         * A
         * @param delta the &delta; parameter for the default stopping criterion
         * @param check {@code true} if self-adjointedness of both matrix and
         * preconditioner should be checked
         */
        State(const Real_Linear_Operator a, const Real_Linear_Operator m, const Real_Vector b, const bool goodb, const double shift, const double delta, const bool check) 
        {
            this.a = a;
            this.m = m;
            this.b = b;
            this.x_l = Array_Real_Vector(b.get_dimension());
            this.goodb = goodb;
            this.shift = shift;
            this.mb = m == NULL ? b : m.operate(b);
            this.has_converged = false;
            this.check = check;
            this.delta = delta;
        }

        /**
         * Performs a symmetry check on the specified linear operator, and an
         * exception in case this check fails. Given a linear operator L, and a
         * vector x, this method checks that
         * x' &middot; L &middot; y = y' &middot; L &middot; x
         * (within a given accuracy), where y = L &middot; x.
         * @param x the candidate vector x
         * @param y the candidate vector y = L &middot; x
         * @param z the vector z = L &middot; y
         *
         * @ when the test fails
         */
        private static void check_symmetry(const Real_Vector x, const Real_Vector y, const Real_Vector z)
        {
            const double s = y.dot_product(y);
            const double t = x.dot_product(z);
            const double epsa = (s + MACH_PREC) * CBRT_MACH_PREC;
            if (std::abs(s - t) > epsa) 
            {
                throw std::exception("not implemented");
                //throw (hipparchus::exception::Localized_Core_Formats_Type::NON_SELF_ADJOINT_OPERATOR);
            }
        }

        /**
         * Throws a {@link } with
         * appropriate context.
         * @ in any circumstances
         */
        private static void throw_npdlo_exception()  
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NON_POSITIVE_DEFINITE_OPERATOR);
        }

        /**
         * A clone of the BLAS {@code DAXPY} function, which carries out the
         * operation y &larr; a &middot; x + y. This is for internal use only: no
         * dimension checks are provided.
         *
         * @param a the scalar by which {@code x} is to be multiplied
         * @param x the vector to be added to {@code y}
         * @param y the vector to be incremented
         */
        private static void daxpy(const double& a, const Real_Vector x, const Real_Vector y) 
        {
            const int n = x.get_dimension();
            for (int i{}; i < n; i++) 
            {
                y.set_entry(i, a * x.get_entry(i) + y.get_entry(i));
            }
        }

        /**
         * A BLAS-like function, for the operation z &larr; a &middot; x + b
         * &middot; y + z. This is for internal use only: no dimension checks are
         * provided.
         *
         * @param a the scalar by which {@code x} is to be multiplied
         * @param x the first vector to be added to {@code z}
         * @param b the scalar by which {@code y} is to be multiplied
         * @param y the second vector to be added to {@code z}
         * @param z the vector to be incremented
         */
        private static void daxpbypz(const double& a, const Real_Vector x, const double b, const Real_Vector y, const Real_Vector z) 
        {
            const int n = z.get_dimension();
            for (int i{}; i < n; i++) 
            {
                const double zi;
                zi = a * x.get_entry(i) + b * y.get_entry(i) + z.get_entry(i);
                z.set_entry(i, zi);
            }
        }

        /**
         * <p>
         * Move to the CG point if it seems better. In this version of SYMMLQ, * the convergence tests involve only cgnorm, so we're unlikely to stop
         * at an LQ point, except if the iteration limit interferes.
         * </p>
         * <p>
         * Additional upudates are also carried out in case {@code goodb} is set
         * to {@code true}.
         * </p>
         *
         * @param x the vector to be updated with the refined value of x_l
         */
         void refine_solution(const Real_Vector x) 
         {
            const int n = this.x_l.get_dimension();
            if (lqnorm < cgnorm) 
            {
                if (!goodb) 
                {
                    x.set_sub_vector(0, this.x_l);
                }
else 
                {
                    const double step = bstep / beta1;
                    for (int i{}; i < n; i++) 
                    {
                        const double bi = mb.get_entry(i);
                        const double xi = this.x_l.get_entry(i);
                        x.set_entry(i, xi + step * bi);
                    }
                }
            }
else 
            {
                const double& anorm = std::sqrt(tnorm);
                const double diag = gbar == 0. ? anorm * MACH_PREC : gbar;
                const double zbar = gamma_zeta / diag;
                const double step = (bstep + snprod * zbar) / beta1;
                // ynorm = std::sqrt(ynorm2 + zbar * zbar);
                if (!goodb) 
                {
                    for (int i{}; i < n; i++) 
                    {
                        const double xi = this.x_l.get_entry(i);
                        const double wi = wbar.get_entry(i);
                        x.set_entry(i, xi + zbar * wi);
                    }
                }
else 
                {
                    for (int i{}; i < n; i++) 
                    {
                        const double xi = this.x_l.get_entry(i);
                        const double wi = wbar.get_entry(i);
                        const double bi = mb.get_entry(i);
                        x.set_entry(i, xi + zbar * wi + step * bi);
                    }
                }
            }
        }

        /**
         * Performs the initial phase of the SYMMLQ algorithm. On return, the
         * value of the state variables of {@code this} object correspond to k =
         * 1.
         */
         void init() 
         {
            this.x_l.set(0.);
            /*
             * Set up y for the first Lanczos vector. y and beta1 will be zero
             * if b = 0.
             */
            this.r1 = this.b.copy();
            this.y = this.m == NULL ? this.b.copy() : this.m.operate(this.r1);
            if ((this.m != NULL) && this.check) 
            {
                check_symmetry(this.r1, this.y, this.m.operate(this.y));
            }

            this.beta1 = this.r1.dot_product(this.y);
            if (this.beta1 < 0.) 
            {
                throw_npdlo_exception();
            }
            if (this.beta1 == 0.) 
            {
                /* If b = 0 exactly, stop with x = 0. */
                this.b_is_null = true;
                return;
            }
            this.b_is_null = false;
            this.beta1 = std::sqrt(this.beta1);
            /* At this point
             *   r1 = b, *   y = M * b, *   beta1 = beta[1].
             */
            const Real_Vector v = this.y.map_multiply(1. / this.beta1);
            this.y = this.a.operate(v);
            if (this.check) 
            {
                check_symmetry(v, this.y, this.a.operate(this.y));
            }
            /*
             * Set up y for the second Lanczos vector. y and beta will be zero
             * or very small if b is an eigenvector.
             */
            daxpy(-this.shift, v, this.y);
            const double& alpha = v.dot_product(this.y);
            daxpy(-alpha / this.beta1, this.r1, this.y);
            /*
             * At this point
             *   alpha = alpha[1]
             *   y     = beta[2] * M^(-1) * P' * v[2]
             */
            /* Make sure r2 will be orthogonal to the first v. */
            const double vty = v.dot_product(this.y);
            const double vtv = v.dot_product(v);
            daxpy(-vty / vtv, v, this.y);
            this.r2 = this.y.copy();
            if (this.m != NULL) 
            {
                this.y = this.m.operate(this.r2);
            }
            this.oldb = this.beta1;
            this.beta = this.r2.dot_product(this.y);
            if (this.beta < 0.) 
            {
                throw_npdlo_exception();
            }
            this.beta = std::sqrt(this.beta);
            /*
             * At this point
             *   oldb = beta[1]
             *   beta = beta[2]
             *   y  = beta[2] * P' * v[2]
             *   r2 = beta[2] * M^(-1) * P' * v[2]
             */
            this.cgnorm = this.beta1;
            this.gbar = alpha;
            this.dbar = this.beta;
            this.gamma_zeta = this.beta1;
            this.minus_eps_zeta = 0.;
            this.bstep = 0.;
            this.snprod = 1.;
            this.tnorm = alpha * alpha + this.beta * this.beta;
            this.ynorm2 = 0.;
            this.gmax = std::abs(alpha) + MACH_PREC;
            this.gmin = this.gmax;

            if (this.goodb) 
            {
                this.wbar = Array_Real_Vector(this.a.get_row_dimension());
                this.wbar.set(0.);
            }
else 
            {
                this.wbar = v;
            }
            update_norms();
        }

        /**
         * Performs the next iteration of the algorithm. The iteration count
         * should be incremented prior to calling this method. On return, the
         * value of the state variables of {@code this} object correspond to the
         * current iteration count {@code k}.
         */
        void update() 
        {
            const Real_Vector v = y.map_multiply(1. / beta);
            y = a.operate(v);
            daxpbypz(-shift, v, -beta / oldb, r1, y);
            const double& alpha = v.dot_product(y);
            /*
             * At this point
             *   v     = P' * v[k], *   y     = (A - shift * I) * P' * v[k] - beta[k] * M^(-1) * P' * v[k-1], *   alpha = v'[k] * P * (A - shift * I) * P' * v[k]
             *           - beta[k] * v[k]' * P * M^(-1) * P' * v[k-1]
             *         = v'[k] * P * (A - shift * I) * P' * v[k]
             *           - beta[k] * v[k]' * v[k-1]
             *         = alpha[k].
             */
            daxpy(-alpha / beta, r2, y);
            /*
             * At this point
             *   y = (A - shift * I) * P' * v[k] - alpha[k] * M^(-1) * P' * v[k]
             *       - beta[k] * M^(-1) * P' * v[k-1]
             *     = M^(-1) * P' * (P * (A - shift * I) * P' * v[k] -alpha[k] * v[k]
             *       - beta[k] * v[k-1])
             *     = beta[k+1] * M^(-1) * P' * v[k+1], * from Paige and Saunders (1975), equation (3.2).
             *
             * WATCH-IT: the two following lines work only because y is no longer
             * updated up to the end of the present iteration, and is
             * reinitialized at the beginning of the next iteration.
             */
            r1 = r2;
            r2 = y;
            if (m != NULL) 
            {
                y = m.operate(r2);
            }
            oldb = beta;
            beta = r2.dot_product(y);
            if (beta < 0.) 
            {
                throw_npdlo_exception();
            }
            beta = std::sqrt(beta);
            /*
             * At this point
             *   r1 = beta[k] * M^(-1) * P' * v[k], *   r2 = beta[k+1] * M^(-1) * P' * v[k+1], *   y  = beta[k+1] * P' * v[k+1], *   oldb = beta[k], *   beta = beta[k+1].
             */
            tnorm += alpha * alpha + oldb * oldb + beta * beta;
            /*
             * Compute the next plane rotation for Q. See Paige and Saunders
             * (1975), equation (5.6), with
             *   gamma = gamma[k-1], *   c     = c[k-1], *   s     = s[k-1].
             */
            const double gamma = std::sqrt(gbar * gbar + oldb * oldb);
            const double c = gbar / gamma;
            const double s = oldb / gamma;
            /*
             * The relations
             *   gbar[k] = s[k-1] * (-c[k-2] * beta[k]) - c[k-1] * alpha[k]
             *           = s[k-1] * dbar[k] - c[k-1] * alpha[k], *   delta[k] = c[k-1] * dbar[k] + s[k-1] * alpha[k], * are not stated in Paige and Saunders (1975), but can be retrieved
             * by expanding the (k, k-1) and (k, k) coefficients of the matrix in
             * equation (5.5).
             */
            const double deltak = c * dbar + s * alpha;
            gbar = s * dbar - c * alpha;
            const double eps = s * beta;
            dbar = -c * beta;
            const double zeta = gamma_zeta / gamma;
            /*
             * At this point
             *   gbar   = gbar[k]
             *   deltak = delta[k]
             *   eps    = eps[k+1]
             *   dbar   = dbar[k+1]
             *   zeta   = zeta[k-1]
             */
            const double zeta_c = zeta * c;
            const double zeta_s = zeta * s;
            const int n = x_l.get_dimension();
            for (int i{}; i < n; i++) 
            {
                const double xi = x_l.get_entry(i);
                const double vi = v.get_entry(i);
                const double wi = wbar.get_entry(i);
                x_l.set_entry(i, xi + wi * zeta_c + vi * zeta_s);
                wbar.set_entry(i, wi * s - vi * c);
            }
            /*
             * At this point
             *   x = x_l[k-1], *   ptwbar = P' wbar[k], * see Paige and Saunders (1975), equations (5.9) and (5.10).
             */
            bstep += snprod * c * zeta;
            snprod *= s;
            gmax = std::max(gmax, gamma);
            gmin = std::min(gmin, gamma);
            ynorm2 += zeta * zeta;
            gamma_zeta = minus_eps_zeta - deltak * zeta;
            minus_eps_zeta = -eps * zeta;
            /*
             * At this point
             *   snprod       = s[1] * ... * s[k-1], *   gmax         = max(|alpha[1]|, gamma[1], ..., gamma[k-1]), *   gmin         = min(|alpha[1]|, gamma[1], ..., gamma[k-1]), *   ynorm2       = zeta[1]^2 + ... + zeta[k-1]^2, *   gamma_zeta    = gamma[k] * zeta[k], *   minus_eps_zeta = -eps[k+1] * zeta[k-1].
             * The relation for gamma_zeta can be retrieved from Paige and
             * Saunders (1975), equation (5.4a), last line of the vector
             * gbar[k] * zbar[k] = -eps[k] * zeta[k-2] - delta[k] * zeta[k-1].
             */
            update_norms();
        }

        /**
         * Computes the norms of the residuals, and checks for convergence.
         * Updates {@link #lqnorm} and {@link #cgnorm}.
         */
        private void update_norms() 
        {
            const double& anorm = std::sqrt(tnorm);
            const double ynorm = std::sqrt(ynorm2);
            const double epsa = anorm * MACH_PREC;
            const double epsx = anorm * ynorm * MACH_PREC;
            const double epsr = anorm * ynorm * delta;
            const double diag = gbar == 0. ? epsa : gbar;
            lqnorm = std::sqrt(gamma_zeta * gamma_zeta +
                                   minus_eps_zeta * minus_eps_zeta);
            const double qrnorm = snprod * beta1;
            cgnorm = qrnorm * beta / std::abs(diag);

            /*
             * Estimate cond(A). In this version we look at the diagonals of L
             * in the factorization of the tridiagonal matrix, T = L * Q.
             * Sometimes, T[k] can be misleadingly ill-conditioned when T[k+1]
             * is not, so we must be careful not to overestimate acond.
             */
            const double& acond;
            if (lqnorm <= cgnorm) 
            {
                acond = gmax / gmin;
            }
            else 
            {
                acond = gmax / std::min(gmin, std::abs(diag));
            }
            if (acond * MACH_PREC >= 0.1) 
            {
                throw std::exception("not implemented");
                //throw (hipparchus::exception::Localized_Core_Formats_Type::ILL_CONDITIONED_OPERATOR, acond);
            }
            if (beta1 <= epsx) 
            {
                /*
                 * x has converged to an eigenvector of A corresponding to the
                 * eigenvalue shift.
                 */
                throw std::exception("not implemented");
                //throw (hipparchus::exception::Localized_Core_Formats_Type::SINGULAR_OPERATOR);
            }
            rnorm = std::min(cgnorm, lqnorm);
            has_converged = (cgnorm <= epsx) || (cgnorm <= epsr);
        }

        /**
         * Returns {@code true} if the default stopping criterion is fulfilled.
         *
         * @return {@code true} if convergence of the iterations has occurred
         */
        bool has_converged() 
        {
            return has_converged;
        }

        /**
         * Returns {@code true} if the right-hand side vector is zero exactly.
         *
         * @return the bool value of {@code b == 0}
         */
        bool b_equals_null_vector() 
        {
            return b_is_null;
        }

        /**
         * Returns {@code true} if {@code beta} is essentially zero. This method
         * is used to check for early stop of the iterations.
         *
         * @return {@code true} if {@code beta < }{@link #MACH_PREC}
         */
        bool beta_equals_zero() 
        {
            return beta < MACH_PREC;
        }

        /**
         * Returns the norm of the updated, preconditioned residual.
         *
         * @return the norm of the residual, ||P * r||
         */
        double get_norm_of_residual() 
        {
            return rnorm;
        }
    }

    /**
     * Creates a instance of this class, with <a href="#stopcrit">default
     * stopping criterion</a>. Note that setting {@code check} to {@code true}
     * entails an extra matrix-vector product in the initial phase.
     *
     * @param max_iterations the maximum number of iterations
     * @param delta the &delta; parameter for the default stopping criterion
     * @param check {@code true} if self-adjointedness of both matrix and
     * preconditioner should be checked
     */
    public Symm_LQ(const int max_iterations, const double delta, const bool check) 
    {
        super(max_iterations);
        this.delta = delta;
        this.check = check;
    }

    /**
     * Creates a instance of this class, with <a href="#stopcrit">default
     * stopping criterion</a> and custom iteration manager. Note that setting
     * {@code check} to {@code true} entails an extra matrix-vector product in
     * the initial phase.
     *
     * @param manager the custom iteration manager
     * @param delta the &delta; parameter for the default stopping criterion
     * @param check {@code true} if self-adjointedness of both matrix and
     * preconditioner should be checked
     */
    public Symm_LQ(const Iteration_Manager manager, const double delta, const bool check) 
    {
        super(manager);
        this.delta = delta;
        this.check = check;
    }

    /**
     * Returns {@code true} if symmetry of the matrix, and symmetry as well as
     * positive definiteness of the preconditioner should be checked.
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
     * @ if {@link #should_check()} is
     * {@code true}, and {@code a} or {@code m} is not self-adjoint
     * @ if {@code m} is not
     * positive definite
     * @ if {@code a} is ill-conditioned
     */
    //override
    public Real_Vector solve(const Real_Linear_Operator a, const Real_Linear_Operator m, const Real_Vector b) throws
        ,  
        {
        //Math_Utils::check_not_null(a);
        const Real_Vector x = Array_Real_Vector(a.get_column_dimension());
        return solve_in_place(a, m, b, x, false, 0.);
    }

    /**
     * Returns an estimate of the solution to the linear system (A - shift
     * &middot; I) &middot; x = b.
     * <p>
     * If the solution x is expected to contain a large multiple of {@code b}
     * (as in Rayleigh-quotient iteration), then better precision may be
     * achieved with {@code goodb} set to {@code true}; this however requires an
     * extra call to the preconditioner.
     * </p>
     * <p>
     * {@code shift} should be zero if the system A &middot; x = b is to be
     * solved. Otherwise, it could be an approximation to an eigenvalue of A, * such as the Rayleigh quotient b<sup>T</sup> &middot; A &middot; b /
     * (b<sup>T</sup> &middot; b) corresponding to the vector b. If b is
     * sufficiently like an eigenvector corresponding to an eigenvalue near
     * shift, then the computed x may have very large components. When
     * normalized, x may be closer to an eigenvector than b.
     * </p>
     *
     * @param a the linear operator A of the system
     * @param m the preconditioner, M (can be {@code NULL})
     * @param b the right-hand side vector
     * @param goodb usually {@code false}, except if {@code x} is expected to
     * contain a large multiple of {@code b}
     * @param shift the amount to be subtracted to all diagonal elements of A
     * @return a reference to {@code x} (shallow copy)
     * @ if one of the parameters is {@code NULL}
     * @ if {@code a} or {@code m} is not square
     * @ if {@code m} or {@code b} have dimensions
     * inconsistent with {@code a}
     * @Math_Illegal_State_Exception at exhaustion of the iteration count, * unless a custom
     * {@link org.hipparchus.util.Incrementor.Max_Count_Exceeded_Callback callback}
     * has been set at construction of the {@link Iteration_Manager}
     * @ if {@link #should_check()} is
     * {@code true}, and {@code a} or {@code m} is not self-adjoint
     * @ if {@code m} is not
     * positive definite
     * @ if {@code a} is ill-conditioned
     */
    public Real_Vector solve(const Real_Linear_Operator a, const Real_Linear_Operator m, const Real_Vector b, const bool goodb, const double shift)
         
        {
        //Math_Utils::check_not_null(a);
        const Real_Vector x = Array_Real_Vector(a.get_column_dimension());
        return solve_in_place(a, m, b, x, goodb, shift);
    }

    /**
     * {@inherit_doc}
     *
     * @param x not meaningful in this implementation; should not be considered
     * as an initial guess (<a href="#initguess">more</a>)
     * @ if {@link #should_check()} is
     * {@code true}, and {@code a} or {@code m} is not self-adjoint
     * @ if {@code m} is not positive
     * definite
     * @ if {@code a} is ill-conditioned
     */
    //override
    public Real_Vector solve(const Real_Linear_Operator a, const Real_Linear_Operator m, const Real_Vector b, const Real_Vector x)
        ,  
        {
        //Math_Utils::check_not_null(x);
        return solve_in_place(a, m, b, x.copy(), false, 0.);
    }

    /**
     * {@inherit_doc}
     *
     * @ if {@link #should_check()} is
     * {@code true}, and {@code a} is not self-adjoint
     * @ if {@code a} is ill-conditioned
     */
    //override
    public Real_Vector solve(const Real_Linear_Operator a, const Real_Vector b)
        ,  
        {
        //Math_Utils::check_not_null(a);
        const Real_Vector x = Array_Real_Vector(a.get_column_dimension());
        x.set(0.);
        return solve_in_place(a, NULL, b, x, false, 0.);
    }

    /**
     * Returns the solution to the system (A - shift &middot; I) &middot; x = b.
     * <p>
     * If the solution x is expected to contain a large multiple of {@code b}
     * (as in Rayleigh-quotient iteration), then better precision may be
     * achieved with {@code goodb} set to {@code true}.
     * </p>
     * <p>
     * {@code shift} should be zero if the system A &middot; x = b is to be
     * solved. Otherwise, it could be an approximation to an eigenvalue of A, * such as the Rayleigh quotient b<sup>T</sup> &middot; A &middot; b /
     * (b<sup>T</sup> &middot; b) corresponding to the vector b. If b is
     * sufficiently like an eigenvector corresponding to an eigenvalue near
     * shift, then the computed x may have very large components. When
     * normalized, x may be closer to an eigenvector than b.
     * </p>
     *
     * @param a the linear operator A of the system
     * @param b the right-hand side vector
     * @param goodb usually {@code false}, except if {@code x} is expected to
     * contain a large multiple of {@code b}
     * @param shift the amount to be subtracted to all diagonal elements of A
     * @return a reference to {@code x}
     * @ if one of the parameters is {@code NULL}
     * @ if {@code a} is not square
     * @ if {@code b} has dimensions
     * inconsistent with {@code a}
     * @Math_Illegal_State_Exception at exhaustion of the iteration count, * unless a custom
     * {@link org.hipparchus.util.Incrementor.Max_Count_Exceeded_Callback callback}
     * has been set at construction of the {@link Iteration_Manager}
     * @ if {@link #should_check()} is
     * {@code true}, and {@code a} is not self-adjoint
     * @ if {@code a} is ill-conditioned
     */
    public Real_Vector solve(const Real_Linear_Operator a, const Real_Vector b, const bool goodb, const double shift)
         
        {
        //Math_Utils::check_not_null(a);
        const Real_Vector x = Array_Real_Vector(a.get_column_dimension());
        return solve_in_place(a, NULL, b, x, goodb, shift);
    }

    /**
     * {@inherit_doc}
     *
     * @param x not meaningful in this implementation; should not be considered
     * as an initial guess (<a href="#initguess">more</a>)
     * @ if {@link #should_check()} is
     * {@code true}, and {@code a} is not self-adjoint
     * @ if {@code a} is ill-conditioned
     */
    //override
    public Real_Vector solve(const Real_Linear_Operator a, const Real_Vector b, const Real_Vector x) ,  
    {
        //Math_Utils::check_not_null(x);
        return solve_in_place(a, NULL, b, x.copy(), false, 0.);
    }

    /**
     * {@inherit_doc}
     *
     * @param x the vector to be updated with the solution; {@code x} should
     * not be considered as an initial guess (<a href="#initguess">more</a>)
     * @ if {@link #should_check()} is
     * {@code true}, and {@code a} or {@code m} is not self-adjoint
     * @ if {@code m} is not
     * positive definite
     * @ if {@code a} is ill-conditioned
     */
    //override
    public Real_Vector solve_in_place(const Real_Linear_Operator a, const Real_Linear_Operator m, const Real_Vector b, const Real_Vector x)
        ,  
        {
        return solve_in_place(a, m, b, x, false, 0.);
    }

    /**
     * Returns an estimate of the solution to the linear system (A - shift
     * &middot; I) &middot; x = b. The solution is computed in-place.
     * <p>
     * If the solution x is expected to contain a large multiple of {@code b}
     * (as in Rayleigh-quotient iteration), then better precision may be
     * achieved with {@code goodb} set to {@code true}; this however requires an
     * extra call to the preconditioner.
     * </p>
     * <p>
     * {@code shift} should be zero if the system A &middot; x = b is to be
     * solved. Otherwise, it could be an approximation to an eigenvalue of A, * such as the Rayleigh quotient b<sup>T</sup> &middot; A &middot; b /
     * (b<sup>T</sup> &middot; b) corresponding to the vector b. If b is
     * sufficiently like an eigenvector corresponding to an eigenvalue near
     * shift, then the computed x may have very large components. When
     * normalized, x may be closer to an eigenvector than b.
     * </p>
     *
     * @param a the linear operator A of the system
     * @param m the preconditioner, M (can be {@code NULL})
     * @param b the right-hand side vector
     * @param x the vector to be updated with the solution; {@code x} should
     * not be considered as an initial guess (<a href="#initguess">more</a>)
     * @param goodb usually {@code false}, except if {@code x} is expected to
     * contain a large multiple of {@code b}
     * @param shift the amount to be subtracted to all diagonal elements of A
     * @return a reference to {@code x} (shallow copy).
     * @ if one of the parameters is {@code NULL}
     * @ if {@code a} or {@code m} is not square
     * @ if {@code m}, {@code b} or {@code x}
     * have dimensions inconsistent with {@code a}.
     * @Math_Illegal_State_Exception at exhaustion of the iteration count, * unless a custom
     * {@link org.hipparchus.util.Incrementor.Max_Count_Exceeded_Callback callback}
     * has been set at construction of the {@link Iteration_Manager}
     * @ if {@link #should_check()} is
     * {@code true}, and {@code a} or {@code m} is not self-adjoint
     * @ if {@code m} is not positive definite
     * @ if {@code a} is ill-conditioned
     */
    public Real_Vector solve_in_place(const Real_Linear_Operator a, const Real_Linear_Operator m, const Real_Vector b, const Real_Vector x, const bool goodb, const double shift)
         
        {
        check_parameters(a, m, b, x);

        const Iteration_Manager manager = get_iteration_manager();
        /* Initialization counts as an iteration. */
        manager.reset_iteration_count();
        manager.increment_iteration_count();

        const State state;
        state = State(a, m, b, goodb, shift, delta, check);
        state.init();
        state.refine_solution(x);
        Iterative_Linear_SolverEvent event;
        event = Default_iterativeLinearSolverEvent(this, manager.get_iterations(), x, b, state.get_norm_of_residual());
        if (state.b_equals_null_vector()) 
        {
            /* If b = 0 exactly, stop with x = 0. */
            manager.fire_termination_event(event);
            return x;
        }
        /* Cause termination if beta is essentially zero. */
        const bool early_stop;
        early_stop = state.beta_equals_zero() || state.has_converged();
        manager.fire_initialization_event(event);
        if (!early_stop) 
        {
            do 
            {
                manager.increment_iteration_count();
                event = Default_iterativeLinearSolverEvent(this, manager.get_iterations(), x, b, state.get_norm_of_residual());
                manager.fire_iteration_started_event(event);
                state.update();
                state.refine_solution(x);
                event = Default_iterativeLinearSolverEvent(this, manager.get_iterations(), x, b, state.get_norm_of_residual());
                manager.fire_iteration_performed_event(event);
            } while (!state.has_converged());
        }
        event = Default_iterativeLinearSolverEvent(this, manager.get_iterations(), x, b, state.get_norm_of_residual());
        manager.fire_termination_event(event);
        return x;
    }

    /**
     * {@inherit_doc}
     *
     * @param x the vector to be updated with the solution; {@code x} should
     * not be considered as an initial guess (<a href="#initguess">more</a>)
     * @ if {@link #should_check()} is
     * {@code true}, and {@code a} is not self-adjoint
     * @ if {@code a} is ill-conditioned
     */
    //override
    public Real_Vector solve_in_place(const Real_Linear_Operator a, const Real_Vector b, const Real_Vector x) ,  
    {
        return solve_in_place(a, NULL, b, x, false, 0.);
    }
}


