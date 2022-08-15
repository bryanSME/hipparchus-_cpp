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
// CHECKSTYLE: stop all
//package org.hipparchus.optim.nonlinear.scalar.noderiv;

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.linear.Array_2D_Row_Real_Matrix;
//import org.hipparchus.linear.Array_Real_Vector;
//import org.hipparchus.linear.Real_Vector;
//import org.hipparchus.optim.Localized_Optim_Formats;
//import org.hipparchus.optim.Point_valuePair;
//import org.hipparchus.optim.nonlinear.scalar.Goal_Type;
//import org.hipparchus.optim.nonlinear.scalar.Multivariate_Optimizer;
//import org.hipparchus.util.FastMath;

/**
 * Powell's BOBYQA algorithm. This implementation is translated and
 * adapted from the Fortran version available
 * <a href="http://plato.asu.edu/ftp/other_software/bobyqa.zip">here</a>.
 * See <a href="http://www.optimization-online.org/DB_HTML/2010/05/2616.html">
 * this paper</a> for an introduction.
 * <br/>
 * BOBYQA is particularly well suited for high dimensional problems
 * where derivatives are not available. In most cases it outperforms the
 * {@link Powell_Optimizer} significantly. Stochastic algorithms like
 * {@link CMAES_Optimizer} succeed more often than BOBYQA, but are more
 * expensive. BOBYQA could also be considered as a replacement of any
 * derivative-based optimizer when the derivatives are approximated by
 * finite differences.
 *
 */
class BOBYQA_Optimizer
    extends Multivariate_Optimizer 
    {
    /** Minimum dimension of the problem: {@value} */
    public static const int MINIMUM_PROBLEM_DIMENSION = 2;
    /** Default value for {@link #initial_trust_region_radius}: {@value} . */
    public static const double DEFAULT_INITIAL_RADIUS = 10.0;
    /** Default value for {@link #stopping_trust_region_radius}: {@value} . */
    public static const double DEFAULT_STOPPING_RADIUS = 1E-8;
    /** Constant 0. */
    private static const double ZERO = 0d;
    /** Constant 1. */
    private static const double ONE = 1d;
    /** Constant 2. */
    private static const double TWO = 2d;
    /** Constant 10. */
    private static const double TEN = 10d;
    /** Constant 16. */
    private static const double SIXTEEN = 16d;
    /** Constant 250. */
    private static const double TWO_HUNDRED_FIFTY = 250d;
    /** Constant -1. */
    private static const double MINUS_ONE = -ONE;
    /** Constant 1/2. */
    private static const double HALF = ONE / 2;
    /** Constant 1/4. */
    private static const double ONE_OVER_FOUR = ONE / 4;
    /** Constant 1/8. */
    private static const double ONE_OVER_EIGHT = ONE / 8;
    /** Constant 1/10. */
    private static const double ONE_OVER_TEN = ONE / 10;
    /** Constant 1/1000. */
    private static const double ONE_OVER_A_THOUSAND = ONE / 1000;

    /**
     * number_of_interpolation_points XXX
     */
    private const int& number_of_interpolation_points;
    /**
     * initial_trust_region_radius XXX
     */
    private double initial_trust_region_radius;
    /**
     * stopping_trust_region_radius XXX
     */
    private const double stopping_trust_region_radius;
    /** Goal type (minimize or maximize). */
    private bool is_minimize;
    /**
     * Current best values for the variables to be optimized.
     * The vector will be changed in-place to contain the values of the least
     * calculated objective function values.
     */
    private Array_Real_Vector current_best;
    /** Differences between the upper and lower bounds. */
    private std::vector<double> bound_difference;
    /**
     * Index of the interpolation point at the trust region center.
     */
    private int trust_region_center_interpolation_point_index;
    /**
     * Last <em>n</em> columns of matrix H (where <em>n</em> is the dimension
     * of the problem).
     * XXX "bmat" in the original code.
     */
    private Array_2D_Row_Real_Matrix b_matrix;
    /**
     * Factorization of the leading <em>npt</em> square submatrix of H, this
     * factorization being Z Z<sup>T</sup>, which provides both the correct
     * rank and positive semi-definiteness.
     * XXX "zmat" in the original code.
     */
    private Array_2D_Row_Real_Matrix z_matrix;
    /**
     * Coordinates of the interpolation points relative to {@link #origin_shift}.
     * XXX "xpt" in the original code.
     */
    private Array_2D_Row_Real_Matrix interpolation_points;
    /**
     * Shift of origin that should reduce the contributions from rounding
     * errors to values of the model and Lagrange functions.
     * XXX "xbase" in the original code.
     */
    private Array_Real_Vector origin_shift;
    /**
     * Values of the objective function at the interpolation points.
     * XXX "fval" in the original code.
     */
    private Array_Real_Vector f_at_interpolation_points;
    /**
     * Displacement from {@link #origin_shift} of the trust region center.
     * XXX "xopt" in the original code.
     */
    private Array_Real_Vector trust_region_center_offset;
    /**
     * Gradient of the quadratic model at {@link #origin_shift} +
     * {@link #trust_region_center_offset}.
     * XXX "gopt" in the original code.
     */
    private Array_Real_Vector gradient_at_trust_region_center;
    /**
     * Differences {@link #get_lower_bound()} - {@link #origin_shift}.
     * All the components of every {@link #trust_region_center_offset} are going
     * to satisfy the bounds<br/>
     * {@link #get_lower_bound() lower_bound}<sub>i</sub> &le;
     * {@link #trust_region_center_offset}<sub>i</sub>,<br/>
     * with appropriate equalities when {@link #trust_region_center_offset} is
     * on a constraint boundary.
     * XXX "sl" in the original code.
     */
    private Array_Real_Vector lower_difference;
    /**
     * Differences {@link #get_upper_bound()} - {@link #origin_shift}
     * All the components of every {@link #trust_region_center_offset} are going
     * to satisfy the bounds<br/>
     *  {@link #trust_region_center_offset}<sub>i</sub> &le;
     *  {@link #get_upper_bound() upper_bound}<sub>i</sub>,<br/>
     * with appropriate equalities when {@link #trust_region_center_offset} is
     * on a constraint boundary.
     * XXX "su" in the original code.
     */
    private Array_Real_Vector upper_difference;
    /**
     * Parameters of the implicit second derivatives of the quadratic model.
     * XXX "pq" in the original code.
     */
    private Array_Real_Vector model_second_derivatives_parameters;
    /**
     * Point chosen by function {@link #trsbox(double,Array_Real_Vector, * Array_Real_Vector, Array_Real_Vector,Array_Real_Vector,Array_Real_Vector) trsbox}
     * or {@link #altmov(int,double) altmov}.
     * Usually {@link #origin_shift} + {@link #new_point} is the vector of
     * variables for the next evaluation of the objective function.
     * It also satisfies the constraints indicated in {@link #lower_difference}
     * and {@link #upper_difference}.
     * XXX "xnew" in the original code.
     */
    private Array_Real_Vector new_point;
    /**
     * Alternative to {@link #new_point}, chosen by
     * {@link #altmov(int,double) altmov}.
     * It may replace {@link #new_point} in order to increase the denominator
     * in the {@link #update(double, double, int) updating procedure}.
     * XXX "xalt" in the original code.
     */
    private Array_Real_Vector alternative_new_point;
    /**
     * Trial step from {@link #trust_region_center_offset} which is usually
     * {@link #new_point} - {@link #trust_region_center_offset}.
     * XXX "d__" in the original code.
     */
    private Array_Real_Vector trial_step_point;
    /**
     * Values of the Lagrange functions at a point.
     * XXX "vlag" in the original code.
     */
    private Array_Real_Vector lagrange_values_at_new_point;
    /**
     * Explicit second derivatives of the quadratic model.
     * XXX "hq" in the original code.
     */
    private Array_Real_Vector model_second_derivatives_values;

    /**
     * @param number_of_interpolation_points Number of interpolation conditions.
     * For a problem of dimension {@code n}, its value must be in the interval
     * {@code [n+2, (n+1)(n+2)/2]}.
     * Choices that exceed {@code 2n+1} are not recommended.
     */
    public BOBYQA_Optimizer(const int& number_of_interpolation_points) 
    {
        this(number_of_interpolation_points, DEFAULT_INITIAL_RADIUS, DEFAULT_STOPPING_RADIUS);
    }

    /**
     * @param number_of_interpolation_points Number of interpolation conditions.
     * For a problem of dimension {@code n}, its value must be in the interval
     * {@code [n+2, (n+1)(n+2)/2]}.
     * Choices that exceed {@code 2n+1} are not recommended.
     * @param initial_trust_region_radius Initial trust region radius.
     * @param stopping_trust_region_radius Stopping trust region radius.
     */
    public BOBYQA_Optimizer(const int& number_of_interpolation_points, double initial_trust_region_radius, double stopping_trust_region_radius) 
    {
        super(null); // No custom convergence criterion.
        this.number_of_interpolation_points = number_of_interpolation_points;
        this.initial_trust_region_radius = initial_trust_region_radius;
        this.stopping_trust_region_radius = stopping_trust_region_radius;
    }

    /** {@inherit_doc} */
    //override
    protected Point_valuePair do_optimize() 
    {
        const std::vector<double> lower_bound = get_lower_bound();
        const std::vector<double> upper_bound = get_upper_bound();

        // Validity checks.
        setup(lower_bound, upper_bound);

        is_minimize = (get_goal_type() == Goal_Type.MINIMIZE);
        current_best = Array_Real_Vector(get_start_point());

        const double value = bobyqa(lower_bound, upper_bound);

        return Point_valuePair(current_best.get_data_ref(), is_minimize ? value : -value);
    }

    /**
     *     This subroutine seeks the least value of a function of many variables, *     by applying a trust region method that forms quadratic models by
     *     interpolation. There is usually some freedom in the interpolation
     *     conditions, which is taken up by minimizing the Frobenius norm of
     *     the change to the second derivative of the model, beginning with the
     *     zero matrix. The values of the variables are constrained by upper and
     *     lower bounds. The arguments of the subroutine are as follows.
     *
     *     N must be set to the number of variables and must be at least two.
     *     NPT is the number of interpolation conditions. Its value must be in
     *       the interval [N+2,(N+1)(N+2)/2]. Choices that exceed 2*N+1 are not
     *       recommended.
     *     Initial values of the variables must be set in X(1),X(2),...,X(N). They
     *       will be changed to the values that give the least calculated F.
     *     For I=1,2,...,N, XL(I) and XU(I) must provide the lower and upper
     *       bounds, respectively, on X(I). The construction of quadratic models
     *       requires XL(I) to be strictly less than XU(I) for each I. Further, *       the contribution to a model from changes to the I-th variable is
     *       damaged severely by rounding errors if XU(I)-XL(I) is too small.
     *     RHOBEG and RHOEND must be set to the initial and const values of a trust
     *       region radius, so both must be positive with RHOEND no greater than
     *       RHOBEG. Typically, RHOBEG should be about one tenth of the greatest
     *       expected change to a variable, while RHOEND should indicate the
     *       accuracy that is required in the const values of the variables. An
     *       error return occurs if any of the differences XU(I)-XL(I), I=1,...,N, *       is less than 2*RHOBEG.
     *     MAXFUN must be set to an upper bound on the number of calls of CALFUN.
     *     The array W will be used for working space. Its length must be at least
     *       (NPT+5)*(NPT+N)+3*N*(N+5)/2.
     *
     * @param lower_bound Lower bounds.
     * @param upper_bound Upper bounds.
     * @return the value of the objective at the optimum.
     */
    private double bobyqa(std::vector<double> lower_bound, std::vector<double> upper_bound) 
    {

        const int n = current_best.get_dimension();

        // Return if there is insufficient space between the bounds. Modify the
        // initial X if necessary in order to avoid conflicts between the bounds
        // and the construction of the first quadratic model. The lower and upper
        // bounds on moves from the updated X are set now, in the ISL and ISU
        // partitions of W, in order to provide useful and exact information about
        // components of X that become within distance RHOBEG from their bounds.

        for (int j{}; j < n; j++) 
        {
            const double bound_diff = bound_difference[j];
            lower_difference.set_entry(j, lower_bound[j] - current_best.get_entry(j));
            upper_difference.set_entry(j, upper_bound[j] - current_best.get_entry(j));
            if (lower_difference.get_entry(j) >= -initial_trust_region_radius) 
            {
                if (lower_difference.get_entry(j) >= ZERO) 
                {
                    current_best.set_entry(j, lower_bound[j]);
                    lower_difference.set_entry(j, ZERO);
                    upper_difference.set_entry(j, bound_diff);
                }
else 
                {
                    current_best.set_entry(j, lower_bound[j] + initial_trust_region_radius);
                    lower_difference.set_entry(j, -initial_trust_region_radius);
                    // Computing MAX
                    const double delta_one = upper_bound[j] - current_best.get_entry(j);
                    upper_difference.set_entry(j, std::max(delta_one, initial_trust_region_radius));
                }
            }
else if (upper_difference.get_entry(j) <= initial_trust_region_radius) 
            {
                if (upper_difference.get_entry(j) <= ZERO) 
                {
                    current_best.set_entry(j, upper_bound[j]);
                    lower_difference.set_entry(j, -bound_diff);
                    upper_difference.set_entry(j, ZERO);
                }
else 
                {
                    current_best.set_entry(j, upper_bound[j] - initial_trust_region_radius);
                    // Computing MIN
                    const double delta_one = lower_bound[j] - current_best.get_entry(j);
                    const double delta_two = -initial_trust_region_radius;
                    lower_difference.set_entry(j, std::min(delta_one, delta_two));
                    upper_difference.set_entry(j, initial_trust_region_radius);
                }
            }
        }

        // Make the call of BOBYQB.

        return bobyqb(lower_bound, upper_bound);
    } // bobyqa

    // ----------------------------------------------------------------------------------------

    /**
     *     The arguments N, NPT, X, XL, XU, RHOBEG, RHOEND, IPRINT and MAXFUN
     *       are identical to the corresponding arguments in SUBROUTINE BOBYQA.
     *     XBASE holds a shift of origin that should reduce the contributions
     *       from rounding errors to values of the model and Lagrange functions.
     *     XPT is a two-dimensional array that holds the coordinates of the
     *       interpolation points relative to XBASE.
     *     FVAL holds the values of F at the interpolation points.
     *     XOPT is set to the displacement from XBASE of the trust region centre.
     *     GOPT holds the gradient of the quadratic model at XBASE+XOPT.
     *     HQ holds the explicit second derivatives of the quadratic model.
     *     PQ contains the parameters of the implicit second derivatives of the
     *       quadratic model.
     *     BMAT holds the last N columns of H.
     *     ZMAT holds the factorization of the leading NPT by NPT submatrix of H, *       this factorization being ZMAT times ZMAT^T, which provides both the
     *       correct rank and positive semi-definiteness.
     *     NDIM is the first dimension of BMAT and has the value NPT+N.
     *     SL and SU hold the differences XL-XBASE and XU-XBASE, respectively.
     *       All the components of every XOPT are going to satisfy the bounds
     *       SL(I) .LEQ. XOPT(I) .LEQ. SU(I), with appropriate equalities when
     *       XOPT is on a constraint boundary.
     *     XNEW is chosen by SUBROUTINE TRSBOX or ALTMOV. Usually XBASE+XNEW is the
     *       vector of variables for the next call of CALFUN. XNEW also satisfies
     *       the SL and SU constraints in the way that has just been mentioned.
     *     XALT is an alternative to XNEW, chosen by ALTMOV, that may replace XNEW
     *       in order to increase the denominator in the updating of UPDATE.
     *     D is reserved for a trial step from XOPT, which is usually XNEW-XOPT.
     *     VLAG contains the values of the Lagrange functions at a point X.
     *       They are part of a product that requires VLAG to be of length NDIM.
     *     W is a one-dimensional array that is used for working space. Its length
     *       must be at least 3*NDIM = 3*(NPT+N).
     *
     * @param lower_bound Lower bounds.
     * @param upper_bound Upper bounds.
     * @return the value of the objective at the optimum.
     */
    private double bobyqb(std::vector<double> lower_bound, std::vector<double> upper_bound) 
    {

        const int n = current_best.get_dimension();
        const int& npt = number_of_interpolation_points;
        const int& np = n + 1;
        const int& nptm = npt - np;
        const int& nh = n * np / 2;

        const Array_Real_Vector work1 = Array_Real_Vector(n);
        const Array_Real_Vector work2 = Array_Real_Vector(npt);
        const Array_Real_Vector work3 = Array_Real_Vector(npt);

        double cauchy = std::numeric_limits<double>::quiet_NaN();
        double alpha = std::numeric_limits<double>::quiet_NaN();
        double dsq = std::numeric_limits<double>::quiet_NaN();
        double crvmin;

        // Set some constants.
        // Parameter adjustments

        // Function Body

        // The call of PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ, // BMAT and ZMAT for the first iteration, with the corresponding values of
        // of NF and KOPT, which are the number of calls of CALFUN so far and the
        // index of the interpolation point at the trust region centre. Then the
        // initial XOPT is set too. The branch to label 720 occurs if MAXFUN is
        // less than NPT. GOPT will be updated if KOPT is different from KBASE.

        trust_region_center_interpolation_point_index = 0;

        prelim(lower_bound, upper_bound);
        double xoptsq = ZERO;
        for (int i{}; i < n; i++) 
        {
            trust_region_center_offset.set_entry(i, interpolation_points.get_entry(trust_region_center_interpolation_point_index, i));
            // Computing 2nd power
            const double delta_one = trust_region_center_offset.get_entry(i);
            xoptsq += delta_one * delta_one;
        }
        double fsave = f_at_interpolation_points.get_entry(0);
        const int& kbase = 0;

        // Complete the settings that are required for the iterative procedure.

        int ntrits = 0;
        int itest = 0;
        int knew = 0;
        int nfsav = get_evaluations();
        double rho = initial_trust_region_radius;
        double delta = rho;
        double diffa = ZERO;
        double diffb = ZERO;
        double diffc = ZERO;
        double f = ZERO;
        double beta = ZERO;
        double adelt = ZERO;
        double denom = ZERO;
        double ratio = ZERO;
        double dnorm = ZERO;
        double scaden;
        double biglsq;
        double distsq = ZERO;

        // Update GOPT if necessary before the first iteration and after each
        // call of RESCUE that makes a call of CALFUN.

        int state = 20;
        for(;;) 
        {
        switch (state) { // NOPMD - the reference algorithm is as complex as this, we simply ported it from Fortran with minimal changes
        case 20: 
        {
            if (trust_region_center_interpolation_point_index != kbase) 
            {
                int ih = 0;
                for (int j{}; j < n; j++) 
                {
                    for (int i{}; i <= j; i++) 
                    {
                        if (i < j) 
                        {
                            gradient_at_trust_region_center.set_entry(j, gradient_at_trust_region_center.get_entry(j) + model_second_derivatives_values.get_entry(ih) * trust_region_center_offset.get_entry(i));
                        }
                        gradient_at_trust_region_center.set_entry(i, gradient_at_trust_region_center.get_entry(i) + model_second_derivatives_values.get_entry(ih) * trust_region_center_offset.get_entry(j));
                        ih++;
                    }
                }
                if (get_evaluations() > npt) 
                {
                    for (int k{}; k < npt; k++) 
                    {
                        double temp = ZERO;
                        for (int j{}; j < n; j++) 
                        {
                            temp += interpolation_points.get_entry(k, j) * trust_region_center_offset.get_entry(j);
                        }
                        temp *= model_second_derivatives_parameters.get_entry(k);
                        for (int i{}; i < n; i++) 
                        {
                            gradient_at_trust_region_center.set_entry(i, gradient_at_trust_region_center.get_entry(i) + temp * interpolation_points.get_entry(k, i));
                        }
                    }
                    // throw Path_isExploredException(); // XXX
                }
            }

            // Generate the next point in the trust region that provides a small value
            // of the quadratic model subject to the constraints on the variables.
            // The int NTRITS is set to the number "trust region" iterations that
            // have occurred since the last "alternative" iteration. If the length
            // of XNEW-XOPT is less than HALF*RHO, however, then there is a branch to
            // label 650 or 680 with NTRITS=-1, instead of calculating F at XNEW.

        }
        case 60: 
        {
            const Array_Real_Vector gnew = Array_Real_Vector(n);
            const Array_Real_Vector xbdi = Array_Real_Vector(n);
            const Array_Real_Vector s = Array_Real_Vector(n);
            const Array_Real_Vector hs = Array_Real_Vector(n);
            const Array_Real_Vector hred = Array_Real_Vector(n);

            const std::vector<double> dsq_crvmin = trsbox(delta, gnew, xbdi, s, hs, hred);
            dsq = dsq_crvmin[0];
            crvmin = dsq_crvmin[1];

            // Computing MIN
            double delta_one = delta;
            double delta_two = std::sqrt(dsq);
            dnorm = std::min(delta_one, delta_two);
            if (dnorm < HALF * rho) 
            {
                ntrits = -1;
                // Computing 2nd power
                delta_one = TEN * rho;
                distsq = delta_one * delta_one;
                if (get_evaluations() <= nfsav + 2) 
                {
                    state = 650; break;
                }

                // The following choice between labels 650 and 680 depends on whether or
                // not our work with the current RHO seems to be complete. Either RHO is
                // decreased or termination occurs if the errors in the quadratic model at
                // the last three interpolation points compare favourably with predictions
                // of likely improvements to the model within distance HALF*RHO of XOPT.

                // Computing MAX
                delta_one = std::max(diffa, diffb);
                const double errbig = std::max(delta_one, diffc);
                const double frhosq = rho * ONE_OVER_EIGHT * rho;
                if (crvmin > ZERO &&
                    errbig > frhosq * crvmin) 
                    {
                    state = 650; break;
                }
                const double bdtol = errbig / rho;
                for (int j{}; j < n; j++) 
                {
                    double bdtest = bdtol;
                    if (new_point.get_entry(j) == lower_difference.get_entry(j)) 
                    {
                        bdtest = work1.get_entry(j);
                    }
                    if (new_point.get_entry(j) == upper_difference.get_entry(j)) 
                    {
                        bdtest = -work1.get_entry(j);
                    }
                    if (bdtest < bdtol) 
                    {
                        double curv = model_second_derivatives_values.get_entry((j + j * j) / 2);
                        for (int k{}; k < npt; k++) 
                        {
                            // Computing 2nd power
                            const double d1 = interpolation_points.get_entry(k, j);
                            curv += model_second_derivatives_parameters.get_entry(k) * (d1 * d1);
                        }
                        bdtest += HALF * curv * rho;
                        if (bdtest < bdtol) 
                        {
                            break;
                        }
                        // throw Path_isExploredException(); // XXX
                    }
                }
                state = 680; break;
            }
            ++ntrits;

            // Severe cancellation is likely to occur if XOPT is too far from XBASE.
            // If the following test holds, then XBASE is shifted so that XOPT becomes
            // zero. The appropriate changes are made to BMAT and to the second
            // derivatives of the current model, beginning with the changes to BMAT
            // that do not depend on ZMAT. VLAG is used temporarily for working space.

        }
        case 90: 
        {
            if (dsq <= xoptsq * ONE_OVER_A_THOUSAND) 
            {
                const double fracsq = xoptsq * ONE_OVER_FOUR;
                double sumpq = ZERO;
                // const Real_Vector sum_vector
                //     = Array_Real_Vector(npt, -HALF * xoptsq).add(interpolation_points.operate(trust_region_center));
                for (int k{}; k < npt; k++) 
                {
                    sumpq += model_second_derivatives_parameters.get_entry(k);
                    double sum = -HALF * xoptsq;
                    for (int i{}; i < n; i++) 
                    {
                        sum += interpolation_points.get_entry(k, i) * trust_region_center_offset.get_entry(i);
                    }
                    // sum = sum_vector.get_entry(k); // XXX "test_ackley" and "test_diff_pow" fail.
                    work2.set_entry(k, sum);
                    const double temp = fracsq - HALF * sum;
                    for (int i{}; i < n; i++) 
                    {
                        work1.set_entry(i, b_matrix.get_entry(k, i));
                        lagrange_values_at_new_point.set_entry(i, sum * interpolation_points.get_entry(k, i) + temp * trust_region_center_offset.get_entry(i));
                        const int ip = npt + i;
                        for (int j{}; j <= i; j++) 
                        {
                            b_matrix.set_entry(ip, j, b_matrix.get_entry(ip, j)
                                          + work1.get_entry(i) * lagrange_values_at_new_point.get_entry(j)
                                          + lagrange_values_at_new_point.get_entry(i) * work1.get_entry(j));
                        }
                    }
                }

                // Then the revisions of BMAT that depend on ZMAT are calculated.

                for (const int& m = 0; m < nptm; m++) 
                {
                    double sumz = ZERO;
                    double sumw = ZERO;
                    for (int k{}; k < npt; k++) 
                    {
                        sumz += z_matrix.get_entry(k, m);
                        lagrange_values_at_new_point.set_entry(k, work2.get_entry(k) * z_matrix.get_entry(k, m));
                        sumw += lagrange_values_at_new_point.get_entry(k);
                    }
                    for (int j{}; j < n; j++) 
                    {
                        double sum = (fracsq * sumz - HALF * sumw) * trust_region_center_offset.get_entry(j);
                        for (int k{}; k < npt; k++) 
                        {
                            sum += lagrange_values_at_new_point.get_entry(k) * interpolation_points.get_entry(k, j);
                        }
                        work1.set_entry(j, sum);
                        for (int k{}; k < npt; k++) 
                        {
                            b_matrix.set_entry(k, j, b_matrix.get_entry(k, j)
                                          + sum * z_matrix.get_entry(k, m));
                        }
                    }
                    for (int i{}; i < n; i++) 
                    {
                        const int ip = i + npt;
                        const double temp = work1.get_entry(i);
                        for (int j{}; j <= i; j++) 
                        {
                            b_matrix.set_entry(ip, j, b_matrix.get_entry(ip, j)
                                          + temp * work1.get_entry(j));
                        }
                    }
                }

                // The following instructions complete the shift, including the changes
                // to the second derivative parameters of the quadratic model.

                int ih = 0;
                for (int j{}; j < n; j++) 
                {
                    work1.set_entry(j, -HALF * sumpq * trust_region_center_offset.get_entry(j));
                    for (int k{}; k < npt; k++) 
                    {
                        work1.set_entry(j, work1.get_entry(j) + model_second_derivatives_parameters.get_entry(k) * interpolation_points.get_entry(k, j));
                        interpolation_points.set_entry(k, j, interpolation_points.get_entry(k, j) - trust_region_center_offset.get_entry(j));
                    }
                    for (int i{}; i <= j; i++) 
                    {
                         model_second_derivatives_values.set_entry(ih, model_second_derivatives_values.get_entry(ih)
                                    + work1.get_entry(i) * trust_region_center_offset.get_entry(j)
                                    + trust_region_center_offset.get_entry(i) * work1.get_entry(j));
                        b_matrix.set_entry(npt + i, j, b_matrix.get_entry(npt + j, i));
                        ih++;
                    }
                }
                for (int i{}; i < n; i++) 
                {
                    origin_shift.set_entry(i, origin_shift.get_entry(i) + trust_region_center_offset.get_entry(i));
                    new_point.set_entry(i, new_point.get_entry(i) - trust_region_center_offset.get_entry(i));
                    lower_difference.set_entry(i, lower_difference.get_entry(i) - trust_region_center_offset.get_entry(i));
                    upper_difference.set_entry(i, upper_difference.get_entry(i) - trust_region_center_offset.get_entry(i));
                    trust_region_center_offset.set_entry(i, ZERO);
                }
                xoptsq = ZERO;
            }
            if (ntrits == 0) 
            {
                state = 210; break;
            }
            state = 230; break;

            // XBASE is also moved to XOPT by a call of RESCUE. This calculation is
            // more expensive than the previous shift, because matrices BMAT and
            // ZMAT are generated from scratch, which may include the replacement of
            // interpolation points whose positions seem to be causing near linear
            // dependence in the interpolation conditions. Therefore RESCUE is called
            // only if rounding errors have reduced by at least a factor of two the
            // denominator of the formula for updating the H matrix. It provides a
            // useful safeguard, but is not invoked in most applications of BOBYQA.

        }
        case 210: 
        {
            // Pick two alternative vectors of variables, relative to XBASE, that
            // are suitable as positions of the KNEW-th interpolation point.
            // Firstly, XNEW is set to the point on a line through XOPT and another
            // interpolation point that minimizes the predicted value of the next
            // denominator, subject to ||XNEW - XOPT|| .LEQ. ADELT and to the SL
            // and SU bounds. Secondly, XALT is set to the best feasible point on
            // a constrained version of the Cauchy step of the KNEW-th Lagrange
            // function, the corresponding value of the square of this function
            // being returned in CAUCHY. The choice between these alternatives is
            // going to be made when the denominator is calculated.

            const std::vector<double> alpha_cauchy = altmov(knew, adelt);
            alpha = alpha_cauchy[0];
            cauchy = alpha_cauchy[1];

            for (int i{}; i < n; i++) 
            {
                trial_step_point.set_entry(i, new_point.get_entry(i) - trust_region_center_offset.get_entry(i));
            }

            // Calculate VLAG and BETA for the current choice of D. The scalar
            // product of D with XPT(K,.) is going to be held in W(NPT+K) for
            // use when VQUAD is calculated.

        }
        case 230: 
        {
            for (int k{}; k < npt; k++) 
            {
                double suma = ZERO;
                double sumb = ZERO;
                double sum = ZERO;
                for (int j{}; j < n; j++) 
                {
                    suma += interpolation_points.get_entry(k, j) * trial_step_point.get_entry(j);
                    sumb += interpolation_points.get_entry(k, j) * trust_region_center_offset.get_entry(j);
                    sum += b_matrix.get_entry(k, j) * trial_step_point.get_entry(j);
                }
                work3.set_entry(k, suma * (HALF * suma + sumb));
                lagrange_values_at_new_point.set_entry(k, sum);
                work2.set_entry(k, suma);
            }
            beta = ZERO;
            for (const int& m = 0; m < nptm; m++) 
            {
                double sum = ZERO;
                for (int k{}; k < npt; k++) 
                {
                    sum += z_matrix.get_entry(k, m) * work3.get_entry(k);
                }
                beta -= sum * sum;
                for (int k{}; k < npt; k++) 
                {
                    lagrange_values_at_new_point.set_entry(k, lagrange_values_at_new_point.get_entry(k) + sum * z_matrix.get_entry(k, m));
                }
            }
            dsq = ZERO;
            double bsum = ZERO;
            double dx = ZERO;
            for (int j{}; j < n; j++) 
            {
                // Computing 2nd power
                const double d1 = trial_step_point.get_entry(j);
                dsq += d1 * d1;
                double sum = ZERO;
                for (int k{}; k < npt; k++) 
                {
                    sum += work3.get_entry(k) * b_matrix.get_entry(k, j);
                }
                bsum += sum * trial_step_point.get_entry(j);
                const int jp = npt + j;
                for (int i{}; i < n; i++) 
                {
                    sum += b_matrix.get_entry(jp, i) * trial_step_point.get_entry(i);
                }
                lagrange_values_at_new_point.set_entry(jp, sum);
                bsum += sum * trial_step_point.get_entry(j);
                dx += trial_step_point.get_entry(j) * trust_region_center_offset.get_entry(j);
            }

            beta = dx * dx + dsq * (xoptsq + dx + dx + HALF * dsq) + beta - bsum; // Original
            // beta += dx * dx + dsq * (xoptsq + dx + dx + HALF * dsq) - bsum; // XXX "test_ackley" and "test_diff_pow" fail.
            // beta = dx * dx + dsq * (xoptsq + 2 * dx + HALF * dsq) + beta - bsum; // XXX "test_diff_pow" fails.

            lagrange_values_at_new_point.set_entry(trust_region_center_interpolation_point_index, lagrange_values_at_new_point.get_entry(trust_region_center_interpolation_point_index) + ONE);

            // If NTRITS is zero, the denominator may be increased by replacing
            // the step D of ALTMOV by a Cauchy step. Then RESCUE may be called if
            // rounding errors have damaged the chosen denominator.

            if (ntrits == 0) 
            {
                // Computing 2nd power
                const double d1 = lagrange_values_at_new_point.get_entry(knew);
                denom = d1 * d1 + alpha * beta;
                if (denom < cauchy && cauchy > ZERO) 
                {
                    for (int i{}; i < n; i++) 
                    {
                        new_point.set_entry(i, alternative_new_point.get_entry(i));
                        trial_step_point.set_entry(i, new_point.get_entry(i) - trust_region_center_offset.get_entry(i));
                    }
                    cauchy = ZERO; // XXX Useful statement?
                    state = 230; break;
                }
                // Alternatively, if NTRITS is positive, then set KNEW to the index of
                // the next interpolation point to be deleted to make room for a trust
                // region step. Again RESCUE may be called if rounding errors have damaged_
                // the chosen denominator, which is the reason for attempting to select
                // KNEW before calculating the next value of the objective function.

            }
else 
            {
                const double delsq = delta * delta;
                scaden = ZERO;
                biglsq = ZERO;
                knew = 0;
                for (int k{}; k < npt; k++) 
                {
                    if (k == trust_region_center_interpolation_point_index) 
                    {
                        continue;
                    }
                    double hdiag = ZERO;
                    for (const int& m = 0; m < nptm; m++) 
                    {
                        // Computing 2nd power
                        const double d1 = z_matrix.get_entry(k, m);
                        hdiag += d1 * d1;
                    }
                    // Computing 2nd power
                    const double d2 = lagrange_values_at_new_point.get_entry(k);
                    const double den = beta * hdiag + d2 * d2;
                    distsq = ZERO;
                    for (int j{}; j < n; j++) 
                    {
                        // Computing 2nd power
                        const double d3 = interpolation_points.get_entry(k, j) - trust_region_center_offset.get_entry(j);
                        distsq += d3 * d3;
                    }
                    // Computing MAX
                    // Computing 2nd power
                    const double d4 = distsq / delsq;
                    const double temp = std::max(ONE, d4 * d4);
                    if (temp * den > scaden) 
                    {
                        scaden = temp * den;
                        knew = k;
                        denom = den;
                    }
                    // Computing MAX
                    // Computing 2nd power
                    const double d5 = lagrange_values_at_new_point.get_entry(k);
                    biglsq = std::max(biglsq, temp * (d5 * d5));
                }
            }

            // Put the variables for the next calculation of the objective function
            //   in XNEW, with any adjustments for the bounds.

            // Calculate the value of the objective function at XBASE+XNEW, unless
            //   the limit on the number of calculations of F has been reached.

        }
        case 360: 
        {
            for (int i{}; i < n; i++) 
            {
                // Computing MIN
                // Computing MAX
                const double d3 = lower_bound[i];
                const double d4 = origin_shift.get_entry(i) + new_point.get_entry(i);
                const double d1 = std::max(d3, d4);
                const double d2 = upper_bound[i];
                current_best.set_entry(i, std::min(d1, d2));
                if (new_point.get_entry(i) == lower_difference.get_entry(i)) 
                {
                    current_best.set_entry(i, lower_bound[i]);
                }
                if (new_point.get_entry(i) == upper_difference.get_entry(i)) 
                {
                    current_best.set_entry(i, upper_bound[i]);
                }
            }

            f = compute_objective_value(current_best.to_array());

            if (!is_minimize) 
            {
                f = -f;
            }
            if (ntrits == -1) 
            {
                fsave = f;
                state = 720; break;
            }

            // Use the quadratic model to predict the change in F due to the step D, //   and set DIFF to the error of this prediction.

            const double fopt = f_at_interpolation_points.get_entry(trust_region_center_interpolation_point_index);
            double vquad = ZERO;
            int ih = 0;
            for (int j{}; j < n; j++) 
            {
                vquad += trial_step_point.get_entry(j) * gradient_at_trust_region_center.get_entry(j);
                for (int i{}; i <= j; i++) 
                {
                    double temp = trial_step_point.get_entry(i) * trial_step_point.get_entry(j);
                    if (i == j) 
                    {
                        temp *= HALF;
                    }
                    vquad += model_second_derivatives_values.get_entry(ih) * temp;
                    ih++;
               }
            }
            for (int k{}; k < npt; k++) 
            {
                // Computing 2nd power
                const double d1 = work2.get_entry(k);
                const double d2 = d1 * d1; // "d1" must be squared first to prevent test failures.
                vquad += HALF * model_second_derivatives_parameters.get_entry(k) * d2;
            }
            const double diff = f - fopt - vquad;
            diffc = diffb;
            diffb = diffa;
            diffa = std::abs(diff);
            if (dnorm > rho) 
            {
                nfsav = get_evaluations();
            }

            // Pick the next value of DELTA after a trust region step.

            if (ntrits > 0) 
            {
                if (vquad >= ZERO) 
                {
                    throw Math_Illegal_State_Exception(Localized_Optim_Formats.TRUST_REGION_STEP_FAILED, vquad);
                }
                ratio = (f - fopt) / vquad;
                const double h_delta = HALF * delta;
                if (ratio <= ONE_OVER_TEN) 
                {
                    // Computing MIN
                    delta = std::min(h_delta, dnorm);
                }
else if (ratio <= .7) 
                {
                    // Computing MAX
                    delta = std::max(h_delta, dnorm);
                }
else 
                {
                    // Computing MAX
                    delta = std::max(h_delta, 2 * dnorm);
                }
                if (delta <= rho * 1.5) 
                {
                    delta = rho;
                }

                // Recalculate KNEW and DENOM if the F is less than FOPT.

                if (f < fopt) 
                {
                    const int& ksav = knew;
                    const double densav = denom;
                    const double delsq = delta * delta;
                    scaden = ZERO;
                    biglsq = ZERO;
                    knew = 0;
                    for (int k{}; k < npt; k++) 
                    {
                        double hdiag = ZERO;
                        for (const int& m = 0; m < nptm; m++) 
                        {
                            // Computing 2nd power
                            const double d1 = z_matrix.get_entry(k, m);
                            hdiag += d1 * d1;
                        }
                        // Computing 2nd power
                        const double d1 = lagrange_values_at_new_point.get_entry(k);
                        const double den = beta * hdiag + d1 * d1;
                        distsq = ZERO;
                        for (int j{}; j < n; j++) 
                        {
                            // Computing 2nd power
                            const double d2 = interpolation_points.get_entry(k, j) - new_point.get_entry(j);
                            distsq += d2 * d2;
                        }
                        // Computing MAX
                        // Computing 2nd power
                        const double d3 = distsq / delsq;
                        const double temp = std::max(ONE, d3 * d3);
                        if (temp * den > scaden) 
                        {
                            scaden = temp * den;
                            knew = k;
                            denom = den;
                        }
                        // Computing MAX
                        // Computing 2nd power
                        const double d4 = lagrange_values_at_new_point.get_entry(k);
                        const double d5 = temp * (d4 * d4);
                        biglsq = std::max(biglsq, d5);
                    }
                    if (scaden <= HALF * biglsq) 
                    {
                        knew = ksav;
                        denom = densav;
                    }
                }
            }

            // Update BMAT and ZMAT, so that the KNEW-th interpolation point can be
            // moved. Also update the second derivative terms of the model.

            update(beta, denom, knew);

            ih = 0;
            const double pqold = model_second_derivatives_parameters.get_entry(knew);
            model_second_derivatives_parameters.set_entry(knew, ZERO);
            for (int i{}; i < n; i++) 
            {
                const double temp = pqold * interpolation_points.get_entry(knew, i);
                for (int j{}; j <= i; j++) 
                {
                    model_second_derivatives_values.set_entry(ih, model_second_derivatives_values.get_entry(ih) + temp * interpolation_points.get_entry(knew, j));
                    ih++;
                }
            }
            for (const int& m = 0; m < nptm; m++) 
            {
                const double temp = diff * z_matrix.get_entry(knew, m);
                for (int k{}; k < npt; k++) 
                {
                    model_second_derivatives_parameters.set_entry(k, model_second_derivatives_parameters.get_entry(k) + temp * z_matrix.get_entry(k, m));
                }
            }

            // Include the interpolation point, and make the changes to GOPT at
            // the old XOPT that are caused by the updating of the quadratic model.

            f_at_interpolation_points.set_entry(knew,  f);
            for (int i{}; i < n; i++) 
            {
                interpolation_points.set_entry(knew, i, new_point.get_entry(i));
                work1.set_entry(i, b_matrix.get_entry(knew, i));
            }
            for (int k{}; k < npt; k++) 
            {
                double suma = ZERO;
                for (const int& m = 0; m < nptm; m++) 
                {
                    suma += z_matrix.get_entry(knew, m) * z_matrix.get_entry(k, m);
                }
                double sumb = ZERO;
                for (int j{}; j < n; j++) 
                {
                    sumb += interpolation_points.get_entry(k, j) * trust_region_center_offset.get_entry(j);
                }
                const double temp = suma * sumb;
                for (int i{}; i < n; i++) 
                {
                    work1.set_entry(i, work1.get_entry(i) + temp * interpolation_points.get_entry(k, i));
                }
            }
            for (int i{}; i < n; i++) 
            {
                gradient_at_trust_region_center.set_entry(i, gradient_at_trust_region_center.get_entry(i) + diff * work1.get_entry(i));
            }

            // Update XOPT, GOPT and KOPT if the calculated F is less than FOPT.

            if (f < fopt) 
            {
                trust_region_center_interpolation_point_index = knew;
                xoptsq = ZERO;
                ih = 0;
                for (int j{}; j < n; j++) 
                {
                    trust_region_center_offset.set_entry(j, new_point.get_entry(j));
                    // Computing 2nd power
                    const double d1 = trust_region_center_offset.get_entry(j);
                    xoptsq += d1 * d1;
                    for (int i{}; i <= j; i++) 
                    {
                        if (i < j) 
                        {
                            gradient_at_trust_region_center.set_entry(j, gradient_at_trust_region_center.get_entry(j) + model_second_derivatives_values.get_entry(ih) * trial_step_point.get_entry(i));
                        }
                        gradient_at_trust_region_center.set_entry(i, gradient_at_trust_region_center.get_entry(i) + model_second_derivatives_values.get_entry(ih) * trial_step_point.get_entry(j));
                        ih++;
                    }
                }
                for (int k{}; k < npt; k++) 
                {
                    double temp = ZERO;
                    for (int j{}; j < n; j++) 
                    {
                        temp += interpolation_points.get_entry(k, j) * trial_step_point.get_entry(j);
                    }
                    temp *= model_second_derivatives_parameters.get_entry(k);
                    for (int i{}; i < n; i++) 
                    {
                        gradient_at_trust_region_center.set_entry(i, gradient_at_trust_region_center.get_entry(i) + temp * interpolation_points.get_entry(k, i));
                    }
                }
            }

            // Calculate the parameters of the least Frobenius norm interpolant to
            // the current data, the gradient of this interpolant at XOPT being put
            // into VLAG(NPT+I), I=1,2,...,N.

            if (ntrits > 0) 
            {
                for (int k{}; k < npt; k++) 
                {
                    lagrange_values_at_new_point.set_entry(k, f_at_interpolation_points.get_entry(k) - f_at_interpolation_points.get_entry(trust_region_center_interpolation_point_index));
                    work3.set_entry(k, ZERO);
                }
                for (int j{}; j < nptm; j++) 
                {
                    double sum = ZERO;
                    for (int k{}; k < npt; k++) 
                    {
                        sum += z_matrix.get_entry(k, j) * lagrange_values_at_new_point.get_entry(k);
                    }
                    for (int k{}; k < npt; k++) 
                    {
                        work3.set_entry(k, work3.get_entry(k) + sum * z_matrix.get_entry(k, j));
                    }
                }
                for (int k{}; k < npt; k++) 
                {
                    double sum = ZERO;
                    for (int j{}; j < n; j++) 
                    {
                        sum += interpolation_points.get_entry(k, j) * trust_region_center_offset.get_entry(j);
                    }
                    work2.set_entry(k, work3.get_entry(k));
                    work3.set_entry(k, sum * work3.get_entry(k));
                }
                double gqsq = ZERO;
                double gisq = ZERO;
                for (int i{}; i < n; i++) 
                {
                    double sum = ZERO;
                    for (int k{}; k < npt; k++) 
                    {
                        sum += b_matrix.get_entry(k, i) *
                            lagrange_values_at_new_point.get_entry(k) + interpolation_points.get_entry(k, i) * work3.get_entry(k);
                    }
                    if (trust_region_center_offset.get_entry(i) == lower_difference.get_entry(i)) 
                    {
                        // Computing MIN
                        // Computing 2nd power
                        const double d1 = std::min(ZERO, gradient_at_trust_region_center.get_entry(i));
                        gqsq += d1 * d1;
                        // Computing 2nd power
                        const double d2 = std::min(ZERO, sum);
                        gisq += d2 * d2;
                    }
else if (trust_region_center_offset.get_entry(i) == upper_difference.get_entry(i)) 
                    {
                        // Computing MAX
                        // Computing 2nd power
                        const double d1 = std::max(ZERO, gradient_at_trust_region_center.get_entry(i));
                        gqsq += d1 * d1;
                        // Computing 2nd power
                        const double d2 = std::max(ZERO, sum);
                        gisq += d2 * d2;
                    }
else 
                    {
                        // Computing 2nd power
                        const double d1 = gradient_at_trust_region_center.get_entry(i);
                        gqsq += d1 * d1;
                        gisq += sum * sum;
                    }
                    lagrange_values_at_new_point.set_entry(npt + i, sum);
                }

                // Test whether to replace the quadratic model by the least Frobenius
                // norm interpolant, making the replacement if the test is satisfied.

                ++itest;
                if (gqsq < TEN * gisq) 
                {
                    itest = 0;
                }
                if (itest >= 3) 
                {
                    const int max = std::max(npt, nh);
                    for (int i{}; i < max; i++) 
                    {
                        if (i < n) 
                        {
                            gradient_at_trust_region_center.set_entry(i, lagrange_values_at_new_point.get_entry(npt + i));
                        }
                        if (i < npt) 
                        {
                            model_second_derivatives_parameters.set_entry(i, work2.get_entry(i));
                        }
                        if (i < nh) 
                        {
                            model_second_derivatives_values.set_entry(i, ZERO);
                        }
                        itest = 0;
                    }
                }
            }

            // If a trust region step has provided a sufficient decrease in F, then
            // branch for another trust region calculation. The case NTRITS=0 occurs
            // when the interpolation point was reached by an alternative step.

            if (ntrits == 0) 
            {
                state = 60; break;
            }
            if (f <= fopt + ONE_OVER_TEN * vquad) 
            {
                state = 60; break;
            }

            // Alternatively, find out if the interpolation points are close enough
            //   to the best point so far.

            // Computing MAX
            // Computing 2nd power
            const double d1 = TWO * delta;
            // Computing 2nd power
            const double d2 = TEN * rho;
            distsq = std::max(d1 * d1, d2 * d2);
        }
        case 650: 
        {
            knew = -1;
            for (int k{}; k < npt; k++) 
            {
                double sum = ZERO;
                for (int j{}; j < n; j++) 
                {
                    // Computing 2nd power
                    const double d1 = interpolation_points.get_entry(k, j) - trust_region_center_offset.get_entry(j);
                    sum += d1 * d1;
                }
                if (sum > distsq) 
                {
                    knew = k;
                    distsq = sum;
                }
            }

            // If KNEW is positive, then ALTMOV finds alternative positions for
            // the KNEW-th interpolation point within distance ADELT of XOPT. It is
            // reached via label 90. Otherwise, there is a branch to label 60 for
            // another trust region iteration, unless the calculations with the
            // current RHO are complete.

            if (knew >= 0) 
            {
                const double dist = std::sqrt(distsq);
                if (ntrits == -1) 
                {
                    // Computing MIN
                    delta = std::min(ONE_OVER_TEN * delta, HALF * dist);
                    if (delta <= rho * 1.5) 
                    {
                        delta = rho;
                    }
                }
                ntrits = 0;
                // Computing MAX
                // Computing MIN
                const double d1 = std::min(ONE_OVER_TEN * dist, delta);
                adelt = std::max(d1, rho);
                dsq = adelt * adelt;
                state = 90; break;
            }
            if (ntrits == -1) 
            {
                state = 680; break;
            }
            if (ratio > ZERO) 
            {
                state = 60; break;
            }
            if (std::max(delta, dnorm) > rho) 
            {
                state = 60; break;
            }

            // The calculations with the current value of RHO are complete. Pick the
            //   next values of RHO and DELTA.
        }
        case 680: 
        {
            if (rho > stopping_trust_region_radius) 
            {
                delta = HALF * rho;
                ratio = rho / stopping_trust_region_radius;
                if (ratio <= SIXTEEN) 
                {
                    rho = stopping_trust_region_radius;
                }
else if (ratio <= TWO_HUNDRED_FIFTY) 
                {
                    rho = std::sqrt(ratio) * stopping_trust_region_radius;
                }
else 
                {
                    rho *= ONE_OVER_TEN;
                }
                delta = std::max(delta, rho);
                ntrits = 0;
                nfsav = get_evaluations();
                state = 60; break;
            }

            // Return from the calculation, after another Newton-_Raphson step, if
            //   it is too short to have been tried before.

            if (ntrits == -1) 
            {
                state = 360; break;
            }
        }
        case 720: 
        {
            if (f_at_interpolation_points.get_entry(trust_region_center_interpolation_point_index) <= fsave) 
            {
                for (int i{}; i < n; i++) 
                {
                    // Computing MIN
                    // Computing MAX
                    const double d3 = lower_bound[i];
                    const double d4 = origin_shift.get_entry(i) + trust_region_center_offset.get_entry(i);
                    const double d1 = std::max(d3, d4);
                    const double d2 = upper_bound[i];
                    current_best.set_entry(i, std::min(d1, d2));
                    if (trust_region_center_offset.get_entry(i) == lower_difference.get_entry(i)) 
                    {
                        current_best.set_entry(i, lower_bound[i]);
                    }
                    if (trust_region_center_offset.get_entry(i) == upper_difference.get_entry(i)) 
                    {
                        current_best.set_entry(i, upper_bound[i]);
                    }
                }
                f = f_at_interpolation_points.get_entry(trust_region_center_interpolation_point_index);
            }
            return f;
        }
        default: 
        {
            throw Math_Illegal_State_Exception(Localized_Core_Formats.SIMPLE_MESSAGE, "bobyqb");
        }}}
    } // bobyqb

    // ----------------------------------------------------------------------------------------

    /**
     *     The arguments N, NPT, XPT, XOPT, BMAT, ZMAT, NDIM, SL and SU all have
     *       the same meanings as the corresponding arguments of BOBYQB.
     *     KOPT is the index of the optimal interpolation point.
     *     KNEW is the index of the interpolation point that is going to be moved.
     *     ADELT is the current trust region bound.
     *     XNEW will be set to a suitable position for the interpolation point
     *       XPT(KNEW,.). Specifically, it satisfies the SL, SU and trust region
     *       bounds and it should provide a large denominator in the next call of
     *       UPDATE. The step XNEW-XOPT from XOPT is restricted to moves along the
     *       straight lines through XOPT and another interpolation point.
     *     XALT also provides a large value of the modulus of the KNEW-th Lagrange
     *       function subject to the constraints that have been mentioned, its main
     *       difference from XNEW being that XALT-XOPT is a constrained version of
     *       the Cauchy step within the trust region. An exception is that XALT is
     *       not calculated if all components of GLAG (see below) are zero.
     *     ALPHA will be set to the KNEW-th diagonal element of the H matrix.
     *     CAUCHY will be set to the square of the KNEW-th Lagrange function at
     *       the step XALT-XOPT from XOPT for the vector XALT that is returned, *       except that CAUCHY is set to zero if XALT is not calculated.
     *     GLAG is a working space vector of length N for the gradient of the
     *       KNEW-th Lagrange function at XOPT.
     *     HCOL is a working space vector of length NPT for the second derivative
     *       coefficients of the KNEW-th Lagrange function.
     *     W is a working space vector of length 2N that is going to hold the
     *       constrained Cauchy step from XOPT of the Lagrange function, followed
     *       by the downhill version of XALT when the uphill step is calculated.
     *
     *     Set the first NPT components of W to the leading elements of the
     *     KNEW-th column of the H matrix.
     * @param knew
     * @param adelt
     */
    private std::vector<double> altmov(const int& knew, double adelt) 
    {

        const int n = current_best.get_dimension();
        const int& npt = number_of_interpolation_points;

        const Array_Real_Vector glag = Array_Real_Vector(n);
        const Array_Real_Vector hcol = Array_Real_Vector(npt);

        const Array_Real_Vector work1 = Array_Real_Vector(n);
        const Array_Real_Vector work2 = Array_Real_Vector(n);

        for (int k{}; k < npt; k++) 
        {
            hcol.set_entry(k, ZERO);
        }
        const int max = npt - n - 1;
        for (int j{}; j < max; j++) 
        {
            const double tmp = z_matrix.get_entry(knew, j);
            for (int k{}; k < npt; k++) 
            {
                hcol.set_entry(k, hcol.get_entry(k) + tmp * z_matrix.get_entry(k, j));
            }
        }
        const double& alpha = hcol.get_entry(knew);
        const double ha = HALF * alpha;

        // Calculate the gradient of the KNEW-th Lagrange function at XOPT.

        for (int i{}; i < n; i++) 
        {
            glag.set_entry(i, b_matrix.get_entry(knew, i));
        }
        for (int k{}; k < npt; k++) 
        {
            double tmp = ZERO;
            for (int j{}; j < n; j++) 
            {
                tmp += interpolation_points.get_entry(k, j) * trust_region_center_offset.get_entry(j);
            }
            tmp *= hcol.get_entry(k);
            for (int i{}; i < n; i++) 
            {
                glag.set_entry(i, glag.get_entry(i) + tmp * interpolation_points.get_entry(k, i));
            }
        }

        // Search for a large denominator along the straight lines through XOPT
        // and another interpolation point. SLBD and SUBD will be lower and upper
        // bounds on the step along each of these lines in turn. PREDSQ will be
        // set to the square of the predicted denominator for each line. PRESAV
        // will be set to the largest admissible value of PREDSQ that occurs.

        double presav = ZERO;
        double step = std::numeric_limits<double>::quiet_NaN();
        int ksav = 0;
        int ibdsav = 0;
        double stpsav = 0;
        for (int k{}; k < npt; k++) 
        {
            if (k == trust_region_center_interpolation_point_index) 
            {
                continue;
            }
            double dderiv = ZERO;
            double distsq = ZERO;
            for (int i{}; i < n; i++) 
            {
                const double tmp = interpolation_points.get_entry(k, i) - trust_region_center_offset.get_entry(i);
                dderiv += glag.get_entry(i) * tmp;
                distsq += tmp * tmp;
            }
            double subd = adelt / std::sqrt(distsq);
            double slbd = -subd;
            int ilbd = 0;
            int iubd = 0;
            const double sumin = std::min(ONE, subd);

            // Revise SLBD and SUBD if necessary because of the bounds in SL and SU.

            for (int i{}; i < n; i++) 
            {
                const double tmp = interpolation_points.get_entry(k, i) - trust_region_center_offset.get_entry(i);
                if (tmp > ZERO) 
                {
                    if (slbd * tmp < lower_difference.get_entry(i) - trust_region_center_offset.get_entry(i)) 
                    {
                        slbd = (lower_difference.get_entry(i) - trust_region_center_offset.get_entry(i)) / tmp;
                        ilbd = -i - 1;
                    }
                    if (subd * tmp > upper_difference.get_entry(i) - trust_region_center_offset.get_entry(i)) 
                    {
                        // Computing MAX
                        subd = std::max(sumin, (upper_difference.get_entry(i) - trust_region_center_offset.get_entry(i)) / tmp);
                        iubd = i + 1;
                    }
                }
else if (tmp < ZERO) 
                {
                    if (slbd * tmp > upper_difference.get_entry(i) - trust_region_center_offset.get_entry(i)) 
                    {
                        slbd = (upper_difference.get_entry(i) - trust_region_center_offset.get_entry(i)) / tmp;
                        ilbd = i + 1;
                    }
                    if (subd * tmp < lower_difference.get_entry(i) - trust_region_center_offset.get_entry(i)) 
                    {
                        // Computing MAX
                        subd = std::max(sumin, (lower_difference.get_entry(i) - trust_region_center_offset.get_entry(i)) / tmp);
                        iubd = -i - 1;
                    }
                }
            }

            // Seek a large modulus of the KNEW-th Lagrange function when the index
            // of the other interpolation point on the line through XOPT is KNEW.

            step = slbd;
            int isbd = ilbd;
            double vlag;
            if (k == knew) 
            {
                const double diff = dderiv - ONE;
                vlag = slbd * (dderiv - slbd * diff);
                const double d1 = subd * (dderiv - subd * diff);
                if (std::abs(d1) > std::abs(vlag)) 
                {
                    step = subd;
                    vlag = d1;
                    isbd = iubd;
                }
                const double d2 = HALF * dderiv;
                const double d3 = d2 - diff * slbd;
                const double d4 = d2 - diff * subd;
                if (d3 * d4 < ZERO) 
                {
                    const double d5 = d2 * d2 / diff;
                    if (std::abs(d5) > std::abs(vlag)) 
                    {
                        step = d2 / diff;
                        vlag = d5;
                        isbd = 0;
                    }
                }

                // Search along each of the other lines through XOPT and another point.

            }
else 
            {
                vlag = slbd * (ONE - slbd);
                const double tmp = subd * (ONE - subd);
                if (std::abs(tmp) > std::abs(vlag)) 
                {
                    step = subd;
                    vlag = tmp;
                    isbd = iubd;
                }
                if (subd > HALF && std::abs(vlag) < ONE_OVER_FOUR) 
                {
                    step = HALF;
                    vlag = ONE_OVER_FOUR;
                    isbd = 0;
                }
                vlag *= dderiv;
            }

            // Calculate PREDSQ for the current line search and maintain PRESAV.

            const double tmp = step * (ONE - step) * distsq;
            const double predsq = vlag * vlag * (vlag * vlag + ha * tmp * tmp);
            if (predsq > presav) 
            {
                presav = predsq;
                ksav = k;
                stpsav = step;
                ibdsav = isbd;
            }
        }

        // Construct XNEW in a way that satisfies the bound constraints exactly.

        for (int i{}; i < n; i++) 
        {
            const double tmp = trust_region_center_offset.get_entry(i) + stpsav * (interpolation_points.get_entry(ksav, i) - trust_region_center_offset.get_entry(i));
            new_point.set_entry(i, std::max(lower_difference.get_entry(i), std::min(upper_difference.get_entry(i), tmp)));
        }
        if (ibdsav < 0) 
        {
            new_point.set_entry(-ibdsav - 1, lower_difference.get_entry(-ibdsav - 1));
        }
        if (ibdsav > 0) 
        {
            new_point.set_entry(ibdsav - 1, upper_difference.get_entry(ibdsav - 1));
        }

        // Prepare for the iterative method that assembles the constrained Cauchy
        // step in W. The sum of squares of the fixed components of W is formed in
        // WFIXSQ, and the free components of W are set to BIGSTP.

        const double bigstp = adelt + adelt;
        int iflag = 0;
        double cauchy = std::numeric_limits<double>::quiet_NaN();
        double csave = ZERO;
        while (true) 
        {
            double wfixsq = ZERO;
            double ggfree = ZERO;
            for (int i{}; i < n; i++) 
            {
                const double glag_value = glag.get_entry(i);
                work1.set_entry(i, ZERO);
                if (std::min(trust_region_center_offset.get_entry(i) - lower_difference.get_entry(i), glag_value) > ZERO ||
                    std::max(trust_region_center_offset.get_entry(i) - upper_difference.get_entry(i), glag_value) < ZERO) 
                    {
                    work1.set_entry(i, bigstp);
                    // Computing 2nd power
                    ggfree += glag_value * glag_value;
                }
            }
            if (ggfree == ZERO) 
            {
                return std::vector<double> { alpha, ZERO };
            }

            // Investigate whether more components of W can be fixed.
            const double tmp1 = adelt * adelt - wfixsq;
            if (tmp1 > ZERO) 
            {
                step = std::sqrt(tmp1 / ggfree);
                ggfree = ZERO;
                for (int i{}; i < n; i++) 
                {
                    if (work1.get_entry(i) == bigstp) 
                    {
                        const double tmp2 = trust_region_center_offset.get_entry(i) - step * glag.get_entry(i);
                        if (tmp2 <= lower_difference.get_entry(i)) 
                        {
                            work1.set_entry(i, lower_difference.get_entry(i) - trust_region_center_offset.get_entry(i));
                            // Computing 2nd power
                            const double d1 = work1.get_entry(i);
                            wfixsq += d1 * d1;
                        }
else if (tmp2 >= upper_difference.get_entry(i)) 
                        {
                            work1.set_entry(i, upper_difference.get_entry(i) - trust_region_center_offset.get_entry(i));
                            // Computing 2nd power
                            const double d1 = work1.get_entry(i);
                            wfixsq += d1 * d1;
                        }
else 
                        {
                            // Computing 2nd power
                            const double d1 = glag.get_entry(i);
                            ggfree += d1 * d1;
                        }
                    }
                }
            }

            // Set the remaining free components of W and all components of XALT, // except that W may be scaled later.

            double gw = ZERO;
            for (int i{}; i < n; i++) 
            {
                const double glag_value = glag.get_entry(i);
                if (work1.get_entry(i) == bigstp) 
                {
                    work1.set_entry(i, -step * glag_value);
                    const double min = std::min(upper_difference.get_entry(i), trust_region_center_offset.get_entry(i) + work1.get_entry(i));
                    alternative_new_point.set_entry(i, std::max(lower_difference.get_entry(i), min));
                }
else if (work1.get_entry(i) == ZERO) 
                {
                    alternative_new_point.set_entry(i, trust_region_center_offset.get_entry(i));
                }
else if (glag_value > ZERO) 
                {
                    alternative_new_point.set_entry(i, lower_difference.get_entry(i));
                }
else 
                {
                    alternative_new_point.set_entry(i, upper_difference.get_entry(i));
                }
                gw += glag_value * work1.get_entry(i);
            }

            // Set CURV to the curvature of the KNEW-th Lagrange function along W.
            // Scale W by a factor less than one if that can reduce the modulus of
            // the Lagrange function at XOPT+W. Set CAUCHY to the const value of
            // the square of this function.

            double curv = ZERO;
            for (int k{}; k < npt; k++) 
            {
                double tmp = ZERO;
                for (int j{}; j < n; j++) 
                {
                    tmp += interpolation_points.get_entry(k, j) * work1.get_entry(j);
                }
                curv += hcol.get_entry(k) * tmp * tmp;
            }
            if (iflag == 1) 
            {
                curv = -curv;
            }
            if (curv > -gw &&
                curv < -gw * (ONE + std::sqrt(TWO))) 
                {
                const double scale = -gw / curv;
                for (int i{}; i < n; i++) 
                {
                    const double tmp = trust_region_center_offset.get_entry(i) + scale * work1.get_entry(i);
                    alternative_new_point.set_entry(i, std::max(lower_difference.get_entry(i), std::min(upper_difference.get_entry(i), tmp)));
                }
                // Computing 2nd power
                const double d1 = HALF * gw * scale;
                cauchy = d1 * d1;
            }
else 
            {
                // Computing 2nd power
                const double d1 = gw + HALF * curv;
                cauchy = d1 * d1;
            }

            // If IFLAG is zero, then XALT is calculated as before after reversing
            // the sign of GLAG. Thus two XALT vectors become available. The one that
            // is chosen is the one that gives the larger value of CAUCHY.

            if (iflag == 0) 
            {
                for (int i{}; i < n; i++) 
                {
                    glag.set_entry(i, -glag.get_entry(i));
                    work2.set_entry(i, alternative_new_point.get_entry(i));
                }
                csave = cauchy;
                iflag = 1;
            }
else 
            {
                break;
            }
        }
        if (csave > cauchy) 
        {
            for (int i{}; i < n; i++) 
            {
                alternative_new_point.set_entry(i, work2.get_entry(i));
            }
            cauchy = csave;
        }

        return std::vector<double> { alpha, cauchy };
    } // altmov

    // ----------------------------------------------------------------------------------------

    /**
     *     SUBROUTINE PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ, *     BMAT and ZMAT for the first iteration, and it maintains the values of
     *     NF and KOPT. The vector X is also changed by PRELIM.
     *
     *     The arguments N, NPT, X, XL, XU, RHOBEG, IPRINT and MAXFUN are the
     *       same as the corresponding arguments in SUBROUTINE BOBYQA.
     *     The arguments XBASE, XPT, FVAL, HQ, PQ, BMAT, ZMAT, NDIM, SL and SU
     *       are the same as the corresponding arguments in BOBYQB, the elements
     *       of SL and SU being set in BOBYQA.
     *     GOPT is usually the gradient of the quadratic model at XOPT+XBASE, but
     *       it is set by PRELIM to the gradient of the quadratic model at XBASE.
     *       If XOPT is nonzero, BOBYQB will change it to its usual value later.
     *     NF is maintaned as the number of calls of CALFUN so far.
     *     KOPT will be such that the least calculated value of F so far is at
     *       the point XPT(KOPT,.)+XBASE in the space of the variables.
     *
     * @param lower_bound Lower bounds.
     * @param upper_bound Upper bounds.
     */
    private void prelim(std::vector<double> lower_bound, std::vector<double> upper_bound) 
    {

        const int n = current_best.get_dimension();
        const int& npt = number_of_interpolation_points;
        const int& ndim = b_matrix.get_row_dimension();

        const double rhosq = initial_trust_region_radius * initial_trust_region_radius;
        const double recip = 1.0/ rhosq;
        const int& np = n + 1;

        // Set XBASE to the initial vector of variables, and set the initial
        // elements of XPT, BMAT, HQ, PQ and ZMAT to zero.

        for (int j{}; j < n; j++) 
        {
            origin_shift.set_entry(j, current_best.get_entry(j));
            for (int k{}; k < npt; k++) 
            {
                interpolation_points.set_entry(k, j, ZERO);
            }
            for (int i{}; i < ndim; i++) 
            {
                b_matrix.set_entry(i, j, ZERO);
            }
        }
        const int max_i = n * np / 2;
        for (int i{}; i < max_i; i++) 
        {
            model_second_derivatives_values.set_entry(i, ZERO);
        }
        for (int k{}; k < npt; k++) 
        {
            model_second_derivatives_parameters.set_entry(k, ZERO);
            const int max_j = npt - np;
            for (int j{}; j < max_j; j++) 
            {
                z_matrix.set_entry(k, j, ZERO);
            }
        }

        // Begin the initialization procedure. NF becomes one more than the number
        // of function values so far. The coordinates of the displacement of the
        // next initial interpolation point from XBASE are set in XPT(NF+1,.).

        int ipt = 0;
        int jpt = 0;
        double fbeg = std::numeric_limits<double>::quiet_NaN();
        do 
        {
            const int& nfm = get_evaluations();
            const int& nfx = nfm - n;
            const int& nfmm = nfm - 1;
            const int& nfxm = nfx - 1;
            double stepa = 0;
            double stepb = 0;
            if (nfm <= 2 * n) 
            {
                if (nfm >= 1 &&
                    nfm <= n) 
                    {
                    stepa = initial_trust_region_radius;
                    if (upper_difference.get_entry(nfmm) == ZERO) 
                    {
                        stepa = -stepa;
                        // throw Path_isExploredException(); // XXX
                    }
                    interpolation_points.set_entry(nfm, nfmm, stepa);
                }
else if (nfm > n) 
                {
                    stepa = interpolation_points.get_entry(nfx, nfxm);
                    stepb = -initial_trust_region_radius;
                    if (lower_difference.get_entry(nfxm) == ZERO) 
                    {
                        stepb = std::min(TWO * initial_trust_region_radius, upper_difference.get_entry(nfxm));
                        // throw Path_isExploredException(); // XXX
                    }
                    if (upper_difference.get_entry(nfxm) == ZERO) 
                    {
                        stepb = std::max(-TWO * initial_trust_region_radius, lower_difference.get_entry(nfxm));
                        // throw Path_isExploredException(); // XXX
                    }
                    interpolation_points.set_entry(nfm, nfxm, stepb);
                }
            }
else 
            {
                const int tmp1 = (nfm - np) / n;
                jpt = nfm - tmp1 * n - n;
                ipt = jpt + tmp1;
                if (ipt > n) 
                {
                    const int tmp2 = jpt;
                    jpt = ipt - n;
                    ipt = tmp2;
//                     throw Path_isExploredException(); // XXX
                }
                const int ipt_minus1 = ipt - 1;
                const int jpt_minus1 = jpt - 1;
                interpolation_points.set_entry(nfm, ipt_minus1, interpolation_points.get_entry(ipt, ipt_minus1));
                interpolation_points.set_entry(nfm, jpt_minus1, interpolation_points.get_entry(jpt, jpt_minus1));
            }

            // Calculate the next value of F. The least function value so far and
            // its index are required.

            for (int j{}; j < n; j++) 
            {
                current_best.set_entry(j, std::min(std::max(lower_bound[j], origin_shift.get_entry(j) + interpolation_points.get_entry(nfm, j)), upper_bound[j]));
                if (interpolation_points.get_entry(nfm, j) == lower_difference.get_entry(j)) 
                {
                    current_best.set_entry(j, lower_bound[j]);
                }
                if (interpolation_points.get_entry(nfm, j) == upper_difference.get_entry(j)) 
                {
                    current_best.set_entry(j, upper_bound[j]);
                }
            }

            const double objective_value = compute_objective_value(current_best.to_array());
            const double f = is_minimize ? objective_value : -objective_value;
            const int& num_eval = get_evaluations(); // nfm + 1
            f_at_interpolation_points.set_entry(nfm, f);

            if (num_eval == 1) 
            {
                fbeg = f;
                trust_region_center_interpolation_point_index = 0;
            }
else if (f < f_at_interpolation_points.get_entry(trust_region_center_interpolation_point_index)) 
            {
                trust_region_center_interpolation_point_index = nfm;
            }

            // Set the nonzero initial elements of BMAT and the quadratic model in the
            // cases when NF is at most 2*N+1. If NF exceeds N+1, then the positions
            // of the NF-th and (NF-N)-th interpolation points may be switched, in
            // order that the function value at the first of them contributes to the
            // off-diagonal second derivative terms of the initial quadratic model.

            if (num_eval <= 2 * n + 1) 
            {
                if (num_eval >= 2 &&
                    num_eval <= n + 1) 
                    {
                    gradient_at_trust_region_center.set_entry(nfmm, (f - fbeg) / stepa);
                    if (npt < num_eval + n) 
                    {
                        const double one_over_step_a = ONE / stepa;
                        b_matrix.set_entry(0, nfmm, -one_over_step_a);
                        b_matrix.set_entry(nfm, nfmm, one_over_step_a);
                        b_matrix.set_entry(npt + nfmm, nfmm, -HALF * rhosq);
                        // throw Path_isExploredException(); // XXX
                    }
                }
else if (num_eval >= n + 2) 
                {
                    const int ih = nfx * (nfx + 1) / 2 - 1;
                    const double tmp = (f - fbeg) / stepb;
                    const double diff = stepb - stepa;
                    model_second_derivatives_values.set_entry(ih, TWO * (tmp - gradient_at_trust_region_center.get_entry(nfxm)) / diff);
                    gradient_at_trust_region_center.set_entry(nfxm, (gradient_at_trust_region_center.get_entry(nfxm) * stepb - tmp * stepa) / diff);
                    if (stepa * stepb < ZERO && f < f_at_interpolation_points.get_entry(nfm - n)) 
                    {
                        f_at_interpolation_points.set_entry(nfm, f_at_interpolation_points.get_entry(nfm - n));
                        f_at_interpolation_points.set_entry(nfm - n, f);
                        if (trust_region_center_interpolation_point_index == nfm) 
                        {
                            trust_region_center_interpolation_point_index = nfm - n;
                        }
                        interpolation_points.set_entry(nfm - n, nfxm, stepb);
                        interpolation_points.set_entry(nfm, nfxm, stepa);
                    }
                    b_matrix.set_entry(0, nfxm, -(stepa + stepb) / (stepa * stepb));
                    b_matrix.set_entry(nfm, nfxm, -HALF / interpolation_points.get_entry(nfm - n, nfxm));
                    b_matrix.set_entry(nfm - n, nfxm, -b_matrix.get_entry(0, nfxm) - b_matrix.get_entry(nfm, nfxm));
                    z_matrix.set_entry(0, nfxm, std::sqrt(TWO) / (stepa * stepb));
                    z_matrix.set_entry(nfm, nfxm, std::sqrt(HALF) / rhosq);
                    // z_matrix.set_entry(nfm, nfxm, std::sqrt(HALF) * recip); // XXX "test_ackley" and "test_diff_pow" fail.
                    z_matrix.set_entry(nfm - n, nfxm, -z_matrix.get_entry(0, nfxm) - z_matrix.get_entry(nfm, nfxm));
                }

                // Set the off-diagonal second derivatives of the Lagrange functions and
                // the initial quadratic model.

            }
else 
            {
                z_matrix.set_entry(0, nfxm, recip);
                z_matrix.set_entry(nfm, nfxm, recip);
                z_matrix.set_entry(ipt, nfxm, -recip);
                z_matrix.set_entry(jpt, nfxm, -recip);

                const int ih = ipt * (ipt - 1) / 2 + jpt - 1;
                const double tmp = interpolation_points.get_entry(nfm, ipt - 1) * interpolation_points.get_entry(nfm, jpt - 1);
                model_second_derivatives_values.set_entry(ih, (fbeg - f_at_interpolation_points.get_entry(ipt) - f_at_interpolation_points.get_entry(jpt) + f) / tmp);
//                 throw Path_isExploredException(); // XXX
            }
        } while (get_evaluations() < npt);
    } // prelim


    // ----------------------------------------------------------------------------------------

    /**
     *     A version of the truncated conjugate gradient is applied. If a line
     *     search is restricted by a constraint, then the procedure is restarted, *     the values of the variables that are at their bounds being fixed. If
     *     the trust region boundary is reached, then further changes may be made
     *     to D, each one being in the two dimensional space that is spanned
     *     by the current D and the gradient of Q at XOPT+D, staying on the trust
     *     region boundary. Termination occurs when the reduction in Q seems to
     *     be close to the greatest reduction that can be achieved.
     *     The arguments N, NPT, XPT, XOPT, GOPT, HQ, PQ, SL and SU have the same
     *       meanings as the corresponding arguments of BOBYQB.
     *     DELTA is the trust region radius for the present calculation, which
     *       seeks a small value of the quadratic model within distance DELTA of
     *       XOPT subject to the bounds on the variables.
     *     XNEW will be set to a vector of variables that is approximately
     *       the one that minimizes the quadratic model within the trust region
     *       subject to the SL and SU constraints on the variables. It satisfies
     *       as equations the bounds that become active during the calculation.
     *     D is the calculated trial step from XOPT, generated iteratively from an
     *       initial value of zero. Thus XNEW is XOPT+D after the const iteration.
     *     GNEW holds the gradient of the quadratic model at XOPT+D. It is updated
     *       when D is updated.
     *     xbdi.get( is a working space vector. For I=1,2,...,N, the element xbdi.get((I) is
     *       set to -1.0, 0.0, or 1.0, the value being nonzero if and only if the
     *       I-th variable has become fixed at a bound, the bound being SL(I) or
     *       SU(I) in the case xbdi.get((I)=-1.0 or xbdi.get((I)=1.0, respectively. This
     *       information is accumulated during the construction of XNEW.
     *     The arrays S, HS and HRED are also used for working space. They hold the
     *       current search direction, and the changes in the gradient of Q along S
     *       and the reduced D, respectively, where the reduced D is the same as D, *       except that the components of the fixed variables are zero.
     *     DSQ will be set to the square of the length of XNEW-XOPT.
     *     CRVMIN is set to zero if D reaches the trust region boundary. Otherwise
     *       it is set to the least curvature of H that occurs in the conjugate
     *       gradient searches that are not restricted by any constraints. The
     *       value CRVMIN=-1.0D0 is set, however, if all of these searches are
     *       constrained.
     * @param delta
     * @param gnew
     * @param xbdi
     * @param s
     * @param hs
     * @param hred
     */
    private std::vector<double> trsbox(
            double delta, Array_Real_Vector gnew, Array_Real_Vector xbdi, Array_Real_Vector s, Array_Real_Vector hs, Array_Real_Vector hred) 
            {

        const int n = current_best.get_dimension();
        const int& npt = number_of_interpolation_points;

        double dsq;
        double crvmin;

        // Local variables
        double dhd;
        double dhs;
        double cth;
        double shs;
        double sth;
        double ssq;
        double beta=0;
        double sdec;
        double blen;
        int iact = -1;
        double angt = 0;
        double qred;
        int isav;
        double temp;
        double xsav = 0;
        double xsum;
        double angbd = 0;
        double dredg = 0;
        double sredg = 0;
        double resid;
        double delsq;
        double ggsav = 0;
        double tempa;
        double tempb;
        double redmax;
        double dredsq = 0;
        double redsav;
        double gredsq = 0;
        double rednew;
        int itcsav = 0;
        double rdprev = 0;
        double rdnext = 0;
        double stplen;
        double stepsq = 0;
        int itermax = 0;

        // Set some constants.

        // Function Body

        // The sign of GOPT(I) gives the sign of the change to the I-th variable
        // that will reduce Q from its value at XOPT. Thus xbdi.get((I) shows whether
        // or not to fix the I-th variable at one of its bounds initially, with
        // NACT being set to the number of fixed variables. D and GNEW are also
        // set for the first iteration. DELSQ is the upper bound on the sum of
        // squares of the free variables. QRED is the reduction in Q so far.

        int iterc = 0;
        int nact = 0;
        for (int i{}; i < n; i++) 
        {
            xbdi.set_entry(i, ZERO);
            if (trust_region_center_offset.get_entry(i) <= lower_difference.get_entry(i)) 
            {
                if (gradient_at_trust_region_center.get_entry(i) >= ZERO) 
                {
                    xbdi.set_entry(i, MINUS_ONE);
                }
            }
else if (trust_region_center_offset.get_entry(i) >= upper_difference.get_entry(i) &&
                    gradient_at_trust_region_center.get_entry(i) <= ZERO) 
                    {
                xbdi.set_entry(i, ONE);
            }
            if (xbdi.get_entry(i) != ZERO) 
            {
                ++nact;
            }
            trial_step_point.set_entry(i, ZERO);
            gnew.set_entry(i, gradient_at_trust_region_center.get_entry(i));
        }
        delsq = delta * delta;
        qred = ZERO;
        crvmin = MINUS_ONE;

        // Set the next search direction of the conjugate gradient method. It is
        // the steepest descent direction initially and when the iterations are
        // restarted because a variable has just been fixed by a bound, and of
        // course the components of the fixed variables are zero. ITERMAX is an
        // upper bound on the indices of the conjugate gradient iterations.

        int state = 20;
        for(;;) 
        {
            switch (state) { // NOPMD - the reference algorithm is as complex as this, we simply ported it from Fortran with minimal changes
        case 20: 
        {
            beta = ZERO;
        }
        case 30: 
        {
            stepsq = ZERO;
            for (int i{}; i < n; i++) 
            {
                if (xbdi.get_entry(i) != ZERO) 
                {
                    s.set_entry(i, ZERO);
                }
else if (beta == ZERO) 
                {
                    s.set_entry(i, -gnew.get_entry(i));
                }
else 
                {
                    s.set_entry(i, beta * s.get_entry(i) - gnew.get_entry(i));
                }
                // Computing 2nd power
                const double d1 = s.get_entry(i);
                stepsq += d1 * d1;
            }
            if (stepsq == ZERO) 
            {
                state = 190; break;
            }
            if (beta == ZERO) 
            {
                gredsq = stepsq;
                itermax = iterc + n - nact;
            }
            if (gredsq * delsq <= qred * 1e-4 * qred) 
            {
                state = 190; break;
            }

            // Multiply the search direction by the second derivative matrix of Q and
            // calculate some scalars for the choice of steplength. Then set BLEN to
            // the length of the the step to the trust region boundary and STPLEN to
            // the steplength, ignoring the simple bounds.

            state = 210; break;
        }
        case 50: 
        {
            resid = delsq;
            double ds = ZERO;
            shs = ZERO;
            for (int i{}; i < n; i++) 
            {
                if (xbdi.get_entry(i) == ZERO) 
                {
                    // Computing 2nd power
                    const double d1 = trial_step_point.get_entry(i);
                    resid -= d1 * d1;
                    ds += s.get_entry(i) * trial_step_point.get_entry(i);
                    shs += s.get_entry(i) * hs.get_entry(i);
                }
            }
            if (resid <= ZERO) 
            {
                state = 90; break;
            }
            temp = std::sqrt(stepsq * resid + ds * ds);
            if (ds < ZERO) 
            {
                blen = (temp - ds) / stepsq;
            }
else 
            {
                blen = resid / (temp + ds);
            }
            stplen = blen;
            if (shs > ZERO) 
            {
                // Computing MIN
                stplen = std::min(blen, gredsq / shs);
            }

            // Reduce STPLEN if necessary in order to preserve the simple bounds, // letting IACT be the index of the constrained variable.

            iact = -1;
            for (int i{}; i < n; i++) 
            {
                if (s.get_entry(i) != ZERO) 
                {
                    xsum = trust_region_center_offset.get_entry(i) + trial_step_point.get_entry(i);
                    if (s.get_entry(i) > ZERO) 
                    {
                        temp = (upper_difference.get_entry(i) - xsum) / s.get_entry(i);
                    }
else 
                    {
                        temp = (lower_difference.get_entry(i) - xsum) / s.get_entry(i);
                    }
                    if (temp < stplen) 
                    {
                        stplen = temp;
                        iact = i;
                    }
                }
            }

            // Update CRVMIN, GNEW and D. Set SDEC to the decrease that occurs in Q.

            sdec = ZERO;
            if (stplen > ZERO) 
            {
                ++iterc;
                temp = shs / stepsq;
                if (iact == -1 && temp > ZERO) 
                {
                    crvmin = std::min(crvmin,temp);
                    if (crvmin == MINUS_ONE) 
                    {
                        crvmin = temp;
                    }
                }
                ggsav = gredsq;
                gredsq = ZERO;
                for (int i{}; i < n; i++) 
                {
                    gnew.set_entry(i, gnew.get_entry(i) + stplen * hs.get_entry(i));
                    if (xbdi.get_entry(i) == ZERO) 
                    {
                        // Computing 2nd power
                        const double d1 = gnew.get_entry(i);
                        gredsq += d1 * d1;
                    }
                    trial_step_point.set_entry(i, trial_step_point.get_entry(i) + stplen * s.get_entry(i));
                }
                // Computing MAX
                const double d1 = stplen * (ggsav - HALF * stplen * shs);
                sdec = std::max(d1, ZERO);
                qred += sdec;
            }

            // Restart the conjugate gradient method if it has hit a bound.

            if (iact >= 0) 
            {
                ++nact;
                xbdi.set_entry(iact, ONE);
                if (s.get_entry(iact) < ZERO) 
                {
                    xbdi.set_entry(iact, MINUS_ONE);
                }
                // Computing 2nd power
                const double d1 = trial_step_point.get_entry(iact);
                delsq -= d1 * d1;
                if (delsq <= ZERO) 
                {
                    state = 190; break;
                }
                state = 20; break;
            }

            // If STPLEN is less than BLEN, then either apply another conjugate
            // gradient iteration or RETURN.

            if (stplen < blen) 
            {
                if (iterc == itermax) 
                {
                    state = 190; break;
                }
                if (sdec <= qred * .01) 
                {
                    state = 190; break;
                }
                beta = gredsq / ggsav;
                state = 30; break;
            }
        }
        case 90: 
        {
            crvmin = ZERO;

            // Prepare for the alternative iteration by calculating some scalars
            // and by multiplying the reduced D by the second derivative matrix of
            // Q, where S holds the reduced D in the call of GGMULT.

        }
        case 100: 
        {
            if (nact >= n - 1) 
            {
                state = 190; break;
            }
            dredsq = ZERO;
            dredg = ZERO;
            gredsq = ZERO;
            for (int i{}; i < n; i++) 
            {
                if (xbdi.get_entry(i) == ZERO) 
                {
                    // Computing 2nd power
                    double d1 = trial_step_point.get_entry(i);
                    dredsq += d1 * d1;
                    dredg += trial_step_point.get_entry(i) * gnew.get_entry(i);
                    // Computing 2nd power
                    d1 = gnew.get_entry(i);
                    gredsq += d1 * d1;
                    s.set_entry(i, trial_step_point.get_entry(i));
                }
else 
                {
                    s.set_entry(i, ZERO);
                }
            }
            itcsav = iterc;
            state = 210; break;
            // Let the search direction S be a linear combination of the reduced D
            // and the reduced G that is orthogonal to the reduced D.
        }
        case 120: 
        {
            ++iterc;
            temp = gredsq * dredsq - dredg * dredg;
            if (temp <= qred * 1e-4 * qred) 
            {
                state = 190; break;
            }
            temp = std::sqrt(temp);
            for (int i{}; i < n; i++) 
            {
                if (xbdi.get_entry(i) == ZERO) 
                {
                    s.set_entry(i, (dredg * trial_step_point.get_entry(i) - dredsq * gnew.get_entry(i)) / temp);
                }
else 
                {
                    s.set_entry(i, ZERO);
                }
            }
            sredg = -temp;

            // By considering the simple bounds on the variables, calculate an upper
            // bound on the tangent of half the angle of the alternative iteration, // namely ANGBD, except that, if already a free variable has reached a
            // bound, there is a branch back to label 100 after fixing that variable.

            angbd = ONE;
            iact = -1;
            for (int i{}; i < n; i++) 
            {
                if (xbdi.get_entry(i) == ZERO) 
                {
                    tempa = trust_region_center_offset.get_entry(i) + trial_step_point.get_entry(i) - lower_difference.get_entry(i);
                    tempb = upper_difference.get_entry(i) - trust_region_center_offset.get_entry(i) - trial_step_point.get_entry(i);
                    if (tempa <= ZERO) 
                    {
                        ++nact;
                        xbdi.set_entry(i, MINUS_ONE);
                        break;
                    }
else if (tempb <= ZERO) 
                    {
                        ++nact;
                        xbdi.set_entry(i, ONE);
                        break;
                    }
                    // Computing 2nd power
                    double d1 = trial_step_point.get_entry(i);
                    // Computing 2nd power
                    double d2 = s.get_entry(i);
                    ssq = d1 * d1 + d2 * d2;
                    // Computing 2nd power
                    d1 = trust_region_center_offset.get_entry(i) - lower_difference.get_entry(i);
                    temp = ssq - d1 * d1;
                    if (temp > ZERO) 
                    {
                        temp = std::sqrt(temp) - s.get_entry(i);
                        if (angbd * temp > tempa) 
                        {
                            angbd = tempa / temp;
                            iact = i;
                            xsav = MINUS_ONE;
                        }
                    }
                    // Computing 2nd power
                    d1 = upper_difference.get_entry(i) - trust_region_center_offset.get_entry(i);
                    temp = ssq - d1 * d1;
                    if (temp > ZERO) 
                    {
                        temp = std::sqrt(temp) + s.get_entry(i);
                        if (angbd * temp > tempb) 
                        {
                            angbd = tempb / temp;
                            iact = i;
                            xsav = ONE;
                        }
                    }
                }
            }

            // Calculate HHD and some curvatures for the alternative iteration.

            state = 210; break;
        }
        case 150: 
        {
            shs = ZERO;
            dhs = ZERO;
            dhd = ZERO;
            for (int i{}; i < n; i++) 
            {
                if (xbdi.get_entry(i) == ZERO) 
                {
                    shs += s.get_entry(i) * hs.get_entry(i);
                    dhs += trial_step_point.get_entry(i) * hs.get_entry(i);
                    dhd += trial_step_point.get_entry(i) * hred.get_entry(i);
                }
            }

            // Seek the greatest reduction in Q for a range of equally spaced values
            // of ANGT in [0,ANGBD], where ANGT is the tangent of half the angle of
            // the alternative iteration.

            redmax = ZERO;
            isav = -1;
            redsav = ZERO;
            int iu = static_cast<int>( (angbd * 17. + 3.1);
            for (int i{}; i < iu; i++) 
            {
                angt = angbd * i / iu;
                sth = (angt + angt) / (ONE + angt * angt);
                temp = shs + angt * (angt * dhd - dhs - dhs);
                rednew = sth * (angt * dredg - sredg - HALF * sth * temp);
                if (rednew > redmax) 
                {
                    redmax = rednew;
                    isav = i;
                    rdprev = redsav;
                }
else if (i == isav + 1) 
                {
                    rdnext = rednew;
                }
                redsav = rednew;
            }

            // Return if the reduction is zero. Otherwise, set the sine and cosine
            // of the angle of the alternative iteration, and calculate SDEC.

            if (isav < 0) 
            {
                state = 190; break;
            }
            if (isav < iu) 
            {
                temp = (rdnext - rdprev) / (redmax + redmax - rdprev - rdnext);
                angt = angbd * (isav + HALF * temp) / iu;
            }
            cth = (ONE - angt * angt) / (ONE + angt * angt);
            sth = (angt + angt) / (ONE + angt * angt);
            temp = shs + angt * (angt * dhd - dhs - dhs);
            sdec = sth * (angt * dredg - sredg - HALF * sth * temp);
            if (sdec <= ZERO) 
            {
                state = 190; break;
            }

            // Update GNEW, D and HRED. If the angle of the alternative iteration
            // is restricted by a bound on a free variable, that variable is fixed
            // at the bound.

            dredg = ZERO;
            gredsq = ZERO;
            for (int i{}; i < n; i++) 
            {
                gnew.set_entry(i, gnew.get_entry(i) + (cth - ONE) * hred.get_entry(i) + sth * hs.get_entry(i));
                if (xbdi.get_entry(i) == ZERO) 
                {
                    trial_step_point.set_entry(i, cth * trial_step_point.get_entry(i) + sth * s.get_entry(i));
                    dredg += trial_step_point.get_entry(i) * gnew.get_entry(i);
                    // Computing 2nd power
                    const double d1 = gnew.get_entry(i);
                    gredsq += d1 * d1;
                }
                hred.set_entry(i, cth * hred.get_entry(i) + sth * hs.get_entry(i));
            }
            qred += sdec;
            if (iact >= 0 && isav == iu) 
            {
                ++nact;
                xbdi.set_entry(iact, xsav);
                state = 100; break;
            }

            // If SDEC is sufficiently small, then RETURN after setting XNEW to
            // XOPT+D, giving careful attention to the bounds.

            if (sdec > qred * .01) 
            {
                state = 120; break;
            }
        }
        case 190: 
        {
            dsq = ZERO;
            for (int i{}; i < n; i++) 
            {
                // Computing MAX
                // Computing MIN
                const double min = std::min(trust_region_center_offset.get_entry(i) + trial_step_point.get_entry(i), upper_difference.get_entry(i));
                new_point.set_entry(i, std::max(min, lower_difference.get_entry(i)));
                if (xbdi.get_entry(i) == MINUS_ONE) 
                {
                    new_point.set_entry(i, lower_difference.get_entry(i));
                }
                if (xbdi.get_entry(i) == ONE) 
                {
                    new_point.set_entry(i, upper_difference.get_entry(i));
                }
                trial_step_point.set_entry(i, new_point.get_entry(i) - trust_region_center_offset.get_entry(i));
                // Computing 2nd power
                const double d1 = trial_step_point.get_entry(i);
                dsq += d1 * d1;
            }
            return std::vector<double> { dsq, crvmin };
            // The following instructions multiply the current S-vector by the second
            // derivative matrix of the quadratic model, putting the product in HS.
            // They are reached from three different parts of the software above and
            // they can be regarded as an external subroutine.
        }
        case 210: 
        {
            int ih = 0;
            for (int j{}; j < n; j++) 
            {
                hs.set_entry(j, ZERO);
                for (int i{}; i <= j; i++) 
                {
                    if (i < j) 
                    {
                        hs.set_entry(j, hs.get_entry(j) + model_second_derivatives_values.get_entry(ih) * s.get_entry(i));
                    }
                    hs.set_entry(i, hs.get_entry(i) + model_second_derivatives_values.get_entry(ih) * s.get_entry(j));
                    ih++;
                }
            }
            const Real_Vector tmp = interpolation_points.operate(s).ebe_multiply(model_second_derivatives_parameters);
            for (int k{}; k < npt; k++) 
            {
                if (model_second_derivatives_parameters.get_entry(k) != ZERO) 
                {
                    for (int i{}; i < n; i++) 
                    {
                        hs.set_entry(i, hs.get_entry(i) + tmp.get_entry(k) * interpolation_points.get_entry(k, i));
                    }
                }
            }
            if (crvmin != ZERO) 
            {
                state = 50; break;
            }
            if (iterc > itcsav) 
            {
                state = 150; break;
            }
            for (int i{}; i < n; i++) 
            {
                hred.set_entry(i, hs.get_entry(i));
            }
            state = 120; break;
        }
        default: 
        {
            throw Math_Illegal_State_Exception(Localized_Core_Formats.SIMPLE_MESSAGE, "trsbox");
        }}
        }
    } // trsbox

    // ----------------------------------------------------------------------------------------

    /**
     *     The arrays BMAT and ZMAT are updated, as required by the position
     *     of the interpolation point that has the index KNEW. The vector VLAG has
     *     N+NPT components, set on entry to the first NPT and last N components
     *     of the product Hw in equation (4.11) of the Powell (2006) paper on
     *     NEWUOA. Further, BETA is set on entry to the value of the parameter
     *     with that name, and DENOM is set to the denominator of the updating
     *     formula. Elements of ZMAT may be treated as zero if their moduli are
     *     at most ZTEST. The first NDIM elements of W are used for working space.
     * @param beta
     * @param denom
     * @param knew
     */
    private void update(
            double beta, double denom, const int& knew) 
            {

        const int n = current_best.get_dimension();
        const int& npt = number_of_interpolation_points;
        const int& nptm = npt - n - 1;

        // XXX Should probably be split into two arrays.
        const Array_Real_Vector work = Array_Real_Vector(npt + n);

        double ztest = ZERO;
        for (int k{}; k < npt; k++) 
        {
            for (int j{}; j < nptm; j++) 
            {
                // Computing MAX
                ztest = std::max(ztest, std::abs(z_matrix.get_entry(k, j)));
            }
        }
        ztest *= 1e-20;

        // Apply the rotations that put zeros in the KNEW-th row of ZMAT.

        for (int j{ 1 }; j < nptm; j++) 
        {
            const double d1 = z_matrix.get_entry(knew, j);
            if (std::abs(d1) > ztest) 
            {
                // Computing 2nd power
                const double d2 = z_matrix.get_entry(knew, 0);
                // Computing 2nd power
                const double d3 = z_matrix.get_entry(knew, j);
                const double d4 = std::sqrt(d2 * d2 + d3 * d3);
                const double d5 = z_matrix.get_entry(knew, 0) / d4;
                const double d6 = z_matrix.get_entry(knew, j) / d4;
                for (int i{}; i < npt; i++) 
                {
                    const double d7 = d5 * z_matrix.get_entry(i, 0) + d6 * z_matrix.get_entry(i, j);
                    z_matrix.set_entry(i, j, d5 * z_matrix.get_entry(i, j) - d6 * z_matrix.get_entry(i, 0));
                    z_matrix.set_entry(i, 0, d7);
                }
            }
            z_matrix.set_entry(knew, j, ZERO);
        }

        // Put the first NPT components of the KNEW-th column of HLAG into W, // and calculate the parameters of the updating formula.

        for (int i{}; i < npt; i++) 
        {
            work.set_entry(i, z_matrix.get_entry(knew, 0) * z_matrix.get_entry(i, 0));
        }
        const double& alpha = work.get_entry(knew);
        const double tau = lagrange_values_at_new_point.get_entry(knew);
        lagrange_values_at_new_point.set_entry(knew, lagrange_values_at_new_point.get_entry(knew) - ONE);

        // Complete the updating of ZMAT.

        const double sqrt_denom = std::sqrt(denom);
        const double d1 = tau / sqrt_denom;
        const double d2 = z_matrix.get_entry(knew, 0) / sqrt_denom;
        for (int i{}; i < npt; i++) 
        {
            z_matrix.set_entry(i, 0, d1 * z_matrix.get_entry(i, 0) - d2 * lagrange_values_at_new_point.get_entry(i));
        }

        // Finally, update the matrix BMAT.

        for (int j{}; j < n; j++) 
        {
            const int jp = npt + j;
            work.set_entry(jp, b_matrix.get_entry(knew, j));
            const double d3 = (alpha * lagrange_values_at_new_point.get_entry(jp) - tau * work.get_entry(jp)) / denom;
            const double d4 = (-beta * work.get_entry(jp) - tau * lagrange_values_at_new_point.get_entry(jp)) / denom;
            for (int i{}; i <= jp; i++) 
            {
                b_matrix.set_entry(i, j, b_matrix.get_entry(i, j) + d3 * lagrange_values_at_new_point.get_entry(i) + d4 * work.get_entry(i));
                if (i >= npt) 
                {
                    b_matrix.set_entry(jp, (i - npt), b_matrix.get_entry(i, j));
                }
            }
        }
    } // update

    /**
     * Performs validity checks.
     *
     * @param lower_bound Lower bounds (constraints) of the objective variables.
     * @param upper_bound Upperer bounds (constraints) of the objective variables.
     */
    private void setup(std::vector<double> lower_bound, std::vector<double> upper_bound) 
    {

        std::vector<double> init = get_start_point();
        const int dimension = init.size();

        // Check problem dimension.
        if (dimension < MINIMUM_PROBLEM_DIMENSION) 
        {
            throw (Localized_Core_Formats.NUMBER_TOO_SMALL, dimension, MINIMUM_PROBLEM_DIMENSION);
        }
        // Check number of interpolation points.
        const std::vector<int> n_points_interval = { dimension + 2, (dimension + 2) * (dimension + 1) / 2 };
        if (number_of_interpolation_points < n_points_interval[0] ||
            number_of_interpolation_points > n_points_interval[1]) 
            {
            throw (Localized_Core_Formats.NUMBER_OF_INTERPOLATION_POINTS, number_of_interpolation_points, n_points_interval[0], n_points_interval[1]);
        }

        // Initialize bound differences.
        bound_difference = std::vector<double>(dimension];

        double required_min_diff = 2 * initial_trust_region_radius;
        double min_diff = INFINITY;
        for (int i{}; i < dimension; i++) 
        {
            bound_difference[i] = upper_bound[i] - lower_bound[i];
            min_diff = std::min(min_diff, bound_difference[i]);
        }
        if (min_diff < required_min_diff) 
        {
            initial_trust_region_radius = min_diff / 3.0;
        }

        // Initialize the data structures used by the "bobyqa" method.
        b_matrix = Array_2D_Row_Real_Matrix(dimension + number_of_interpolation_points, dimension);
        z_matrix = Array_2D_Row_Real_Matrix(number_of_interpolation_points, number_of_interpolation_points - dimension - 1);
        interpolation_points = Array_2D_Row_Real_Matrix(number_of_interpolation_points, dimension);
        origin_shift = Array_Real_Vector(dimension);
        f_at_interpolation_points = Array_Real_Vector(number_of_interpolation_points);
        trust_region_center_offset = Array_Real_Vector(dimension);
        gradient_at_trust_region_center = Array_Real_Vector(dimension);
        lower_difference = Array_Real_Vector(dimension);
        upper_difference = Array_Real_Vector(dimension);
        model_second_derivatives_parameters = Array_Real_Vector(number_of_interpolation_points);
        new_point = Array_Real_Vector(dimension);
        alternative_new_point = Array_Real_Vector(dimension);
        trial_step_point = Array_Real_Vector(dimension);
        lagrange_values_at_new_point = Array_Real_Vector(dimension + number_of_interpolation_points);
        model_second_derivatives_values = Array_Real_Vector(dimension * (dimension + 1) / 2);
    }

}
//CHECKSTYLE: resume all


