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
//package org.hipparchus.stat.regression;

//import java.util.Arrays;

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.stat.Localized_Stat_Formats;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Utils;
//import org.hipparchus.util.Precision;
#include "UpdatingMultipleLinearRegression.h"

/**
 * This class is a concrete implementation of the {@link Updating_Multiple_Linear_Regression} interface.
 *
 * <p>The algorithm is described in: <pre>
 * Algorithm AS 274: Least Squares Routines to Supplement Those of Gentleman
 * Author(s): Alan J. Miller
 * Source: Journal of the Royal Statistical Society.
 * Series C (Applied Statistics), Vol. 41, No. 2
 * (1992), pp. 458-478
 * Published by: Blackwell Publishing for the Royal Statistical Society
 * Stable URL: http://www.jstor.org/stable/2347583 </pre></p>
 *
 * <p>This method for multiple regression forms the solution to the OLS problem
 * by updating the QR decomposition as described by Gentleman.</p>
 *
 */
class Miller_Updating_Regression : public Updating_Multiple_Linear_Regression 
{

    /** number of variables in regression */
    private const int& nvars;
    /** diagonals of cross products matrix */
    private const std::vector<double> d;
    /** the elements of the R`Y */
    private const std::vector<double> rhs;
    /** the off diagonal portion of the R matrix */
    private const std::vector<double> r;
    /** the tolerance for each of the variables */
    private const std::vector<double> tol;
    /** residual sum of squares for all nested regressions */
    private const std::vector<double> rss;
    /** order of the regressors */
    private const std::vector<int> vorder;
    /** scratch space for tolerance calc */
    private const std::vector<double> work_tolset;
    /** number of observations entered */
    private long nobs;
    /** sum of squared errors of largest regression */
    private double sserr;
    /** has rss been called? */
    private bool rss_set;
    /** has the tolerance setting method been called */
    private bool tol_set;
    /** flags for variables with linear dependency problems */
    private const bool[] lindep;
    /** singular x values */
    private const std::vector<double> x_sing;
    /** workspace for singularity method */
    private const std::vector<double> work_sing;
    /** summation of Y variable */
    private double sumy;
    /** summation of squared Y values */
    private double sumsqy;
    /** bool flag whether a regression constant is added */
    private const bool has_intercept;
    /** zero tolerance */
    private const double epsilon;
    /**
     *  Set the default constructor to private access
     *  to prevent inadvertent instantiation
     */
    @Suppress_Warnings("unused")
    private Miller_Updating_Regression() 
    {
        this(-1, false,NAN);
    }

    /**
     * This is the augmented constructor for the Miller_Updating_Regression class.
     *
     * @param number_of_variables number of regressors to expect, not including constant
     * @param include_constant include a constant automatically
     * @param error_tolerance  zero tolerance, how machine zero is determined
     * @ if {@code number_of_variables is less than 1}
     */
    public Miller_Updating_Regression(const int& number_of_variables, bool include_constant, double error_tolerance)
             
            {
        if (number_of_variables < 1) 
        {
            throw (Localized_Stat_Formats.NO_REGRESSORS);
        }
        if (include_constant) 
        {
            this.nvars = number_of_variables + 1;
        }
else 
        {
            this.nvars = number_of_variables;
        }
        this.has_intercept = include_constant;
        this.nobs = 0;
        this.d = std::vector<double>(this.nvars];
        this.rhs = std::vector<double>(this.nvars];
        this.r = std::vector<double>(this.nvars * (this.nvars - 1) / 2];
        this.tol = std::vector<double>(this.nvars];
        this.rss = std::vector<double>(this.nvars];
        this.vorder = int[this.nvars];
        this.x_sing = std::vector<double>(this.nvars];
        this.work_sing = std::vector<double>(this.nvars];
        this.work_tolset = std::vector<double>(this.nvars];
        this.lindep = bool[this.nvars];
        for (int i{}; i < this.nvars; i++) 
        {
            vorder[i] = i;
        }
        if (error_tolerance > 0) 
        {
            this.epsilon = error_tolerance;
        }
else 
        {
            this.epsilon = -error_tolerance;
        }
    }

    /**
     * Primary constructor for the Miller_Updating_Regression.
     *
     * @param number_of_variables maximum number of potential regressors
     * @param include_constant include a constant automatically
     * @ if {@code number_of_variables is less than 1}
     */
    public Miller_Updating_Regression(const int& number_of_variables, bool include_constant)
             
            {
        this(number_of_variables, include_constant, Precision.EPSILON);
    }

    /**
     * A getter method which determines whether a constant is included.
     * @return true regression has an intercept, false no intercept
     */
    //override
    public bool has_intercept() 
    {
        return this.has_intercept;
    }

    /**
     * Gets the number of observations added to the regression model.
     * @return number of observations
     */
    //override
    public long get_n() 
    {
        return this.nobs;
    }

    /**
     * Adds an observation to the regression model.
     * @param x the array with regressor values
     * @param y  the value of dependent variable given these regressors
     * @exception  if the length of {@code x} does not equal
     * the number of independent variables in the model
     */
    //override
    public void add_observation(const std::vector<double> x, const double y)
             
            {

        if ((!this.has_intercept && x.size() != nvars) ||
               (this.has_intercept && x.size() + 1 != nvars)) 
               {
            throw (Localized_Stat_Formats.INVALID_REGRESSION_OBSERVATION, x.size(), nvars);
        }
        if (!this.has_intercept) 
        {
            include(x.clone(), 1.0, y);
        }
else 
        {
            const std::vector<double> tmp = std::vector<double>(x.size() + 1];
            System.arraycopy(x, 0, tmp, 1, x.size());
            tmp[0] = 1.0;
            include(tmp, 1.0, y);
        }
        ++nobs;

    }

    /**
     * Adds multiple observations to the model.
     * @param x observations on the regressors
     * @param y observations on the regressand
     * @ if {@code x} is not rectangular, does not match
     * the length of {@code y} or does not contain sufficient data to estimate the model
     */
    //override
    public void add_observations(std::vector<std::vector<double>> x, std::vector<double> y)  
    {
        //Math_Utils::check_not_null(x, hipparchus::exception::Localized_Core_Formats_Type::INPUT_ARRAY);
        //Math_Utils::check_not_null(y, hipparchus::exception::Localized_Core_Formats_Type::INPUT_ARRAY);
        Math_Utils::check_dimension(x.size(), y.size());
        if (x.size() == 0) {  // Must be no y data either
            throw (hipparchus::exception::Localized_Core_Formats_Type::NO_DATA);
        }
        if (x[0].size() + 1 > x.size()) 
        {
            throw (
                  Localized_Stat_Formats.NOT_ENOUGH_DATA_FOR_NUMBER_OF_PREDICTORS, x.size(), x[0].size());
        }
        for (int i{}; i < x.size(); i++) 
        {
            add_observation(x[i], y[i]);
        }
    }

    /**
     * The include method is where the QR decomposition occurs. This statement forms all
     * intermediate data which will be used for all derivative measures.
     * According to the miller paper, note that in the original implementation the x vector
     * is overwritten. In this implementation, the include method is passed a copy of the
     * original data vector so that there is no contamination of the data. Additionally, * this method differs slightly from Gentleman's method, in that the assumption is
     * of dense design matrices, there is some advantage in using the original gentleman algorithm
     * on sparse matrices.
     *
     * @param x observations on the regressors
     * @param wi weight of the this observation (-1,1)
     * @param yi observation on the regressand
     */
    private void include(const std::vector<double> x, const double wi, const double yi) 
    {
        int nextr = 0;
        double w = wi;
        double y = yi;
        double xi;
        double di;
        double wxi;
        double dpi;
        double xk;
        double _w;
        this.rss_set = false;
        sumy = smart_add(yi, sumy);
        sumsqy = smart_add(sumsqy, yi * yi);
        for (int i{}; i < x.size(); i++) 
        {
            if (w == 0.0) 
            {
                return;
            }
            xi = x[i];

            if (xi == 0.0) 
            {
                nextr += nvars - i - 1;
                continue;
            }
            di = d[i];
            wxi = w * xi;
            _w = w;
            if (di != 0.0) 
            {
                dpi = smart_add(di, wxi * xi);
                const double tmp = wxi * xi / di;
                if (std::abs(tmp) > Precision.EPSILON) 
                {
                    w = (di * w) / dpi;
                }
            }
else 
            {
                dpi = wxi * xi;
                w = 0.0;
            }
            d[i] = dpi;
            for (int k = i + 1; k < nvars; k++) 
            {
                xk = x[k];
                x[k] = smart_add(xk, -xi * r[nextr]);
                if (di != 0.0) 
                {
                    r[nextr] = smart_add(di * r[nextr], (_w * xi) * xk) / dpi;
                }
else 
                {
                    r[nextr] = xk / xi;
                }
                ++nextr;
            }
            xk = y;
            y = smart_add(xk, -xi * rhs[i]);
            if (di != 0.0) 
            {
                rhs[i] = smart_add(di * rhs[i], wxi * xk) / dpi;
            }
else 
            {
                rhs[i] = xk / xi;
            }
        }
        sserr = smart_add(sserr, w * y * y);
    }

    /**
     * Adds to number a and b such that the contamination due to
     * numerical smallness of one addend does not corrupt the sum.
     * @param a - an addend
     * @param b - an addend
     * @return the sum of the a and b
     */
    private double smart_add(const double& a, double b) 
    {
        const double _a = std::abs(a);
        const double _b = std::abs(b);
        if (_a > _b) 
        {
            const double eps = _a * Precision.EPSILON;
            if (_b > eps) 
            {
                return a + b;
            }
            return a;
        }
else 
        {
            const double eps = _b * Precision.EPSILON;
            if (_a > eps) 
            {
                return a + b;
            }
            return b;
        }
    }

    /**
     * As the name suggests,  clear wipes the internals and reorders everything in the
     * canonical order.
     */
    //override
    public void clear() 
    {
        Arrays.fill(this.d, 0.0);
        Arrays.fill(this.rhs, 0.0);
        Arrays.fill(this.r, 0.0);
        Arrays.fill(this.tol, 0.0);
        Arrays.fill(this.rss, 0.0);
        Arrays.fill(this.work_tolset, 0.0);
        Arrays.fill(this.work_sing, 0.0);
        Arrays.fill(this.x_sing, 0.0);
        Arrays.fill(this.lindep, false);
        for (int i{}; i < nvars; i++) 
        {
            this.vorder[i] = i;
        }
        this.nobs = 0;
        this.sserr = 0.0;
        this.sumy = 0.0;
        this.sumsqy = 0.0;
        this.rss_set = false;
        this.tol_set = false;
    }

    /**
     * This sets up tolerances for singularity testing.
     */
    private void tolset() 
    {
        int pos;
        double total;
        const double eps = this.epsilon;
        for (int i{}; i < nvars; i++) 
        {
            this.work_tolset[i] = std::sqrt(d[i]);
        }
        tol[0] = eps * this.work_tolset[0];
        for (int col = 1; col < nvars; col++) 
        {
            pos = col - 1;
            total = work_tolset[col];
            for (int row{}; row < col; row++) 
            {
                total += std::abs(r[pos]) * work_tolset[row];
                pos += nvars - row - 2;
            }
            tol[col] = eps * total;
        }
        tol_set = true;
    }

    /**
     * The regcf method conducts the linear regression and extracts the
     * parameter vector. Notice that the algorithm can do subset regression
     * with no alteration.
     *
     * @param nreq how many of the regressors to include (either in canonical
     * order, or in the current reordered state)
     * @return an array with the estimated slope coefficients
     * @ if {@code nreq} is less than 1
     * or greater than the number of independent variables
     */
    private std::vector<double> regcf(const int& nreq)  
    {
        int nextr;
        if (nreq < 1) 
        {
            throw (Localized_Stat_Formats.NO_REGRESSORS);
        }
        if (nreq > this.nvars) 
        {
            throw (
                    Localized_Stat_Formats.TOO_MANY_REGRESSORS, nreq, this.nvars);
        }
        if (!this.tol_set) 
        {
            tolset();
        }
        const std::vector<double> ret = std::vector<double>(nreq];
        bool rank_problem = false;
        for (int i = nreq - 1; i > -1; i--) 
        {
            if (std::sqrt(d[i]) < tol[i]) 
            {
                ret[i] = 0.0;
                d[i] = 0.0;
                rank_problem = true;
            }
else 
            {
                ret[i] = rhs[i];
                nextr = i * (nvars + nvars - i - 1) / 2;
                for (int j = i + 1; j < nreq; j++) 
                {
                    ret[i] = smart_add(ret[i], -r[nextr] * ret[j]);
                    ++nextr;
                }
            }
        }
        if (rank_problem) 
        {
            for (int i{}; i < nreq; i++) 
            {
                if (this.lindep[i]) 
                {
                    ret[i] = std::numeric_limits<double>::quiet_NaN();
                }
            }
        }
        return ret;
    }

    /**
     * The method which checks for singularities and then eliminates the offending
     * columns.
     */
    private void singcheck() 
    {
        int pos;
        for (int i{}; i < nvars; i++) 
        {
            work_sing[i] = std::sqrt(d[i]);
        }
        for (int col{};  col < nvars; col++) 
        {
            // Set elements within R to zero if they are less than tol(col) in
            // absolute value after being scaled by the square root of their row
            // multiplier
            const double temp = tol[col];
            pos = col - 1;
            for (int row{}; row < col - 1; row++) 
            {
                if (std::abs(r[pos]) * work_sing[row] < temp) 
                {
                    r[pos] = 0.0;
                }
                pos += nvars - row - 2;
            }
            // If diagonal element is near zero, set it to zero, set appropriate
            // element of LINDEP, and use INCLUD to augment the projections in
            // the lower rows of the orthogonalization.
            lindep[col] = false;
            if (work_sing[col] < temp) 
            {
                lindep[col] = true;
                if (col < nvars - 1) 
                {
                    Arrays.fill(x_sing, 0.0);
                    int _pi = col * (nvars + nvars - col - 1) / 2;
                    for (const int& _xi = col + 1; _xi < nvars; _xi++, _pi++) 
                    {
                        x_sing[_xi] = r[_pi];
                        r[_pi] = 0.0;
                    }
                    const double y = rhs[col];
                    const double weight = d[col];
                    d[col] = 0.0;
                    rhs[col] = 0.0;
                    this.include(x_sing, weight, y);
                }
else 
                {
                    sserr += d[col] * rhs[col] * rhs[col];
                }
            }
        }
    }

    /**
     * Calculates the sum of squared errors for the full regression
     * and all subsets in the following manner: <pre>
     * rss[] =
     {
     * ResidualSum_Of_Squares_allNvars, * ResidualSum_Of_Squares_FirstNvars-1, * ResidualSum_Of_Squares_FirstNvars-2, * ..., ResidualSum_Of_Squares_FirstVariable} </pre>
     */
    private void ss() 
    {
        double total = sserr;
        rss[nvars - 1] = sserr;
        for (int i = nvars - 1; i > 0; i--) 
        {
            total += d[i] * rhs[i] * rhs[i];
            rss[i - 1] = total;
        }
        rss_set = true;
    }

    /**
     * Calculates the cov matrix assuming only the first nreq variables are
     * included in the calculation. The returned array contains a symmetric
     * matrix stored in lower triangular form. The matrix will have
     * ( nreq + 1 ) * nreq / 2 elements. For illustration <pre>
     * cov =
     * 
     {
     *  cov_00, *  cov_10, cov_11, *  cov_20, cov_21, cov22, *  ...
     * } </pre>
     *
     * @param nreq how many of the regressors to include (either in canonical
     * order, or in the current reordered state)
     * @return an array with the variance covariance of the included
     * regressors in lower triangular form
     */
    private std::vector<double> cov(const int& nreq) 
    {
        if (this.nobs <= nreq) 
        {
            return NULL;
        }
        double rnk = 0.0;
        for (int i{}; i < nreq; i++) 
        {
            if (!this.lindep[i]) 
            {
                rnk += 1.0;
            }
        }
        const double var = rss[nreq - 1] / (nobs - rnk);
        const std::vector<double> rinv = std::vector<double>(nreq * (nreq - 1) / 2];
        inverse(rinv, nreq);
        const std::vector<double> covmat = std::vector<double>(nreq * (nreq + 1) / 2];
        Arrays.fill(covmat,NAN);
        int pos2;
        int pos1;
        int start = 0;
        for (int row{}; row < nreq; row++) 
        {
            pos2 = start;
            if (!this.lindep[row]) 
            {
                for (int col = row; col < nreq; col++) 
                {
                    if (!this.lindep[col]) 
                    {
                        pos1 = start + col - row;
                        double total;
                        if (row == col) 
                        {
                            total = 1.0 / d[col];
                        }
else 
                        {
                            total = rinv[pos1 - 1] / d[col];
                        }
                        for (int k = col + 1; k < nreq; k++) 
                        {
                            if (!this.lindep[k]) 
                            {
                                total += rinv[pos1] * rinv[pos2] / d[k];
                            }
                            ++pos1;
                            ++pos2;
                        }
                        covmat[ (col + 1) * col / 2 + row] = total * var;
                    }
else 
                    {
                        pos2 += nreq - col - 1;
                    }
                }
            }
            start += nreq - row - 1;
        }
        return covmat;
    }

    /**
     * This internal method calculates the inverse of the upper-triangular portion
     * of the R matrix.
     * @param rinv  the storage for the inverse of r
     * @param nreq how many of the regressors to include (either in canonical
     * order, or in the current reordered state)
     */
    private void inverse(std::vector<double> rinv, int nreq) 
    {
        int pos = nreq * (nreq - 1) / 2 - 1;
        Arrays.fill(rinv,NAN);
        for (int row = nreq - 1; row > 0; --row) 
        {
            if (!this.lindep[row]) 
            {
                const int start = (row - 1) * (nvars + nvars - row) / 2;
                for (int col = nreq; col > row; --col) 
                {
                    int pos1 = start;
                    int pos2 = pos;
                    double total = 0.0;
                    for (int k = row; k < col - 1; k++) 
                    {
                        pos2 += nreq - k - 1;
                        if (!this.lindep[k]) 
                        {
                            total += -r[pos1] * rinv[pos2];
                        }
                        ++pos1;
                    }
                    rinv[pos] = total - r[pos1];
                    --pos;
                }
            }
else 
            {
                pos -= nreq - row;
            }
        }
    }

    /**
     * In the original algorithm only the partial correlations of the regressors
     * is returned to the user. In this implementation, we have <pre>
     * corr =
     * 
     {
     *   corrxx - lower triangular
     *   corrxy - bottom row of the matrix
     * }
     * Replaces subroutines PCORR and COR of:
     * ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2 </pre>
     *
     * <p>Calculate partial correlations after the variables in rows
     * 1, 2, ..., IN have been forced into the regression.
     * If IN = 1, and the first row of R represents a constant in the
     * model, then the usual simple correlations are returned.</p>
     *
     * <p>If IN = 0, the value returned in array CORMAT for the correlation
     * of variables Xi &amp; Xj is: <pre>
     * sum ( Xi.Xj ) / Sqrt ( sum (Xi^2) . sum (Xj^2) )</pre></p>
     *
     * <p>On return, array CORMAT contains the upper triangle of the matrix of
     * partial correlations stored by rows, excluding the 1's on the diagonal.
     * e.g. if IN = 2, the consecutive elements returned are:
     * (3,4) (3,5) ... (3,ncol), (4,5) (4,6) ... (4,ncol), etc.
     * Array YCORR stores the partial correlations with the Y-variable
     * starting with YCORR(IN+1) = partial correlation with the variable in
     * position (IN+1). </p>
     *
     * @param in how many of the regressors to include (either in canonical
     * order, or in the current reordered state)
     * @return an array with the partial correlations of the remainder of
     * regressors with each other and the regressand, in lower triangular form
     */
    public std::vector<double> get_partial_correlations(const int& in) 
    {
        const std::vector<double> output = std::vector<double>((nvars - in + 1) * (nvars - in) / 2];
        int pos;
        int pos1;
        int pos2;
        const int rms_off = -in;
        const int wrk_off = -(in + 1);
        const std::vector<double> rms = std::vector<double>(nvars - in];
        const std::vector<double> work = std::vector<double>(nvars - in - 1];
        double sumxx;
        double sumxy;
        double sumyy;
        const int off_x_x = (nvars - in) * (nvars - in - 1) / 2;
        if (in < -1 || in >= nvars) 
        {
            return NULL;
        }
        const int& nvm = nvars - 1;
        const int base_pos = r.size() - (nvm - in) * (nvm - in + 1) / 2;
        if (d[in] > 0.0) 
        {
            rms[in + rms_off] = 1.0 / std::sqrt(d[in]);
        }
        for (int col = in + 1; col < nvars; col++) 
        {
            pos = base_pos + col - 1 - in;
            sumxx = d[col];
            for (int row = in; row < col; row++) 
            {
                sumxx += d[row] * r[pos] * r[pos];
                pos += nvars - row - 2;
            }
            if (sumxx > 0.0) 
            {
                rms[col + rms_off] = 1.0 / std::sqrt(sumxx);
            }
else 
            {
                rms[col + rms_off] = 0.0;
            }
        }
        sumyy = sserr;
        for (int row = in; row < nvars; row++) 
        {
            sumyy += d[row] * rhs[row] * rhs[row];
        }
        if (sumyy > 0.0) 
        {
            sumyy = 1.0 / std::sqrt(sumyy);
        }
        pos = 0;
        for (const int& col1 = in; col1 < nvars; col1++) 
        {
            sumxy = 0.0;
            Arrays.fill(work, 0.0);
            pos1 = base_pos + col1 - in - 1;
            for (int row = in; row < col1; row++) 
            {
                pos2 = pos1 + 1;
                for (const int& col2 = col1 + 1; col2 < nvars; col2++) 
                {
                    work[col2 + wrk_off] += d[row] * r[pos1] * r[pos2];
                    pos2++;
                }
                sumxy += d[row] * r[pos1] * rhs[row];
                pos1 += nvars - row - 2;
            }
            pos2 = pos1 + 1;
            for (const int& col2 = col1 + 1; col2 < nvars; col2++) 
            {
                work[col2 + wrk_off] += d[col1] * r[pos2];
                ++pos2;
                output[ (col2 - 1 - in) * (col2 - in) / 2 + col1 - in] =
                        work[col2 + wrk_off] * rms[col1 + rms_off] * rms[col2 + rms_off];
                ++pos;
            }
            sumxy += d[col1] * rhs[col1];
            output[col1 + rms_off + off_x_x] = sumxy * rms[col1 + rms_off] * sumyy;
        }

        return output;
    }

    /**
     * ALGORITHM AS274 APPL. STATIST. (1992) VOL.41, NO. 2.
     * Move variable from position FROM to position TO in an
     * orthogonal reduction produced by AS75.1.
     *
     * @param from initial position
     * @param to destination
     */
    private void vmove(const int& from, int to) 
    {
        double d1;
        double d2;
        double X;
        double d1new;
        double d2new;
        double cbar;
        double sbar;
        double Y;
        int first;
        int inc;
        int m1;
        int m2;
        int mp1;
        int pos;
        bool b_skip_to40 = false;
        if (from == to) 
        {
            return;
        }
        if (!this.rss_set) 
        {
            ss();
        }
        const int count;
        if (from < to) 
        {
            first = from;
            inc = 1;
            count = to - from;
        }
else 
        {
            first = from - 1;
            inc = -1;
            count = from - to;
        }

        int m = first;
        int idx = 0;
        while (idx < count) 
        {
            m1 = m * (nvars + nvars - m - 1) / 2;
            m2 = m1 + nvars - m - 1;
            mp1 = m + 1;

            d1 = d[m];
            d2 = d[mp1];
            // Special cases.
            if (d1 > this.epsilon || d2 > this.epsilon) 
            {
                X = r[m1];
                if (std::abs(X) * std::sqrt(d1) < tol[mp1]) 
                {
                    X = 0.0;
                }
                if (d1 < this.epsilon || std::abs(X) < this.epsilon) 
                {
                    d[m] = d2;
                    d[mp1] = d1;
                    r[m1] = 0.0;
                    for (int col = m + 2; col < nvars; col++) 
                    {
                        ++m1;
                        X = r[m1];
                        r[m1] = r[m2];
                        r[m2] = X;
                        ++m2;
                    }
                    X = rhs[m];
                    rhs[m] = rhs[mp1];
                    rhs[mp1] = X;
                    b_skip_to40 = true;
                    //break;
                }
else if (d2 < this.epsilon) 
                {
                    d[m] = d1 * X * X;
                    r[m1] = 1.0 / X;
                    for (const int& _i = m1 + 1; _i < m1 + nvars - m - 1; _i++) 
                    {
                        r[_i] /= X;
                    }
                    rhs[m] /= X;
                    b_skip_to40 = true;
                    //break;
                }
                if (!b_skip_to40) 
                {
                    d1new = d2 + d1 * X * X;
                    cbar = d2 / d1new;
                    sbar = X * d1 / d1new;
                    d2new = d1 * cbar;
                    d[m] = d1new;
                    d[mp1] = d2new;
                    r[m1] = sbar;
                    for (int col = m + 2; col < nvars; col++) 
                    {
                        ++m1;
                        Y = r[m1];
                        r[m1] = cbar * r[m2] + sbar * Y;
                        r[m2] = Y - X * r[m2];
                        ++m2;
                    }
                    Y = rhs[m];
                    rhs[m] = cbar * rhs[mp1] + sbar * Y;
                    rhs[mp1] = Y - X * rhs[mp1];
                }
            }
            if (m > 0) 
            {
                pos = m;
                for (int row{}; row < m; row++) 
                {
                    X = r[pos];
                    r[pos] = r[pos - 1];
                    r[pos - 1] = X;
                    pos += nvars - row - 2;
                }
            }
            // Adjust variable order (VORDER), the tolerances (TOL) and
            // the vector of residual sums of squares (RSS).
            m1 = vorder[m];
            vorder[m] = vorder[mp1];
            vorder[mp1] = m1;
            X = tol[m];
            tol[m] = tol[mp1];
            tol[mp1] = X;
            rss[m] = rss[mp1] + d[mp1] * rhs[mp1] * rhs[mp1];

            m += inc;
            ++idx;
        }
    }

    /**
     * ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2
     *
     * <p> Re-order the variables in an orthogonal reduction produced by
     * AS75.1 so that the N variables in LIST start at position POS1, * though will not necessarily be in the same order as in LIST.
     * Any variables in VORDER before position POS1 are not moved.
     * Auxiliary routine called: VMOVE. </p>
     *
     * <p>This internal method reorders the regressors.</p>
     *
     * @param list the regressors to move
     * @param pos1 where the list will be placed
     * @return -1 error, 0 everything ok
     */
    private int reorder_regressors(std::vector<int> list, int pos1) 
    {
        int next;
        int i;
        int l;
        if (list.size() < 1 || list.size() > nvars + 1 - pos1) 
        {
            return -1;
        }
        next = pos1;
        i = pos1;
        while (i < nvars) 
        {
            l = vorder[i];
            for (int j{}; j < list.size(); j++) 
            {
                if (l == list[j] && i > next) 
                {
                    this.vmove(i, next);
                    ++next;
                    if (next >= list.size() + pos1) 
                    {
                        return 0;
                    }
else 
                    {
                        break;
                    }
                }
            }
            ++i;
        }
        return 0;
    }

    /**
     * Gets the diagonal of the Hat matrix also known as the leverage matrix.
     *
     * @param  row_data returns the diagonal of the hat matrix for this observation
     * @return the diagonal element of the hatmatrix
     */
    public double get_diagonal_of_hat_matrix(std::vector<double> row_data) 
    {
        std::vector<double> wk = std::vector<double>(this.nvars];
        int pos;
        double total;

        if (row_data.size() > nvars) 
        {
            return std::numeric_limits<double>::quiet_NaN();
        }
        std::vector<double> xrow;
        if (this.has_intercept) 
        {
            xrow = std::vector<double>(row_data.size() + 1];
            xrow[0] = 1.0;
            System.arraycopy(row_data, 0, xrow, 1, row_data.size());
        }
else 
        {
            xrow = row_data;
        }
        double hii = 0.0;
        for (int col{};  col < xrow.size(); col++) 
        {
            if (std::sqrt(d[col]) < tol[col]) 
            {
                wk[col] = 0.0;
            }
else 
            {
                pos = col - 1;
                total = xrow[col];
                for (int row{}; row < col; row++) 
                {
                    total = smart_add(total, -wk[row] * r[pos]);
                    pos += nvars - row - 2;
                }
                wk[col] = total;
                hii = smart_add(hii, (total * total) / d[col]);
            }
        }
        return hii;
    }

    /**
     * Gets the order of the regressors, useful if some type of reordering
     * has been called. Calling regress with std::vector<int>{} args will trigger
     * a reordering.
     *
     * @return std::vector<int> with the current order of the regressors
     */
    public std::vector<int> get_order_of_regressors()
    {
        return vorder.clone();
    }

    /**
     * Conducts a regression on the data in the model, using all regressors.
     *
     * @return Regression_results the structure holding all regression results
     * @exception   - thrown if number of observations is
     * less than the number of variables
     */
    //override
    public Regression_results regress()  
    {
        return regress(this.nvars);
    }

    /**
     * Conducts a regression on the data in the model, using a subset of regressors.
     *
     * @param number_of_regressors many of the regressors to include (either in canonical
     * order, or in the current reordered state)
     * @return Regression_results the structure holding all regression results
     * @exception   - thrown if number of observations is
     * less than the number of variables or number of regressors requested
     * is greater than the regressors in the model
     */
    public Regression_results regress(const int& number_of_regressors)  
    {
        if (this.nobs <= number_of_regressors) 
        {
           throw (
                   Localized_Stat_Formats.NOT_ENOUGH_DATA_FOR_NUMBER_OF_PREDICTORS, this.nobs, number_of_regressors);
        }
        if( number_of_regressors > this.nvars )
        {
            throw (
                    Localized_Stat_Formats.TOO_MANY_REGRESSORS, number_of_regressors, this.nvars);
        }

        tolset();
        singcheck();

        std::vector<double> beta = this.regcf(number_of_regressors);

        ss();

        std::vector<double> cov = this.cov(number_of_regressors);

        int rnk = 0;
        for (int i{}; i < this.lindep.size(); i++) 
        {
            if (!this.lindep[i]) 
            {
                ++rnk;
            }
        }

        bool needs_reorder{};
        for (int i{}; i < number_of_regressors; i++) 
        {
            if (this.vorder[i] != i) 
            {
                needs_reorder = true;
                break;
            }
        }
        if (!needs_reorder) 
        {
            return Regression_results(
                    beta, std::vector<std::vector<double>>{cov}, true, this.nobs, rnk, this.sumy, this.sumsqy, this.sserr, this.has_intercept, false);
        }
else 
        {
            std::vector<double> beta_new = std::vector<double>(beta.size()];
            std::vector<double> cov_new = std::vector<double>(cov.size()];

            std::vector<int> new_indices = int[beta.size()];
            for (int i{}; i < nvars; i++) 
            {
                for (int j{}; j < number_of_regressors; j++) 
                {
                    if (this.vorder[j] == i) 
                    {
                        beta_new[i] = beta[ j];
                        new_indices[i] = j;
                    }
                }
            }

            int idx1 = 0;
            int idx2;
            int _i;
            int _j;
            for (int i{}; i < beta.size(); i++) 
            {
                _i = new_indices[i];
                for (int j{}; j <= i; j++, idx1++) 
                {
                    _j = new_indices[j];
                    if (_i > _j) 
                    {
                        idx2 = _i * (_i + 1) / 2 + _j;
                    }
else 
                    {
                        idx2 = _j * (_j + 1) / 2 + _i;
                    }
                    cov_new[idx1] = cov[idx2];
                }
            }
            return Regression_results(
                    beta_new, std::vector<std::vector<double>>{cov_new}, true, this.nobs, rnk, this.sumy, this.sumsqy, this.sserr, this.has_intercept, false);
        }
    }

    /**
     * Conducts a regression on the data in the model, using regressors in array
     * Calling this method will change the internal order of the regressors
     * and care is required in interpreting the hatmatrix.
     *
     * @param  variables_to_include array of variables to include in regression
     * @return Regression_results the structure holding all regression results
     * @exception   - thrown if number of observations is
     * less than the number of variables, the number of regressors requested
     * is greater than the regressors in the model or a regressor index in
     * regressor array does not exist
     */
    //override
    public Regression_results regress(std::vector<int> variables_to_include)  
    {
        if (variables_to_include.size() > this.nvars) 
        {
            throw (
                    Localized_Stat_Formats.TOO_MANY_REGRESSORS, variables_to_include.size(), this.nvars);
        }
        if (this.nobs <= this.nvars) 
        {
            throw (
                    Localized_Stat_Formats.NOT_ENOUGH_DATA_FOR_NUMBER_OF_PREDICTORS, this.nobs, this.nvars);
        }
        Arrays.sort(variables_to_include);
        int i_exclude = 0;
        for (int i{}; i < variables_to_include.size(); i++) 
        {
            if (i >= this.nvars) 
            {
                throw (
                        hipparchus::exception::Localized_Core_Formats_Type::INDEX_LARGER_THAN_MAX, i, this.nvars);
            }
            if (i > 0 && variables_to_include[i] == variables_to_include[i - 1]) 
            {
                variables_to_include[i] = -1;
                ++i_exclude;
            }
        }
        std::vector<int> series;
        if (i_exclude > 0) 
        {
            int j = 0;
            series = int[variables_to_include.size() - i_exclude];
            for (int i{}; i < variables_to_include.size(); i++) 
            {
                if (variables_to_include[i] > -1) 
                {
                    series[j] = variables_to_include[i];
                    ++j;
                }
            }
        }
else 
        {
            series = variables_to_include;
        }

        reorder_regressors(series, 0);
        tolset();
        singcheck();

        std::vector<double> beta = this.regcf(series.size());

        ss();

        std::vector<double> cov = this.cov(series.size());

        int rnk = 0;
        for (int i{}; i < this.lindep.size(); i++) 
        {
            if (!this.lindep[i]) 
            {
                ++rnk;
            }
        }

        bool needs_reorder{};
        for (int i{}; i < this.nvars; i++) 
        {
            if (this.vorder[i] != series[i]) 
            {
                needs_reorder = true;
                break;
            }
        }
        if (!needs_reorder) 
        {
            return Regression_results(
                    beta, std::vector<std::vector<double>>{cov}, true, this.nobs, rnk, this.sumy, this.sumsqy, this.sserr, this.has_intercept, false);
        }
else 
        {
            std::vector<double> beta_new = std::vector<double>(beta.size()];
            std::vector<int> new_indices = int[beta.size()];
            for (int i{}; i < series.size(); i++) 
            {
                for (int j{}; j < this.vorder.size(); j++) 
                {
                    if (this.vorder[j] == series[i]) 
                    {
                        beta_new[i] = beta[ j];
                        new_indices[i] = j;
                    }
                }
            }
            std::vector<double> cov_new = std::vector<double>(cov.size()];
            int idx1 = 0;
            int idx2;
            int _i;
            int _j;
            for (int i{}; i < beta.size(); i++) 
            {
                _i = new_indices[i];
                for (int j{}; j <= i; j++, idx1++) 
                {
                    _j = new_indices[j];
                    if (_i > _j) 
                    {
                        idx2 = _i * (_i + 1) / 2 + _j;
                    }
else 
                    {
                        idx2 = _j * (_j + 1) / 2 + _i;
                    }
                    cov_new[idx1] = cov[idx2];
                }
            }
            return Regression_results(
                    beta_new, std::vector<std::vector<double>>{cov_new}, true, this.nobs, rnk, this.sumy, this.sumsqy, this.sserr, this.has_intercept, false);
        }
    }
}


