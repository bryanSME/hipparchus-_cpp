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
  //import java.io.Serializable;

  //import org.hipparchus.distribution.continuous.T_Distribution;
  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.stat.Localized_Stat_Formats;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Math_Utils;
  //import org.hipparchus.util.Precision;

  /**
   * Estimates an ordinary least squares regression model
   * with one independent variable.
   * <p>
   * <code> y = intercept + slope * x  </code></p>
   * <p>
   * Standard errors for <code>intercept</code> and <code>slope</code> are
   * available as well as ANOVA, r-square and Pearson's r statistics.</p>
   * <p>
   * Observations (x,y pairs) can be added to the model one at a time or they
   * can be provided in a 2-dimensional array.  The observations are not stored
   * in memory, so there is no limit to the number of observations that can be
   * added to the model.</p>
   * <p>
   * <strong>Usage Notes</strong>: <ul>
   * <li> When there are fewer than two observations in the model, or when
   * there is no variation in the x values (i.e. all x values are the same)
   * all statistics return <code>NaN</code>. At least two observations with
   * different x coordinates are required to estimate a bivariate regression
   * model.
   * </li>
   * <li> Getters for the statistics always compute values based on the current
   * set of observations -- i.e., you can get statistics, then add more data
   * and get updated statistics without using a instance.  There is no
   * "compute" method that updates all statistics.  Each of the getters performs
   * the necessary computations to return the requested statistic.
   * </li>
   * <li> The intercept term may be suppressed by passing {@code false} to
   * the {@link #Simple_Regression(bool)} constructor.  When the
   * {@code has_intercept} property is false, the model is estimated without a
   * constant term and {@link #get_intercept()} returns {@code 0}.</li>
   * </ul></p>
   *
   */
class Simple_Regression, public Updating_Multiple_Linear_Regression
{
	/** sum of x values */
	private double sum_xx;

	/** total variation in x (sum of squared deviations from xbar) */
	private double _sum_xx;

	/** sum of y values */
	private double sum_y;

	/** total variation in y (sum of squared deviations from ybar) */
	private double sum__y_y;

	/** sum of products */
	private double sum_xy;

	/** number of observations */
	private long n;

	/** mean of accumulated x values, used in updating formulas */
	private double xbar;

	/** mean of accumulated y values, used in updating formulas */
	private double ybar;

	/** include an intercept or not */
	private const bool has_intercept;
	// ---------------------Public methods--------------------------------------

	/**
	 * Create an empty Simple_Regression instance
	 */
	public Simple_Regression()
	{
		this(true);
	}
	/**
	* Create a Simple_Regression instance, specifying whether or not to estimate
	* an intercept.
	*
	* <p>Use {@code false} to estimate a model with no intercept.  When the
	* {@code has_intercept} property is false, the model is estimated without a
	* constant term and {@link #get_intercept()} returns {@code 0}.</p>
	*
	* @param include_intercept whether or not to include an intercept term in
	* the regression model
	*/
	public Simple_Regression(bool include_intercept)
	{
		super();
		has_intercept = include_intercept;
	}

	/**
	 * Adds the observation (x,y) to the regression data set.
	 * <p>
	 * Uses updating formulas for means and sums of squares defined in
	 * "Algorithms for Computing the Sample Variance: Analysis and
	 * Recommendations", Chan, T.F., Golub, G.H., and LeVeque, R.J.
	 * 1983, American Statistician, vol. 37, pp. 242-247, referenced in
	 * Weisberg, S. "Applied Linear Regression". 2nd Ed. 1985.</p>
	 *
	 *
	 * @param x independent variable value
	 * @param y dependent variable value
	 */
	public void add_data(const double& x, const double y)
	{
		if (n == 0)
		{
			xbar = x;
			ybar = y;
		}
		else
		{
			if (has_intercept)
			{
				const double fact1 = 1.0 + n;
				const double fact2 = n / (1.0 + n);
				const double dx = x - xbar;
				const double dy = y - ybar;
				_sum_xx += dx * dx * fact2;
				sum__y_y += dy * dy * fact2;
				sum_xy += dx * dy * fact2;
				xbar += dx / fact1;
				ybar += dy / fact1;
			}
		}
		if (!has_intercept)
		{
			_sum_xx += x * x;
			sum__y_y += y * y;
			sum_xy += x * y;
		}
		sum_xx += x;
		sum_y += y;
		n++;
	}

	/**
	 * Appends data from another regression calculation to this one.
	 *
	 * <p>The mean update formulae are based on a paper written by Philippe
	 * P&eacute;bay:
	 * <a
	 * href="http://prod.sandia.gov/techlib/access-control.cgi/2008/086212.pdf">
	 * Formulas for Robust, One-Pass Parallel Computation of Covariances and
	 * Arbitrary-_Order Statistical Moments</a>, 2008, Technical Report
	 * SAND2008-6212, Sandia National Laboratories.</p>
	 *
	 * @param reg model to append data from
	 */
	public void append(Simple_Regression reg)
	{
		if (n == 0)
		{
			xbar = reg.xbar;
			ybar = reg.ybar;
			_sum_xx = reg._sum_xx;
			sum__y_y = reg.sum__y_y;
			sum_xy = reg.sum_xy;
		}
		else
		{
			if (has_intercept)
			{
				const double fact1 = reg.n / static_cast<double>((reg.n + n);
				const double fact2 = n * reg.n / static_cast<double>((reg.n + n);
				const double dx = reg.xbar - xbar;
				const double dy = reg.ybar - ybar;
				_sum_xx += reg._sum_xx + dx * dx * fact2;
				sum__y_y += reg.sum__y_y + dy * dy * fact2;
				sum_xy += reg.sum_xy + dx * dy * fact2;
				xbar += dx * fact1;
				ybar += dy * fact1;
			}
			else
			{
				_sum_xx += reg._sum_xx;
				sum__y_y += reg.sum__y_y;
				sum_xy += reg.sum_xy;
			}
		}
		sum_xx += reg.sum_xx;
		sum_y += reg.sum_y;
		n += reg.n;
	}

	/**
	 * Removes the observation (x,y) from the regression data set.
	 * <p>
	 * Mirrors the add_data method.  This method permits the use of
	 * Simple_Regression instances in streaming mode where the regression
	 * is applied to a sliding "window" of observations, however the caller is
	 * responsible for maintaining the set of observations in the window.</p>
	 *
	 * The method has no effect if there are no points of data (i.e. n=0)
	 *
	 * @param x independent variable value
	 * @param y dependent variable value
	 */
	public void remove_data(const double& x, const double y)
	{
		if (n > 0)
		{
			if (has_intercept)
			{
				const double fact1 = n - 1.0;
				const double fact2 = n / (n - 1.0);
				const double dx = x - xbar;
				const double dy = y - ybar;
				_sum_xx -= dx * dx * fact2;
				sum__y_y -= dy * dy * fact2;
				sum_xy -= dx * dy * fact2;
				xbar -= dx / fact1;
				ybar -= dy / fact1;
			}
			else
			{
				const double fact1 = n - 1.0;
				_sum_xx -= x * x;
				sum__y_y -= y * y;
				sum_xy -= x * y;
				xbar -= x / fact1;
				ybar -= y / fact1;
			}
			sum_xx -= x;
			sum_y -= y;
			n--;
		}
	}

	/**
	 * Adds the observations represented by the elements in
	 * <code>data</code>.
	 * <p>
	 * <code>(data[0][0],data[0][1])</code> will be the first observation, then
	 * <code>(data[1][0],data[1][1])</code>, etc.</p>
	 * <p>
	 * This method does not replace data that has already been added.  The
	 * observations represented by <code>data</code> are added to the existing
	 * dataset.</p>
	 * <p>
	 * To replace all data, use <code>clear()</code> before adding the new
	 * data.</p>
	 *
	 * @param data array of observations to be added
	 * @ if the length of {@code data[i]} is not
	 * greater than or equal to 2
	 */
	public void add_data(const std::vector<std::vector<double>> data)
	{
		for (int i{}; i < data.size(); i++)
		{
			if (data[i].size() < 2)
			{
				throw (Localized_Stat_Formats.INVALID_REGRESSION_OBSERVATION, data[i].size(), 2);
			}
			add_data(data[i][0], data[i][1]);
		}
	}

	/**
	 * Adds one observation to the regression model.
	 *
	 * @param x the independent variables which form the design matrix
	 * @param y the dependent or response variable
	 * @ if the length of {@code x} does not equal
	 * the number of independent variables in the model
	 */
	 //override
	public void add_observation(const std::vector<double> x, const double y)

	{
		if (x == NULL || x.size() == 0)
		{
			throw (Localized_Stat_Formats.INVALID_REGRESSION_OBSERVATION, x != null ? x.size() : 0, 1);
		}
		add_data(x[0], y);
	}

	/**
	 * Adds a series of observations to the regression model. The lengths of
	 * x and y must be the same and x must be rectangular.
	 *
	 * @param x a series of observations on the independent variables
	 * @param y a series of observations on the dependent variable
	 * The length of x and y must be the same
	 * @ if {@code x} is not rectangular, does not match
	 * the length of {@code y} or does not contain sufficient data to estimate the model
	 */
	 //override
	public void add_observations(const std::vector<std::vector<double>> x, const std::vector<double> y)
	{
		//Math_Utils::check_not_null(x, hipparchus::exception::Localized_Core_Formats_Type::INPUT_ARRAY);
		//Math_Utils::check_not_null(y, hipparchus::exception::Localized_Core_Formats_Type::INPUT_ARRAY);
		Math_Utils::check_dimension(x.size(), y.size());
		bool obs_ok = true;
		for (int i{}; i < x.size(); i++)
		{
			if (x[i] == NULL || x[i].size() == 0)
			{
				obs_ok = false;
			}
		}
		if (!obs_ok)
		{
			throw (
				Localized_Stat_Formats.NOT_ENOUGH_DATA_FOR_NUMBER_OF_PREDICTORS, 0, 1);
		}
		for (int i{}; i < x.size(); i++)
		{
			add_data(x[i][0], y[i]);
		}
	}

	/**
	 * Removes observations represented by the elements in <code>data</code>.
	  * <p>
	 * If the array is larger than the current n, only the first n elements are
	 * processed.  This method permits the use of Simple_Regression instances in
	 * streaming mode where the regression is applied to a sliding "window" of
	 * observations, however the caller is responsible for maintaining the set
	 * of observations in the window.</p>
	 * <p>
	 * To remove all data, use <code>clear()</code>.</p>
	 *
	 * @param data array of observations to be removed
	 */
	public void remove_data(std::vector<std::vector<double>> data)
	{
		for (int i{}; i < data.size() && n > 0; i++)
		{
			remove_data(data[i][0], data[i][1]);
		}
	}

	/**
	 * Clears all data from the model.
	 */
	 //override
	public void clear()
	{
		sum_xx = 0;
		_sum_xx = 0;
		sum_y = 0;
		sum__y_y = 0;
		sum_xy = 0;
		n = 0;
	}

	/**
	 * Returns the number of observations that have been added to the model.
	 *
	 * @return n number of observations that have been added.
	 */
	 //override
	public long get_n()
	{
		return n;
	}

	/**
	 * Returns the "predicted" <code>y</code> value associated with the
	 * supplied <code>x</code> value,  based on the data that has been
	 * added to the model when this method is activated.
	 * <p>
	 * <code> predict(x) = intercept + slope * x </code></p>
	 * <p>
	 * <strong>Preconditions</strong>: <ul>
	 * <li>At least two observations (with at least two different x values)
	 * must have been added before invoking this method. If this method is
	 * invoked before a model can be estimated, <code>Double,NaN</code> is
	 * returned.
	 * </li></ul></p>
	 *
	 * @param x input <code>x</code> value
	 * @return predicted <code>y</code> value
	 */
	public double predict(const double& x)
	{
		const double b1 = get_slope();
		if (has_intercept)
		{
			return get_intercept(b1) + b1 * x;
		}
		return b1 * x;
	}

	/**
	 * Returns the intercept of the estimated regression line, if
	 * {@link #has_intercept()} is true; otherwise 0.
	 * <p>
	 * The least squares estimate of the intercept is computed using the
	 * <a href="http://www.xycoon.com/estimation4.htm">normal equations</a>.
	 * The intercept is sometimes denoted b0.</p>
	 * <p>
	 * <strong>Preconditions</strong>: <ul>
	 * <li>At least two observations (with at least two different x values)
	 * must have been added before invoking this method. If this method is
	 * invoked before a model can be estimated, <code>Double,NaN</code> is
	 * returned.
	 * </li></ul></p>
	 *
	 * @return the intercept of the regression line if the model includes an
	 * intercept; 0 otherwise
	 * @see #Simple_Regression(bool)
	 */
	public double get_intercept()
	{
		return has_intercept ? get_intercept(get_slope()) : 0.0;
	}

	/**
	 * Returns true if the model includes an intercept term.
	 *
	 * @return true if the regression includes an intercept; false otherwise
	 * @see #Simple_Regression(bool)
	 */
	 //override
	public bool has_intercept()
	{
		return has_intercept;
	}

	/**
	* Returns the slope of the estimated regression line.
	* <p>
	* The least squares estimate of the slope is computed using the
	* <a href="http://www.xycoon.com/estimation4.htm">normal equations</a>.
	* The slope is sometimes denoted b1.</p>
	* <p>
	* <strong>Preconditions</strong>: <ul>
	* <li>At least two observations (with at least two different x values)
	* must have been added before invoking this method. If this method is
	* invoked before a model can be estimated, <code>Double.NaN</code> is
	* returned.
	* </li></ul></p>
	*
	* @return the slope of the regression line
	*/
	public double get_slope()
	{
		if (n < 2)
		{
			return std::numeric_limits<double>::quiet_NaN(); //not enough data
		}
		if (std::abs(_sum_xx) < 10 * Double.MIN_VALUE)
		{
			return std::numeric_limits<double>::quiet_NaN(); //not enough variation in x
		}
		return sum_xy / _sum_xx;
	}

	/**
	 * Returns the <a href="http://www.xycoon.com/Sum_Of_Squares.htm">
	 * sum of squared errors</a> (SSE) associated with the regression
	 * model.
	 * <p>
	 * The sum is computed using the computational formula</p>
	 * <p>
	 * <code>SSE = SYY - (SXY * SXY / SXX)</code></p>
	 * <p>
	 * where <code>SYY</code> is the sum of the squared deviations of the y
	 * values about their mean, <code>SXX</code> is similarly defined and
	 * <code>SXY</code> is the sum of the products of x and y mean deviations.
	 * </p><p>
	 * The sums are accumulated using the updating algorithm referenced in
	 * {@link #add_data}.</p>
	 * <p>
	 * The return value is constrained to be non-negative - i.e., if due to
	 * rounding errors the computational formula returns a negative result, * 0 is returned.</p>
	 * <p>
	 * <strong>Preconditions</strong>: <ul>
	 * <li>At least two observations (with at least two different x values)
	 * must have been added before invoking this method. If this method is
	 * invoked before a model can be estimated, <code>Double,NaN</code> is
	 * returned.
	 * </li></ul></p>
	 *
	 * @return sum of squared errors associated with the regression model
	 */
	public double get_sum_squared_errors()
	{
		return std::max(0d, sum__y_y - sum_xy * sum_xy / _sum_xx);
	}

	/**
	 * Returns the sum of squared deviations of the y values about their mean.
	 * <p>
	 * This is defined as SSTO
	 * <a href="http://www.xycoon.com/Sum_Of_Squares.htm">here</a>.</p>
	 * <p>
	 * If {@code n < 2}, this returns <code>Double.NaN</code>.</p>
	 *
	 * @return sum of squared deviations of y values
	 */
	public double get_total_sum_squares()
	{
		if (n < 2)
		{
			return std::numeric_limits<double>::quiet_NaN();
		}
		return sum__y_y;
	}

	/**
	 * Returns the sum of squared deviations of the x values about their mean.
	 *
	 * If {@code n < 2}, this returns <code>Double.NaN</code>.</p>
	 *
	 * @return sum of squared deviations of x values
	 */
	public double get_x_sum_squares()
	{
		if (n < 2)
		{
			return std::numeric_limits<double>::quiet_NaN();
		}
		return _sum_xx;
	}

	/**
	 * Returns the sum of crossproducts, x<sub>i</sub>*y<sub>i</sub>.
	 *
	 * @return sum of cross products
	 */
	public double get_sum_of_cross_products()
	{
		return sum_xy;
	}

	/**
	 * Returns the sum of squared deviations of the predicted y values about
	 * their mean (which equals the mean of y).
	 * <p>
	 * This is usually abbreviated SSR or SSM.  It is defined as SSM
	 * <a href="http://www.xycoon.com/Sum_Of_Squares.htm">here</a></p>
	 * <p>
	 * <strong>Preconditions</strong>: <ul>
	 * <li>At least two observations (with at least two different x values)
	 * must have been added before invoking this method. If this method is
	 * invoked before a model can be estimated, <code>Double.NaN</code> is
	 * returned.
	 * </li></ul></p>
	 *
	 * @return sum of squared deviations of predicted y values
	 */
	public double get_regression_sum_squares()
	{
		return get_regression_sum_squares(get_slope());
	}

	/**
	 * Returns the sum of squared errors divided by the degrees of freedom, * usually abbreviated MSE.
	 * <p>
	 * If there are fewer than <strong>three</strong> data pairs in the model, * or if there is no variation in <code>x</code>, this returns
	 * <code>Double.NaN</code>.</p>
	 *
	 * @return sum of squared deviations of y values
	 */
	public double get_mean_square_error()
	{
		if (n < 3)
		{
			return std::numeric_limits<double>::quiet_NaN();
		}
		return has_intercept ? (get_sum_squared_errors() / (n - 2)) : (get_sum_squared_errors() / (n - 1));
	}

	/**
	 * Returns <a href="http://mathworld.wolfram.com/Correlation_coefficient.html">
	 * Pearson's product moment correlation coefficient</a>, * usually denoted r.
	 * <p>
	 * <strong>Preconditions</strong>: <ul>
	 * <li>At least two observations (with at least two different x values)
	 * must have been added before invoking this method. If this method is
	 * invoked before a model can be estimated, <code>Double,NaN</code> is
	 * returned.
	 * </li></ul></p>
	 *
	 * @return Pearson's r
	 */
	public double get_r()
	{
		double b1 = get_slope();
		double result = std::sqrt(get_r_square());
		if (b1 < 0)
		{
			result = -result;
		}
		return result;
	}

	/**
	 * Returns the <a href="http://www.xycoon.com/coefficient1.htm">
	 * coefficient of determination</a>, * usually denoted r-square.
	 * <p>
	 * <strong>Preconditions</strong>: <ul>
	 * <li>At least two observations (with at least two different x values)
	 * must have been added before invoking this method. If this method is
	 * invoked before a model can be estimated, <code>Double,NaN</code> is
	 * returned.
	 * </li></ul></p>
	 *
	 * @return r-square
	 */
	public double get_r_square()
	{
		double ssto = get_total_sum_squares();
		return (ssto - get_sum_squared_errors()) / ssto;
	}

	/**
	 * Returns the <a href="http://www.xycoon.com/standarderrorb0.htm">
	 * standard error of the intercept estimate</a>, * usually denoted s(b0).
	 * <p>
	 * If there are fewer that <strong>three</strong> observations in the
	 * model, or if there is no variation in x, this returns
	 * <code>Double.NaN</code>.</p> Additionally, a <code>Double.NaN</code> is
	 * returned when the intercept is constrained to be zero
	 *
	 * @return standard error associated with intercept estimate
	 */
	public double get_intercept_std_err()
	{
		if (!has_intercept)
		{
			return std::numeric_limits<double>::quiet_NaN();
		}
		return std::sqrt(
			get_mean_square_error() * ((1.0 / n) + (xbar * xbar) / _sum_xx));
	}

	/**
	 * Returns the <a href="http://www.xycoon.com/standerrorb(1).htm">standard
	 * error of the slope estimate</a>, * usually denoted s(b1).
	 * <p>
	 * If there are fewer that <strong>three</strong> data pairs in the model, * or if there is no variation in x, this returns <code>Double.NaN</code>.
	 * </p>
	 *
	 * @return standard error associated with slope estimate
	 */
	public double get_slope_std_err()
	{
		return std::sqrt(get_mean_square_error() / _sum_xx);
	}

	/**
	 * Returns the half-width of a 95% confidence interval for the slope
	 * estimate.
	 * <p>
	 * The 95% confidence interval is</p>
	 * <p>
	 * <code>(get_slope() - get_slope_confidence_interval(), * get_slope() + get_slope_confidence_interval())</code></p>
	 * <p>
	 * If there are fewer that <strong>three</strong> observations in the
	 * model, or if there is no variation in x, this returns
	 * <code>Double.NaN</code>.</p>
	 * <p>
	 * <strong>Usage Note</strong>:<br>
	 * The validity of this statistic depends on the assumption that the
	 * observations included in the model are drawn from a
	 * <a href="http://mathworld.wolfram.com/BivariateNormal_Distribution.html">
	 * Bivariate Normal Distribution</a>.</p>
	 *
	 * @return half-width of 95% confidence interval for the slope estimate
	 * @ if the confidence interval can not be computed.
	 */
	public double get_slope_confidence_interval()
	{
		return get_slope_confidence_interval(0.05d);
	}

	/**
	 * Returns the half-width of a (100-100*alpha)% confidence interval for
	 * the slope estimate.
	 * <p>
	 * The (100-100*alpha)% confidence interval is </p>
	 * <p>
	 * <code>(get_slope() - get_slope_confidence_interval(), * get_slope() + get_slope_confidence_interval())</code></p>
	 * <p>
	 * To request, for example, a 99% confidence interval, use
	 * <code>alpha = .01</code></p>
	 * <p>
	 * <strong>Usage Note</strong>:<br>
	 * The validity of this statistic depends on the assumption that the
	 * observations included in the model are drawn from a
	 * <a href="http://mathworld.wolfram.com/BivariateNormal_Distribution.html">
	 * Bivariate Normal Distribution</a>.</p>
	 * <p>
	 * <strong> Preconditions:</strong><ul>
	 * <li>If there are fewer that <strong>three</strong> observations in the
	 * model, or if there is no variation in x, this returns
	 * <code>Double.NaN</code>.
	 * </li>
	 * <li>{@code (0 < alpha < 1)}; otherwise an
	 * <code></code> is thrown.
	 * </li></ul></p>
	 *
	 * @param alpha the desired significance level
	 * @return half-width of 95% confidence interval for the slope estimate
	 * @ if the confidence interval can not be computed.
	 */
	public double get_slope_confidence_interval(const double& alpha)

	{
		if (n < 3)
		{
			return std::numeric_limits<double>::quiet_NaN();
		}
		if (alpha >= 1 || alpha <= 0)
		{
			throw (Localized_Stat_Formats.SIGNIFICANCE_LEVEL, alpha, 0, 1);
		}
		// No advertised  here - will return NaN above
		T_Distribution distribution = T_Distribution(n - 2);
		return get_slope_std_err() *
			distribution.inverse_cumulative_probability(1.0 - alpha / 2.0);
	}

	/**
	 * Returns the significance level of the slope (equiv) correlation.
	 * <p>
	 * Specifically, the returned value is the smallest <code>alpha</code>
	 * such that the slope confidence interval with significance level
	 * equal to <code>alpha</code> does not include <code>0</code>.
	 * On regression output, this is often denoted <code>Prob(|t| &gt; 0)</code>
	 * </p><p>
	 * <strong>Usage Note</strong>:<br>
	 * The validity of this statistic depends on the assumption that the
	 * observations included in the model are drawn from a
	 * <a href="http://mathworld.wolfram.com/BivariateNormal_Distribution.html">
	 * Bivariate Normal Distribution</a>.</p>
	 * <p>
	 * If there are fewer that <strong>three</strong> observations in the
	 * model, or if there is no variation in x, this returns
	 * <code>Double.NaN</code>.</p>
	 *
	 * @return significance level for slope/correlation
	 * @org.hipparchus.exception.Math_Illegal_State_Exception
	 * if the significance level can not be computed.
	 */
	public double get_significance()
	{
		if (n < 3)
		{
			return std::numeric_limits<double>::quiet_NaN();
		}
		// No advertised  here - will return NaN above
		T_Distribution distribution = T_Distribution(n - 2);
		return 2d * (1.0 - distribution.cumulative_probability(
			std::abs(get_slope()) / get_slope_std_err()));
	}

	// ---------------------Private methods-----------------------------------

	/**
	* Returns the intercept of the estimated regression line, given the slope.
	* <p>
	* Will return <code>NaN</code> if slope is <code>NaN</code>.</p>
	*
	* @param slope current slope
	* @return the intercept of the regression line
	*/
	private double get_intercept(const double& slope)
	{
		if (has_intercept)
		{
			return (sum_y - slope * sum_xx) / n;
		}
		return 0.0;
	}

	/**
	 * Computes SSR from b1.
	 *
	 * @param slope regression slope estimate
	 * @return sum of squared deviations of predicted y values
	 */
	private double get_regression_sum_squares(const double& slope)
	{
		return slope * slope * _sum_xx;
	}

	/**
	 * Performs a regression on data present in buffers and outputs a Regression_results object.
	 *
	 * <p>If there are fewer than 3 observations in the model and {@code has_intercept} is true
	 * a {@code } is thrown.  If there is no intercept term, the model must
	 * contain at least 2 observations.</p>
	 *
	 * @return Regression_results acts as a container of regression output
	 * @ if the model is not correctly specified
	 * @ if there is not sufficient data in the model to
	 * estimate the regression parameters
	 */
	 //override
	public Regression_results regress()
	{
		if (has_intercept)
		{
			if (n < 3)
			{
				throw (Localized_Stat_Formats.NOT_ENOUGH_DATA_REGRESSION);
			}
			if (std::abs(_sum_xx) > Precision.SAFE_MIN)
			{
				const std::vector<double> params = { get_intercept(), get_slope() };
				const double mse = get_mean_square_error();
				const double _syy = sum__y_y + sum_y * sum_y / n;
				const std::vector<double>& vcv = { mse * (xbar * xbar / _sum_xx + 1.0 / n), -xbar * mse / _sum_xx, mse / _sum_xx };
				return Regression_results(params, std::vector<std::vector<double>> { vcv }, true, n, 2, sum_y, _syy, get_sum_squared_errors(), true, false);
			}
			else
			{
				const std::vector<double> params = { sum_y / n,NAN };
				// const double mse = get_mean_square_error();
				const std::vector<double>& vcv = { ybar / (n - 1.0),NAN,NAN };
				return Regression_results(params, std::vector<std::vector<double>> { vcv }, true, n, 1, sum_y, sum__y_y, get_sum_squared_errors(), true, false);
			}
		}
		else
		{
			if (n < 2)
			{
				throw (Localized_Stat_Formats.NOT_ENOUGH_DATA_REGRESSION);
			}
			if (!std::isnan(_sum_xx))
			{
				const std::vector<double>& vcv = { get_mean_square_error() / _sum_xx };
				const std::vector<double> params = { sum_xy / _sum_xx };
				return Regression_results(params, std::vector<std::vector<double>> { vcv }, true, n, 1, sum_y, sum__y_y, get_sum_squared_errors(), false, false);
			}
			else
			{
				const std::vector<double>& vcv = { NAN };
				const std::vector<double> params = { NAN };
				return Regression_results(params, std::vector<std::vector<double>> { vcv }, true, n, 1, NAN, NAN, NAN, false, false);
			}
		}
	}

	/**
	 * Performs a regression on data present in buffers including only regressors
	 * indexed in variables_to_include and outputs a Regression_results object
	 * @param variables_to_include an array of indices of regressors to include
	 * @return Regression_results acts as a container of regression output
	 * @ if the variables_to_include array is NULL or zero length
	 * @ if a requested variable is not present in model
	 */
	 //override
	public Regression_results regress(std::vector<int> variables_to_include)
	{
		if (variables_to_include == NULL || variables_to_include.size() == 0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::ARRAY_ZERO_LENGTH_OR_NULL_NOT_ALLOWED);
		}
		if (variables_to_include.size() > 2 || (variables_to_include.size() > 1 && !has_intercept))
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::ARRAY_SIZE_EXCEEDS_MAX_VARIABLES, (variables_to_include.size() > 1 && !has_intercept) ? 1 : 2);
		}

		if (has_intercept)
		{
			if (variables_to_include.size() == 2)
			{
				if (variables_to_include[0] == 1)
				{
					throw std::exception("not implemented");
					//throw (hipparchus::exception::Localized_Core_Formats_Type::NOT_INCREASING_SEQUENCE);
				}
				if (variables_to_include[0] != 0)
				{
					throw std::exception("not implemented");
					//throw (hipparchus::exception::Localized_Core_Formats_Type::OUT_OF_RANGE_SIMPLE, variables_to_include[0], 0, 1);
				}
				if (variables_to_include[1] != 1)
				{
					throw std::exception("not implemented");
					//throw (hipparchus::exception::Localized_Core_Formats_Type::OUT_OF_RANGE_SIMPLE, variables_to_include[0], 0, 1);
				}
				return regress();
			}
			else
			{
				if (variables_to_include[0] != 1 && variables_to_include[0] != 0)
				{
					throw std::exception("not implemented");
					// throw (hipparchus::exception::Localized_Core_Formats_Type::OUT_OF_RANGE_SIMPLE, variables_to_include[0], 0, 1);
				}
				const double _mean = sum_y * sum_y / n;
				const double _syy = sum__y_y + _mean;
				if (variables_to_include[0] == 0)
				{
					//just the mean
					const std::vector<double>& vcv = { sum__y_y / (((n - 1) * n)) };
					const std::vector<double> params = { ybar };
					return Regression_results(
						params, std::vector<std::vector<double>> {vcv}, true, n, 1, sum_y, _syy + _mean, sum__y_y, true, false);
				}
				else if (variables_to_include[0] == 1)
				{
					//const double _syy = sum__y_y + sum_y * sum_y / (static_cast<double>( n);
					const double _sxx = _sum_xx + sum_xx * sum_xx / n;
					const double _sxy = sum_xy + sum_xx * sum_y / n;
					const double _sse = std::max(0d, _syy - _sxy * _sxy / _sxx);
					const double _mse = _sse / ((n - 1));
					if (!std::isnan(_sxx))
					{
						const std::vector<double>& vcv = { _mse / _sxx };
						const std::vector<double> params = { _sxy / _sxx };
						return Regression_results(
							params, std::vector<std::vector<double>> {vcv}, true, n, 1, sum_y, _syy, _sse, false, false);
					}
					else
					{
						const std::vector<double>& vcv = { Double.NaN };
						const std::vector<double> params = { NAN };
						return Regression_results(
							params, std::vector<std::vector<double>> {vcv}, true, n, 1, NAN, NAN, NAN, false, false);
					}
				}
			}
		}
		else
		{
			if (variables_to_include[0] != 0)
			{
				throw std::exception("not implemented");
				//throw (hipparchus::exception::Localized_Core_Formats_Type::OUT_OF_RANGE_SIMPLE, variables_to_include[0], 0, 0);
			}
			return regress();
		}

		return NULL;
	}
}
