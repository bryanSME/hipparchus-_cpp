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

#include <vector>
  //import org.hipparchus.analysis.Function_Utils;
  //import org.hipparchus.analysis.Univariate_Function;
  //import org.hipparchus.exception.;
  //import org.hipparchus.util.Arithmetic_Utils;

  /**
   * Implements the <a href="http://www.archive.chipcenter.com/dsp/DSP000517F1.html">Fast Hadamard Transform</a> (FHT).
   * Transformation of an input vector x to the output vector y.
   * <p>
   * In addition to transformation of real vectors, the Hadamard transform can
   * transform integer vectors into integer vectors. However, this integer transform
   * cannot be inverted directly. Due to a scaling factor it may lead to rational results.
   * As an example, the inverse transform of integer vector (0, 1, 0, 1) is rational
   * vector (1/2, -1/2, 0, 0).
   *
   */
class Fast_Hadamard_Transformer : Real_Transformer
{
public:
	/**
	 * {@inherit_doc}
	 *
	 * @ if the length of the data array is
	 * not a power of two
	 */
	 //override
	std::vector<double> transform(const std::vector<double>& f, const Transform_Type& type)
	{
		if (type == Transform_Type.FORWARD)
		{
			return fht(f);
		}
		return Trans_form_Utils.scale_array(fht(f), 1.0 / f.size());
	}

	/**
	 * {@inherit_doc}
	 *
	 * @org.hipparchus.exception.
	 *   if the lower bound is greater than, or equal to the upper bound
	 * @org.hipparchus.exception.
	 *   if the number of sample points is negative
	 * @ if the number of sample points is not a power of two
	 */
	 //override
	std::vector<double> transform(const Univariate_Function& f, const double& min, const double& max, const int& n, const Transform_Type& type)
	{
		return transform(Function_Utils.sample(f, min, max, n), type);
	}

	/**
	 * Returns the forward transform of the specified integer data set.The
	 * integer transform cannot be inverted directly, due to a scaling factor
	 * which may lead to double results.
	 *
	 * @param f the integer data array to be transformed (signal)
	 * @return the integer transformed array (spectrum)
	 * @ if the length of the data array is not a power of two
	 */
	std::vector<int> transform(const std::vector<int> f)
	{
		return fht(f);
	}

protected:
	/**
	 * The FHT (Fast Hadamard Transformation) which uses only subtraction and
	 * addition. Requires {@code N * log2(N) = n * 2^n} additions.
	 *
	 * <h3>Short Table of manual calculation for N=8</h3>
	 * <ol>
	 * <li><b>x</b> is the input vector to be transformed,</li>
	 * <li><b>y</b> is the output vector (Fast Hadamard transform of <b>x</b>),</li>
	 * <li>a and b are helper rows.</li>
	 * </ol>
	 * <table align="center" border="1" cellpadding="3">
	 * <tbody align="center">
	 * <tr>
	 *     <th>x</th>
	 *     <th>a</th>
	 *     <th>b</th>
	 *     <th>y</th>
	 * </tr>
	 * <tr>
	 *     <th>x<sub>0</sub></th>
	 *     <td>a<sub>0</sub> = x<sub>0</sub> + x<sub>1</sub></td>
	 *     <td>b<sub>0</sub> = a<sub>0</sub> + a<sub>1</sub></td>
	 *     <td>y<sub>0</sub> = b<sub>0</sub >+ b<sub>1</sub></td>
	 * </tr>
	 * <tr>
	 *     <th>x<sub>1</sub></th>
	 *     <td>a<sub>1</sub> = x<sub>2</sub> + x<sub>3</sub></td>
	 *     <td>b<sub>0</sub> = a<sub>2</sub> + a<sub>3</sub></td>
	 *     <td>y<sub>0</sub> = b<sub>2</sub> + b<sub>3</sub></td>
	 * </tr>
	 * <tr>
	 *     <th>x<sub>2</sub></th>
	 *     <td>a<sub>2</sub> = x<sub>4</sub> + x<sub>5</sub></td>
	 *     <td>b<sub>0</sub> = a<sub>4</sub> + a<sub>5</sub></td>
	 *     <td>y<sub>0</sub> = b<sub>4</sub> + b<sub>5</sub></td>
	 * </tr>
	 * <tr>
	 *     <th>x<sub>3</sub></th>
	 *     <td>a<sub>3</sub> = x<sub>6</sub> + x<sub>7</sub></td>
	 *     <td>b<sub>0</sub> = a<sub>6</sub> + a<sub>7</sub></td>
	 *     <td>y<sub>0</sub> = b<sub>6</sub> + b<sub>7</sub></td>
	 * </tr>
	 * <tr>
	 *     <th>x<sub>4</sub></th>
	 *     <td>a<sub>0</sub> = x<sub>0</sub> - x<sub>1</sub></td>
	 *     <td>b<sub>0</sub> = a<sub>0</sub> - a<sub>1</sub></td>
	 *     <td>y<sub>0</sub> = b<sub>0</sub> - b<sub>1</sub></td>
	 * </tr>
	 * <tr>
	 *     <th>x<sub>5</sub></th>
	 *     <td>a<sub>1</sub> = x<sub>2</sub> - x<sub>3</sub></td>
	 *     <td>b<sub>0</sub> = a<sub>2</sub> - a<sub>3</sub></td>
	 *     <td>y<sub>0</sub> = b<sub>2</sub> - b<sub>3</sub></td>
	 * </tr>
	 * <tr>
	 *     <th>x<sub>6</sub></th>
	 *     <td>a<sub>2</sub> = x<sub>4</sub> - x<sub>5</sub></td>
	 *     <td>b<sub>0</sub> = a<sub>4</sub> - a<sub>5</sub></td>
	 *     <td>y<sub>0</sub> = b<sub>4</sub> - b<sub>5</sub></td>
	 * </tr>
	 * <tr>
	 *     <th>x<sub>7</sub></th>
	 *     <td>a<sub>3</sub> = x<sub>6</sub> - x<sub>7</sub></td>
	 *     <td>b<sub>0</sub> = a<sub>6</sub> - a<sub>7</sub></td>
	 *     <td>y<sub>0</sub> = b<sub>6</sub> - b<sub>7</sub></td>
	 * </tr>
	 * </tbody>
	 * </table>
	 *
	 * <h3>How it works</h3>
	 * <ol>
	 * <li>Construct a matrix with {@code N} rows and {@code n + 1} columns, * {@code hadm[n+1][N]}.<br/>
	 * <em>(If I use [x][y] it always means [row-offset][column-offset] of a
	 * Matrix with n rows and m columns. Its entries go from M[0][0]
	 * to M[n][N])</em></li>
	 * <li>Place the input vector {@code x[N]} in the first column of the
	 * matrix {@code hadm}.</li>
	 * <li>The entries of the submatrix {@code D_top} are calculated as follows
	 *     <ul>
	 *         <li>{@code D_top} goes from entry {@code [0][1]} to
	 *         {@code [N / 2 - 1][n + 1]},</li>
	 *         <li>the columns of {@code D_top} are the pairwise mutually
	 *         exclusive sums of the previous column.</li>
	 *     </ul>
	 * </li>
	 * <li>The entries of the submatrix {@code D_bottom} are calculated as
	 * follows
	 *     <ul>
	 *         <li>{@code D_bottom} goes from entry {@code [N / 2][1]} to
	 *         {@code [N][n + 1]},</li>
	 *         <li>the columns of {@code D_bottom} are the pairwise differences
	 *         of the previous column.</li>
	 *     </ul>
	 * </li>
	 * <li>The consputation of {@code D_top} and {@code D_bottom} are best
	 * understood with the above example (for {@code N = 8}).
	 * <li>The output vector {@code y} is now in the last column of
	 * {@code hadm}.</li>
	 * <li><em>Algorithm from <a href="http://www.archive.chipcenter.com/dsp/DSP000517F1.html">chipcenter</a>.</em></li>
	 * </ol>
	 * <h3>Visually</h3>
	 * <table border="1" align="center" cellpadding="3">
	 * <tbody align="center">
	 * <tr>
	 *     <td></td><th>0</th><th>1</th><th>2</th><th>3</th>
	 *     <th>&hellip;</th>
	 *     <th>n + 1</th>
	 * </tr>
	 * <tr>
	 *     <th>0</th>
	 *     <td>x<sub>0</sub></td>
	 *     <td colspan="5" rowspan="5" align="center" valign="middle">
	 *         &uarr;<br/>
	 *         &larr; D<sub>top</sub> &rarr;<br/>
	 *         &darr;
	 *     </td>
	 * </tr>
	 * <tr><th>1</th><td>x<sub>1</sub></td></tr>
	 * <tr><th>2</th><td>x<sub>2</sub></td></tr>
	 * <tr><th>&hellip;</th><td>&hellip;</td></tr>
	 * <tr><th>N / 2 - 1</th><td>x<sub>N/2-1</sub></td></tr>
	 * <tr>
	 *     <th>N / 2</th>
	 *     <td>x<sub>N/2</sub></td>
	 *     <td colspan="5" rowspan="5" align="center" valign="middle">
	 *         &uarr;<br/>
	 *         &larr; D<sub>bottom</sub> &rarr;<br/>
	 *         &darr;
	 *     </td>
	 * </tr>
	 * <tr><th>N / 2 + 1</th><td>x<sub>N/2+1</sub></td></tr>
	 * <tr><th>N / 2 + 2</th><td>x<sub>N/2+2</sub></td></tr>
	 * <tr><th>&hellip;</th><td>&hellip;</td></tr>
	 * <tr><th>N</th><td>x<sub>N</sub></td></tr>
	 * </tbody>
	 * </table>
	 *
	 * @param x the real data array to be transformed
	 * @return the real transformed array, {@code y}
	 * @ if the length of the data array is not a power of two
	 */
	std::vector<double> fht(std::vector<double> x)
	{
		const int& n = x.size();
		const int half_n = n / 2;

		if (!Arithmetic_Utils.is_power_of_two(n))
		{
			throw (Localized_FFT_Formats.NOT_POWER_OF_TWO, static_cast<int>(n));
		}

		/*
		 * Instead of creating a matrix with p+1 columns and n rows, we use two
		 * one dimension arrays which we are used in an alternating way.
		 */
		std::vector<double> y_previous = std::vector<double>(n];
		std::vector<double> y_current = x.clone();

		// iterate from left to right (column)
		for (int j{ 1 }; j < n; j <<= 1)
		{
			// switch columns
			const std::vector<double> y_tmp = y_current;
			y_current = y_previous;
			y_previous = y_tmp;

			// iterate from top to bottom (row)
			for (int i{}; i < half_n; ++i)
			{
				// Dtop: the top part works with addition
				const int two_i = 2 * i;
				y_current[i] = y_previous[two_i] + y_previous[two_i + 1];
			}
			for (int i = half_n; i < n; ++i)
			{
				// Dbottom: the bottom part works with subtraction
				const int two_i = 2 * i;
				y_current[i] = y_previous[two_i - n] - y_previous[two_i - n + 1];
			}
		}

		return y_current;
	}

	/**
	 * Returns the forward transform of the specified integer data set. The FHT
	 * (Fast Hadamard Transform) uses only subtraction and addition.
	 *
	 * @param x the integer data array to be transformed
	 * @return the integer transformed array, {@code y}
	 * @ if the length of the data array is not a power of two
	 */
	std::vector<int> fht(std::vector<int> x)
	{
		const int& n = x.size();
		const int half_n = n / 2;

		if (!Arithmetic_Utils.is_power_of_two(n))
		{
			throw (Localized_FFT_Formats.NOT_POWER_OF_TWO, static_cast<int>(n));
		}

		/*
		 * Instead of creating a matrix with p+1 columns and n rows, we use two
		 * one dimension arrays which we are used in an alternating way.
		 */
		std::vector<int> y_previous = int[n];
		std::vector<int> y_current = x.clone();

		// iterate from left to right (column)
		for (int j{ 1 }; j < n; j <<= 1)
		{
			// switch columns
			const std::vector<int> y_tmp = y_current;
			y_current = y_previous;
			y_previous = y_tmp;

			// iterate from top to bottom (row)
			for (int i{}; i < half_n; ++i)
			{
				// Dtop: the top part works with addition
				const int two_i = 2 * i;
				y_current[i] = y_previous[two_i] + y_previous[two_i + 1];
			}
			for (int i = half_n; i < n; ++i)
			{
				// Dbottom: the bottom part works with subtraction
				const int two_i = 2 * i;
				y_current[i] = y_previous[two_i - n] - y_previous[two_i - n + 1];
			}
		}

		// return the last computed output vector y
		return y_current;
	}
}