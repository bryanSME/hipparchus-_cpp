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
#include <cmath>

  //import org.hipparchus.analysis.Function_Utils;
  //import org.hipparchus.analysis.Univariate_Function;
  //import org.hipparchus.complex.std::complex<double>;
  //import org.hipparchus.exception.;
  //import org.hipparchus.util.Arithmetic_Utils;
  //import org.hipparchus.util.FastMath;

  /**
   * Implements the Fast Sine Transform for transformation of one-dimensional real
   * data sets. For reference, see James S. Walker, <em>Fast Fourier
   * Transforms</em>, chapter 3 (ISBN 0849371635).
   * <p>
   * There are several variants of the discrete sine transform. The present
   * implementation corresponds to DST-I, with various normalization conventions, * which are specified by the parameter {@link Dst_Normalization}.
   * <strong>It should be noted that regardless to the convention, the first
   * element of the dataset to be transformed must be zero.</strong>
   * <p>
   * DST-I is equivalent to DFT of an <em>odd extension</em> of the data series.
   * More precisely, if x<sub>0</sub>, &hellip;, x<sub>N-1</sub> is the data set
   * to be sine transformed, the extended data set x<sub>0</sub><sup>&#35;</sup>, * &hellip;, x<sub>2N-1</sub><sup>&#35;</sup> is defined as follows
   * <ul>
   * <li>x<sub>0</sub><sup>&#35;</sup> = x<sub>0</sub> = 0,</li>
   * <li>x<sub>k</sub><sup>&#35;</sup> = x<sub>k</sub> if 1 &le; k &lt; N,</li>
   * <li>x<sub>N</sub><sup>&#35;</sup> = 0,</li>
   * <li>x<sub>k</sub><sup>&#35;</sup> = -x<sub>2N-k</sub> if N + 1 &le; k &lt;
   * 2N.</li>
   * </ul>
   * <p>
   * Then, the standard DST-I y<sub>0</sub>, &hellip;, y<sub>N-1</sub> of the real
   * data set x<sub>0</sub>, &hellip;, x<sub>N-1</sub> is equal to <em>half</em>
   * of i (the pure imaginary number) times the N first elements of the DFT of the
   * extended data set x<sub>0</sub><sup>&#35;</sup>, &hellip;, * x<sub>2N-1</sub><sup>&#35;</sup> <br />
   * y<sub>n</sub> = (i / 2) &sum;<sub>k=0</sub><sup>2N-1</sup>
   * x<sub>k</sub><sup>&#35;</sup> exp[-2&pi;i nk / (2N)]
   * &nbsp;&nbsp;&nbsp;&nbsp;k = 0, &hellip;, N-1.
   * <p>
   * The present implementation of the discrete sine transform as a fast sine
   * transform requires the length of the data to be a power of two. Besides, * it implicitly assumes that the sampled function is odd. In particular, the
   * first element of the data set must be 0, which is enforced in
   * {@link #transform(Univariate_Function, double, double, int, Transform_Type)}, * after sampling.
   *
   */
class Fast_Sine_Transformer : Real_Transformer
{

private:
	/** The type of DST to be performed. */
	private Dst_Normalization my_normalization;

public:
	/**
	 * Creates a instance of this class, with various normalization conventions.
	 *
	 * @param normalization the type of normalization to be applied to the transformed data
	 */
	Fast_Sine_Transformer(const Dst_Normalization& normalization)
	{
		my_normalization = normalization;
	}

	/**
	 * {@inherit_doc}
	 *
	 * The first element of the specified data set is required to be {@code 0}.
	 *
	 * @ if the length of the data array is
	 *   not a power of two, or the first element of the data array is not zero
	 */
	//override
	std::vector<double> transform(const std::vector<double>& f, const Transform_Type& type)
	{
		if (normalization == Dst_Normalization.ORTHOGONAL_DST_I)
		{
			const double s = std::sqrt(2.0 / f.size());
			return Trans_form_Utils.scale_array(fst(f), s);
		}
		if (type == Transform_Type.FORWARD)
		{
			return fst(f);
		}
		const double s = 2.0 / f.size();
		return Trans_form_Utils.scale_array(fst(f), s);
	}

	/**
	 * {@inherit_doc}
	 *
	 * This implementation enforces {@code f(x) = 0.0} at {@code x = 0.0}.
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
		auto data = Function_Utils.sample(f, min, max, n);
		data[0] = 0.0;
		return transform(data, type);
	}

protected:
	/**
	 * Perform the FST algorithm (including inverse). The first element of the
	 * data set is required to be {@code 0}.
	 *
	 * @param f the real data array to be transformed
	 * @return the real transformed array
	 * @ if the length of the data array is
	 *   not a power of two, or the first element of the data array is not zero
	 */
	std::vector<double> fst(const std::vector<double>& f)
	{
		auto transformed = std::vector<double>(f.size()];

		if (!Arithmetic_Utils.is_power_of_two(f.size()))
		{
			throw (
				Localized_FFT_Formats.NOT_POWER_OF_TWO_CONSIDER_PADDING, static_cast<int>(f.size()));
		}
		if (f[0] != 0.0)
		{
			throw (
				Localized_FFT_Formats.FIRST_ELEMENT_NOT_ZERO, static_cast<double>(f[0]));
		}
		const int n = f.size();
		if (n == 1) {       // trivial case
			transformed[0] = 0.0;
			return transformed;
		}

		// construct a array and perform FFT on it
		auto x = std::vector<double>(n);
		x[0] = 0.0;
		x[n >> 1] = 2.0 * f[n >> 1];
		for (int i{ 1 }; i < (n >> 1); i++)
		{
			const double& a = std::sin(i * std::numbers::pi / n) * (f[i] + f[n - i]);
			const double b = 0.5 * (f[i] - f[n - i]);
			x[i] = a + b;
			x[n - i] = a - b;
		}
		Fast_Fourier_Transformer transformer;
		transformer = Fast_Fourier_Transformer(Dft_Normalization.STANDARD);
		std::vector<std::complex<double>> y = transformer.transform(x, Transform_Type.FORWARD);

		// reconstruct the FST result for the original array
		transformed[0] = 0.0;
		transformed[1] = 0.5 * y[0].get_real();
		for (int i{ 1 }; i < (n >> 1); i++)
		{
			transformed[2 * i] = -y[i].get_imaginary();
			transformed[2 * i + 1] = y[i].get_real() + transformed[2 * i - 1];
		}

		return transformed;
	}
};