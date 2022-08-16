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
#include <numbers>

  //import org.hipparchus.analysis.Function_Utils;
  //import org.hipparchus.analysis.Univariate_Function;
  //import org.hipparchus.complex.std::complex<double>;
  //import org.hipparchus.exception.;
  //import org.hipparchus.util.Arithmetic_Utils;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Sin_Cos;

/**
* Implements the Fast Cosine Transform for transformation of one-dimensional
* real data sets. For reference, see James S. Walker, <em>Fast Fourier
* Transforms</em>, chapter 3 (ISBN 0849371635).
* <p>
* There are several variants of the discrete cosine transform. The present
* implementation corresponds to DCT-I, with various normalization conventions, * which are specified by the parameter {@link Dct_Normalization}.
* <p>
* DCT-I is equivalent to DFT of an <em>even extension</em> of the data series.
* More precisely, if x<sub>0</sub>, &hellip;, x<sub>N-1</sub> is the data set
* to be cosine transformed, the extended data set
* x<sub>0</sub><sup>&#35;</sup>, &hellip;, x<sub>2N-3</sub><sup>&#35;</sup>
* is defined as follows
* <ul>
* <li>x<sub>k</sub><sup>&#35;</sup> = x<sub>k</sub> if 0 &le; k &lt; N,</li>
* <li>x<sub>k</sub><sup>&#35;</sup> = x<sub>2N-2-k</sub>
* if N &le; k &lt; 2N - 2.</li>
* </ul>
* <p>
* Then, the standard DCT-I y<sub>0</sub>, &hellip;, y<sub>N-1</sub> of the real
* data set x<sub>0</sub>, &hellip;, x<sub>N-1</sub> is equal to <em>half</em>
* of the N first elements of the DFT of the extended data set
* x<sub>0</sub><sup>&#35;</sup>, &hellip;, x<sub>2N-3</sub><sup>&#35;</sup>
* <br/>
* y<sub>n</sub> = (1 / 2) &sum;<sub>k=0</sub><sup>2N-3</sup>
* x<sub>k</sub><sup>&#35;</sup> exp[-2&pi;i nk / (2N - 2)]
* &nbsp;&nbsp;&nbsp;&nbsp;k = 0, &hellip;, N-1.
* <p>
* The present implementation of the discrete cosine transform as a fast cosine
* transform requires the length of the data set to be a power of two plus one
* (N&nbsp;=&nbsp;2<sup>n</sup>&nbsp;+&nbsp;1). Besides, it implicitly assumes
* that the sampled function is even.
*
*/
class Fast_Cosine_Transformer : Real_Transformer
{
private:
	/** The type of DCT to be performed. */
	const Dct_Normalization my_normalization;

public:
	/**
	 * Creates a instance of this class, with various normalization
	 * conventions.
	 *
	 * @param normalization the type of normalization to be applied to the
	 * transformed data
	 */
	Fast_Cosine_Transformer(const Dct_Normalization& normalization)
	{
		my_normalization = normalization;
	}

	/**
	 * {@inherit_doc}
	 *
	 * @ if the length of the data array is
	 * not a power of two plus one
	 */
	 //override
	std::vector<double> transform(const std::vector<double>& f, const Transform_Type& type)
	{
		if (type == Transform_Type.FORWARD)
		{
			if (my_normalization == Dct_Normalization.ORTHOGONAL_DCT_I)
			{
				const double s = std::sqrt(2.0 / (f.size() - 1));
				return Trans_form_Utils.scale_array(fct(f), s);
			}
			return fct(f);
		}
		const double s2 = 2.0 / (f.size() - 1);
		double s1;
		if (my_normalization == Dct_Normalization.ORTHOGONAL_DCT_I)
		{
			s1 = std::sqrt(s2);
		}
		else
		{
			s1 = s2;
		}
		return Trans_form_Utils.scale_array(fct(f), s1);
	}

	/**
	 * {@inherit_doc}
	 *
	 * @org.hipparchus.exception.
	 * if the lower bound is greater than, or equal to the upper bound
	 * @org.hipparchus.exception.
	 * if the number of sample points is negative
	 * @ if the number of sample points is
	 * not a power of two plus one
	 */
	 //override
	std::vector<double> transform(const Univariate_Function& f, const double& min, const double& max, const int& n, const Transform_Type& type)
	{
		const std::vector<double> data = Function_Utils.sample(f, min, max, n);
		return transform(data, type);
	}

protected:
	/**
	 * Perform the FCT algorithm (including inverse).
	 *
	 * @param f the real data array to be transformed
	 * @return the real transformed array
	 * @ if the length of the data array is
	 * not a power of two plus one
	 */
	std::vector<double> fct(std::vector<double> f)
	{
		auto transformed = std::vector<double>(f.size()];

		const int n = f.size() - 1;
		if (!Arithmetic_Utils.is_power_of_two(n))
		{
			throw (Localized_FFT_Formats.NOT_POWER_OF_TWO_PLUS_ONE, static_cast<int>(f.size()));
		}
		if (n == 1) {       // trivial case
			transformed[0] = 0.5 * (f[0] + f[1]);
			transformed[1] = 0.5 * (f[0] - f[1]);
			return transformed;
		}

		// construct a array and perform FFT on it
		auto x = std::vector<double>(n);
		x[0] = 0.5 * (f[0] + f[n]);
		x[n >> 1] = f[n >> 1];
		// temporary variable for transformed[1]
		double t1 = 0.5 * (f[0] - f[n]);
		for (int i{ 1 }; i < (n >> 1); i++)
		{
			const Sin_Cos sc = Sin_Cos(i * std::numbers::pi / n);
			const double& a = 0.5 * (f[i] + f[n - i]);
			const double b = sc.sin() * (f[i] - f[n - i]);
			const double c = sc.cos() * (f[i] - f[n - i]);
			x[i] = a - b;
			x[n - i] = a + b;
			t1 += c;
		}
		Fast_Fourier_Transformer transformer;
		transformer = Fast_Fourier_Transformer(Dft_Normalization.STANDARD);
		std::vector<std::complex<double>> y = transformer.transform(x, Transform_Type.FORWARD);

		// reconstruct the FCT result for the original array
		transformed[0] = y[0].get_real();
		transformed[1] = t1;
		for (int i{ 1 }; i < (n >> 1); i++)
		{
			transformed[2 * i] = y[i].get_real();
			transformed[2 * i + 1] = transformed[2 * i - 1] - y[i].get_imaginary();
		}
		transformed[n] = y[n >> 1].get_real();

		return transformed;
	}
};