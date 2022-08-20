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
  //package org.hipparchus.transform;

  //import java.io.Serializable;

  //import org.hipparchus.analysis.Function_Utils;
  //import org.hipparchus.analysis.Univariate_Function;
  //import org.hipparchus.complex.std::complex<double>;
  //import org.hipparchus.exception.;
  //import org.hipparchus.exception.Math_Runtime_Exception;
  //import org.hipparchus.util.Arithmetic_Utils;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Math_Arrays;
  //import org.hipparchus.util.Math_Utils;
#include <vector>
#include <cmath>
#include <complex>
#include <assert.h>
#include "DftNormalization.h"
#include "TransformType.h"
#include "../../core/util/MathArrays.h"
#include "../../core/util/MathUtils.h"
#include "../../core/util/ArithmeticUtils.h"
#include "../../core/complex/Complex.h"
#include "../../core/analysis/UnivariateFunction.h"

  /**
   * Implements the Fast Fourier Transform for transformation of one-dimensional
   * real or complex data sets. For reference, see <em>Applied Numerical Linear
   * Algebra</em>, ISBN 0898713897, chapter 6.
   * <p>
   * There are several variants of the discrete Fourier transform, with various
   * normalization conventions, which are specified by the parameter
   * {@link Dft_Normalization}.
   * <p>
   * The current implementation of the discrete Fourier transform as a fast
   * Fourier transform requires the length of the data set to be a power of 2.
   * This greatly simplifies and speeds up the code. Users can pad the data with
   * zeros to meet this requirement. There are other flavors of FFT, for
   * reference, see S. Winograd, * <i>On computing the discrete Fourier transform</i>, Mathematics of
   * Computation, 32 (1978), 175 - 199.
   *
   * @see Dft_Normalization
   */
class Fast_Fourier_Transformer
{
private:
	/**
	 * {@code W_SUB_N_R[i]} is the real part of
	 * {@code exp(- 2 * i * pi / n)}:
	 * {@code W_SUB_N_R[i] = cos(2 * pi/ n)}, where {@code n = 2^i}.
	 */
	static const std::vector<double> W_SUB_N_R = 
	{
		0x1.0p0, -0x1.0p0, 0x1.1a62633145c07p-54, 0x1.6a09e667f3bcdp-1,
		0x1.d906bcf328d46p-1, 0x1.f6297cff75cbp-1, 0x1.fd88da3d12526p-1, 0x1.ff621e3796d7ep-1,
		0x1.ffd886084cd0dp-1, 0x1.fff62169b92dbp-1, 0x1.fffd8858e8a92p-1, 0x1.ffff621621d02p-1,
		0x1.ffffd88586ee6p-1, 0x1.fffff62161a34p-1, 0x1.fffffd8858675p-1, 0x1.ffffff621619cp-1,
		0x1.ffffffd885867p-1, 0x1.fffffff62161ap-1, 0x1.fffffffd88586p-1, 0x1.ffffffff62162p-1,
		0x1.ffffffffd8858p-1, 0x1.fffffffff6216p-1, 0x1.fffffffffd886p-1, 0x1.ffffffffff621p-1,
		0x1.ffffffffffd88p-1, 0x1.fffffffffff62p-1, 0x1.fffffffffffd9p-1, 0x1.ffffffffffff6p-1,
		0x1.ffffffffffffep-1, 0x1.fffffffffffffp-1, 0x1.0p0, 0x1.0p0,
		0x1.0p0, 0x1.0p0, 0x1.0p0, 0x1.0p0,
		0x1.0p0, 0x1.0p0, 0x1.0p0, 0x1.0p0,
		0x1.0p0, 0x1.0p0, 0x1.0p0, 0x1.0p0,
		0x1.0p0, 0x1.0p0, 0x1.0p0, 0x1.0p0,
		0x1.0p0, 0x1.0p0, 0x1.0p0, 0x1.0p0,
		0x1.0p0, 0x1.0p0, 0x1.0p0, 0x1.0p0,
		0x1.0p0, 0x1.0p0, 0x1.0p0, 0x1.0p0,
		0x1.0p0, 0x1.0p0, 0x1.0p0
	};

	/**
	 * {@code W_SUB_N_I[i]} is the imaginary part of
	 * {@code exp(- 2 * i * pi / n)}:
	 * {@code W_SUB_N_I[i] = -sin(2 * pi/ n)}, where {@code n = 2^i}.
	 */
	static const std::vector<double> W_SUB_N_I =
	{ 0x1.1a62633145c07p-52, -0x1.1a62633145c07p-53, -0x1.0p0, -0x1.6a09e667f3bccp-1
	, -0x1.87de2a6aea963p-2, -0x1.8f8b83c69a60ap-3, -0x1.917a6bc29b42cp-4, -0x1.91f65f10dd814p-5
	, -0x1.92155f7a3667ep-6, -0x1.921d1fcdec784p-7, -0x1.921f0fe670071p-8, -0x1.921f8becca4bap-9
	, -0x1.921faaee6472dp-10, -0x1.921fb2aecb36p-11, -0x1.921fb49ee4ea6p-12, -0x1.921fb51aeb57bp-13
	, -0x1.921fb539ecf31p-14, -0x1.921fb541ad59ep-15, -0x1.921fb5439d73ap-16, -0x1.921fb544197ap-17
	, -0x1.921fb544387bap-18, -0x1.921fb544403c1p-19, -0x1.921fb544422c2p-20, -0x1.921fb54442a83p-21
	, -0x1.921fb54442c73p-22, -0x1.921fb54442cefp-23, -0x1.921fb54442d0ep-24, -0x1.921fb54442d15p-25
	, -0x1.921fb54442d17p-26, -0x1.921fb54442d18p-27, -0x1.921fb54442d18p-28, -0x1.921fb54442d18p-29
	, -0x1.921fb54442d18p-30, -0x1.921fb54442d18p-31, -0x1.921fb54442d18p-32, -0x1.921fb54442d18p-33
	, -0x1.921fb54442d18p-34, -0x1.921fb54442d18p-35, -0x1.921fb54442d18p-36, -0x1.921fb54442d18p-37
	, -0x1.921fb54442d18p-38, -0x1.921fb54442d18p-39, -0x1.921fb54442d18p-40, -0x1.921fb54442d18p-41
	, -0x1.921fb54442d18p-42, -0x1.921fb54442d18p-43, -0x1.921fb54442d18p-44, -0x1.921fb54442d18p-45
	, -0x1.921fb54442d18p-46, -0x1.921fb54442d18p-47, -0x1.921fb54442d18p-48, -0x1.921fb54442d18p-49
	, -0x1.921fb54442d18p-50, -0x1.921fb54442d18p-51, -0x1.921fb54442d18p-52, -0x1.921fb54442d18p-53
	, -0x1.921fb54442d18p-54, -0x1.921fb54442d18p-55, -0x1.921fb54442d18p-56, -0x1.921fb54442d18p-57
	, -0x1.921fb54442d18p-58, -0x1.921fb54442d18p-59, -0x1.921fb54442d18p-60 };

	/** The type of DFT to be performed. */
	const Dft_Normalization my_normalization;

	/**
	 * Performs identical index bit reversal shuffles on two arrays of identical
	 * size. Each element in the array is swapped with another element based on
	 * the bit-reversal of the index. For example, in an array with length 16, * item at binary index 0011 (decimal 3) would be swapped with the item at
	 * binary index 1100 (decimal 12).
	 *
	 * @param a the first array to be shuffled
	 * @param b the second array to be shuffled
	 */
	static void bit_reversal_shuffle2(std::vector<double> a, std::vector<double> b)
	{
		const int n = a.size();
		assert(b.size() == n);
		const int half_of_n = n >> 1;

		int j = 0;
		for (int i{}; i < n; i++)
		{
			if (i < j)
			{
				// swap indices i & j
				double temp = a[i];
				a[i] = a[j];
				a[j] = temp;

				temp = b[i];
				b[i] = b[j];
				b[j] = temp;
			}

			int k = half_of_n;
			while (k <= j && k > 0)
			{
				j -= k;
				k >>= 1;
			}
			j += k;
		}
	}

	/**
	 * Applies the proper normalization to the specified transformed data.
	 *
	 * @param data_r_i the unscaled transformed data
	 * @param normalization the normalization to be applied
	 * @param type the type of transform (forward, inverse) which resulted in the specified data
	 */
	static void normalize_transformed_data(const std::vector<std::vector<double>> data_r_i, const Dft_Normalization normalization, const Transform_Type& type)
	{
		auto data_r = data_r_i[0];
		auto data_i = data_r_i[1];
		const int n = data_r.size();
		assert(data_i.size() == n);

		switch (normalization)
		{
		case STANDARD:
			if (type == Transform_Type::INVERSE)
			{
				const double scale_factor = 1.0 / n;
				for (int i{}; i < n; i++)
				{
					data_r[i] *= scale_factor;
					data_i[i] *= scale_factor;
				}
			}
			break;
		case UNITARY:
			const double scale_factor = 1.0 / std::sqrt(n);
			for (int i{}; i < n; i++)
			{
				data_r[i] *= scale_factor;
				data_i[i] *= scale_factor;
			}
			break;
		default:
			// This should never occur in normal conditions. However this
			// clause has been added as a safeguard if other types of
			// normalizations are ever implemented, and the corresponding
			// test is forgotten in the present switch.
			throw std::exception("not implemented");
			//throw Math_Runtime_Exception.create_internal_error();
		}
	}

public:
	/**
		* Creates a instance of this class, with various normalization
		* conventions.
		*
		* @param normalization the type of normalization to be applied to the
		* transformed data
		*/
	Fast_Fourier_Transformer(const Dft_Normalization& normalization)
		:
		my_normalization{ normalization }
	{
	}

	/**
	 * Computes the standard transform of the specified complex data. The
	 * computation is done in place. The input data is laid out as follows
	 * <ul>
	 *   <li>{@code data_r_i[0][i]} is the real part of the {@code i}-th data point,</li>
	 *   <li>{@code data_r_i[1][i]} is the imaginary part of the {@code i}-th data point.</li>
	 * </ul>
	 *
	 * @param data_r_i the two dimensional array of real and imaginary parts of the data
	 * @param normalization the normalization to be applied to the transformed data
	 * @param type the type of transform (forward, inverse) to be performed
	 * @ if the number of rows of the specified
	 *   array is not two, or the array is not rectangular
	 * @ if the number of data points is not
	 *   a power of two
	 */
	static void transform_in_place(const std::vector<std::vector<double>>& data_r_i, const Dft_Normalization& normalization, const Transform_Type& type)
	{
		Math_Utils::check_dimension(data_r_i.size(), 2);
		auto data_r = data_r_i[0];
		auto data_i = data_r_i[1];
		Math_Arrays::check_equal_length(data_r, data_i);

		const int n = data_r.size();
		if (!Arithmetic_Utils::is_power_of_two(n))
		{
			throw std::exception("not implemented");
			//throw (Localized_FFT_Formats.NOT_POWER_OF_TWO_CONSIDER_PADDING, static_cast<int>(n));
		}

		if (n == 1)
		{
			return;
		}
		else if (n == 2)
		{
			const double src_r0 = data_r[0];
			const double src_i0 = data_i[0];
			const double src_r1 = data_r[1];
			const double src_i1 = data_i[1];

			// X_0 = x_0 + x_1
			data_r[0] = src_r0 + src_r1;
			data_i[0] = src_i0 + src_i1;
			// X_1 = x_0 - x_1
			data_r[1] = src_r0 - src_r1;
			data_i[1] = src_i0 - src_i1;

			normalize_transformed_data(data_r_i, normalization, type);
			return;
		}

		bit_reversal_shuffle2(data_r, data_i);

		// Do 4-term DFT.
		if (type == Transform_Type::INVERSE)
		{
			for (int i0 = 0; i0 < n; i0 += 4)
			{
				const int i1 = i0 + 1;
				const int i2 = i0 + 2;
				const int i3 = i0 + 3;

				const double src_r0 = data_r[i0];
				const double src_i0 = data_i[i0];
				const double src_r1 = data_r[i2];
				const double src_i1 = data_i[i2];
				const double src_r2 = data_r[i1];
				const double src_i2 = data_i[i1];
				const double src_r3 = data_r[i3];
				const double src_i3 = data_i[i3];

				// 4-term DFT
				// X_0 = x_0 + x_1 + x_2 + x_3
				data_r[i0] = src_r0 + src_r1 + src_r2 + src_r3;
				data_i[i0] = src_i0 + src_i1 + src_i2 + src_i3;
				// X_1 = x_0 - x_2 + j * (x_3 - x_1)
				data_r[i1] = src_r0 - src_r2 + (src_i3 - src_i1);
				data_i[i1] = src_i0 - src_i2 + (src_r1 - src_r3);
				// X_2 = x_0 - x_1 + x_2 - x_3
				data_r[i2] = src_r0 - src_r1 + src_r2 - src_r3;
				data_i[i2] = src_i0 - src_i1 + src_i2 - src_i3;
				// X_3 = x_0 - x_2 + j * (x_1 - x_3)
				data_r[i3] = src_r0 - src_r2 + (src_i1 - src_i3);
				data_i[i3] = src_i0 - src_i2 + (src_r3 - src_r1);
			}
		}
		else
		{
			for (int i0 = 0; i0 < n; i0 += 4)
			{
				const int i1 = i0 + 1;
				const int i2 = i0 + 2;
				const int i3 = i0 + 3;

				const double src_r0 = data_r[i0];
				const double src_i0 = data_i[i0];
				const double src_r1 = data_r[i2];
				const double src_i1 = data_i[i2];
				const double src_r2 = data_r[i1];
				const double src_i2 = data_i[i1];
				const double src_r3 = data_r[i3];
				const double src_i3 = data_i[i3];

				// 4-term DFT
				// X_0 = x_0 + x_1 + x_2 + x_3
				data_r[i0] = src_r0 + src_r1 + src_r2 + src_r3;
				data_i[i0] = src_i0 + src_i1 + src_i2 + src_i3;
				// X_1 = x_0 - x_2 + j * (x_3 - x_1)
				data_r[i1] = src_r0 - src_r2 + (src_i1 - src_i3);
				data_i[i1] = src_i0 - src_i2 + (src_r3 - src_r1);
				// X_2 = x_0 - x_1 + x_2 - x_3
				data_r[i2] = src_r0 - src_r1 + src_r2 - src_r3;
				data_i[i2] = src_i0 - src_i1 + src_i2 - src_i3;
				// X_3 = x_0 - x_2 + j * (x_1 - x_3)
				data_r[i3] = src_r0 - src_r2 + (src_i3 - src_i1);
				data_i[i3] = src_i0 - src_i2 + (src_r1 - src_r3);
			}
		}

		int last_n0{ 4 };
		int last_log_n0{ 2 };
		while (last_n0 < n)
		{
			int n0 = last_n0 << 1;
			int log_n0 = last_log_n0 + 1;
			double wSubN0R = W_SUB_N_R[log_n0];
			double wSubN0I = W_SUB_N_I[log_n0];
			if (type == Transform_Type::INVERSE)
			{
				wSubN0I = -wSubN0I;
			}

			// Combine even/odd transforms of size last_n0 into a transform of size N0 (last_n0 * 2).
			for (int dest_even_start_index{}; dest_even_start_index < n; dest_even_start_index += n0)
			{
				int dest_odd_start_index = dest_even_start_index + last_n0;

				double w_sub_n0_to_r_r = 1;
				double w_sub_n0_to_r_i = 0;

				for (int r{}; r < last_n0; r++)
				{
					double grR = data_r[dest_even_start_index + r];
					double grI = data_i[dest_even_start_index + r];
					double hr_r = data_r[dest_odd_start_index + r];
					double hr_i = data_i[dest_odd_start_index + r];

					// dest[dest_even_start_index + r] = Gr + WsubN0ToR * Hr
					data_r[dest_even_start_index + r] = grR + w_sub_n0_to_r_r * hr_r - w_sub_n0_to_r_i * hr_i;
					data_i[dest_even_start_index + r] = grI + w_sub_n0_to_r_r * hr_i + w_sub_n0_to_r_i * hr_r;
					// dest[dest_odd_start_index + r] = Gr - WsubN0ToR * Hr
					data_r[dest_odd_start_index + r] = grR - (w_sub_n0_to_r_r * hr_r - w_sub_n0_to_r_i * hr_i);
					data_i[dest_odd_start_index + r] = grI - (w_sub_n0_to_r_r * hr_i + w_sub_n0_to_r_i * hr_r);

					// WsubN0ToR *= WsubN0R
					double next_wsub_n0_to_r_r = w_sub_n0_to_r_r * wSubN0R - w_sub_n0_to_r_i * wSubN0I;
					double next_wsub_n0_to_r_i = w_sub_n0_to_r_r * wSubN0I + w_sub_n0_to_r_i * wSubN0R;
					w_sub_n0_to_r_r = next_wsub_n0_to_r_r;
					w_sub_n0_to_r_i = next_wsub_n0_to_r_i;
				}
			}

			last_n0 = n0;
			last_log_n0 = log_n0;
		}

		normalize_transformed_data(data_r_i, normalization, type);
	}

	/**
	 * Returns the (forward, inverse) transform of the specified real data set.
	 *
	 * @param f the real data array to be transformed
	 * @param type the type of transform (forward, inverse) to be performed
	 * @return the complex transformed array
	 * @ if the length of the data array is not a power of two
	 */
	std::vector<std::complex<double>> transform(const std::vector<double>& f, const Transform_Type& type)
	{
		const std::vector<std::vector<double>> data_r_i = { f, std::vector<double>(f.size()) };
		transform_in_place(data_r_i, my_normalization, type);
		return Trans_form_Utils::create_complex_array(data_r_i);
	}

	/**
	 * Returns the (forward, inverse) transform of the specified real function, * sampled on the specified interval.
	 *
	 * @param f the function to be sampled and transformed
	 * @param min the (inclusive) lower bound for the interval
	 * @param max the (exclusive) upper bound for the interval
	 * @param n the number of sample points
	 * @param type the type of transform (forward, inverse) to be performed
	 * @return the complex transformed array
	 * @org.hipparchus.exception.
	 *   if the lower bound is greater than, or equal to the upper bound
	 * @org.hipparchus.exception.
	 *   if the number of sample points {@code n} is negative
	 * @ if the number of sample points
	 *   {@code n} is not a power of two
	 */
	std::vector<std::complex<double>> transform(const Univariate_Function& f, const double& min, const double& max, const int& n, const Transform_Type& type)
	{
		const std::vector<double> data = Function_Utils::sample(f, min, max, n);
		return transform(data, type);
	}

	/**
	 * Returns the (forward, inverse) transform of the specified complex data set.
	 *
	 * @param f the complex data array to be transformed
	 * @param type the type of transform (forward, inverse) to be performed
	 * @return the complex transformed array
	 * @ if the length of the data array is not a power of two
	 */
	std::vector<std::complex<double>> transform(const std::vector<std::complex<double>>& f, const Transform_Type& type)
	{
		const std::vector<std::vector<double>> data_r_i = Transform_Utils::create_real_imaginary_array(f);

		transform_in_place(data_r_i, my_normalization, type);

		return Trans_form_Utils::create_complex_array(data_r_i);
	}
};