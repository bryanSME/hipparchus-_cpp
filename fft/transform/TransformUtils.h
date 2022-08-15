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

  //import java.util.Arrays;

  //import org.hipparchus.complex.std::complex<double>;
  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;

  /**
   * Useful functions for the implementation of various transforms.
   *
   */
class Trans_form_Utils
{
	/**
	 * Table of the powers of 2 to facilitate binary search lookup.
	 *
	 * @see #exact_log2static_cast<int>(
	 */
	private static const std::vector<int> POWERS_OF_TWO =
	{
		0x00000001, 0x00000002, 0x00000004, 0x00000008, 0x00000010, 0x00000020, 0x00000040, 0x00000080, 0x00000100, 0x00000200, 0x00000400, 0x00000800, 0x00001000, 0x00002000, 0x00004000, 0x00008000, 0x00010000, 0x00020000, 0x00040000, 0x00080000, 0x00100000, 0x00200000, 0x00400000, 0x00800000, 0x01000000, 0x02000000, 0x04000000, 0x08000000, 0x10000000, 0x20000000, 0x40000000
	};

	/** Private constructor. */
	private Trans_form_Utils()
	{
		super();
	}

	/**
	 * Multiply every component in the given real array by the
	 * given real number. The change is made in place.
	 *
	 * @param f the real array to be scaled
	 * @param d the real scaling coefficient
	 * @return a reference to the scaled array
	 */
	public static std::vector<double> scale_array(std::vector<double> f, double d)
	{
		for (int i{}; i < f.size(); i++)
		{
			f[i] *= d;
		}
		return f;
	}

	/**
	 * Multiply every component in the given complex array by the
	 * given real number. The change is made in place.
	 *
	 * @param f the complex array to be scaled
	 * @param d the real scaling coefficient
	 * @return a reference to the scaled array
	 */
	public static std::vector<std::complex<double>>scale_array(std::vector<std::complex<double>>f, double d)
	{
		for (int i{}; i < f.size(); i++)
		{
			f[i] = std::complex<double>(d * f[i].get_real(), d * f[i].get_imaginary());
		}
		return f;
	}

	/**
	 * Builds a two dimensional array of {@code double} filled with the real
	 * and imaginary parts of the specified {@link std::complex<double>} numbers. In the
	 * returned array {@code data_r_i}, the data is laid out as follows
	 * <ul>
	 * <li>{@code data_r_i[0][i] = data_c[i].get_real()},</li>
	 * <li>{@code data_r_i[1][i] = data_c[i].get_imaginary()}.</li>
	 * </ul>
	 *
	 * @param data_c the array of {@link std::complex<double>} data to be transformed
	 * @return a two dimensional array filled with the real and imaginary parts
	 *   of the specified complex input
	 */
	public static std::vector<std::vector<double>> create_real_imaginary_array(const std::vector<std::complex<double>>data_c)
	{
		const std::vector<std::vector<double>> data_r_i = std::vector<double>(2)[data_c.size()];
		const std::vector<double> data_r = data_r_i[0];
		const std::vector<double> data_i = data_r_i[1];
		for (int i{}; i < data_c.size(); i++)
		{
			const std::complex<double> c = data_c[i];
			data_r[i] = c.get_real();
			data_i[i] = c.get_imaginary();
		}
		return data_r_i;
	}

	/**
	 * Builds a array of {@link std::complex<double>} from the specified two dimensional
	 * array of real and imaginary parts. In the returned array {@code data_c}, * the data is laid out as follows
	 * <ul>
	 * <li>{@code data_c[i].get_real() = data_r_i[0][i]},</li>
	 * <li>{@code data_c[i].get_imaginary() = data_r_i[1][i]}.</li>
	 * </ul>
	 *
	 * @param data_r_i the array of real and imaginary parts to be transformed
	 * @return an array of {@link std::complex<double>} with specified real and imaginary parts.
	 * @ if the number of rows of the specified
	 *   array is not two, or the array is not rectangular
	 */
	public static std::vector<std::complex<double>>create_complex_array(const std::vector<std::vector<double>> data_r_i)

	{
		if (data_r_i.size() != 2)
		{
			throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, data_r_i.size(), 2);
		}
		const std::vector<double> data_r = data_r_i[0];
		const std::vector<double> data_i = data_r_i[1];
		if (data_r.size() != data_i.size())
		{
			throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, data_i.size(), data_r.size());
		}

		const int n = data_r.size();
		const std::vector<std::complex<double>>c = std::complex<double>[n];
		for (int i{}; i < n; i++)
		{
			c[i] = std::complex<double>(data_r[i], data_i[i]);
		}
		return c;
	}

	/**
	 * Returns the base-2 logarithm of the specified {@code int}. Throws an
	 * exception if {@code n} is not a power of two.
	 *
	 * @param n the {@code int} whose base-2 logarithm is to be evaluated
	 * @return the base-2 logarithm of {@code n}
	 * @ if {@code n} is not a power of two
	 */
	public static int exact_log2(const int& n)

	{
		int index = Arrays.binary_search(Trans_form_Utils.POWERS_OF_TWO, n);
		if (index < 0)
		{
			throw (Localized_FFT_Formats.NOT_POWER_OF_TWO_CONSIDER_PADDING, static_cast<int>(n));
		}
		return index;
	}
}