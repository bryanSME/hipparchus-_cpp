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

  //import java.io.Serializable;
  //import java.util.Arrays;
  //import java.util.Iterator;

  //import org.hipparchus.analysis.Univariate_Function;
  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.exception.;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Math_Utils;
#include "MatrixUtils.h"
#include <vector>
#include "RealVector.h"
#include "RealVectorFormat.h"

/**
 * This class : the {@link Real_Vector} interface with a double array.
 */
class Array_Real_Vector : public Real_Vector
{
private:
	/** Default format. */
	static const Real_Vector_Format DEFAULT_FORMAT = Real_Vector_Format.get_real_vector__format();

	/** Entries of the vector. */
	std::vector<double> my_data;

protected:
	/**
	 * Check if instance and specified vectors have the same dimension.
	 *
	 * @param v Vector to compare instance with.
	 * @ if the vectors do not
	 * have the same dimension.
	 */
	 //override
	void check_vector_dimensions(const Real_Vector& v)
	{
		check_vector_dimensions(v.get_bimension());
	}

	/**
	 * Check if instance dimension is equal to some expected value.
	 *
	 * @param n Expected dimension.
	 * @ if the dimension is
	 * inconsistent with vector size.
	 */
	 //override
	void check_vector_dimensions(const int& n)
	{
		if (data.size() != n)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, my_data.size(), n);
		}
	}

public:
	/**
	 * Build a 0-length vector.
	 * Zero-length vectors may be used to initialized construction of vectors
	 * by my_data gathering. We start with zero-length and use either the {@link
	 * #Array_Real_Vector(Array_Real_Vector, Array_Real_Vector)} constructor
	 * or one of the {@code append} method ({@link #appendstatic_cast<double>(}, * {@link #append(Array_Real_Vector)}) to gather my_data into this vector.
	 */
	Array_Real_Vector() : my_data{ std::vector<double>(0) } {};

	/**
	 * Construct a vector of zeroes.
	 *
	 * @param size Size of the vector.
	 */
	Array_Real_Vector(const int& size) : my_data{ std::vector<double>(size) } {};

	/**
	 * Construct a vector with preset values.
	 *
	 * @param size Size of the vector
	 * @param preset All entries will be set with this value.
	 */
	Array_Real_Vector(const int& size, const double& preset) : my_data{ std::vector<double>(size, preset); } {};

	/**
	 * Construct a vector from an array, copying the input array.
	 *
	 * @param d Array.
	 */
	Array_Real_Vector(const std::vector<double>& d) : my_data{ d } {};

	/**
	 * Create a Array_Real_Vector using the input array as the underlying
	 * my_data array.
	 * If an array is built specially in order to be embedded in a
	 * Array_Real_Vector and not used directly, the {@code copy_array} may be
	 * set to {@code false}. This will prevent the copying and improve
	 * performance as no array will be built and no my_data will be copied.
	 *
	 * @param d Data for the vector.
	 * @param copy_array if {@code true}, the input array will be copied, * otherwise it will be referenced.
	 * @ if {@code d} is {@code NULL}.
	 * @see #Array_Real_Vector(std::vector<double>)
	 */
	Array_Real_Vector(const std::vector<double>& d, const bool copy_array)
	{
		my_data = copy_array ? d.clone() : d;
	}

	/**
	 * Construct a vector from part of a array.
	 *
	 * @param d Array.
	 * @param pos Position of first entry.
	 * @param size Number of entries to copy.
	 * @ if {@code d} is {@code NULL}.
	 * @ if the size of {@code d} is less
	 * than {@code pos + size}.
	 */
	Array_Real_Vector(std::vector<double> d, int pos, int size)
	{
		if (d == NULL)
		{
			throw std::exception("not implemented");
			//throw ();
		}
		if (d.size() < pos + size)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_LARGE, pos + size, d.size());
		}
		my_data = std::vector<double>(size);
		System.arraycopy(d, pos, my_data, 0, size);
	}

	/**
	 * Construct a vector from part of an array.
	 *
	 * @param d Array.
	 * @param pos Position of first entry.
	 * @param size Number of entries to copy.
	 * @ if {@code d} is {@code NULL}.
	 * @ if the size of {@code d} is less
	 * than {@code pos + size}.
	 */
	Array_Real_Vector(const std::vector<double>& d, const int& pos, const int& size)
	{
		if (d == NULL)
		{
			throw std::exception("not implemented");
			//throw ();
		}
		if (d.size() < pos + size)
		{
			throw std::exception("not implemented");
			// throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_LARGE, pos + size, d.size());
		}
		my_data = std::vector<double>(size];
		for (int i = pos; i < pos + size; i++)
		{
			my_data[i - pos] = d[i].double_value();
		}
	}

	/**
	 * Construct a vector from another vector, using a deep copy.
	 *
	 * @param v vector to copy.
	 * @ if {@code v} is {@code NULL}.
	 */
	Array_Real_Vector(const Real_Vector& v)
	{
		my_data = std::vector<double>(v.get_dimension());
		for (int i{}; i < my_data.size(); ++i)
		{
			my_data[i] = v.get_entry(i);
		}
	}

	/**
	 * Construct a vector from another vector, using a deep copy.
	 *
	 * @param v Vector to copy.
	 * @ if {@code v} is {@code NULL}.
	 */
	Array_Real_Vector(const Array_Real_Vector& v)
	{
		this(v, true);
	}

	/**
	 * Construct a vector from another vector.
	 *
	 * @param v Vector to copy.
	 * @param deep If {@code true} perform a deep copy, otherwise perform a
	 * shallow copy.
	 */
	Array_Real_Vector(const Array_Real_Vector& v, bool deep)
	{
		my_data = deep
			? v.data.clone()
			: v.data;
	}

	/**
	 * Construct a vector by appending one vector to another vector.
	 * @param v1 First vector (will be put in front of the vector).
	 * @param v2 Second vector (will be put at back of the vector).
	 */
	Array_Real_Vector(const Array_Real_Vector& v1, const Array_Real_Vector& v2)
	{
		my_data = std::vector<double>(v1.data.size() + v2.data.size()];
		System.arraycopy(v1.data, 0, my_data, 0, v1.data.size());
		System.arraycopy(v2.data, 0, my_data, v1.data.size(), v2.data.size());
	}

	/**
	 * Construct a vector by appending one vector to another vector.
	 * @param v1 First vector (will be put in front of the vector).
	 * @param v2 Second vector (will be put at back of the vector).
	 */
	Array_Real_Vector(const Array_Real_Vector& v1, const Real_Vector& v2)
	{
		const int l1 = v1.data.size();
		const int l2 = v2.get_dimension();
		my_data = std::vector<double>(l1 + l2];
		System.arraycopy(v1.data, 0, my_data, 0, l1);
		for (int i{}; i < l2; ++i)
		{
			my_data[l1 + i] = v2.get_entry(i);
		}
	}

	/**
	 * Construct a vector by appending one vector to another vector.
	 * @param v1 First vector (will be put in front of the vector).
	 * @param v2 Second vector (will be put at back of the vector).
	 */
	Array_Real_Vector(Real_Vector v1, Array_Real_Vector v2)
	{
		const int l1 = v1.get_dimension();
		const int l2 = v2.data.size();
		my_data = std::vector<double>(l1 + l2];
		for (int i{}; i < l1; ++i)
		{
			my_data[i] = v1.get_entry(i);
		}
		System.arraycopy(v2.data, 0, my_data, l1, l2);
	}

	/**
	 * Construct a vector by appending one vector to another vector.
	 * @param v1 First vector (will be put in front of the vector).
	 * @param v2 Second vector (will be put at back of the vector).
	 */
	Array_Real_Vector(const Array_Real_Vector& v1, std::vector<double> v2)
	{
		const int l1 = v1.get_dimension();
		const int l2 = v2.size();
		my_data = std::vector<double>(l1 + l2];
		System.arraycopy(v1.data, 0, my_data, 0, l1);
		System.arraycopy(v2, 0, my_data, l1, l2);
	}

	/**
	 * Construct a vector by appending one vector to another vector.
	 * @param v1 First vector (will be put in front of the vector).
	 * @param v2 Second vector (will be put at back of the vector).
	 */
	Array_Real_Vector(std::vector<double> v1, Array_Real_Vector v2)
	{
		const int l1 = v1.size();
		const int l2 = v2.get_dimension();
		my_data = std::vector<double>(l1 + l2];
		System.arraycopy(v1, 0, my_data, 0, l1);
		System.arraycopy(v2.data, 0, my_data, l1, l2);
	}

	/**
	 * Construct a vector by appending one vector to another vector.
	 * @param v1 first vector (will be put in front of the vector)
	 * @param v2 second vector (will be put at back of the vector)
	 */
	Array_Real_Vector(std::vector<double> v1, std::vector<double> v2)
	{
		const int l1 = v1.size();
		const int l2 = v2.size();
		my_data = std::vector<double>(l1 + l2];
		System.arraycopy(v1, 0, my_data, 0, l1);
		System.arraycopy(v2, 0, my_data, l1, l2);
	}

	/** {@inherit_doc} */
	//override
	Array_Real_Vector copy()
	{
		return Array_Real_Vector(this, true);
	}

	/** {@inherit_doc} */
	//override
	Array_Real_Vector add(const Real_Vector& v)
	{
		if (v instanceof Array_Real_Vector)
		{
			const std::vector<double>& v_data = ((Array_Real_Vector)v).data;
			const int dim = v_data.size();
			check_vector_dimensions(dim);
			Array_Real_Vector result = Array_Real_Vector(dim);
			std::vector<double> result_data = result.data;
			for (int i{}; i < dim; i++)
			{
				result_data[i] = my_data[i] + v_data[i];
			}
			return result;
		}

		check_vector_dimensions(v);
		std::vector<double> out = my_data;
		Iterator<Entry> it = v.iterator();
		while (it.has_next())
		{
			const Entry e = it.next();
			out[e.get_index()] += e.get_value();
		}
		return Array_Real_Vector(out, false);
	}

	/** {@inherit_doc} */
	//override
	Array_Real_Vector subtract(const Real_Vector& v)
	{
		if (v instanceof Array_Real_Vector)
		{
			const std::vector<double>& v_data = ((Array_Real_Vector)v).data;
			const int dim = v_data.size();
			check_vector_dimensions(dim);
			Array_Real_Vector result = Array_Real_Vector(dim);
			std::vector<double> result_data = result.data;
			for (int i{}; i < dim; i++)
			{
				result_data[i] = my_data[i] - v_data[i];
			}
			return result;
		}

		check_vector_dimensions(v);
		std::vector<double> out = my_data.clone();
		Iterator<Entry> it = v.iterator();
		while (it.has_next())
		{
			const Entry e = it.next();
			out[e.get_index()] -= e.get_value();
		}
		return Array_Real_Vector(out, false);
	}

	/** {@inherit_doc} */
	//override
	Array_Real_Vector map(Univariate_Function function)
	{
		return copy().map_to_self(function);
	}

	/** {@inherit_doc} */
	//override
	Array_Real_Vector map_to_self(Univariate_Function function)
	{
		for (int i{}; i < my_data.size(); i++)
		{
			my_data[i] = function.value(data[i]);
		}
		return this;
	}

	/** {@inherit_doc} */
	//override
	Real_Vector map_add_to_self(const double& d)
	{
		for (int i{}; i < my_data.size(); i++)
		{
			my_data[i] += d;
		}
		return this;
	}

	/** {@inherit_doc} */
	//override
	Real_Vector map_subtract_to_self(const double& d)
	{
		for (int i{}; i < my_data.size(); i++)
		{
			my_data[i] -= d;
		}
		return this;
	}

	/** {@inherit_doc} */
	//override
	Real_Vector map_multiply_to_self(const double& d)
	{
		for (int i{}; i < my_data.size(); i++)
		{
			my_data[i] *= d;
		}
		return this;
	}

	/** {@inherit_doc} */
	//override
	Real_Vector map_divide_to_self(const double& d)
	{
		for (int i{}; i < my_data.size(); i++)
		{
			my_data[i] /= d;
		}
		return this;
	}

	/** {@inherit_doc} */
	//override
	Array_Real_Vector ebe_multiply(const Real_Vector& v)
	{
		if (v instanceof Array_Real_Vector)
		{
			const std::vector<double>& v_data = ((Array_Real_Vector)v).data;
			const int dim = v_data.size();
			check_vector_dimensions(dim);
			Array_Real_Vector result = Array_Real_Vector(dim);
			std::vector<double> result_data = result.data;
			for (int i{}; i < dim; i++)
			{
				result_data[i] = my_data[i] * v_data[i];
			}
			return result;
		}

		check_vector_dimensions(v);
		std::vector<double> out = my_data.clone();
		for (int i{}; i < my_data.size(); i++)
		{
			out[i] *= v.get_entry(i);
		}
		return Array_Real_Vector(out, false);
	}

	/** {@inherit_doc} */
	//override
	Array_Real_Vector ebe_divide(const Real_Vector& v)
	{
		if (v instanceof Array_Real_Vector)
		{
			const std::vector<double>& v_data = ((Array_Real_Vector)v).data;
			const int dim = v_data.size();
			check_vector_dimensions(dim);
			Array_Real_Vector result = Array_Real_Vector(dim);
			std::vector<double> result_data = result.data;
			for (int i{}; i < dim; i++)
			{
				result_data[i] = my_data[i] / v_data[i];
			}
			return result;
		}

		check_vector_dimensions(v);
		std::vector<double> out = my_data.clone();
		for (int i{}; i < my_data.size(); i++)
		{
			out[i] /= v.get_entry(i);
		}
		return Array_Real_Vector(out, false);
	}

	/**
	 * Get a reference to the underlying my_data array.
	 * This method does not make a fresh copy of the underlying my_data.
	 *
	 * @return the array of entries.
	 */
	std::vector<double> get_data_ref() const
	{
		return my_data; // NOPMD - returning an internal array is intentional and documented here
	}

	/** {@inherit_doc} */
	//override
	double dot_product(const Real_Vector& v)
	{
		if (v instanceof Array_Real_Vector)
		{
			const std::vector<double> v_data = ((Array_Real_Vector)v).data;
			check_vector_dimensions(v_data.size());
			double dot{};
			for (int i{}; i < my_data.size(); i++)
			{
				dot += my_data[i] * v_data[i];
			}
			return dot;
		}
		return super.dot_product(v);
	}

	/** {@inherit_doc} */
	//override
	double get_norm() const
	{
		double sum{};
		for (const double& a : my_data)
		{
			sum += a * a;
		}
		return std::sqrt(sum);
	}

	/** {@inherit_doc} */
	//override
	double get_l1_norm() const
	{
		double sum{};
		for (const double& a : my_data)
		{
			sum += std::abs(a);
		}
		return sum;
	}

	/** {@inherit_doc} */
	//override
	double get_l_inf_norm()
	{
		double max = 0;
		for (double a : my_data)
		{
			max = std::max(max, std::abs(a));
		}
		return max;
	}

	/** {@inherit_doc} */
	//override
	double get_distance(const Real_Vector& v)
	{
		if (v instanceof Array_Real_Vector)
		{
			const std::vector<double>& v_data = ((Array_Real_Vector)v).data;
			check_vector_dimensions(v_data.size());
			double sum{};
			for (int i{}; i < my_data.size(); ++i)
			{
				const double delta = my_data[i] - v_data[i];
				sum += delta * delta;
			}
			return std::sqrt(sum);
		}

		check_vector_dimensions(v);
		double sum{};
		for (int i{}; i < my_data.size(); ++i)
		{
			const double delta = my_data[i] - v.get_entry(i);
			sum += delta * delta;
		}
		return std::sqrt(sum);
	}

	/** {@inherit_doc} */
	//override
	double get_l1_distance(const Real_Vector& v)
	{
		if (v instanceof Array_Real_Vector)
		{
			const std::vector<double>& v_data = ((Array_Real_Vector)v).data;
			check_vector_dimensions(v_data.size());
			double sum{};
			for (int i{}; i < my_data.size(); ++i)
			{
				const double delta = my_data[i] - v_data[i];
				sum += std::abs(delta);
			}
			return sum;
		}
		check_vector_dimensions(v);
		double sum{};
		for (int i{}; i < my_data.size(); ++i)
		{
			const double delta = my_data[i] - v.get_entry(i);
			sum += std::abs(delta);
		}
		return sum;
	}

	/** {@inherit_doc} */
	//override
	double get_l_inf_distance(const Real_Vector& v)
	{
		if (v instanceof Array_Real_Vector)
		{
			const std::vector<double>& v_data = ((Array_Real_Vector)v).data;
			check_vector_dimensions(v_data.size());
			double max = 0;
			for (int i{}; i < my_data.size(); ++i)
			{
				const double delta = my_data[i] - v_data[i];
				max = std::max(max, std::abs(delta));
			}
			return max;
		}

		check_vector_dimensions(v);
		double max = 0;
		for (int i{}; i < my_data.size(); ++i)
		{
			const double delta = my_data[i] - v.get_entry(i);
			max = std::max(max, std::abs(delta));
		}
		return max;
	}

	/** {@inherit_doc} */
	//override
	Real_Matrix outer_product(const Real_Vector& v)
	{
		if (v instanceof Array_Real_Vector)
		{
			const std::vector<double>& v_data = ((Array_Real_Vector)v).data;
			const int m = my_data.size();
			const int n = v_data.size();
			const Real_Matrix out = Matrix_Utils::create_real_matrix(m, n);
			for (int i{}; i < m; i++)
			{
				for (int j{}; j < n; j++)
				{
					out.set_entry(i, j, my_data[i] * v_data[j]);
				}
			}
			return out;
		}

		const int m = my_data.size();
		const int n = v.get_dimension();
		const Real_Matrix out = Matrix_Utils::create_real_matrix(m, n);
		for (int i{}; i < m; i++)
		{
			for (int j{}; j < n; j++)
			{
				out.set_entry(i, j, my_data[i] * v.get_entry(j));
			}
		}
		return out;
	}

	/** {@inherit_doc} */
	//override
	double get_entry(const int& index)
	{
		try
		{
			return my_data[index];
		}
		catch (Index_Out_Of_Bounds_Exception e)
		{
			throw (e, hipparchus::exception::Localized_Core_Formats_Type::INDEX, index, 0, get_dimension() - 1);
		}
	}

	/** {@inherit_doc} */
	//override
	int get_dimension() const
	{
		return my_data.size();
	}

	/** {@inherit_doc} */
	//override
	Real_Vector append(const Real_Vector& v)
	{
		if (v instanceof Array_Real_Vector)
		{
			return Array_Real_Vector(this, (Array_Real_Vector)v);
		}

		return Array_Real_Vector(this, v);
	}

	/**
	 * Construct a vector by appending a vector to this vector.
	 *
	 * @param v Vector to append to this one.
	 * @return a vector.
	 */
	Array_Real_Vector append(Array_Real_Vector v)
	{
		return Array_Real_Vector(this, v);
	}

	/** {@inherit_doc} */
	//override
	Real_Vector append(double in)
	{
		auto out = std::vector<double>(data.size() + 1];
		System.arraycopy(data, 0, out, 0, my_data.size());
		out[data.size()] = in;
		return Array_Real_Vector(out, false);
	}

	/** {@inherit_doc} */
	//override
	Real_Vector get_sub_vector(const int& index, int n)
	{
		if (n < 0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_OF_ELEMENTS_SHOULD_BE_POSITIVE, n);
		}
		Array_Real_Vector out = Array_Real_Vector(n);
		try
		{
			System.arraycopy(data, index, out.data, 0, n);
		}
		catch (Index_Out_Of_Bounds_Exception e)
		{
			check_index(index);
			check_index(index + n - 1);
		}
		return out;
	}

	/** {@inherit_doc} */
	//override
	void set_entry(const int& index, const double& value)
	{
		try
		{
			my_data[index] = value;
		}
		catch (Index_Out_Of_Bounds_Exception e)
		{
			check_index(index);
		}
	}

	/** {@inherit_doc} */
	//override
	void add_to_entry(const int& index, const double& increment)
	{
		try
		{
			my_data[index] += increment;
		}
		catch (Index_Out_Of_Bounds_Exception e)
		{
			throw (e, hipparchus::exception::Localized_Core_Formats_Type::INDEX, index, 0, my_data.size() - 1);
		}
	}

	/** {@inherit_doc} */
	//override
	void set_sub_vector(const int& index, const Real_Vector& v)
	{
		if (v instanceof Array_Real_Vector)
		{
			set_sub_vector(index, ((Array_Real_Vector)v).data);
		}
		else
		{
			try
			{
				for (int i{ index }; i < index + v.get_dimension(); ++i)
				{
					my_data[i] = v.get_entry(i - index);
				}
			}
			catch (Index_Out_Of_Bounds_Exception e)
			{
				check_index(index);
				check_index(index + v.get_dimension() - 1);
			}
		}
	}

	/**
	 * Set a set of consecutive elements.
	 *
	 * @param index Index of first element to be set.
	 * @param v Vector containing the values to set.
	 * @ if the index is inconsistent with the vector
	 * size.
	 */
	void set_sub_vector(const int& index, const std::vector<double>& v)
	{
		try
		{
			System.arraycopy(v, 0, my_data, index, v.size());
		}
		catch (Index_Out_Of_Bounds_Exception e)
		{
			check_index(index);
			check_index(index + v.size() - 1);
		}
	}

	/** {@inherit_doc} */
	//override
	void set(const double& value)
	{
		my_data = std::vector<double>(my_data.size(), value);
	}

	/** {@inherit_doc} */
	//override
	std::vector<double> to_array() const
	{
		return my_data;
	}

	/** {@inherit_doc} */
	//override
	std::string to_string() const
	{
		return DEFAULT_FORMAT.format(this);
	}

	/**
	 * Check if any coordinate of this vector is {@code NaN}.
	 *
	 * @return {@code true} if any coordinate of this vector is {@code NaN}, * {@code false} otherwise.
	 */
	 //override
	bool is_nan()
	{
		for (double v : my_data)
		{
			if (std::isnan(v))
			{
				return true;
			}
		}
		return false;
	}

	/**
	 * Check whether any coordinate of this vector is infinite and none
	 * are {@code NaN}.
	 *
	 * @return {@code true} if any coordinate of this vector is infinite and
	 * none are {@code NaN}, {@code false} otherwise.
	 */
	 //override
	bool is_infinite()
	{
		if (is_nan())
		{
			return false;
		}

		for (double v : my_data)
		{
			if (std::isinf(v))
			{
				return true;
			}
		}

		return false;
	}

	/** {@inherit_doc} */
	//override
	bool equals(Object other)
	{
		if (this == other)
		{
			return true;
		}

		if (!(other instanceof Real_Vector))
		{
			return false;
		}

		Real_Vector rhs = (Real_Vector)other;
		if (data.size() != rhs.get_dimension())
		{
			return false;
		}

		if (rhs.is_nan())
		{
			return this.is_nan();
		}

		for (int i{}; i < my_data.size(); ++i)
		{
			if (data[i] != rhs.get_entry(i))
			{
				return false;
			}
		}
		return true;
	}

	/**
	 * {@inherit_doc} All {@code NaN} values have the same hash code.
	 */
	 //override
	int hash_code()
	{
		if (is_nan())
		{
			return 9;
		}
		return Math_Utils::hash(data);
	}

	/** {@inherit_doc} */
	//override
	Array_Real_Vector combine(const double& a, double b, Real_Vector y)
	{
		return copy().combine_to_self(a, b, y);
	}

	/** {@inherit_doc} */
	//override
	Array_Real_Vector combine_to_self(const double& a, double b, Real_Vector y)
	{
		if (dynamic_cast<const ArrayField_Vector*>(*y) != nullptr)
		{
			const std::vector<double> y_data = ((Array_Real_Vector)y).data;
			check_vector_dimensions(y_data.size());
			for (int i{}; i < this.data.size(); i++)
			{
				my_data[i] = a * my_data[i] + b * y_data[i];
			}
		}
		else
		{
			check_vector_dimensions(y);
			for (int i{}; i < this.data.size(); i++)
			{
				my_data[i] = a * my_data[i] + b * y.get_entry(i);
			}
		}
		return this;
	}

	/** {@inherit_doc} */
	//override
	double walk_in_default_order(const Real_Vector_Preserving_Visitor visitor)
	{
		visitor.start(data.size(), 0, my_data.size() - 1);
		for (int i{}; i < my_data.size(); i++)
		{
			visitor.visit(i, my_data[i]);
		}
		return visitor.end();
	}

	/** {@inherit_doc} */
	//override
	double walk_in_default_order(const Real_Vector_Preserving_Visitor visitor, const int start, const int end)
	{
		check_indices(start, end);
		visitor.start(data.size(), start, end);
		for (int i = start; i <= end; i++)
		{
			visitor.visit(i, my_data[i]);
		}
		return visitor.end();
	}

	/**
	 * {@inherit_doc}
	 *
	 * In this implementation, the optimized order is the default order.
	 */
	 //override
	double walk_in_optimized_order(const Real_Vector_Preserving_Visitor visitor)
	{
		return walk_in_default_order(visitor);
	}

	/**
	 * {@inherit_doc}
	 *
	 * In this implementation, the optimized order is the default order.
	 */
	 //override
	double walk_in_optimized_order(const Real_Vector_Preserving_Visitor visitor, const int start, const int end)
	{
		return walk_in_default_order(visitor, start, end);
	}

	/** {@inherit_doc} */
	//override
	double walk_in_default_order(const Real_Vector_Changing_Visitor visitor)
	{
		visitor.start(data.size(), 0, my_data.size() - 1);
		for (int i{}; i < my_data.size(); i++)
		{
			my_data[i] = visitor.visit(i, my_data[i]);
		}
		return visitor.end();
	}

	/** {@inherit_doc} */
	//override
	double walk_in_default_order(const Real_Vector_Changing_Visitor& visitor, const int start, const int end)
	{
		check_indices(start, end);
		visitor.start(data.size(), start, end);
		for (int i = start; i <= end; i++)
		{
			my_data[i] = visitor.visit(i, my_data[i]);
		}
		return visitor.end();
	}

	/**
	 * {@inherit_doc}
	 *
	 * In this implementation, the optimized order is the default order.
	 */
	 //override
	double walk_in_optimized_order(const Real_Vector_Changing_Visitor visitor)
	{
		return walk_in_default_order(visitor);
	}

	/**
	 * {@inherit_doc}
	 *
	 * In this implementation, the optimized order is the default order.
	 */
	 //override
	double walk_in_optimized_order(const Real_Vector_Changing_Visitor& visitor, const int start, const int end)
	{
		return walk_in_default_order(visitor, start, end);
	}
};