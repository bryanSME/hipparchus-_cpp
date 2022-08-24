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
  //package org.hipparchus.analysis.differentiation;

#include <cmath>
#include "../../util/FieldSinCos.h"
//import java.util.Array_list;
//import java.util.Arrays;
//import java.util.List;
//import java.util.concurrent.atomic.Atomic_Reference;

//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.Field;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Runtime_Exception;
//import org.hipparchus.util.Combinatorics_Utils;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Field_Sin_Cos;
//import org.hipparchus.util.Field_Sinh_Cosh;
//import org.hipparchus.util.Math_Arrays;
//import org.hipparchus.util.Math_Utils;
//import org.hipparchus.util.Sin_Cos;
//import org.hipparchus.util.Sinh_Cosh;
#include <vector>

/** Class holding "compiled" computation rules for derivative structures.
 * <p>This class : the computation rules described in Dan Kalman's paper <a
 * href="http://www1.american.edu/cas/mathstat/People/kalman/pdffiles/mmgautodiff.pdf">Doubly
 * Recursive Multivariate Automatic Differentiation</a>, Mathematics Magazine, vol. 75, * no. 3, June 2002. However, in order to avoid performances bottlenecks, the recursive
 * rules are "compiled" once in an unfold form. This class does this recursion unrolling
 * and stores the computation rules as simple loops with pre-computed indirection arrays.</p>
 * <p>
 * This class maps all derivative computation into single dimension arrays that hold the
 * value and partial derivatives. The class does not hold these arrays, which remains under
 * the responsibility of the caller. For each combination of number of free parameters and
 * derivation order, only one compiler is necessary, and this compiler will be used to
 * perform computations on all arrays provided to it, which can represent hundreds or
 * thousands of different parameters kept together with all their partial derivatives.
 * </p>
 * <p>
 * The arrays on which compilers operate contain only the partial derivatives together
 * with the 0<sup>th</sup> derivative, i.e. the value. The partial derivatives are stored in
 * a compiler-specific order, which can be retrieved using methods {@link
 * #get_partial_derivative_index(int...) get_partial_derivative_index} and {@link
 * #get_partial_derivative_ordersstatic_cast<int>(}. The value is guaranteed to be stored as the first element
 * (i.e. the {@link #get_partial_derivative_index(int...) get_partial_derivative_index} method returns
 * 0 when called with 0 for all derivation orders and {@link #get_partial_derivative_ordersstatic_cast<int>(
 * get_partial_derivative_orders} returns an array filled with 0 when called with 0 as the index).
 * </p>
 * <p>
 * Note that the ordering changes with number of parameters and derivation order. For example
 * given 2 parameters x and y, df/dy is stored at index 2 when derivation order is set to 1 (in
 * this case the array has three elements: f, df/dx and df/dy). If derivation order is set to
 * 2, then df/dy will be stored at index 3 (in this case the array has six elements: f, df/dx, * d\xc2\xb2f/dxdx, df/dy, d\xc2\xb2f/dxdy and d\xc2\xb2f/dydy).
 * </p>
 * <p>
 * Given this structure, users can perform some simple operations like adding, subtracting
 * or multiplying constants and negating the elements by themselves, knowing if they want to
 * mutate their array or create a array. These simple operations are not provided by
 * the compiler. The compiler provides only the more complex operations between several arrays.
 * </p>
 * <p>This class is mainly used as the engine for scalar variable {@link Derivative_Structure}.
 * It can also be used directly to hold several variables in arrays for more complex data
 * structures. User can for example store a vector of n variables depending on three x, y
 * and z free parameters in one array as follows:</p> <pre>
 *   // parameter 0 is x, parameter 1 is y, parameter 2 is z
 *   int parameters = 3;
 *   DS_Compiler compiler = DS_Compiler.get_compiler(parameters, order);
 *   int size = compiler.get_size();
 *
 *   // pack all elements in a single array
 *   std::vector<double> array = std::vector<double>(n * size];
 *   for (int i{}; i &lt; n; ++i)
 {
 *
 *     // we know value is guaranteed to be the first element
 *     array[i * size] = v[i];
 *
 *     // we don't know where first derivatives are stored, so we ask the compiler
 *     array[i * size + compiler.get_partial_derivative_index(1, 0, 0) = dvOnDx[i][0];
 *     array[i * size + compiler.get_partial_derivative_index(0, 1, 0) = dvOnDy[i][0];
 *     array[i * size + compiler.get_partial_derivative_index(0, 0, 1) = dvOnDz[i][0];
 *
 *     // we let all higher order derivatives set to 0
 *
 *   }
 * </pre>
 * <p>Then in another function, user can perform some operations on all elements stored
 * in the single array, such as a simple product of all variables:</p> <pre>
 *   // compute the product of all elements
 *   std::vector<double> product = std::vector<double>(size];
 *   prod[0] = 1.0;
 *   for (int i{}; i &lt; n; ++i)
 {
 *     std::vector<double> tmp = product.clone();
 *     compiler.multiply(tmp, 0, array, i * size, product, 0);
 *   }
 *
 *   // value
 *   double p = product[0];
 *
 *   // first derivatives
 *   double dPdX = product[compiler.get_partial_derivative_index(1, 0, 0)];
 *   double dPdY = product[compiler.get_partial_derivative_index(0, 1, 0)];
 *   double dPdZ = product[compiler.get_partial_derivative_index(0, 0, 1)];
 *
 *   // cross derivatives (assuming order was at least 2)
 *   double dPdXdX = product[compiler.get_partial_derivative_index(2, 0, 0)];
 *   double dPdXdY = product[compiler.get_partial_derivative_index(1, 1, 0)];
 *   double dPdXdZ = product[compiler.get_partial_derivative_index(1, 0, 1)];
 *   double dPdYdY = product[compiler.get_partial_derivative_index(0, 2, 0)];
 *   double dPdYdZ = product[compiler.get_partial_derivative_index(0, 1, 1)];
 *   double dPdZdZ = product[compiler.get_partial_derivative_index(0, 0, 2)];
 * </pre>
 * @see Derivative_Structure
 * @see Field_Derivative_Structure
 */
class DS_Compiler
{
private:
	/** Array of all compilers created so far. */
	static Atomic_Reference<std::vector<std::vector<DS_Compiler>>> compilers;

	/** Number of free parameters. */
	const int parameters;

	/** Derivation order. */
	const int order;

	/** Number of partial derivatives (including the single 0 order derivative element). */
	const std::vector<std::vector<int>> my_sizes;

	/** Indirection array for partial derivatives. */
	const std::vector<std::vector<int>> my_derivatives_indirection;

	/** Indirection array of the lower derivative elements. */
	const std::vector<int> my_lower_indirection;

	/** Indirection arrays for multiplication. */
	const std::vector<std::vector<int>>[] my_mult_indirection;

	/** Indirection arrays for function composition. */
	const std::vector<std::vector<int>>[] my_comp_indirection;

	/** Private constructor, reserved for the factory method {@link #get_compiler(int, int)}.
	 * @param parameters number of free parameters
	 * @param order derivation order
	 * @param value_compiler compiler for the value part
	 * @param derivative_compiler compiler for the derivative part
	 * @ if order is too large
	 */
	DS_Compiler(const int& parameters, const int& order, const DS_Compiler& value_compiler, const DS_Compiler& derivative_compiler)
	{
		my_parameters = parameters;
		my_order = order;
		my_sizes = compile_sizes(parameters, order, value_compiler);
		my_derivatives_indirection = compile_derivatives_indirection(parameters, order, value_compiler, derivative_compiler);
		my_lower_indirection = compile_lower_indirection(parameters, order, value_compiler, derivative_compiler);
		my_mult_indirection = compile_multiplication_indirection(parameters, order, value_compiler, derivative_compiler, lower_indirection);
		my_comp_indirection = compile_composition_indirection(parameters, order, value_compiler, derivative_compiler, sizes, derivatives_indirection);
	}

	/** Compile the sizes array.
	 * @param parameters number of free parameters
	 * @param order derivation order
	 * @param value_compiler compiler for the value part
	 * @return sizes array
	 */
	static std::vector<std::vector<int>> compile_sizes(const int& parameters, const int& order, const DS_Compiler& value_compiler)
	{
		auto sizes = std::vector < std::vector<int>(parameters + 1, std::vector<int>(order + 1));
		if (parameters == 0)
		{
			Arrays.fill(sizes[0], 1);
		}
		else
		{
			System.arraycopy(value_compiler.sizes, 0, sizes, 0, parameters);
			sizes[parameters][0] = 1;
			for (int i{}; i < order; ++i)
			{
				sizes[parameters][i + 1] = sizes[parameters][i] + sizes[parameters - 1][i + 1];
			}
		}

		return sizes;
	}

	/** Compile the derivatives indirection array.
	 * @param parameters number of free parameters
	 * @param order derivation order
	 * @param value_compiler compiler for the value part
	 * @param derivative_compiler compiler for the derivative part
	 * @return derivatives indirection array
	 */
	static std::vector<std::vector<int>> compile_derivatives_indirection(const int& parameters, const int& order, const DS_Compiler& value_compiler, const DS_Compiler& derivative_compiler)
	{
		if (parameters == 0 || order == 0)
		{
			return std::vector < std::vector<int>(1, std::vector<int>(parameters));
		}

		const int v_size = value_compiler.derivatives_indirection.size();
		const int d_size = derivative_compiler.derivatives_indirection.size();
		auto derivatives_indirection = std::vector < std::vector<int>(v_size + d_size, std::vector<int>(parameters));

		// set up the indices for the value part
		for (int i{}; i < v_size; ++i)
		{
			// copy the first indices, the last one remaining set to 0
			System.arraycopy(value_compiler.derivatives_indirection[i], 0, derivatives_indirection[i], 0, parameters - 1);
		}

		// set up the indices for the derivative part
		for (int i{}; i < d_size; ++i)
		{
			// copy the indices
			System.arraycopy(derivative_compiler.derivatives_indirection[i], 0, derivatives_indirection[v_size + i], 0, parameters);

			// increment the derivation order for the last parameter
			derivatives_indirection[v_size + i][parameters - 1]++;
		}

		return derivatives_indirection;
	}

	/** Compile the lower derivatives indirection array.
	 * <p>
	 * This indirection array contains the indices of all elements
	 * except derivatives for last derivation order.
	 * </p>
	 * @param parameters number of free parameters
	 * @param order derivation order
	 * @param value_compiler compiler for the value part
	 * @param derivative_compiler compiler for the derivative part
	 * @return lower derivatives indirection array
	 */
	static std::vector<int> compile_lower_indirection(const int& parameters, const int& order, const DS_Compiler& value_compiler, const DS_Compiler& derivative_compiler)
	{
		if (parameters == 0 || order <= 1)
		{
			return std::vector<int> { 0 };
		}

		// this is an implementation of definition 6 in Dan Kalman's paper.
		const int v_size = value_compiler.lower_indirection.size();
		const int d_size = derivative_compiler.lower_indirection.size();
		const std::vector<int> lower_indirection = int[v_size + d_size];
		System.arraycopy(value_compiler.lower_indirection, 0, lower_indirection, 0, v_size);
		for (int i{}; i < d_size; ++i)
		{
			lower_indirection[v_size + i] = value_compiler.get_size() + derivative_compiler.lower_indirection[i];
		}

		return lower_indirection;
	}

	/** Compile the multiplication indirection array.
	 * <p>
	 * This indirection array contains the indices of all pairs of elements
	 * involved when computing a multiplication. This allows a straightforward
	 * loop-based multiplication (see {@link #multiply(std::vector<double>, int, std::vector<double>, int, std::vector<double>, int)}).
	 * </p>
	 * @param parameters number of free parameters
	 * @param order derivation order
	 * @param value_compiler compiler for the value part
	 * @param derivative_compiler compiler for the derivative part
	 * @param lower_indirection lower derivatives indirection array
	 * @return multiplication indirection array
	 */
	static std::vector<std::vector<int>> compile_multiplication_indirection(const int& parameters, const int& order, const DS_Compiler& value_compiler, const DS_Compiler& derivative_compiler, const std::vector<int> lower_indirection)
	{
		if (parameters == 0 || order == 0)
		{
			return std::vector<std::vector<int>>[] { { { 1, 0, 0 } } };
		}

		// this is an implementation of definition 3 in Dan Kalman's paper.
		const int v_size = value_compiler.mult_indirection.size();
		const int d_size = derivative_compiler.mult_indirection.size();
		const std::vector<std::vector<int>>[] mult_indirection = int[v_size + d_size][][];

		System.arraycopy(value_compiler.mult_indirection, 0, mult_indirection, 0, v_size);

		for (int i{}; i < d_size; ++i)
		{
			const std::vector<std::vector<int>> d_row = derivative_compiler.mult_indirection[i];
			auto row = std::vector<int>(d_row.size() * 2);
			for (int j{}; j < d_row.size(); ++j)
			{
				row.add(new std::vector<int>{ d_row[j][0], lower_indirection[d_row[j][1]], v_size + d_row[j][2] });
				row.add(new std::vector<int>{ d_row[j][0], v_size + d_row[j][1], lower_indirection[d_row[j][2]] });
			}

			// combine terms with similar derivation orders
			const List<std::vector<int>> combined = Array_list<>(row.size());
			for (int j{}; j < row.size(); ++j)
			{
				auto term_j = row.get(j);
				if (term_j[0] > 0)
				{
					for (int k{ j + 1 }; k < row.size(); ++k)
					{
						const std::vector<int> term_k = row.get(k);
						if (term_j[1] == term_k[1] && term_j[2] == term_k[2])
						{
							// combine term_j and term_k
							term_j[0] += term_k[0];
							// make sure we will skip term_k later on in the outer loop
							term_k[0] = 0;
						}
					}
					combined.add(term_j);
				}
			}

			mult_indirection[v_size + i] = combined.to_array(new int[0][]);
		}

		return mult_indirection;
	}

	/** Compile the function composition indirection array.
	 * <p>
	 * This indirection array contains the indices of all sets of elements
	 * involved when computing a composition. This allows a straightforward
	 * loop-based composition (see {@link #compose(std::vector<double>, int, std::vector<double>, std::vector<double>, int)}).
	 * </p>
	 * @param parameters number of free parameters
	 * @param order derivation order
	 * @param value_compiler compiler for the value part
	 * @param derivative_compiler compiler for the derivative part
	 * @param sizes sizes array
	 * @param derivatives_indirection derivatives indirection array
	 * @return multiplication indirection array
	 * @ if order is too large
	 */
	static std::vector<std::vector<int>>[] compile_composition_indirection(const int& parameters, const int& order, const DS_Compiler& value_compiler, const DS_Compiler& derivative_compiler, const std::vector<std::vector<int>>& sizes, const std::vector<std::vector<int>>& derivatives_indirection)
	{
		if (parameters == 0 || order == 0)
		{
			return std::vector<std::vector<int>>[] { { { 1, 0 } } };
		}

		const int v_size = value_compiler.comp_indirection.size();
		const int d_size = derivative_compiler.comp_indirection.size();
		auto comp_indirection = std::vector < std::vector < std::vector<int>(v_size + d_size);

		// the composition rules from the value part can be reused as is
		System.arraycopy(value_compiler.comp_indirection, 0, comp_indirection, 0, v_size);

		// the composition rules for the derivative part are deduced by
		// differentiation the rules from the underlying compiler once
		// with respect to the parameter this compiler handles and the
		// underlying one did not handle
		for (int i{}; i < d_size; ++i)
		{
			auto row = std::vector<int>();
			for (const auto& term : derivative_compiler.comp_indirection[i])
			{
				// handle term p * f_k(g(x)) * g_l1(x) * g_l2(x) * ... * g_lp(x)

				// derive the first factor in the term: f_k with respect to parameter
				auto derived_term_f = std::vector<int>(term.size() + 1);
				derived_term_f[0] = term[0];     // p
				derived_term_f[1] = term[1] + 1; // f_(k+1)
				std::vector<int> orders = int[parameters];
				orders[parameters - 1] = 1;
				derived_term_f[term.size()] = get_partial_derivative_index(parameters, order, sizes, orders);  // g_1
				for (int j{ 2 }; j < term.size(); ++j)
				{
					// convert the indices as the mapping for the current order
					// is different from the mapping with one less order
					derived_term_f[j] = convert_index(term[j], parameters, derivative_compiler.derivatives_indirection, parameters, order, sizes);
				}
				Arrays.sort(derived_term_f, 2, derived_term_f.size());
				row.add(derived_term_f);

				// derive the various g_l
				for (const int l{ 2 }; l < term.size(); ++l)
				{
					auto derived_term_g = std::vector<int>(term.size());
					derived_term_g[0] = term[0];
					derived_term_g[1] = term[1];
					for (int j = 2; j < term.size(); ++j)
					{
						// convert the indices as the mapping for the current order
						// is different from the mapping with one less order
						derived_term_g[j] = convert_index(term[j], parameters, derivative_compiler.derivatives_indirection, parameters, order, sizes);
						if (j == l)
						{
							// derive this term
							System.arraycopy(derivatives_indirection[derived_term_g[j]], 0, orders, 0, parameters);
							orders[parameters - 1]++;
							derived_term_g[j] = get_partial_derivative_index(parameters, order, sizes, orders);
						}
					}
					Arrays.sort(derived_term_g, 2, derived_term_g.size());
					row.add(derived_term_g);
				}
			}

			// combine terms with similar derivation orders
			const auto combined = std::vector<int>(row.size());
			for (int j{}; j < row.size(); ++j)
			{
				const auto term_j = row.get(j);
				if (term_j[0] > 0)
				{
					for (int k = j + 1; k < row.size(); ++k)
					{
						const auto term_k = row.get(k);
						bool equals = term_j.size() == term_k.size();
						for (int l{ 1 }; equals && l < term_j.size(); ++l)
						{
							equals &= term_j[l] == term_k[l];
						}
						if (equals)
						{
							// combine term_j and term_k
							term_j[0] += term_k[0];
							// make sure we will skip term_k later on in the outer loop
							term_k[0] = 0;
						}
					}
					combined.add(term_j);
				}
			}

			comp_indirection[v_size + i] = combined.to_array(new int[0][]);
		}

		return comp_indirection;
	}

	/** Get the index of a partial derivative in an array.
	 * @param parameters number of free parameters
	 * @param order derivation order
	 * @param sizes sizes array
	 * @param orders derivation orders with respect to each parameter
	 * (the length of this array must match the number of parameters)
	 * @return index of the partial derivative
	 * @exception  if sum of derivation orders is larger
	 * than the instance limits
	 */
	static int get_partial_derivative_index(const int& parameters, const int& order, const std::vector<std::vector<int>>& sizes, const int& ... orders)
	{
		// the value is obtained by diving into the recursive Dan Kalman's structure
		// this is theorem 2 of his paper, with recursion replaced by iteration
		int index{};
		int m{ order };
		int orders_sum{};
		for (int i{ parameters - 1 }; i >= 0; --i)
		{
			// derivative order for current free parameter
			int derivative_order = orders[i];

			// safety check
			orders_sum += derivative_order;
			if (orders_sum > order)
			{
				throw std::exception("not implmented");
				//throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_LARGE, orders_sum, order);
			}

			while (derivative_order > 0)
			{
				--derivative_order;
				// as long as we differentiate according to current free parameter, // we have to skip the value part and dive into the derivative part
				// so we add the size of the value part to the base index
				index += sizes[i][m--];
			}
		}

		return index;
	}

	/** Convert an index from one (parameters, order) structure to another.
	 * @param index index of a partial derivative in source derivative structure
	 * @param src_p number of free parameters in source derivative structure
	 * @param src_derivatives_indirection derivatives indirection array for the source
	 * derivative structure
	 * @param dest_p number of free parameters in destination derivative structure
	 * @param dest_o derivation order in destination derivative structure
	 * @param dest_sizes sizes array for the destination derivative structure
	 * @return index of the partial derivative with the <em>same</em> characteristics
	 * in destination derivative structure
	 * @ if order is too large
	 */
	static int convert_index(const int& index, const int& src_p, const std::vector<std::vector<int>>& src_derivatives_indirection, const int& dest_p, const int& dest_o, const std::vector<std::vector<int>>& dest_sizes)
	{
		auto orders = std::vector<int>(dest_p);
		System.arraycopy(src_derivatives_indirection[index], 0, orders, 0, std::min(src_p, dest_p));
		return get_partial_derivative_index(dest_p, dest_o, dest_sizes, orders);
	}

public:
	/** Get the compiler for number of free parameters and order.
	 * @param parameters number of free parameters
	 * @param order derivation order
	 * @return cached rules set
	 * @ if order is too large
	 */
	static DS_Compiler get_compiler(const int& parameters, const int& order)
	{
		// get the cached compilers
		const auto cache = compilers.get();
		if (cache != NULL && cache.size() > parameters && cache[parameters].size() > order && cache[parameters][order] != NULL)
		{
			// the compiler has already been created
			return cache[parameters][order];
		}

		// we need to create more compilers
		const int max_parameters = std::max(parameters, cache == NULL ? 0 : cache.size());
		const int max_order = std::max(order, cache == NULL ? 0 : cache[0].size());
		const auto new_cache = std::vector<std::vector<DS_Compiler>(max_parameters + 1, std::vector<DS_Compiler>(max_order + 1));

		if (cache != NULL)
		{
			// preserve the already created compilers
			for (int i{}; i < cache.size(); ++i)
			{
				System.arraycopy(cache[i], 0, new_cache[i], 0, cache[i].size());
			}
		}

		// create the array in increasing diagonal order
		for (const int& diag{}; diag <= parameters + order; ++diag)
		{
			for (const int& o = std::max(0, diag - parameters); o <= std::min(order, diag); ++o)
			{
				const int p = diag - o;
				if (new_cache[p][o] == NULL)
				{
					const DS_Compiler value_compiler = p == 0
						? NULL
						: new_cache[p - 1][o];
					const DS_Compiler derivative_compiler = o == 0
						? NULL
						: new_cache[p][o - 1];
					new_cache[p][o] = DS_Compiler(p, o, value_compiler, derivative_compiler);
				}
			}
		}

		// atomically reset the cached compilers array
		compilers.compare_and_set(cache, new_cache);

		return new_cache[parameters][order];
	}

	/** Get the index of a partial derivative in the array.
	 * <p>
	 * If all orders are set to 0, then the 0<sup>th</sup> order derivative
	 * is returned, which is the value of the function.
	 * </p>
	 * <p>The indices of derivatives are between 0 and {@link #get_size() get_size()} - 1.
	 * Their specific order is fixed for a given compiler, but otherwise not
	 * publicly specified. There are however some simple cases which have guaranteed
	 * indices:
	 * </p>
	 * <ul>
	 *   <li>the index of 0<sup>th</sup> order derivative is always 0</li>
	 *   <li>if there is only 1 {@link #get_free_parameters() free parameter}, then the
	 *   derivatives are sorted in increasing derivation order (i.e. f at index 0, df/dp
	 *   at index 1, d<sup>2</sup>f/dp<sup>2</sup> at index 2 ...
	 *   d<sup>k</sup>f/dp<sup>k</sup> at index k),</li>
	 *   <li>if the {@link #get_order() derivation order} is 1, then the derivatives
	 *   are sorted in increasing free parameter order (i.e. f at index 0, df/dx<sub>1</sub>
	 *   at index 1, df/dx<sub>2</sub> at index 2 ... df/dx<sub>k</sub> at index k),</li>
	 *   <li>all other cases are not publicly specified</li>
	 * </ul>
	 * <p>
	 * This method is the inverse of method {@link #get_partial_derivative_ordersstatic_cast<int>(}
	 * </p>
	 * @param orders derivation orders with respect to each parameter
	 * @return index of the partial derivative
	 * @exception  if the numbers of parameters does not
	 * match the instance
	 * @exception  if sum of derivation orders is larger
	 * than the instance limits
	 * @see #get_partial_derivative_ordersstatic_cast<int>(
	 */
	int get_partial_derivative_index(const int& ... orders)
	{
		// safety check
		Math_Utils::check_dimension(orders.size(), get_free_parameters());
		return get_partial_derivative_index(parameters, order, sizes, orders);
	}

	/** Get the derivation orders for a specific index in the array.
	 * <p>
	 * This method is the inverse of {@link #get_partial_derivative_index(int...)}.
	 * </p>
	 * @param index of the partial derivative
	 * @return orders derivation orders with respect to each parameter
	 * @see #get_partial_derivative_index(int...)
	 */
	std::vector<int> get_partial_derivative_orders(const int& index)
	{
		return derivatives_indirection[index];
	}

	/** Get the number of free parameters.
	 * @return number of free parameters
	 */
	int get_free_parameters() const
	{
		return my_parameters;
	}

	/** Get the derivation order.
	 * @return derivation order
	 */
	int get_order() const
	{
		return my_order;
	}

	/** Get the array size required for holding partial derivatives data.
	 * <p>
	 * This number includes the single 0 order derivative element, which is
	 * guaranteed to be stored in the first element of the array.
	 * </p>
	 * @return array size required for holding partial derivatives data
	 */
	int get_size() const
	{
		return sizes[parameters][order];
	}

	/** Compute linear combination.
	 * The derivative structure built will be a1 * ds1 + a2 * ds2
	 * @param a1 first scale factor
	 * @param c1 first base (unscaled) component
	 * @param offset1 offset of first operand in its array
	 * @param a2 second scale factor
	 * @param c2 second base (unscaled) component
	 * @param offset2 offset of second operand in its array
	 * @param result array where result must be stored (it may be
	 * one of the input arrays)
	 * @param result_offset offset of the result in its array
	 */
	void linear_combination(const double& a1, const std::vector<double>& c1, const int& offset1, const double& a2, const std::vector<double>& c2, const int& offset2, const std::vector<double>& result, const int& result_offset)
	{
		for (int i{}; i < get_size(); ++i)
		{
			result[result_offset + i] = Math_Arrays::linear_combination(a1, c1[offset1 + i], a2, c2[offset2 + i]);
		}
	}

	/** Compute linear combination.
	 * The derivative structure built will be a1 * ds1 + a2 * ds2
	 * @param a1 first scale factor
	 * @param c1 first base (unscaled) component
	 * @param offset1 offset of first operand in its array
	 * @param a2 second scale factor
	 * @param c2 second base (unscaled) component
	 * @param offset2 offset of second operand in its array
	 * @param result array where result must be stored (it may be
	 * one of the input arrays)
	 * @param result_offset offset of the result in its array
	 * @param <T> the type of the function parameters and value
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	void linear_combination(const T& a1, const std::vector<T>& c1, const int& offset1, const T& a2, const std::vector<T>& c2, const int& offset2, std::vector<T>& result, const int& result_offset)
	{
		for (int i{}; i < get_size(); ++i)
		{
			result[result_offset + i] = a1.linear_combination(a1, c1[offset1 + i], a2, c2[offset2 + i]);
		}
	}

	/** Compute linear combination.
	 * The derivative structure built will be a1 * ds1 + a2 * ds2
	 * @param a1 first scale factor
	 * @param c1 first base (unscaled) component
	 * @param offset1 offset of first operand in its array
	 * @param a2 second scale factor
	 * @param c2 second base (unscaled) component
	 * @param offset2 offset of second operand in its array
	 * @param result array where result must be stored (it may be
	 * one of the input arrays)
	 * @param result_offset offset of the result in its array
	 * @param <T> the type of the function parameters and value
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	void linear_combination(const double& a1, const std::vector<T>& c1, const int& offset1, const double& a2, const std::vector<T>& c2, const int& offset2, std::vector<T>& result, const int& result_offset)
	{
		for (int i{}; i < get_size(); ++i)
		{
			result[result_offset + i] = c1[offset1].linear_combination(a1, c1[offset1 + i], a2, c2[offset2 + i]);
		}
	}

	/** Compute linear combination.
	 * The derivative structure built will be a1 * ds1 + a2 * ds2 + a3 * ds3 + a4 * ds4
	 * @param a1 first scale factor
	 * @param c1 first base (unscaled) component
	 * @param offset1 offset of first operand in its array
	 * @param a2 second scale factor
	 * @param c2 second base (unscaled) component
	 * @param offset2 offset of second operand in its array
	 * @param a3 third scale factor
	 * @param c3 third base (unscaled) component
	 * @param offset3 offset of third operand in its array
	 * @param result array where result must be stored (it may be
	 * one of the input arrays)
	 * @param result_offset offset of the result in its array
	 */
	void linear_combination(const double& a1, const std::vector<double>& c1, const int& offset1, const double& a2, const std::vector<double>& c2, const int& offset2, const double& a3, const std::vector<double> c3, const int& offset3, std::vector<double>& result, const int& result_offset)
	{
		for (int i{}; i < get_size(); ++i)
		{
			result[result_offset + i] = Math_Arrays::linear_combination(a1, c1[offset1 + i], a2, c2[offset2 + i], a3, c3[offset3 + i]);
		}
	}

	/** Compute linear combination.
	 * The derivative structure built will be a1 * ds1 + a2 * ds2 + a3 * ds3 + a4 * ds4
	 * @param a1 first scale factor
	 * @param c1 first base (unscaled) component
	 * @param offset1 offset of first operand in its array
	 * @param a2 second scale factor
	 * @param c2 second base (unscaled) component
	 * @param offset2 offset of second operand in its array
	 * @param a3 third scale factor
	 * @param c3 third base (unscaled) component
	 * @param offset3 offset of third operand in its array
	 * @param result array where result must be stored (it may be
	 * one of the input arrays)
	 * @param result_offset offset of the result in its array
	 * @param <T> the type of the function parameters and value
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	void linear_combination(const T& a1, const std::vector<T>& c1, const int& offset1, const T& a2, const std::vector<T>& c2, const int& offset2, const T& a3, const std::vector<T>& c3, const int& offset3, std::vector<T>& result, const int& result_offset)
	{
		for (int i{}; i < get_size(); ++i)
		{
			result[result_offset + i] = a1.linear_combination(a1, c1[offset1 + i], a2, c2[offset2 + i], a3, c3[offset3 + i]);
		}
	}

	/** Compute linear combination.
	 * The derivative structure built will be a1 * ds1 + a2 * ds2 + a3 * ds3 + a4 * ds4
	 * @param a1 first scale factor
	 * @param c1 first base (unscaled) component
	 * @param offset1 offset of first operand in its array
	 * @param a2 second scale factor
	 * @param c2 second base (unscaled) component
	 * @param offset2 offset of second operand in its array
	 * @param a3 third scale factor
	 * @param c3 third base (unscaled) component
	 * @param offset3 offset of third operand in its array
	 * @param result array where result must be stored (it may be
	 * one of the input arrays)
	 * @param result_offset offset of the result in its array
	 * @param <T> the type of the function parameters and value
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	void linear_combination(const double& a1, const std::vector<T>& c1, const int& offset1, const double& a2, const std::vector<T>& c2, const int& offset2, const double& a3, const std::vector<T>& c3, const int& offset3, std::vector<T>& result, const int& result_offset)
	{
		for (int i{}; i < get_size(); ++i)
		{
			result[result_offset + i] = c1[offset1].linear_combination(a1, c1[offset1 + i], a2, c2[offset2 + i], a3, c3[offset3 + i]);
		}
	}

	/** Compute linear combination.
	 * The derivative structure built will be a1 * ds1 + a2 * ds2 + a3 * ds3 + a4 * ds4
	 * @param a1 first scale factor
	 * @param c1 first base (unscaled) component
	 * @param offset1 offset of first operand in its array
	 * @param a2 second scale factor
	 * @param c2 second base (unscaled) component
	 * @param offset2 offset of second operand in its array
	 * @param a3 third scale factor
	 * @param c3 third base (unscaled) component
	 * @param offset3 offset of third operand in its array
	 * @param a4 fourth scale factor
	 * @param c4 fourth base (unscaled) component
	 * @param offset4 offset of fourth operand in its array
	 * @param result array where result must be stored (it may be
	 * one of the input arrays)
	 * @param result_offset offset of the result in its array
	 */
	void linear_combination(const double& a1, const std::vector<double>& c1, const int& offset1, const double& a2, const std::vector<double>& c2, const int& offset2, const double& a3, const std::vector<double>& c3, const int& offset3, const double& a4, const std::vector<double> c4, const int& offset4, std::vector<double>& result, const int& result_offset)
	{
		for (int i{}; i < get_size(); ++i)
		{
			result[result_offset + i] = Math_Arrays::linear_combination(a1, c1[offset1 + i], a2, c2[offset2 + i], a3, c3[offset3 + i], a4, c4[offset4 + i]);
		}
	}

	/** Compute linear combination.
	 * The derivative structure built will be a1 * ds1 + a2 * ds2 + a3 * ds3 + a4 * ds4
	 * @param a1 first scale factor
	 * @param c1 first base (unscaled) component
	 * @param offset1 offset of first operand in its array
	 * @param a2 second scale factor
	 * @param c2 second base (unscaled) component
	 * @param offset2 offset of second operand in its array
	 * @param a3 third scale factor
	 * @param c3 third base (unscaled) component
	 * @param offset3 offset of third operand in its array
	 * @param a4 fourth scale factor
	 * @param c4 fourth base (unscaled) component
	 * @param offset4 offset of fourth operand in its array
	 * @param result array where result must be stored (it may be
	 * one of the input arrays)
	 * @param result_offset offset of the result in its array
	 * @param <T> the type of the function parameters and value
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	void linear_combination(const T& a1, const std::vector<T>& c1, const int& offset1, const T& a2, const std::vector<T>& c2, const int& offset2, const T& a3, const std::vector<T>& c3, const int& offset3, const T& a4, const std::vector<T>& c4, const int& offset4, std::vector<T>& result, const int& result_offset)
	{
		for (int i{}; i < get_size(); ++i)
		{
			result[result_offset + i] =
				a1.linear_combination(a1, c1[offset1 + i], a2, c2[offset2 + i], a3, c3[offset3 + i], a4, c4[offset4 + i]);
		}
	}

	/** Compute linear combination.
	 * The derivative structure built will be a1 * ds1 + a2 * ds2 + a3 * ds3 + a4 * ds4
	 * @param a1 first scale factor
	 * @param c1 first base (unscaled) component
	 * @param offset1 offset of first operand in its array
	 * @param a2 second scale factor
	 * @param c2 second base (unscaled) component
	 * @param offset2 offset of second operand in its array
	 * @param a3 third scale factor
	 * @param c3 third base (unscaled) component
	 * @param offset3 offset of third operand in its array
	 * @param a4 fourth scale factor
	 * @param c4 fourth base (unscaled) component
	 * @param offset4 offset of fourth operand in its array
	 * @param result array where result must be stored (it may be
	 * one of the input arrays)
	 * @param result_offset offset of the result in its array
	 * @param <T> the type of the function parameters and value
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	void linear_combination(const double& a1, const std::vector<T>& c1, const int& offset1, const double& a2, const std::vector<T>& c2, const int& offset2, const double& a3, const std::vector<T>& c3, const int& offset3, const double& a4, const std::vector<T>& c4, const int& offset4, std::vector<T>& result, const int& result_offset)
	{
		for (int i{}; i < get_size(); ++i)
		{
			result[result_offset + i] = c1[offset1].linear_combination(a1, c1[offset1 + i], a2, c2[offset2 + i], a3, c3[offset3 + i], a4, c4[offset4 + i]);
		}
	}

	/** Perform addition of two derivative structures.
	 * @param lhs array holding left hand side of addition
	 * @param lhs_offset offset of the left hand side in its array
	 * @param rhs array right hand side of addition
	 * @param rhs_offset offset of the right hand side in its array
	 * @param result array where result must be stored (it may be
	 * one of the input arrays)
	 * @param result_offset offset of the result in its array
	 */
	void add(const std::vector<double>& lhs, const int& lhs_offset, const std::vector<double>& rhs, const int& rhs_offset, std::vector<double>& result, const int& result_offset)
	{
		for (int i{}; i < get_size(); ++i)
		{
			result[result_offset + i] = lhs[lhs_offset + i] + rhs[rhs_offset + i];
		}
	}

	/** Perform addition of two derivative structures.
	 * @param lhs array holding left hand side of addition
	 * @param lhs_offset offset of the left hand side in its array
	 * @param rhs array right hand side of addition
	 * @param rhs_offset offset of the right hand side in its array
	 * @param result array where result must be stored (it may be
	 * one of the input arrays)
	 * @param result_offset offset of the result in its array
	 * @param <T> the type of the function parameters and value
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	void add(const std::vector<T>& lhs, const int& lhs_offset, const std::vector<T>& rhs, const int& rhs_offset, std::vector<T>& result, const int& result_offset)
	{
		for (int i{}; i < get_size(); ++i)
		{
			result[result_offset + i] = lhs[lhs_offset + i].add(rhs[rhs_offset + i]);
		}
	}

	/** Perform subtraction of two derivative structures.
	 * @param lhs array holding left hand side of subtraction
	 * @param lhs_offset offset of the left hand side in its array
	 * @param rhs array right hand side of subtraction
	 * @param rhs_offset offset of the right hand side in its array
	 * @param result array where result must be stored (it may be
	 * one of the input arrays)
	 * @param result_offset offset of the result in its array
	 */
	void subtract(const std::vector<double>& lhs, const int& lhs_offset, const std::vector<double>& rhs, const int& rhs_offset, std::vector<double>& result, const int& result_offset)
	{
		for (int i{}; i < get_size(); ++i)
		{
			result[result_offset + i] = lhs[lhs_offset + i] - rhs[rhs_offset + i];
		}
	}

	/** Perform subtraction of two derivative structures.
	 * @param lhs array holding left hand side of subtraction
	 * @param lhs_offset offset of the left hand side in its array
	 * @param rhs array right hand side of subtraction
	 * @param rhs_offset offset of the right hand side in its array
	 * @param result array where result must be stored (it may be
	 * one of the input arrays)
	 * @param result_offset offset of the result in its array
	 * @param <T> the type of the function parameters and value
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	void subtract(const std::vector<T>& lhs, const int& lhs_offset, const std::vector<T>& rhs, const int& rhs_offset, std::vector<T>& result, const int& result_offset)
	{
		for (int i{}; i < get_size(); ++i)
		{
			result[result_offset + i] = lhs[lhs_offset + i].subtract(rhs[rhs_offset + i]);
		}
	}

	/** Perform multiplication of two derivative structures.
	  * @param lhs array holding left hand side of multiplication
	  * @param lhs_offset offset of the left hand side in its array
	  * @param rhs array right hand side of multiplication
	  * @param rhs_offset offset of the right hand side in its array
	  * @param result array where result must be stored (for
	  * multiplication the result array <em>cannot</em> be one of
	  * the input arrays)
	  * @param result_offset offset of the result in its array
	  */
	void multiply(const std::vector<double>& lhs, const int& lhs_offset, const std::vector<double>& rhs, const int& rhs_offset, std::vector<double>& result, const int& result_offset)
	{
		for (int i{}; i < mult_indirection.size(); ++i)
		{
			const std::vector<std::vector<int>> mapping_i = mult_indirection[i];
			double r{};
			for (int j{}; j < mapping_i.size(); ++j)
			{
				r += mapping_i[j][0] *
					lhs[lhs_offset + mapping_i[j][1]] *
					rhs[rhs_offset + mapping_i[j][2]];
			}
			result[result_offset + i] = r;
		}
	}

	/** Perform multiplication of two derivative structures.
	 * @param lhs array holding left hand side of multiplication
	 * @param lhs_offset offset of the left hand side in its array
	 * @param rhs array right hand side of multiplication
	 * @param rhs_offset offset of the right hand side in its array
	 * @param result array where result must be stored (for
	 * multiplication the result array <em>cannot</em> be one of
	 * the input arrays)
	 * @param result_offset offset of the result in its array
	 * @param <T> the type of the function parameters and value
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	void multiply(const std::vector<T>& lhs, const int& lhs_offset, const std::vector<T>& rhs, const int& rhs_offset, std::vector<T>& result, const int& result_offset)
	{
		T zero = lhs[lhs_offset].get_field().get_zero();
		for (int i{}; i < mult_indirection.size(); ++i)
		{
			const auto mapping_i = mult_indirection[i];
			T r = zero;
			for (int j{}; j < mapping_i.size(); ++j)
			{
				r = r.add(lhs[lhs_offset + mapping_i[j][1]].
					multiply(rhs[rhs_offset + mapping_i[j][2]]).
					multiply(mapping_i[j][0]));
			}
			result[result_offset + i] = r;
		}
	}

	/** Perform division of two derivative structures.
	 * @param lhs array holding left hand side of division
	 * @param lhs_offset offset of the left hand side in its array
	 * @param rhs array right hand side of division
	 * @param rhs_offset offset of the right hand side in its array
	 * @param result array where result must be stored (for
	 * division the result array <em>cannot</em> be one of
	 * the input arrays)
	 * @param result_offset offset of the result in its array
	 */
	void divide(const std::vector<double>& lhs, const int& lhs_offset, const std::vector<double>& rhs, const int& rhs_offset, std::vector<double>& result, const int& result_offset)
	{
		const auto reciprocal = std::vector<double>(get_size()];
		pow(rhs, lhs_offset, -1, reciprocal, 0);
		multiply(lhs, lhs_offset, reciprocal, 0, result, result_offset);
	}

	/** Perform division of two derivative structures.
	 * @param lhs array holding left hand side of division
	 * @param lhs_offset offset of the left hand side in its array
	 * @param rhs array right hand side of division
	 * @param rhs_offset offset of the right hand side in its array
	 * @param result array where result must be stored (for
	 * division the result array <em>cannot</em> be one of
	 * the input arrays)
	 * @param result_offset offset of the result in its array
	 * @param <T> the type of the function parameters and value
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	void divide(const std::vector<T>& lhs, const int& lhs_offset, const std::vector<T>& rhs, const int& rhs_offset, std::vector<T>& result, const int& result_offset)
	{
		const std::vector<T> reciprocal = Math_Arrays::build_array(lhs[lhs_offset].get_field(), get_size());
		pow(rhs, lhs_offset, -1, reciprocal, 0);
		multiply(lhs, lhs_offset, reciprocal, 0, result, result_offset);
	}

	/** Perform remainder of two derivative structures.
	 * @param lhs array holding left hand side of remainder
	 * @param lhs_offset offset of the left hand side in its array
	 * @param rhs array right hand side of remainder
	 * @param rhs_offset offset of the right hand side in its array
	 * @param result array where result must be stored (it may be
	 * one of the input arrays)
	 * @param result_offset offset of the result in its array
	 */
	void remainder(const std::vector<double>& lhs, const int& lhs_offset, const std::vector<double>& rhs, const int& rhs_offset, std::vector<double>& result, const int& result_offset)
	{
		// compute k such that lhs % rhs = lhs - k rhs
		const double rem = std::remainder(lhs[lhs_offset], rhs[rhs_offset]);
		const double k = std::rint((lhs[lhs_offset] - rem) / rhs[rhs_offset]);

		// set up value
		result[result_offset] = rem;

		// set up partial derivatives
		for (int i{ 1 }; i < get_size(); ++i)
		{
			result[result_offset + i] = lhs[lhs_offset + i] - k * rhs[rhs_offset + i];
		}
	}

	/** Perform remainder of two derivative structures.
	 * @param lhs array holding left hand side of remainder
	 * @param lhs_offset offset of the left hand side in its array
	 * @param rhs array right hand side of remainder
	 * @param rhs_offset offset of the right hand side in its array
	 * @param result array where result must be stored (it may be
	 * one of the input arrays)
	 * @param result_offset offset of the result in its array
	 * @param <T> the type of the function parameters and value
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	void remainder(const std::vector<T>& lhs, const int& lhs_offset, const std::vector<T>& rhs, const int& rhs_offset, std::vector<T>& result, const int& result_offset)
	{
		// compute k such that lhs % rhs = lhs - k rhs
		const T rem = lhs[lhs_offset].remainder(rhs[rhs_offset]);
		const double k = std::rint((lhs[lhs_offset].get_real() - rem.get_real()) / rhs[rhs_offset].get_real());

		// set up value
		result[result_offset] = rem;

		// set up partial derivatives
		for (int i{ 1 }; i < get_size(); ++i)
		{
			result[result_offset + i] = lhs[lhs_offset + i].subtract(rhs[rhs_offset + i].multiply(k));
		}
	}

	/** Compute power of a double to a derivative structure.
	 * @param a number to exponentiate
	 * @param operand array holding the power
	 * @param operand_offset offset of the power in its array
	 * @param result array where result must be stored (for
	 * power the result array <em>cannot</em> be the input
	 * array)
	 * @param result_offset offset of the result in its array
	 */
	void pow(const double& a, const std::vector<double>& operand, const int& operand_offset, std::vector<double>& result, const int& result_offset)
	{
		// create the function value and derivatives
		// [a^x, ln(a) a^x, ln(a)^2 a^x,, ln(a)^3 a^x, ... ]
		auto function = std::vector<double>(1 + order);
		if (a == 0)
		{
			if (operand[operand_offset] == 0)
			{
				function[0] = 1;
				double infinity = INFINITY;
				for (int i{ 1 }; i < function.size(); ++i)
				{
					infinity = -infinity;
					function[i] = infinity;
				}
			}
			else if (operand[operand_offset] < 0)
			{
				Arrays.fill(function, NAN);
			}
		}
		else
		{
			function[0] = std::pow(a, operand[operand_offset]);
			const double ln_a = std::log(a);
			for (int i{ 1 }; i < function.size(); ++i)
			{
				function[i] = ln_a * function[i - 1];
			}
		}

		// apply function composition
		compose(operand, operand_offset, function, result, result_offset);
	}

	/** Compute power of a double to a derivative structure.
	 * @param a number to exponentiate
	 * @param operand array holding the power
	 * @param operand_offset offset of the power in its array
	 * @param result array where result must be stored (for
	 * power the result array <em>cannot</em> be the input
	 * array)
	 * @param result_offset offset of the result in its array
	 * @param <T> the type of the function parameters and value
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	void pow(const double& a, const std::vector<T>& operand, const int& operand_offset, std::vector<T>& result, const int& result_offset)
	{
		const T zero = operand[operand_offset].get_field().get_zero();

		// create the function value and derivatives
		// [a^x, ln(a) a^x, ln(a)^2 a^x,, ln(a)^3 a^x, ... ]
		const auto function = Math_Arrays::build_array(operand[operand_offset].get_field(), 1 + order);
		if (a == 0)
		{
			if (operand[operand_offset].get_real() == 0)
			{
				function[0] = zero.add(1);
				T infinity = zero.add(INFINITY);
				for (int i{ 1 }; i < function.size(); ++i)
				{
					infinity = infinity.negate();
					function[i] = infinity;
				}
			}
			else if (operand[operand_offset].get_real() < 0)
			{
				Arrays.fill(function, zero.add(Double.NaN));
			}
		}
		else
		{
			function[0] = zero.add(a).pow(operand[operand_offset]);
			const double ln_a = std::log(a);
			for (int i{ 1 }; i < function.size(); ++i)
			{
				function[i] = function[i - 1].multiply(ln_a);
			}
		}

		// apply function composition
		compose(operand, operand_offset, function, result, result_offset);
	}

	/** Compute power of a derivative structure.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param p power to apply
	 * @param result array where result must be stored (for
	 * power the result array <em>cannot</em> be the input
	 * array)
	 * @param result_offset offset of the result in its array
	 */
	void pow(const std::vector<double>& operand, const int& operand_offset, const double p, std::vector<double>& result, const int& result_offset)
	{
		if (p == 0)
		{
			// special case, x^0 = 1 for all x
			result[result_offset] = 1.0;
			Arrays.fill(result, result_offset + 1, result_offset + get_size(), 0);
			return;
		}

		if (operand[operand_offset] == 0)
		{
			// special case, 0^p = 0 for all p
			Arrays.fill(result, result_offset, result_offset + get_size(), 0);
			return;
		}

		// create the function value and derivatives
		// [x^p, px^(p-1), p(p-1)x^(p-2), ... ]
		std::vector<double> function = std::vector<double>(1 + order);
		double xk = std::pow(operand[operand_offset], p - order);
		for (int i = order; i > 0; --i)
		{
			function[i] = xk;
			xk *= operand[operand_offset];
		}
		function[0] = xk;
		double coefficient = p;
		for (int i{ 1 }; i <= order; ++i)
		{
			function[i] *= coefficient;
			coefficient *= p - i;
		}

		// apply function composition
		compose(operand, operand_offset, function, result, result_offset);
	}

	/** Compute power of a derivative structure.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param p power to apply
	 * @param result array where result must be stored (for
	 * power the result array <em>cannot</em> be the input
	 * array)
	 * @param result_offset offset of the result in its array
	 * @param <T> the type of the function parameters and value
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	void pow(const std::vector<T>& operand, const int& operand_offset, const double p, std::vector<T>& result, const int& result_offset)
	{
		const Field<T> field = operand[operand_offset].get_field();

		if (p == 0)
		{
			// special case, x^0 = 1 for all x
			result[result_offset] = field.get_one();
			Arrays.fill(result, result_offset + 1, result_offset + get_size(), field.get_zero());
			return;
		}

		if (operand[operand_offset].get_real() == 0)
		{
			// special case, 0^p = 0 for all p
			Arrays.fill(result, result_offset, result_offset + get_size(), field.get_zero());
			return;
		}

		// create the function value and derivatives
		// [x^p, px^(p-1), p(p-1)x^(p-2), ... ]
		std::vector<T> function = Math_Arrays::build_array(field, 1 + order);
		T xk = operand[operand_offset].pow(p - order);
		for (int i = order; i > 0; --i)
		{
			function[i] = xk;
			xk = xk.multiply(operand[operand_offset]);
		}
		function[0] = xk;
		double coefficient = p;
		for (int i{ 1 }; i <= order; ++i)
		{
			function[i] = function[i].multiply(coefficient);
			coefficient *= p - i;
		}

		// apply function composition
		compose(operand, operand_offset, function, result, result_offset);
	}

	/** Compute integer power of a derivative structure.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param n power to apply
	 * @param result array where result must be stored (for
	 * power the result array <em>cannot</em> be the input
	 * array)
	 * @param result_offset offset of the result in its array
	 */
	void pow(const std::vector<double>& operand, const int& operand_offset, const int& n, std::vector<double>& result, const int& result_offset)
	{
		if (n == 0)
		{
			// special case, x^0 = 1 for all x
			result[result_offset] = 1.0;
			Arrays.fill(result, result_offset + 1, result_offset + get_size(), 0);
			return;
		}

		// create the power function value and derivatives
		// [x^n, nx^(n-1), n(n-1)x^(n-2), ... ]
		auto function = std::vector<double>(1 + order);

		if (n > 0)
		{
			// strictly positive power
			const int max_order = std::min(order, n);
			double xk = std::pow(operand[operand_offset], n - max_order);
			for (int i = max_order; i > 0; --i)
			{
				function[i] = xk;
				xk *= operand[operand_offset];
			}
			function[0] = xk;
		}
		else
		{
			// strictly negative power
			const double inv = 1.0 / operand[operand_offset];
			double xk = std::pow(inv, -n);
			for (int i{}; i <= order; ++i)
			{
				function[i] = xk;
				xk *= inv;
			}
		}

		double coefficient = n;
		for (int i{ 1 }; i <= order; ++i)
		{
			function[i] *= coefficient;
			coefficient *= n - i;
		}

		// apply function composition
		compose(operand, operand_offset, function, result, result_offset);
	}

	/** Compute integer power of a derivative structure.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param n power to apply
	 * @param result array where result must be stored (for
	 * power the result array <em>cannot</em> be the input
	 * array)
	 * @param result_offset offset of the result in its array
	 * @param <T> the type of the function parameters and value
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	void pow(const std::vector<T>& operand, const int& operand_offset, const int& n, std::vector<T>& result, const int& result_offset)
	{
		const Field<T> field = operand[operand_offset].get_field();

		if (n == 0)
		{
			// special case, x^0 = 1 for all x
			result[result_offset] = field.get_one();
			Arrays.fill(result, result_offset + 1, result_offset + get_size(), field.get_zero());
			return;
		}

		// create the power function value and derivatives
		// [x^n, nx^(n-1), n(n-1)x^(n-2), ... ]
		std::vector<T> function = Math_Arrays::build_array(field, 1 + order);

		if (n > 0)
		{
			// strictly positive power
			const int max_order = std::min(order, n);
			T xk = operand[operand_offset].pow(n - max_order);
			for (int i = max_order; i > 0; --i)
			{
				function[i] = xk;
				xk = xk.multiply(operand[operand_offset]);
			}
			function[0] = xk;
		}
		else
		{
			// strictly negative power
			const T inv = operand[operand_offset].reciprocal();
			T xk = inv.pow(-n);
			for (int i{}; i <= order; ++i)
			{
				function[i] = xk;
				xk = xk.multiply(inv);
			}
		}

		double coefficient = n;
		for (int i{ 1 }; i <= order; ++i)
		{
			function[i] = function[i].multiply(coefficient);
			coefficient *= n - i;
		}

		// apply function composition
		compose(operand, operand_offset, function, result, result_offset);
	}

	/** Compute power of a derivative structure.
	 * @param x array holding the base
	 * @param x_offset offset of the base in its array
	 * @param y array holding the exponent
	 * @param y_offset offset of the exponent in its array
	 * @param result array where result must be stored (for
	 * power the result array <em>cannot</em> be the input
	 * array)
	 * @param result_offset offset of the result in its array
	 */
	void pow(const std::vector<double> x, const int& x_offset, const std::vector<double>& y, const int& y_offset, std::vector<double>& result, const int& result_offset)
	{
		const std::vector<double> log_x = std::vector<double>(get_size()];
		log(x, x_offset, log_x, 0);
		const std::vector<double> y_log_x = std::vector<double>(get_size()];
		multiply(log_x, 0, y, y_offset, y_log_x, 0);
		exp(y_log_x, 0, result, result_offset);
	}

	/** Compute power of a derivative structure.
	 * @param x array holding the base
	 * @param x_offset offset of the base in its array
	 * @param y array holding the exponent
	 * @param y_offset offset of the exponent in its array
	 * @param result array where result must be stored (for
	 * power the result array <em>cannot</em> be the input
	 * array)
	 * @param result_offset offset of the result in its array
	 * @param <T> the type of the function parameters and value
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	void pow(const std::vector<T> x, const int& x_offset, const std::vector<T>& y, const int& y_offset, std::vector<T>& result, const int& result_offset)
	{
		const std::vector<T> log_x = Math_Arrays::build_array(x[x_offset].get_field(), get_size());
		log(x, x_offset, log_x, 0);
		const std::vector<T> y_log_x = Math_Arrays::build_array(x[x_offset].get_field(), get_size());
		multiply(log_x, 0, y, y_offset, y_log_x, 0);
		exp(y_log_x, 0, result, result_offset);
	}

	/** Compute n<sup>th</sup> root of a derivative structure.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param n order of the root
	 * @param result array where result must be stored (for
	 * n<sup>th</sup> root the result array <em>cannot</em> be the input
	 * array)
	 * @param result_offset offset of the result in its array
	 */
	void root_n(const std::vector<double>& operand, const int& operand_offset, const int& n, std::vector<double>& result, const int& result_offset)
	{
		// create the function value and derivatives
		// [x^(1/n), (1/n)x^((1/n)-1), (1-n)/n^2x^((1/n)-2), ... ]
		std::vector<double> function = std::vector<double>(1 + order);
		double xk;
		if (n == 2)
		{
			function[0] = std::sqrt(operand[operand_offset]);
			xk = 0.5 / function[0];
		}
		else if (n == 3)
		{
			function[0] = std::cbrt(operand[operand_offset]);
			xk = 1.0 / (3.0 * function[0] * function[0]);
		}
		else
		{
			function[0] = std::pow(operand[operand_offset], 1.0 / n);
			xk = 1.0 / (n * std::pow(function[0], n - 1));
		}
		const double n_reciprocal = 1.0 / n;
		const double x_reciprocal = 1.0 / operand[operand_offset];
		for (int i{ 1 }; i <= order; ++i)
		{
			function[i] = xk;
			xk *= x_reciprocal * (n_reciprocal - i);
		}

		// apply function composition
		compose(operand, operand_offset, function, result, result_offset);
	}

	/** Compute n<sup>th</sup> root of a derivative structure.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param n order of the root
	 * @param result array where result must be stored (for
	 * n<sup>th</sup> root the result array <em>cannot</em> be the input
	 * array)
	 * @param result_offset offset of the result in its array
	 * @param <T> the type of the function parameters and value
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	void root_n(const std::vector<T>& operand, const int& operand_offset, const int& n, std::vector<T>& result, const int& result_offset)
	{
		const Field<T> field = operand[operand_offset].get_field();

		// create the function value and derivatives
		// [x^(1/n), (1/n)x^((1/n)-1), (1-n)/n^2x^((1/n)-2), ... ]
		std::vector<T> function = Math_Arrays::build_array(field, 1 + order);
		T xk;
		if (n == 2)
		{
			function[0] = operand[operand_offset].sqrt();
			xk = function[0].add(function[0]).reciprocal();
		}
		else if (n == 3)
		{
			function[0] = operand[operand_offset].cbrt();
			xk = function[0].multiply(function[0]).multiply(3).reciprocal();
		}
		else
		{
			function[0] = operand[operand_offset].pow(1.0 / n);
			xk = function[0].pow(n - 1).multiply(n).reciprocal();
		}
		const double n_reciprocal = 1.0 / n;
		const T      x_reciprocal = operand[operand_offset].reciprocal();
		for (int i{ 1 }; i <= order; ++i)
		{
			function[i] = xk;
			xk = xk.multiply(x_reciprocal.multiply(n_reciprocal - i));
		}

		// apply function composition
		compose(operand, operand_offset, function, result, result_offset);
	}

	/** Compute exponential of a derivative structure.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param result array where result must be stored (for
	 * exponential the result array <em>cannot</em> be the input
	 * array)
	 * @param result_offset offset of the result in its array
	 */
	void exp(const std::vector<double>& operand, const int& operand_offset, std::vector<double>& result, const int& result_offset)
	{
		// create the function value and derivatives
		std::vector<double> function = std::vector<double>(1 + order);
		Arrays.fill(function, std::exp(operand[operand_offset]));

		// apply function composition
		compose(operand, operand_offset, function, result, result_offset);
	}

	/** Compute exponential of a derivative structure.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param result array where result must be stored (for
	 * exponential the result array <em>cannot</em> be the input
	 * array)
	 * @param result_offset offset of the result in its array
	 * @param <T> the type of the function parameters and value
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	void exp(const std::vector<T>& operand, const int& operand_offset, std::vector<T>& result, const int& result_offset)
	{
		const Field<T> field = operand[operand_offset].get_field();

		// create the function value and derivatives
		std::vector<T> function = Math_Arrays::build_array(field, 1 + order);
		Arrays.fill(function, operand[operand_offset].exp());

		// apply function composition
		compose(operand, operand_offset, function, result, result_offset);
	}

	/** Compute exp(x) - 1 of a derivative structure.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param result array where result must be stored (for
	 * exponential the result array <em>cannot</em> be the input
	 * array)
	 * @param result_offset offset of the result in its array
	 */
	void expm1(const std::vector<double>& operand, const int& operand_offset, std::vector<double>& result, const int& result_offset)
	{
		// create the function value and derivatives
		std::vector<double> function = std::vector<double>(1 + order);
		function[0] = std::expm1(operand[operand_offset]);
		Arrays.fill(function, 1, 1 + order, std::exp(operand[operand_offset]));

		// apply function composition
		compose(operand, operand_offset, function, result, result_offset);
	}

	/** Compute exp(x) - 1 of a derivative structure.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param result array where result must be stored (for
	 * exponential the result array <em>cannot</em> be the input
	 * array)
	 * @param result_offset offset of the result in its array
	 * @param <T> the type of the function parameters and value
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	void expm1(const std::vector<T>& operand, const int& operand_offset, std::vector<T>& result, const int& result_offset)
	{
		const Field<T> field = operand[operand_offset].get_field();

		// create the function value and derivatives
		std::vector<T> function = Math_Arrays::build_array(field, 1 + order);
		function[0] = operand[operand_offset].expm1();
		Arrays.fill(function, 1, 1 + order, operand[operand_offset].exp());

		// apply function composition
		compose(operand, operand_offset, function, result, result_offset);
	}

	/** Compute natural logarithm of a derivative structure.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param result array where result must be stored (for
	 * logarithm the result array <em>cannot</em> be the input
	 * array)
	 * @param result_offset offset of the result in its array
	 */
	void log(const std::vector<double>& operand, const int& operand_offset, std::vector<double>& result, const int& result_offset)
	{
		// create the function value and derivatives
		std::vector<double> function = std::vector<double>(1 + order);
		function[0] = std::log(operand[operand_offset]);
		if (order > 0)
		{
			double inv = 1.0 / operand[operand_offset];
			double xk = inv;
			for (int i{ 1 }; i <= order; ++i)
			{
				function[i] = xk;
				xk *= -i * inv;
			}
		}

		// apply function composition
		compose(operand, operand_offset, function, result, result_offset);
	}

	/** Compute natural logarithm of a derivative structure.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param result array where result must be stored (for
	 * logarithm the result array <em>cannot</em> be the input
	 * array)
	 * @param result_offset offset of the result in its array
	 * @param <T> the type of the function parameters and value
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	void log(const std::vector<T>& operand, const int& operand_offset, std::vector<T>& result, const int& result_offset)
	{
		const Field<T> field = operand[operand_offset].get_field();

		// create the function value and derivatives
		std::vector<T> function = Math_Arrays::build_array(field, 1 + order);
		function[0] = operand[operand_offset].log();
		if (order > 0)
		{
			T inv = operand[operand_offset].reciprocal();
			T xk = inv;
			for (int i{ 1 }; i <= order; ++i)
			{
				function[i] = xk;
				xk = xk.multiply(inv.multiply(-i));
			}
		}

		// apply function composition
		compose(operand, operand_offset, function, result, result_offset);
	}

	/** Computes shifted logarithm of a derivative structure.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param result array where result must be stored (for
	 * shifted logarithm the result array <em>cannot</em> be the input array)
	 * @param result_offset offset of the result in its array
	 */
	void log1p(const std::vector<double>& operand, const int& operand_offset, std::vector<double>& result, const int& result_offset)
	{
		// create the function value and derivatives
		auto function = std::vector<double>(1 + order);
		function[0] = std::log1p(operand[operand_offset]);
		if (order > 0)
		{
			double inv = 1.0 / (1.0 + operand[operand_offset]);
			double xk = inv;
			for (int i{ 1 }; i <= order; ++i)
			{
				function[i] = xk;
				xk *= -i * inv;
			}
		}

		// apply function composition
		compose(operand, operand_offset, function, result, result_offset);
	}

	/** Computes shifted logarithm of a derivative structure.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param result array where result must be stored (for
	 * shifted logarithm the result array <em>cannot</em> be the input array)
	 * @param result_offset offset of the result in its array
	 * @param <T> the type of the function parameters and value
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	void log1p(const std::vector<T>& operand, const int& operand_offset, std::vector<T>& result, const int& result_offset)
	{
		const Field<T> field = operand[operand_offset].get_field();

		// create the function value and derivatives
		std::vector<T> function = Math_Arrays::build_array(field, 1 + order);
		function[0] = operand[operand_offset].log1p();
		if (order > 0)
		{
			T inv = operand[operand_offset].add(1).reciprocal();
			T xk = inv;
			for (int i{ 1 }; i <= order; ++i)
			{
				function[i] = xk;
				xk = xk.multiply(inv.multiply(-i));
			}
		}

		// apply function composition
		compose(operand, operand_offset, function, result, result_offset);
	}

	/** Computes base 10 logarithm of a derivative structure.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param result array where result must be stored (for
	 * base 10 logarithm the result array <em>cannot</em> be the input array)
	 * @param result_offset offset of the result in its array
	 */
	void log10(const std::vector<double>& operand, const int& operand_offset, std::vector<double>& result, const int& result_offset)
	{
		// create the function value and derivatives
		std::vector<double> function = std::vector<double>(1 + order);
		function[0] = std::log10(operand[operand_offset]);
		if (order > 0)
		{
			double inv = 1.0 / operand[operand_offset];
			double xk = inv / std::log(10.0);
			for (int i{ 1 }; i <= order; ++i)
			{
				function[i] = xk;
				xk *= -i * inv;
			}
		}

		// apply function composition
		compose(operand, operand_offset, function, result, result_offset);
	}

	/** Computes base 10 logarithm of a derivative structure.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param result array where result must be stored (for
	 * base 10 logarithm the result array <em>cannot</em> be the input array)
	 * @param result_offset offset of the result in its array
	 * @param <T> the type of the function parameters and value
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	void log10(const std::vector<T>& operand, const int& operand_offset, std::vector<T>& result, const int& result_offset)
	{
		const Field<T> field = operand[operand_offset].get_field();

		// create the function value and derivatives
		std::vector<T> function = Math_Arrays::build_array(field, 1 + order);
		function[0] = operand[operand_offset].log10();
		if (order > 0)
		{
			T inv = operand[operand_offset].reciprocal();
			T xk = inv.multiply(1.0 / std::log(10.0));
			for (int i{ 1 }; i <= order; ++i)
			{
				function[i] = xk;
				xk = xk.multiply(inv.multiply(-i));
			}
		}

		// apply function composition
		compose(operand, operand_offset, function, result, result_offset);
	}

	/** Compute cosine of a derivative structure.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param result array where result must be stored (for
	 * cosine the result array <em>cannot</em> be the input
	 * array)
	 * @param result_offset offset of the result in its array
	 */
	void cos(const std::vector<double>& operand, const int& operand_offset, std::vector<double>& result, const int& result_offset)
	{
		// create the function value and derivatives
		std::vector<double> function = std::vector<double>(1 + order);
		const Sin_Cos sin_cos = Sin_Cos(operand[operand_offset]);
		function[0] = sin_cos.cos();
		if (order > 0)
		{
			function[1] = -sin_cos.sin();
			for (int i{ 2 }; i <= order; ++i)
			{
				function[i] = -function[i - 2];
			}
		}

		// apply function composition
		compose(operand, operand_offset, function, result, result_offset);
	}

	/** Compute cosine of a derivative structure.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param result array where result must be stored (for
	 * cosine the result array <em>cannot</em> be the input
	 * array)
	 * @param result_offset offset of the result in its array
	 * @param <T> the type of the function parameters and value
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	void cos(const std::vector<T>& operand, const int& operand_offset, std::vector<T>& result, const int& result_offset)
	{
		const Field<T> field = operand[operand_offset].get_field();

		// create the function value and derivatives
		std::vector<T> function = Math_Arrays::build_array(field, 1 + order);
		const Field_Sin_Cos<T> sin_cos = Sin_Cos(operand[operand_offset]);
		function[0] = sin_cos.cos();
		if (order > 0)
		{
			function[1] = sin_cos.sin().negate();
			if (order > 1)
			{
				function[2] = sin_cos.cos().negate();
				if (order > 2)
				{
					function[3] = sin_cos.sin();
					for (int i = 4; i <= order; ++i)
					{
						function[i] = function[i - 4];
					}
				}
			}
		}

		// apply function composition
		compose(operand, operand_offset, function, result, result_offset);
	}

	/** Compute sine of a derivative structure.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param result array where result must be stored (for
	 * sine the result array <em>cannot</em> be the input
	 * array)
	 * @param result_offset offset of the result in its array
	 */
	void sin(const std::vector<double>& operand, const int& operand_offset, std::vector<double>& result, const int& result_offset)
	{
		// create the function value and derivatives
		std::vector<double> function = std::vector<double>(1 + order);
		const Sin_Cos sin_cos = Sin_Cos(operand[operand_offset]);
		function[0] = sin_cos.sin();
		if (order > 0)
		{
			function[1] = sin_cos.cos();
			for (int i{ 2 }; i <= order; ++i)
			{
				function[i] = -function[i - 2];
			}
		}

		// apply function composition
		compose(operand, operand_offset, function, result, result_offset);
	}

	/** Compute sine of a derivative structure.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param result array where result must be stored (for
	 * sine the result array <em>cannot</em> be the input
	 * array)
	 * @param result_offset offset of the result in its array
	 * @param <T> the type of the function parameters and value
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	void sin(const std::vector<T>& operand, const int& operand_offset, std::vector<T>& result, const int& result_offset)
	{
		const Field<T> field = operand[operand_offset].get_field();

		// create the function value and derivatives
		std::vector<T> function = Math_Arrays::build_array(field, 1 + order);
		const Field_Sin_Cos<T> sin_cos = Sin_Cos(operand[operand_offset]);
		function[0] = sin_cos.sin();
		if (order > 0)
		{
			function[1] = sin_cos.cos();
			if (order > 1)
			{
				function[2] = sin_cos.sin().negate();
				if (order > 2)
				{
					function[3] = sin_cos.cos().negate();
					for (int i = 4; i <= order; ++i)
					{
						function[i] = function[i - 4];
					}
				}
			}
		}

		// apply function composition
		compose(operand, operand_offset, function, result, result_offset);
	}

	/** Compute combined sine and cosine of a derivative structure.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param sin array where sine must be stored (for
	 * sine the result array <em>cannot</em> be the input
	 * array)
	 * @param sin_offset offset of the result in its array
	 * @param cos array where cosine must be stored (for
	 * cosine the result array <em>cannot</em> be the input
	 * array)
	 * @param cos_offset offset of the result in its array
	 * @since 1.4
	 */
	void sin_cos(const std::vector<double>& operand, const int& operand_offset, const std::vector<double> sin, const int& sin_offset, const std::vector<double> cos, const int& cos_offset)
	{
		// create the function value and derivatives
		std::vector<double> function_sin = std::vector<double>(1 + order);
		std::vector<double> function_cos = std::vector<double>(1 + order);
		const Sin_Cos sin_cos = Sin_Cos(operand[operand_offset]);
		function_sin[0] = sin_cos.sin();
		function_cos[0] = sin_cos.cos();
		if (order > 0)
		{
			function_sin[1] = sin_cos.cos();
			function_cos[1] = -sin_cos.sin();
			for (int i{ 2 }; i <= order; ++i)
			{
				function_sin[i] = -function_sin[i - 2];
				function_cos[i] = -function_cos[i - 2];
			}
		}

		// apply function composition
		compose(operand, operand_offset, function_sin, sin, sin_offset);
		compose(operand, operand_offset, function_cos, cos, cos_offset);
	}

	/** Compute combined sine and cosine of a derivative structure.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param sin array where sine must be stored (for
	 * sine the result array <em>cannot</em> be the input
	 * array)
	 * @param sin_offset offset of the result in its array
	 * @param cos array where cosine must be stored (for
	 * cosine the result array <em>cannot</em> be the input
	 * array)
	 * @param cos_offset offset of the result in its array
	 * @param <T> the type of the function parameters and value
	 * @since 1.4
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	void sin_cos(const std::vector<T>& operand, const int& operand_offset, const std::vector<T> sin, const int& sin_offset, const std::vector<T> cos, const int& cos_offset)
	{
		const Field<T> field = operand[operand_offset].get_field();

		// create the function value and derivatives
		std::vector<T> function_sin = Math_Arrays::build_array(field, 1 + order);
		std::vector<T> function_cos = Math_Arrays::build_array(field, 1 + order);
		const Field_Sin_Cos<T> sin_cos = Sin_Cos(operand[operand_offset]);
		function_cos[0] = sin_cos.cos();
		if (order > 0)
		{
			function_cos[1] = sin_cos.sin().negate();
			if (order > 1)
			{
				function_cos[2] = sin_cos.cos().negate();
				if (order > 2)
				{
					function_cos[3] = sin_cos.sin();
					for (int i = 4; i <= order; ++i)
					{
						function_cos[i] = function_cos[i - 4];
					}
				}
			}
		}
		function_sin[0] = sin_cos.sin();
		System.arraycopy(function_cos, 0, function_sin, 1, order);

		// apply function composition
		compose(operand, operand_offset, function_sin, sin, sin_offset);
		compose(operand, operand_offset, function_cos, cos, cos_offset);
	}

	/** Compute tangent of a derivative structure.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param result array where result must be stored (for
	 * tangent the result array <em>cannot</em> be the input
	 * array)
	 * @param result_offset offset of the result in its array
	 */
	void tan(const std::vector<double>& operand, const int& operand_offset, std::vector<double>& result, const int& result_offset)
	{
		// create the function value and derivatives
		const auto function = std::vector<double>(1 + order);
		const double t = std::tan(operand[operand_offset]);
		function[0] = t;

		if (order > 0)
		{
			// the nth order derivative of tan has the form:
			// dn(tan(x)/dxn = P_n(tan(x))
			// where P_n(t) is a degree n+1 polynomial with same parity as n+1
			// P_0(t) = t, P_1(t) = 1 + t^2, P_2(t) = 2 t (1 + t^2) ...
			// the general recurrence relation for P_n is:
			// P_n(x) = (1+t^2) P_(n-1)'(t)
			// as per polynomial parity, we can store coefficients of both P_(n-1) and P_n in the same array
			auto p = std::vector<double>(order + 2];
			p[1] = 1;
			const double t2 = t * t;
			for (const int n{ 1 }; n <= order; ++n)
			{
				// update and evaluate polynomial P_n(t)
				double v{};
				p[n + 1] = n * p[n];
				for (int k = n + 1; k >= 0; k -= 2)
				{
					v = v * t2 + p[k];
					if (k > 2)
					{
						p[k - 2] = (k - 1) * p[k - 1] + (k - 3) * p[k - 3];
					}
					else if (k == 2)
					{
						p[0] = p[1];
					}
				}
				if ((n & 0x1) == 0)
				{
					v *= t;
				}

				function[n] = v;
			}
		}

		// apply function composition
		compose(operand, operand_offset, function, result, result_offset);
	}

	/** Compute tangent of a derivative structure.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param result array where result must be stored (for
	 * tangent the result array <em>cannot</em> be the input
	 * array)
	 * @param result_offset offset of the result in its array
	 * @param <T> the type of the function parameters and value
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	void tan(const std::vector<T>& operand, const int& operand_offset, std::vector<T>& result, const int& result_offset)
	{
		const Field<T> field = operand[operand_offset].get_field();

		// create the function value and derivatives
		std::vector<T> function = Math_Arrays::build_array(field, 1 + order);
		const T t = operand[operand_offset].tan();
		function[0] = t;

		if (order > 0)
		{
			// the nth order derivative of tan has the form:
			// dn(tan(x)/dxn = P_n(tan(x))
			// where P_n(t) is a degree n+1 polynomial with same parity as n+1
			// P_0(t) = t, P_1(t) = 1 + t^2, P_2(t) = 2 t (1 + t^2) ...
			// the general recurrence relation for P_n is:
			// P_n(x) = (1+t^2) P_(n-1)'(t)
			// as per polynomial parity, we can store coefficients of both P_(n-1) and P_n in the same array
			const std::vector<T> p = Math_Arrays::build_array(field, order + 2);
			p[1] = field.get_one();
			const T t2 = t.multiply(t);
			for (const int n{ 1 }; n <= order; ++n)
			{
				// update and evaluate polynomial P_n(t)
				T v = field.get_zero();
				p[n + 1] = p[n].multiply(n);
				for (int k = n + 1; k >= 0; k -= 2)
				{
					v = v.multiply(t2).add(p[k]);
					if (k > 2)
					{
						p[k - 2] = p[k - 1].multiply(k - 1).add(p[k - 3].multiply(k - 3));
					}
					else if (k == 2)
					{
						p[0] = p[1];
					}
				}
				if ((n & 0x1) == 0)
				{
					v = v.multiply(t);
				}

				function[n] = v;
			}
		}

		// apply function composition
		compose(operand, operand_offset, function, result, result_offset);
	}

	/** Compute arc cosine of a derivative structure.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param result array where result must be stored (for
	 * arc cosine the result array <em>cannot</em> be the input
	 * array)
	 * @param result_offset offset of the result in its array
	 */
	void acos(const std::vector<double>& operand, const int& operand_offset, std::vector<double>& result, const int& result_offset)
	{
		// create the function value and derivatives
		std::vector<double> function = std::vector<double>(1 + order);
		const double x = operand[operand_offset];
		function[0] = std::acos(x);
		if (order > 0)
		{
			// the nth order derivative of acos has the form:
			// dn(acos(x)/dxn = P_n(x) / [1 - x^2]^((2n-1)/2)
			// where P_n(x) is a degree n-1 polynomial with same parity as n-1
			// P_1(x) = -1, P_2(x) = -x, P_3(x) = -2x^2 - 1 ...
			// the general recurrence relation for P_n is:
			// P_n(x) = (1-x^2) P_(n-1)'(x) + (2n-3) x P_(n-1)(x)
			// as per polynomial parity, we can store coefficients of both P_(n-1) and P_n in the same array
			auto p = std::vector<double>(order);
			p[0] = -1;
			const double x2 = x * x;
			const double f = 1.0 / (1 - x2);
			double coeff = std::sqrt(f);
			function[1] = coeff * p[0];
			for (const int& n = 2; n <= order; ++n)
			{
				// update and evaluate polynomial P_n(x)
				double v{};
				p[n - 1] = (n - 1) * p[n - 2];
				for (int k{ n - 1 }; k >= 0; k -= 2)
				{
					v = v * x2 + p[k];
					if (k > 2)
					{
						p[k - 2] = (k - 1) * p[k - 1] + (2 * n - k) * p[k - 3];
					}
					else if (k == 2)
					{
						p[0] = p[1];
					}
				}
				if ((n & 0x1) == 0)
				{
					v *= x;
				}

				coeff *= f;
				function[n] = coeff * v;
			}
		}

		// apply function composition
		compose(operand, operand_offset, function, result, result_offset);
	}

	/** Compute arc cosine of a derivative structure.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param result array where result must be stored (for
	 * arc cosine the result array <em>cannot</em> be the input
	 * array)
	 * @param result_offset offset of the result in its array
	 * @param <T> the type of the function parameters and value
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	void acos(const std::vector<T>& operand, const int& operand_offset, std::vector<T>& result, const int& result_offset)
	{
		const Field<T> field = operand[operand_offset].get_field();

		// create the function value and derivatives
		std::vector<T> function = Math_Arrays::build_array(field, 1 + order);
		const T x = operand[operand_offset];
		function[0] = x.acos();
		if (order > 0)
		{
			// the nth order derivative of acos has the form:
			// dn(acos(x)/dxn = P_n(x) / [1 - x^2]^((2n-1)/2)
			// where P_n(x) is a degree n-1 polynomial with same parity as n-1
			// P_1(x) = -1, P_2(x) = -x, P_3(x) = -2x^2 - 1 ...
			// the general recurrence relation for P_n is:
			// P_n(x) = (1-x^2) P_(n-1)'(x) + (2n-3) x P_(n-1)(x)
			// as per polynomial parity, we can store coefficients of both P_(n-1) and P_n in the same array
			const std::vector<T> p = Math_Arrays::build_array(field, order);
			p[0] = field.get_one().negate();
			const T x2 = x.multiply(x);
			const T f = x2.subtract(1).negate().reciprocal();
			T coeff = f.sqrt();
			function[1] = coeff.multiply(p[0]);
			for (const int& n = 2; n <= order; ++n)
			{
				// update and evaluate polynomial P_n(x)
				T v = field.get_zero();
				p[n - 1] = p[n - 2].multiply(n - 1);
				for (int k{ n - 1 }; k >= 0; k -= 2)
				{
					v = v.multiply(x2).add(p[k]);
					if (k > 2)
					{
						p[k - 2] = p[k - 1].multiply(k - 1).add(p[k - 3].multiply(2 * n - k));
					}
					else if (k == 2)
					{
						p[0] = p[1];
					}
				}
				if ((n & 0x1) == 0)
				{
					v = v.multiply(x);
				}

				coeff = coeff.multiply(f);
				function[n] = coeff.multiply(v);
			}
		}

		// apply function composition
		compose(operand, operand_offset, function, result, result_offset);
	}

	/** Compute arc sine of a derivative structure.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param result array where result must be stored (for
	 * arc sine the result array <em>cannot</em> be the input
	 * array)
	 * @param result_offset offset of the result in its array
	 */
	void asin(const std::vector<double>& operand, const int& operand_offset, std::vector<double>& result, const int& result_offset)
	{
		// create the function value and derivatives
		auto function = std::vector<double>(1 + order);
		const double x = operand[operand_offset];
		function[0] = std::asin(x);
		if (order > 0)
		{
			// the nth order derivative of asin has the form:
			// dn(asin(x)/dxn = P_n(x) / [1 - x^2]^((2n-1)/2)
			// where P_n(x) is a degree n-1 polynomial with same parity as n-1
			// P_1(x) = 1, P_2(x) = x, P_3(x) = 2x^2 + 1 ...
			// the general recurrence relation for P_n is:
			// P_n(x) = (1-x^2) P_(n-1)'(x) + (2n-3) x P_(n-1)(x)
			// as per polynomial parity, we can store coefficients of both P_(n-1) and P_n in the same array
			auto p = std::vector<double>(order);
			p[0] = 1;
			const double x2 = x * x;
			const double f = 1.0 / (1 - x2);
			double coeff = std::sqrt(f);
			function[1] = coeff * p[0];
			for (const int& n = 2; n <= order; ++n)
			{
				// update and evaluate polynomial P_n(x)
				double v{};
				p[n - 1] = (n - 1) * p[n - 2];
				for (int k{ n - 1 }; k >= 0; k -= 2)
				{
					v = v * x2 + p[k];
					if (k > 2)
					{
						p[k - 2] = (k - 1) * p[k - 1] + (2 * n - k) * p[k - 3];
					}
					else if (k == 2)
					{
						p[0] = p[1];
					}
				}
				if ((n & 0x1) == 0)
				{
					v *= x;
				}

				coeff *= f;
				function[n] = coeff * v;
			}
		}

		// apply function composition
		compose(operand, operand_offset, function, result, result_offset);
	}

	/** Compute arc sine of a derivative structure.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param result array where result must be stored (for
	 * arc sine the result array <em>cannot</em> be the input
	 * array)
	 * @param result_offset offset of the result in its array
	 * @param <T> the type of the function parameters and value
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	void asin(const std::vector<T>& operand, const int& operand_offset, std::vector<T>& result, const int& result_offset)
	{
		const Field<T> field = operand[operand_offset].get_field();

		// create the function value and derivatives
		std::vector<T> function = Math_Arrays::build_array(field, 1 + order);
		const T x = operand[operand_offset];
		function[0] = x.asin();
		if (order > 0)
		{
			// the nth order derivative of asin has the form:
			// dn(asin(x)/dxn = P_n(x) / [1 - x^2]^((2n-1)/2)
			// where P_n(x) is a degree n-1 polynomial with same parity as n-1
			// P_1(x) = 1, P_2(x) = x, P_3(x) = 2x^2 + 1 ...
			// the general recurrence relation for P_n is:
			// P_n(x) = (1-x^2) P_(n-1)'(x) + (2n-3) x P_(n-1)(x)
			// as per polynomial parity, we can store coefficients of both P_(n-1) and P_n in the same array
			const std::vector<T> p = Math_Arrays::build_array(field, order);
			p[0] = field.get_one();
			const T x2 = x.multiply(x);
			const T f = x2.subtract(1).negate().reciprocal();
			T coeff = f.sqrt();
			function[1] = coeff.multiply(p[0]);
			for (const int& n = 2; n <= order; ++n)
			{
				// update and evaluate polynomial P_n(x)
				T v = field.get_zero();
				p[n - 1] = p[n - 2].multiply(n - 1);
				for (int k{ n - 1 }; k >= 0; k -= 2)
				{
					v = v.multiply(x2).add(p[k]);
					if (k > 2)
					{
						p[k - 2] = p[k - 1].multiply(k - 1).add(p[k - 3].multiply(2 * n - k));
					}
					else if (k == 2)
					{
						p[0] = p[1];
					}
				}
				if ((n & 0x1) == 0)
				{
					v = v.multiply(x);
				}

				coeff = coeff.multiply(f);
				function[n] = coeff.multiply(v);
			}
		}

		// apply function composition
		compose(operand, operand_offset, function, result, result_offset);
	}

	/** Compute arc tangent of a derivative structure.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param result array where result must be stored (for
	 * arc tangent the result array <em>cannot</em> be the input
	 * array)
	 * @param result_offset offset of the result in its array
	 */
	void atan(const std::vector<double>& operand, const int& operand_offset, std::vector<double>& result, const int& result_offset)
	{
		// create the function value and derivatives
		std::vector<double> function = std::vector<double>(1 + order);
		const double x = operand[operand_offset];
		function[0] = std::atan(x);
		if (order > 0)
		{
			// the nth order derivative of atan has the form:
			// dn(atan(x)/dxn = Q_n(x) / (1 + x^2)^n
			// where Q_n(x) is a degree n-1 polynomial with same parity as n-1
			// Q_1(x) = 1, Q_2(x) = -2x, Q_3(x) = 6x^2 - 2 ...
			// the general recurrence relation for Q_n is:
			// Q_n(x) = (1+x^2) Q_(n-1)'(x) - 2(n-1) x Q_(n-1)(x)
			// as per polynomial parity, we can store coefficients of both Q_(n-1) and Q_n in the same array
			auto q = std::vector<double>(order);
			q[0] = 1;
			const double x2 = x * x;
			const double f = 1.0 / (1 + x2);
			double coeff = f;
			function[1] = coeff * q[0];
			for (const int& n = 2; n <= order; ++n)
			{
				// update and evaluate polynomial Q_n(x)
				double v{};
				q[n - 1] = -n * q[n - 2];
				for (int k{ n - 1 }; k >= 0; k -= 2)
				{
					v = v * x2 + q[k];
					if (k > 2)
					{
						q[k - 2] = (k - 1) * q[k - 1] + (k - 1 - 2 * n) * q[k - 3];
					}
					else if (k == 2)
					{
						q[0] = q[1];
					}
				}
				if ((n & 0x1) == 0)
				{
					v *= x;
				}

				coeff *= f;
				function[n] = coeff * v;
			}
		}

		// apply function composition
		compose(operand, operand_offset, function, result, result_offset);
	}

	/** Compute arc tangent of a derivative structure.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param result array where result must be stored (for
	 * arc tangent the result array <em>cannot</em> be the input
	 * array)
	 * @param result_offset offset of the result in its array
	 * @param <T> the type of the function parameters and value
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	void atan(const std::vector<T>& operand, const int& operand_offset, std::vector<T>& result, const int& result_offset)
	{
		const Field<T> field = operand[operand_offset].get_field();

		// create the function value and derivatives
		std::vector<T> function = Math_Arrays::build_array(field, 1 + order);
		const T x = operand[operand_offset];
		function[0] = x.atan();
		if (order > 0)
		{
			// the nth order derivative of atan has the form:
			// dn(atan(x)/dxn = Q_n(x) / (1 + x^2)^n
			// where Q_n(x) is a degree n-1 polynomial with same parity as n-1
			// Q_1(x) = 1, Q_2(x) = -2x, Q_3(x) = 6x^2 - 2 ...
			// the general recurrence relation for Q_n is:
			// Q_n(x) = (1+x^2) Q_(n-1)'(x) - 2(n-1) x Q_(n-1)(x)
			// as per polynomial parity, we can store coefficients of both Q_(n-1) and Q_n in the same array
			const std::vector<T> q = Math_Arrays::build_array(field, order);
			q[0] = field.get_one();
			const T x2 = x.multiply(x);
			const T f = x2.add(1).reciprocal();
			T coeff = f;
			function[1] = coeff.multiply(q[0]);
			for (const int n{ 2 }; n <= order; ++n)
			{
				// update and evaluate polynomial Q_n(x)
				T v = field.get_zero();
				q[n - 1] = q[n - 2].multiply(-n);
				for (int k{ n - 1 }; k >= 0; k -= 2)
				{
					v = v.multiply(x2).add(q[k]);
					if (k > 2)
					{
						q[k - 2] = q[k - 1].multiply(k - 1).add(q[k - 3].multiply(k - 1 - 2 * n));
					}
					else if (k == 2)
					{
						q[0] = q[1];
					}
				}
				if ((n & 0x1) == 0)
				{
					v = v.multiply(x);
				}

				coeff = coeff.multiply(f);
				function[n] = coeff.multiply(v);
			}
		}

		// apply function composition
		compose(operand, operand_offset, function, result, result_offset);
	}

	/** Compute two arguments arc tangent of a derivative structure.
	 * @param y array holding the first operand
	 * @param y_offset offset of the first operand in its array
	 * @param x array holding the second operand
	 * @param x_offset offset of the second operand in its array
	 * @param result array where result must be stored (for
	 * two arguments arc tangent the result array <em>cannot</em>
	 * be the input array)
	 * @param result_offset offset of the result in its array
	 */
	void atan2(const std::vector<double>& y, const int& y_offset, const std::vector<double>& x, const int& x_offset, std::vector<double>& result, const int& result_offset)
	{
		// compute r = sqrt(x^2+y^2)
		auto tmp1 = std::vector<double>(get_size());
		multiply(x, x_offset, x, x_offset, tmp1, 0);      // x^2
		auto tmp2 = std::vector<double>(get_size());
		multiply(y, y_offset, y, y_offset, tmp2, 0);      // y^2
		add(tmp1, 0, tmp2, 0, tmp2, 0);                 // x^2 + y^2
		root_n(tmp2, 0, 2, tmp1, 0);                     // r = sqrt(x^2 + y^2)

		if (x[x_offset] >= 0)
		{
			// compute atan2(y, x) = 2 atan(y / (r + x))
			add(tmp1, 0, x, x_offset, tmp2, 0);          // r + x
			divide(y, y_offset, tmp2, 0, tmp1, 0);       // y /(r + x)
			atan(tmp1, 0, tmp2, 0);                     // atan(y / (r + x))
			for (int i{}; i < tmp2.size(); ++i)
			{
				result[result_offset + i] = 2 * tmp2[i]; // 2 * atan(y / (r + x))
			}
		}
		else
		{
			// compute atan2(y, x) = +/- pi - 2 atan(y / (r - x))
			subtract(tmp1, 0, x, x_offset, tmp2, 0);     // r - x
			divide(y, y_offset, tmp2, 0, tmp1, 0);       // y /(r - x)
			atan(tmp1, 0, tmp2, 0);                     // atan(y / (r - x))
			result[result_offset] = ((tmp2[0] <= 0) ? -std::numbers::pi : std::numbers::pi) - 2 * tmp2[0]; // +/-pi - 2 * atan(y / (r - x))
			for (int i{ 1 }; i < tmp2.size(); ++i)
			{
				result[result_offset + i] = -2 * tmp2[i]; // +/-pi - 2 * atan(y / (r - x))
			}
		}

		// fix value to take special cases (+0/+0, +0/-0, -0/+0, -0/-0, +/-infinity) correctly
		result[result_offset] = std::atan2(y[y_offset], x[x_offset]);
	}

	/** Compute two arguments arc tangent of a derivative structure.
	 * @param y array holding the first operand
	 * @param y_offset offset of the first operand in its array
	 * @param x array holding the second operand
	 * @param x_offset offset of the second operand in its array
	 * @param result array where result must be stored (for
	 * two arguments arc tangent the result array <em>cannot</em>
	 * be the input array)
	 * @param result_offset offset of the result in its array
	 * @param <T> the type of the function parameters and value
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	void atan2(const std::vector<T>& y, const int& y_offset, const std::vector<T>& x, const int& x_offset, std::vector<T>& result, const int& result_offset)
	{
		const Field<T> field = y[y_offset].get_field();

		// compute r = sqrt(x^2+y^2)
		std::vector<T> tmp1 = Math_Arrays::build_array(field, get_size());
		multiply(x, x_offset, x, x_offset, tmp1, 0);      // x^2
		std::vector<T> tmp2 = Math_Arrays::build_array(field, get_size());
		multiply(y, y_offset, y, y_offset, tmp2, 0);      // y^2
		add(tmp1, 0, tmp2, 0, tmp2, 0);                 // x^2 + y^2
		root_n(tmp2, 0, 2, tmp1, 0);                     // r = sqrt(x^2 + y^2)

		if (x[x_offset].get_real() >= 0)
		{
			// compute atan2(y, x) = 2 atan(y / (r + x))
			add(tmp1, 0, x, x_offset, tmp2, 0);          // r + x
			divide(y, y_offset, tmp2, 0, tmp1, 0);       // y /(r + x)
			atan(tmp1, 0, tmp2, 0);                     // atan(y / (r + x))
			for (int i{}; i < tmp2.size(); ++i)
			{
				result[result_offset + i] = tmp2[i].add(tmp2[i]); // 2 * atan(y / (r + x))
			}
		}
		else
		{
			// compute atan2(y, x) = +/- pi - 2 atan(y / (r - x))
			subtract(tmp1, 0, x, x_offset, tmp2, 0);     // r - x
			divide(y, y_offset, tmp2, 0, tmp1, 0);       // y /(r - x)
			atan(tmp1, 0, tmp2, 0);                     // atan(y / (r - x))
			result[result_offset] = tmp2[0].add(tmp2[0]).negate().
				add((tmp2[0].get_real() <= 0) ? -std::numbers::pi : std::numbers::pi); // +/-pi - 2 * atan(y / (r - x))
			for (int i{ 1 }; i < tmp2.size(); ++i)
			{
				result[result_offset + i] = tmp2[i].add(tmp2[i]).negate(); // +/-pi - 2 * atan(y / (r - x))
			}
		}

		// fix value to take special cases (+0/+0, +0/-0, -0/+0, -0/-0, +/-infinity) correctly
		result[result_offset] = y[y_offset].atan2(x[x_offset]);
	}

	/** Compute hyperbolic cosine of a derivative structure.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param result array where result must be stored (for
	 * hyperbolic cosine the result array <em>cannot</em> be the input
	 * array)
	 * @param result_offset offset of the result in its array
	 */
	void cosh(const std::vector<double>& operand, const int& operand_offset, std::vector<double>& result, const int& result_offset)
	{
		// create the function value and derivatives
		auto function = std::vector<double>(1 + order);
		function[0] = std::cosh(operand[operand_offset]);
		if (order > 0)
		{
			function[1] = std::sinh(operand[operand_offset]);
			for (int i{ 2 }; i <= order; ++i)
			{
				function[i] = function[i - 2];
			}
		}

		// apply function composition
		compose(operand, operand_offset, function, result, result_offset);
	}

	/** Compute hyperbolic cosine of a derivative structure.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param result array where result must be stored (for
	 * hyperbolic cosine the result array <em>cannot</em> be the input
	 * array)
	 * @param result_offset offset of the result in its array
	 * @param <T> the type of the function parameters and value
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	void cosh(const std::vector<T>& operand, const int& operand_offset, std::vector<T>& result, const int& result_offset)
	{
		const Field<T> field = operand[operand_offset].get_field();

		// create the function value and derivatives
		std::vector<T> function = Math_Arrays::build_array(field, 1 + order);
		function[0] = operand[operand_offset].cosh();
		if (order > 0)
		{
			function[1] = operand[operand_offset].sinh();
			for (int i{ 2 }; i <= order; ++i)
			{
				function[i] = function[i - 2];
			}
		}

		// apply function composition
		compose(operand, operand_offset, function, result, result_offset);
	}

	/** Compute hyperbolic sine of a derivative structure.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param result array where result must be stored (for
	 * hyperbolic sine the result array <em>cannot</em> be the input
	 * array)
	 * @param result_offset offset of the result in its array
	 */
	void sinh(const std::vector<double>& operand, const int& operand_offset, std::vector<double>& result, const int& result_offset)
	{
		// create the function value and derivatives
		auto function = std::vector<double>(1 + order);
		function[0] = std::sinh(operand[operand_offset]);
		if (order > 0)
		{
			function[1] = std::cosh(operand[operand_offset]);
			for (int i{ 2 }; i <= order; ++i)
			{
				function[i] = function[i - 2];
			}
		}

		// apply function composition
		compose(operand, operand_offset, function, result, result_offset);
	}

	/** Compute hyperbolic sine of a derivative structure.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param result array where result must be stored (for
	 * hyperbolic sine the result array <em>cannot</em> be the input
	 * array)
	 * @param result_offset offset of the result in its array
	 * @param <T> the type of the function parameters and value
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	void sinh(const std::vector<T>& operand, const int& operand_offset, std::vector<T>& result, const int& result_offset)
	{
		const Field<T> field = operand[operand_offset].get_field();

		// create the function value and derivatives
		std::vector<T> function = Math_Arrays::build_array(field, 1 + order);
		function[0] = operand[operand_offset].sinh();
		if (order > 0)
		{
			function[1] = operand[operand_offset].cosh();
			for (int i{ 2 }; i <= order; ++i)
			{
				function[i] = function[i - 2];
			}
		}

		// apply function composition
		compose(operand, operand_offset, function, result, result_offset);
	}

	/** Compute combined hyperbolic sine and cosine of a derivative structure.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param sinh array where hyperbolic sine must be stored (for
	 * sine the result array <em>cannot</em> be the input
	 * array)
	 * @param sinh_offset offset of the result in its array
	 * @param cosh array where hyperbolic <em>cannot</em> be the input
	 * array)
	 * @param cosh_offset offset of the result in its array
	 * @since 2.0
	 */
	void sinh_cosh(const std::vector<double>& operand, const int& operand_offset, const std::vector<double>& sinh, const int& sinh_offset, const std::vector<double> cosh, const int& cosh_offset)
	{
		// create the function value and derivatives
		auto function_sinh = std::vector<double>(1 + order);
		auto function_cosh = std::vector<double>(1 + order);
		const auto sinh_cosh = std::sinh_cosh(operand[operand_offset]);
		function_sinh[0] = sinh_cosh.sinh();
		function_cosh[0] = sinh_cosh.cosh();
		if (order > 0)
		{
			function_sinh[1] = sinh_cosh.cosh();
			function_cosh[1] = sinh_cosh.sinh();
			for (int i{ 2 }; i <= order; ++i)
			{
				function_sinh[i] = function_sinh[i - 2];
				function_cosh[i] = function_cosh[i - 2];
			}
		}

		// apply function composition
		compose(operand, operand_offset, function_sinh, sinh, sinh_offset);
		compose(operand, operand_offset, function_cosh, cosh, cosh_offset);
	}

	/** Compute combined hyperbolic sine and cosine of a derivative structure.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param sinh array where hyperbolic sine must be stored (for
	 * sine the result array <em>cannot</em> be the input
	 * array)
	 * @param sinh_offset offset of the result in its array
	 * @param cosh array where hyperbolic cosine must be stored (for
	 * cosine the result array <em>cannot</em> be the input
	 * array)
	 * @param cosh_offset offset of the result in its array
	 * @param <T> the type of the function parameters and value
	 * @since 1.4
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	void sinh_cosh(const std::vector<T>& operand, const int& operand_offset, const std::vector<T>& sinh, const int& sinh_offset, const std::vector<T> cosh, const int& cosh_offset)
	{
		const Field<T> field = operand[operand_offset].get_field();

		// create the function value and derivatives
		std::vector<T> function_sinh = Math_Arrays::build_array(field, 1 + order);
		std::vector<T> function_cosh = Math_Arrays::build_array(field, 1 + order);
		const Field_Sinh_Cosh<T> sinh_cosh = std::sinh_cosh(operand[operand_offset]);
		function_sinh[0] = sinh_cosh.sinh();
		function_cosh[0] = sinh_cosh.cosh();
		for (int i{ 1 }; i <= order; ++i)
		{
			function_sinh[i] = function_cosh[i - 1];
			function_cosh[i] = function_sinh[i - 1];
		}

		// apply function composition
		compose(operand, operand_offset, function_sinh, sinh, sinh_offset);
		compose(operand, operand_offset, function_cosh, cosh, cosh_offset);
	}

	/** Compute hyperbolic tangent of a derivative structure.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param result array where result must be stored (for
	 * hyperbolic tangent the result array <em>cannot</em> be the input
	 * array)
	 * @param result_offset offset of the result in its array
	 */
	void tanh(const std::vector<double>& operand, const int& operand_offset, std::vector<double>& result, const int& result_offset)
	{
		// create the function value and derivatives
		const auto function = std::vector<double>(1 + order);
		const double t = std::tanh(operand[operand_offset]);
		function[0] = t;

		if (order > 0)
		{
			// the nth order derivative of tanh has the form:
			// dn(tanh(x)/dxn = P_n(tanh(x))
			// where P_n(t) is a degree n+1 polynomial with same parity as n+1
			// P_0(t) = t, P_1(t) = 1 - t^2, P_2(t) = -2 t (1 - t^2) ...
			// the general recurrence relation for P_n is:
			// P_n(x) = (1-t^2) P_(n-1)'(t)
			// as per polynomial parity, we can store coefficients of both P_(n-1) and P_n in the same array
			auto p = std::vector<double>(order + 2);
			p[1] = 1;
			const double t2{ t * t };
			for (const int n{ 1 }; n <= order; ++n)
			{
				// update and evaluate polynomial P_n(t)
				double v{};
				p[n + 1] = -n * p[n];
				for (int k = n + 1; k >= 0; k -= 2)
				{
					v = v * t2 + p[k];
					if (k > 2)
					{
						p[k - 2] = (k - 1) * p[k - 1] - (k - 3) * p[k - 3];
					}
					else if (k == 2)
					{
						p[0] = p[1];
					}
				}
				if ((n & 0x1) == 0)
				{
					v *= t;
				}

				function[n] = v;
			}
		}

		// apply function composition
		compose(operand, operand_offset, function, result, result_offset);
	}

	/** Compute hyperbolic tangent of a derivative structure.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param result array where result must be stored (for
	 * hyperbolic tangent the result array <em>cannot</em> be the input
	 * array)
	 * @param result_offset offset of the result in its array
	 * @param <T> the type of the function parameters and value
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	void tanh(const std::vector<T>& operand, const int& operand_offset, std::vector<T>& result, const int& result_offset)
	{
		const Field<T> field = operand[operand_offset].get_field();

		// create the function value and derivatives
		std::vector<T> function = Math_Arrays::build_array(field, 1 + order);
		const T t = operand[operand_offset].tanh();
		function[0] = t;

		if (order > 0)
		{
			// the nth order derivative of tanh has the form:
			// dn(tanh(x)/dxn = P_n(tanh(x))
			// where P_n(t) is a degree n+1 polynomial with same parity as n+1
			// P_0(t) = t, P_1(t) = 1 - t^2, P_2(t) = -2 t (1 - t^2) ...
			// the general recurrence relation for P_n is:
			// P_n(x) = (1-t^2) P_(n-1)'(t)
			// as per polynomial parity, we can store coefficients of both P_(n-1) and P_n in the same array
			const std::vector<T> p = Math_Arrays::build_array(field, order + 2);
			p[1] = field.get_one();
			const T t2 = t.multiply(t);
			for (const int n{ 1 }; n <= order; ++n)
			{
				// update and evaluate polynomial P_n(t)
				T v = field.get_zero();
				p[n + 1] = p[n].multiply(-n);
				for (int k = n + 1; k >= 0; k -= 2)
				{
					v = v.multiply(t2).add(p[k]);
					if (k > 2)
					{
						p[k - 2] = p[k - 1].multiply(k - 1).subtract(p[k - 3].multiply(k - 3));
					}
					else if (k == 2)
					{
						p[0] = p[1];
					}
				}
				if ((n & 0x1) == 0)
				{
					v = v.multiply(t);
				}

				function[n] = v;
			}
		}

		// apply function composition
		compose(operand, operand_offset, function, result, result_offset);
	}

	/** Compute inverse hyperbolic cosine of a derivative structure.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param result array where result must be stored (for
	 * inverse hyperbolic cosine the result array <em>cannot</em> be the input
	 * array)
	 * @param result_offset offset of the result in its array
	 */
	void acosh(const std::vector<double>& operand, const int& operand_offset, std::vector<double>& result, const int& result_offset)
	{
		// create the function value and derivatives
		std::vector<double> function = std::vector<double>(1 + order);
		const double x = operand[operand_offset];
		function[0] = std::acosh(x);
		if (order > 0)
		{
			// the nth order derivative of acosh has the form:
			// dn(acosh(x)/dxn = P_n(x) / [x^2 - 1]^((2n-1)/2)
			// where P_n(x) is a degree n-1 polynomial with same parity as n-1
			// P_1(x) = 1, P_2(x) = -x, P_3(x) = 2x^2 + 1 ...
			// the general recurrence relation for P_n is:
			// P_n(x) = (x^2-1) P_(n-1)'(x) - (2n-3) x P_(n-1)(x)
			// as per polynomial parity, we can store coefficients of both P_(n-1) and P_n in the same array
			auto p = std::vector<double>(order);
			p[0] = 1;
			const double x2 = x * x;
			const double f = 1.0 / (x2 - 1);
			double coeff = std::sqrt(f);
			function[1] = coeff * p[0];
			for (const int& n = 2; n <= order; ++n)
			{
				// update and evaluate polynomial P_n(x)
				double v{};
				p[n - 1] = (1 - n) * p[n - 2];
				for (int k{ n - 1 }; k >= 0; k -= 2)
				{
					v = v * x2 + p[k];
					if (k > 2)
					{
						p[k - 2] = (1 - k) * p[k - 1] + (k - 2 * n) * p[k - 3];
					}
					else if (k == 2)
					{
						p[0] = -p[1];
					}
				}
				if ((n & 0x1) == 0)
				{
					v *= x;
				}

				coeff *= f;
				function[n] = coeff * v;
			}
		}

		// apply function composition
		compose(operand, operand_offset, function, result, result_offset);
	}

	/** Compute inverse hyperbolic cosine of a derivative structure.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param result array where result must be stored (for
	 * inverse hyperbolic cosine the result array <em>cannot</em> be the input
	 * array)
	 * @param result_offset offset of the result in its array
	 * @param <T> the type of the function parameters and value
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	void acosh(const std::vector<T>& operand, const int& operand_offset, std::vector<T>& result, const int& result_offset)
	{
		const Field<T> field = operand[operand_offset].get_field();

		// create the function value and derivatives
		std::vector<T> function = Math_Arrays::build_array(field, 1 + order);
		const T x = operand[operand_offset];
		function[0] = x.acosh();
		if (order > 0)
		{
			// the nth order derivative of acosh has the form:
			// dn(acosh(x)/dxn = P_n(x) / [x^2 - 1]^((2n-1)/2)
			// where P_n(x) is a degree n-1 polynomial with same parity as n-1
			// P_1(x) = 1, P_2(x) = -x, P_3(x) = 2x^2 + 1 ...
			// the general recurrence relation for P_n is:
			// P_n(x) = (x^2-1) P_(n-1)'(x) - (2n-3) x P_(n-1)(x)
			// as per polynomial parity, we can store coefficients of both P_(n-1) and P_n in the same array
			const std::vector<T> p = Math_Arrays::build_array(field, order);
			p[0] = field.get_one();
			const T x2 = x.multiply(x);
			const T f = x2.subtract(1).reciprocal();
			T coeff = f.sqrt();
			function[1] = coeff.multiply(p[0]);
			for (const int& n = 2; n <= order; ++n)
			{
				// update and evaluate polynomial P_n(x)
				T v = field.get_zero();
				p[n - 1] = p[n - 2].multiply(1 - n);
				for (int k{ n - 1 }; k >= 0; k -= 2)
				{
					v = v.multiply(x2).add(p[k]);
					if (k > 2)
					{
						p[k - 2] = p[k - 1].multiply(1 - k).add(p[k - 3].multiply(k - 2 * n));
					}
					else if (k == 2)
					{
						p[0] = p[1].negate();
					}
				}
				if ((n & 0x1) == 0)
				{
					v = v.multiply(x);
				}

				coeff = coeff.multiply(f);
				function[n] = coeff.multiply(v);
			}
		}

		// apply function composition
		compose(operand, operand_offset, function, result, result_offset);
	}

	/** Compute inverse hyperbolic sine of a derivative structure.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param result array where result must be stored (for
	 * inverse hyperbolic sine the result array <em>cannot</em> be the input
	 * array)
	 * @param result_offset offset of the result in its array
	 */
	void asinh(const std::vector<double>& operand, const int& operand_offset, std::vector<double>& result, const int& result_offset)
	{
		// create the function value and derivatives
		std::vector<double> function = std::vector<double>(1 + order);
		const double x = operand[operand_offset];
		function[0] = std::asinh(x);
		if (order > 0)
		{
			// the nth order derivative of asinh has the form:
			// dn(asinh(x)/dxn = P_n(x) / [x^2 + 1]^((2n-1)/2)
			// where P_n(x) is a degree n-1 polynomial with same parity as n-1
			// P_1(x) = 1, P_2(x) = -x, P_3(x) = 2x^2 - 1 ...
			// the general recurrence relation for P_n is:
			// P_n(x) = (x^2+1) P_(n-1)'(x) - (2n-3) x P_(n-1)(x)
			// as per polynomial parity, we can store coefficients of both P_(n-1) and P_n in the same array
			auto p = std::vector<double>(order);
			p[0] = 1;
			const double x2 = x * x;
			const double f = 1.0 / (1 + x2);
			double coeff = std::sqrt(f);
			function[1] = coeff * p[0];
			for (const int& n = 2; n <= order; ++n)
			{
				// update and evaluate polynomial P_n(x)
				double v{};
				p[n - 1] = (1 - n) * p[n - 2];
				for (int k{ n - 1 }; k >= 0; k -= 2)
				{
					v = v * x2 + p[k];
					if (k > 2)
					{
						p[k - 2] = (k - 1) * p[k - 1] + (k - 2 * n) * p[k - 3];
					}
					else if (k == 2)
					{
						p[0] = p[1];
					}
				}
				if ((n & 0x1) == 0)
				{
					v *= x;
				}

				coeff *= f;
				function[n] = coeff * v;
			}
		}

		// apply function composition
		compose(operand, operand_offset, function, result, result_offset);
	}

	/** Compute inverse hyperbolic sine of a derivative structure.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param result array where result must be stored (for
	 * inverse hyperbolic sine the result array <em>cannot</em> be the input
	 * array)
	 * @param result_offset offset of the result in its array
	 * @param <T> the type of the function parameters and value
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	void asinh(const std::vector<T>& operand, const int& operand_offset, std::vector<T>& result, const int& result_offset)
	{
		const Field<T> field = operand[operand_offset].get_field();

		// create the function value and derivatives
		std::vector<T> function = Math_Arrays::build_array(field, 1 + order);
		const T x = operand[operand_offset];
		function[0] = x.asinh();
		if (order > 0)
		{
			// the nth order derivative of asinh has the form:
			// dn(asinh(x)/dxn = P_n(x) / [x^2 + 1]^((2n-1)/2)
			// where P_n(x) is a degree n-1 polynomial with same parity as n-1
			// P_1(x) = 1, P_2(x) = -x, P_3(x) = 2x^2 - 1 ...
			// the general recurrence relation for P_n is:
			// P_n(x) = (x^2+1) P_(n-1)'(x) - (2n-3) x P_(n-1)(x)
			// as per polynomial parity, we can store coefficients of both P_(n-1) and P_n in the same array
			const std::vector<T> p = Math_Arrays::build_array(field, order);
			p[0] = field.get_one();
			const T x2 = x.multiply(x);
			const T f = x2.add(1).reciprocal();
			T coeff = f.sqrt();
			function[1] = coeff.multiply(p[0]);
			for (const int& n = 2; n <= order; ++n)
			{
				// update and evaluate polynomial P_n(x)
				T v = field.get_zero();
				p[n - 1] = p[n - 2].multiply(1 - n);
				for (int k{ n - 1 }; k >= 0; k -= 2)
				{
					v = v.multiply(x2).add(p[k]);
					if (k > 2)
					{
						p[k - 2] = p[k - 1].multiply(k - 1).add(p[k - 3].multiply(k - 2 * n));
					}
					else if (k == 2)
					{
						p[0] = p[1];
					}
				}
				if ((n & 0x1) == 0)
				{
					v = v.multiply(x);
				}

				coeff = coeff.multiply(f);
				function[n] = coeff.multiply(v);
			}
		}

		// apply function composition
		compose(operand, operand_offset, function, result, result_offset);
	}

	/** Compute inverse hyperbolic tangent of a derivative structure.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param result array where result must be stored (for
	 * inverse hyperbolic tangent the result array <em>cannot</em> be the input
	 * array)
	 * @param result_offset offset of the result in its array
	 */
	void atanh(const std::vector<double>& operand, const int& operand_offset, std::vector<double>& result, const int& result_offset)
	{
		// create the function value and derivatives
		std::vector<double> function = std::vector<double>(1 + order);
		const double x = operand[operand_offset];
		function[0] = std::atanh(x);
		if (order > 0)
		{
			// the nth order derivative of atanh has the form:
			// dn(atanh(x)/dxn = Q_n(x) / (1 - x^2)^n
			// where Q_n(x) is a degree n-1 polynomial with same parity as n-1
			// Q_1(x) = 1, Q_2(x) = 2x, Q_3(x) = 6x^2 + 2 ...
			// the general recurrence relation for Q_n is:
			// Q_n(x) = (1-x^2) Q_(n-1)'(x) + 2(n-1) x Q_(n-1)(x)
			// as per polynomial parity, we can store coefficients of both Q_(n-1) and Q_n in the same array
			const std::vector<double> q = std::vector<double>(order);
			q[0] = 1;
			const double x2 = x * x;
			const double f = 1.0 / (1 - x2);
			double coeff = f;
			function[1] = coeff * q[0];
			for (const int& n = 2; n <= order; ++n)
			{
				// update and evaluate polynomial Q_n(x)
				double v{};
				q[n - 1] = n * q[n - 2];
				for (int k{ n - 1 }; k >= 0; k -= 2)
				{
					v = v * x2 + q[k];
					if (k > 2)
					{
						q[k - 2] = (k - 1) * q[k - 1] + (2 * n - k + 1) * q[k - 3];
					}
					else if (k == 2)
					{
						q[0] = q[1];
					}
				}
				if ((n & 0x1) == 0)
				{
					v *= x;
				}

				coeff *= f;
				function[n] = coeff * v;
			}
		}

		// apply function composition
		compose(operand, operand_offset, function, result, result_offset);
	}

	/** Compute inverse hyperbolic tangent of a derivative structure.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param result array where result must be stored (for
	 * inverse hyperbolic tangent the result array <em>cannot</em> be the input
	 * array)
	 * @param result_offset offset of the result in its array
	 * @param <T> the type of the function parameters and value
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	void atanh(const std::vector<T>& operand, const int& operand_offset, std::vector<T>& result, const int& result_offset)
	{
		const Field<T> field = operand[operand_offset].get_field();

		// create the function value and derivatives
		std::vector<T> function = Math_Arrays::build_array(field, 1 + order);
		const T x = operand[operand_offset];
		function[0] = x.atanh();
		if (order > 0)
		{
			// the nth order derivative of atanh has the form:
			// dn(atanh(x)/dxn = Q_n(x) / (1 - x^2)^n
			// where Q_n(x) is a degree n-1 polynomial with same parity as n-1
			// Q_1(x) = 1, Q_2(x) = 2x, Q_3(x) = 6x^2 + 2 ...
			// the general recurrence relation for Q_n is:
			// Q_n(x) = (1-x^2) Q_(n-1)'(x) + 2(n-1) x Q_(n-1)(x)
			// as per polynomial parity, we can store coefficients of both Q_(n-1) and Q_n in the same array
			const std::vector<T> q = Math_Arrays::build_array(field, order);
			q[0] = field.get_one();
			const T x2 = x.multiply(x);
			const T f = x2.subtract(1).negate().reciprocal();
			T coeff = f;
			function[1] = coeff.multiply(q[0]);
			for (const int& n = 2; n <= order; ++n)
			{
				// update and evaluate polynomial Q_n(x)
				T v = field.get_zero();
				q[n - 1] = q[n - 2].multiply(n);
				for (int k{ n - 1 }; k >= 0; k -= 2)
				{
					v = v.multiply(x2).add(q[k]);
					if (k > 2)
					{
						q[k - 2] = q[k - 1].multiply(k - 1).add(q[k - 3].multiply(2 * n - k + 1));
					}
					else if (k == 2)
					{
						q[0] = q[1];
					}
				}
				if ((n & 0x1) == 0)
				{
					v = v.multiply(x);
				}

				coeff = coeff.multiply(f);
				function[n] = coeff.multiply(v);
			}
		}

		// apply function composition
		compose(operand, operand_offset, function, result, result_offset);
	}

	/** Compute composition of a derivative structure by a function.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param f array of value and derivatives of the function at
	 * the current point (i.e. at {@code operand[operand_offset]}).
	 * @param result array where result must be stored (for
	 * composition the result array <em>cannot</em> be the input
	 * array)
	 * @param result_offset offset of the result in its array
	 */
	void compose(const std::vector<double>& operand, const int& operand_offset, const std::vector<double>& f, std::vector<double>& result, const int& result_offset)
	{
		for (int i{}; i < comp_indirection.size(); ++i)
		{
			const std::vector<std::vector<int>> mapping_i = comp_indirection[i];
			double r{};
			for (int j{}; j < mapping_i.size(); ++j)
			{
				const std::vector<int> mapping_i_j = mapping_i[j];
				double product = mapping_i_j[0] * f[mapping_i_j[1]];
				for (int k = 2; k < mapping_i_j.size(); ++k)
				{
					product *= operand[operand_offset + mapping_i_j[k]];
				}
				r += product;
			}
			result[result_offset + i] = r;
		}
	}

	/** Compute composition of a derivative structure by a function.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param f array of value and derivatives of the function at
	 * the current point (i.e. at {@code operand[operand_offset]}).
	 * @param result array where result must be stored (for
	 * composition the result array <em>cannot</em> be the input
	 * array)
	 * @param result_offset offset of the result in its array
	 * @param <T> the type of the function parameters and value
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	void compose(const std::vector<T>& operand, const int& operand_offset, const std::vector<T> f, std::vector<T>& result, const int& result_offset)
	{
		const T zero = f[0].get_field().get_zero();
		for (int i{}; i < comp_indirection.size(); ++i)
		{
			const std::vector<std::vector<int>> mapping_i = comp_indirection[i];
			T r = zero;
			for (int j{}; j < mapping_i.size(); ++j)
			{
				const std::vector<int> mapping_i_j = mapping_i[j];
				T product = f[mapping_i_j[1]].multiply(mapping_i_j[0]);
				for (int k = 2; k < mapping_i_j.size(); ++k)
				{
					product = product.multiply(operand[operand_offset + mapping_i_j[k]]);
				}
				r = r.add(product);
			}
			result[result_offset + i] = r;
		}
	}

	/** Compute composition of a derivative structure by a function.
	 * @param operand array holding the operand
	 * @param operand_offset offset of the operand in its array
	 * @param f array of value and derivatives of the function at
	 * the current point (i.e. at {@code operand[operand_offset]}).
	 * @param result array where result must be stored (for
	 * composition the result array <em>cannot</em> be the input
	 * array)
	 * @param result_offset offset of the result in its array
	 * @param <T> the type of the function parameters and value
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	void compose(const std::vector<T>& operand, const int& operand_offset, const std::vector<double>& f, std::vector<T>& result, const int& result_offset)
	{
		const T zero = operand[operand_offset].get_field().get_zero();
		for (int i{}; i < comp_indirection.size(); ++i)
		{
			const std::vector<std::vector<int>> mapping_i = comp_indirection[i];
			T r = zero;
			for (int j{}; j < mapping_i.size(); ++j)
			{
				const std::vector<int> mapping_i_j = mapping_i[j];
				T product = zero.add(f[mapping_i_j[1]] * mapping_i_j[0]);
				for (int k = 2; k < mapping_i_j.size(); ++k)
				{
					product = product.multiply(operand[operand_offset + mapping_i_j[k]]);
				}
				r = r.add(product);
			}
			result[result_offset + i] = r;
		}
	}

	/** Evaluate Taylor expansion of a derivative structure.
	 * @param ds array holding the derivative structure
	 * @param ds_offset offset of the derivative structure in its array
	 * @param delta parameters offsets (&Delta;x, &Delta;y, ...)
	 * @return value of the Taylor expansion at x + &Delta;x, y + &Delta;y, ...
	 * @Math_Runtime_Exception if factorials becomes too large
	 */
	double taylor(const std::vector<double> ds, const int& ds_offset, const double ... delta)
		Math_Runtime_Exception
	{
		double value = 0;
		for (int i = get_size() - 1; i >= 0; --i)
		{
			const std::vector<int> orders = derivatives_indirection[i];
			double term = ds[ds_offset + i];
			for (int k{}; k < orders.size(); ++k)
			{
				if (orders[k] > 0)
				{
					term *= std::pow(delta[k], orders[k]) /
						Combinatorics_Utils.factorial(orders[k]);
				}
			}
			value += term;
		}
		return value;
	}

	/** Evaluate Taylor expansion of a derivative structure.
	 * @param ds array holding the derivative structure
	 * @param ds_offset offset of the derivative structure in its array
	 * @param delta parameters offsets (&Delta;x, &Delta;y, ...)
	 * @return value of the Taylor expansion at x + &Delta;x, y + &Delta;y, ...
	 * @Math_Runtime_Exception if factorials becomes too large
	 * @param <T> the type of the function parameters and value
	 */
	 //@Safe_Varargs
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	const  T taylor(const std::vector<T> ds, const int& ds_offset, const T ... delta)
		Math_Runtime_Exception
	{
		const Field<T> field = ds[ds_offset].get_field();
		T value = field.get_zero();
		for (int i = get_size() - 1; i >= 0; --i)
		{
			const std::vector<int> orders = derivatives_indirection[i];
			T term = ds[ds_offset + i];
			for (int k{}; k < orders.size(); ++k)
			{
				if (orders[k] > 0)
				{
					term = term.multiply(delta[k].pow(orders[k]).
						divide(Combinatorics_Utils.factorial(orders[k])));
				}
			}
			value = value.add(term);
		}
		return value;
	}

	/** Evaluate Taylor expansion of a derivative structure.
	 * @param ds array holding the derivative structure
	 * @param ds_offset offset of the derivative structure in its array
	 * @param delta parameters offsets (&Delta;x, &Delta;y, ...)
	 * @return value of the Taylor expansion at x + &Delta;x, y + &Delta;y, ...
	 * @Math_Runtime_Exception if factorials becomes too large
	 * @param <T> the type of the function parameters and value
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	T taylor(const std::vector<T> ds, const int& ds_offset, const double ... delta)
		Math_Runtime_Exception
	{
		const Field<T> field = ds[ds_offset].get_field();
		T value = field.get_zero();
		for (int i = get_size() - 1; i >= 0; --i)
		{
			const std::vector<int> orders = derivatives_indirection[i];
			T term = ds[ds_offset + i];
			for (int k{}; k < orders.size(); ++k)
			{
				if (orders[k] > 0)
				{
					term = term.multiply(field.get_zero().new_instance(delta[k]).pow(orders[k]).
						divide(Combinatorics_Utils.factorial(orders[k])));
				}
			}
			value = value.add(term);
		}
		return value;
	}

	/** Check rules set compatibility.
	 * @param compiler other compiler to check against instance
	 * @exception  if number of free parameters or orders are inconsistent
	 */
	void check_compatibility(const DS_Compiler compiler)
	{
		Math_Utils::check_dimension(parameters, compiler.parameters);
		Math_Utils::check_dimension(order, compiler.order);
	}
};