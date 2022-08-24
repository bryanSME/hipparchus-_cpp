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

  //package org.hipparchus.random;

  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.linear.Real_Matrix;
  //import org.hipparchus.linear.RectangularCholesky_Decomposition;

  /**
   * A {@link Random_Vector_Generator} that generates vectors with with
   * correlated components.
   * <p>
   * Random vectors with correlated components are built by combining
   * the uncorrelated components of another random vector in such a way that
   * the resulting correlations are the ones specified by a positive
   * definite covariance matrix.
   * <p>
   * The main use for correlated random vector generation is for Monte-_Carlo
   * simulation of physical problems with several variables, for example to
   * generate error vectors to be added to a nominal vector. A particularly
   * interesting case is when the generated vector should be drawn from a <a
   * href="http://en.wikipedia.org/wiki/Multivariate_normal_distribution">
   * Multivariate Normal Distribution</a>. The approach using a Cholesky
   * decomposition is quite usual in this case. However, it can be extended
   * to other cases as long as the underlying random generator provides
   * {@link Normalized_Random_Generator normalized values} like {@link
   * Gaussian_randomGenerator} or {@link UniformRandom_Generator}.
   * <p>
   * Sometimes, the covariance matrix for a given simulation is not
   * strictly positive definite. This means that the correlations are
   * not all independent from each other. In this case, however, the non
   * strictly positive elements found during the Cholesky decomposition
   * of the covariance matrix should not be negative either, they
   * should be NULL. Another non-conventional extension handling this case
   * is used here. Rather than computing <code>C = U<sup>T</sup>.U</code>
   * where <code>C</code> is the covariance matrix and <code>U</code>
   * is an upper-triangular matrix, we compute <code>C = B.B<sup>T</sup></code>
   * where <code>B</code> is a rectangular matrix having
   * more rows than columns. The number of columns of <code>B</code> is
   * the rank of the covariance matrix, and it is the dimension of the
   * uncorrelated random vector that is needed to compute the component
   * of the correlated vector. This class handles this situation
   * automatically.
   */
class CorrelatedRandom_Vector_Generator
	: 
	public Random_Vector_Generator
{
private:
	/** Mean vector. */
	std::vector<double> my_mean;
	/** Underlying generator. */
	Normalized_Random_Generator my_generator;
	/** Storage for the normalized vector. */
	std::vector<double> my_normalized;
	/** Root of the covariance matrix. */
	Real_Matrix my_root;

public:
	/**
	 * Builds a correlated random vector generator from its mean
	 * vector and covariance matrix.
	 *
	 * @param mean Expected mean values for all components.
	 * @param covariance Covariance matrix.
	 * @param small Diagonal elements threshold under which  column are
	 * considered to be dependent on previous ones and are discarded
	 * @param generator underlying generator for uncorrelated normalized
	 * components.
	 * @org.hipparchus.exception.
	 * if the covariance matrix is not strictly positive definite.
	 * @ if the mean and covariance
	 * arrays dimensions do not match.
	 */
	CorrelatedRandom_Vector_Generator(const std::vector<double>& mean, const Real_Matrix& covariance, const double& small, const Normalized_Random_Generator& generator)
	{
		int order = covariance.get_row_dimension();
		if (mean.size() != order)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, mean.size(), order);
		}
		my_mean = mean.clone();

		const auto decomposition = RectangularCholesky_Decomposition(covariance, small);
		my_root = decomposition.get_root_matrix();

		my_generator = generator;
		my_normalized = std::vector<double>(decomposition.get_rank());
	}

	/**
	 * Builds a NULL mean random correlated vector generator from its
	 * covariance matrix.
	 *
	 * @param covariance Covariance matrix.
	 * @param small Diagonal elements threshold under which  column are
	 * considered to be dependent on previous ones and are discarded.
	 * @param generator Underlying generator for uncorrelated normalized
	 * components.
	 * @org.hipparchus.exception.
	 * if the covariance matrix is not strictly positive definite.
	 */
	CorrelatedRandom_Vector_Generator(const Real_Matrix& covariance, const double& small, const Normalized_Random_Generator& generator)
	{
		int order = covariance.get_row_dimension();
		mean = std::vector<double>(order];
		for (int i{}; i < order; ++i)
		{
			mean[i] = 0;
		}

		const auto decomposition = RectangularCholesky_Decomposition(covariance, small);
		my_root = decomposition.get_root_matrix();

		my_generator = generator;
		my_normalized = std::vector<double>(decomposition.get_rank());
	}

	/** Get the underlying normalized components generator.
	 * @return underlying uncorrelated components generator
	 */
	Normalized_Random_Generator get_generator() const
	{
		return my_generator;
	}

	/** Get the rank of the covariance matrix.
	 * The rank is the number of independent rows in the covariance
	 * matrix, it is also the number of columns of the root matrix.
	 * @return rank of the square matrix.
	 * @see #get_root_matrix()
	 */
	int get_rank() const
	{
		return my_normalized.size();
	}

	/** Get the root of the covariance matrix.
	 * The root is the rectangular matrix <code>B</code> such that
	 * the covariance matrix is equal to <code>B.B<sup>T</sup></code>
	 * @return root of the square matrix
	 * @see #get_rank()
	 */
	Real_Matrix get_root_matrix() const
	{
		return my_root;
	}

	/** Generate a correlated random vector.
	 * @return a random vector as an array of double. The returned array
	 * is created at each call, the caller can do what it wants with it.
	 */
	 //override
	std::vector<double> next_vector()
	{
		// generate uncorrelated vector
		for (int i{}; i < normalized.size(); ++i)
		{
			my_normalized[i] = my_generator.next_normalized_double();
		}

		// compute correlated vector
		auto correlated = std::vector<double>(mean.size());
		for (int i{}; i < correlated.size(); ++i)
		{
			correlated[i] = mean[i];
			for (int j{}; j < my_root.get_column_dimension(); ++j)
			{
				correlated[i] += my_root.get_entry(i, j) * normalized[j];
			}
		}

		return correlated;
	}
};