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
#include "MatrixUtils.h"
//import java.util.Arrays;

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.util.FastMath;
#include <vector>
#include "RealMatrix.h"

/**
 * Class transforming a symmetrical matrix to tridiagonal shape.
 * <p>A symmetrical m &times; m matrix A can be written as the product of three matrices:
 * A = Q &times; T &times; Q<sup>T</sup> with Q an orthogonal matrix and T a symmetrical
 * tridiagonal matrix. Both Q and T are m &times; m matrices.</p>
 * <p>This implementation only uses the upper part of the matrix, the part below the
 * diagonal is not accessed at all.</p>
 * <p>Transformation to tridiagonal shape is often not a goal by itself, but it is
 * an intermediate step in more general decomposition algorithms like {@link
 * Eigen_Decomposition eigen decomposition}. This class is therefore intended for internal
 * use by the library and is not public. As a consequence of this explicitly limited scope, * many methods directly returns references to internal arrays, not copies.</p>
 */
class Tri_Diagonal_Transformer
{
private:
	/** Householder vectors. */
	std::vector<std::vector<double>> my_householder_vectors;
	/** Main diagonal. */
	std::vector<double> my_main;
	/** Secondary diagonal. */
	std::vector<double> my_secondary;
	/** Cached value of Q. */
	Real_Matrix my_cached_q;
	/** Cached value of Qt. */
	Real_Matrix my_cached_qt;
	/** Cached value of T. */
	Real_Matrix my_cached_t;

	/**
	 * Transform original matrix to tridiagonal form.
	 * <p>Transformation is done using Householder transforms.</p>
	 */
	void transform()
	{
		const int m = my_householder_vectors.size();
		const std::vector<double> z = std::vector<double>(m];
		for (int k{}; k < m - 1; k++)
		{
			//zero-out a row and a column simultaneously
			auto h_k = my_householder_vectors[k];
			my_main[k] = h_k[k];
			double x_norm_sqr{};
			for (int j{ k + 1 }; j < m; ++j)
			{
				const double c = h_k[j];
				x_norm_sqr += c * c;
			}
			const double& a = (h_k[k + 1] > 0) ? -std::sqrt(x_norm_sqr) : std::sqrt(x_norm_sqr);
			my_secondary[k] = a;
			if (a != 0.0)
			{
				// apply Householder transform from left and right simultaneously

				h_k[k + 1] -= a;
				const double beta = -1 / (a * h_k[k + 1]);

				// compute a = beta A v, where v is the Householder vector
				// this loop is written in such a way
				//   1) only the upper triangular part of the matrix is accessed
				//   2) access is cache-friendly for a matrix stored in rows
				Arrays.fill(z, k + 1, m, 0);
				for (int i = k + 1; i < m; ++i)
				{
					const std::vector<double> h_i = my_householder_vectors[i];
					const double h_k_i = h_k[i];
					double zI = h_i[i] * h_k_i;
					for (int j = i + 1; j < m; ++j)
					{
						const double h_i_j = h_i[j];
						zI += h_i_j * h_k[j];
						z[j] += h_i_j * h_k_i;
					}
					z[i] = beta * (z[i] + zI);
				}

				// compute gamma = beta vT z / 2
				double gamma = 0;
				for (int i = k + 1; i < m; ++i)
				{
					gamma += z[i] * h_k[i];
				}
				gamma *= beta / 2;

				// compute z = z - gamma v
				for (int i = k + 1; i < m; ++i)
				{
					z[i] -= gamma * h_k[i];
				}

				// update matrix: A = A - v zT - z vT
				// only the upper triangular part of the matrix is updated
				for (int i = k + 1; i < m; ++i)
				{
					auto h_i = my_householder_vectors[i];
					for (int j = i; j < m; ++j)
					{
						h_i[j] -= h_k[i] * z[j] + z[i] * h_k[j];
					}
				}
			}
		}
		my_main[m - 1] = my_householder_vectors[m - 1][m - 1];
	}

public:
	/**
	 * Build the transformation to tridiagonal shape of a symmetrical matrix.
	 * <p>The specified matrix is assumed to be symmetrical without any check.
	 * Only the upper triangular part of the matrix is used.</p>
	 *
	 * @param matrix Symmetrical matrix to transform.
	 * @ if the matrix is not square.
	 */
	Tri_Diagonal_Transformer(Real_Matrix matrix)
	{
		if (!matrix.is_square())
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NON_SQUARE_MATRIX, matrix.get_row_dimension(), matrix.get_column_dimension());
		}

		const int m = matrix.get_row_dimension();
		my_householder_vectors = matrix.get_data();
		my_main = std::vector<double>(m);
		my_secondary = std::vector<double>(m - 1);
		my_cached_q = NULL;
		my_cached_qt = NULL;
		my_cached_t = NULL;

		// transform matrix
		transform();
	}

	/**
	 * Returns the matrix Q of the transform.
	 * <p>Q is an orthogonal matrix, i.e. its transpose is also its inverse.</p>
	 * @return the Q matrix
	 */
	Real_Matrix get_q()
	{
		if (my_cached_q == NULL)
		{
			my_cached_q = get_q_t().transpose();
		}
		return my_cached_q;
	}

	/**
	 * Returns the transpose of the matrix Q of the transform.
	 * <p>Q is an orthogonal matrix, i.e. its transpose is also its inverse.</p>
	 * @return the Q matrix
	 */
	Real_Matrix get_q_t()
	{
		if (my_cached_qt == NULL)
		{
			const int m = my_householder_vectors.size();
			auto qta = std::vector<std::vector<double>>(m, std::vector<double>(m));

			// build up first part of the matrix by applying Householder transforms
			for (int k = m - 1; k >= 1; --k)
			{
				const std::vector<double> h_k = my_householder_vectors[k - 1];
				qta[k][k] = 1;
				if (h_k[k] != 0.0)
				{
					const double inv = 1.0 / (my_secondary[k - 1] * h_k[k]);
					double beta = 1.0 / my_secondary[k - 1];
					qta[k][k] = 1 + beta * h_k[k];
					for (int i = k + 1; i < m; ++i)
					{
						qta[k][i] = beta * h_k[i];
					}
					for (int j{ k + 1 }; j < m; ++j)
					{
						beta = 0;
						for (int i = k + 1; i < m; ++i)
						{
							beta += qta[j][i] * h_k[i];
						}
						beta *= inv;
						qta[j][k] = beta * h_k[k];
						for (int i = k + 1; i < m; ++i)
						{
							qta[j][i] += beta * h_k[i];
						}
					}
				}
			}
			qta[0][0] = 1;
			my_cached_qt = Matrix_Utils::create_real_matrix(qta);
		}

		// return the cached matrix
		return my_cached_qt;
	}

	/**
	 * Returns the tridiagonal matrix T of the transform.
	 * @return the T matrix
	 */
	Real_Matrix get_t()
	{
		if (my_cached_t == NULL)
		{
			const int m = my_main.size();
			auto ta = std::vector<std::vector<double>>(m, std::vector<double>(m));
			for (int i{}; i < m; ++i)
			{
				ta[i][i] = my_main[i];
				if (i > 0)
				{
					ta[i][i - 1] = my_secondary[i - 1];
				}
				if (i < my_main.size() - 1)
				{
					ta[i][i + 1] = my_secondary[i];
				}
			}
			my_cached_t = Matrix_Utils::create_real_matrix(ta);
		}

		// return the cached matrix
		return my_cached_t;
	}

	/**
	 * Get the Householder vectors of the transform.
	 * <p>Note that since this class is only intended for internal use, * it returns directly a reference to its internal arrays, not a copy.</p>
	 * @return the main diagonal elements of the B matrix
	 */
	std::vector<std::vector<double>> get_householder_vectors_ref() const
	{
		return my_householder_vectors; // NOPMD - returning an internal array is intentional and documented here
	}

	/**
	 * Get the main diagonal elements of the matrix T of the transform.
	 * <p>Note that since this class is only intended for internal use, * it returns directly a reference to its internal arrays, not a copy.</p>
	 * @return the main diagonal elements of the T matrix
	 */
	std::vector<double> get_main_diagonal_ref() const
	{
		return my_main; // NOPMD - returning an internal array is intentional and documented here
	}

	/**
	 * Get the secondary diagonal elements of the matrix T of the transform.
	 * <p>Note that since this class is only intended for internal use, * it returns directly a reference to its internal arrays, not a copy.</p>
	 * @return the secondary diagonal elements of the T matrix
	 */
	std::vector<double> get_secondary_diagonal_ref() const
	{
		return my_secondary; // NOPMD - returning an internal array is intentional and documented here
	}
};