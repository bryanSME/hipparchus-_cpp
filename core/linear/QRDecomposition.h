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

  //import java.util.Arrays;

  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.util.FastMath;
#include <vector>
#include "RealMatrix.h"
#include "MatrixUtils.h"

/**
 * Calculates the QR-decomposition of a matrix.
 * <p>The QR-decomposition of a matrix A consists of two matrices Q and R
 * that satisfy: A = QR, Q is orthogonal (Q<sup>T</sup>Q = I), and R is
 * upper triangular. If A is m&times;n, Q is m&times;m and R m&times;n.</p>
 * <p>This class compute the decomposition using Householder reflectors.</p>
 * <p>For efficiency purposes, the decomposition in packed form is transposed.
 * This allows inner loop to iterate inside rows, which is much more cache-efficient
 * in Java.</p>
 * <p>This class is based on the class with similar name from the
 * <a href="http://math.nist.gov/javanumerics/jama/">JAMA</a> library, with the
 * following changes:</p>
 * <ul>
 *   <li>a {@link #get_q_t() get_q_t} method has been added,</li>
 *   <li>the {@code solve} and {@code is_full_rank} methods have been replaced
 *   by a {@link #get_solver() get_solver} method and the equivalent methods
 *   provided by the returned {@link Decomposition_Solver}.</li>
 * </ul>
 *
 * @see <a href="http://mathworld.wolfram.com/QR_Decomposition.html">MathWorld</a>
 * @see <a href="http://en.wikipedia.org/wiki/QR_decomposition">Wikipedia</a>
 *
 */
class QR_Decomposition
{
private:
	/**
	 * A packed TRANSPOSED representation of the QR decomposition.
	 * <p>The elements BELOW the diagonal are the elements of the UPPER triangular
	 * matrix R, and the rows ABOVE the diagonal are the Householder reflector vectors
	 * from which an explicit form of Q can be recomputed if desired.</p>
	 */
	std::vector<std::vector<double>> my_qrt;
	/** The diagonal elements of R. */
	std::vector<double> my_r_diag;
	/** Cached value of Q. */
	Real_Matrix my_cached_q;
	/** Cached value of QT. */
	Real_Matrix my_cached_q_t;
	/** Cached value of R. */
	Real_Matrix my_cached_r;
	/** Cached value of H. */
	Real_Matrix my_cached_h;
	/** Singularity my_threshold. */
	const double my_threshold;

protected:
	/** Decompose matrix.
	 * @param matrix transposed matrix
	 */
	void decompose(automatrix)
	{
		for (const int& minor = 0; minor < std::min(matrix.size(), matrix[0].size()); minor++)
		{
			perform_householder_reflection(minor, matrix);
		}
	}

	/** Perform Householder reflection for a minor A(minor, minor) of A.
	 * @param minor minor index
	 * @param matrix transposed matrix
	 */
	void perform_householder_reflection(const int& minor, automatrix)
	{
		auto qrt_minor = matrix[minor];

		/*
		 * Let x be the first column of the minor, and a^2 = |x|^2.
		 * x will be in the positions qr[minor][minor] through qr[m][minor].
		 * The first column of the transformed minor will be (a,0,0,..)'
		 * The sign of a is chosen to be opposite to the sign of the first
		 * component of x. Let's find a:
		 */
		double x_norm_sqr{};
		for (int row = minor; row < qrt_minor.size(); row++)
		{
			const double c = qrt_minor[row];
			x_norm_sqr += c * c;
		}
		const double a = (qrt_minor[minor] > 0) ? -std::sqrt(x_norm_sqr) : std::sqrt(x_norm_sqr);
		my_r_diag[minor] = a;

		if (a != 0.0)
		{
			/*
			 * Calculate the normalized reflection vector v and transform
			 * the first column. We know the norm of v beforehand: v = x-ae
			 * so |v|^2 = <x-ae,x-ae> = <x,x>-2a<x,e>+a^2<e,e> =
			 * a^2+a^2-2a<x,e> = 2a*(a - <x,e>).
			 * Here <x, e> is now qr[minor][minor].
			 * v = x-ae is stored in the column at qr:
			 */
			qrt_minor[minor] -= a; // now |v|^2 = -2a*(qr[minor][minor])

			/*
			 * Transform the rest of the columns of the minor:
			 * They will be transformed by the matrix H = I-2vv'/|v|^2.
			 * If x is a column vector of the minor, then
			 * Hx = (I-2vv'/|v|^2)x = x-2vv'x/|v|^2 = x - 2<x,v>/|v|^2 v.
			 * Therefore the transformation is easily calculated by
			 * subtracting the column vector (2<x,v>/|v|^2)v from x.
			 *
			 * Let 2<x,v>/|v|^2 = alpha. From above we have
			 * |v|^2 = -2a*(qr[minor][minor]), so
			 * alpha = -<x,v>/(a*qr[minor][minor])
			 */
			for (int col{ minor + 1 }; col < matrix.size(); col++)
			{
				auto qrt_col = matrix[col];
				double alpha{};
				for (int row = minor; row < qrt_col.size(); row++)
				{
					alpha -= qrt_col[row] * qrt_minor[row];
				}
				alpha /= a * qrt_minor[minor];

				// Subtract the column vector alpha*v from x.
				for (int row = minor; row < qrt_col.size(); row++)
				{
					qrt_col[row] -= alpha * qrt_minor[row];
				}
			}
		}
	}

public:
	/**
	 * Calculates the QR-decomposition of the given matrix.
	 * The singularity my_threshold defaults to zero.
	 *
	 * @param matrix The matrix to decompose.
	 *
	 * @see #QR_Decomposition(Real_Matrix,double)
	 */
	QR_Decomposition(Real_Matrix matrix) : my_threshold{ 0.0 }
	{
		QR_Decomposition(matrix, 0.0);
	}

	/**
	 * Calculates the QR-decomposition of the given matrix.
	 *
	 * @param matrix The matrix to decompose.
	 * @param my_threshold Singularity my_threshold.
	 */
	QR_Decomposition(const Real_Matrix& matrix, const double& threshold) : my_threshold{ threshold }
	{
		const int m = matrix.get_row_dimension();
		const int n = matrix.get_column_dimension();
		my_qrt = matrix.transpose().get_data();
		my_r_diag = std::vector<double>(std::min(m, n));
		my_cached_q = NULL;
		my_cached_q_t = NULL;
		my_cached_r = NULL;
		my_cached_h = NULL;

		decompose(my_qrt);
	}

	/**
	 * Returns the matrix R of the decomposition.
	 * <p>R is an upper-triangular matrix</p>
	 * @return the R matrix
	 */
	Real_Matrix get_r()
	{
		if (my_cached_r == NULL)
		{
			// R is supposed to be m x n
			const int n = my_qrt.size();
			const int m = my_qrt[0].size();
			auto ra = std::vector<std::vector<double>>(m, std::vector<double>(n));
			// copy the diagonal from my_r_diag and the upper triangle of qr
			for (int row = std::min(m, n) - 1; row >= 0; row--)
			{
				ra[row][row] = my_r_diag[row];
				for (int col = row + 1; col < n; col++)
				{
					ra[row][col] = my_qrt[col][row];
				}
			}
			my_cached_r = Matrix_Utils::create_real_matrix(ra);
		}

		// return the cached matrix
		return my_cached_r;
	}

	/**
	 * Returns the matrix Q of the decomposition.
	 * <p>Q is an orthogonal matrix</p>
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
	 * Returns the transpose of the matrix Q of the decomposition.
	 * <p>Q is an orthogonal matrix</p>
	 * @return the transpose of the Q matrix, Q<sup>T</sup>
	 */
	Real_Matrix get_q_t()
	{
		if (my_cached_q_t == NULL)
		{
			// QT is supposed to be m x m
			const int n = my_qrt.size();
			const int m = my_qrt[0].size();
			auto qta = std::vector<std::vector<double>>(m, std::vector<double>(m));

			/*
			 * Q = Q1 Q2 ... Q_m, so Q is formed by first constructing Q_m and then
			 * applying the Householder transformations Q_(m-1),Q_(m-2),...,Q1 in
			 * succession to the result
			 */
			for (int minor{ m - 1 }; minor >= std::min(m, n); minor--)
			{
				qta[minor][minor] = 1.0;
			}

			for (int minor = std::min(m, n) - 1; minor >= 0; minor--)
			{
				const std::vector<double> qrt_minor = my_qrt[minor];
				qta[minor][minor] = 1.0;
				if (qrt_minor[minor] != 0.0)
				{
					for (int col = minor; col < m; col++)
					{
						double alpha = 0;
						for (int row = minor; row < m; row++)
						{
							alpha -= qta[col][row] * qrt_minor[row];
						}
						alpha /= my_r_diag[minor] * qrt_minor[minor];

						for (int row = minor; row < m; row++)
						{
							qta[col][row] += -alpha * qrt_minor[row];
						}
					}
				}
			}
			my_cached_q_t = Matrix_Utils::create_real_matrix(qta);
		}

		// return the cached matrix
		return my_cached_q_t;
	}

	/**
	 * Returns the Householder reflector vectors.
	 * <p>H is a lower trapezoidal matrix whose columns represent
	 * each successive Householder reflector vector. This matrix is used
	 * to compute Q.</p>
	 * @return a matrix containing the Householder reflector vectors
	 */
	Real_Matrix get_h()
	{
		if (my_cached_h == NULL)
		{
			const int n = my_qrt.size();
			const int m = my_qrt[0].size();
			auto ha = std::vector<std::vector<double>>(m, std::vector<double>(n));
			for (int i{}; i < m; ++i)
			{
				for (int j{}; j < std::min(i + 1, n); ++j)
				{
					ha[i][j] = my_qrt[j][i] / -my_r_diag[j];
				}
			}
			my_cached_h = Matrix_Utils::create_real_matrix(ha);
		}

		// return the cached matrix
		return my_cached_h;
	}

	/**
	 * Get a solver for finding the A &times; X = B solution in least square sense.
	 * <p>
	 * Least Square sense means a solver can be computed for an overdetermined system, * (i.e. a system with more equations than unknowns, which corresponds to a tall A
	 * matrix with more rows than columns). In any case, if the matrix is singular
	 * within the tolerance set at {@link QR_Decomposition#QR_Decomposition(Real_Matrix, * double) construction}, an error will be triggered when
	 * the {@link Decomposition_Solver#solve(Real_Vector) solve} method will be called.
	 * </p>
	 * @return a solver
	 */
	Decomposition_Solver get_solver()
	{
		return Solver();
	}

	/** Specialized solver. */
	private class Solver : Decomposition_Solver
	{
		/** {@inherit_doc} */
		//override
		public bool is_non_singular()
		{
			return !check_singular(my_r_diag, my_threshold, false);
		}

		/** {@inherit_doc} */
		//override
		public Real_Vector solve(Real_Vector b)
		{
			const int n = my_qrt.size();
			const int m = my_qrt[0].size();
			if (b.get_dimension() != m)
			{
				throw std::exception("not implemented");
				//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, b.get_dimension(), m);
			}
			check_singular(my_r_diag, my_threshold, true);

			const std::vector<double> x = std::vector<double>(n];
			const std::vector<double> y = b.to_array();

			// apply Householder transforms to solve Q.y = b
			for (const int& minor = 0; minor < std::min(m, n); minor++)
			{
				const std::vector<double> qrt_minor = my_qrt[minor];
				double dot_product = 0;
				for (int row = minor; row < m; row++)
				{
					dot_product += y[row] * qrt_minor[row];
				}
				dot_product /= my_r_diag[minor] * qrt_minor[minor];

				for (int row = minor; row < m; row++)
				{
					y[row] += dot_product * qrt_minor[row];
				}
			}

			// solve triangular system R.x = y
			for (int row = my_r_diag.size() - 1; row >= 0; --row)
			{
				y[row] /= my_r_diag[row];
				const double y_row = y[row];
				const std::vector<double> qrt_row = my_qrt[row];
				x[row] = y_row;
				for (int i{}; i < row; i++)
				{
					y[i] -= y_row * qrt_row[i];
				}
			}

			return Array_Real_Vector(x, false);
		}

		/** {@inherit_doc} */
		//override
		public Real_Matrix solve(Real_Matrix b)
		{
			const int n = my_qrt.size();
			const int m = my_qrt[0].size();
			if (b.get_row_dimension() != m)
			{
				throw std::exception("not implemented");
				//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, b.get_row_dimension(), m);
			}
			check_singular(my_r_diag, my_threshold, true);

			const int columns = b.get_column_dimension();
			const int block_size = Block_Real_Matrix.BLOCK_SIZE;
			const int c_blocks = (columns + block_size - 1) / block_size;
			const autox_blocks = Block_Real_Matrix.create_blocks_layout(n, columns);
			const autoy = std::vector<double>(b.get_row_dimension()][block_size];
			const std::vector<double>   alpha = std::vector<double>(block_size];

			for (const int& k_block = 0; k_block < c_blocks; ++k_block)
			{
				const int& k_start = k_block * block_size;
				const int& k_end = std::min(k_start + block_size, columns);
				const int& k_width = k_end - k_start;

				// get the right hand side vector
				b.copy_sub_matrix(0, m - 1, k_start, k_end - 1, y);

				// apply Householder transforms to solve Q.y = b
				for (const int& minor = 0; minor < std::min(m, n); minor++)
				{
					const std::vector<double> qrt_minor = my_qrt[minor];
					const double factor = 1.0 / (my_r_diag[minor] * qrt_minor[minor]);

					Arrays.fill(alpha, 0, k_width, 0.0);
					for (int row = minor; row < m; ++row)
					{
						const double   d = qrt_minor[row];
						const std::vector<double> y_row = y[row];
						for (int k{}; k < k_width; ++k)
						{
							alpha[k] += d * y_row[k];
						}
					}
					for (int k{}; k < k_width; ++k)
					{
						alpha[k] *= factor;
					}

					for (int row = minor; row < m; ++row)
					{
						const double   d = qrt_minor[row];
						const std::vector<double> y_row = y[row];
						for (int k{}; k < k_width; ++k)
						{
							y_row[k] += alpha[k] * d;
						}
					}
				}

				// solve triangular system R.x = y
				for (int j = my_r_diag.size() - 1; j >= 0; --j)
				{
					const int      j_block = j / block_size;
					const int      j_start = j_block * block_size;
					const double   factor = 1.0 / my_r_diag[j];
					const std::vector<double> yJ = y[j];
					const std::vector<double> x_block = x_blocks[j_block * c_blocks + k_block];
					int index = (j - j_start) * k_width;
					for (int k{}; k < k_width; ++k)
					{
						yJ[k] *= factor;
						x_block[index++] = yJ[k];
					}

					const std::vector<double> qrtJ = my_qrt[j];
					for (int i{}; i < j; ++i)
					{
						const double rIJ = qrtJ[i];
						const std::vector<double> y_i = y[i];
						for (int k{}; k < k_width; ++k)
						{
							y_i[k] -= yJ[k] * rIJ;
						}
					}
				}
			}

			return Block_Real_Matrix(n, columns, x_blocks, false);
		}

		/**
		 * {@inherit_doc}
		 * @ if the decomposed matrix is singular.
		 */
		 //override
		public Real_Matrix get_inverse()
		{
			return solve(Matrix_Utils::create_real_identity_matrix(my_qrt[0].size()));
		}

		/**
		 * Check singularity.
		 *
		 * @param diag Diagonal elements of the R matrix.
		 * @param min Singularity my_threshold.
		 * @param raise Whether to raise a {@link }
		 * if any element of the diagonal fails the check.
		 * @return {@code true} if any element of the diagonal is smaller
		 * or equal to {@code min}.
		 * @ if the matrix is singular and
		 * {@code raise} is {@code true}.
		 */
		private bool check_singular(std::vector<double> diag, const double& min, bool raise)
		{
			const int len = diag.size();
			for (int i{}; i < len; i++)
			{
				const double d = diag[i];
				if (std::abs(d) <= min)
				{
					if (raise)
					{
						throw std::exception("not implemented");
						//throw (hipparchus::exception::Localized_Core_Formats_Type::SINGULAR_MATRIX);
					}
					return true;
				}
			}
			return false;
		}

		/** {@inherit_doc} */
		//override
		public int get_row_dimension()
		{
			return my_qrt[0].size();
		}

		/** {@inherit_doc} */
		//override
		public int get_column_dimension()
		{
			return my_qrt.size();
		}
	}
};