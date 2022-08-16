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
#include <vector>
#include "MatrixUtils.h"
#include "QRDecomposer.h"
#include "RealMatrix.h"
//import org.hipparchus.util.FastMath;

/**
 * Calculates the rank-revealing QR-decomposition of a matrix, with column pivoting.
 * <p>The rank-revealing QR-decomposition of a matrix A consists of three matrices Q, * R and P such that AP=QR.  Q is orthogonal (Q<sup>T</sup>Q = I), and R is upper triangular.
 * If A is m&times;n, Q is m&times;m and R is m&times;n and P is n&times;n.</p>
 * <p>QR decomposition with column pivoting produces a rank-revealing QR
 * decomposition and the {@link #get_rankstatic_cast<double>(} method may be used to return the rank of the
 * input matrix A.</p>
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
class RRQR_Decomposition : public QR_Decomposition
{
private:
	/** An array to record the column pivoting for later creation of P. */
	std::vector<int> p;

	/** Cached value of P. */
	Real_Matrix cached_p;

public:
	/**
	 * Calculates the QR-decomposition of the given matrix.
	 * The singularity threshold defaults to zero.
	 *
	 * @param matrix The matrix to decompose.
	 *
	 * @see #RRQR_Decomposition(Real_Matrix, double)
	 */
	RRQR_Decomposition(Real_Matrix matrix)
	{
		this(matrix, 0.0);
	}

	/**
	  * Calculates the QR-decomposition of the given matrix.
	  *
	  * @param matrix The matrix to decompose.
	  * @param threshold Singularity threshold.
	  * @see #RRQR_Decomposition(Real_Matrix)
	  */
	RRQR_Decomposition(Real_Matrix matrix, double threshold)
	{
		super(matrix, threshold);
	}

	/**
	 * Returns the pivot matrix, P, used in the QR Decomposition of matrix A such that AP = QR.
	 *
	 * If no pivoting is used in this decomposition then P is equal to the identity matrix.
	 *
	 * @return a permutation matrix.
	 */
	Real_Matrix get_p()
	{
		if (cached_p == NULL)
		{
			int n = p.size();
			cached_p = Matrix_Utils::create_real_matrix(n, n);
			for (int i{}; i < n; i++)
			{
				cached_p.set_entry(p[i], i, 1);
			}
		}
		return cached_p;
	}

	/**
	 * Return the effective numerical matrix rank.
	 * <p>The effective numerical rank is the number of non-negligible
	 * singular values.</p>
	 * <p>This implementation looks at Frobenius norms of the sequence of
	 * bottom right submatrices.  When a large fall in norm is seen, * the rank is returned. The drop is computed as:</p>
	 * <pre>
	 *   (this_norm/last_norm) * r_norm &lt; drop_threshold
	 * </pre>
	 * <p>
	 * where this_norm is the Frobenius norm of the current submatrix, * last_norm is the Frobenius norm of the previous submatrix, * r_norm is is the Frobenius norm of the complete matrix
	 * </p>
	 *
	 * @param drop_threshold threshold triggering rank computation
	 * @return effective numerical matrix rank
	 */
	int get_rank(const double drop_threshold)
	{
		Real_Matrix r = get_r();
		int rows = r.get_row_dimension();
		int columns = r.get_column_dimension();
		int rank = 1;
		double last_norm = r.get_frobenius_norm();
		double r_norm = last_norm;
		while (rank < std::min(rows, columns))
		{
			double this_norm = r.get_sub_matrix(rank, rows - 1, rank, columns - 1).get_frobenius_norm();
			if (this_norm == 0 || (this_norm / last_norm) * r_norm < drop_threshold)
			{
				break;
			}
			last_norm = this_norm;
			rank++;
		}
		return rank;
	}

	/**
	 * Get a solver for finding the A &times; X = B solution in least square sense.
	 * <p>
	 * Least Square sense means a solver can be computed for an overdetermined system, * (i.e. a system with more equations than unknowns, which corresponds to a tall A
	 * matrix with more rows than columns). In any case, if the matrix is singular
	 * within the tolerance set at {@link RRQR_Decomposition#RRQR_Decomposition(Real_Matrix, * double) construction}, an error will be triggered when
	 * the {@link Decomposition_Solver#solve(Real_Vector) solve} method will be called.
	 * </p>
	 * @return a solver
	 */
	 //override
	Decomposition_Solver get_solver()
	{
		return Solver(super.get_solver(), this.get_p());
	}

	/** Specialized solver. */
	static class Solver : public Decomposition_Solver
	{
		/** Upper level solver. */
		private const Decomposition_Solver upper;

		/** A permutation matrix for the pivots used in the QR decomposition */
		private Real_Matrix p;

		/**
		 * Build a solver from decomposed matrix.
		 *
		 * @param upper upper level solver.
		 * @param p permutation matrix
		 */
		private Solver(const Decomposition_Solver upper, const Real_Matrix p)
		{
			this.upper = upper;
			this.p = p;
		}

		/** {@inherit_doc} */
		//override
		public bool is_non_singular()
		{
			return upper.is_non_singular();
		}

		/** {@inherit_doc} */
		//override
		public Real_Vector solve(Real_Vector b)
		{
			return p.operate(upper.solve(b));
		}

		/** {@inherit_doc} */
		//override
		public Real_Matrix solve(Real_Matrix b)
		{
			return p.multiply(upper.solve(b));
		}

		/**
		 * {@inherit_doc}
		 * @org.hipparchus.exception.
		 * if the decomposed matrix is singular.
		 */
		 //override
		public Real_Matrix get_inverse()
		{
			return solve(Matrix_Utils::create_real_identity_matrix(p.get_row_dimension()));
		}
		/** {@inherit_doc} */
		//override
		public int get_row_dimension()
		{
			return upper.get_row_dimension();
		}

		/** {@inherit_doc} */
		//override
		public int get_column_dimension()
		{
			return upper.get_column_dimension();
		}
	}
}
