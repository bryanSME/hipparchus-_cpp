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
#include "RealVector.h"
#include "RealMatrix.h"

  /**
   * Interface handling decomposition algorithms that can solve A &times; X = B.
   * <p>
   * Decomposition algorithms decompose an A matrix as a product of several specific
   * matrices from which they can solve A &times; X = B in least squares sense: they find X
   * such that ||A &times; X - B|| is minimal.
   * <p>
   * Some solvers like {@link LU_Decomposition} can only find the solution for
   * square matrices and when the solution is an exact linear solution, i.e. when
   * ||A &times; X - B|| is exactly 0. Other solvers can also find solutions
   * with non-square matrix A and with non-null minimal norm. If an exact linear
   * solution exists it is also the minimal norm solution.
   *
   */
class Decomposition_Solver
{
	/**
	 * Solve the linear equation A &times; X = B for matrices A.
	 * <p>
	 * The A matrix is implicit, it is provided by the underlying
	 * decomposition algorithm.
	 *
	 * @param b right-hand side of the equation A &times; X = B
	 * @return a vector X that minimizes the two norm of A &times; X - B
	 * @org.hipparchus.exception.
	 * if the matrices dimensions do not match.
	 * @ if the decomposed matrix is singular.
	 */
	virtual Real_Vector solve(Real_Vector b) = 0;

	/**
	 * Solve the linear equation A &times; X = B for matrices A.
	 * <p>
	 * The A matrix is implicit, it is provided by the underlying
	 * decomposition algorithm.
	 *
	 * @param b right-hand side of the equation A &times; X = B
	 * @return a matrix X that minimizes the two norm of A &times; X - B
	 * @org.hipparchus.exception.
	 * if the matrices dimensions do not match.
	 * @ if the decomposed matrix is singular.
	 */
	virtual Real_Matrix solve(Real_Matrix b) = 0;

	/**
	 * Check if the decomposed matrix is non-singular.
	 * @return true if the decomposed matrix is non-singular.
	 */
	virtual bool is_non_singular() = 0;

	/**
	 * Get the <a href="http://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_pseudoinverse">pseudo-inverse</a>
	 * of the decomposed matrix.
	 * <p>
	 * <em>This is equal to the inverse  of the decomposed matrix, if such an inverse exists.</em>
	 * <p>
	 * If no such inverse exists, then the result has properties that resemble that of an inverse.
	 * <p>
	 * In particular, in this case, if the decomposed matrix is A, then the system of equations
	 * \( A x = b \) may have no solutions, or many. If it has no solutions, then the pseudo-inverse
	 * \( A^+ \) gives the "closest" solution \( z = A^+ b \), meaning \( \left \| A z - b \right \|_2 \)
	 * is minimized. If there are many solutions, then \( z = A^+ b \) is the smallest solution, * meaning \( \left \| z \right \|_2 \) is minimized.
	 * <p>
	 * Note however that some decompositions cannot compute a pseudo-inverse for all matrices.
	 * For example, the {@link LU_Decomposition} is not defined for non-square matrices to begin
	 * with. The {@link QR_Decomposition} can operate on non-square matrices, but will throw
	 * {@link } if the decomposed matrix is singular. Refer to the javadoc
	 * of specific decomposition implementations for more details.
	 *
	 * @return pseudo-inverse matrix (which is the inverse, if it exists), * if the decomposition can pseudo-invert the decomposed matrix
	 * @ if the decomposed matrix is singular and the decomposition
	 * can not compute a pseudo-inverse
	 */
	virtual Real_Matrix get_inverse() = 0;

	/**
	 * Returns the number of rows in the matrix.
	 *
	 * @return row_dimension
	 * @since 2.0
	 */
	virtual int get_row_dimension() = 0;

	/**
	 * Returns the number of columns in the matrix.
	 *
	 * @return column_dimension
	 * @since 2.0
	 */
	virtual int get_column_dimension() = 0;
};