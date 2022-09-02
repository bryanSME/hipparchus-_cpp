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

  //import org.hipparchus.complex.std::complex<double>;
  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.exception.Math_Runtime_Exception;
  //import org.hipparchus.exception.Math_Illegal_State_Exception;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Precision;
#include <vector>
#include "MatrixUtils.h"
#include "Solver.h"

/**
 * Calculates the eigen decomposition of a real matrix.
 * <p>
 * The eigen decomposition of matrix A is a set of two matrices:
 * V and D such that A = V &times; D &times; V<sup>T</sup>.
 * A, V and D are all m &times; m matrices.
 * <p>
 * This class is similar in spirit to the {@code Eigenvalue_Decomposition}
 * class from the <a href="http://math.nist.gov/javanumerics/jama/">JAMA</a>
 * library, with the following changes:
 * <ul>
 *   <li>a {@link #get_v_t() get_vt} method has been added,</li>
 *   <li>two {@link #get_real_eigenvaluestatic_cast<int>( get_real_eigenvalue} and
 *       {@link #get_imag_eigenvaluestatic_cast<int>( get_imag_eigenvalue} methods to pick up a
 *       single eigenvalue have been added,</li>
 *   <li>a {@link #get_eigenvectorstatic_cast<int>( get_eigenvector} method to pick up a
 *       single eigenvector has been added,</li>
 *   <li>a {@link #get_determinant() get_determinant} method has been added.</li>
 *   <li>a {@link #get_solver() get_solver} method has been added.</li>
 * </ul>
 * <p>
 * As of 3.1, this class supports general real matrices (both symmetric and non-symmetric):
 * <p>
 * If A is symmetric, then A = V*D*V' where the eigenvalue matrix D is diagonal
 * and the eigenvector matrix V is orthogonal, i.e.
 * {@code A = V.multiply(D.multiply(V.transpose()))} and
 * {@code V.multiply(V.transpose())} equals the identity matrix.
 * </p>
 * <p>
 * If A is not symmetric, then the eigenvalue matrix D is block diagonal with the real
 * eigenvalues in 1-by-1 blocks and any complex eigenvalues, lambda + i*mu, in 2-by-2
 * blocks:
 * <pre>
 *    [lambda, mu    ]
 *    [   -mu, lambda]
 * </pre>
 * The columns of V represent the eigenvectors in the sense that {@code A*V = V*D}, * i.e. A.multiply(V) equals V.multiply(D).
 * The matrix V may be badly conditioned, or even singular, so the validity of the
 * equation {@code A = V*D*inverse(V)} depends upon the condition of V.
 * <p>
 * This implementation is based on the paper by A. Drubrulle, R.S. Martin and
 * J.H. Wilkinson "The Implicit QL Algorithm" in Wilksinson and Reinsch (1971)
 * Handbook for automatic computation, vol. 2, Linear algebra, Springer-Verlag, * New-York.
 *
 * @see <a href="http://mathworld.wolfram.com/Eigen_Decomposition.html">MathWorld</a>
 * @see <a href="http://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix">Wikipedia</a>
 */
class Eigen_Decomposition
{
private:
	/** Default epsilon value to use for internal epsilon **/
	static constexpr double DEFAULT_EPSILON{ 1e-12 };
	/** Maximum number of iterations accepted in the implicit QL transformation */
	static const std::byte MAX_ITER{ 30 };
	/** Internally used epsilon criteria. */
	const double my_epsilon;
	/** Main diagonal of the tridiagonal matrix. */
	std::vector<double> my_main;
	/** Secondary diagonal of the tridiagonal matrix. */
	std::vector<double> my_secondary;
	/**
	 * Transformer to tridiagonal (may be NULL if matrix is already
	 * tridiagonal).
	 */
	Tri_Diagonal_Transformer my_transformer;
	/** Real part of the real_eigenvalues. */
	std::vector<double> my_real_eigenvalues;
	/** Imaginary part of the real_eigenvalues. */
	std::vector<double> my_imag_eigenvalues;
	/** Eigenvectors. */
	std::vector<Array_Real_Vector> my_eigenvectors;
	/** Cached value of V. */
	Real_Matrix my_cached_v;
	/** Cached value of D. */
	Real_Matrix my_cached_d;
	/** Cached value of Vt. */
	Real_Matrix my_cached_vt;
	/** Whether the matrix is symmetric. */
	const bool my_is_symmetric;

public:
	/**
	 * Calculates the eigen decomposition of the given real matrix.
	 * <p>
	 * Supports decomposition of a general matrix since 3.1.
	 *
	 * @param matrix Matrix to decompose.
	 * @Math_Illegal_State_Exception if the algorithm fails to converge.
	 * @Math_Runtime_Exception if the decomposition of a general matrix
	 * results in a matrix with zero norm
	 */
	Eigen_Decomposition(const Real_Matrix& matrix)
	{
		Eigen_Decomposition(matrix, DEFAULT_EPSILON);
	}

	/**
	 * Calculates the eigen decomposition of the given real matrix.
	 * <p>
	 * Supports decomposition of a general matrix since 3.1.
	 *
	 * @param matrix Matrix to decompose.
	 * @param epsilon Epsilon used for internal tests (e.g. is singular, eigenvalue ratio, etc.)
	 * @Math_Illegal_State_Exception if the algorithm fails to converge.
	 * @Math_Runtime_Exception if the decomposition of a general matrix
	 * results in a matrix with zero norm
	 */
	Eigen_Decomposition(const Real_Matrix& matrix, const double& epsilon)
		:
		my_epsilon{ epsilon }
	{
		const double sym_tol = 10 * matrix.get_row_dimension() * matrix.get_column_dimension() * Precision.EPSILON;
		my_is_symmetric = Matrix_Utils::is_symmetric(matrix, sym_tol);
		if (is_symmetric)
		{
			transform_to_tridiagonal(matrix);
			find_eigen_vectors(transformer.get_q().get_data());
		}
		else
		{
			const Schur_Transformer t = transform_to_schur(matrix);
			find_eigen_vectors_from_schur(t);
		}
	}

	/**
	 * Calculates the eigen decomposition of the symmetric tridiagonal
	 * matrix.  The Householder matrix is assumed to be the identity matrix.
	 *
	 * @param main Main diagonal of the symmetric tridiagonal form.
	 * @param secondary Secondary of the tridiagonal form.
	 * @Math_Illegal_State_Exception if the algorithm fails to converge.
	 */
	Eigen_Decomposition(const std::vector<double>& main, const std::vector<double>& secondary)
	{
		Eigen_Decomposition(main, secondary, DEFAULT_EPSILON);
	}

	/**
	 * Calculates the eigen decomposition of the symmetric tridiagonal
	 * matrix.  The Householder matrix is assumed to be the identity matrix.
	 *
	 * @param main Main diagonal of the symmetric tridiagonal form.
	 * @param secondary Secondary of the tridiagonal form.
	 * @param epsilon Epsilon used for internal tests (e.g. is singular, eigenvalue ratio, etc.)
	 * @Math_Illegal_State_Exception if the algorithm fails to converge.
	 */
	Eigen_Decomposition(const std::vector<double>& main, const std::vector<double>& secondary, const double& epsilon)
		:
		my_epsilon{ epsilon },
		my_is_symmetric{ true },
		my_main{ main },
		my_secondary{ secondary },
		my_transformer{ NULL },
	{
		const int size = main.size();
		auto z = std::vector<std::vector<double>>(size, std::vector<double>(size));
		for (int i{}; i < size; i++)
		{
			z[i][i] = 1.0;
		}
		find_eigen_vectors(z);
	}

	/**
	 * Gets the matrix V of the decomposition.
	 * V is an orthogonal matrix, i.e. its transpose is also its inverse.
	 * The columns of V are the eigenvectors of the original matrix.
	 * No assumption is made about the orientation of the system axes formed
	 * by the columns of V (e.g. in a 3-dimension space, V can form a left-
	 * or right-handed system).
	 *
	 * @return the V matrix.
	 */
	Real_Matrix get_v()
	{
		if (cached_v == NULL)
		{
			const int m = eigenvectors.size();
			cached_v = Matrix_Utils::create_real_matrix(m, m);
			for (int k{}; k < m; ++k)
			{
				cached_v.set_column_vector(k, eigenvectors[k]);
			}
		}
		// return the cached matrix
		return cached_v;
	}

	/**
	 * Gets the block diagonal matrix D of the decomposition.
	 * D is a block diagonal matrix.
	 * Real eigenvalues are on the diagonal while complex values are on
	 * 2x2 blocks { {real +imaginary}, {-imaginary, real} }.
	 *
	 * @return the D matrix.
	 *
	 * @see #get_real_eigenvalues()
	 * @see #get_imag_eigenvalues()
	 */
	Real_Matrix get_d()
	{
		if (cached_d == NULL)
		{
			// cache the matrix for subsequent calls
			cached_d = Matrix_Utils::create_real_matrix(real_eigenvalues.size(), real_eigenvalues.size());
			for (int i{}; i < real_eigenvalues.size(); ++i)
			{
				cached_d.set_entry(i, i, real_eigenvalues[i]);
			}

			for (int i{}; i < imag_eigenvalues.size(); i++)
			{
				if (Precision.compare_to(imag_eigenvalues[i], 0.0, epsilon) > 0)
				{
					cached_d.set_entry(i, i + 1, imag_eigenvalues[i]);
				}
				else if (Precision.compare_to(imag_eigenvalues[i], 0.0, epsilon) < 0)
				{
					cached_d.set_entry(i, i - 1, imag_eigenvalues[i]);
				}
			}
		}
		return cached_d;
	}

	/**
	 * Get's the value for epsilon which is used for internal tests (e.g. is singular, eigenvalue ratio, etc.)
	 *
	 * @return the epsilon value.
	 */
	double get_epsilon() { return epsilon; }

	/**
	 * Gets the transpose of the matrix V of the decomposition.
	 * V is an orthogonal matrix, i.e. its transpose is also its inverse.
	 * The columns of V are the eigenvectors of the original matrix.
	 * No assumption is made about the orientation of the system axes formed
	 * by the columns of V (e.g. in a 3-dimension space, V can form a left-
	 * or right-handed system).
	 *
	 * @return the transpose of the V matrix.
	 */
	Real_Matrix get_v_t()
	{
		if (cached_vt == NULL)
		{
			const int m = eigenvectors.size();
			cached_vt = Matrix_Utils::create_real_matrix(m, m);
			for (int k{}; k < m; ++k)
			{
				cached_vt.set_row_vector(k, eigenvectors[k]);
			}
		}

		// return the cached matrix
		return cached_vt;
	}

	/**
	 * Returns whether the calculated eigen values are complex or real.
	 * <p>The method performs a zero check for each element of the
	 * {@link #get_imag_eigenvalues()} array and returns {@code true} if any
	 * element is not equal to zero.
	 *
	 * @return {@code true} if the eigen values are complex, {@code false} otherwise
	 */
	bool has_complex_eigenvalues()
	{
		for (int i{}; i < imag_eigenvalues.size(); i++)
		{
			if (!Precision::equals(imag_eigenvalues[i], 0.0, epsilon))
			{
				return true;
			}
		}
		return false;
	}

	/**
	 * Gets a copy of the real parts of the eigenvalues of the original matrix.
	 *
	 * @return a copy of the real parts of the eigenvalues of the original matrix.
	 *
	 * @see #get_d()
	 * @see #get_real_eigenvaluestatic_cast<int>(
	 * @see #get_imag_eigenvalues()
	 */
	std::vector<double> get_real_eigenvalues() const
	{
		return my_real_eigenvalues;
	}

	/**
	 * Returns the real part of the i<sup>th</sup> eigenvalue of the original
	 * matrix.
	 *
	 * @param i index of the eigenvalue (counting from 0)
	 * @return real part of the i<sup>th</sup> eigenvalue of the original
	 * matrix.
	 *
	 * @see #get_d()
	 * @see #get_real_eigenvalues()
	 * @see #get_imag_eigenvaluestatic_cast<int>(
	 */
	double get_real_eigenvalue(const int& i) const
	{
		return my_real_eigenvalues[i];
	}

	/**
	 * Gets a copy of the imaginary parts of the eigenvalues of the original
	 * matrix.
	 *
	 * @return a copy of the imaginary parts of the eigenvalues of the original
	 * matrix.
	 *
	 * @see #get_d()
	 * @see #get_imag_eigenvaluestatic_cast<int>(
	 * @see #get_real_eigenvalues()
	 */
	std::vector<double> get_imag_eigenvalues() const
	{
		return my_imag_eigenvalues;
	}

	/**
	 * Gets the imaginary part of the i<sup>th</sup> eigenvalue of the original
	 * matrix.
	 *
	 * @param i Index of the eigenvalue (counting from 0).
	 * @return the imaginary part of the i<sup>th</sup> eigenvalue of the original
	 * matrix.
	 *
	 * @see #get_d()
	 * @see #get_imag_eigenvalues()
	 * @see #get_real_eigenvaluestatic_cast<int>(
	 */
	double get_imag_eigenvalue(const int& i) const
	{
		return my_imag_eigenvalues[i];
	}

	/**
	 * Gets a copy of the i<sup>th</sup> eigenvector of the original matrix.
	 *
	 * @param i Index of the eigenvector (counting from 0).
	 * @return a copy of the i<sup>th</sup> eigenvector of the original matrix.
	 * @see #get_d()
	 */
	Real_Vector get_eigenvector(const int& i) const
	{
		return my_eigenvectors[i];
	}

	/**
	 * Computes the determinant of the matrix.
	 *
	 * @return the determinant of the matrix.
	 */
	double get_determinant()
	{
		double determinant{ 1 };
		for (const auto& lambda : my_real_eigenvalues)
		{
			determinant *= lambda;
		}
		return determinant;
	}

	/**
	 * Computes the square-root of the matrix.
	 * This implementation assumes that the matrix is symmetric and positive
	 * definite.
	 *
	 * @return the square-root of the matrix.
	 * @Math_Runtime_Exception if the matrix is not
	 * symmetric or not positive definite.
	 */
	Real_Matrix get_square_root()
	{
		if (!is_symmetric)
		{
			throw std::exception("not implemented");
			//throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::UNSUPPORTED_OPERATION);
		}

		auto sqrt_eigen_values = std::vector<double>(my_real_eigenvalues.size());
		for (int i{}; i < my_real_eigenvalues.size(); i++)
		{
			const double eigen{ my_real_eigenvalues[i] };
			if (eigen <= 0)
			{
				throw std::exception("not implemented");
				//throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::UNSUPPORTED_OPERATION);
			}
			sqrt_eigen_values[i] = std::sqrt(eigen);
		}
		const Real_Matrix sqrt_eigen = Matrix_Utils::create_real_diagonal_matrix(sqrt_eigen_values);
		const Real_Matrix v = get_v();
		const Real_Matrix vT = get_v_t();

		return v.multiply(sqrt_eigen).multiply(vT);
	}

	/**
	 * Gets a solver for finding the A &times; X = B solution in exact
	 * linear sense.
	 * <p>
	 * sin_ce 3.1, eigen decomposition of a general matrix is supported, * but the {@link Decomposition_Solver} only supports real eigenvalues.
	 *
	 * @return a solver
	 * @Math_Runtime_Exception if the decomposition resulted in
	 * complex eigenvalues
	 */
	Decomposition_Solver get_solver()
	{
		if (has_complex_eigenvalues())
		{
			throw std::exception("not implemented");
			//throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::UNSUPPORTED_OPERATION);
		}
		return Solver();
	}
};