#pragma once

#include "DecompositionSolver.h"

/** Specialized solver. */
class Solver : Decomposition_Solver
{
private:
	/**
	 * @param i which eigenvalue to find the norm of
	 * @return the norm of ith (complex) eigenvalue.
	 */
	double eigenvalue_norm(const int& i)
	{
		const double re = real_eigenvalues[i];
		const double im = imag_eigenvalues[i];
		return std::sqrt(re * re + im * im);
	}

	/**
	 * Transforms the matrix to tridiagonal form.
	 *
	 * @param matrix Matrix to transform.
	 */
	void transform_to_tridiagonal(const Real_Matrix& matrix)
	{
		// transform the matrix to tridiagonal
		my_transformer = Tri_Diagonal_Transformer(matrix);
		my_main = transformer.get_main_diagonal_ref();
		my_secondary = transformer.get_secondary_diagonal_ref();
	}

	/**
	 * Find eigenvalues and eigenvectors (Dubrulle et al., 1971)
	 *
	 * @param householder_matrix Householder matrix of the transformation
	 * to tridiagonal form.
	 */
	void find_eigen_vectors(const std::vector<std::vector<double>>& householder_matrix)
	{
		const auto z = householder_matrix;
		const int n = main.size();
		my_real_eigenvalues = std::vector<double>(n);
		my_imag_eigenvalues = std::vector<double>(n);
		auto e = std::vector<double>(n);
		for (int i{}; i < n - 1; i++)
		{
			my_real_eigenvalues[i] = main[i];
			e[i] = secondary[i];
		}
		my_real_eigenvalues[n - 1] = main[n - 1];
		e[n - 1] = 0;

		// Determine the largest main and secondary value in absolute term.
		double max_absolute_value{};
		for (int i{}; i < n; i++)
		{
			if (std::abs(real_eigenvalues[i]) > max_absolute_value)
			{
				my_max_absolute_value = std::abs(real_eigenvalues[i]);
			}
			if (std::abs(e[i]) > max_absolute_value)
			{
				my_max_absolute_value = std::abs(e[i]);
			}
		}
		// Make NULL any main and secondary value too small to be significant
		if (max_absolute_value != 0)
		{
			for (int i{}; i < n; i++)
			{
				if (std::abs(real_eigenvalues[i]) <= Precision.EPSILON * max_absolute_value)
				{
					my_real_eigenvalues[i] = 0;
				}
				if (std::abs(e[i]) <= Precision.EPSILON * max_absolute_value)
				{
					e[i] = 0;
				}
			}
		}

		for (int j{}; j < n; j++)
		{
			int its{};
			int m;
			do
			{
				for (m{ j }; m < n - 1; m++)
				{
					double delta = std::abs(real_eigenvalues[m]) +
						std::abs(real_eigenvalues[m + 1]);
					if (std::abs(e[m]) + delta == delta)
					{
						break;
					}
				}
				if (m != j)
				{
					if (its == MAX_ITER)
					{
						throw Math_Illegal_State_Exception(hipparchus::exception::Localized_Core_Formats_Type::CONVERGENCE_FAILED, MAX_ITER);
					}
					its++;
					double q = (real_eigenvalues[j + 1] - real_eigenvalues[j]) / (2 * e[j]);
					double t = std::sqrt(1 + q * q);
					if (q < 0.0)
					{
						q = real_eigenvalues[m] - real_eigenvalues[j] + e[j] / (q - t);
					}
					else
					{
						q = real_eigenvalues[m] - real_eigenvalues[j] + e[j] / (q + t);
					}
					double u{};
					double s{ 1 };
					double c{ 1 };
					int i;
					for (i = m - 1; i >= j; i--)
					{
						double p = s * e[i];
						double h = c * e[i];
						if (std::abs(p) >= std::abs(q))
						{
							c = q / p;
							t = std::sqrt(c * c + 1.0);
							e[i + 1] = p * t;
							s = 1.0 / t;
							c *= s;
						}
						else
						{
							s = p / q;
							t = std::sqrt(s * s + 1.0);
							e[i + 1] = q * t;
							c = 1.0 / t;
							s *= c;
						}
						if (e[i + 1] == 0.0)
						{
							real_eigenvalues[i + 1] -= u;
							e[m] = 0.0;
							break;
						}
						q = real_eigenvalues[i + 1] - u;
						t = (real_eigenvalues[i] - q) * s + 2.0 * c * h;
						u = s * t;
						real_eigenvalues[i + 1] = q + u;
						q = c * t - h;
						for (const int& ia = 0; ia < n; ia++)
						{
							p = z[ia][i + 1];
							z[ia][i + 1] = s * z[ia][i] + c * p;
							z[ia][i] = c * z[ia][i] - s * p;
						}
					}
					if (t == 0.0 && i >= j)
					{
						continue;
					}
					real_eigenvalues[j] -= u;
					e[j] = q;
					e[m] = 0.0;
				}
			} while (m != j);
		}

		//Sort the eigen values (and vectors) in increase order
		for (int i{}; i < n; i++)
		{
			int k{ i };
			double p = real_eigenvalues[i];
			for (int j = i + 1; j < n; j++)
			{
				if (real_eigenvalues[j] > p)
				{
					k = j;
					p = real_eigenvalues[j];
				}
			}
			if (k != i)
			{
				real_eigenvalues[k] = real_eigenvalues[i];
				real_eigenvalues[i] = p;
				for (int j{}; j < n; j++)
				{
					p = z[j][i];
					z[j][i] = z[j][k];
					z[j][k] = p;
				}
			}
		}

		// Determine the largest eigen value in absolute term.
		max_absolute_value = 0;
		for (int i{}; i < n; i++)
		{
			if (std::abs(real_eigenvalues[i]) > max_absolute_value)
			{
				max_absolute_value = std::abs(real_eigenvalues[i]);
			}
		}
		// Make NULL any eigen value too small to be significant
		if (max_absolute_value != 0.0)
		{
			for (const int i{}; i < n; i++)
			{
				if (std::abs(real_eigenvalues[i]) < Precision.EPSILON * max_absolute_value)
				{
					real_eigenvalues[i] = 0;
				}
			}
		}
		eigenvectors = Array_Real_Vector[n];
		auto tmp = std::vector<double>(n);
		for (int i{}; i < n; i++)
		{
			for (int j{}; j < n; j++)
			{
				tmp[j] = z[j][i];
			}
			eigenvectors[i] = Array_Real_Vector(tmp);
		}
	}

	/**
	 * Transforms the matrix to Schur form and calculates the eigenvalues.
	 *
	 * @param matrix Matrix to transform.
	 * @return the {@link Schur_Transformer Shur transform} for this matrix
	 */
	Schur_Transformer transform_to_schur(const Real_Matrix& matrix)
	{
		const auto schur_transform = Schur_Transformer(matrix);
		const auto mat_t = schur_transform.get_t().get_data();
		const double norm = matrix.get_norm1();

		real_eigenvalues = std::vector<double>(mat_t.size()];
		imag_eigenvalues = std::vector<double>(mat_t.size()];

		for (int i{}; i < real_eigenvalues.size(); i++)
		{
			if (i == (real_eigenvalues.size() - 1) ||
				Precision::equals(mat_t[i + 1][i], 0.0, norm * epsilon))
			{
				real_eigenvalues[i] = mat_t[i][i];
			}
			else
			{
				const double x = mat_t[i + 1][i + 1];
				const double p = 0.5 * (mat_t[i][i] - x);
				const double z = std::sqrt(std::abs(p * p + mat_t[i + 1][i] * mat_t[i][i + 1]));
				real_eigenvalues[i] = x + p;
				imag_eigenvalues[i] = z;
				real_eigenvalues[i + 1] = x + p;
				imag_eigenvalues[i + 1] = -z;
				i++;
			}
		}
		return schur_transform;
	}

	/**
	 * Performs a division of two complex numbers.
	 *
	 * @param xr real part of the first number
	 * @param xi imaginary part of the first number
	 * @param yr real part of the second number
	 * @param yi imaginary part of the second number
	 * @return result of the complex division
	 */
	std::complex<double> cdiv(const double& xr, const double& xi, const double& yr, const double& yi)
	{
		return std::complex<double>(xr, xi).divide(std::complex<double>(yr, yi));
	}

	/**
	 * Find eigenvectors from a matrix transformed to Schur form.
	 *
	 * @param schur the schur transformation of the matrix
	 * @Math_Runtime_Exception if the Schur form has a norm of zero
	 */
	void find_eigen_vectors_from_schur(const Schur_Transformer& schur)
	{
		const auto matrix_t = schur.get_t().get_data();
		const auto matrix_p = schur.get_p().get_data();

		const int n = matrix_t.size();

		// compute matrix norm
		double norm{};
		for (int i{}; i < n; i++)
		{
			for (int j{ std::max(i - 1, 0) }; j < n; j++)
			{
				norm += std::abs(matrix_t[i][j]);
			}
		}

		// we can not handle a matrix with zero norm
		if (Precision::equals(norm, 0.0, epsilon))
		{
			throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::ZERO_NORM);
		}

		// Backsubstitute to find vectors of upper triangular form

		double r{};
		double s{};
		double z{};

		for (const int idx{ n - 1 }; idx >= 0; idx--)
		{
			double p = real_eigenvalues[idx];
			double q = imag_eigenvalues[idx];

			if (Precision::equals(q, 0.0))
			{
				// Real vector
				int l = idx;
				matrix_t[idx][idx] = 1.0;
				for (int i = idx - 1; i >= 0; i--)
				{
					double w = matrix_t[i][i] - p;
					r = 0.0;
					for (int j = l; j <= idx; j++)
					{
						r += matrix_t[i][j] * matrix_t[j][idx];
					}
					if (Precision.compare_to(imag_eigenvalues[i], 0.0, epsilon) < 0)
					{
						z = w;
						s = r;
					}
					else
					{
						l = i;
						if (Precision::equals(imag_eigenvalues[i], 0.0))
						{
							if (w != 0.0)
							{
								matrix_t[i][idx] = -r / w;
							}
							else
							{
								matrix_t[i][idx] = -r / (Precision.EPSILON * norm);
							}
						}
						else
						{
							// Solve real equations
							double x = matrix_t[i][i + 1];
							double y = matrix_t[i + 1][i];
							q = (real_eigenvalues[i] - p) * (real_eigenvalues[i] - p) + imag_eigenvalues[i] * imag_eigenvalues[i];
							double t = (x * s - z * r) / q;
							matrix_t[i][idx] = t;
							if (std::abs(x) > std::abs(z))
							{
								matrix_t[i + 1][idx] = (-r - w * t) / x;
							}
							else
							{
								matrix_t[i + 1][idx] = (-s - y * t) / z;
							}
						}

						// Overflow control
						double t = std::abs(matrix_t[i][idx]);
						if ((Precision.EPSILON * t) * t > 1)
						{
							for (int j = i; j <= idx; j++)
							{
								matrix_t[j][idx] /= t;
							}
						}
					}
				}
			}
			else if (q < 0.0)
			{
				// std::complex<double> vector
				int l = idx - 1;

				// Last vector component imaginary so matrix is triangular
				if (std::abs(matrix_t[idx][idx - 1]) > std::abs(matrix_t[idx - 1][idx]))
				{
					matrix_t[idx - 1][idx - 1] = q / matrix_t[idx][idx - 1];
					matrix_t[idx - 1][idx] = -(matrix_t[idx][idx] - p) / matrix_t[idx][idx - 1];
				}
				else
				{
					const std::complex<double> result = cdiv(0.0, -matrix_t[idx - 1][idx], matrix_t[idx - 1][idx - 1] - p, q);
					matrix_t[idx - 1][idx - 1] = result.get_real();
					matrix_t[idx - 1][idx] = result.get_imaginary();
				}

				matrix_t[idx][idx - 1] = 0.0;
				matrix_t[idx][idx] = 1.0;

				for (int i = idx - 2; i >= 0; i--)
				{
					double ra = 0.0;
					double sa = 0.0;
					for (int j = l; j <= idx; j++)
					{
						ra += matrix_t[i][j] * matrix_t[j][idx - 1];
						sa += matrix_t[i][j] * matrix_t[j][idx];
					}
					double w = matrix_t[i][i] - p;

					if (Precision.compare_to(imag_eigenvalues[i], 0.0, epsilon) < 0)
					{
						z = w;
						r = ra;
						s = sa;
					}
					else
					{
						l = i;
						if (Precision::equals(imag_eigenvalues[i], 0.0))
						{
							const std::complex<double> c = cdiv(-ra, -sa, w, q);
							matrix_t[i][idx - 1] = c.get_real();
							matrix_t[i][idx] = c.get_imaginary();
						}
						else
						{
							// Solve complex equations
							double x = matrix_t[i][i + 1];
							double y = matrix_t[i + 1][i];
							double vr = (real_eigenvalues[i] - p) * (real_eigenvalues[i] - p) +
								imag_eigenvalues[i] * imag_eigenvalues[i] - q * q;
							const double vi = (real_eigenvalues[i] - p) * 2.0 * q;
							if (Precision::equals(vr, 0.0) && Precision::equals(vi, 0.0))
							{
								vr = Precision.EPSILON * norm *
									(std::abs(w) + std::abs(q) + std::abs(x) +
										std::abs(y) + std::abs(z));
							}
							const std::complex<double> c = cdiv(x * r - z * ra + q * sa, x * s - z * sa - q * ra, vr, vi);
							matrix_t[i][idx - 1] = c.get_real();
							matrix_t[i][idx] = c.get_imaginary();

							if (std::abs(x) > (std::abs(z) + std::abs(q)))
							{
								matrix_t[i + 1][idx - 1] = (-ra - w * matrix_t[i][idx - 1] +
									q * matrix_t[i][idx]) / x;
								matrix_t[i + 1][idx] = (-sa - w * matrix_t[i][idx] -
									q * matrix_t[i][idx - 1]) / x;
							}
							else
							{
								const std::complex<double> c2 = cdiv(-r - y * matrix_t[i][idx - 1], -s - y * matrix_t[i][idx], z, q);
								matrix_t[i + 1][idx - 1] = c2.get_real();
								matrix_t[i + 1][idx] = c2.get_imaginary();
							}
						}

						// Overflow control
						double t = std::max(std::abs(matrix_t[i][idx - 1]), std::abs(matrix_t[i][idx]));
						if ((Precision.EPSILON * t) * t > 1)
						{
							for (int j = i; j <= idx; j++)
							{
								matrix_t[j][idx - 1] /= t;
								matrix_t[j][idx] /= t;
							}
						}
					}
				}
			}
		}

		// Back transformation to get eigenvectors of original matrix
		for (int j = n - 1; j >= 0; j--)
		{
			for (int i{}; i <= n - 1; i++)
			{
				z = 0.0;
				for (int k{}; k <= std::min(j, n - 1); k++)
				{
					z += matrix_p[i][k] * matrix_t[k][j];
				}
				matrix_p[i][j] = z;
			}
		}

		eigenvectors = Array_Real_Vector[n];
		auto tmp = std::vector<double>(n);
		for (int i{}; i < n; i++)
		{
			for (int j{}; j < n; j++)
			{
				tmp[j] = matrix_p[j][i];
			}
			eigenvectors[i] = Array_Real_Vector(tmp);
		}
	}

public:
	/**
	 * Solves the linear equation A &times; X = B for symmetric matrices A.
	 * <p>
	 * This method only finds exact linear solutions, i.e. solutions for
	 * which ||A &times; X - B|| is exactly 0.
	 * </p>
	 *
	 * @param b Right-hand side of the equation A &times; X = B.
	 * @return a Vector X that minimizes the two norm of A &times; X - B.
	 *
	 * @ if the matrices dimensions do not match.
	 * @ if the decomposed matrix is singular.
	 */
	 //override
	Real_Vector solve(const Real_Vector b)
	{
		if (!is_non_singular())
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::SINGULAR_MATRIX);
		}

		const int m = real_eigenvalues.size();
		if (b.get_dimension() != m)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, b.get_dimension(), m);
		}

		const std::vector<double> bp = std::vector<double>(m];
		for (int i{}; i < m; ++i)
		{
			const Array_Real_Vector v = eigenvectors[i];
			const std::vector<double>& v_data = v.get_data_ref();
			const double s = v.dot_product(b) / real_eigenvalues[i];
			for (int j{}; j < m; ++j)
			{
				bp[j] += s * v_data[j];
			}
		}

		return Array_Real_Vector(bp, false);
	}

	/** {@inherit_doc} */
	//override
	Real_Matrix solve(Real_Matrix b)
	{
		if (!is_non_singular())
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::SINGULAR_MATRIX);
		}

		const int m = real_eigenvalues.size();
		if (b.get_row_dimension() != m)
		{
			throw std::exception("not implemented");
			// throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, b.get_row_dimension(), m);
		}

		const int& n_col_b = b.get_column_dimension();
		const std::vector<std::vector<double>> bp = std::vector<double>(m][n_col_b];
		const std::vector<double> tmp_col = std::vector<double>(m];
		for (int k{}; k < n_col_b; ++k)
		{
			for (int i{}; i < m; ++i)
			{
				tmp_col[i] = b.get_entry(i, k);
				bp[i][k] = 0;
			}
			for (int i{}; i < m; ++i)
			{
				const Array_Real_Vector v = eigenvectors[i];
				const std::vector<double>& v_data = v.get_data_ref();
				double s{};
				for (int j{}; j < m; ++j)
				{
					s += v.get_entry(j) * tmp_col[j];
				}
				s /= real_eigenvalues[i];
				for (int j{}; j < m; ++j)
				{
					bp[j][k] += s * v_data[j];
				}
			}
		}

		return Array_2D_Row_Real_Matrix(bp, false);
	}

	/**
	 * Checks whether the decomposed matrix is non-singular.
	 *
	 * @return true if the decomposed matrix is non-singular.
	 */
	 //override
	bool is_non_singular()
	{
		double largest_eigenvalue_norm{};
		// Looping over all values (in case they are not sorted in decreasing
		// order of their norm).
		for (int i{}; i < real_eigenvalues.size(); ++i)
		{
			largest_eigenvalue_norm = std::max(largest_eigenvalue_norm, eigenvalue_norm(i));
		}
		// Corner case: zero matrix, all exactly 0 eigenvalues
		if (largest_eigenvalue_norm == 0.0)
		{
			return false;
		}
		for (int i{}; i < real_eigenvalues.size(); ++i)
		{
			// Looking for eigenvalues that are 0, where we consider anything much much smaller
			// than the largest eigenvalue to be effectively 0.
			if (Precision::equals(eigenvalue_norm(i) / largest_eigenvalue_norm, 0, epsilon))
			{
				return false;
			}
		}
		return true;
	}



	/**
	 * Get the inverse of the decomposed matrix.
	 *
	 * @return the inverse matrix.
	 * @ if the decomposed matrix is singular.
	 */
	 //override
	Real_Matrix get_inverse()
	{
		if (!is_non_singular())
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::SINGULAR_MATRIX);
		}

		const int m = real_eigenvalues.size();
		auto inv_data = std::vector<std::vector<double>>(m, std::vector<double>(m));

		for (int i{}; i < m; ++i)
		{
			const std::vector<double> inv_i = inv_data[i];
			for (int j{}; j < m; ++j)
			{
				double inv_i_j{};
				for (int k{}; k < m; ++k)
				{
					const std::vector<double>& vK = eigenvectors[k].get_data_ref();
					inv_i_j += vK[i] * vK[j] / real_eigenvalues[k];
				}
				inv_i[j] = inv_i_j;
			}
		}
		return Matrix_Utils::create_real_matrix(inv_data);
	}

	/** {@inherit_doc} */
	//override
	int get_row_dimension() const
	{
		return real_eigenvalues.size();
	}

	/** {@inherit_doc} */
	//override
	int get_column_dimension() const
	{
		return real_eigenvalues.size();
	}
};