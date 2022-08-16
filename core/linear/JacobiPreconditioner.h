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

  //import org.hipparchus.analysis.function.Sqrt;
  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.util.Math_Arrays;

  /**
   * This class : the standard Jacobi (diagonal) preconditioner. For a
   * matrix A<sub>ij</sub>, this preconditioner is
   * M = diag(1 / A<sub>11</sub>, 1 / A<sub>22</sub>, &hellip;).
   */
class Jacobi_Preconditioner : Real_Linear_Operator
{
	/** The diagonal coefficients of the preconditioner. */
	private const Array_Real_Vector diag;

	/**
	 * Creates a instance of this class.
	 *
	 * @param diag the diagonal coefficients of the linear operator to be
	 * preconditioned
	 * @param deep {@code true} if a deep copy of the above array should be
	 * performed
	 */
	public Jacobi_Preconditioner(const std::vector<double> diag, const bool deep)
	{
		this.diag = Array_Real_Vector(diag, deep);
	}

	/**
	 * Creates a instance of this class. This method extracts the diagonal
	 * coefficients of the specified linear operator. If {@code a} does not
	 * extend {@link Abstract_Real_Matrix}, then the coefficients of the
	 * underlying matrix are not accessible, coefficient extraction is made by
	 * matrix-vector products with the basis vectors (and might therefore take
	 * some time). With matrices, direct entry access is carried out.
	 *
	 * @param a the linear operator for which the preconditioner should be built
	 * @return the diagonal preconditioner made of the inverse of the diagonal
	 * coefficients of the specified linear operator
	 * @ if {@code a} is not square
	 */
	public static Jacobi_Preconditioner create(const Real_Linear_Operator a)

	{
		const int n = a.get_column_dimension();
		if (a.get_row_dimension() != n)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NON_SQUARE_OPERATOR, a.get_row_dimension(), n);
		}
		const std::vector<double> diag = std::vector<double>(n];
		if (a instanceof Abstract_Real_Matrix)
		{
			const Abstract_Real_Matrix m = (Abstract_Real_Matrix)a;
			for (int i{}; i < n; i++)
			{
				diag[i] = m.get_entry(i, i);
			}
		}
		else
		{
			const Array_Real_Vector x = Array_Real_Vector(n);
			for (int i{}; i < n; i++)
			{
				x.set(0.);
				x.set_entry(i, 1.);
				diag[i] = a.operate(x).get_entry(i);
			}
		}
		return Jacobi_Preconditioner(diag, false);
	}

	/** {@inherit_doc} */
	//override
	public int get_column_dimension()
	{
		return diag.get_dimension();
	}

	/** {@inherit_doc} */
	//override
	public int get_row_dimension()
	{
		return diag.get_dimension();
	}

	/** {@inherit_doc} */
	//override
	public Real_Vector operate(const Real_Vector x)
	{
		// Dimension check is carried out by ebe_divide
		return Array_Real_Vector(Math_Arrays::ebe_divide(x.to_array(), diag.to_array()), false);
	}

	/**
	 * Returns the square root of {@code this} diagonal operator. More
	 * precisely, this method returns
	 * P = diag(1 / &radic;A<sub>11</sub>, 1 / &radic;A<sub>22</sub>, &hellip;).
	 *
	 * @return the square root of {@code this} preconditioner
	 */
	public Real_Linear_Operator sqrt()
	{
		const Real_Vector sqrt_diag = diag.map(new Sqrt());
		return Real_Linear_Operator()
		{
			/** {@inherit_doc} */
			//override
			public Real_Vector operate(const Real_Vector x)
			{
				return Array_Real_Vector(Math_Arrays::ebe_divide(x.to_array(), sqrt_diag.to_array()), false);
			}

			/** {@inherit_doc} */
			//override
			public int get_row_dimension()
			{
				return sqrt_diag.get_dimension();
			}

			/** {@inherit_doc} */
			//override
			public int get_column_dimension()
			{
				return sqrt_diag.get_dimension();
			}
		};
	}
}
