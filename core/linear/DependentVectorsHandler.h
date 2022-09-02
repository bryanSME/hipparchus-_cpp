#pragma once
/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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
 //package org.hipparchus.linear;
#include "MatrixUtils.h"
//import java.util.List;
#include <type_traits>
#include "../CalculusFieldElement.hpp"

//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.Field;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;

/** Enumerate to specify how dependent vectors should be handled in
 * {@link Matrix_Utils#orthonormalize(List, double, Dependent_Vectors_Handler)} and
 * {@link Matrix_Utils#orthonormalize(Field, List, Calculus_Field_Element, Dependent_Vectors_Handler)}.
 * @since 2.1
 */
enum Dependent_Vectors_Handler
{
	/** Generate a {@link } if dependent vectors are found. */
	GENERATE_EXCEPTION
	{
		/** {@inherit_doc} */
		//override
		public int manage_dependent(const int index, const List<Real_Vector> basis)
		{
			// generate exception, dependent vectors are forbidden with this settings
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::ZERO_NORM);
		}

		/** {@inherit_doc} */
		//override
		template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
		public  int manage_dependent(const Field<T> field, const int index, const List<Field_Vector<T>> basis)
		{
			// generate exception, dependent vectors are forbidden with this settings
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::ZERO_NORM);
		}
	},
	/** Replace dependent vectors by vectors with norm 0.
	 * <p>
	 * This behavior matches the Wolfram language API. It keeps the
	 * number of output vectors equal to the number of input vectors.
	 * The only two norms output vectors can have are 0 and 1.
	 * </p>
	 */
	ADD_ZERO_VECTOR
	{
		/** {@inherit_doc} */
		//override
		public int manage_dependent(const int index, const List<Real_Vector> basis)
		{
			// add a zero vector, preserving output vector size (and dropping its normalization property)
			basis.set(index, Matrix_Utils::create_real__vector(basis.get(index).get_dimension()));
			return index + 1;
		}

		/** {@inherit_doc} */
		//override
		template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
		public  int manage_dependent(const Field<T> field, const int index, const List<Field_Vector<T>> basis)
		{
			// add a zero vector, preserving output vector size (and dropping its normalization property)
			basis.set(index, Matrix_Utils::create_field_vector(field, basis.get(index).get_dimension()));
			return index + 1;
		}
	},
	/** Ignore dependent vectors.
	 * <p>
	 * This behavior ensures the output vectors form an orthonormal
	 * basis, i.e. all vectors are independent and they all have norm 1.
	 * The number of output vectors may be smaller than the number of
	 * input vectors, this number corresponds to the dimension of the
	 * span of the input vectors.
	 * </p>
	 */
	REDUCE_BASE_TO_SPAN
	{
		/** {@inherit_doc} */
		//override
		public int manage_dependent(const int index, const List<Real_Vector> basis)
		{
		// remove dependent vector
		basis.remove(index);
		return index;
	}

	/** {@inherit_doc} */
	//override
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
		public  int manage_dependent(const Field<T> field, const int index, const List<Field_Vector<T>> basis)
		{
			// remove dependent vector
			basis.remove(index);
			return index;
		}
	};

	/** Manage a dependent vector.
	 * @param index of the vector in the basis
	 * @param basis placeholder for basis vectors
	 * @return next index to manage
	 */
	public virtual int manage_dependent(const int& index, List<Real_Vector> basis);

	/** Manage a dependent vector.
	 * @param <T> type of the vectors components
	 * @param field field to which the vectors belong
	 * @param index of the vector in the basis
	 * @param basis placeholder for basis vectors
	 * @return next index to manage
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public virtual  int manage_dependent(Field<T> field, int index, List<Field_Vector<T>> basis);
};