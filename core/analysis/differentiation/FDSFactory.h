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
 //package org.hipparchus.analysis.differentiation;

 //import org.hipparchus.Calculus_Field_Element;
 //import org.hipparchus.Field;
 //import org.hipparchus.exception.Localized_Core_Formats;
 //import org.hipparchus.exception.;
 //import org.hipparchus.util.Math_Arrays;

 /** Factory for {@link Field_Derivative_Structure}.
  * <p>This class is a factory for {@link Field_Derivative_Structure} instances.</p>
  * <p>Instances of this class are guaranteed to be immutable.</p>
  * @see Field_Derivative_Structure
  * @param <T> the type of the function parameters and value
  */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
class FDS_Factory
{
	/** Compiler for the current dimensions. */
	private const DS_Compiler compiler;

	/** Field the value and parameters of the function belongs to. */
	private const Field<T> value_field;

	/** Field the {@link Field_Derivative_Structure} instances belong to. */
	private const Derivative_Field<T> derivative_field;

	/** Simple constructor.
	 * @param value_field field for the function parameters and value
	 * @param parameters number of free parameters
	 * @param order derivation order
	 */
	public FDS_Factory(const Field<T> value_field, const int parameters, const int order)
	{
		this.compiler = DS_Compiler.get_compiler(parameters, order);
		this.value_field = value_field;
		this.derivative_field = Derivative_Field<>(constant(value_field.get_zero()), constant(value_field.get_one()), constant(value_field.get_zero().get_pi()));
	}

	/** Get the {@link Field} the value and parameters of the function belongs to.
	 * @return {@link Field} the value and parameters of the function belongs to
	 */
	public Field<T> get_value_field()
	{
		return value_field;
	}

	/** Get the {@link Field} the {@link Field_Derivative_Structure} instances belong to.
	 * @return {@link Field} the {@link Field_Derivative_Structure} instances belong to
	 */
	public Derivative_Field<T> get_derivative_field()
	{
		return derivative_field;
	}

	/** Build a {@link Field_Derivative_Structure} representing a constant value.
	 * @param value value of the constant
	 * @return a {@link Field_Derivative_Structure} representing a constant value
	 */
	public Field_Derivative_Structure<T> constant(double value)
	{
		return constant(value_field.get_zero().add(value));
	}

	/** Build a {@link Field_Derivative_Structure} representing a constant value.
	 * @param value value of the constant
	 * @return a {@link Field_Derivative_Structure} representing a constant value
	 */
	public Field_Derivative_Structure<T> constant(const T value)
	{
		const Field_Derivative_Structure<T> fds = Field_Derivative_Structure<>(this);
		fds.set_derivative_component(0, value);
		return fds;
	}

	/** Build a {@link Field_Derivative_Structure} representing a variable.
	 * <p>Instances built using this method are considered
	 * to be the free variables with respect to which differentials
	 * are computed. As such, their differential with respect to
	 * themselves is +1.</p>
	 * @param index index of the variable (from 0 to
	 * {@link #get_compiler()}.{@link DS_Compiler#get_free_parameters() get_free_parameters()} - 1)
	 * @param value value of the variable
	 * @return a {@link Field_Derivative_Structure} representing a variable
	 * @exception  if index if greater or
	 * equal to {@link #get_compiler()}.{@link DS_Compiler#get_free_parameters() get_free_parameters()}.
	 */
	public Field_Derivative_Structure<T> variable(const int index, const T value)

	{
		if (index >= get_compiler().get_free_parameters())
		{
			throw std::exception("not implmented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_LARGE_BOUND_EXCLUDED, index, get_compiler().get_free_parameters());
		}

		const Field_Derivative_Structure<T> fds = Field_Derivative_Structure<>(this);
		fds.set_derivative_component(0, value);

		if (get_compiler().get_order() > 0)
		{
			// the derivative of the variable with respect to itself is 1.
			fds.set_derivative_component(DS_Compiler.get_compiler(index, get_compiler().get_order()).get_size(), value_field.get_one());
		}

		return fds;
	}

	/** Build a {@link Field_Derivative_Structure} representing a variable.
	 * <p>Instances built using this method are considered
	 * to be the free variables with respect to which differentials
	 * are computed. As such, their differential with respect to
	 * themselves is +1.</p>
	 * @param index index of the variable (from 0 to
	 * {@link #get_compiler()}.{@link DS_Compiler#get_free_parameters() get_free_parameters()} - 1)
	 * @param value value of the variable
	 * @return a {@link Field_Derivative_Structure} representing a variable
	 * @exception  if index if greater or
	 * equal to {@link #get_compiler()}.{@link DS_Compiler#get_free_parameters() get_free_parameters()}.
	 */
	public Field_Derivative_Structure<T> variable(const int index, const double value)

	{
		if (index >= get_compiler().get_free_parameters())
		{
			throw std::exception("not implmented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_LARGE_BOUND_EXCLUDED, index, get_compiler().get_free_parameters());
		}

		const Field_Derivative_Structure<T> fds = Field_Derivative_Structure<>(this);
		fds.set_derivative_component(0, value_field.get_zero().new_instance(value));

		if (get_compiler().get_order() > 0)
		{
			// the derivative of the variable with respect to itself is 1.
			fds.set_derivative_component(DS_Compiler.get_compiler(index, get_compiler().get_order()).get_size(), value_field.get_one());
		}

		return fds;
	}

	/** Build a {@link Field_Derivative_Structure} from all its derivatives.
	 * @param derivatives derivatives sorted according to
	 * {@link DS_Compiler#get_partial_derivative_index(int...)}
	 * @return  {@link Field_Derivative_Structure} with specified derivatives
	 * @exception  if derivatives array does not match the
	 * {@link DS_Compiler#get_size() size} expected by the compiler
	 * @exception  if order is too large
	 * @see Field_Derivative_Structure#get_all_derivatives()
	 */
	 //@Safe_Varargs
	public const Field_Derivative_Structure<T> build(const T ... derivatives)
	{
		const std::vector<T> data = build_array();
		if (derivatives.size() != data.size())
		{
			throw std::exception("not implmented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, derivatives.size(), data.size());
		}
		System.arraycopy(derivatives, 0, data, 0, data.size());

		return Field_Derivative_Structure<>(this, data);
	}

	/** Build a {@link Field_Derivative_Structure} from all its derivatives.
	 * @param derivatives derivatives sorted according to
	 * {@link DS_Compiler#get_partial_derivative_index(int...)}
	 * @return  {@link Field_Derivative_Structure} with specified derivatives
	 * @exception  if derivatives array does not match the
	 * {@link DS_Compiler#get_size() size} expected by the compiler
	 * @exception  if order is too large
	 * @see Field_Derivative_Structure#get_all_derivatives()
	 */
	public Field_Derivative_Structure<T> build(const double ... derivatives)

	{
		const std::vector<T> data = build_array();
		if (derivatives.size() != data.size())
		{
			throw std::exception("not implmented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, derivatives.size(), data.size());
		}
		for (int i{}; i < data.size(); ++i)
		{
			data[i] = value_field.get_zero().add(derivatives[i]);
		}

		return Field_Derivative_Structure<>(this, data);
	}

	/** Build a {@link Field_Derivative_Structure} with an uninitialized array.
	 * <p>This method is intended only for Field_Derivative_Structure internal use.</p>
	 * @return a {@link Field_Derivative_Structure} with an uninitialized array
	 */
	Field_Derivative_Structure<T> build()
	{
		return Field_Derivative_Structure<>(this);
	}

	/** Build an uninitialized array for derivatives data.
	 * @return uninitialized array for derivatives data
	 */
	private std::vector<T> build_array()
	{
		return Math_Arrays::build_array(value_field, compiler.get_size());
	}

	/** Get the compiler for the current dimensions.
	 * @return compiler for the current dimensions
	 */
	public DS_Compiler get_compiler()
	{
		return compiler;
	}

	/** Check rules set compatibility.
	 * @param factory other factory field to check against instance
	 * @exception  if number of free parameters or orders are inconsistent
	 */
	void check_compatibility(const FDS_Factory<T> factory)
	{
		compiler.check_compatibility(factory.compiler);
	}

	/** Field for {link Field_Derivative_Structure} instances.
	 * @param <T> the type of the function parameters and value
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public static class Derivative_Field : Field<Field_Derivative_Structure<T>>
	{
		/** Constant function evaluating to 0.0. */
		private const Field_Derivative_Structure<T> zero;

		/** Constant function evaluating to 1.0. */
		private const Field_Derivative_Structure<T> one;

		/** Constant function evaluating to \xcf\x80. */
		private const Field_Derivative_Structure<T> pi;

		/** Simple constructor.
		 * @param zero constant function evaluating to 0.0
		 * @param one constant function evaluating to 1.0
		 * @param pi constant function evaluating to \xcf\x80
		 */
		Derivative_Field(const Field_Derivative_Structure<T> zero, const Field_Derivative_Structure<T> one, const Field_Derivative_Structure<T> pi)
		{
			this.zero = zero;
			this.one = one;
			this.pi = pi;
		}

		/** {@inherit_doc} */
		//override
		public Field_Derivative_Structure<T> get_zero()
		{
			return zero;
		}

		/** {@inherit_doc} */
		//override
		public Field_Derivative_Structure<T> get_one()
		{
			return one;
		}

		/** Get the Archimedes constant \xcf\x80.
		 * <p>
		 * Archimedes constant is the ratio of a circle's circumference to its diameter.
		 * </p>
		 * @return Archimedes constant \xcf\x80
		 * @since 2.0
		 */
		public Field_Derivative_Structure<T> get_pi()
		{
			return pi;
		}

		/** {@inherit_doc} */
		//@Suppress_Warnings("unchecked")
		//override
		public Class<Field_Derivative_Structure<T>> get_runtime_class()
		{
			return (Class<Field_Derivative_Structure<T>>) zero.get_class();
		}

		/** {@inherit_doc} */
		//override
		public bool equals(const Object& other)
		{
			if (this == other)
			{
				return true;
			}
			if (dynamic_cast<const Derivative_Field*>(*other) != nullptr)
			{
				FDS_Factory<T> lhs_factory = zero.get_factory();
				FDS_Factory< ? > rhs_factory = ((Derivative_Field< ? >) other).zero.get_factory();
				return lhs_factory.compiler == rhs_factory.compiler &&
					lhs_factory.value_field.equals(rhs_factory.value_field);
			}
			return false;
		}

		/** {@inherit_doc} */
		//override
		public int hash_code()
		{
			const DS_Compiler compiler = zero.get_factory().get_compiler();
			return 0x58d35de8 ^ (compiler.get_free_parameters() << 16 & compiler.get_order());
		}
	}
}
