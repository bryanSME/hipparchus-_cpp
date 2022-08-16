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
 //package org.hipparchus.util;

 //import java.util.Arrays;

 //import org.hipparchus.Calculus_Field_Element;
 //import org.hipparchus.Field;
 //import org.hipparchus.exception.;
#include <vector>
#include <type_traits>
#include "../CalculusFieldElement.hpp"

/**
 * This class allows to perform the same computation of all components of a Tuple at once.
 * @param <T> the type of the field elements
 * @since 1.2
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
class Field_Tuple : public Calculus_Field_Element<Field_Tuple<T>>
{
private:
	/** Components of the tuple. */
	const std::vector<T> my_values;

	/** Field the instance belongs to. */
	const Field_Tuple_Field<T> my_field;

	/** Creates a instance from its components.
	 * @param field field the instance belongs to
	 * @param x components of the tuple (beware, it is <em>not</em> copied, it is shared with caller)
	 */
	Field_Tuple(const Field_Tuple_Field<T> field, const std::vector<T> x) { // NOPMD - storing user-supplied array is intentional and documented here
		this.values = x;
		this.field = field;
	}

public:
	/** Creates a instance from its components.
	 * @param x components of the tuple
	 */
	 //@Safe_Varargs
	Field_Tuple(const T... x)
	{
		this(new Field_Tuple_Field<>(x[0].get_field(), x.size()), x.clone());
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> new_instance(const double value)
	{
		////@Suppress_Warnings("unchecked")
		const Field_Tuple<T> t = Field_Tuple<>(field, (std::vector<T>) Math_Arrays::build_array(field, values.size()));
		Arrays.fill(t.values, value);
		return t;
	}

	/** Get the dimension of the tuple.
	 * @return dimension of the tuple
	 */
	int get_dimension() const
	{
		return my_values.size();
	}

	/** Get one component of the tuple.
	 * @param index index of the component, between 0 and {@link #get_dimension() get_dimension()} - 1
	 * @return value of the component
	 */
	T get_component(const int& index) const
	{
		return my_values[index];
	}

	/** Get all components of the tuple.
	 * @return all components
	 */
	std::vector<T> get_components() const
	{
		return my_values;
	}

	/** {@inherit_doc} */
	//override
	Field<Field_Tuple<T>> get_field() const
	{
		return my_field;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> add(const Field_Tuple<T> a)
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].add(a.values[i]);
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> subtract(const Field_Tuple<T> a)
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].subtract(a.values[i]);
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> negate()
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].negate();
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> multiply(const Field_Tuple<T> a)
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].multiply(a.values[i]);
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> multiply(const int& n)
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].multiply(n);
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> divide(const Field_Tuple<T> a)
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].divide(a.values[i]);
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> reciprocal()
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].reciprocal();
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	bool equals(const Object& obj)
	{
		if (dynamic_cast<const Field_Tuple*>(*obj) != nullptr)
		{
			////@Suppress_Warnings("unchecked")
			const Field_Tuple<T> that = (Field_Tuple<T>) obj;
			if (get_dimension() == that.get_dimension())
			{
				bool equals = true;
				for (int i{}; i < values.size(); ++i)
				{
					equals &= values[i].equals(that.values[i]);
				}
				return equals;
			}
		}
		return false;
	}

	/** {@inherit_doc} */
	//override
	int hash_code()
	{
		return  0x58f61de5 + Arrays.hash_code(values);
	}

	/** {@inherit_doc} */
	//override
	double get_real()
	{
		return my_values[0].get_real();
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> add(const double& a)
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].add(a);
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> subtract(const double& a)
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].subtract(a);
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> multiply(const double& a)
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].multiply(a);
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> divide(const double& a)
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].divide(a);
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> remainder(const double& a)
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].remainder(a);
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> remainder(const Field_Tuple<T> a)
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].remainder(a.values[i]);
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> abs()
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].abs();
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> ceil()
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].ceil();
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> floor()
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].floor();
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> rint()
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].rint();
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> sign()
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].sign();
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> copy_sign(const Field_Tuple<T> sign)
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].copy_sign(sign.values[i]);
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> copy_sign(const double sign)
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].copy_sign(sign);
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> scalb(const int& n)
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].scalb(n);
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> ulp()
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].ulp();
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> hypot(const Field_Tuple<T> y)
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].hypot(y.values[i]);
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> sqrt()
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].sqrt();
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> cbrt()
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].cbrt();
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> root_n(const int& n)
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].root_n(n);
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> pow(const double& p)
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].pow(p);
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> pow(const int& n)
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].pow(n);
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> pow(const Field_Tuple<T> e)
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].pow(e.values[i]);
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> exp()
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].exp();
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> expm1()
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].expm1();
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> log()
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].log();
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> log1p()
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].log1p();
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> log10()
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].log10();
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> cos()
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].cos();
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> sin()
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].sin();
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Sin_Cos<Field_Tuple<T>> sin_cos()
	{
		const Field_Tuple<T> sin = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		const Field_Tuple<T> cos = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			const Field_Sin_Cos<T> sc = Sin_Cos(values[i]);
			sin.values[i] = sc.sin();
			cos.values[i] = sc.cos();
		}
		return Field_Sin_Cos<>(sin, cos);
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> tan()
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].tan();
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> acos()
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].acos();
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> asin()
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].asin();
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> atan()
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].atan();
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> atan2(const Field_Tuple<T> x)
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].atan2(x.values[i]);
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> cosh()
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].cosh();
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> sinh()
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].sinh();
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Sinh_Cosh<Field_Tuple<T>> sinh_cosh()
	{
		const Field_Tuple<T> sinh = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		const Field_Tuple<T> cosh = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			const Field_Sinh_Cosh<T> sc = std::sinh_cosh(values[i]);
			sinh.values[i] = sc.sinh();
			cosh.values[i] = sc.cosh();
		}
		return Field_Sinh_Cosh<>(sinh, cosh);
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> tanh()
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].tanh();
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> acosh()
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].acosh();
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> asinh()
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].asinh();
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> atanh()
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].atanh();
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> to_degrees()
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].to_degrees();
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> to_radians()
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = values[i].to_radians();
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> linear_combination(const Field_Tuple<T>[] a, const Field_Tuple<T>[] b)
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		Math_Utils::check_dimension(a.size(), b.size());
		const std::vector<T> aT = Math_Arrays::build_array(values[0].get_field(), a.size());
		const std::vector<T> b_t = Math_Arrays::build_array(values[0].get_field(), b.size());
		for (int i{}; i < values.size(); ++i)
		{
			for (int j{}; j < a.size(); ++j)
			{
				aT[j] = a[j].values[i];
				b_t[j] = b[j].values[i];
			}
			result.values[i] = aT[0].linear_combination(aT, b_t);
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> linear_combination(const std::vector<double> a, const Field_Tuple<T>[] b)
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		Math_Utils::check_dimension(a.size(), b.size());
		const std::vector<T> b_t = Math_Arrays::build_array(values[0].get_field(), b.size());
		for (int i{}; i < values.size(); ++i)
		{
			for (int j{}; j < a.size(); ++j)
			{
				b_t[j] = b[j].values[i];
			}
			result.values[i] = b_t[0].linear_combination(a, b_t);
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> linear_combination(const Field_Tuple<T> a1, const Field_Tuple<T> b1, const Field_Tuple<T> a2, const Field_Tuple<T> b2)
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = a1.values[0].linear_combination(a1.values[i], b1.values[i], a2.values[i], b2.values[i]);
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> linear_combination(const double& a1, const Field_Tuple<T> b1, const double& a2, const Field_Tuple<T> b2)
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = b1.values[0].linear_combination(a1, b1.values[i], a2, b2.values[i]);
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> linear_combination(const Field_Tuple<T> a1, const Field_Tuple<T> b1, const Field_Tuple<T> a2, const Field_Tuple<T> b2, const Field_Tuple<T> a3, const Field_Tuple<T> b3)
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = a1.values[0].linear_combination(a1.values[i], b1.values[i], a2.values[i], b2.values[i], a3.values[i], b3.values[i]);
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> linear_combination(const double& a1, const Field_Tuple<T> b1, const double& a2, const Field_Tuple<T> b2, const double& a3, const Field_Tuple<T> b3)
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = b1.values[0].linear_combination(a1, b1.values[i], a2, b2.values[i], a3, b3.values[i]);
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> linear_combination(const Field_Tuple<T> a1, const Field_Tuple<T> b1, const Field_Tuple<T> a2, const Field_Tuple<T> b2, const Field_Tuple<T> a3, const Field_Tuple<T> b3, const Field_Tuple<T> a4, const Field_Tuple<T> b4)
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = a1.values[0].linear_combination(a1.values[i], b1.values[i], a2.values[i], b2.values[i], a3.values[i], b3.values[i], a4.values[i], b4.values[i]);
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> linear_combination(const double& a1, const Field_Tuple<T> b1, const double& a2, const Field_Tuple<T> b2, const double& a3, const Field_Tuple<T> b3, const double& a4, const Field_Tuple<T> b4)
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		for (int i{}; i < values.size(); ++i)
		{
			result.values[i] = b1.values[0].linear_combination(a1, b1.values[i], a2, b2.values[i], a3, b3.values[i], a4, b4.values[i]);
		}
		return result;
	}

	/** {@inherit_doc} */
	//override
	Field_Tuple<T> get_pi()
	{
		const Field_Tuple<T> result = Field_Tuple<>(field, Math_Arrays::build_array(values[0].get_field(), values.size()));
		Arrays.fill(result.values, values[0].get_pi());
		return result;
	}

	/** Field for {link Field_Tuple} instances.
	 * @param <T> the type of the field elements
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	//private static class Field_Tuple_Field : Field<Field_Tuple<T>>
	//{
	//    /** Constant function evaluating to 0.0. */
	//    private const Field_Tuple<T> zero;

	//    /** Constant function evaluating to 1.0. */
	//    private const Field_Tuple<T> one;

	//    /** Simple constructor.
	//     * @param field field to which the elements belong
	//     * @param dimension dimension of the tuple
	//     */
	//    Field_Tuple_Field(const Field<T> field, const int& dimension)
	//    {
	//        const std::vector<T> zero_data = Math_Arrays::build_array(field, dimension);
	//        Arrays.fill(zero_data, field.get_zero());
	//        const std::vector<T> one_data  = Math_Arrays::build_array(field, dimension);
	//        Arrays.fill(one_data, field.get_one());
	//        this.zero = Field_Tuple<>(this, zero_data);
	//        this.one  = Field_Tuple<>(this, one_data);
	//    }

	//    /** {@inherit_doc} */
	//    //override
	//    public Field_Tuple<T> get_zero()
	//    {
	//        return zero;
	//    }

	//    /** {@inherit_doc} */
	//    //override
	//    public Field_Tuple<T> get_one()
	//    {
	//        return one;
	//    }

	//    /** {@inherit_doc} */
	//    ////@Suppress_Warnings("unchecked")
	//    //override
	//    public Class<Field_Tuple<T>> get_runtime_class()
	//    {
	//        return (Class<Field_Tuple<T>>) zero.get_class();
	//    }

	//    /** {@inherit_doc} */
	//    //override
	//    public bool equals(const Object& other)
	//    {
	//        if (dynamic_cast<const Field_Tuple*>(*other) != nullptr)
	//        {
	//            ////@Suppress_Warnings("unchecked")
	//            const Field_Tuple_Field<T> that = (Field_Tuple_Field<T>) other;
	//            return zero.get_dimension() == that.zero.get_dimension();
	//        }
	//        else
	//        {
	//            return false;
	//        }
	//    }

	//    /** {@inherit_doc} */
	//    //override
	//    public int hash_code()
	//    {
	//        return 0xb4a533e1 ^ zero.get_dimension();
	//    }

	//}
};