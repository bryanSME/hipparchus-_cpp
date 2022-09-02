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

#include <type_traits>
  //import org.hipparchus.Field_Element;

  /**
   * Default implementation of the {@link Field_Matrix_Changing_Visitor} interface.
   * <p>
   * This class is a convenience to create custom visitors without defining all
   * methods. This class provides default implementations that do nothing.
   * </p>
   *
   * @param <T> the type of the field elements
   */
template<typename T, typename std::enable_if<std::is_base_of<Field_Element<T>, T>::value>::type* = nullptr>
class DefaultField_Matrix_Changing_Visitor
	: public Field_Matrix_Changing_Visitor<T>
{
private:
	/** Zero element of the field. */
	const T my_zero;

public:
	/** Build a instance.
	 * @param zero additive identity of the field
	 */
	DefaultField_Matrix_Changing_Visitor(const T& zero) : my_zero{ zero }
	{
	}

	/** {@inherit_doc} */
	//override
	void start(const int& rows, const int& columns, const int& start_row, const int& end_row, const int& start_column, const int& end_column)
	{
	}

	/** {@inherit_doc} */
	//override
	T visit(const int& row, const int& column, T value)
	{
		return value;
	}

	/** {@inherit_doc} */
	//override
	T end() const
	{
		return my_zero;
	}
};