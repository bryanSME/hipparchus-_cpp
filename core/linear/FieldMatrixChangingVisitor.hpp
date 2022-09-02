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

//import org.hipparchus.Field_Element;
#include <type_traits>


  /**
   * Interface defining a visitor for matrix entries.
   *
   * @param <T> the type of the field elements
   */
template<typename T, typename std::enable_if<std::is_base_of<Field_Element, T>::value>::type* = nullptr>
class Field_Matrix_Changing_Visitor
{
	/**
	 * Start visiting a matrix.
	 * <p>This method is called once before any entry of the matrix is visited.</p>
	 * @param rows number of rows of the matrix
	 * @param columns number of columns of the matrix
	 * @param start_row Initial row index
	 * @param end_row Final row index (inclusive)
	 * @param start_column Initial column index
	 * @param end_column Final column index (inclusive)
	 */
	void start(const int& rows, const int& columns, const int& start_row, const int& end_row, const int& start_column, const int& end_column);

	/**
	 * Visit one matrix entry.
	 * @param row row index of the entry
	 * @param column column index of the entry
	 * @param value current value of the entry
	 * @return the value to be set for the entry
	 */
	T visit(const int& row, const int& column, const T& value);

	/**
	 * End visiting a matrix.
	 * <p>This method is called once after all entries of the matrix have been visited.</p>
	 * @return the value that the <code>walk_in_xxx_order</code> must return
	 */
	T end();
};