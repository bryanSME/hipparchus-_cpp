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

  //import org.hipparchus.Field_Element;

  /**
   * This interface defines a visitor for the entries of a vector. Visitors
   * implementing this interface do not alter the entries of the vector being
   * visited.
   *
   * @param <T> the type of the field elements
   */
class Field_VectorPreservingVisitor<T extends Field_Element< ? >>
{
	/**
	 * Start visiting a vector. This method is called once, before any entry
	 * of the vector is visited.
	 *
	 * @param dimension the size of the vector
	 * @param start the index of the first entry to be visited
	 * @param end the index of the last entry to be visited (inclusive)
	 */
	void start(const int& dimension, int start, int end);

	/**
	 * Visit one entry of the vector.
	 *
	 * @param index the index of the entry being visited
	 * @param value the value of the entry being visited
	 */
	void visit(const int& index, T value);

	/**
	 * End visiting a vector. This method is called once, after all entries of
	 * the vector have been visited.
	 *
	 * @return the value returned after visiting all entries
	 */
	T end();
}
