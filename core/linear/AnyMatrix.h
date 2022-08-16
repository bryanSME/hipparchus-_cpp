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

  /**
   * Interface defining very basic matrix operations.
   */
class Any_Matrix
{
	/**
	 * Is this a square matrix?
	 * @return true if the matrix is square (row_dimension = column_dimension)
	 */
	virtual bool is_square() = 0;

	/**
	 * Returns the number of rows in the matrix.
	 *
	 * @return row_dimension
	 */
	virtual int get_row_dimension() = 0;

	/**
	 * Returns the number of columns in the matrix.
	 *
	 * @return column_dimension
	 */
	virtual int get_column_dimension() = 0;
};