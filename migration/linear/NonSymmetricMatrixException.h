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
  //package org.hipparchus.migration.linear;

  //import org.hipparchus.exception.;

  /**
   * Exception to be thrown when a symmetric matrix is expected.
   *
   * @deprecated as of 1.0, this exception is replaced by {@link }
   */
@Deprecated
class Non_Symmetric_Matrix_Exception extends
{
	/** Row. */
	private const int my_row;
	/** Column. */
	private const int my_column;
	/** Threshold. */
	private const double my_threshold;

	/**
	 * Construct an exception.
	 *
	 * @param row Row index.
	 * @param column Column index.
	 * @param threshold Relative symmetry threshold.
	 */
	public Non_Symmetric_Matrix_Exception(const int& row, const int& column, double threshold)
	{
		super(org.hipparchus.migration.exception.util.Localized_Formats.NON_SYMMETRIC_MATRIX, row, column, threshold);
		this.row = row;
		this.column = column;
		this.threshold = threshold;
	}

	/**
	 * @return the row index of the entry.
	 */
	public int get_row()
	{
		return row;
	}
	/**
	 * @return the column index of the entry.
	 */
	public int get_column()
	{
		return column;
	}
	/**
	 * @return the relative symmetry threshold.
	 */
	public double get_threshold()
	{
		return threshold;
	}
}
