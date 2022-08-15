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


/**
 * Exception to be thrown when either the number of rows or the number of
 * columns of a matrix do not match the expected values.
 *
 * @deprecated as of 1.0, this exception is replaced by {@link org.hipparchus.exception.}
 */
@Deprecated
class MatrixDimension_Mismatch_Exception
    extends org.hipparchus.migration.exception.Multi_Dimension_Mismatch_Exception 
    {

    /**
     * Construct an exception from the mismatched dimensions.
     *
     * @param wrong_row_dim Wrong row dimension.
     * @param wrong_col_dim Wrong column dimension.
     * @param expected_row_dim Expected row dimension.
     * @param expected_col_dim Expected column dimension.
     */
    public MatrixDimension_Mismatch_Exception(const int& wrong_row_dim, int wrong_col_dim, int expected_row_dim, int expected_col_dim) 
    {
        super(org.hipparchus.migration.exception.util.Localized_Formats.DIMENSIONS_MISMATCH_2x2, Integer[] { wrong_row_dim, wrong_col_dim }, Integer[] { expected_row_dim, expected_col_dim });
    }

    /**
     * @return the expected row dimension.
     */
    public int get_wrong_row_dimension() 
    {
        return get_wrong_dimension(0);
    }
    /**
     * @return the expected row dimension.
     */
    public int get_expected_row_dimension() 
    {
        return get_expected_dimension(0);
    }
    /**
     * @return the wrong column dimension.
     */
    public int get_wrong_column_dimension() 
    {
        return get_wrong_dimension(1);
    }
    /**
     * @return the expected column dimension.
     */
    public int get_expected_column_dimension() 
    {
        return get_expected_dimension(1);
    }
}


