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
 * Exception to be thrown when a square matrix is expected.
 *
 * @deprecated as of 1.0, this exception is replaced by {@link org.hipparchus.exception.}
 */
@Deprecated
class Non_Square_Matrix_Exception
    extends org.hipparchus.migration.exception.Dimension_Mismatch_Exception 
    {

    /**
     * Construct an exception from the mismatched dimensions.
     *
     * @param wrong Row dimension.
     * @param expected Column dimension.
     */
    public Non_Square_Matrix_Exception(const int& wrong, int expected) 
    {
        super(org.hipparchus.migration.exception.util.Localized_Formats.NON_SQUARE_MATRIX, wrong, expected);
    }
}


