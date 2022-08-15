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
 * Exception to be thrown when a self-adjoint {@link org.hipparchus.linear.Real_Linear_Operator}
 * is expected.
 * sin_ce the coefficients of the matrix are not accessible, the most
 * general definition is used to check that A is not self-adjoint, i.e.
 * there exist x and y such as {@code | x' A y - y' A x | >= eps}, * where {@code eps} is a user-specified tolerance, and {@code x'}
 * denotes the transpose of {@code x}.
 * In the terminology of this exception, {@code A} is the "offending"
 * linear operator, {@code x} and {@code y} are the first and second
 * "offending" vectors, respectively.
 *
 * @deprecated as of 1.0, this exception is replaced by {@link }
 */
@Deprecated
class Non_Self_Adjoint_Operator_Exception  
{

    /** Creates a instance of this class. */
    public Non_Self_Adjoint_Operator_Exception() 
    {
        super(org.hipparchus.migration.exception.util.Localized_Formats.NON_SELF_ADJOINT_OPERATOR);
    }
}


