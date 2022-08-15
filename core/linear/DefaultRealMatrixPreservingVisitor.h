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

/**
 * Default implementation of the {@link Real_Matrix_Preserving_Visitor} interface.
 * <p>
 * This class is a convenience to create custom visitors without defining all
 * methods. This class provides default implementations that do nothing.
 * </p>
 *
 */
class DefaultReal_Matrix_Preserving_Visitor : Real_Matrix_Preserving_Visitor 
{
    /** {@inherit_doc} */
    //override
    public void start(const int& rows, int columns, int start_row, int end_row, int start_column, int end_column) 
    {
    }

    /** {@inherit_doc} */
    //override
    public void visit(const int& row, const int& column, double value) {}

    /** {@inherit_doc} */
    //override
    public double end() 
    {
        return 0;
    }
}

