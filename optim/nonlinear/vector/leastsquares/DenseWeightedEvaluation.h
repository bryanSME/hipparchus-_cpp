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
//package org.hipparchus.optim.nonlinear.vector.leastsquares;

//import org.hipparchus.linear.Real_Matrix;
//import org.hipparchus.linear.Real_Vector;
//import org.hipparchus.optim.nonlinear.vector.leastsquares.Least_Squares_Problem.Evaluation;

/**
 * Applies a dense weight matrix to an evaluation.
 *
 */
class DenseWeighted_Evaluation extends Abstract_Evaluation 
{

    /** the unweighted evaluation */
    private const Evaluation unweighted;
    /** reference to the weight square root matrix */
    private const Real_Matrix weight_sqrt;

    /**
     * Create a weighted evaluation from an unweighted one.
     *
     * @param unweighted the evalutation before weights are applied
     * @param weight_sqrt the matrix square root of the weight matrix
     */
    DenseWeighted_Evaluation(const Evaluation unweighted, const Real_Matrix weight_sqrt) 
    {
        // weight square root is square, n_r=n_c=number of observations
        super(weight_sqrt.get_column_dimension());
        this.unweighted = unweighted;
        this.weight_sqrt = weight_sqrt;
    }

    /* apply weights */

    /** {@inherit_doc} */
    //override
    public Real_Matrix get_jacobian() 
    {
        return weight_sqrt.multiply(this.unweighted.get_jacobian());
    }

    /** {@inherit_doc} */
    //override
    public Real_Vector get_residuals() 
    {
        return this.weight_sqrt.operate(this.unweighted.get_residuals());
    }

    /* delegate */

    /** {@inherit_doc} */
    //override
    public Real_Vector get_point() 
    {
        return unweighted.get_point();
    }

}


