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
//package org.hipparchus.optim.nonlinear.scalar;

//import java.util.Array_list;
//import java.util.Collections;
//import java.util.Comparator;
//import java.util.List;

//import org.hipparchus.exception.;
//import org.hipparchus.exception.Null_Argument_Exception;
//import org.hipparchus.optim.BaseMulti_startMultivariate_Optimizer;
//import org.hipparchus.optim.Point_valuePair;
//import org.hipparchus.random.Random_Vector_Generator;

/**
 * Multi-start optimizer.
 *
 * This class wraps an optimizer in order to use it several times in
 * turn with different starting points (trying to avoid being trapped
 * in a local extremum when looking for a global one).
 *
 */
class Multi_startMultivariate_Optimizer
    extends BaseMulti_startMultivariate_Optimizer<Point_valuePair> 
    {
    /** Underlying optimizer. */
    private const Multivariate_Optimizer optimizer;
    /** Found optima. */
    private const List<Point_valuePair> optima;

    /**
     * Create a multi-start optimizer from a single-start optimizer.
     *
     * @param optimizer Single-start optimizer to wrap.
     * @param starts Number of starts to perform.
     * If {@code starts == 1}, the result will be same as if {@code optimizer}
     * is called directly.
     * @param generator Random vector generator to use for restarts.
     * @Null_Argument_Exception if {@code optimizer} or {@code generator}
     * is {@code NULL}.
     * @ if {@code starts < 1}.
     */
    public Multi_startMultivariate_Optimizer(const Multivariate_Optimizer optimizer, const int starts, const Random_Vector_Generator generator)
        , Null_Argument_Exception 
        {
        super(optimizer, starts, generator);
        this.optimizer = optimizer;
        this.optima   = Array_list<>();
    }

    /**
     * {@inherit_doc}
     */
    //override
    public Point_valuePair[] get_optima() 
    {
        Collections.sort(optima, get_pair_comparator());
        return optima.to_array(new Point_valuePair[0]);
    }

    /**
     * {@inherit_doc}
     */
    //override
    protected void store(Point_valuePair optimum) 
    {
        optima.add(optimum);
    }

    /**
     * {@inherit_doc}
     */
    //override
    protected void clear() 
    {
        optima.clear();
    }

    /**
     * @return a comparator for sorting the optima.
     */
    private Comparator<Point_valuePair> get_pair_comparator() 
    {
        return Comparator<Point_valuePair>() 
        {
            /** {@inherit_doc} */
            //override
            public int compare(const Point_valuePair o1, const Point_valuePair o2) 
            {
                if (o1 == NULL) 
                {
                    return (o2 == NULL) ? 0 : 1;
                }
else if (o2 == NULL) 
                {
                    return -1;
                }
                const double v1 = o1.get_value();
                const double v2 = o2.get_value();
                return (optimizer.get_goal_type() == Goal_Type.MINIMIZE) ?
                    Double.compare(v1, v2) : Double.compare(v2, v1);
            }
        };
    }
}


