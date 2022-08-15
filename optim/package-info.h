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
 * <p>
 *  Generally, optimizers are algorithms that will either
 *  {@link org.hipparchus.optim.nonlinear.scalar.Goal_Type#MINIMIZE minimize} or
 *  {@link org.hipparchus.optim.nonlinear.scalar.Goal_Type#MAXIMIZE maximize}
 *  a scalar function, called the
 *  {@link org.hipparchus.optim.nonlinear.scalar.Objective_Function <em>objective
 *  function</em>}.
 *  <br/>
 *  For some scalar objective functions the gradient can be computed (analytically
 *  or numerically). Algorithms that use this knowledge are defined in the
 *  {@link org.hipparchus.optim.nonlinear.scalar.gradient} //package.
 *  The algorithms that do not need this additional information are located in
 *  the {@link org.hipparchus.optim.nonlinear.scalar.noderiv} //package.
 * </p>
 *
 * <p>
 *  Some problems are solved more efficiently by algorithms that, instead of an
 *  objective function, need access to all the observations.
 *  Such methods are implemented in the fitting module.
 * </p>
 *
 * <p>
 *  This //package provides common functionality for the optimization algorithms.
 *  Abstract classes ({@link org.hipparchus.optim.Base_Optimizer} and
 *  {@link org.hipparchus.optim.BaseMultivariate_Optimizer}) contain
 *  boiler-plate code for storing {@link org.hipparchus.optim.Max_Eval
 *  evaluations} and {@link org.hipparchus.optim.Max_Iter iterations}
 *  counters and a user-defined
 *  {@link org.hipparchus.optim.Convergence_Checker convergence checker}.
 * </p>
 *
 * <p>
 *  For each of the optimizer types, there is a special implementation that
 *  wraps an optimizer instance and provides a "multi-start" feature: it calls
 *  the underlying optimizer several times with different starting points and
 *  returns the best optimum found, or all optima if so desired.
 *  This could be useful to avoid being trapped in a local extremum.
 * </p>
 */
//package org.hipparchus.optim;


