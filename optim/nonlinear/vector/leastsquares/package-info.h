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
 * This //package provides algorithms that minimize the residuals
 * between observations and model values.
 * The {@link org.hipparchus.optim.nonlinear.vector.leastsquares.Least_Squares_Optimizer
 * least-squares optimizers} minimize the distance (called
 * <em>cost</em> or <em>&chi;<sup>2</sup></em>) between model and
 * observations.
 *
 * <br/>
 * Algorithms in this category need access to a <em>problem</em>
 * (represented by a {@link org.hipparchus.optim.nonlinear.vector.leastsquares.Least_Squares_Problem
 * Least_Squares_Problem}).
 * Such a model predicts a set of values which the algorithm tries to match
 * with a set of given set of observed values.
 * <br/>
 * The problem can be created progressively using a {@link
 * org.hipparchus.optim.nonlinear.vector.leastsquares.Least_Squares_Builder builder} or it can
 * be created at once using a {@link org.hipparchus.optim.nonlinear.vector.leastsquares.Least_Squares_Factory
 * factory}.
 */
//package org.hipparchus.optim.nonlinear.vector.leastsquares;


