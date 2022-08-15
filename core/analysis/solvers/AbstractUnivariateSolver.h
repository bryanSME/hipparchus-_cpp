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

#include <exception>
#include "BaseAbstractUnivariateSolver.h"
#include "../UnivariateFunction.h"
#include "UnivariateSolver.h"

/**
 * Base class for solvers.
 *
 */
class Abstract_Univariate_Solver : public Base_Abstract_Univariate_Solver<Univariate_Function>, public Univariate_Solver
{
protected:
    /**
     * Construct a solver with given absolute accuracy.
     *
     * @param absolute_accuracy Maximum absolute error.
     */
    Abstract_Univariate_Solver(const double& absolute_accuracy)
    {
        throw std::exception("not implemented");
        //super(absolute_accuracy);
    }
    /**
     * Construct a solver with given accuracies.
     *
     * @param relative_accuracy Maximum relative error.
     * @param absolute_accuracy Maximum absolute error.
     */
    Abstract_Univariate_Solver(const double& relative_accuracy, const double& absolute_accuracy)
    {
        throw std::exception("not implemented");
        //super(relative_accuracy, absolute_accuracy);
    }
    /**
     * Construct a solver with given accuracies.
     *
     * @param relative_accuracy Maximum relative error.
     * @param absolute_accuracy Maximum absolute error.
     * @param function_value_accuracy Maximum function value error.
     */
    Abstract_Univariate_Solver(const double& relative_accuracy, const double& absolute_accuracy, const double& function_value_accuracy)
    {
        throw std::exception("not implemented");
        //super(relative_accuracy, absolute_accuracy, function_value_accuracy);
    }
};