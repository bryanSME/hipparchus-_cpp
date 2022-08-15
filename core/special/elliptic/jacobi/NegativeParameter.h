#pragma once
/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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
//package org.hipparchus.special.elliptic.jacobi;

#include <cmath>
#include "CopolarD.h"
#include "JacobiElliptic.h"
#include "JacobiEllipticBuilder.hpp"


/** Algorithm for computing the principal Jacobi functions for negative parameter m.
 * <p>
 * The rules for negative parameter change are given in Abramowitz and Stegun, section 16.10.
 * </p>
 * @since 2.0
 */
class Negative_Parameter extends Jacobi_Elliptic 
{
private:
    /** Algorithm to use for the positive parameter. */
    const Jacobi_Elliptic my_algorithm;

    /** Input scaling factor. */
    const double my_input_scale;

    /** output scaling factor. */
    const double my_output_scale;

public:
    /** Simple constructor.
     * @param m parameter of the Jacobi elliptic function (must be negative here)
     */
    Negative_Parameter(const double& m) 
    {
        super(m);
        const auto om_m = 1.0 - m;
        my_algorithm = Jacobi_Elliptic_Builder.build(-m / om_m);
        my_input_scale = std::sqrt(om_m);
        my_output_scale = 1.0 / input_scale;
    }

    /** {@inherit_doc} */
    //override
    Copolar_N values_n(const double& u) 
    {
        const auto trio_d = Copolar_D(algorithm.values_n(u * input_scale));
        return Copolar_N(output_scale * trio_d.sd(), trio_d.cd(), trio_d.nd());
    }

}


