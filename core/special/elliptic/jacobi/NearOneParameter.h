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

#include "../../../util/SinhCosh.h"
#include "CopolarN.h"
#include "JacobiElliptic.h"

/** Algorithm for computing the principal Jacobi functions for parameters slightly below one.
 * <p>
 * The algorithm for evaluating the functions is based on approximation
 * in terms of hyperbolic functions. It is given in Abramowitz and Stegun, * sections 16.15.
 * </p>
 * @since 2.0
 */
class Near_One_Parameter : public Jacobi_Elliptic 
{
private:
    /** Complementary parameter of the Jacobi elliptic function. */
    const double my_m1;

public:
    /** Simple constructor.
     * @param m parameter of the Jacobi elliptic function (must be one or slightly below one here)
     */
    Near_One_Parameter(const double& m) 
    {
        super(m);
        my_m1 = 1.0 - m;
    }

    /** {@inherit_doc} */
    //override
    Copolar_N values_n(const double& u) 
    {
        const Sinh_Cosh sch  = std::sinh_cosh(u);
        const double sech   =  1.0 / sch.cosh();
        const double t      = sch.sinh() * sech;
        const double factor = 0.25 * my_m1 * (sch.sinh() * sch.cosh()  - u) * sech;
        return Copolar_N(t + factor * sech,  // equation 16.15.1
                            sech - factor * t,  // equation 16.15.2
                            sech + factor * t); // equation 16.15.3
    }

}


