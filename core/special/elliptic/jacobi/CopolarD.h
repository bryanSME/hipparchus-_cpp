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

#include "CopolarN.h"

/** Copolar trio with pole at point d in Glaisher’s Notation.
 * <p>
 * This is a container for the three subsidiary Jacobi elliptic functions
 * {@code nd(u|m)}, {@code sd(u|m)}, and {@code cd(u|m)}.
 * </p>
 * @since 2.0
 */
class Copolar_D
{
private:
    /** Value of the nd function. */
    const double my_nd;

    /** Value of the sd function. */
    const double my_sd;

    /** Value of the cd function. */
    const double my_cd;

public:
    /** Simple constructor.
     * @param trio_n copolar trio with pole at point n in Glaisher’s Notation
     */
    Copolar_D(const Copolar_N& trio_n)
        : 
        my_nd{ 1.0 / trio_n.dn() },
        my_sd{ my_nd * trio_n.sn() },
        my_cd{ my_nd * trio_n.cn() }
    {};

    /** Get the value of the nd function.
     * @return nd(u|m)
     */
    double nd() const
    {
        return my_nd;
    }

    /** Get the value of the sd function.
     * @return sd(u|m)
     */
    double sd() const
    {
        return my_sd;
    }

    /** Get the value of the cd function.
     * @return cd(u|m)
     */
    double cd() const
    {
        return my_cd;
    }
};