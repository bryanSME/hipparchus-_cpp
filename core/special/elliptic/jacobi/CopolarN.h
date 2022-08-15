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

/** Copolar trio with pole at point n in Glaisher’s Notation.
 * <p>
 * This is a container for the three principal Jacobi elliptic functions
 * {@code sn(u|m)}, {@code cn(u|m)}, and {@code dn(u|m)}.
 * </p>
 * @since 2.0
 */
class Copolar_N
{
private:
    /** Value of the sn function. */
    const double my_sn;

    /** Value of the cn function. */
    const double my_cn;

    /** Value of the dn function. */
    const double my_dn;

public:
    /** Simple constructor.
     * @param sn value of the sn function
     * @param cn value of the cn function
     * @param dn value of the dn function
     */
    Copolar_N(const double& sn, const double& cn, const double& dn)
        :
        my_sn{ sn },
        my_cn{ cn },
        my_dn{ dn }
    {};


    /** Get the value of the sn function.
     * @return sn(u|m)
     */
    double sn() const
    {
        return my_sn;
    }

    /** Get the value of the cn function.
     * @return cn(u|m)
     */
    double cn() const
    {
        return my_cn;
    }

    /** Get the value of the dn function.
     * @return dn(u|m)
     */
    double dn() const
    {
        return my_dn;
    }
};