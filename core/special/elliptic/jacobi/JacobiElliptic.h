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

//import org.hipparchus.special.elliptic.carlson.Carlson_Elliptic_Integral;
//import org.hipparchus.special.elliptic.legendre.Legendre_Elliptic_Integral;
//import org.hipparchus.util.FastMath;
#include "../carlson/CarlsonEllipticIntegral.h"
#include "../legendre/LegendreEllipticIntegral.h"
#include "../jacobi/CopolarN.h"
#include "../jacobi/CopolarC.h"
#include "../jacobi/CopolarD.h"
#include "../jacobi/CopolarS.h"

/** Algorithm computing Jacobi elliptic functions.
 * @since 2.0
 */
class Jacobi_Elliptic 
{
private:
    /** Parameter of the function. */
    const double my_m;

    /** Evaluate inverse of Jacobi elliptic function ps.
     * <p>
     * Here p, q, r are any permutation of the letters c, d, n.
     * </p>
     * @param x value of Jacobi elliptic function {@code ps(u|m)}
     * @param deltaQP Δ⁡(q, p) = q⁣s²⁡(u|m) - p⁣s²(u|m) (equation 19.5.28 of DLMF)
     * @param deltaRP Δ⁡(r, p) = r⁣s²⁡(u|m) - p⁣s²⁡(u|m) (equation 19.5.28 of DLMF)
     * @return u such that {@code x=ps(u|m)}
     * @since 2.1
     */
    double arcps(const double& x, const double deltaQP, const double deltaRP)
    {
        // see equation 19.25.32 in Digital Library of Mathematical Functions
        // https://dlmf.nist.gov/19.25.E32
        const double x2 = x * x;
        return std::copysign(Carlson_Elliptic_Integral::r_f(x2, x2 + deltaQP, x2 + deltaRP), x);
    }

    /** Evaluate inverse of Jacobi elliptic function sp.
     * <p>
     * Here p, q, r are any permutation of the letters c, d, n.
     * </p>
     * @param x value of Jacobi elliptic function {@code sp(u|m)}
     * @param deltaQP Δ⁡(q, p) = q⁣s²⁡(u|m) - p⁣s²(u|m) (equation 19.5.28 of DLMF)
     * @param deltaRP Δ⁡(r, p) = r⁣s²⁡(u|m) - p⁣s²⁡(u|m) (equation 19.5.28 of DLMF)
     * @return u such that {@code x=sp(u|m)}
     * @since 2.1
     */
    double arcsp(const double& x, const double deltaQP, const double deltaRP)
    {
        // see equation 19.25.33 in Digital Library of Mathematical Functions
        // https://dlmf.nist.gov/19.25.E33
        const double x2 = x * x;
        return x * Carlson_Elliptic_Integral::r_f(1, 1 + deltaQP * x2, 1 + deltaRP * x2);
    }

    /** Evaluate inverse of Jacobi elliptic function pq.
     * <p>
     * Here p, q, r are any permutation of the letters c, d, n.
     * </p>
     * @param x value of Jacobi elliptic function {@code pq(u|m)}
     * @param deltaQP Δ⁡(q, p) = q⁣s²⁡(u|m) - p⁣s²(u|m) (equation 19.5.28 of DLMF)
     * @param deltaRQ Δ⁡(r, q) = r⁣s²⁡(u|m) - q⁣s²⁡(u|m) (equation 19.5.28 of DLMF)
     * @return u such that {@code x=pq(u|m)}
     * @since 2.1
     */
    double arcpq(const double& x, const double deltaQP, const double deltaRQ)
    {
        // see equation 19.25.34 in Digital Library of Mathematical Functions
        // https://dlmf.nist.gov/19.25.E34
        const double x2 = x * x;
        const double w = (1 - x2) / deltaQP;
        const double positive = std::sqrt(w) * Carlson_Elliptic_Integral::r_f(x2, 1, 1 + deltaRQ * w);
        return x < 0
            ? 2 * Legendre_Elliptic_Integral::big_k(get_m()) - positive
            : positive;
    }

protected:
    /** Simple constructor.
     * @param m parameter of the function
     */
    Jacobi_Elliptic(const double m) : my_m{ m } {};

public:
    /** Get the parameter of the function.
     * @return parameter of the function
     */
    double get_m() const
    {
        return my_m;
    }

    /** Evaluate the three principal Jacobi elliptic functions with pole at point n in Glaisher’s Notation.
     * @param u argument of the functions
     * @return copolar trio containing the three principal Jacobi
     * elliptic functions {@code sn(u|m)}, {@code cn(u|m)}, and {@code dn(u|m)}.
     */
    virtual Copolar_N values_n(const double& u);

    /** Evaluate the three subsidiary Jacobi elliptic functions with pole at point s in Glaisher’s Notation.
     * @param u argument of the functions
     * @return copolar trio containing the three subsidiary Jacobi
     * elliptic functions {@code cs(u|m)}, {@code ds(u|m)} and {@code ns(u|m)}.
     */
    Copolar_S values_s(const double& u) 
    {
        return Copolar_S(values_n(u));
    }

    /** Evaluate the three subsidiary Jacobi elliptic functions with pole at point c in Glaisher’s Notation.
     * @param u argument of the functions
     * @return copolar trio containing the three subsidiary Jacobi
     * elliptic functions {@code dc(u|m)}, {@code nc(u|m)}, and {@code sc(u|m)}.
     */
    Copolar_C values_c(const double& u) 
    {
        return Copolar_C(values_n(u));
    }

    /** Evaluate the three subsidiary Jacobi elliptic functions with pole at point d in Glaisher’s Notation.
     * @param u argument of the functions
     * @return copolar trio containing the three subsidiary Jacobi
     * elliptic functions {@code nd(u|m)}, {@code sd(u|m)}, and {@code cd(u|m)}.
     */
    Copolar_D values_d(const double& u) 
    {
        return Copolar_D(values_n(u));
    }

    /** Evaluate inverse of Jacobi elliptic function sn.
     * @param x value of Jacobi elliptic function {@code sn(u|m)}
     * @return u such that {@code x=sn(u|m)}
     * @since 2.1
     */
    double arcsn(const double& x) 
    {
        // p = n, q = c, r = d, see DLMF 19.25.29 for evaluating Δ⁡(q, p) and Δ⁡(r, p)
        return arcsp(x, -1, -get_m());
    }

    /** Evaluate inverse of Jacobi elliptic function cn.
     * @param x value of Jacobi elliptic function {@code cn(u|m)}
     * @return u such that {@code x=cn(u|m)}
     * @since 2.1
     */
    double arccn(const double& x) 
    {
        // p = c, q = n, r = d, see DLMF 19.25.29 for evaluating Δ⁡(q, p) and Δ⁡(r, q)
        return arcpq(x, 1, -get_m());
    }

    /** Evaluate inverse of Jacobi elliptic function dn.
     * @param x value of Jacobi elliptic function {@code dn(u|m)}
     * @return u such that {@code x=dn(u|m)}
     * @since 2.1
     */
    double arcdn(const double& x) 
    {
        // p = d, q = n, r = c, see DLMF 19.25.29 for evaluating Δ⁡(q, p) and Δ⁡(r, q)
        return arcpq(x, get_m(), -1);
    }

    /** Evaluate inverse of Jacobi elliptic function cs.
     * @param x value of Jacobi elliptic function {@code cs(u|m)}
     * @return u such that {@code x=cs(u|m)}
     * @since 2.1
     */
    double arccs(const double& x) 
    {
        // p = c, q = n, r = d, see DLMF 19.25.29 for evaluating Δ⁡(q, p) and Δ⁡(r, p)
        return arcps(x, 1, 1 - get_m());
    }

    /** Evaluate inverse of Jacobi elliptic function ds.
     * @param x value of Jacobi elliptic function {@code ds(u|m)}
     * @return u such that {@code x=ds(u|m)}
     * @since 2.1
     */
    double arcds(const double& x) 
    {
        // p = d, q = c, r = n, see DLMF 19.25.29 for evaluating Δ⁡(q, p) and Δ⁡(r, p)
        return arcps(x, get_m() - 1, get_m());
    }

    /** Evaluate inverse of Jacobi elliptic function ns.
     * @param x value of Jacobi elliptic function {@code ns(u|m)}
     * @return u such that {@code x=ns(u|m)}
     * @since 2.1
     */
    double arcns(const double& x) 
    {
        // p = n, q = c, r = d, see DLMF 19.25.29 for evaluating Δ⁡(q, p) and Δ⁡(r, p)
        return arcps(x, -1, -get_m());
    }

    /** Evaluate inverse of Jacobi elliptic function dc.
     * @param x value of Jacobi elliptic function {@code dc(u|m)}
     * @return u such that {@code x=dc(u|m)}
     * @since 2.1
     */
    double arcdc(const double& x) 
    {
        // p = d, q = c, r = n, see DLMF 19.25.29 for evaluating Δ⁡(q, p) and Δ⁡(r, q)
        return arcpq(x, get_m() - 1, 1);
    }

    /** Evaluate inverse of Jacobi elliptic function nc.
     * @param x value of Jacobi elliptic function {@code nc(u|m)}
     * @return u such that {@code x=nc(u|m)}
     * @since 2.1
     */
    double arcnc(const double& x) 
    {
        // p = n, q = c, r = d, see DLMF 19.25.29 for evaluating Δ⁡(q, p) and Δ⁡(r, q)
        return arcpq(x, -1, 1 - get_m());
    }

    /** Evaluate inverse of Jacobi elliptic function sc.
     * @param x value of Jacobi elliptic function {@code sc(u|m)}
     * @return u such that {@code x=sc(u|m)}
     * @since 2.1
     */
    double arcsc(const double& x) 
    {
        // p = c, q = n, r = d, see DLMF 19.25.29 for evaluating Δ⁡(q, p) and Δ⁡(r, p)
        return arcsp(x, 1, 1 - get_m());
    }

    /** Evaluate inverse of Jacobi elliptic function nd.
     * @param x value of Jacobi elliptic function {@code nd(u|m)}
     * @return u such that {@code x=nd(u|m)}
     * @since 2.1
     */
    double arcnd(const double& x) 
    {
        // p = n, q = d, r = c, see DLMF 19.25.29 for evaluating Δ⁡(q, p) and Δ⁡(r, q)
        return arcpq(x, -get_m(), get_m() - 1);
    }

    /** Evaluate inverse of Jacobi elliptic function sd.
     * @param x value of Jacobi elliptic function {@code sd(u|m)}
     * @return u such that {@code x=sd(u|m)}
     * @since 2.1
     */
    double arcsd(const double& x) 
    {
        // p = d, q = n, r = c, see DLMF 19.25.29 for evaluating Δ⁡(q, p) and Δ⁡(r, p)
        return arcsp(x, get_m(), get_m() - 1);
    }

    /** Evaluate inverse of Jacobi elliptic function cd.
     * @param x value of Jacobi elliptic function {@code cd(u|m)}
     * @return u such that {@code x=cd(u|m)}
     * @since 2.1
     */
    double arccd(const double& x) 
    {
        // p = c, q = d, r = n, see DLMF 19.25.29 for evaluating Δ⁡(q, p) and Δ⁡(r, q)
        return arcpq(x, 1 - get_m(), get_m());
    }
};