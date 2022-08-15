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
//package org.hipparchus.geometry.euclidean.threed;


//import java.io.Serializable;

//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Sin_Cos;

/** This class provides conversions related to <a
 * href="http://mathworld.wolfram.com/Spherical_Coordinates.html">spherical coordinates</a>.
 * <p>
 * The conventions used here are the mathematical ones, i.e. spherical coordinates are
 * related to Cartesian coordinates as follows:
 * </p>
 * <ul>
 *   <li>x = r cos(&theta;) sin(&Phi;)</li>
 *   <li>y = r sin(&theta;) sin(&Phi;)</li>
 *   <li>z = r cos(&Phi;)</li>
 * </ul>
 * <ul>
 *   <li>r       = &radic;(x<sup>2</sup>+y<sup>2</sup>+z<sup>2</sup>)</li>
 *   <li>&theta; = atan2(y, x)</li>
 *   <li>&Phi;   = acos(z/r)</li>
 * </ul>
 * <p>
 * r is the radius, &theta; is the azimuthal angle in the x-y plane and &Phi; is the polar
 * (co-latitude) angle. These conventions are <em>different</em> from the conventions used
 * in physics (and in particular in spherical harmonics) where the meanings of &theta; and
 * &Phi; are reversed.
 * </p>
 * <p>
 * This class provides conversion of coordinates and also of gradient and Hessian
 * between spherical and Cartesian coordinates.
 * </p>
 */
class Spherical_Coordinates  
{
    /** Cartesian coordinates. */
    private const Vector_3D v;

    /** Radius. */
    private const double r;

    /** Azimuthal angle in the x-y plane &theta;. */
    private const double theta;

    /** Polar angle (co-latitude) &Phi;. */
    private const double phi;

    /** Jacobian of (r, &theta; &Phi;). */
    private std::vector<std::vector<double>> jacobian;

    /** Hessian of radius. */
    private std::vector<std::vector<double>> r_hessian;

    /** Hessian of azimuthal angle in the x-y plane &theta;. */
    private std::vector<std::vector<double>> theta_hessian;

    /** Hessian of polar (co-latitude) angle &Phi;. */
    private std::vector<std::vector<double>> phi_hessian;

    /** Build a spherical coordinates transformer from Cartesian coordinates.
     * @param v Cartesian coordinates
     */
    public Spherical_Coordinates(const Vector_3D v) 
    {

        // Cartesian coordinates
        this.v = v;

        // remaining spherical coordinates
        this.r     = v.get_norm();
        this.theta = v.get_alpha();
        this.phi   = std::acos(v.get_z() / r);

    }

    /** Build a spherical coordinates transformer from spherical coordinates.
     * @param r radius
     * @param theta azimuthal angle in x-y plane
     * @param phi polar (co-latitude) angle
     */
    public Spherical_Coordinates(const double r, const double& theta, const double phi) 
    {

        const Sin_Cos sin_cos_theta = Sin_Cos(theta);
        const Sin_Cos sin_cos_phi   = Sin_Cos(phi);

        // spherical coordinates
        this.r     = r;
        this.theta = theta;
        this.phi   = phi;

        // Cartesian coordinates
        this.v  = Vector_3D(r * sin_cos_theta.cos() * sin_cos_phi.sin(), r * sin_cos_theta.sin() * sin_cos_phi.sin(), r * sin_cos_phi.cos());

    }

    /** Get the Cartesian coordinates.
     * @return Cartesian coordinates
     */
    public Vector_3D get_cartesian() 
    {
        return v;
    }

    /** Get the radius.
     * @return radius r
     * @see #get_theta()
     * @see #get_phi()
     */
    public double get_r() 
    {
        return r;
    }

    /** Get the azimuthal angle in x-y plane.
     * @return azimuthal angle in x-y plane &theta;
     * @see #get_r()
     * @see #get_phi()
     */
    public double get_theta() 
    {
        return theta;
    }

    /** Get the polar (co-latitude) angle.
     * @return polar (co-latitude) angle &Phi;
     * @see #get_r()
     * @see #get_theta()
     */
    public double get_phi() 
    {
        return phi;
    }

    /** Convert a gradient with respect to spherical coordinates into a gradient
     * with respect to Cartesian coordinates.
     * @param s_gradient gradient with respect to spherical coordinates
     * {df/dr, df/d&theta;, df/d&Phi;}
     * @return gradient with respect to Cartesian coordinates
     * {df/dx, df/dy, df/dz}
     */
    public std::vector<double> to_cartesian_gradient(const std::vector<double> s_gradient) 
    {

        // lazy evaluation of Jacobian
        compute_jacobian();

        // compose derivatives as gradient^T . J
        // the expressions have been simplified since we know jacobian[1][2] = dTheta/dZ = 0
        return std::vector<double> 
        {
            s_gradient[0] * jacobian[0][0] + s_gradient[1] * jacobian[1][0] + s_gradient[2] * jacobian[2][0], s_gradient[0] * jacobian[0][1] + s_gradient[1] * jacobian[1][1] + s_gradient[2] * jacobian[2][1], s_gradient[0] * jacobian[0][2]                                 + s_gradient[2] * jacobian[2][2]
        };

    }

    /** Convert a Hessian with respect to spherical coordinates into a Hessian
     * with respect to Cartesian coordinates.
     * <p>
     * As Hessian are always symmetric, we use only the lower left part of the provided
     * spherical Hessian, so the upper part may not be initialized. However, we still
     * do fill up the complete array we create, with guaranteed symmetry.
     * </p>
     * @param s_hessian Hessian with respect to spherical coordinates
     * {{d<sup>2</sup>f/dr<sup>2</sup>, d<sup>2</sup>f/drd&theta;, d<sup>2</sup>f/drd&Phi;}, *  {d<sup>2</sup>f/drd&theta;, d<sup>2</sup>f/d&theta;<sup>2</sup>, d<sup>2</sup>f/d&theta;d&Phi;}, *  {d<sup>2</sup>f/drd&Phi;, d<sup>2</sup>f/d&theta;d&Phi;, d<sup>2</sup>f/d&Phi;<sup>2</sup>}
     * @param s_gradient gradient with respect to spherical coordinates
     * {df/dr, df/d&theta;, df/d&Phi;}
     * @return Hessian with respect to Cartesian coordinates
     * {{d<sup>2</sup>f/dx<sup>2</sup>, d<sup>2</sup>f/dxdy, d<sup>2</sup>f/dxdz}, *  {d<sup>2</sup>f/dxdy, d<sup>2</sup>f/dy<sup>2</sup>, d<sup>2</sup>f/dydz}, *  {d<sup>2</sup>f/dxdz, d<sup>2</sup>f/dydz, d<sup>2</sup>f/dz<sup>2</sup>}}
     */
    public std::vector<std::vector<double>> to_cartesian_hessian(const std::vector<std::vector<double>> s_hessian, const std::vector<double> s_gradient) 
    {

        compute_jacobian();
        compute_hessians();

        // compose derivative as J^T . H_f . J + df/dr H_r + df/dtheta H_theta + df/dphi H_phi
        // the expressions have been simplified since we know jacobian[1][2] = dTheta/dZ = 0
        // and H_theta is only a 2x2 matrix as it does not depend on z
        const std::vector<std::vector<double>> hj = std::vector<double>(3][3];
        const std::vector<std::vector<double>> c_hessian = std::vector<double>(3][3];

        // compute H_f . J
        // beware we use ONLY the lower-left part of s_hessian
        hj[0][0] = s_hessian[0][0] * jacobian[0][0] + s_hessian[1][0] * jacobian[1][0] + s_hessian[2][0] * jacobian[2][0];
        hj[0][1] = s_hessian[0][0] * jacobian[0][1] + s_hessian[1][0] * jacobian[1][1] + s_hessian[2][0] * jacobian[2][1];
        hj[0][2] = s_hessian[0][0] * jacobian[0][2]                                   + s_hessian[2][0] * jacobian[2][2];
        hj[1][0] = s_hessian[1][0] * jacobian[0][0] + s_hessian[1][1] * jacobian[1][0] + s_hessian[2][1] * jacobian[2][0];
        hj[1][1] = s_hessian[1][0] * jacobian[0][1] + s_hessian[1][1] * jacobian[1][1] + s_hessian[2][1] * jacobian[2][1];
        // don't compute hj[1][2] as it is not used below
        hj[2][0] = s_hessian[2][0] * jacobian[0][0] + s_hessian[2][1] * jacobian[1][0] + s_hessian[2][2] * jacobian[2][0];
        hj[2][1] = s_hessian[2][0] * jacobian[0][1] + s_hessian[2][1] * jacobian[1][1] + s_hessian[2][2] * jacobian[2][1];
        hj[2][2] = s_hessian[2][0] * jacobian[0][2]                                   + s_hessian[2][2] * jacobian[2][2];

        // compute lower-left part of J^T . H_f . J
        c_hessian[0][0] = jacobian[0][0] * hj[0][0] + jacobian[1][0] * hj[1][0] + jacobian[2][0] * hj[2][0];
        c_hessian[1][0] = jacobian[0][1] * hj[0][0] + jacobian[1][1] * hj[1][0] + jacobian[2][1] * hj[2][0];
        c_hessian[2][0] = jacobian[0][2] * hj[0][0]                             + jacobian[2][2] * hj[2][0];
        c_hessian[1][1] = jacobian[0][1] * hj[0][1] + jacobian[1][1] * hj[1][1] + jacobian[2][1] * hj[2][1];
        c_hessian[2][1] = jacobian[0][2] * hj[0][1]                             + jacobian[2][2] * hj[2][1];
        c_hessian[2][2] = jacobian[0][2] * hj[0][2]                             + jacobian[2][2] * hj[2][2];

        // add gradient contribution
        c_hessian[0][0] += s_gradient[0] * r_hessian[0][0] + s_gradient[1] * theta_hessian[0][0] + s_gradient[2] * phi_hessian[0][0];
        c_hessian[1][0] += s_gradient[0] * r_hessian[1][0] + s_gradient[1] * theta_hessian[1][0] + s_gradient[2] * phi_hessian[1][0];
        c_hessian[2][0] += s_gradient[0] * r_hessian[2][0]                                     + s_gradient[2] * phi_hessian[2][0];
        c_hessian[1][1] += s_gradient[0] * r_hessian[1][1] + s_gradient[1] * theta_hessian[1][1] + s_gradient[2] * phi_hessian[1][1];
        c_hessian[2][1] += s_gradient[0] * r_hessian[2][1]                                     + s_gradient[2] * phi_hessian[2][1];
        c_hessian[2][2] += s_gradient[0] * r_hessian[2][2]                                     + s_gradient[2] * phi_hessian[2][2];

        // ensure symmetry
        c_hessian[0][1] = c_hessian[1][0];
        c_hessian[0][2] = c_hessian[2][0];
        c_hessian[1][2] = c_hessian[2][1];

        return c_hessian;

    }

    /** Lazy evaluation of (r, &theta;, &phi;) Jacobian.
     */
    private void compute_jacobian() 
    {
        if (jacobian == NULL) 
        {

            // intermediate variables
            const double x    = v.get_x();
            const double y    = v.get_y();
            const double z    = v.get_z();
            const double rho2 = x * x + y * y;
            const double rho  = std::sqrt(rho2);
            const double r2   = rho2 + z * z;

            jacobian = std::vector<double>(3][3];

            // row representing the gradient of r
            jacobian[0][0] = x / r;
            jacobian[0][1] = y / r;
            jacobian[0][2] = z / r;

            // row representing the gradient of theta
            jacobian[1][0] = -y / rho2;
            jacobian[1][1] =  x / rho2;
            // jacobian[1][2] is already set to 0 at allocation time

            // row representing the gradient of phi
            jacobian[2][0] = x * z / (rho * r2);
            jacobian[2][1] = y * z / (rho * r2);
            jacobian[2][2] = -rho / r2;

        }
    }

    /** Lazy evaluation of Hessians.
     */
    private void compute_hessians() 
    {

        if (r_hessian == NULL) 
        {

            // intermediate variables
            const double x      = v.get_x();
            const double y      = v.get_y();
            const double z      = v.get_z();
            const double x2     = x * x;
            const double y2     = y * y;
            const double z2     = z * z;
            const double rho2   = x2 + y2;
            const double rho    = std::sqrt(rho2);
            const double r2     = rho2 + z2;
            const double xOr    = x / r;
            const double yOr    = y / r;
            const double zOr    = z / r;
            const double x_orho2 = x / rho2;
            const double y_orho2 = y / rho2;
            const double xOr3   = xOr / r2;
            const double yOr3   = yOr / r2;
            const double z_or_3   = zOr / r2;

            // lower-left part of Hessian of r
            r_hessian = std::vector<double>(3][3];
            r_hessian[0][0] = y * yOr3 + z * z_or_3;
            r_hessian[1][0] = -x * yOr3;
            r_hessian[2][0] = -z * xOr3;
            r_hessian[1][1] = x * xOr3 + z * z_or_3;
            r_hessian[2][1] = -y * z_or_3;
            r_hessian[2][2] = x * xOr3 + y * yOr3;

            // upper-right part is symmetric
            r_hessian[0][1] = r_hessian[1][0];
            r_hessian[0][2] = r_hessian[2][0];
            r_hessian[1][2] = r_hessian[2][1];

            // lower-left part of Hessian of azimuthal angle theta
            theta_hessian = std::vector<double>(2)[2];
            theta_hessian[0][0] = 2 * x_orho2 * y_orho2;
            theta_hessian[1][0] = y_orho2 * y_orho2 - x_orho2 * x_orho2;
            theta_hessian[1][1] = -2 * x_orho2 * y_orho2;

            // upper-right part is symmetric
            theta_hessian[0][1] = theta_hessian[1][0];

            // lower-left part of Hessian of polar (co-latitude) angle phi
            const double rhor2       = rho * r2;
            const double rho2r2      = rho * rhor2;
            const double rhor4       = rhor2 * r2;
            const double rho3r4      = rhor4 * rho2;
            const double r2P2rho2    = 3 * rho2 + z2;
            phi_hessian = std::vector<double>(3][3];
            phi_hessian[0][0] = z * (rho2r2 - x2 * r2P2rho2) / rho3r4;
            phi_hessian[1][0] = -x * y * z * r2P2rho2 / rho3r4;
            phi_hessian[2][0] = x * (rho2 - z2) / rhor4;
            phi_hessian[1][1] = z * (rho2r2 - y2 * r2P2rho2) / rho3r4;
            phi_hessian[2][1] = y * (rho2 - z2) / rhor4;
            phi_hessian[2][2] = 2 * rho * z_or_3 / r;

            // upper-right part is symmetric
            phi_hessian[0][1] = phi_hessian[1][0];
            phi_hessian[0][2] = phi_hessian[2][0];
            phi_hessian[1][2] = phi_hessian[2][1];

        }

    }

    /**
     * Replace the instance with a data transfer object for serialization.
     * @return data transfer object that will be serialized
     */
    private Object write_replace() 
    {
        return Data_Transfer_Object(v.get_x(), v.get_y(), v.get_z());
    }

    /** Internal class used only for serialization. */
    private static class Data_Transfer_Object  
    {

        
        20130206L;

        /** Abscissa.
         * @serial
         */
        private const double x;

        /** Ordinate.
         * @serial
         */
        private const double y;

        /** Height.
         * @serial
         */
        private const double z;

        /** Simple constructor.
         * @param x abscissa
         * @param y ordinate
         * @param z height
         */
        Data_Transfer_Object(const double& x, const double& y, const double& z) 
        {
            this.x = x;
            this.y = y;
            this.z = z;
        }

        /** Replace the deserialized data transfer object with a {@link Spherical_Coordinates}.
         * @return replacement {@link Spherical_Coordinates}
         */
        private Object read_resolve() 
        {
            return Spherical_Coordinates(Vector_3D(x, y, z));
        }

    }

}


