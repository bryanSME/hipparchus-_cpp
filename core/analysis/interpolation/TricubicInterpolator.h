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
//package org.hipparchus.analysis.interpolation;

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.util.Math_Arrays;
//import org.hipparchus.util.Math_Utils;

/**
 * Generates a tricubic interpolating function.
 *
 */
class Tricubic_Interpolator
    : Trivariate_Grid_Interpolator 
    {
    /**
     * {@inherit_doc}
     */
    //override
    public Tricubic_Interpolating_Function interpolate(const std::vector<double>& xval, const std::vector<double>& yval, const std::vector<double> zval, const std::vector<std::vector<double>>[] fval)
         
        {
        if (xval.size() == 0 || yval.size() == 0 || zval.size() == 0 || fval.size() == 0) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::NO_DATA);
        }
        if (xval.size() != fval.size()) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, xval.size(), fval.size());
        }

        Math_Arrays::check_order(xval);
        Math_Arrays::check_order(yval);
        Math_Arrays::check_order(zval);

        const int x_len = xval.size();
        const int y_len = yval.size();
        const int z_len = zval.size();

        // Approximation to the partial derivatives using finite differences.
        const std::vector<std::vector<double>>[] dFdX = std::vector<double>(x_len][y_len][z_len];
        const std::vector<std::vector<double>>[] d_fd_y = std::vector<double>(x_len][y_len][z_len];
        const std::vector<std::vector<double>>[] d_fd_z = std::vector<double>(x_len][y_len][z_len];
        const std::vector<std::vector<double>>[] d2FdXdY = std::vector<double>(x_len][y_len][z_len];
        const std::vector<std::vector<double>>[] d2_fd_xd_z = std::vector<double>(x_len][y_len][z_len];
        const std::vector<std::vector<double>>[] d2_fd_yd_z = std::vector<double>(x_len][y_len][z_len];
        const std::vector<std::vector<double>>[] d3_fd_xd_yd_z = std::vector<double>(x_len][y_len][z_len];

        for (int i{ 1 }; i < x_len - 1; i++) 
        {
            Math_Utils::check_dimension(yval.size(), fval[i].size());

            const int& nI = i + 1;
            const int pI = i - 1;

            const double nX = xval[nI];
            const double pX = xval[pI];

            const double delta_x = nX - pX;

            for (int j{ 1 }; j < y_len - 1; j++) 
            {
                Math_Utils::check_dimension(zval.size(), fval[i][j].size());

                const int& nJ = j + 1;
                const int pJ = j - 1;

                const double nY = yval[nJ];
                const double pY = yval[pJ];

                const double delta_y = nY - pY;
                const double delta_x_y = delta_x * delta_y;

                for (int k{ 1 }; k < z_len - 1; k++) 
                {
                    const int& nK = k + 1;
                    const int pK = k - 1;

                    const double nZ = zval[nK];
                    const double pZ = zval[pK];

                    const double delta_z = nZ - pZ;

                    dFdX[i][j][k] = (fval[nI][j][k] - fval[pI][j][k]) / delta_x;
                    d_fd_y[i][j][k] = (fval[i][nJ][k] - fval[i][pJ][k]) / delta_y;
                    d_fd_z[i][j][k] = (fval[i][j][nK] - fval[i][j][pK]) / delta_z;

                    const double delta_x_z = delta_x * delta_z;
                    const double delta_y_z = delta_y * delta_z;

                    d2FdXdY[i][j][k] = (fval[nI][nJ][k] - fval[nI][pJ][k] - fval[pI][nJ][k] + fval[pI][pJ][k]) / delta_x_y;
                    d2_fd_xd_z[i][j][k] = (fval[nI][j][nK] - fval[nI][j][pK] - fval[pI][j][nK] + fval[pI][j][pK]) / delta_x_z;
                    d2_fd_yd_z[i][j][k] = (fval[i][nJ][nK] - fval[i][nJ][pK] - fval[i][pJ][nK] + fval[i][pJ][pK]) / delta_y_z;

                    const double delta_x_y_z = delta_x_y * delta_z;

                    d3_fd_xd_yd_z[i][j][k] = (fval[nI][nJ][nK] - fval[nI][pJ][nK] -
                                          fval[pI][nJ][nK] + fval[pI][pJ][nK] -
                                          fval[nI][nJ][pK] + fval[nI][pJ][pK] +
                                          fval[pI][nJ][pK] - fval[pI][pJ][pK]) / delta_x_y_z;
                }
            }
        }

        // Create the interpolating function.
        return Tricubic_Interpolating_Function(xval, yval, zval, fval, dFdX, d_fd_y, d_fd_z, d2FdXdY, d2_fd_xd_z, d2_fd_yd_z, d3_fd_xd_yd_z) 
        {
            /** {@inherit_doc} */
            //override
            public bool is_valid_point(const double& x, const double& y, const double& z) 
            {
                if (x < xval[1] ||
                    x > xval[xval.size() - 2] ||
                    y < yval[1] ||
                    y > yval[yval.size() - 2] ||
                    z < zval[1] ||
                    z > zval[zval.size() - 2]) 
                    {
                    return false;
                }
else 
                {
                    return true;
                }
            }
        };
    }
}


