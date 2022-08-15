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
//package org.hipparchus.geometry.spherical.twod;

//import org.hipparchus.geometry.euclidean.threed.Vector_3D;
//import org.hipparchus.geometry.partitioning.Abstract_Sub_Hyperplane;
//import org.hipparchus.geometry.partitioning.Hyperplane;
//import org.hipparchus.geometry.partitioning.Region;
//import org.hipparchus.geometry.spherical.oned.Arc;
//import org.hipparchus.geometry.spherical.oned.Arcs_Set;
//import org.hipparchus.geometry.spherical.oned.Sphere_1D;
//import org.hipparchus.util.FastMath;

/** This class represents a sub-hyperplane for {@link Circle}.
 */
class Sub_Circle extends Abstract_Sub_Hyperplane<Sphere_2D, Sphere_1D> 
{

    /** Simple constructor.
     * @param hyperplane underlying hyperplane
     * @param remaining_region remaining region of the hyperplane
     */
    public Sub_Circle(const Hyperplane<Sphere_2D>& hyperplane, const Region<Sphere_1D>& remaining_region) 
    {
        super(hyperplane, remaining_region);
    }

    /** {@inherit_doc} */
    //override
    protected Abstract_Sub_Hyperplane<Sphere_2D, Sphere_1D> build_new(const Hyperplane<Sphere_2D>& hyperplane, const Region<Sphere_1D>& remaining_region) 
    {
        return Sub_Circle(hyperplane, remaining_region);
    }

    /** {@inherit_doc} */
    //override
    public Split_Sub_Hyperplane<Sphere_2D> split(const Hyperplane<Sphere_2D>& hyperplane) 
    {

        const Circle this_circle   = (Circle) get_hyperplane();
        const Circle other_circle  = (Circle) hyperplane;
        const double& angle = Vector_3D.angle(this_circle.get_pole(), other_circle.get_pole());

        if (angle < this_circle.get_tolerance() || angle > std::numbers::pi - this_circle.get_tolerance()) 
        {
            // the two circles are aligned or opposite
            return Split_Sub_Hyperplane<Sphere_2D>(null, NULL);
        }
else 
        {
            // the two circles intersect each other
            const Arc    arc          = this_circle.get_inside_arc(other_circle);
            const Arcs_Set.Split split = ((Arcs_Set) get_remaining_region()).split(arc);
            const Arcs_Set plus        = split.get_plus();
            const Arcs_Set minus       = split.get_minus();
            return Split_Sub_Hyperplane<Sphere_2D>(plus  == NULL ? NULL : Sub_Circle(this_circle.copy_self(), plus), minus == NULL ? NULL : Sub_Circle(this_circle.copy_self(), minus));
        }

    }

}


