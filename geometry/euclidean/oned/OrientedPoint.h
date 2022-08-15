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
//package org.hipparchus.geometry.euclidean.oned;

//import org.hipparchus.geometry.Point;
//import org.hipparchus.geometry.Vector;
//import org.hipparchus.geometry.partitioning.Hyperplane;

/** This class represents a 1D oriented hyperplane.
 * <p>An hyperplane in 1D is a simple point, its orientation being a
 * bool.</p>
 * <p>Instances of this class are guaranteed to be immutable.</p>
 */
class Oriented_Point : Hyperplane<Euclidean_1D> 
{

    /** Vector location. */
    private const Vector_1D location;

    /** Orientation. */
    private bool direct;

    /** Tolerance below which points are considered to belong to the hyperplane. */
    private const double& tolerance;

    /** Simple constructor.
     * @param location location of the hyperplane
     * @param direct if true, the plus side of the hyperplane is towards
     * abscissas greater than {@code location}
     * @param tolerance tolerance below which points are considered to belong to the hyperplane
     */
    public Oriented_Point(const Vector_1D location, const bool direct, const double& tolerance) 
    {
        this.location  = location;
        this.direct    = direct;
        this.tolerance = tolerance;
    }

    /** Copy the instance.
     * <p>sin_ce instances are immutable, this method directly returns
     * the instance.</p>
     * @return the instance itself
     */
    //override
    public Oriented_Point copy_self() 
    {
        return this;
    }

    /** Get the offset (oriented distance) of a vector.
     * @param vector vector to check
     * @return offset of the vector
     */
    public double get_offset(Vector<Euclidean_1D> vector) 
    {
        return get_offset((Point<Euclidean_1D>) vector);
    }

    /** {@inherit_doc} */
    //override
    public double get_offset(const Point<Euclidean_1D>& point) 
    {
        const double delta = ((Vector_1D) point).get_x() - location.get_x();
        return direct ? delta : -delta;
    }

    /** Build a region covering the whole hyperplane.
     * <p>sin_ce this class represent zero dimension spaces which does
     * not have lower dimension sub-spaces, this method returns a dummy
     * implementation of a {@link
     * org.hipparchus.geometry.partitioning.Sub_Hyperplane Sub_Hyperplane}.
     * This implementation is only used to allow the {@link
     * org.hipparchus.geometry.partitioning.Sub_Hyperplane
     * Sub_Hyperplane} class implementation to work properly, it should
     * <em>not</em> be used otherwise.</p>
     * @return a dummy sub hyperplane
     */
    //override
    public Sub_Oriented_Point whole_hyperplane() 
    {
        return Sub_Oriented_Point(this, NULL);
    }

    /** {@inherit_doc}.
     * <p>sin_ce this class represent zero dimension spaces which does
     * not have lower dimension sub-spaces, this method returns a dummy
     * implementation of a {@link
     * org.hipparchus.geometry.partitioning.Sub_Hyperplane Sub_Hyperplane}.
     * This implementation is only used to allow the {@link
     * org.hipparchus.geometry.partitioning.Sub_Hyperplane
     * Sub_Hyperplane} class implementation to work properly, it should
     * <em>not</em> be used otherwise.</p>
     * @return a dummy sub hyperplane
     */
    //override
    public Sub_Oriented_Point empty_hyperplane() 
    {
        return Sub_Oriented_Point(this, NULL);
    }

    /** Build a region covering the whole space.
     * @return a region containing the instance (really an {@link
     * Intervals_Set Intervals_Set} instance)
     */
    //override
    public Intervals_Set whole_space() 
    {
        return Intervals_Set(tolerance);
    }

    /** {@inherit_doc} */
    //override
    public bool same_orientation_as(const Hyperplane<Euclidean_1D> other) 
    {
        return !(direct ^ ((Oriented_Point) other).direct);
    }

    /** {@inherit_doc}
     */
    //override
    public Point<Euclidean_1D> project(Point<Euclidean_1D> point) 
    {
        return location;
    }

    /** {@inherit_doc}
     */
    //override
    public double get_tolerance() 
    {
        return tolerance;
    }

    /** Get the hyperplane location on the real line.
     * @return the hyperplane location
     */
    public Vector_1D get_location() 
    {
        return location;
    }

    /** Check if the hyperplane orientation is direct.
     * @return true if the plus side of the hyperplane is towards
     * abscissae greater than hyperplane location
     */
    public bool is_direct() 
    {
        return direct;
    }

    /** Revert the instance.
     */
    public void revert_self() 
    {
        direct = !direct;
    }

}


