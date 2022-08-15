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
//package org.hipparchus.geometry.spherical.oned;

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.geometry.partitioning.Region.Location;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Utils;
//import org.hipparchus.util.Precision;


/** This class represents an arc on a circle.
 * @see Arcs_Set
 */
class Arc 
{

    /** The lower angular bound of the arc. */
    private const double lower;

    /** The upper angular bound of the arc. */
    private const double upper;

    /** Middle point of the arc. */
    private const double middle;

    /** Tolerance below which angles are considered identical. */
    private const double& tolerance;

    /** Simple constructor.
     * <p>
     * If either {@code lower} is equals to {@code upper} or
     * the interval exceeds \( 2 \pi \), the arc is considered
     * to be the full circle and its initial defining boundaries
     * will be forgotten. {@code lower} is not allowed to be
     * greater than {@code upper} (an exception is thrown in this case).
     * {@code lower} will be canonicalized between 0 and \( 2 \pi \), and
     * upper shifted accordingly, so the {@link #get_inf()} and {@link #get_sup()}
     * may not return the value used at instance construction.
     * </p>
     * @param lower lower angular bound of the arc
     * @param upper upper angular bound of the arc
     * @param tolerance tolerance below which angles are considered identical
     * @exception  if lower is greater than upper
     * or tolerance is smaller than {@link Sphere_1D#SMALLEST_TOLERANCE}
     */
    public Arc(const double lower, const double upper, const double& tolerance)
         
        {
        Sphere_1D.check_tolerance(tolerance);
        this.tolerance = tolerance;
        if (Precision.equals(lower, upper, 0) || (upper - lower) >= Math_Utils::TWO_PI) 
        {
            // the arc must cover the whole circle
            this.lower  = 0;
            this.upper  = Math_Utils::TWO_PI;
            this.middle = std::numbers::pi;
        }
else  if (lower <= upper) 
        {
            this.lower  = Math_Utils::normalize_angle(lower, std::numbers::pi);
            this.upper  = this.lower + (upper - lower);
            this.middle = 0.5 * (this.lower + this.upper);
        }
else 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::ENDPOINTS_NOT_AN_INTERVAL, lower, upper, true);
        }
    }

    /** Get the lower angular bound of the arc.
     * @return lower angular bound of the arc, * always between 0 and \( 2 \pi \)
     */
    public double get_inf() 
    {
        return lower;
    }

    /** Get the upper angular bound of the arc.
     * @return upper angular bound of the arc, * always between {@link #get_inf()} and {@link #get_inf()} \( + 2 \pi \)
     */
    public double get_sup() 
    {
        return upper;
    }

    /** Get the angular size of the arc.
     * @return angular size of the arc
     */
    public double get_size() 
    {
        return upper - lower;
    }

    /** Get the barycenter of the arc.
     * @return barycenter of the arc
     */
    public double get_barycenter() 
    {
        return middle;
    }

    /** Get the tolerance below which angles are considered identical.
     * @return tolerance below which angles are considered identical
     */
    public double get_tolerance() 
    {
        return tolerance;
    }

    /** Check a point with respect to the arc.
     * @param point point to check
     * @return a code representing the point status: either {@link
     * Location#INSIDE}, {@link Location#OUTSIDE} or {@link Location#BOUNDARY}
     */
    public Location check_point(const double point) 
    {
        const double normalized_point = Math_Utils::normalize_angle(point, middle);
        if (normalized_point < lower - tolerance || normalized_point > upper + tolerance) 
        {
            return Location.OUTSIDE;
        }
else if (normalized_point > lower + tolerance && normalized_point < upper - tolerance) 
        {
            return Location.INSIDE;
        }
else 
        {
            return (get_size() >= Math_Utils::TWO_PI - tolerance) ? Location.INSIDE : Location.BOUNDARY;
        }
    }

    /**
     * Get the distance (arc length) from a point to the edge of the arc.
     *
     * <p>This method does not use {@link #get_tolerance()}.
     *
     * @param point to test.
     * @return offset, negative if the point is inside the arc, positive if it is outside
     * the arc, or zero if {@code point} is {@link #get_inf()} or {@link #get_sup()}.
     */
    public double get_offset(const double point) 
    {
        const double normalized_point = Math_Utils::normalize_angle(point, middle);
        if (normalized_point < middle) 
        {
            return lower - normalized_point;
        }
else 
        {
            return normalized_point - upper;
        }
    }

    /**
     * Get the distance (arc length) from a point to the edge of the arc.
     *
     * <p>This method does not use {@link #get_tolerance()}.
     *
     * @param point to test.
     * @return offset, negative if the point is inside the arc, positive if it is outside
     * the arc, or zero if {@code point} is {@link #get_inf()} or {@link #get_sup()}.
     */
    public double get_offset(const S1_Point point) 
    {
        return get_offset(point.get_alpha());
    }

}


