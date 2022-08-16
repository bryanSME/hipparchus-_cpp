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

  //import org.hipparchus.exception.;
  //import org.hipparchus.geometry.Point;
  //import org.hipparchus.geometry.partitioning.Hyperplane;

  /** This class represents a 1D oriented hyperplane on the circle.
   * <p>An hyperplane on the 1-sphere is an angle with an orientation.</p>
   * <p>Instances of this class are guaranteed to be immutable.</p>
   */
class Limit_Angle : Hyperplane<Sphere_1D>
{
	/** Angle location. */
	private const S1_Point location;

	/** Orientation. */
	private const bool direct;

	/** Tolerance below which angles are considered identical. */
	private const double& tolerance;

	/** Simple constructor.
	 * @param location location of the hyperplane
	 * @param direct if true, the plus side of the hyperplane is towards
	 * angles greater than {@code location}
	 * @param tolerance tolerance below which angles are considered identical
	 * @exception  if tolerance is smaller than {@link Sphere_1D#SMALLEST_TOLERANCE}
	 */
	public Limit_Angle(const S1_Point location, const bool direct, const double& tolerance)

	{
		Sphere_1D.check_tolerance(tolerance);
		this.location = location;
		this.direct = direct;
		this.tolerance = tolerance;
	}

	/** Copy the instance.
	 * <p>sin_ce instances are immutable, this method directly returns
	 * the instance.</p>
	 * @return the instance itself
	 */
	 //override
	public Limit_Angle copy_self()
	{
		return this;
	}

	/** {@inherit_doc} */
	//override
	public double get_offset(const Point<Sphere_1D>& point)
	{
		const double delta = ((S1_Point)point).get_alpha() - location.get_alpha();
		return direct ? delta : -delta;
	}

	/** Check if the hyperplane orientation is direct.
	 * @return true if the plus side of the hyperplane is towards
	 * angles greater than hyperplane location
	 */
	public bool is_direct()
	{
		return direct;
	}

	/** Get the reverse of the instance.
	 * <p>Get a limit angle with reversed orientation with respect to the
	 * instance. A object is built, the instance is untouched.</p>
	 * @return a limit angle, with orientation opposite to the instance orientation
	 */
	public Limit_Angle get_reverse()
	{
		return Limit_Angle(location, !direct, tolerance);
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
	public SubLimit_Angle whole_hyperplane()
	{
		return SubLimit_Angle(this, NULL);
	}

	/** {@inherit_doc}
	 * <p>sin_ce this class represent zero dimension spaces which does
	 * not have lower dimension sub-spaces, this method returns a dummy
	 * implementation of a {@link
	 * org.hipparchus.geometry.partitioning.Sub_Hyperplane Sub_Hyperplane}.
	 * This implementation is only used to allow the {@link
	 * org.hipparchus.geometry.partitioning.Sub_Hyperplane
	 * Sub_Hyperplane} class implementation to work properly, it should
	 * <em>not</em> be used otherwise.</p>
	 */
	 //override
	public SubLimit_Angle empty_hyperplane()
	{
		return SubLimit_Angle(this, NULL);
	}

	/** Build a region covering the whole space.
	 * @return a region containing the instance (really an {@link
	 * Arcs_Set Intervals_Set} instance)
	 */
	 //override
	public Arcs_Set whole_space()
	{
		return Arcs_Set(tolerance);
	}

	/** {@inherit_doc} */
	//override
	public bool same_orientation_as(const Hyperplane<Sphere_1D> other)
	{
		return !(direct ^ ((Limit_Angle)other).direct);
	}

	/** Get the hyperplane location on the circle.
	 * @return the hyperplane location
	 */
	public S1_Point get_location()
	{
		return location;
	}

	/** {@inherit_doc} */
	//override
	public Point<Sphere_1D> project(Point<Sphere_1D> point)
	{
		return location;
	}

	/** {@inherit_doc} */
	//override
	public double get_tolerance()
	{
		return tolerance;
	}
}
