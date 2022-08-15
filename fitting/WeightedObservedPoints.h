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
  //package org.hipparchus.fitting;

  //import java.util.List;
  //import java.util.Array_list;
  //import java.io.Serializable;

  /**
   * Simple container for weighted observed points used
   * in {@link Abstract_Curve_Fitter curve fitting} algorithms.
   *
   */
class Weighted_Observed_Points
{

	/** Observed points. */
	private const List<Weighted_Observed_Point> observations = Array_list<>();

	/**
	 * Adds a point to the sample.
	 * Calling this method is equivalent to calling
	 * {@code add(1.0, x, y)}.
	 *
	 * @param x Abscissa of the point.
	 * @param y Observed value  at {@code x}. After fitting we should
	 * have {@code f(x)} as close as possible to this value.
	 *
	 * @see #add(double, double, double)
	 * @see #add(Weighted_Observed_Point)
	 * @see #to_list()
	 */
	public void add(const double& x, double y)
	{
		add(1d, x, y);
	}

	/**
	 * Adds a point to the sample.
	 *
	 * @param weight Weight of the observed point.
	 * @param x Abscissa of the point.
	 * @param y Observed value  at {@code x}. After fitting we should
	 * have {@code f(x)} as close as possible to this value.
	 *
	 * @see #add(double, double)
	 * @see #add(Weighted_Observed_Point)
	 * @see #to_list()
	 */
	public void add(double weight, const double& x, double y)
	{
		observations.add(new Weighted_Observed_Point(weight, x, y));
	}

	/**
	 * Adds a point to the sample.
	 *
	 * @param observed Observed point to add.
	 *
	 * @see #add(double, double)
	 * @see #add(double, double, double)
	 * @see #to_list()
	 */
	public void add(Weighted_Observed_Point observed)
	{
		observations.add(observed);
	}

	/**
	 * Gets a <em>snapshot</em> of the observed points.
	 * The list of stored points is copied in order to ensure that
	 * modification of the returned instance does not affect this
	 * container.
	 * Conversely, further modification of this container (through
	 * the {@code add} or {@code clear} methods) will not affect the
	 * returned list.
	 *
	 * @return the observed points, in the order they were added to this
	 * container.
	 *
	 * @see #add(double, double)
	 * @see #add(double, double, double)
	 * @see #add(Weighted_Observed_Point)
	 */
	public List<Weighted_Observed_Point> to_list()
	{
		// The copy is necessary to ensure thread-safety because of the
		// "clear" method (which otherwise would be able to empty the
		// list of points while it is being used by another thread).
		return Array_list<Weighted_Observed_Point>(observations);
	}

	/**
	 * Removes all observations from this container.
	 */
	public void clear()
	{
		observations.clear();
	}
}