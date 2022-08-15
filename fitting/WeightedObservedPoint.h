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

  //import java.io.Serializable;

  /**
   * This class is a simple container for weighted observed point in
   * {@link Abstract_Curve_Fitter curve fitting}.
   * <p>Instances of this class are guaranteed to be immutable.</p>
   */
class Weighted_Observed_Point
{
	/** Weight of the measurement in the fitting process. */
	private const double weight;
	/** Abscissa of the point. */
	private const double x;
	/** Observed value of the function at x. */
	private const double y;

	/**
	 * Simple constructor.
	 *
	 * @param weight Weight of the measurement in the fitting process.
	 * @param x Abscissa of the measurement.
	 * @param y Ordinate of the measurement.
	 */
	public Weighted_Observed_Point(const double weight, const double& x, const double y)
	{
		this.weight = weight;
		this.x = x;
		this.y = y;
	}

	/**
	 * Gets the weight of the measurement in the fitting process.
	 *
	 * @return the weight of the measurement in the fitting process.
	 */
	public double get_weight()
	{
		return weight;
	}

	/**
	 * Gets the abscissa of the point.
	 *
	 * @return the abscissa of the point.
	 */
	public double get_x()
	{
		return x;
	}

	/**
	 * Gets the observed value of the function at x.
	 *
	 * @return the observed value of the function at x.
	 */
	public double get_y()
	{
		return y;
	}
}