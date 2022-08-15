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

#include <type_traits>
#include <vector>
#include "Clusterable.h"

  /**
   * Cluster holding a set of {@link Clusterable} points.
   * @param <T> the type of points that can be clustered
   */
template<typename T, typename std::enable_if<std::is_base_of<Clusterable, T>::value>::type* = nullptr>
class Cluster
{
private:
	/** The points contained in this cluster. */
	const std::vector<T> points;

private:
	/**
	 * Build a cluster centered at a specified point.
	 */
	Cluster() = default;

	/**
	 * Add a point to this cluster.
	 * @param point point to add
	 */
	void add_point(const T& point)
	{
		my_points.add(point);
	}

	/**
	 * Get the points contained in the cluster.
	 * @return points contained in the cluster
	 */
	public List<T> get_points()
	{
		return points;
	}
};