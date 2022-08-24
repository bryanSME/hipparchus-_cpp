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
  //package org.hipparchus.clustering;

  /**
   * A Cluster used by centroid-based clustering algorithms.
   * <p>
   * Defines additionally a cluster center which may not necessarily be a member
   * of the original data set.
   *
   * @param <T> the type of points that can be clustered
   */
template<typename T, typename std::enable_if<std::is_base_of<Clusterable, T>::value>::type* = nullptr>
class Centroid_Cluster : public Cluster<T>
{
private:
	/** Center of the cluster. */
	Clusterable my_center;

public:
	/**
	 * Build a cluster centered at a specified point.
	 * @param center the point which is to be the center of this cluster
	 */
	Centroid_Cluster(const Clusterable& center)
	{
		Cluster<T>();
		my_center = center;
	}

	/**
	 * Get the point chosen to be the center of this cluster.
	 * @return chosen cluster center
	 */
	Clusterable get_center() const
	{
		return my_center;
	}
};