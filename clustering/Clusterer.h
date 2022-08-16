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

#include "Clusterable.h"

  //import org.hipparchus.clustering.distance.Distance_Measure;
  //import org.hipparchus.exception.;
  //import org.hipparchus.exception.Math_Illegal_State_Exception;

  /**
   * Base class for clustering algorithms.
   *
   * @param <T> the type of points that can be clustered
   */
template<typename T, typename std::enable_if<std::is_base_of<Clusterable, T>::value>::type* = nullptr>
class Clusterer
{
private:
	/** The distance measure to use. */
	Distance_Measure my_measure;

protected:
	/**
	 * Build a clusterer with the given {@link Distance_Measure}.
	 *
	 * @param measure the distance measure to use
	 */
	Clusterer(const Distance_Measure& measure)
	{
		my_measure = measure;
	}

	/**
	 * Calculates the distance between two {@link Clusterable} instances
	 * with the configured {@link Distance_Measure}.
	 *
	 * @param p1 the first clusterable
	 * @param p2 the second clusterable
	 * @return the distance between the two clusterables
	 */
	protected double distance(const Clusterable p1, const Clusterable p2)
	{
		return measure.compute(p1.get_point(), p2.get_point());
	}

public:
	/**
	 * Perform a cluster analysis on the given set of {@link Clusterable} instances.
	 *
	 * @param points the set of {@link Clusterable} instances
	 * @return a {@link List} of clusters
	 * @ if points are NULL or the number of
	 *   data points is not compatible with this clusterer
	 * @Math_Illegal_State_Exception if the algorithm has not yet converged after
	 *   the maximum number of iterations has been exceeded
	 */
	virtual List< ? extends Cluster<T>> cluster(Collection<T> points) = 0;

	/**
	 * Returns the {@link Distance_Measure} instance used by this clusterer.
	 *
	 * @return the distance measure
	 */
	Distance_Measure get_distance_measure() const
	{
		return my_measure;
	}
};