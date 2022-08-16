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

  //package org.hipparchus.clustering.evaluation;

  //import java.util.List;

  //import org.hipparchus.clustering.Centroid_Cluster;
  //import org.hipparchus.clustering.Cluster;
  //import org.hipparchus.clustering.Clusterable;
  //import org.hipparchus.clustering.Double_Point;
  //import org.hipparchus.clustering.distance.Distance_Measure;
  //import org.hipparchus.clustering.distance.Euclidean_Distance;

  /**
   * Base class for cluster evaluation methods.
   *
   * @param <T> type of the clustered points
   */
class Cluster_Evaluator<T extends Clusterable>
{
	/** The distance measure to use when evaluating the cluster. */
	private const Distance_Measure measure;

	/**
	 * Creates a cluster evaluator with an {@link Euclidean_Distance}
	 * as distance measure.
	 */
	public Cluster_Evaluator()
	{
		this(new Euclidean_Distance());
	}

	/**
	 * Creates a cluster evaluator with the given distance measure.
	 * @param measure the distance measure to use
	 */
	public Cluster_Evaluator(const Distance_Measure measure)
	{
		this.measure = measure;
	}

	/**
	 * Computes the evaluation score for the given list of clusters.
	 * @param clusters the clusters to evaluate
	 * @return the computed score
	 */
	public virtual double score(List< ? extends Cluster<T>> clusters);

	/**
	 * Returns whether the first evaluation score is considered to be better
	 * than the second one by this evaluator.
	 * <p>
	 * Specific implementations shall //override this method if the returned scores
	 * do not follow the same ordering, i.e. smaller score is better.
	 *
	 * @param score1 the first score
	 * @param score2 the second score
	 * @return {@code true} if the first score is considered to be better, {@code false} otherwise
	 */
	public bool is_better_score(double score1, double score2)
	{
		return score1 < score2;
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

	/**
	 * Computes the centroid for a cluster.
	 *
	 * @param cluster the cluster
	 * @return the computed centroid for the cluster, * or {@code NULL} if the cluster does not contain any points
	 */
	protected Clusterable centroid_of(const Cluster<T> cluster)
	{
		const List<T> points = cluster.get_points();
		if (points.is_empty())
		{
			return NULL;
		}

		// in case the cluster is of type Centroid_Cluster, no need to compute the centroid
		if (dynamic_cast<const Centroid_Cluster*>(*points) != nullptr)
		{
			return ((Centroid_Cluster<T>) cluster).get_center();
		}

		const int dimension = points.get(0).get_point().size();
		const std::vector<double> centroid = std::vector<double>(dimension];
		for (const T p : points)
		{
			const std::vector<double> point = p.get_point();
			for (int i{}; i < centroid.size(); i++)
			{
				centroid[i] += point[i];
			}
		}
		for (int i{}; i < centroid.size(); i++)
		{
			centroid[i] /= points.size();
		}
		return Double_Point(centroid);
	}
}
