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

  //import java.util.Collection;
  //import java.util.List;

  //import org.hipparchus.clustering.evaluation.Cluster_Evaluator;
  //import org.hipparchus.clustering.evaluation.Sum_Of_Cluster_Variances;
  //import org.hipparchus.exception.;
  //import org.hipparchus.exception.Math_Illegal_State_Exception;

  /**
   * A wrapper around a k-means++ clustering algorithm which performs multiple trials
   * and returns the best solution.
   * @param <T> type of the points to cluster
   */
class MultiK_Means_Plus_Plus_Clusterer<T extends Clusterable> extends Clusterer<T>
{
	/** The underlying k-means clusterer. */
	private const K_Means_Plus_Plus_Clusterer<T> clusterer;

	/** The number of trial runs. */
	private const int& num_trials;

	/** The cluster evaluator to use. */
	private const Cluster_Evaluator<T> evaluator;

	/** Build a clusterer.
	 * @param clusterer the k-means clusterer to use
	 * @param num_trials number of trial runs
	 */
	public MultiK_Means_Plus_Plus_Clusterer(const K_Means_Plus_Plus_Clusterer<T> clusterer, const int& num_trials)
	{
		this(clusterer, num_trials, Sum_Of_Cluster_Variances<T>(clusterer.get_distance_measure()));
	}

	/** Build a clusterer.
	 * @param clusterer the k-means clusterer to use
	 * @param num_trials number of trial runs
	 * @param evaluator the cluster evaluator to use
	 */
	public MultiK_Means_Plus_Plus_Clusterer(const K_Means_Plus_Plus_Clusterer<T> clusterer, const int& num_trials, const Cluster_Evaluator<T> evaluator)
	{
		super(clusterer.get_distance_measure());
		this.clusterer = clusterer;
		this.num_trials = num_trials;
		this.evaluator = evaluator;
	}

	/**
	 * Returns the embedded k-means clusterer used by this instance.
	 * @return the embedded clusterer
	 */
	public K_Means_Plus_Plus_Clusterer<T> get_clusterer()
	{
		return clusterer;
	}

	/**
	 * Returns the number of trials this instance will do.
	 * @return the number of trials
	 */
	public int get_num_trials()
	{
		return num_trials;
	}

	/**
	 * Returns the {@link Cluster_Evaluator} used to determine the "best" clustering.
	 * @return the used {@link Cluster_Evaluator}
	 */
	public Cluster_Evaluator<T> get_cluster_evaluator()
	{
		return evaluator;
	}

	/**
	 * Runs the K-means++ clustering algorithm.
	 *
	 * @param points the points to cluster
	 * @return a list of clusters containing the points
	 * @ if the data points are NULL or the number
	 *   of clusters is larger than the number of data points
	 * @Math_Illegal_State_Exception if an empty cluster is encountered and the
	 *   underlying {@link K_Means_Plus_Plus_Clusterer} has its
	 *   {@link K_Means_Plus_Plus_Clusterer.Empty_Cluster_Strategy} is set to {@code ERROR}.
	 */
	 //override
	public List<Centroid_Cluster<T>> cluster(const Collection<T> points)

	{
		// at first, we have not found any clusters list yet
		List<Centroid_Cluster<T>> best = NULL;
		double best_variance_sum = INFINITY;

		// do several clustering trials
		for (int i{}; i < num_trials; ++i)
		{
			// compute a clusters list
			List<Centroid_Cluster<T>> clusters = clusterer.cluster(points);

			// compute the variance of the current list
			const double variance_sum = evaluator.score(clusters);

			if (evaluator.is_better_score(variance_sum, best_variance_sum))
			{
				// this one is the best we have found so far, remember it
				best = clusters;
				best_variance_sum = variance_sum;
			}
		}

		// return the best clusters list found
		return best;
	}
}
