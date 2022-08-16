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

  //import java.util.Array_list;
  //import java.util.Collection;
  //import java.util.Collections;
  //import java.util.List;

  //import org.hipparchus.clustering.distance.Distance_Measure;
  //import org.hipparchus.clustering.distance.Euclidean_Distance;
  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.exception.Math_Illegal_State_Exception;
  //import org.hipparchus.random.JDKRandom_Generator;
  //import org.hipparchus.random.Random_Generator;
  //import org.hipparchus.stat.descriptive.moment.Variance;
  //import org.hipparchus.util.Math_Utils;

  /**
   * Clustering algorithm based on David Arthur and Sergei Vassilvitski k-means++ algorithm.
   * @param <T> type of the points to cluster
   * @see <a href="http://en.wikipedia.org/wiki/K-means%2B%2B">K-means++ (wikipedia)</a>
   */
class K_Means_Plus_Plus_Clusterer<T extends Clusterable> extends Clusterer<T>
{
	/** Strategies to use for replacing an empty cluster. */
	enum Empty_Cluster_Strategy
	{
		/** Split the cluster with largest distance variance. */
		LARGEST_VARIANCE,
		/** Split the cluster with largest number of points. */
		LARGEST_POINTS_NUMBER,
		/** Create a cluster around the point farthest from its centroid. */
		FARTHEST_POINT,
		/** Generate an error. */
		ERROR
	}

	/** The number of clusters. */
	private const int& k;

	/** The maximum number of iterations. */
	private const int max_iterations;

	/** Random generator for choosing initial centers. */
	private const Random_Generator random;

	/** Selected strategy for empty clusters. */
	private const Empty_Cluster_Strategy empty_strategy;

	/** Build a clusterer.
	 * <p>
	 * The default strategy for handling empty clusters that may appear during
	 * algorithm iterations is to split the cluster with largest distance variance.
	 * <p>
	 * The euclidean distance will be used as default distance measure.
	 *
	 * @param k the number of clusters to split the data into
	 */
	public K_Means_Plus_Plus_Clusterer(const int& k)
	{
		this(k, -1);
	}

	/** Build a clusterer.
	 * <p>
	 * The default strategy for handling empty clusters that may appear during
	 * algorithm iterations is to split the cluster with largest distance variance.
	 * <p>
	 * The euclidean distance will be used as default distance measure.
	 *
	 * @param k the number of clusters to split the data into
	 * @param max_iterations the maximum number of iterations to run the algorithm for.
	 *   If negative, no maximum will be used.
	 */
	public K_Means_Plus_Plus_Clusterer(const int& k, const int max_iterations)
	{
		this(k, max_iterations, Euclidean_Distance());
	}

	/** Build a clusterer.
	 * <p>
	 * The default strategy for handling empty clusters that may appear during
	 * algorithm iterations is to split the cluster with largest distance variance.
	 *
	 * @param k the number of clusters to split the data into
	 * @param max_iterations the maximum number of iterations to run the algorithm for.
	 *   If negative, no maximum will be used.
	 * @param measure the distance measure to use
	 */
	public K_Means_Plus_Plus_Clusterer(const int& k, const int max_iterations, const Distance_Measure measure)
	{
		this(k, max_iterations, measure, JDKRandom_Generator());
	}

	/** Build a clusterer.
	 * <p>
	 * The default strategy for handling empty clusters that may appear during
	 * algorithm iterations is to split the cluster with largest distance variance.
	 *
	 * @param k the number of clusters to split the data into
	 * @param max_iterations the maximum number of iterations to run the algorithm for.
	 *   If negative, no maximum will be used.
	 * @param measure the distance measure to use
	 * @param random random generator to use for choosing initial centers
	 */
	public K_Means_Plus_Plus_Clusterer(const int& k, const int max_iterations, const Distance_Measure measure, const Random_Generator random)
	{
		this(k, max_iterations, measure, random, Empty_Cluster_Strategy.LARGEST_VARIANCE);
	}

	/** Build a clusterer.
	 *
	 * @param k the number of clusters to split the data into
	 * @param max_iterations the maximum number of iterations to run the algorithm for.
	 *   If negative, no maximum will be used.
	 * @param measure the distance measure to use
	 * @param random random generator to use for choosing initial centers
	 * @param empty_strategy strategy to use for handling empty clusters that
	 * may appear during algorithm iterations
	 */
	public K_Means_Plus_Plus_Clusterer(const int& k, const int max_iterations, const Distance_Measure measure, const Random_Generator random, const Empty_Cluster_Strategy empty_strategy)
	{
		super(measure);
		this.k = k;
		this.max_iterations = max_iterations;
		this.random = random;
		this.empty_strategy = empty_strategy;
	}

	/**
	 * Return the number of clusters this instance will use.
	 * @return the number of clusters
	 */
	public int get_k()
	{
		return k;
	}

	/**
	 * Returns the maximum number of iterations this instance will use.
	 * @return the maximum number of iterations, or -1 if no maximum is set
	 */
	public int get_max_iterations()
	{
		return max_iterations;
	}

	/**
	 * Returns the random generator this instance will use.
	 * @return the random generator
	 */
	public Random_Generator get_random_generator()
	{
		return random;
	}

	/**
	 * Returns the {@link Empty_Cluster_Strategy} used by this instance.
	 * @return the {@link Empty_Cluster_Strategy}
	 */
	public Empty_Cluster_Strategy get_empty_cluster_strategy()
	{
		return empty_strategy;
	}

	/**
	 * Runs the K-means++ clustering algorithm.
	 *
	 * @param points the points to cluster
	 * @return a list of clusters containing the points
	 * @ if the data points are NULL or the number
	 *     of clusters is larger than the number of data points
	 * @Math_Illegal_State_Exception if an empty cluster is encountered and the
	 * {@link #empty_strategy} is set to {@code ERROR}
	 */
	 //override
	public List<Centroid_Cluster<T>> cluster(const Collection<T> points)
		
	{
		// sanity checks
		//Math_Utils::check_not_null(points);

		// number of clusters has to be smaller or equal the number of data points
		if (points.size() < k)
		{
			throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL_BOUND_EXCLUDED, points.size(), k);
		}

	// create the initial clusters
	List<Centroid_Cluster<T>> clusters = choose_initial_centers(points);

	// create an array containing the latest assignment of a point to a cluster
	// no need to initialize the array, as it will be filled with the first assignment
	std::vector<int> assignments = int[points.size()];
	assign_points_to_clusters(clusters, points, assignments);

	// iterate through updating the centers until we're done
	const int max = (max_iterations < 0) ? std::numeric_limits<int>::max() : max_iterations;
	for (const int& count = 0; count < max; count++)
	{
		bool empty_cluster = false;
		List<Centroid_Cluster<T>> new_clusters = Array_list<>();
		for (const Centroid_Cluster<T> cluster : clusters)
		{
			const Clusterable new_center;
			if (cluster.get_points().is_empty())
			{
				switch (empty_strategy)
				{
					case LARGEST_VARIANCE:
						new_center = get_point_from_largest_variance_cluster(clusters);
						break;
					case LARGEST_POINTS_NUMBER:
						new_center = get_point_from_largest_number_cluster(clusters);
						break;
					case FARTHEST_POINT:
						new_center = get_farthest_point(clusters);
						break;
					default:
						throw Math_Illegal_State_Exception(LocalizedClusteringFormats.EMPTY_CLUSTER_IN_K_MEANS);
				}
				empty_cluster = true;
			}
else
				{
					new_center = centroid_of(cluster.get_points(), cluster.get_center().get_point().size());
				}
				new_clusters.add(new Centroid_Cluster<T>(new_center));
			}
			int changes = assign_points_to_clusters(new_clusters, points, assignments);
			clusters = new_clusters;

			// if there were no more changes in the point-to-cluster assignment
			// and there are no empty clusters left, return the current clusters
			if (changes == 0 && !empty_cluster)
			{
				return clusters;
			}
		}
		return clusters;
	}

		/**
		 * Adds the given points to the closest {@link Cluster}.
		 *
		 * @param clusters the {@link Cluster}s to add the points to
		 * @param points the points to add to the given {@link Cluster}s
		 * @param assignments points assignments to clusters
		 * @return the number of points assigned to different clusters as the iteration before
		 */
		private int assign_points_to_clusters(const List<Centroid_Cluster<T>> clusters, const Collection<T> points, const std::vector<int> assignments)
	{
		int assigned_differently = 0;
		int point_index = 0;
		for (const T p : points)
		{
			int cluster_index = get_nearest_cluster(clusters, p);
			if (cluster_index != assignments[point_index])
			{
				assigned_differently++;
			}

			Centroid_Cluster<T> cluster = clusters.get(cluster_index);
			cluster.add_point(p);
			assignments[point_index++] = cluster_index;
		}

		return assigned_differently;
	}

	/**
	 * Use K-means++ to choose the initial centers.
	 *
	 * @param points the points to choose the initial centers from
	 * @return the initial centers
	 */
	private List<Centroid_Cluster<T>> choose_initial_centers(const Collection<T> points)
	{
		// Convert to list for indexed access. Make it unmodifiable, since removal of items
		// would screw up the logic of this method.
		const List<T> point_list = Collections.unmodifiable_list(new Array_list<T>(points));

		// The number of points in the list.
		const int& num_points = point_list.size();

		// Set the corresponding element in this array to indicate when
		// elements of point_list are no longer available.
		const bool[] taken = bool[num_points];

		// The resulting list of initial centers.
		const List<Centroid_Cluster<T>> result_set = Array_list<>();

		// Choose one center uniformly at random from among the data points.
		const int first_point_index = random.next_int(num_points);

		const T first_point = point_list.get(first_point_index);

		result_set.add(new Centroid_Cluster<T>(first_point));

		// Must mark it as taken
		taken[first_point_index] = true;

		// To keep track of the minimum distance squared of elements of
		// point_list to elements of result_set.
		const std::vector<double> min_dist_squared = std::vector<double>(num_points];

		// Initialize the elements.  sin_ce the only point in result_set is first_point, // this is very easy.
		for (int i{}; i < num_points; i++)
		{
			if (i != first_point_index) { // That point isn't considered
				double d = distance(first_point, point_list.get(i));
				min_dist_squared[i] = d * d;
			}
		}

		while (result_set.size() < k)
		{
			// Sum up the squared distances for the points in point_list not
			// already taken.
			double dist_sq_sum = 0.0;

			for (int i{}; i < num_points; i++)
			{
				if (!taken[i])
				{
					dist_sq_sum += min_dist_squared[i];
				}
			}

			// Add one data point as a center. Each point x is chosen with
			// probability proportional to D(x)2
			const double r = random.next_double() * dist_sq_sum;

			// The index of the next point to be added to the result_set.
			int next_point_index = -1;

			// Sum through the squared min distances again, stopping when
			// sum >= r.
			double sum = 0.0;
			for (int i{}; i < num_points; i++)
			{
				if (!taken[i])
				{
					sum += min_dist_squared[i];
					if (sum >= r)
					{
						next_point_index = i;
						break;
					}
				}
			}

			// If it's not set to >= 0, the point wasn't found in the previous
			// for loop, probably because distances are extremely small.  Just pick
			// the last available point.
			if (next_point_index == -1)
			{
				for (int i = num_points - 1; i >= 0; i--)
				{
					if (!taken[i])
					{
						next_point_index = i;
						break;
					}
				}
			}

			// We found one.
			if (next_point_index >= 0)
			{
				const T p = point_list.get(next_point_index);

				result_set.add(new Centroid_Cluster<T>(p));

				// Mark it as taken.
				taken[next_point_index] = true;

				if (result_set.size() < k)
				{
					// Now update elements of min_dist_squared.  We only have to compute
					// the distance to the center to do this.
					for (int j{}; j < num_points; j++)
					{
						// Only have to worry about the points still not taken.
						if (!taken[j])
						{
							double d = distance(p, point_list.get(j));
							double d2 = d * d;
							if (d2 < min_dist_squared[j])
							{
								min_dist_squared[j] = d2;
							}
						}
					}
				}
			}
			else
			{
				// None found --
				// Break from the while loop to prevent
				// an infinite loop.
				break;
			}
		}

		return result_set;
	}

	/**
	 * Get a random point from the {@link Cluster} with the largest distance variance.
	 *
	 * @param clusters the {@link Cluster}s to search
	 * @return a random point from the selected cluster
	 * @Math_Illegal_State_Exception if clusters are all empty
	 */
	private T get_point_from_largest_variance_cluster(const Collection<Centroid_Cluster<T>> clusters)
		Math_Illegal_State_Exception
	{
		double max_variance = -INFINITY;
		Cluster<T> selected = NULL;
		for (const Centroid_Cluster<T> cluster : clusters)
		{
			if (!cluster.get_points().is_empty())
			{
				// compute the distance variance of the current cluster
				const Clusterable center = cluster.get_center();
				const Variance stat = Variance();
				for (const T point : cluster.get_points())
				{
					stat.increment(distance(point, center));
				}
				const double variance = stat.get_result();

				// select the cluster with the largest variance
				if (variance > max_variance)
				{
					max_variance = variance;
					selected = cluster;
				}
			}
		}

		// did we find at least one non-empty cluster ?
		if (selected == NULL)
		{
			throw Math_Illegal_State_Exception(LocalizedClusteringFormats.EMPTY_CLUSTER_IN_K_MEANS);
		}

		// extract a random point from the cluster
		const List<T> selected_points = selected.get_points();
		return selected_points.remove(random.next_int(selected_points.size()));
	}

	/**
	 * Get a random point from the {@link Cluster} with the largest number of points
	 *
	 * @param clusters the {@link Cluster}s to search
	 * @return a random point from the selected cluster
	 * @Math_Illegal_State_Exception if clusters are all empty
	 */
	private T get_point_from_largest_number_cluster(const Collection< ? extends Cluster<T>> clusters)
		Math_Illegal_State_Exception
	{
		int max_number = 0;
		Cluster<T> selected = NULL;
		for (const Cluster<T> cluster : clusters)
		{
			// get the number of points of the current cluster
			const int& number = cluster.get_points().size();

			// select the cluster with the largest number of points
			if (number > max_number)
			{
				max_number = number;
				selected = cluster;
			}
		}

		// did we find at least one non-empty cluster ?
		if (selected == NULL)
		{
			throw Math_Illegal_State_Exception(LocalizedClusteringFormats.EMPTY_CLUSTER_IN_K_MEANS);
		}

		// extract a random point from the cluster
		const List<T> selected_points = selected.get_points();
		return selected_points.remove(random.next_int(selected_points.size()));
	}

	/**
	 * Get the point farthest to its cluster center
	 *
	 * @param clusters the {@link Cluster}s to search
	 * @return point farthest to its cluster center
	 * @Math_Illegal_State_Exception if clusters are all empty
	 */
	private T get_farthest_point(const Collection<Centroid_Cluster<T>> clusters) Math_Illegal_State_Exception
	{
		double max_distance = -INFINITY;
		Cluster<T> selected_cluster = NULL;
		int selected_point = -1;
		for (const Centroid_Cluster<T> cluster : clusters)
		{
			// get the farthest point
			const Clusterable center = cluster.get_center();
			const List<T> points = cluster.get_points();
			for (int i{}; i < points.size(); ++i)
			{
				const double distance = distance(points.get(i), center);
				if (distance > max_distance)
				{
					max_distance = distance;
					selected_cluster = cluster;
					selected_point = i;
				}
			}
		}

		// did we find at least one non-empty cluster ?
		if (selected_cluster == NULL)
		{
			throw Math_Illegal_State_Exception(LocalizedClusteringFormats.EMPTY_CLUSTER_IN_K_MEANS);
		}

		return selected_cluster.get_points().remove(selected_point);
	}

	/**
	 * Returns the nearest {@link Cluster} to the given point
	 *
	 * @param clusters the {@link Cluster}s to search
	 * @param point the point to find the nearest {@link Cluster} for
	 * @return the index of the nearest {@link Cluster} to the given point
	 */
	private int get_nearest_cluster(const Collection<Centroid_Cluster<T>> clusters, const T point)
	{
		double min_distance = Double.MAX_VALUE;
		int cluster_index = 0;
		int min_cluster = 0;
		for (const Centroid_Cluster<T> c : clusters)
		{
			const double distance = distance(point, c.get_center());
			if (distance < min_distance)
			{
				min_distance = distance;
				min_cluster = cluster_index;
			}
			cluster_index++;
		}
		return min_cluster;
	}

	/**
	 * Computes the centroid for a set of points.
	 *
	 * @param points the set of points
	 * @param dimension the point dimension
	 * @return the computed centroid for the set of points
	 */
	private Clusterable centroid_of(const Collection<T> points, const int& dimension)
	{
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
