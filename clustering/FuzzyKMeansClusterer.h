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
  //import org.hipparchus.linear.Matrix_Utils;
  //import org.hipparchus.linear.Real_Matrix;
  //import org.hipparchus.random.JDKRandom_Generator;
  //import org.hipparchus.random.Random_Generator;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Math_Arrays;
  //import org.hipparchus.util.Math_Utils;
#include <type_traits>
#include <cmath>
#include <vector>
#include "../core/linear/MatrixUtils.h"
#include "Clusterable.h"
/**
 * Fuzzy K-_Means clustering algorithm.
 * <p>
 * The Fuzzy K-_Means algorithm is a variation of the classical K-_Means algorithm, with the
 * major difference that a single data point is not uniquely assigned to a single cluster.
 * Instead, each point i has a set of weights u<sub>ij</sub> which indicate the degree of membership
 * to the cluster j.
 * <p>
 * The algorithm then tries to minimize the objective function:
 * <pre>
 * J = &#8721;<sub>i=1..C</sub>&#8721;<sub>k=1..N</sub> u<sub>ik</sub><sup>m</sup>d<sub>ik</sub><sup>2</sup>
 * </pre>
 * with d<sub>ik</sub> being the distance between data point i and the cluster center k.
 * <p>
 * The algorithm requires two parameters:
 * <ul>
 *   <li>k: the number of clusters
 *   <li>fuzziness: determines the level of cluster fuzziness, larger values lead to fuzzier clusters
 * </ul>
 * Additional, optional parameters:
 * <ul>
 *   <li>max_iterations: the maximum number of iterations
 *   <li>epsilon: the convergence criteria, default is 1e-3
 * </ul>
 * <p>
 * The fuzzy variant of the K-_Means algorithm is more robust with regard to the selection
 * of the initial cluster centers.
 *
 * @param <T> type of the points to cluster
 */
template<typename T, typename std::enable_if<std::is_base_of<Clusterable, T>::value>::type* = nullptr>
class Fuzzy_K_Means_Clusterer : public Clusterer<T>
{
private:
	/** The default value for the convergence criteria. */
	static constexpr double DEFAULT_EPSILON{ 1e-3 };

	/** The number of clusters. */
	const int my_k;

	/** The maximum number of iterations. */
	const int my_max_iterations;

	/** The fuzziness factor. */
	const double my_fuzziness;

	/** The convergence criteria. */
	const double my_epsilon;

	/** Random generator for choosing initial centers. */
	const Random_Generator my_random;

	/** The membership matrix. */
	std::vector<std::vector<double>> my_membership_matrix;

	/** The list of points used in the last call to {@link #cluster(Collection)}. */
	List<T> points;

	/** The list of clusters resulting from the last call to {@link #cluster(Collection)}. */
	List<Centroid_Cluster<T>> my_clusters;

	/**
	 * Update the cluster centers.
	 */
	void update_cluster_centers()
	{
		int j{ 0 };
		const List<Centroid_Cluster<T>> new_clusters = Array_list<>(k);
		for (const Centroid_Cluster<T> cluster : clusters)
		{
			const Clusterable center = cluster.get_center();
			int i{};
			std::vector<double> arr = std::vector<double>(center.get_point().size()];
			double sum{};
			for (const T& point : points)
			{
				const double u = std::pow(membership_matrix[i][j], fuzziness);
				const std::vector<double> point_arr = point.get_point();
				for (const int& idx = 0; idx < arr.size(); idx++)
				{
					arr[idx] += u * point_arr[idx];
				}
				sum += u;
				i++;
			}
			Math_Arrays::scale_in_place(1.0 / sum, arr);
			new_clusters.add(new Centroid_Cluster<T>(new Double_Point(arr)));
			j++;
		}
		clusters.clear();
		clusters = new_clusters;
	}

	/**
	 * Updates the membership matrix and assigns the points to the cluster with
	 * the highest membership.
	 */
	void update_membership_matrix()
	{
		for (int i{}; i < points.size(); i++)
		{
			const T point = points.get(i);
			double max_membership = Double.MIN_VALUE;
			int new_cluster = -1;
			for (int j{}; j < clusters.size(); j++)
			{
				double sum{};
				const double dist_a = std::abs(distance(point, clusters.get(j).get_center()));

				if (dist_a != 0.0)
				{
					for (const Centroid_Cluster<T> c : clusters)
					{
						const double dist_b = std::abs(distance(point, c.get_center()));
						if (dist_b == 0.0)
						{
							sum = INFINITY;
							break;
						}
						sum += std::pow(dist_a / dist_b, 2.0 / (fuzziness - 1.0));
					}
				}

				double membership;
				if (sum == 0.0)
				{
					membership = 1.0;
				}
				else if (std::isinf(sum)
				{
					membership = 0.0;
				}
				else
				{
					membership = 1.0 / sum;
				}
				membership_matrix[i][j] = membership;

				if (membership_matrix[i][j] > max_membership)
				{
					max_membership = membership_matrix[i][j];
					new_cluster = j;
				}
			}
			clusters.get(new_cluster).add_point(point);
		}
	}

	/**
	 * Initialize the membership matrix with random values.
	 */
	void initialize_membership_matrix()
	{
		for (int i{}; i < points.size(); i++)
		{
			for (int j{}; j < k; j++)
			{
				membership_matrix[i][j] = random.next_double();
			}
			membership_matrix[i] = Math_Arrays::normalize_array(membership_matrix[i], 1.0);
		}
	}

	/**
	 * Calculate the maximum element-by-element change of the membership matrix
	 * for the current iteration.
	 *
	 * @param matrix the membership matrix of the previous iteration
	 * @return the maximum membership matrix change
	 */
	double calculate_max_membership_change(const std::vector<std::vector<double>>& matrix)
	{
		double max_membership = 0.0;
		for (int i{}; i < points.size(); i++)
		{
			for (int j{}; j < clusters.size(); j++)
			{
				double v = std::abs(membership_matrix[i][j] - matrix[i][j]);
				max_membership = std::max(v, max_membership);
			}
		}
		return max_membership;
	}

	/**
	 * Copy the membership matrix into the provided matrix.
	 *
	 * @param matrix the place to store the membership matrix
	 */
	void save_membership_matrix(const std::vector<std::vector<double>>& matrix)
	{
		for (int i{}; i < points.size(); i++)
		{
			System.arraycopy(membership_matrix[i], 0, matrix[i], 0, clusters.size());
		}
	}

public:
	/**
	 * Creates a instance of a Fuzzy_K_Means_Clusterer.
	 * <p>
	 * The euclidean distance will be used as default distance measure.
	 *
	 * @param k the number of clusters to split the data into
	 * @param fuzziness the fuzziness factor, must be &gt; 1.0
	 * @ if {@code fuzziness <= 1.0}
	 */
	Fuzzy_K_Means_Clusterer(const int& k, const double fuzziness)
	{
		this(k, fuzziness, -1, Euclidean_Distance());
	}

	/**
	 * Creates a instance of a Fuzzy_K_Means_Clusterer.
	 *
	 * @param k the number of clusters to split the data into
	 * @param fuzziness the fuzziness factor, must be &gt; 1.0
	 * @param max_iterations the maximum number of iterations to run the algorithm for.
	 *   If negative, no maximum will be used.
	 * @param measure the distance measure to use
	 * @ if {@code fuzziness <= 1.0}
	 */
	Fuzzy_K_Means_Clusterer(const int& k, const double fuzziness, const int max_iterations, const Distance_Measure measure)

	{
		this(k, fuzziness, max_iterations, measure, DEFAULT_EPSILON, JDKRandom_Generator());
	}

	/**
	 * Creates a instance of a Fuzzy_K_Means_Clusterer.
	 *
	 * @param k the number of clusters to split the data into
	 * @param fuzziness the fuzziness factor, must be &gt; 1.0
	 * @param max_iterations the maximum number of iterations to run the algorithm for.
	 *   If negative, no maximum will be used.
	 * @param measure the distance measure to use
	 * @param epsilon the convergence criteria (default is 1e-3)
	 * @param random random generator to use for choosing initial centers
	 * @ if {@code fuzziness <= 1.0}
	 */
	Fuzzy_K_Means_Clusterer(const int& k, const double& fuzziness, const int& max_iterations, const Distance_Measure& measure, const double& epsilon, const Random_Generator& random)
	{
		super(measure);

		if (fuzziness <= 1.0)
		{
			throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL_BOUND_EXCLUDED, fuzziness, 1.0);
		}
		this.k = k;
		this.fuzziness = fuzziness;
		this.max_iterations = max_iterations;
		this.epsilon = epsilon;
		this.random = random;

		this.membership_matrix = NULL;
		this.points = NULL;
		this.clusters = NULL;
	}

	/**
	 * Return the number of clusters this instance will use.
	 * @return the number of clusters
	 */
	int get_k() const
	{
		return my_k;
	}

	/**
	 * Returns the fuzziness factor used by this instance.
	 * @return the fuzziness factor
	 */
	double get_fuzziness() const
	{
		return my_fuzziness;
	}

	/**
	 * Returns the maximum number of iterations this instance will use.
	 * @return the maximum number of iterations, or -1 if no maximum is set
	 */
	int get_max_iterations() const
	{
		return my_max_iterations;
	}

	/**
	 * Returns the convergence criteria used by this instance.
	 * @return the convergence criteria
	 */
	double get_epsilon() const
	{
		return my_epsilon;
	}

	/**
	 * Returns the random generator this instance will use.
	 * @return the random generator
	 */
	Random_Generator get_random_generator() const
	{
		return my_random;
	}

	/**
	 * Returns the {@code nxk} membership matrix, where {@code n} is the number
	 * of data points and {@code k} the number of clusters.
	 * <p>
	 * The element U<sub>i,j</sub> represents the membership value for data point {@code i}
	 * to cluster {@code j}.
	 *
	 * @return the membership matrix
	 * @Math_Illegal_State_Exception if {@link #cluster(Collection)} has not been called before
	 */
	Real_Matrix get_membership_matrix()
	{
		if (membership_matrix == NULL)
		{
			throw std::exception("not implemented");
			//throw Math_Illegal_State_Exception(Localized_Core_Formats::ILLEGAL_STATE);
		}
		return Matrix_Utils::create_real_matrix(membership_matrix);
	}

	/**
	 * Returns an unmodifiable list of the data points used in the last
	 * call to {@link #cluster(Collection)}.
	 * @return the list of data points, or {@code NULL} if {@link #cluster(Collection)} has
	 *   not been called before.
	 */
	std::vector<T> get_data_points() const
	{
		return my_points;
	}

	/**
	 * Returns the list of clusters resulting from the last call to {@link #cluster(Collection)}.
	 * @return the list of clusters, or {@code NULL} if {@link #cluster(Collection)} has
	 *   not been called before.
	 */
	std::vector<Centroid_Cluster<T>> get_clusters() const
	{
		return my_clusters;
	}

	/**
	 * Get the value of the objective function.
	 * @return the objective function evaluation as double value
	 * @Math_Illegal_State_Exception if {@link #cluster(Collection)} has not been called before
	 */
	double get_objective_function_value()
	{
		if (my_points == NULL || my_clusters == NULL)
		{
			throw Math_Illegal_State_Exception(hipparchus::exception::Localized_Core_Formats_Type::ILLEGAL_STATE);
		}

		int i{};
		double obj_function{};
		for (const T& point : my_points)
		{
			int j{};
			for (const Centroid_Cluster<T>& cluster : my_clusters)
			{
				const double dist = distance(point, cluster.get_center());
				obj_function += (dist * dist) * std::pow(membership_matrix[i][j], fuzziness);
				j++;
			}
			i++;
		}
		return obj_function;
	}

	/**
	 * Performs Fuzzy K-_Means cluster analysis.
	 *
	 * @param data_points the points to cluster
	 * @return the list of clusters
	 * @ if the data points are NULL or the number
	 *     of clusters is larger than the number of data points
	 */
	 //override
	std::vector<Centroid_Cluster<T>> > cluster(const Collection<T>& data_points)
	{
		// sanity checks
		//Math_Utils::check_not_null(data_points);

		const int size{ data_points.size() };

		// number of clusters has to be smaller or equal the number of data points
		if (size < k)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL_BOUND_EXCLUDED, size, k);
		}

		// copy the input collection to an unmodifiable list with indexed access
		points = Collections.unmodifiable_list(new Array_list<>(data_points));
		clusters = Array_list<>();
		membership_matrix = std::vector<std::vector<double>>(size, std::vector<double>(k));
		auto old_matrix = std::vector<std::vector<double>>(size, std::vector<double>(k));

		// if no points are provided, return an empty list of clusters
		if (size == 0)
		{
			return clusters;
		}

		initialize_membership_matrix();

		// there is at least one point
		const int point_dimension = points.get(0).get_point().size();
		for (int i{}; i < k; i++)
		{
			clusters.add(new Centroid_Cluster<T>(new Double_Point(std::vector<double>(point_dimension])));
		}

		int iteration = 0;
		const int max = (max_iterations < 0) ? std::numeric_limits<int>::max() : max_iterations;
		double difference;

		do
		{
			save_membership_matrix(old_matrix);
			update_cluster_centers();
			update_membership_matrix();
			difference = calculate_max_membership_change(old_matrix);
		} while (difference > epsilon && ++iteration < max);

		return clusters;
	}
};
