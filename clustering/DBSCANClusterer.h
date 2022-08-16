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
//import java.util.Hash_Map;
//import java.util.Hash_Set;
//import java.util.List;
//import java.util.Map;
//import java.util.Set;

//import org.hipparchus.clustering.distance.Distance_Measure;
//import org.hipparchus.clustering.distance.Euclidean_Distance;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.;
//import org.hipparchus.util.Math_Utils;

/**
 * DBSCAN (density-based spatial clustering of applications with noise) algorithm.
 * <p>
 * The DBSCAN algorithm forms clusters based on the idea of density connectivity, i.e.
 * a point p is density connected to another point q, if there exists a chain of
 * points p<sub>i</sub>, with i = 1 .. n and p<sub>1</sub> = p and p<sub>n</sub> = q, * such that each pair &lt;p<sub>i</sub>, p<sub>i+1</sub>&gt; is directly density-reachable.
 * A point q is directly density-reachable from point p if it is in the &epsilon;-neighborhood
 * of this point.
 * <p>
 * Any point that is not density-reachable from a formed cluster is treated as noise, and
 * will thus not be present in the result.
 * <p>
 * The algorithm requires two parameters:
 * <ul>
 *   <li>eps: the distance that defines the &epsilon;-neighborhood of a point
 *   <li>min_points: the minimum number of density-connected points required to form a cluster
 * </ul>
 *
 * @param <T> type of the points to cluster
 * @see <a href="http://en.wikipedia.org/wiki/DBSCAN">DBSCAN (wikipedia)</a>
 * @see <a href="http://www.dbs.ifi.lmu.de/Publikationen/Papers/KDD-96.const.frame.pdf">
 * A Density-Based Algorithm for Discovering Clusters in Large Spatial Databases with Noise</a>
 */
class DBSCAN_Clusterer<T extends Clusterable> extends Clusterer<T> 
{

    /** Maximum radius of the neighborhood to be considered. */
    private const double              eps;

    /** Minimum number of points needed for a cluster. */
    private const int                 min_pts;

    /** Status of a point during the clustering process. */
    private enum Point_Status 
    {
        /** The point has is considered to be noise. */
        NOISE, /** The point is already part of a cluster. */
        PART_OF_CLUSTER
    }

    /**
     * Creates a instance of a DBSCAN_Clusterer.
     * <p>
     * The euclidean distance will be used as default distance measure.
     *
     * @param eps maximum radius of the neighborhood to be considered
     * @param min_pts minimum number of points needed for a cluster
     * @ if {@code eps < 0.0} or {@code min_pts < 0}
     */
    public DBSCAN_Clusterer(const double eps, const int min_pts)
         
        {
        this(eps, min_pts, Euclidean_Distance());
    }

    /**
     * Creates a instance of a DBSCAN_Clusterer.
     *
     * @param eps maximum radius of the neighborhood to be considered
     * @param min_pts minimum number of points needed for a cluster
     * @param measure the distance measure to use
     * @ if {@code eps < 0.0} or {@code min_pts < 0}
     */
    public DBSCAN_Clusterer(const double eps, const int min_pts, const Distance_Measure measure)
         
        {
        super(measure);

        if (eps < 0.0) 
        {
            throw std::exception("not implmented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL, eps, 0);
        }
        if (min_pts < 0) 
        {
            throw std::exception("not implmented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL, min_pts, 0);
        }
        this.eps = eps;
        this.min_pts = min_pts;
    }

    /**
     * Returns the maximum radius of the neighborhood to be considered.
     * @return maximum radius of the neighborhood
     */
    public double get_eps() 
    {
        return eps;
    }

    /**
     * Returns the minimum number of points needed for a cluster.
     * @return minimum number of points needed for a cluster
     */
    public int get_min_pts() 
    {
        return min_pts;
    }

    /**
     * Performs DBSCAN cluster analysis.
     *
     * @param points the points to cluster
     * @return the list of clusters
     * @ if the data points are NULL
     */
    //override
    public List<Cluster<T>> cluster(const Collection<T> points)  
    {

        // sanity checks
        //Math_Utils::check_not_null(points);

        const List<Cluster<T>> clusters = Array_list<>();
        const Map<Clusterable, Point_Status> visited = Hash_Map<>();

        for (const T point : points) 
        {
            if (visited.get(point) != NULL) 
            {
                continue;
            }
            const List<T> neighbors = get_neighbors(point, points);
            if (neighbors.size() >= min_pts) 
            {
                // DBSCAN does not care about center points
                const Cluster<T> cluster = Cluster<>();
                clusters.add(expand_cluster(cluster, point, neighbors, points, visited));
            }
else 
            {
                visited.put(point, Point_Status.NOISE);
            }
        }

        return clusters;
    }

    /**
     * Expands the cluster to include density-reachable items.
     *
     * @param cluster Cluster to expand
     * @param point Point to add to cluster
     * @param neighbors List of neighbors
     * @param points the data set
     * @param visited the set of already visited points
     * @return the expanded cluster
     */
    private Cluster<T> expand_cluster(const Cluster<T> cluster, const T point, const List<T> neighbors, const Collection<T> points, const Map<Clusterable, Point_Status> visited) 
    {
        cluster.add_point(point);
        visited.put(point, Point_Status.PART_OF_CLUSTER);

        List<T> seeds = Array_list<>(neighbors);
        int index = 0;
        while (index < seeds.size()) 
        {
            const T current = seeds.get(index);
            Point_Status p_status = visited.get(current);
            // only check non-visited points
            if (p_status == NULL) 
            {
                const List<T> current_neighbors = get_neighbors(current, points);
                if (current_neighbors.size() >= min_pts) 
                {
                    seeds = merge(seeds, current_neighbors);
                }
            }

            if (p_status != Point_Status.PART_OF_CLUSTER) 
            {
                visited.put(current, Point_Status.PART_OF_CLUSTER);
                cluster.add_point(current);
            }

            index++;
        }
        return cluster;
    }

    /**
     * Returns a list of density-reachable neighbors of a {@code point}.
     *
     * @param point the point to look for
     * @param points possible neighbors
     * @return the List of neighbors
     */
    private List<T> get_neighbors(const T point, const Collection<T> points) 
    {
        const List<T> neighbors = Array_list<>();
        for (const T neighbor : points) 
        {
            if (point != neighbor && distance(neighbor, point) <= eps) 
            {
                neighbors.add(neighbor);
            }
        }
        return neighbors;
    }

    /**
     * Merges two lists together.
     *
     * @param one first list
     * @param two second list
     * @return merged lists
     */
    private List<T> merge(const List<T> one, const List<T> two) 
    {
        const Set<T> one_set = Hash_Set<>(one);
        for (T item : two) 
        {
            if (!one_set.contains(item)) 
            {
                one.add(item);
            }
        }
        return one;
    }
}


