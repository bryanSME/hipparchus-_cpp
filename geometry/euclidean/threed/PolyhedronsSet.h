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
//package org.hipparchus.geometry.euclidean.threed;

//import java.util.Array_list;
//import java.util.Arrays;
//import java.util.Collection;
//import java.util.List;

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Runtime_Exception;
//import org.hipparchus.geometry.Localized_Geometry_Formats;
//import org.hipparchus.geometry.Point;
//import org.hipparchus.geometry.euclidean.oned.Euclidean_1D;
//import org.hipparchus.geometry.euclidean.twod.Euclidean_2D;
//import org.hipparchus.geometry.euclidean.twod.Polygons_Set;
//import org.hipparchus.geometry.euclidean.twod.Sub_Line;
//import org.hipparchus.geometry.euclidean.twod.Vector_2D;
//import org.hipparchus.geometry.partitioning.Abstract_Region;
//import org.hipparchus.geometry.partitioning.BSP_Tree;
//import org.hipparchus.geometry.partitioning.BSP_Tree_Visitor;
//import org.hipparchus.geometry.partitioning.Boundary_Attribute;
//import org.hipparchus.geometry.partitioning.Hyperplane;
//import org.hipparchus.geometry.partitioning.Region;
//import org.hipparchus.geometry.partitioning.Region_Factory;
//import org.hipparchus.geometry.partitioning.Sub_Hyperplane;
//import org.hipparchus.geometry.partitioning.Transform;
//import org.hipparchus.util.FastMath;

/** This class represents a 3D region: a set of polyhedrons.
 */
class Polyhedrons_Set extends Abstract_Region<Euclidean_3D, Euclidean_2D> 
{

    /** Build a polyhedrons set representing the whole real line.
     * @param tolerance tolerance below which points are considered identical
     */
    public Polyhedrons_Set(const double& tolerance) 
    {
        super(tolerance);
    }

    /** Build a polyhedrons set from a BSP tree.
     * <p>The leaf nodes of the BSP tree <em>must</em> have a
     * {@code Boolean} attribute representing the inside status of
     * the corresponding cell (true for inside cells, false for outside
     * cells). In order to avoid building too many small objects, it is
     * recommended to use the predefined constants
     * {@code Boolean.TRUE} and {@code Boolean.FALSE}</p>
     * <p>
     * This constructor is aimed at expert use, as building the tree may
     * be a difficult task. It is not intended for general use and for
     * performances reasons does not check thoroughly its input, as this would
     * require walking the full tree each time. Failing to provide a tree with
     * the proper attributes, <em>will</em> therefore generate problems like
     * {@link Null_Pointer_Exception} or {@link Class_Cast_Exception} only later on.
     * This limitation is known and explains why this constructor is for expert
     * use only. The caller does have the responsibility to provided correct arguments.
     * </p>
     * @param tree inside/outside BSP tree representing the region
     * @param tolerance tolerance below which points are considered identical
     */
    public Polyhedrons_Set(const BSP_Tree<Euclidean_3D> tree, const double& tolerance) 
    {
        super(tree, tolerance);
    }

    /** Build a polyhedrons set from a Boundary RE_Presentation (B-rep) specified by sub-hyperplanes.
     * <p>The boundary is provided as a collection of {@link
     * Sub_Hyperplane sub-hyperplanes}. Each sub-hyperplane has the
     * interior part of the region on its minus side and the exterior on
     * its plus side.</p>
     * <p>The boundary elements can be in any order, and can form
     * several non-connected sets (like for example polyhedrons with holes
     * or a set of disjoint polyhedrons considered as a whole). In
     * fact, the elements do not even need to be connected together
     * (their topological connections are not used here). However, if the
     * boundary does not really separate an inside open from an outside
     * open (open having here its topological meaning), then subsequent
     * calls to the {@link Region#check_point(Point) check_point} method will
     * not be meaningful anymore.</p>
     * <p>If the boundary is empty, the region will represent the whole
     * space.</p>
     * @param boundary collection of boundary elements, as a
     * collection of {@link Sub_Hyperplane Sub_Hyperplane} objects
     * @param tolerance tolerance below which points are considered identical
     */
    public Polyhedrons_Set(const Collection<Sub_Hyperplane<Euclidean_3D>> boundary, const double& tolerance) 
    {
        super(boundary, tolerance);
    }

    /** Build a polyhedrons set from a Boundary RE_Presentation (B-rep) specified by connected vertices.
     * <p>
     * The boundary is provided as a list of vertices and a list of facets.
     * Each facet is specified as an integer array containing the arrays vertices
     * indices in the vertices list. Each facet normal is oriented by right hand
     * rule to the facet vertices list.
     * </p>
     * <p>
     * Some basic sanity checks are performed but not everything is thoroughly
     * assessed, so it remains under caller responsibility to ensure the vertices
     * and facets are consistent and properly define a polyhedrons set.
     * </p>
     * @param vertices list of polyhedrons set vertices
     * @param facets list of facets, as vertices indices in the vertices list
     * @param tolerance tolerance below which points are considered identical
     * @exception  if some basic sanity checks fail
     */
    public Polyhedrons_Set(const List<Vector_3D> vertices, const List<std::vector<int>> facets, const double& tolerance) 
    {
        super(build_boundary(vertices, facets, tolerance), tolerance);
    }

    /** Build a polyhedrons set from a Boundary RE_Presentation (B-rep) specified by connected vertices.
     * <p>
     * Some basic sanity checks are performed but not everything is thoroughly
     * assessed, so it remains under caller responsibility to ensure the vertices
     * and facets are consistent and properly define a polyhedrons set.
     * </p>
     * @param brep Boundary RE_Presentation of the polyhedron to build
     * @param tolerance tolerance below which points are considered identical
     * @exception  if some basic sanity checks fail
     * @since 1.2
     */
    public Polyhedrons_Set(const B_Rep brep, const double& tolerance) 
    {
        super(build_boundary(brep.get_vertices(), brep.get_facets(), tolerance), tolerance);
    }

    /** Build a parallellepipedic box.
     * @param x_min low bound along the x direction
     * @param x_max high bound along the x direction
     * @param y_min low bound along the y direction
     * @param y_max high bound along the y direction
     * @param z_min low bound along the z direction
     * @param z_max high bound along the z direction
     * @param tolerance tolerance below which points are considered identical
     */
    public Polyhedrons_Set(const double x_min, const double x_max, const double y_min, const double y_max, const double z_min, const double z_max, const double& tolerance) 
    {
        super(build_boundary(x_min, x_max, y_min, y_max, z_min, z_max, tolerance), tolerance);
    }

    /** Build a parallellepipedic box boundary.
     * @param x_min low bound along the x direction
     * @param x_max high bound along the x direction
     * @param y_min low bound along the y direction
     * @param y_max high bound along the y direction
     * @param z_min low bound along the z direction
     * @param z_max high bound along the z direction
     * @param tolerance tolerance below which points are considered identical
     * @return boundary tree
     */
    private static BSP_Tree<Euclidean_3D> build_boundary(const double x_min, const double x_max, const double y_min, const double y_max, const double z_min, const double z_max, const double& tolerance) 
    {
        if ((x_min >= x_max - tolerance) || (y_min >= y_max - tolerance) || (z_min >= z_max - tolerance)) 
        {
            // too thin box, build an empty polygons set
            return BSP_Tree<Euclidean_3D>(Boolean.FALSE);
        }
        const Plane px_min = Plane(Vector_3D(x_min, 0,    0),   Vector_3D.MINUS_I, tolerance);
        const Plane px_max = Plane(Vector_3D(x_max, 0,    0),   Vector_3D.PLUS_I,  tolerance);
        const Plane py_min = Plane(Vector_3D(0,    y_min, 0),   Vector_3D.MINUS_J, tolerance);
        const Plane py_max = Plane(Vector_3D(0,    y_max, 0),   Vector_3D.PLUS_J,  tolerance);
        const Plane pz_min = Plane(Vector_3D(0,    0,   z_min), Vector_3D.MINUS_K, tolerance);
        const Plane pz_max = Plane(Vector_3D(0,    0,   z_max), Vector_3D.PLUS_K,  tolerance);
        const Region<Euclidean_3D> boundary =
        Region_Factory<Euclidean_3D>().build_convex(px_min, px_max, py_min, py_max, pz_min, pz_max);
        return boundary.get_tree(false);
    }

    /** Build boundary from vertices and facets.
     * @param vertices list of polyhedrons set vertices
     * @param facets list of facets, as vertices indices in the vertices list
     * @param tolerance tolerance below which points are considered identical
     * @return boundary as a list of sub-hyperplanes
     * @exception  if some basic sanity checks fail
     */
    private static List<Sub_Hyperplane<Euclidean_3D>> build_boundary(const List<Vector_3D> vertices, const List<std::vector<int>> facets, const double& tolerance) 
    {

        // check vertices distances
        for (int i{}; i < vertices.size() - 1; ++i) 
        {
            const Vector_3D vi = vertices.get(i);
            for (int j = i + 1; j < vertices.size(); ++j) 
            {
                if (Vector_3D.distance(vi, vertices.get(j)) <= tolerance) 
                {
                    throw (Localized_Geometry_Formats.CLOSE_VERTICES, vi.get_x(), vi.get_y(), vi.get_z());
                }
            }
        }

        // find how vertices are referenced by facets
        const std::vector<std::vector<int>> references = find_references(vertices, facets);

        // find how vertices are linked together by edges along the facets they belong to
        const std::vector<std::vector<int>> successors = successors(vertices, facets, references);

        // check edges orientations
        for (const int& v_a = 0; v_a < vertices.size(); ++v_a) 
        {
            for (const int v_b : successors[v_a]) 
            {

                if (v_b >= 0) 
                {
                    // when facets are properly oriented, if v_b is the successor of v_a on facet f1, // then there must be an adjacent facet f2 where v_a is the successor of v_b
                    bool found = false;
                    for (const int v : successors[v_b]) 
                    {
                        found = found || (v == v_a);
                    }
                    if (!found) 
                    {
                        const Vector_3D start = vertices.get(v_a);
                        const Vector_3D end   = vertices.get(v_b);
                        throw (Localized_Geometry_Formats.EDGE_CONNECTED_TO_ONE_FACET, start.get_x(), start.get_y(), start.get_z(), end.get_x(),   end.get_y(),   end.get_z());
                    }
                }
            }
        }

        const List<Sub_Hyperplane<Euclidean_3D>> boundary = Array_list<>();

        for (const std::vector<int> facet : facets) 
        {

            // define facet plane from the first 3 points
            Plane plane = Plane(vertices.get(facet[0]), vertices.get(facet[1]), vertices.get(facet[2]), tolerance);

            // check all points are in the plane
            const Vector_2D[] two_2_points = Vector_2D[facet.size()];
            for (int i = 0 ; i < facet.size(); ++i) 
            {
                const Vector_3D v = vertices.get(facet[i]);
                if (!plane.contains(v)) 
                {
                    throw (Localized_Geometry_Formats.OUT_OF_PLANE, v.get_x(), v.get_y(), v.get_z());
                }
                two_2_points[i] = plane.to_sub_space(v);
            }

            // create the polygonal facet
            boundary.add(new Sub_Plane(plane, Polygons_Set(tolerance, two_2_points)));

        }

        return boundary;

    }

    /** Find the facets that reference each edges.
     * @param vertices list of polyhedrons set vertices
     * @param facets list of facets, as vertices indices in the vertices list
     * @return references array such that r[v][k] = f for some k if facet f contains vertex v
     * @exception  if some facets have fewer than 3 vertices
     */
    private static std::vector<std::vector<int>> find_references(const List<Vector_3D> vertices, const List<std::vector<int>> facets) 
    {

        // find the maximum number of facets a vertex belongs to
        const std::vector<int> nb_facets = int[vertices.size()];
        int max_facets  = 0;
        for (const std::vector<int> facet : facets) 
        {
            if (facet.size() < 3) 
            {
                throw std::exception("not implemented");
                // throw (hipparchus::exception::Localized_Core_Formats_Type::WRONG_NUMBER_OF_POINTS, 3, facet.size(), true);
            }
            for (const int index : facet) 
            {
                max_facets = std::max(max_facets, ++nb_facets[index]);
            }
        }

        // set up the references array
        const std::vector<std::vector<int>> references = int[vertices.size()][max_facets];
        for (std::vector<int> r : references) 
        {
            Arrays.fill(r, -1);
        }
        for (const int& f = 0; f < facets.size(); ++f) 
        {
            for (const int v : facets.get(f)) 
            {
                // vertex v is referenced by facet f
                int k = 0;
                while (k < max_facets && references[v][k] >= 0) 
                {
                    ++k;
                }
                references[v][k] = f;
            }
        }

        return references;

    }

    /** Find the successors of all vertices among all facets they belong to.
     * @param vertices list of polyhedrons set vertices
     * @param facets list of facets, as vertices indices in the vertices list
     * @param references facets references array
     * @return indices of vertices that follow vertex v in some facet (the array
     * may contain extra entries at the end, set to negative indices)
     * @exception  if the same vertex appears more than
     * once in the successors list (which means one facet orientation is wrong)
     */
    private static std::vector<std::vector<int>> successors(const List<Vector_3D> vertices, const List<std::vector<int>> facets, const std::vector<std::vector<int>> references) 
    {

        // create an array large enough
        const std::vector<std::vector<int>> successors = int[vertices.size()][references[0].size()];
        for (const std::vector<int> s : successors) 
        {
            Arrays.fill(s, -1);
        }

        for (const int& v = 0; v < vertices.size(); ++v) 
        {
            for (int k{}; k < successors[v].size() && references[v][k] >= 0; ++k) 
            {

                // look for vertex v
                const std::vector<int> facet = facets.get(references[v][k]);
                int i = 0;
                while (i < facet.size() && facet[i] != v) 
                {
                    ++i;
                }

                // we have found vertex v, we deduce its successor on current facet
                successors[v][k] = facet[(i + 1) % facet.size()];
                for (const int& l = 0; l < k; ++l) 
                {
                    if (successors[v][l] == successors[v][k]) 
                    {
                        const Vector_3D start = vertices.get(v);
                        const Vector_3D end   = vertices.get(successors[v][k]);
                        throw (Localized_Geometry_Formats.FACET_ORIENTATION_MISMATCH, start.get_x(), start.get_y(), start.get_z(), end.get_x(),   end.get_y(),   end.get_z());
                    }
                }

            }
        }

        return successors;

    }

    /** {@inherit_doc} */
    //override
    public Polyhedrons_Set build_new(const BSP_Tree<Euclidean_3D> tree) 
    {
        return Polyhedrons_Set(tree, get_tolerance());
    }

    /** Get the boundary representation of the instance.
     * <p>
     * The boundary representation can be extracted <em>only</em> from
     * bounded polyhedrons sets. If the polyhedrons set is unbounded, * a {@link Math_Runtime_Exception} will be thrown.
     * </p>
     * <p>
     * The boundary representation extracted is not minimal, as for
     * example canonical facets may be split into several smaller
     * independent sub-facets sharing the same plane and connected by
     * their edges.
     * </p>
     * <p>
     * As the {@link B_Rep B-Rep} representation does not support
     * facets with several boundary loops (for example facets with
     * holes), an exception is triggered when attempting to extract
     * B-Rep from such complex polyhedrons sets.
     * </p>
     * @return boundary representation of the instance
     * @exception Math_Runtime_Exception if polyhedrons is unbounded
     * @since 1.2
     */
    public B_Rep get_b_rep() Math_Runtime_Exception 
    {
        B_RepExtractor extractor = B_RepExtractor(get_tolerance());
        get_tree(true).visit(extractor);
        return extractor.get_b_rep();
    }

    /** Visitor extracting B_Rep. */
    private static class B_RepExtractor : BSP_Tree_Visitor<Euclidean_3D> 
    {

        /** Tolerance for vertices identification. */
        private const double& tolerance;

        /** Extracted vertices. */
        private const List<Vector_3D> vertices;

        /** Extracted facets. */
        private const List<std::vector<int>> facets;

        /** Simple constructor.
         * @param tolerance tolerance for vertices identification
         */
        B_RepExtractor(const double& tolerance) 
        {
            this.tolerance = tolerance;
            this.vertices  = Array_list<>();
            this.facets    = Array_list<>();
        }

        /** Get the B_Rep.
         * @return extracted B_Rep
         */
        public B_Rep get_b_rep() 
        {
            return B_Rep(vertices, facets);
        }

        /** {@inherit_doc} */
        //override
        public Order visit_order(const BSP_Tree<Euclidean_3D> node) 
        {
            return Order.MINUS_SUB_PLUS;
        }

        /** {@inherit_doc} */
        //override
        public void visit_internal_node(const BSP_Tree<Euclidean_3D> node) 
        {
            //@Suppress_Warnings("unchecked")
            const Boundary_Attribute<Euclidean_3D> attribute =
                (Boundary_Attribute<Euclidean_3D>) node.get_attribute();
            if (attribute.get_plus_outside() != NULL) 
            {
                add_contribution(attribute.get_plus_outside(), false);
            }
            if (attribute.get_plus_inside() != NULL) 
            {
                add_contribution(attribute.get_plus_inside(), true);
            }
        }

        /** {@inherit_doc} */
        //override
        public void visit_leaf_node(const BSP_Tree<Euclidean_3D> node) 
        {
        }

        /** Add he contribution of a boundary facet.
         * @param facet boundary facet
         * @param reversed if true, the facet has the inside on its plus side
         * @exception Math_Runtime_Exception if facet is unbounded
         */
        private void add_contribution(const Sub_Hyperplane<Euclidean_3D> facet, const bool reversed)
            Math_Runtime_Exception 
            {

            const Plane plane = (Plane) facet.get_hyperplane();
            const Polygons_Set polygon = (Polygons_Set) ((Sub_Plane) facet).get_remaining_region();
            const Vector_2D[][] loops_2d = polygon.get_vertices();
            if (loops_2d.size() == 0) 
            {
                throw Math_Runtime_Exception(Localized_Geometry_Formats.OUTLINE_BOUNDARY_LOOP_OPEN);
            }
else if (loops_2d.size() > 1) 
            {
                throw Math_Runtime_Exception(Localized_Geometry_Formats.FACET_WITH_SEVERAL_BOUNDARY_LOOPS);
            }
else 
            {
                for (const Vector_2D[] loop_2d : polygon.get_vertices()) 
                {
                    const std::vector<int> loop_3d = int[loop_2d.size()];
                    for (int i{}; i < loop_2d.size() ; ++i) 
                    {
                        if (loop_2d[i] == NULL) 
                        {
                            throw Math_Runtime_Exception(Localized_Geometry_Formats.OUTLINE_BOUNDARY_LOOP_OPEN);
                        }
                        loop_3d[reversed ? loop_2d.size() - 1 - i : i] = get_vertex_index(plane.to_space(loop_2d[i]));
                    }
                    facets.add(loop_3d);
                }
            }

        }

        /** Get the index of a vertex.
         * @param vertex vertex as a 3D point
         * @return index of the vertex
         */
        private int get_vertex_index(const Vector_3D vertex) 
        {

            for (int i{}; i < vertices.size(); ++i) 
            {
                if (Vector_3D.distance(vertex, vertices.get(i)) <= tolerance) 
                {
                    // the vertex is already known
                    return i;
                }
            }

            // the vertex is a one, add it
            vertices.add(vertex);
            return vertices.size() - 1;

        }

    }

    /** {@inherit_doc} */
    //override
    protected void compute_geometrical_properties() 
    {

        // compute the contribution of all boundary facets
        get_tree(true).visit(new Facets_Contribution_Visitor());

        if (get_size() < 0) 
        {
            // the polyhedrons set as a finite outside
            // surrounded by an infinite inside
            set_size(INFINITY);
            set_barycenter((Point<Euclidean_3D>) Vector_3D.NaN);
        }
else 
        {
            // the polyhedrons set is finite, apply the remaining scaling factors
            set_size(get_size() / 3.0);
            set_barycenter((Point<Euclidean_3D>) Vector_3D(1.0 / (4 * get_size()), (Vector_3D) get_barycenter()));
        }

    }

    /** Visitor computing geometrical properties. */
    private class Facets_Contribution_Visitor : BSP_Tree_Visitor<Euclidean_3D> 
    {

        /** Simple constructor. */
        Facets_Contribution_Visitor() 
        {
            set_size(0);
            set_barycenter((Point<Euclidean_3D>) Vector_3D(0, 0, 0));
        }

        /** {@inherit_doc} */
        //override
        public Order visit_order(const BSP_Tree<Euclidean_3D> node) 
        {
            return Order.MINUS_SUB_PLUS;
        }

        /** {@inherit_doc} */
        //override
        public void visit_internal_node(const BSP_Tree<Euclidean_3D> node) 
        {
            //@Suppress_Warnings("unchecked")
            const Boundary_Attribute<Euclidean_3D> attribute =
                (Boundary_Attribute<Euclidean_3D>) node.get_attribute();
            if (attribute.get_plus_outside() != NULL) 
            {
                add_contribution(attribute.get_plus_outside(), false);
            }
            if (attribute.get_plus_inside() != NULL) 
            {
                add_contribution(attribute.get_plus_inside(), true);
            }
        }

        /** {@inherit_doc} */
        //override
        public void visit_leaf_node(const BSP_Tree<Euclidean_3D> node) 
        {
        }

        /** Add he contribution of a boundary facet.
         * @param facet boundary facet
         * @param reversed if true, the facet has the inside on its plus side
         */
        private void add_contribution(const Sub_Hyperplane<Euclidean_3D> facet, const bool reversed) 
        {

            const Region<Euclidean_2D> polygon = ((Sub_Plane) facet).get_remaining_region();
            const double& area    = polygon.get_size();

            if (std::isinf(area)) 
            {
                set_size(INFINITY);
                set_barycenter((Point<Euclidean_3D>) Vector_3D.NaN);
            }
else 
            {

                const Plane    plane  = (Plane) facet.get_hyperplane();
                const Vector_3D facet_b = plane.to_space(polygon.get_barycenter());
                double   scaled = area * facet_b.dot_product(plane.get_normal());
                if (reversed) 
                {
                    scaled = -scaled;
                }

                set_size(get_size() + scaled);
                set_barycenter((Point<Euclidean_3D>) Vector_3D(1.0, (Vector_3D) get_barycenter(), scaled, facet_b));

            }

        }

    }

    /** Get the first sub-hyperplane crossed by a semi-infinite line.
     * @param point start point of the part of the line considered
     * @param line line to consider (contains point)
     * @return the first sub-hyperplane crossed by the line after the
     * given point, or NULL if the line does not intersect any
     * sub-hyperplane
     */
    public Sub_Hyperplane<Euclidean_3D> first_intersection(const Vector_3D point, const Line& line) 
    {
        return recurse_first_intersection(get_tree(true), point, line);
    }

    /** Get the first sub-hyperplane crossed by a semi-infinite line.
     * @param node current node
     * @param point start point of the part of the line considered
     * @param line line to consider (contains point)
     * @return the first sub-hyperplane crossed by the line after the
     * given point, or NULL if the line does not intersect any
     * sub-hyperplane
     */
    private Sub_Hyperplane<Euclidean_3D> recurse_first_intersection(const BSP_Tree<Euclidean_3D> node, const Vector_3D point, const Line& line) 
    {

        const Sub_Hyperplane<Euclidean_3D> cut = node.get_cut();
        if (cut == NULL) 
        {
            return NULL;
        }
        const BSP_Tree<Euclidean_3D> minus = node.get_minus();
        const BSP_Tree<Euclidean_3D> plus  = node.get_plus();
        const Plane                plane = (Plane) cut.get_hyperplane();

        // establish search order
        const double offset = plane.get_offset((Point<Euclidean_3D>) point);
        const bool in    = std::abs(offset) < get_tolerance();
        const BSP_Tree<Euclidean_3D> near;
        const BSP_Tree<Euclidean_3D> far;
        if (offset < 0) 
        {
            near = minus;
            far  = plus;
        }
else 
        {
            near = plus;
            far  = minus;
        }

        if (in) 
        {
            // search in the cut hyperplane
            const Sub_Hyperplane<Euclidean_3D> facet = boundary_facet(point, node);
            if (facet != NULL) 
            {
                return facet;
            }
        }

        // search in the near branch
        const Sub_Hyperplane<Euclidean_3D> crossed = recurse_first_intersection(near, point, line);
        if (crossed != NULL) 
        {
            return crossed;
        }

        if (!in) 
        {
            // search in the cut hyperplane
            const Vector_3D hit_3d = plane.intersection(line);
            if (hit_3d != NULL && line.get_abscissa(hit_3d) > line.get_abscissa(point)) 
            {
                const Sub_Hyperplane<Euclidean_3D> facet = boundary_facet(hit_3d, node);
                if (facet != NULL) 
                {
                    return facet;
                }
            }
        }

        // search in the far branch
        return recurse_first_intersection(far, point, line);

    }

    /** Check if a point belongs to the boundary part of a node.
     * @param point point to check
     * @param node node containing the boundary facet to check
     * @return the boundary facet this points belongs to (or NULL if it
     * does not belong to any boundary facet)
     */
    private Sub_Hyperplane<Euclidean_3D> boundary_facet(const Vector_3D point, const BSP_Tree<Euclidean_3D> node) 
    {
        const Vector_2D point_2d = ((Plane) node.get_cut().get_hyperplane()).to_sub_space((Point<Euclidean_3D>) point);
        //@Suppress_Warnings("unchecked")
        const Boundary_Attribute<Euclidean_3D> attribute =
            (Boundary_Attribute<Euclidean_3D>) node.get_attribute();
        if ((attribute.get_plus_outside() != NULL) &&
            (((Sub_Plane) attribute.get_plus_outside()).get_remaining_region().check_point(point_2d) != Location.OUTSIDE)) 
            {
            return attribute.get_plus_outside();
        }
        if ((attribute.get_plus_inside() != NULL) &&
            (((Sub_Plane) attribute.get_plus_inside()).get_remaining_region().check_point(point_2d) != Location.OUTSIDE)) 
            {
            return attribute.get_plus_inside();
        }
        return NULL;
    }

    /** Rotate the region around the specified point.
     * <p>The instance is not modified, a instance is created.</p>
     * @param center rotation center
     * @param rotation vectorial rotation operator
     * @return a instance representing the rotated region
     */
    public Polyhedrons_Set rotate(const Vector_3D center, const Rotation& rotation) 
    {
        return (Polyhedrons_Set) apply_transform(new Rotation_Transform(center, rotation));
    }

    /** 3D rotation as a Transform. */
    private static class Rotation_Transform : Transform<Euclidean_3D, Euclidean_2D> 
    {

        /** Center point of the rotation. */
        private const Vector_3D   center;

        /** Vectorial rotation. */
        private const Rotation   rotation;

        /** Cached original hyperplane. */
        private Plane cached_original;

        /** Cached 2D transform valid inside the cached original hyperplane. */
        private Transform<Euclidean_2D, Euclidean_1D>  cached_transform;

        /** Build a rotation transform.
         * @param center center point of the rotation
         * @param rotation vectorial rotation
         */
        Rotation_Transform(const Vector_3D center, const Rotation& rotation) 
        {
            this.center   = center;
            this.rotation = rotation;
        }

        /** {@inherit_doc} */
        //override
        public Vector_3D apply(const Point<Euclidean_3D> point) 
        {
            const Vector_3D delta = ((Vector_3D) point).subtract(center);
            return Vector_3D(1.0, center, 1.0, rotation.apply_to(delta));
        }

        /** {@inherit_doc} */
        //override
        public Plane apply(const Hyperplane<Euclidean_3D> hyperplane) 
        {
            return ((Plane) hyperplane).rotate(center, rotation);
        }

        /** {@inherit_doc} */
        //override
        public Sub_Hyperplane<Euclidean_2D> apply(const Sub_Hyperplane<Euclidean_2D> sub, const Hyperplane<Euclidean_3D> original, const Hyperplane<Euclidean_3D> transformed) 
        {
            if (original != cached_original) 
            {
                // we have changed hyperplane, reset the in-hyperplane transform

                const Plane    o_plane = (Plane) original;
                const Plane    t_plane = (Plane) transformed;
                const Vector_3D p00    = o_plane.get_origin();
                const Vector_3D p10    = o_plane.to_space((Point<Euclidean_2D>) Vector_2D(1.0, 0.0));
                const Vector_3D p01    = o_plane.to_space((Point<Euclidean_2D>) Vector_2D(0.0, 1.0));
                const Vector_2D t_p00   = t_plane.to_sub_space((Point<Euclidean_3D>) apply(p00));
                const Vector_2D tP10   = t_plane.to_sub_space((Point<Euclidean_3D>) apply(p10));
                const Vector_2D tP01   = t_plane.to_sub_space((Point<Euclidean_3D>) apply(p01));

                cached_original  = (Plane) original;
                cached_transform =
                        org.hipparchus.geometry.euclidean.twod.Line.get_transform(tP10.get_x() - t_p00.get_x(), tP10.get_y() - t_p00.get_y(), tP01.get_x() - t_p00.get_x(), tP01.get_y() - t_p00.get_y(), t_p00.get_x(), t_p00.get_y());

            }
            return ((Sub_Line) sub).apply_transform(cached_transform);
        }

    }

    /** Translate the region by the specified amount.
     * <p>The instance is not modified, a instance is created.</p>
     * @param translation translation to apply
     * @return a instance representing the translated region
     */
    public Polyhedrons_Set translate(const Vector_3D translation) 
    {
        return (Polyhedrons_Set) apply_transform(new Translation_Transform(translation));
    }

    /** 3D translation as a transform. */
    private static class Translation_Transform : Transform<Euclidean_3D, Euclidean_2D> 
    {

        /** Translation vector. */
        private const Vector_3D   translation;

        /** Cached original hyperplane. */
        private Plane cached_original;

        /** Cached 2D transform valid inside the cached original hyperplane. */
        private Transform<Euclidean_2D, Euclidean_1D>  cached_transform;

        /** Build a translation transform.
         * @param translation translation vector
         */
        Translation_Transform(const Vector_3D translation) 
        {
            this.translation = translation;
        }

        /** {@inherit_doc} */
        //override
        public Vector_3D apply(const Point<Euclidean_3D> point) 
        {
            return Vector_3D(1.0, (Vector_3D) point, 1.0, translation);
        }

        /** {@inherit_doc} */
        //override
        public Plane apply(const Hyperplane<Euclidean_3D> hyperplane) 
        {
            return ((Plane) hyperplane).translate(translation);
        }

        /** {@inherit_doc} */
        //override
        public Sub_Hyperplane<Euclidean_2D> apply(const Sub_Hyperplane<Euclidean_2D> sub, const Hyperplane<Euclidean_3D> original, const Hyperplane<Euclidean_3D> transformed) 
        {
            if (original != cached_original) 
            {
                // we have changed hyperplane, reset the in-hyperplane transform

                const Plane   o_plane = (Plane) original;
                const Plane   t_plane = (Plane) transformed;
                const Vector_2D shift  = t_plane.to_sub_space((Point<Euclidean_3D>) apply(o_plane.get_origin()));

                cached_original  = (Plane) original;
                cached_transform =
                        org.hipparchus.geometry.euclidean.twod.Line.get_transform(1, 0, 0, 1, shift.get_x(), shift.get_y());

            }

            return ((Sub_Line) sub).apply_transform(cached_transform);

        }

    }

    /** Container for Boundary RE_Presentation (B-Rep).
     * <p>
     * The boundary is provided as a list of vertices and a list of facets.
     * Each facet is specified as an integer array containing the arrays vertices
     * indices in the vertices list. Each facet normal is oriented by right hand
     * rule to the facet vertices list.
     * </p>
     * @see Polyhedrons_Set#Polyhedrons_Set(BSP_Tree, double)
     * @see Polyhedrons_Set#get_b_rep()
     * @since 1.2
     */
    public static class B_Rep 
    {

        /** List of polyhedrons set vertices. */
        private const List<Vector_3D> vertices;

        /** List of facets, as vertices indices in the vertices list. */
        private const List<std::vector<int>> facets;

        /** Simple constructor.
         * @param vertices list of polyhedrons set vertices
         * @param facets list of facets, as vertices indices in the vertices list
         */
        public B_Rep(const List<Vector_3D> vertices, const List<std::vector<int>> facets) 
        {
            this.vertices = vertices;
            this.facets   = facets;
        }

        /** Get the extracted vertices.
         * @return extracted vertices
         */
        public List<Vector_3D> get_vertices() 
        {
            return vertices;
        }

        /** Get the extracted facets.
         * @return extracted facets
         */
        public List<std::vector<int>> get_facets() 
        {
            return facets;
        }

    }

}


