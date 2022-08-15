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
//package org.hipparchus.geometry.partitioning;
#include "../Space.h"
#include <type_traits>

//import org.hipparchus.geometry.Point;
//import org.hipparchus.geometry.Space;

/** Class holding the result of point projection on region boundary.
 * <p>This class is a simple placeholder, it does not provide any
 * processing methods.</p>
 * <p>Instances of this class are guaranteed to be immutable</p>
 * @param <S> Type of the space.
 * @see Abstract_Region#project_to_boundary(Point)
 */
template<typename S, typename std::enable_if<std::is_base_of<Space, S>::value>::type* = nullptr>
class Boundary_Projection
{

private:
    /** Original point. */
    const Point<S> my_original;

    /** Projected point. */
    const Point<S> my_projected;

    /** Offset of the point with respect to the boundary it is projected on. */
    const double my_offset;

public:
    /** Constructor from raw elements.
     * @param original original point
     * @param projected projected point
     * @param offset offset of the point with respect to the boundary it is projected on
     */
    Boundary_Projection(const Point<S>& original, const Point<S>& projected, const double& offset)
        : my_original{ original }, my_projected{ projected }, my_offset{ offset }
    {};

    /** Get the original point.
     * @return original point
     */
    Point<S> get_original() const
    {
        return my_original;
    }

    /** Projected point.
     * @return projected point, or NULL if there are no boundary
     */
    Point<S> get_projected() const
    {
        return my_projected;
    }

    /** Offset of the point with respect to the boundary it is projected on.
     * <p>
     * The offset with respect to the boundary is negative if the {@link
     * #get_original() original point} is inside the region, and positive otherwise.
     * </p>
     * <p>
     * If there are no boundary, the value is set to either {@code
     * INFINITY} if the region is empty (i.e. all points are
     * outside of the region) or {@code -INFINITY} if the region
     * covers the whole space (i.e. all points are inside of the region).
     * </p>
     * @return offset of the point with respect to the boundary it is projected on
     */
    double get_offset() const
    {
        return my_offset;
    }

};