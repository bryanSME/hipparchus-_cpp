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
//import org.hipparchus.geometry.Space;

/** Class holding boundary attributes.
 * <p>This class is used for the attributes associated with the
 * nodes of region boundary shell trees returned by the {@link
 * Region#get_tree(bool) Region.get_tree(include_boundary_attributes)}
 * when the bool {@code include_boundary_attributes} parameter is
 * set to {@code true}. It contains the parts of the node cut
 * sub-hyperplane that belong to the boundary.</p>
 * <p>This class is a simple placeholder, it does not provide any
 * processing methods.</p>
 * @param <S> Type of the space.
 * @see Region#get_tree
 */
template<typename S, typename std::enable_if<std::is_base_of<Space, S>::value>::type* = nullptr>
class Boundary_Attribute
{

    /** Part of the node cut sub-hyperplane that belongs to the
     * boundary and has the outside of the region on the plus side of
     * its underlying hyperplane (may be NULL).
     */
    private const Sub_Hyperplane<S> plus_outside;

    /** Part of the node cut sub-hyperplane that belongs to the
     * boundary and has the inside of the region on the plus side of
     * its underlying hyperplane (may be NULL).
     */
    private const Sub_Hyperplane<S> plus_inside;

    /** Sub-hyperplanes that were used to split the boundary part. */
    private const Nodes_Set<S> splitters;

    /** Simple constructor.
     * @param plus_outside part of the node cut sub-hyperplane that
     * belongs to the boundary and has the outside of the region on
     * the plus side of its underlying hyperplane (may be NULL)
     * @param plus_inside part of the node cut sub-hyperplane that
     * belongs to the boundary and has the inside of the region on the
     * plus side of its underlying hyperplane (may be NULL)
     * @param splitters sub-hyperplanes that were used to
     * split the boundary part (may be NULL)
     */
    Boundary_Attribute(const Sub_Hyperplane<S> plus_outside, const Sub_Hyperplane<S> plus_inside, const Nodes_Set<S> splitters) 
    {
        this.plus_outside = plus_outside;
        this.plus_inside  = plus_inside;
        this.splitters   = splitters;
    }

    /** Get the part of the node cut sub-hyperplane that belongs to the
     * boundary and has the outside of the region on the plus side of
     * its underlying hyperplane.
     * @return part of the node cut sub-hyperplane that belongs to the
     * boundary and has the outside of the region on the plus side of
     * its underlying hyperplane
     */
    public Sub_Hyperplane<S> get_plus_outside() 
    {
        return plus_outside;
    }

    /** Get the part of the node cut sub-hyperplane that belongs to the
     * boundary and has the inside of the region on the plus side of
     * its underlying hyperplane.
     * @return part of the node cut sub-hyperplane that belongs to the
     * boundary and has the inside of the region on the plus side of
     * its underlying hyperplane
     */
    public Sub_Hyperplane<S> get_plus_inside() 
    {
        return plus_inside;
    }

    /** Get the sub-hyperplanes that were used to split the boundary part.
     * @return sub-hyperplanes that were used to split the boundary part
     */
    public Nodes_Set<S> get_splitters() 
    {
        return splitters;
    }

}


