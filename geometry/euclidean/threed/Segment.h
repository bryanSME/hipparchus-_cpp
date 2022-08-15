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


/** Simple container for a two-points segment.
 */
class Segment 
{

    /** Start point of the segment. */
    private const Vector_3D start;

    /** End point of the segments. */
    private const Vector_3D end;

    /** Line containing the segment. */
    private const Line     line;

    /** Build a segment.
     * @param start start point of the segment
     * @param end end point of the segment
     * @param line line containing the segment
     */
    public Segment(const Vector_3D start, const Vector_3D end, const Line& line) 
    {
        this.start  = start;
        this.end    = end;
        this.line   = line;
    }

    /** Get the start point of the segment.
     * @return start point of the segment
     */
    public Vector_3D get_start() 
    {
        return start;
    }

    /** Get the end point of the segment.
     * @return end point of the segment
     */
    public Vector_3D get_end() 
    {
        return end;
    }

    /** Get the line containing the segment.
     * @return line containing the segment
     */
    public Line get_line() const 
    {
        return line;
    }

}


