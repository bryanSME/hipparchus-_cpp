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
  //package org.hipparchus.geometry.euclidean.twod;

  //import org.hipparchus.geometry.Point;
  //import org.hipparchus.util.FastMath;

  /** Simple container for a two-points segment.
   */
class Segment
{
	/** Start point of the segment. */
	private const Vector_2D start;

	/** End point of the segment. */
	private const Vector_2D end;

	/** Line containing the segment. */
	private const Line     line;

	/** Build a segment.
	 * @param start start point of the segment
	 * @param end end point of the segment
	 * @param tolerance of the line.
	 */
	public Segment(const Vector_2D& start, const Vector_2D& end, const double& tolerance)
	{
		this(start, end, Line(start, end, tolerance));
	}

	/** Build a segment.
	 * @param start start point of the segment
	 * @param end end point of the segment
	 * @param line line containing the segment
	 */
	public Segment(const Vector_2D& start, const Vector_2D& end, const Line& line)
	{
		this.start = start;
		this.end = end;
		this.line = line;
	}

	/** Get the start point of the segment.
	 * @return start point of the segment
	 */
	public Vector_2D get_start()
	{
		return start;
	}

	/** Get the end point of the segment.
	 * @return end point of the segment
	 */
	public Vector_2D get_end() const
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

	/**
	 * Get the length of the line segment.
	 *
	 * @return line segment length.
	 */
	public double get_length()
	{
		return get_end().distance(get_start());
	}

	/** Calculates the shortest distance from a point to this line segment.
	 * <p>
	 * If the perpendicular extension from the point to the line does not
	 * cross in the bounds of the line segment, the shortest distance to
	 * the two end points will be returned.
	 * </p>
	 *
	 * Algorithm adapted from:
	 * <a href="http://www.codeguru.com/forum/printthread.php?s=cc8cf0596231f9a7dba4da6e77c29db3&amp;t=194400&amp;pp=15&amp;page=1">
	 * Thread @ Codeguru</a>
	 *
	 * @param p to check
	 * @return distance between the instance and the point
	 */
	public double distance(const Vector_2D& p)
	{
		const double delta_x = end.get_x() - start.get_x();
		const double delta_y = end.get_y() - start.get_y();

		const double r = ((p.get_x() - start.get_x()) * delta_x + (p.get_y() - start.get_y()) * delta_y) /
			(delta_x * delta_x + delta_y * delta_y);

		// r == 0 => P = start_pt
		// r == 1 => P = end_pt
		// r < 0 => P is on the backward extension of the segment
		// r > 1 => P is on the forward extension of the segment
		// 0 < r < 1 => P is on the segment

		// if point isn't on the line segment, just return the shortest distance to the end points
		if (r < 0 || r > 1)
		{
			const double dist1 = get_start().distance((Point<Euclidean_2D>) p);
			const double dist2 = get_end().distance((Point<Euclidean_2D>) p);

			return std::min(dist1, dist2);
		}
		else
		{
			// find point on line and see if it is in the line segment
			const double px = start.get_x() + r * delta_x;
			const double py = start.get_y() + r * delta_y;

			const Vector_2D inter_pt = Vector_2D(px, py);
			return inter_pt.distance((Point<Euclidean_2D>) p);
		}
	}
}
