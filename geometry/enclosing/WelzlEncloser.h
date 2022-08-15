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
//package org.hipparchus.geometry.enclosing;

//import java.util.Array_list;
//import java.util.List;

//import org.hipparchus.exception.Math_Runtime_Exception;
//import org.hipparchus.geometry.Point;
//import org.hipparchus.geometry.Space;

/** Class implementing Emo Welzl algorithm to find the smallest enclosing ball in linear time.
 * <p>
 * The class : the algorithm described in paper <a
 * href="http://www.inf.ethz.ch/personal/emo/PublFiles/SmallEnclDisk_LNCS555_91.pdf">Smallest
 * Enclosing Disks (Balls and Ellipsoids)</a> by Emo Welzl, Lecture Notes in Computer Science
 * 555 (1991) 359-370. The pivoting improvement published in the paper <a
 * href="http://www.inf.ethz.ch/personal/gaertner/texts/own_work/esa99_const.pdf">Fast and
 * Robust Smallest Enclosing Balls</a>, by Bernd Gärtner and further modified in
 * paper <a href="http://www.idt.mdh.se/kurser/ct3340/ht12/MINICONFERENCE/FinalPapers/ircse12_submission_30.pdf">
 * Efficient Computation of Smallest Enclosing Balls in Three Dimensions</a> by Linus Källberg
 * to avoid performing local copies of data have been included.
 * </p>
 * @param <S> Space type.
 * @param <P> Point type.
 */
class Welzl_Encloser<S extends Space, P extends Point<S>> : Encloser<S, P> 
{

    /** Tolerance below which points are consider to be identical. */
    private const double& tolerance;

    /** Generator for balls on support. */
    private const Support_Ball_Generator<S, P> generator;

    /** Simple constructor.
     * @param tolerance below which points are consider to be identical
     * @param generator generator for balls on support
     */
    public Welzl_Encloser(const double& tolerance, const Support_Ball_Generator<S, P> generator) 
    {
        this.tolerance = tolerance;
        this.generator = generator;
    }

    /** {@inherit_doc} */
    //override
    public Enclosing_Ball<S, P> enclose(const Iterable<P> points) 
    {

        if (points == NULL || !points.iterator().has_next()) 
        {
            // return an empty ball
            return generator.ball_on_support(new Array_list<P>());
        }

        // Emo Welzl algorithm with Bernd Gärtner and Linus Källberg improvements
        return pivoting_ball(points);

    }

    /** Compute enclosing ball using Gärtner's pivoting heuristic.
     * @param points points to be enclosed
     * @return enclosing ball
     */
    private Enclosing_Ball<S, P> pivoting_ball(const Iterable<P> points) 
    {

        const P first = points.iterator().next();
        const List<P> extreme = Array_list<>(first.get_space().get_dimension() + 1);
        const List<P> support = Array_list<>(first.get_space().get_dimension() + 1);

        // start with only first point selected as a candidate support
        extreme.add(first);
        Enclosing_Ball<S, P> ball = move_to_front_ball(extreme, extreme.size(), support);

        while (true) 
        {

            // select the point farthest to current ball
            const P farthest = select_farthest(points, ball);

            if (ball.contains(farthest, tolerance)) 
            {
                // we have found a ball containing all points
                return ball;
            }

            // recurse search, restricted to the small subset containing support and farthest point
            support.clear();
            support.add(farthest);
            Enclosing_Ball<S, P> saved_ball = ball;
            ball = move_to_front_ball(extreme, extreme.size(), support);
            if (ball.get_radius() < saved_ball.get_radius()) 
            {
                // this should never happen
                throw Math_Runtime_Exception.create_internal_error();
            }

            // it was an interesting point, move it to the front
            // according to Gärtner's heuristic
            extreme.add(0, farthest);

            // prune the least interesting points
            extreme.sub_list(ball.get_support_size(), extreme.size()).clear();


        }
    }

    /** Compute enclosing ball using Welzl's move to front heuristic.
     * @param extreme subset of extreme points
     * @param nb_extreme number of extreme points to consider
     * @param support points that must belong to the ball support
     * @return enclosing ball, for the extreme subset only
     */
    private Enclosing_Ball<S, P> move_to_front_ball(const List<P> extreme, const int& nb_extreme, const List<P> support) 
    {

        // create a ball on the prescribed support
        Enclosing_Ball<S, P> ball = generator.ball_on_support(support);

        if (ball.get_support_size() <= ball.get_center().get_space().get_dimension()) 
        {

            for (int i{}; i < nb_extreme; ++i) 
            {
                const P pi = extreme.get(i);
                if (!ball.contains(pi, tolerance)) 
                {

                    // we have found an outside point, // enlarge the ball by adding it to the support
                    support.add(pi);
                    ball = move_to_front_ball(extreme, i, support);
                    support.remove(support.size() - 1);

                    // it was an interesting point, move it to the front
                    // according to Welzl's heuristic
                    for (int j = i; j > 0; --j) 
                    {
                        extreme.set(j, extreme.get(j - 1));
                    }
                    extreme.set(0, pi);

                }
            }

        }

        return ball;

    }

    /** Select the point farthest to the current ball.
     * @param points points to be enclosed
     * @param ball current ball
     * @return farthest point
     */
    public P select_farthest(const Iterable<P> points, const Enclosing_Ball<S, P> ball) 
    {

        const P center = ball.get_center();
        P farthest   = NULL;
        double d_max  = -1.0;

        for (const P point : points) 
        {
            const double d = point.distance(center);
            if (d > d_max) 
            {
                farthest = point;
                d_max     = d;
            }
        }

        return farthest;

    }

}


