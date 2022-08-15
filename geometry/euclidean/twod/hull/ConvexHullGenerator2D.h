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

#include <vector>
#include "../../../hull/ConvexHullGenerator.h"
#include "ConvexHull2D.h"
#include "../../twod/Euclidean2D.h"
#include "../../twod/Vector2D.h"

/**
 * Interface for convex hull generators in the two-dimensional euclidean space.
 *
 */
class Convex_Hull_Generator_2D : public Convex_Hull_Generator<Euclidean_2D, Vector_2D>
{
    /** {@inherit_doc} */
    //override
    Convex_Hull_2D generate(std::vector<Vector_2D> points);

};