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

  //import java.util.Event_Listener;
#include "IterationEvent.h"

/**
 * The listener interface for receiving events occurring in an iterative
 * algorithm.
 */
class Iteration_Listener// : public Event_Listener
{
	/**
	 * Invoked after completion of the initial phase of the iterative algorithm
	 * (prior to the main iteration loop).
	 *
	 * @param e The {@link Iteration_Event} object.
	 */
	virtual void initialization_performed(Iteration_Event e) = 0;

	/**
	 * Invoked each time an iteration is completed (in the main iteration loop).
	 *
	 * @param e The {@link Iteration_Event} object.
	 */
	virtual void iteration_performed(Iteration_Event e) = 0;

	/**
	 * Invoked each time a iteration is completed (in the main iteration
	 * loop).
	 *
	 * @param e The {@link Iteration_Event} object.
	 */
	virtual void iteration_started(Iteration_Event e) = 0;

	/**
	 * Invoked after completion of the operations which occur after breaking out
	 * of the main iteration loop.
	 *
	 * @param e The {@link Iteration_Event} object.
	 */
	virtual void termination_performed(Iteration_Event e) = 0;
};