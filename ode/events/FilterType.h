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

//package org.hipparchus.ode.events;

//import org.hipparchus.exception.Math_Runtime_Exception;

/** Enumerate for {@link Event_Filter filtering events}.
 *
 */

enum Filter_Type 
{

    /** Constant for triggering only decreasing events.
     * <p>When this filter is used, the wrapped {@link ODE_Event_Handler
     * event handler} {@link ODE_Event_Handler#event_occurred(org.hipparchus.ode.ODE_State_And_Derivative, * bool) event_occurred} method will be called <em>only</em> with
     * its {@code increasing} argument set to false.</p>
     */
    TRIGGER_ONLY_DECREASING_EVENTS 
    {
    protected:
        /**  {@inherit_doc} */
        //override
        bool is_triggered_on_increasing() const
        {
            return false;
        }

        /** {@inherit_doc}
         * <p>
         * states scheduling for computing h(t,y) as an altered version of g(t, y)
         * <ul>
         * <li>0 are triggered events for which a zero is produced (here decreasing events)</li>
         * <li>X are ignored events for which zero is masked (here increasing events)</li>
         * </ul>
         * </p>
         * <pre>
         *  g(t)
         *             ___                     ___                     ___
         *            /   \                   /   \                   /   \
         *           /     \                 /     \                 /     \
         *          /  g>0  \               /  g>0  \               /  g>0  \
         *         /         \             /         \             /         \
         *  ----- X --------- 0 --------- X --------- 0 --------- X --------- 0 ---
         *       /             \         /             \         /             \
         *      /               \ g<0   /               \  g<0  /               \ g<0
         *     /                 \     /                 \     /                 \     /
         * ___/                   \___/                   \___/                   \___/
         * </pre>
         * <pre>
         *  h(t,y)) as an alteration of g(t,y)
         *             ___                                 ___         ___
         *    \       /   \                               /   \       /   \
         *     \     /     \ h=+g                        /     \     /     \
         *      \   /       \      h=min(-s,-g,+g)      /       \   /       \
         *       \_/         \                         /         \_/         \
         *  ------ ---------- 0 ----------_---------- 0 --------------------- 0 ---
         *                     \         / \         /                         \
         *   h=max(+s,-g,+g)    \       /   \       /       h=max(+s,-g,+g)     \
         *                       \     /     \     / h=-g                        \     /
         *                        \___/       \___/                               \___/
         * </pre>
         * <p>
         * As shown by the figure above, several expressions are used to compute h, * depending on the current state:
         * <ul>
         *   <li>h = max(+s,-g,+g)</li>
         *   <li>h = +g</li>
         *   <li>h = min(-s,-g,+g)</li>
         *   <li>h = -g</li>
         * </ul>
         * where s is a tiny positive value: {@link org.hipparchus.util.Precision#SAFE_MIN}.
         * </p>
         */
        //override
        Transformer select_transformer(const Transformer previous, const double g, const bool forward) 
        {
            if (forward) 
            {
                switch (previous) 
                {
                    case UNINITIALIZED :
                        // we are initializing the first point
                        if (g > 0) 
                        {
                            // initialize as if previous root (i.e. backward one) was an ignored increasing event
                            return Transformer.MAX;
                        }
                        else if (g < 0) 
                        {
                            // initialize as if previous root (i.e. backward one) was a triggered decreasing event
                            return Transformer.PLUS;
                        }

                            // we are exactly at a root, we don't know if it is an increasing
                            // or a decreasing event, we remain in uninitialized state
                            return Transformer.UNINITIALIZED;
                        
                    case PLUS  :
                        if (g >= 0) 
                        {
                            // we have crossed the zero line on an ignored increasing event, // we must change the transformer
                            return Transformer.MIN;
                        }

                            // we are still in the same status
                            return previous;
                        
                    case MINUS :
                        if (g >= 0) 
                        {
                            // we have crossed the zero line on an ignored increasing event, // we must change the transformer
                            return Transformer.MAX;
                        }

                            // we are still in the same status
                            return previous;
                        
                    case MIN   :
                        if (g <= 0) 
                        {
                            // we have crossed the zero line on a triggered decreasing event, // we must change the transformer
                            return Transformer.MINUS;
                        }

                            // we are still in the same status
                            return previous;
                        
                    case MAX   :
                        if (g <= 0) 
                        {
                            // we have crossed the zero line on a triggered decreasing event, // we must change the transformer
                            return Transformer.PLUS;
                        }

                            // we are still in the same status
                            return previous;
                        
                    default    :
                        // this should never happen
                        throw Math_Runtime_Exception.create_internal_error();
                }
            }
            else 
            {
                switch (previous) 
                {
                    case UNINITIALIZED :
                        // we are initializing the first point
                        if (g > 0) 
                        {
                            // initialize as if previous root (i.e. forward one) was a triggered decreasing event
                            return Transformer.MINUS;
                        }
                        else if (g < 0) 
                        {
                            // initialize as if previous root (i.e. forward one) was an ignored increasing event
                            return Transformer.MIN;
                        }
  
                            // we are exactly at a root, we don't know if it is an increasing
                            // or a decreasing event, we remain in uninitialized state
                            return Transformer.UNINITIALIZED;
                        
                    case PLUS  :
                        if (g <= 0) 
                        {
                            // we have crossed the zero line on an ignored increasing event, // we must change the transformer
                            return Transformer.MAX;
                        }

                            // we are still in the same status
                            return previous;
                        
                    case MINUS :
                        if (g <= 0) 
                        {
                            // we have crossed the zero line on an ignored increasing event, // we must change the transformer
                            return Transformer.MIN;
                        }

                            // we are still in the same status
                            return previous;
                        
                    case MIN   :
                        if (g >= 0) 
                        {
                            // we have crossed the zero line on a triggered decreasing event, // we must change the transformer
                            return Transformer.PLUS;
                        }
                            // we are still in the same status
                            return previous;
                       
                    case MAX   :
                        if (g >= 0) 
                        {
                            // we have crossed the zero line on a triggered decreasing event, // we must change the transformer
                            return Transformer.MINUS;
                        }
   
                            // we are still in the same status
                            return previous;
                        
                    default    :
                        // this should never happen
                        throw Math_Runtime_Exception.create_internal_error();
                }
            }
        }

    }, 
    /** Constant for triggering only increasing events.
     * <p>When this filter is used, the wrapped {@link ODE_Event_Handler
     * event handler} {@link ODE_Event_Handler#event_occurred(org.hipparchus.ode.ODE_State_And_Derivative, * bool) event_occurred} method will be called <em>only</em> with
     * its {@code increasing} argument set to true.</p>
     */
    TRIGGER_ONLY_INCREASING_EVENTS 
    {
    protected:
        /**  {@inherit_doc} */
        //override
        bool is_triggered_on_increasing() 
        {
            return true;
        }

        /** {@inherit_doc}
         * <p>
         * states scheduling for computing h(t,y) as an altered version of g(t, y)
         * <ul>
         * <li>0 are triggered events for which a zero is produced (here increasing events)</li>
         * <li>X are ignored events for which zero is masked (here decreasing events)</li>
         * </ul>
         * </p>
         * <pre>
         *  g(t)
         *             ___                     ___                     ___
         *            /   \                   /   \                   /   \
         *           /     \                 /     \                 /     \
         *          /  g>0  \               /  g>0  \               /  g>0  \
         *         /         \             /         \             /         \
         *  ----- 0 --------- X --------- 0 --------- X --------- 0 --------- X ---
         *       /             \         /             \         /             \
         *      /               \ g<0   /               \  g<0  /               \ g<0
         *     /                 \     /                 \     /                 \     /
         * ___/                   \___/                   \___/                   \___/
         * </pre>
         * <pre>
         *  h(t,y)) as an alteration of g(t,y)
         *                                     ___         ___
         *    \                               /   \       /   \
         *     \ h=-g                        /     \     /     \ h=-g
         *      \      h=min(-s,-g,+g)      /       \   /       \      h=min(-s,-g,+g)
         *       \                         /         \_/         \
         *  ------0 ----------_---------- 0 --------------------- 0 --------- _ ---
         *         \         / \         /                         \         / \
         *          \       /   \       /       h=max(+s,-g,+g)     \       /   \
         *           \     /     \     / h=+g                        \     /     \     /
         *            \___/       \___/                               \___/       \___/
         * </pre>
         * <p>
         * As shown by the figure above, several expressions are used to compute h, * depending on the current state:
         * <ul>
         *   <li>h = max(+s,-g,+g)</li>
         *   <li>h = +g</li>
         *   <li>h = min(-s,-g,+g)</li>
         *   <li>h = -g</li>
         * </ul>
         * where s is a tiny positive value: {@link org.hipparchus.util.Precision#SAFE_MIN}.
         * </p>
         */
        //override
        Transformer select_transformer(const Transformer& previous, const double& g, const bool forward) 
        {
            if (forward) 
            {
                switch (previous) 
                {
                    case UNINITIALIZED :
                        // we are initializing the first point
                        if (g > 0) 
                        {
                            // initialize as if previous root (i.e. backward one) was a triggered increasing event
                            return Transformer.PLUS;
                        }
                        else if (g < 0) 
                        {
                            // initialize as if previous root (i.e. backward one) was an ignored decreasing event
                            return Transformer.MIN;
                        }
      
                            // we are exactly at a root, we don't know if it is an increasing
                            // or a decreasing event, we remain in uninitialized state
                            return Transformer.UNINITIALIZED;
                        
                    case PLUS  :
                        if (g <= 0) 
                        {
                            // we have crossed the zero line on an ignored decreasing event, // we must change the transformer
                            return Transformer.MAX;
                        }
                        // we are still in the same status
                        return previous;
                    case MINUS :
                        if (g <= 0) 
                        {
                            // we have crossed the zero line on an ignored decreasing event, // we must change the transformer
                            return Transformer.MIN;
                        }
                        // we are still in the same status
                        return previous;
                    case MIN   :
                        if (g >= 0) 
                        {
                            // we have crossed the zero line on a triggered increasing event, // we must change the transformer
                            return Transformer.PLUS;
                        }

                            // we are still in the same status
                            return previous;

                    case MAX   :
                        if (g >= 0) 
                        {
                            // we have crossed the zero line on a triggered increasing event, // we must change the transformer
                            return Transformer.MINUS;
                        }

                            // we are still in the same status
                            return previous;

                    default    :
                        // this should never happen
                        throw Math_Runtime_Exception.create_internal_error();
                }
            }
            else 
            {
                switch (previous) 
                {
                    case UNINITIALIZED :
                        // we are initializing the first point
                        if (g > 0) 
                        {
                            // initialize as if previous root (i.e. forward one) was an ignored decreasing event
                            return Transformer.MAX;
                        }
                        else if (g < 0) 
                        {
                            // initialize as if previous root (i.e. forward one) was a triggered increasing event
                            return Transformer.MINUS;
                        }
                            // we are exactly at a root, we don't know if it is an increasing
                            // or a decreasing event, we remain in uninitialized state
                            return Transformer.UNINITIALIZED;
                    case PLUS  :
                        if (g >= 0) 
                        {
                            // we have crossed the zero line on an ignored decreasing event, // we must change the transformer
                            return Transformer.MIN;
                        }

                            // we are still in the same status
                            return previous;
                        
                    case MINUS :
                        if (g >= 0) 
                        {
                            // we have crossed the zero line on an ignored decreasing event, // we must change the transformer
                            return Transformer.MAX;
                        }

                            // we are still in the same status
                            return previous;

                    case MIN   :
                        if (g <= 0) 
                        {
                            // we have crossed the zero line on a triggered increasing event, // we must change the transformer
                            return Transformer.MINUS;
                        }

                            // we are still in the same status
                            return previous;
                        
                    case MAX   :
                        if (g <= 0) 
                        {
                            // we have crossed the zero line on a triggered increasing event, // we must change the transformer
                            return Transformer.PLUS;
                        }

                            // we are still in the same status
                            return previous;
                        
                    default    :
                        // this should never happen
                        throw Math_Runtime_Exception.create_internal_error();
                }
            }
        }

    };

    /** Check if triggered events are increasing events.
     * @return true if triggered events are increasing events
     */
    virtual bool is_triggered_on_increasing();

    /** Get next function transformer in the specified direction.
     * @param previous transformer active on the previous point with respect
     * to integration direction (may be NULL if no previous point is known)
     * @param g current value of the g function
     * @param forward true if integration goes forward
     * @return next transformer transformer
     */
    virtual Transformer select_transformer(Transformer previous, double g, bool forward);

}