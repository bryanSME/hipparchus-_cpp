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
//package org.hipparchus.migration.exception;

//import org.hipparchus.migration.exception.util.Localized_Formats;
//import org.hipparchus.util.Math_Arrays;

/**
 * Exception to be thrown when the a sequence of values is not monotonically
 * increasing or decreasing.
 *
 * @deprecated as of 1.0, this exception is replaced by {@link org.hipparchus.exception.}
 */
//@Deprecated
class Non_Monotonic_Sequence_Exception //: Math_illegalNumberException 
{
    /**
     * Direction (positive for increasing, negative for decreasing).
     */
    private const Math_Arrays.Order_Direction direction;
    /**
     * Whether the sequence must be strictly increasing or decreasing.
     */
    private const bool strict;
    /**
     * Index of the wrong value.
     */
    private const int index;
    /**
     * Previous value.
     */
    private const Number previous;

    /**
     * Construct the exception.
     * This constructor uses default values assuming that the sequence should
     * have been strictly increasing.
     *
     * @param wrong Value that did not match the requirements.
     * @param previous Previous value in the sequence.
     * @param index Index of the value that did not match the requirements.
     */
    public Non_Monotonic_Sequence_Exception(Number wrong, Number previous, int index)
    {
        this(wrong, previous, index, Math_Arrays.Order_Direction::INCREASING, true);
    }

    /**
     * Construct the exception.
     *
     * @param wrong Value that did not match the requirements.
     * @param previous Previous value in the sequence.
     * @param index Index of the value that did not match the requirements.
     * @param direction Strictly positive for a sequence required to be
     * increasing, negative (or zero) for a decreasing sequence.
     * @param strict Whether the sequence must be strictly increasing or
     * decreasing.
     */
    public Non_Monotonic_Sequence_Exception(Number wrong, Number previous, int index, Math_Arrays.Order_Direction direction, bool strict)
    {
        super(direction == Math_Arrays.Order_Direction::INCREASING ?
            (strict ?
                Localized_Formats.NOT_STRICTLY_INCREASING_SEQUENCE :
                Localized_Formats.NOT_INCREASING_SEQUENCE) :
            (strict ?
                Localized_Formats.NOT_STRICTLY_DECREASING_SEQUENCE :
                Localized_Formats.NOT_DECREASING_SEQUENCE), wrong, previous, Integer.value_of(index), Integer.value_of(index - 1));

        this.direction = direction;
        this.strict = strict;
        this.index = index;
        this.previous = previous;
    }

    /**
     * @return the order direction.
     **/
    public Math_Arrays.Order_Direction get_direction()
    {
        return direction;
    }
    /**
     * @return {@code true} is the sequence should be strictly monotonic.
     **/
    public bool get_strict() { // NOPMD - this method name is for a legacy API we cannot change
        return strict;
    }
    /**
     * Get the index of the wrong value.
     *
     * @return the current index.
     */
    public int get_index()
    {
        return index;
    }
    /**
     * @return the previous value.
     */
    public Number get_previous()
    {
        return previous;
    }
};