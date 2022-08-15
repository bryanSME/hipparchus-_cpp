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
//package org.hipparchus.exception;

//import java.util.Missing_Resource_Exception;
//import java.util.Resource_Bundle;
#include <string>
#include <locale>
#include "Localizable.h"
#include <unordered_map>

namespace hipparchus
{
    namespace exception
    {
        namespace _internal
        {
            std::string get_localized_core_format_str(const Localized_Core_Formats_Type& lcf)
            {
                switch (lcf)
                {
                case(ARRAY_SIZE_EXCEEDS_MAX_VARIABLES):
                    return "array size cannot be greater than {0}";
                case(ARRAY_SIZES_SHOULD_HAVE_DIFFERENCE_1):
                    return "array sizes should have difference 1 ({0} != {1} + 1)";
                case(ARRAY_SUMS_TO_ZERO):
                    return "array sums to zero";
                case(AT_LEAST_ONE_COLUMN):
                    return "matrix must have at least one column";
                case(AT_LEAST_ONE_ROW):
                    return "matrix must have at least one row";
                case(BANDWIDTH):
                    return "bandwidth ({0})";
                case(BESSEL_FUNCTION_BAD_ARGUMENT):
                    return "Bessel function of order {0} cannot be computed for x = {1}";
                case(BESSEL_FUNCTION_FAILED_CONVERGENCE):
                    return "Bessel function of order {0} failed to converge for x = {1}";
                case(BINOMIAL_INVALID_PARAMETERS_ORDER):
                    return "must have n >= k for binomial coefficient (n, k), got k = {0}, n = {1}";
                case(BINOMIAL_NEGATIVE_PARAMETER):
                    return "must have n >= 0 for binomial coefficient (n, k), got n = {0}";
                case(CANNOT_COMPUTE_0TH_ROOT_OF_UNITY):
                    return "cannot compute 0-th root of unity, indefinite result";
                case(CANNOT_COMPUTE_BETA_DENSITY_AT_0_FOR_SOME_ALPHA):
                    return "cannot compute beta density at 0 when alpha = {0,number}";
                case(CANNOT_COMPUTE_BETA_DENSITY_AT_1_FOR_SOME_BETA):
                    return "cannot compute beta density at 1 when beta = %.3g";
                case(CANNOT_COMPUTE_NTH_ROOT_FOR_NEGATIVE_N):
                    return "cannot compute nth root for NULL or negative n: {0}";
                case(CANNOT_DISCARD_NEGATIVE_NUMBER_OF_ELEMENTS):
                    return "cannot discard a negative number of elements ({0})";
                case(CANNOT_FORMAT_INSTANCE_AS_COMPLEX):
                    return "cannot format a {0} instance as a complex number";
                case(CANNOT_FORMAT_OBJECT_TO_FRACTION):
                    return "cannot format given object as a fraction number";
                case(CANNOT_SUBSTITUTE_ELEMENT_FROM_EMPTY_ARRAY):
                    return "cannot substitute an element from an empty array";
                case(COLUMN_INDEX):
                    return "column index ({0})" /* keep */;
                case(COMPLEX_CANNOT_BE_CONSIDERED_A_REAL_NUMBER):
                        return "complex number ({0},{1}) cannot be considered to be a real number";
                case(CONSTRAINT):
                        return "constraint" /* keep */;
                case(CONTINUED_FRACTION_INFINITY_DIVERGENCE):
                        return "Continued fraction convergents diverged to +/- infinity for value {0}";
                case(CONTINUED_FRACTION_NAN_DIVERGENCE):
                        return "Continued fraction diverged to NaN for value {0}";
                case(CONTRACTION_CRITERIA_SMALLER_THAN_EXPANSION_FACTOR):
                        return "contraction criteria ({0}) smaller than the expansion factor ({1}). This would lead to a never ending loop of expansion and contraction as a newly expanded internal storage array would immediately satisfy the criteria for contraction.";
                case(CONTRACTION_CRITERIA_SMALLER_THAN_ONE):
                        return "contraction criteria smaller than one ({0}). This would lead to a never ending loop of expansion and contraction as an internal storage array length equal to the number of elements would satisfy the contraction criteria.";
                case(CONVERGENCE_FAILED):
                        return "convergence failed" /* keep */;
                case(CUMULATIVE_PROBABILITY_RETURNED_NAN):
                        return "Cumulative probability function returned NaN for argument {0} p = {1}";
                case(DERIVATION_ORDER_NOT_ALLOWED):
                        return "derivation order {0} is not allowed here";
                case(DIFFERENT_ROWS_LENGTHS):
                        return "some rows have length {0} while others have length {1}";
                case(DIGEST_NOT_INITIALIZED):
                        return "digest not initialized";
                case(DIMENSIONS_MISMATCH_2x2):
                        return "got {0}x{1} but expected {2}x{3}" /* keep */;
                case(DIMENSIONS_MISMATCH):
                        return "inconsistent dimensions: {0} != {1}" /* keep */;
                case(DISCRETE_CUMULATIVE_PROBABILITY_RETURNED_NAN):
                        return "Discrete cumulative probability function returned NaN for argument {0}";
                case(DISTRIBUTION_NOT_LOADED):
                        return "distribution not loaded";
                case(DUPLICATED_ABSCISSA_DIVISION_BY_ZERO):
                        return "duplicated abscissa {0} causes division by zero";
                case(EMPTY_INTERPOLATION_SAMPLE):
                        return "sample for interpolation is empty";
                case(EMPTY_POLYNOMIALS_COEFFICIENTS_ARRAY):
                        return "empty polynomials coefficients array" /* keep */;
                case(EMPTY_SELECTED_COLUMN_INDEX_ARRAY):
                        return "empty selected column index array";
                case(EMPTY_SELECTED_ROW_INDEX_ARRAY):
                        return "empty selected row index array";
                case(ENDPOINTS_NOT_AN_INTERVAL):
                        return "endpoints do not specify an interval: [{0}, {1}]";
                case(EVALUATION):
                        return "evaluation" /* keep */;
                case(EXPANSION_FACTOR_SMALLER_THAN_ONE):
                        return "expansion factor smaller than one ({0})";
                case(FACTORIAL_NEGATIVE_PARAMETER):
                        return "must have n >= 0 for n!, got n = {0}";
                case(FAILED_BRACKETING):
                        return "number of iterations={4}, maximum iterations={5}, initial={6}, lower bound={7}, upper bound={8}, const a value={0}, const b value={1}, f(a)={2}, f(b)={3}";
                case(FAILED_DECOMPOSITION):
                        return "failed decomposition of a {0}x{1} matrix";
                case(FAILED_FRACTION_CONVERSION):
                        return "Unable to convert {0} to fraction after {1} iterations";
                case(FIRST_COLUMNS_NOT_INITIALIZED_YET):
                        return "first {0} columns are not initialized yet";
                case(FIRST_ROWS_NOT_INITIALIZED_YET):
                        return "first {0} rows are not initialized yet";
                case(FRACTION_CONVERSION_OVERFLOW):
                        return "Overflow trying to convert {0} to fraction ({1}/{2})";
                case(GCD_OVERFLOW_32_BITS):
                        return "overflow: gcd({0}, {1}) is 2^31";
                case(GCD_OVERFLOW_64_BITS):
                        return "overflow: gcd({0}, {1}) is 2^63";
                case(ILL_CONDITIONED_OPERATOR):
                        return "condition number {1} is too high ";
                case(INDEX_LARGER_THAN_MAX):
                        return "the index specified: {0} is larger than the current maximal index {1}";
                case(INDEX_NOT_POSITIVE):
                        return "index ({0}) is not positive";
                case(NOT_FINITE_NUMBER):
                        return "{0} is not a finite number" /* keep */;
                case(INFINITE_BOUND):
                        return "interval bounds must be finite";
                case(ARRAY_ELEMENT):
                        return "value {0} at index {1}" /* keep */;
                case(INFINITE_ARRAY_ELEMENT):
                        return "Array contains an infinite element, {0} at index {1}";
                case(INFINITE_VALUE_CONVERSION):
                        return "cannot convert infinite value";
                case(INITIAL_CAPACITY_NOT_POSITIVE):
                        return "initial capacity ({0}) is not positive";
                case(INITIAL_COLUMN_AFTER_FINAL_COLUMN):
                        return "initial column {1} after const column {0}";
                case(INITIAL_ROW_AFTER_FINAL_ROW):
                        return "initial row {1} after const row {0}";
                case(INSUFFICIENT_DATA):
                        return "insufficient data";
                case(INSUFFICIENT_DIMENSION):
                        return "insufficient dimension {0}, must be at least {1}";
                case(DIMENSION):
                        return "dimension ({0})" /* keep */;
                case(INSUFFICIENT_OBSERVED_POINTS_IN_SAMPLE):
                        return "sample contains {0} observed points, at least {1} are required";
                case(INSUFFICIENT_ROWS_AND_COLUMNS):
                        return "insufficient data: only {0} rows and {1} columns.";
                case(INTERNAL_ERROR):
                        return "internal error, please fill a bug report at {0}";
                case(INVALID_MAX_ITERATIONS):
                        return "bad value for maximum iterations number: {0}";
                case(INVALID_ROUNDING_METHOD):
                        return "invalid rounding method {0}, valid methods: {1} ({2}), {3} ({4}), {5} ({6}), {7} ({8}), {9} ({10}), {11} ({12}), {13} ({14}), {15} ({16})";
                case(ITERATIONS):
                        return "iterations" /* keep */;
                case(LCM_OVERFLOW_32_BITS):
                        return "overflow: lcm({0}, {1}) is 2^31";
                case(LCM_OVERFLOW_64_BITS):
                        return "overflow: lcm({0}, {1}) is 2^63";
                case(LOWER_BOUND_NOT_BELOW_UPPER_BOUND):
                        return "lower bound ({0}) must be strictly less than upper bound ({1})" /* keep */;
                case(LOWER_ENDPOINT_ABOVE_UPPER_ENDPOINT):
                        return "lower endpoint ({0}) must be less than or equal to upper endpoint ({1})";
                case(EVALUATIONS):
                        return "evaluations" /* keep */;
                case(MAX_COUNT_EXCEEDED):
                        return "maximal count ({0}) exceeded" /* keep */;
                case(NAN_ELEMENT_AT_INDEX):
                        return "element {0} is NaN";
                case(NAN_VALUE_CONVERSION):
                        return "cannot convert NaN value";
                case(NEGATIVE_COMPLEX_MODULE):
                        return "negative complex module {0}";
                case(NEGATIVE_ELEMENT_AT_INDEX):
                        return "element {0} is negative: {1}";
                case(NUMBER_OF_SUCCESSES):
                        return "number of successes ({0})" /* keep */;
                case(NUMBER_OF_INTERPOLATION_POINTS):
                        return "number of interpolation points ({0})" /* keep */;
                case(NUMBER_OF_TRIALS):
                        return "number of trials ({0})";
                case(ROBUSTNESS_ITERATIONS):
                        return "number of robustness iterations ({0})";
                case(START_POSITION):
                        return "start position ({0})" /* keep */;
                case(NON_CONVERGENT_CONTINUED_FRACTION):
                        return "Continued fraction convergents failed to converge (in less than {0} iterations) for value {1}";
                case(NON_SQUARE_MATRIX):
                        return "non square ({0}x{1}) matrix";
                case(NORM):
                        return "Norm ({0})" /* keep */;
                case(NORMALIZE_INFINITE):
                        return "Cannot normalize to an infinite value";
                case(NORMALIZE_NAN):
                        return "Cannot normalize to NaN";
                case(NOT_DECREASING_SEQUENCE):
                        return "points {3} and {2} are not decreasing ({1} < {0})" /* keep */;
                case(NOT_ENOUGH_POINTS_IN_SPLINE_PARTITION):
                        return "spline partition must have at least {0} points, got {1}";
                case(NOT_INCREASING_SEQUENCE):
                        return "points {3} and {2} are not increasing ({1} > {0})" /* keep */;
                case(NOT_POSITIVE_DEFINITE_MATRIX):
                        return "not positive definite matrix" /* keep */;
                case(NON_POSITIVE_DEFINITE_OPERATOR):
                        return "non positive definite linear operator" /* keep */;
                case(NON_SELF_ADJOINT_OPERATOR):
                        return "non self-adjoint linear operator" /* keep */;
                case(NON_SQUARE_OPERATOR):
                        return "non square ({0}x{1}) linear operator" /* keep */;
                case(DEGREES_OF_FREEDOM):
                        return "degrees of freedom ({0})" /* keep */;
                case(NOT_POSITIVE_EXPONENT):
                        return "invalid exponent {0} (must be positive)";
                case(NUMBER_OF_ELEMENTS_SHOULD_BE_POSITIVE):
                        return "number of elements should be positive ({0})";
                case(BASE):
                        return "base ({0})" /* keep */;
                case(EXPONENT):
                        return "exponent ({0})" /* keep */;
                case(LENGTH):
                        return "length ({0})" /* keep */;
                case(MEAN):
                        return "mean ({0})" /* keep */;
                case(NOT_POSITIVE_NUMBER_OF_SAMPLES):
                        return "number of sample is not positive: {0}";
                case(NUMBER_OF_SAMPLES):
                        return "number of samples ({0})" /* keep */;
                case(PERMUTATION_SIZE):
                        return "permutation size ({0}" /* keep */;
                case(POPULATION_SIZE):
                        return "population size ({0})" /* keep */;
                case(NOT_POSITIVE_SCALE):
                        return "scale must be positive ({0})";
                case(SCALE):
                        return "scale ({0})" /* keep */;
                case(SHAPE):
                        return "shape ({0})" /* keep */;
                case(STANDARD_DEVIATION):
                        return "standard deviation ({0})" /* keep */;
                case(NOT_POSITIVE_WINDOW_SIZE):
                        return "window size must be positive ({0})";
                case(NOT_STRICTLY_DECREASING_SEQUENCE):
                        return "points {3} and {2} are not strictly decreasing ({1} <= {0})" /* keep */;
                case(NOT_STRICTLY_INCREASING_SEQUENCE):
                        return "points {3} and {2} are not strictly increasing ({1} >= {0})" /* keep */;
                case(NON_SYMMETRIC_MATRIX):
                        return "non symmetric matrix: the difference between entries at ({0},{1}) and ({1},{0}) is larger than {2}" /* keep */;
                case(NO_CONVERGENCE_WITH_ANY_START_POINT):
                        return "none of the {0} start points lead to convergence" /* keep */;
                case(NO_DATA):
                        return "no data" /* keep */;
                case(NO_OPTIMUM_COMPUTED_YET):
                        return "no optimum computed yet" /* keep */;
                case(NAN_NOT_ALLOWED):
                        return "NaN is not allowed";
                case(NULL_NOT_ALLOWED):
                        return "null is not allowed" /* keep */;
                case(ARRAY_ZERO_LENGTH_OR_NULL_NOT_ALLOWED):
                        return "a NULL or zero length array not allowed";
                case(DENOMINATOR):
                        return "denominator" /* keep */;
                case(DENOMINATOR_FORMAT):
                        return "denominator format" /* keep */;
                case(FRACTION):
                        return "fraction" /* keep */;
                case(FUNCTION):
                        return "function" /* keep */;
                case(IMAGINARY_FORMAT):
                        return "imaginary format" /* keep */;
                case(INPUT_ARRAY):
                        return "input array" /* keep */;
                case(NUMERATOR):
                        return "numerator" /* keep */;
                case(NUMERATOR_FORMAT):
                        return "numerator format" /* keep */;
                case(REAL_FORMAT):
                        return "real format" /* keep */;
                case(WHOLE_FORMAT):
                        return "whole format" /* keep */;
                case(NUMBER_TOO_LARGE):
                        return "{0} is larger than the maximum ({1})" /* keep */;
                case(NUMBER_TOO_SMALL):
                        return "{0} is smaller than the minimum ({1})" /* keep */;
                case(NUMBER_TOO_LARGE_BOUND_EXCLUDED):
                        return "{0} is larger than, or equal to, the maximum ({1})" /* keep */;
                case(NUMBER_TOO_SMALL_BOUND_EXCLUDED):
                        return "{0} is smaller than, or equal to, the minimum ({1})" /* keep */;
                case(NUMBER_OF_SUCCESS_LARGER_THAN_POPULATION_SIZE):
                        return "number of successes ({0}) must be less than or equal to population size ({1})";
                case(NUMERATOR_OVERFLOW_AFTER_MULTIPLY):
                        return "overflow, numerator too large after multiply: {0}";
                case(OBSERVED_COUNTS_BOTTH_ZERO_FOR_ENTRY):
                        return "observed counts are both zero for entry {0}";
                case(OUT_OF_RANGE_ROOT_OF_UNITY_INDEX):
                        return "out of range root of unity index {0} (must be in [{1};{2}])";
                case(OUT_OF_RANGE):
                        return "out of range" /* keep */;
                case(OUT_OF_RANGE_SIMPLE):
                        return "{0} out of [{1}, {2}] range" /* keep */;
                case(OUT_OF_RANGE_LEFT):
                        return "{0} out of ({1}, {2}] range";
                //case(OVERFLOW):
                //        return "overflow" /* keep */;
                case(OVERFLOW_IN_FRACTION):
                        return "overflow in fraction {0}/{1}, cannot negate";
                case(OVERFLOW_IN_ADDITION):
                        return "overflow in addition: {0} + {1}";
                case(OVERFLOW_IN_SUBTRACTION):
                        return "overflow in subtraction: {0} - {1}";
                case(OVERFLOW_IN_MULTIPLICATION):
                        return "overflow in multiplication: {0} * {1}";
                case(PERMUTATION_EXCEEDS_N):
                        return "permutation size ({0}) exceeds permuation domain ({1})" /* keep */;
                case(POLYNOMIAL):
                        return "polynomial" /* keep */;
                case(ROOTS_OF_UNITY_NOT_COMPUTED_YET):
                        return "roots of unity have not been computed yet";
                case(ROW_INDEX):
                        return "row index ({0})" /* keep */;
                case(NOT_BRACKETING_INTERVAL):
                        return "interval does not bracket a root: f({0,number,##0.################E0}) = {2,number,##0.################E0}, f({1,number,##0.################E0}) = {3,number,##0.################E0}";
                case(START_POINT_NOT_IN_INTERVAL):
                        return "The start point {0} is not in the interval [{1}, {2}]";
                case(SAMPLE_SIZE_EXCEEDS_COLLECTION_SIZE):
                        return "sample size ({0}) exceeds collection size ({1})" /* keep */;
                case(SAMPLE_SIZE_LARGER_THAN_POPULATION_SIZE):
                        return "sample size ({0}) must be less than or equal to population size ({1})";
                case(SIMPLE_MESSAGE):
                        return "{0}";
                case(SINGULAR_MATRIX):
                        return "matrix is singular" /* keep */;
                case(SINGULAR_OPERATOR):
                        return "operator is singular";
                case(SUBARRAY_ENDS_AFTER_ARRAY_END):
                        return "subarray ends after array end";
                case(TOO_LARGE_CUTOFF_SINGULAR_VALUE):
                        return "cutoff singular value is {0}, should be at most {1}";
                case(TOO_MANY_ELEMENTS_TO_DISCARD_FROM_ARRAY):
                        return "cannot discard {0} elements from a {1} elements array";
                case(UNKNOWN_MODE):
                        return "unknown mode {0}, known modes: {1} ({2}), {3} ({4}), {5} ({6}), {7} ({8}), {9} ({10}) and {11} ({12})";
                case(CANNOT_PARSE_AS_TYPE):
                        return "string \"{0}\" unparseable (from position {1}) as an object of type {2}" /* keep */;
                case(CANNOT_PARSE):
                        return "string \"{0}\" unparseable (from position {1})" /* keep */;
                case(UNSUPPORTED_OPERATION):
                        return "unsupported operation" /* keep */;
                case(ARITHMETIC_EXCEPTION):
                        return "arithmetic exception" /* keep */;
                case(ILLEGAL_STATE):
                        return "illegal state" /* keep */;
                case(USER_EXCEPTION):
                        return "exception generated in user code" /* keep */;
                case(URL_CONTAINS_NO_DATA):
                        return "URL {0} contains no data";
                case(VECTOR_MUST_HAVE_AT_LEAST_ONE_ELEMENT):
                        return "vector must have at least one element";
                case(WEIGHT_AT_LEAST_ONE_NON_ZERO):
                        return "weight array must contain at least one non-zero value";
                case(WRONG_NUMBER_OF_POINTS):
                        return "{0} points are required, got only {1}";
                case(NUMBER_OF_POINTS):
                        return "number of points ({0})" /* keep */;
                case(ZERO_DENOMINATOR):
                        return "denominator must be different from 0" /* keep */;
                case(ZERO_DENOMINATOR_IN_FRACTION):
                        return "zero denominator in fraction {0}/{1}";
                case(ZERO_FRACTION_TO_DIVIDE_BY):
                        return "the fraction to divide by must not be zero: {0}/{1}";
                case(ZERO_NORM):
                        return "zero norm";
                case(ZERO_NOT_ALLOWED):
                        return "zero not allowed here";
                }
            }
        };
        enum Localized_Core_Formats_Type
        {
            ARRAY_SIZE_EXCEEDS_MAX_VARIABLES,
            ARRAY_SIZES_SHOULD_HAVE_DIFFERENCE_1,
            ARRAY_SUMS_TO_ZERO,
            AT_LEAST_ONE_COLUMN,
            AT_LEAST_ONE_ROW,
            BANDWIDTH,
            BESSEL_FUNCTION_BAD_ARGUMENT,
            BESSEL_FUNCTION_FAILED_CONVERGENCE,
            BINOMIAL_INVALID_PARAMETERS_ORDER,
            BINOMIAL_NEGATIVE_PARAMETER,
            CANNOT_COMPUTE_0TH_ROOT_OF_UNITY,
            CANNOT_COMPUTE_BETA_DENSITY_AT_0_FOR_SOME_ALPHA,
            CANNOT_COMPUTE_BETA_DENSITY_AT_1_FOR_SOME_BETA,
            CANNOT_COMPUTE_NTH_ROOT_FOR_NEGATIVE_N,
            CANNOT_DISCARD_NEGATIVE_NUMBER_OF_ELEMENTS,
            CANNOT_FORMAT_INSTANCE_AS_COMPLEX,
            CANNOT_FORMAT_OBJECT_TO_FRACTION,
            CANNOT_SUBSTITUTE_ELEMENT_FROM_EMPTY_ARRAY,
            COLUMN_INDEX,
            COMPLEX_CANNOT_BE_CONSIDERED_A_REAL_NUMBER,
            CONSTRAINT,
            CONTINUED_FRACTION_INFINITY_DIVERGENCE,
            CONTINUED_FRACTION_NAN_DIVERGENCE,
            CONTRACTION_CRITERIA_SMALLER_THAN_EXPANSION_FACTOR,
            CONTRACTION_CRITERIA_SMALLER_THAN_ONE,
            CONVERGENCE_FAILED,
            CUMULATIVE_PROBABILITY_RETURNED_NAN,
            DERIVATION_ORDER_NOT_ALLOWED,
            DIFFERENT_ROWS_LENGTHS,
            DIGEST_NOT_INITIALIZED,
            DIMENSIONS_MISMATCH_2x2,
            DIMENSIONS_MISMATCH,
            DISCRETE_CUMULATIVE_PROBABILITY_RETURNED_NAN,
            DISTRIBUTION_NOT_LOADED,
            DUPLICATED_ABSCISSA_DIVISION_BY_ZERO,
            EMPTY_INTERPOLATION_SAMPLE,
            EMPTY_POLYNOMIALS_COEFFICIENTS_ARRAY,
            EMPTY_SELECTED_COLUMN_INDEX_ARRAY,
            EMPTY_SELECTED_ROW_INDEX_ARRAY,
            ENDPOINTS_NOT_AN_INTERVAL,
            EVALUATION,
            EXPANSION_FACTOR_SMALLER_THAN_ONE,
            FACTORIAL_NEGATIVE_PARAMETER,
            FAILED_BRACKETING,
            FAILED_DECOMPOSITION,
            FAILED_FRACTION_CONVERSION,
            FIRST_COLUMNS_NOT_INITIALIZED_YET,
            FIRST_ROWS_NOT_INITIALIZED_YET,
            FRACTION_CONVERSION_OVERFLOW,
            GCD_OVERFLOW_32_BITS,
            GCD_OVERFLOW_64_BITS,
            ILL_CONDITIONED_OPERATOR,
            INDEX_LARGER_THAN_MAX,
            INDEX_NOT_POSITIVE,
            NOT_FINITE_NUMBER,
            INFINITE_BOUND,
            ARRAY_ELEMENT,
            INFINITE_ARRAY_ELEMENT,
            INFINITE_VALUE_CONVERSION,
            INITIAL_CAPACITY_NOT_POSITIVE,
            INITIAL_COLUMN_AFTER_FINAL_COLUMN,
            INITIAL_ROW_AFTER_FINAL_ROW,
            INSUFFICIENT_DATA,
            INSUFFICIENT_DIMENSION,
            DIMENSION,
            INSUFFICIENT_OBSERVED_POINTS_IN_SAMPLE,
            INSUFFICIENT_ROWS_AND_COLUMNS,
            INTERNAL_ERROR,
            INVALID_MAX_ITERATIONS,
            INVALID_ROUNDING_METHOD,
            ITERATIONS,
            LCM_OVERFLOW_32_BITS,
            LCM_OVERFLOW_64_BITS,
            LOWER_BOUND_NOT_BELOW_UPPER_BOUND,
            LOWER_ENDPOINT_ABOVE_UPPER_ENDPOINT,
            EVALUATIONS,
            MAX_COUNT_EXCEEDED,
            NAN_ELEMENT_AT_INDEX,
            NAN_VALUE_CONVERSION,
            NEGATIVE_COMPLEX_MODULE,
            NEGATIVE_ELEMENT_AT_INDEX,
            NUMBER_OF_SUCCESSES,
            NUMBER_OF_INTERPOLATION_POINTS,
            NUMBER_OF_TRIALS,
            ROBUSTNESS_ITERATIONS,
            START_POSITION,
            NON_CONVERGENT_CONTINUED_FRACTION,
            NON_SQUARE_MATRIX,
            NORM,
            NORMALIZE_INFINITE,
            NORMALIZE_NAN,
            NOT_DECREASING_SEQUENCE,
            NOT_ENOUGH_POINTS_IN_SPLINE_PARTITION,
            NOT_INCREASING_SEQUENCE,
            NOT_POSITIVE_DEFINITE_MATRIX,
            NON_POSITIVE_DEFINITE_OPERATOR,
            NON_SELF_ADJOINT_OPERATOR,
            NON_SQUARE_OPERATOR,
            DEGREES_OF_FREEDOM,
            NOT_POSITIVE_EXPONENT,
            NUMBER_OF_ELEMENTS_SHOULD_BE_POSITIVE,
            BASE,
            EXPONENT,
            LENGTH,
            MEAN,
            NOT_POSITIVE_NUMBER_OF_SAMPLES,
            NUMBER_OF_SAMPLES,
            PERMUTATION_SIZE,
            POPULATION_SIZE,
            NOT_POSITIVE_SCALE,
            SCALE,
            SHAPE,
            STANDARD_DEVIATION,
            NOT_POSITIVE_WINDOW_SIZE,
            NOT_STRICTLY_DECREASING_SEQUENCE,
            NOT_STRICTLY_INCREASING_SEQUENCE,
            NON_SYMMETRIC_MATRIX,
            NO_CONVERGENCE_WITH_ANY_START_POINT,
            NO_DATA,
            NO_OPTIMUM_COMPUTED_YET,
            NAN_NOT_ALLOWED,
            NULL_NOT_ALLOWED,
            ARRAY_ZERO_LENGTH_OR_NULL_NOT_ALLOWED,
            DENOMINATOR,
            DENOMINATOR_FORMAT,
            FRACTION,
            FUNCTION,
            IMAGINARY_FORMAT,
            INPUT_ARRAY,
            NUMERATOR,
            NUMERATOR_FORMAT,
            REAL_FORMAT,
            WHOLE_FORMAT,
            NUMBER_TOO_LARGE,
            NUMBER_TOO_SMALL,
            NUMBER_TOO_LARGE_BOUND_EXCLUDED,
            NUMBER_TOO_SMALL_BOUND_EXCLUDED,
            NUMBER_OF_SUCCESS_LARGER_THAN_POPULATION_SIZE,
            NUMERATOR_OVERFLOW_AFTER_MULTIPLY,
            OBSERVED_COUNTS_BOTTH_ZERO_FOR_ENTRY,
            OUT_OF_RANGE_ROOT_OF_UNITY_INDEX,
            OUT_OF_RANGE,
            OUT_OF_RANGE_SIMPLE,
            OUT_OF_RANGE_LEFT,
            //OVERFLOW,
            OVERFLOW_IN_FRACTION,
            OVERFLOW_IN_ADDITION,
            OVERFLOW_IN_SUBTRACTION,
            OVERFLOW_IN_MULTIPLICATION,
            PERMUTATION_EXCEEDS_N,
            POLYNOMIAL,
            ROOTS_OF_UNITY_NOT_COMPUTED_YET,
            ROW_INDEX,
            NOT_BRACKETING_INTERVAL,
            START_POINT_NOT_IN_INTERVAL,
            SAMPLE_SIZE_EXCEEDS_COLLECTION_SIZE,
            SAMPLE_SIZE_LARGER_THAN_POPULATION_SIZE,
            SIMPLE_MESSAGE,
            SINGULAR_MATRIX,
            SINGULAR_OPERATOR,
            SUBARRAY_ENDS_AFTER_ARRAY_END,
            TOO_LARGE_CUTOFF_SINGULAR_VALUE,
            TOO_MANY_ELEMENTS_TO_DISCARD_FROM_ARRAY,
            UNKNOWN_MODE,
            CANNOT_PARSE_AS_TYPE,
            CANNOT_PARSE,
            UNSUPPORTED_OPERATION,
            ARITHMETIC_EXCEPTION,
            ILLEGAL_STATE,
            USER_EXCEPTION,
            URL_CONTAINS_NO_DATA,
            VECTOR_MUST_HAVE_AT_LEAST_ONE_ELEMENT,
            WEIGHT_AT_LEAST_ONE_NON_ZERO,
            WRONG_NUMBER_OF_POINTS,
            NUMBER_OF_POINTS,
            ZERO_DENOMINATOR,
            ZERO_DENOMINATOR_IN_FRACTION,
            ZERO_FRACTION_TO_DIVIDE_BY,
            ZERO_NORM,
            ZERO_NOT_ALLOWED
        };

        /**
         * Enumeration for localized messages formats used in exceptions messages.
         * <p>
         * The constants in this enumeration represent the available
         * formats as localized strings. These formats are intended to be
         * localized using simple properties files, using the constant
         * name as the key and the property value as the message format.
         * The source English format is provided in the constants themselves
         * to serve both as a reminder for developers to understand the parameters
         * needed by each format, as a basis for translators to create
         * localized properties files, and as a default format if some
         * translation is missing.
         * </p>
         */
        //template<typename T, typename std::enable_if<std::is_base_of<Localized_Core_Formats_Type, T>::value>::type* = nullptr>
        class Localized_Core_Formats : public Localizable
        {
        private:
            /** Source English format. */
            std::string my_source_format{};

        public:
            /** Simple constructor.
             * @param source_format source English format to use when no
             * localized version is available
             */
            Localized_Core_Formats(const std::string& source_format)
            {
                my_source_format = source_format;
            }

            /** {@inherit_doc} */
            //override
            std::string get_source_string() const override
            {
                return my_source_format;
            }

            /** {@inherit_doc} */
            //override
            std::string get_localized_string(const std::locale& locale) override
            {
                //try
                //{
                    throw std::exception("LocalizedCoreFormats - not fully implemented");
                    //const std::string path = _internal::get_localized_core_format_str(T);
                    //const std::string path = hipparchus::exception::Localized_Core_Formats_Type::class.get_name().replace_all("\\.", "/");
                    //Resource_Bundle bundle = Resource_Bundle.get_bundle("assets/" + path, locale, UTF8_Control());
                    //if (bundle.get_locale().get_language().equals(locale.get_language()))
                    //{
                    //    const std::string translated = bundle.get_string(name());
                    //    if ((translated != NULL) &&
                    //        (translated.size()() > 0) &&
                    //        (!translated.to_lower_case(locale).contains("missing translation")))
                    //    {
                    //        // the value of the resource is the translated format
                    //        return translated;
                    //    }
                    //}

                //}

                //catch ([[maybe_unused]]Missing_Resource_Exception& mre) // NOPMD
                //{ 
                //    // do nothing here
                //}

                // either the locale is not supported or the resource is unknown
                // don't translate and fall back to using the source format
                return my_source_format;
            }

        };
    };
};