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
//package org.hipparchus.optim.linear;

//import java.util.Array_list;
//import java.util.List;

//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.optim.Localized_Optim_Formats;
//import org.hipparchus.optim.Optimization_data;
//import org.hipparchus.optim.Point_valuePair;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Precision;

/**
 * Solves a linear problem using the "Two-Phase Simplex" method.
 * <p>
 * The {@link Simplex_Solver} supports the following {@link Optimization_data} data provided
 * as arguments to {@link #optimize(Optimization_data...)}:
 * <ul>
 *   <li>objective function: {@link Linear_Objective_Function} - mandatory</li>
 *   <li>linear constraints {@link Linear_ConstraintSet} - mandatory</li>
 *   <li>type of optimization: {@link org.hipparchus.optim.nonlinear.scalar.Goal_Type Goal_Type}
 *    - optional, default: {@link org.hipparchus.optim.nonlinear.scalar.Goal_Type#MINIMIZE MINIMIZE}</li>
 *   <li>whether to allow negative values as solution: {@link Non_Negative_Constraint} - optional, default: true</li>
 *   <li>pivot selection rule: {@link PivotSelection_rule} - optional, default {@link PivotSelection_rule#DANTZIG}</li>
 *   <li>callback for the best solution: {@link Solution_callback} - optional</li>
 *   <li>maximum number of iterations: {@link org.hipparchus.optim.Max_Iter} - optional, default: {@link Integer#MAX_VALUE}</li>
 * </ul>
 * <p>
 * <b>Note:</b> Depending on the problem definition, the default convergence criteria
 * may be too strict, resulting in {@link Math_Illegal_State_Exception} or
 * {@link Math_Illegal_State_Exception}. In such a case it is advised to adjust these
 * criteria with more appropriate values, e.g. relaxing the epsilon value.
 * <p>
 * Default convergence criteria:
 * <ul>
 *   <li>Algorithm convergence: 1e-6</li>
 *   <li>Floating-point comparisons: 10 ulp</li>
 *   <li>Cut-Off value: 1e-10</li>
  * </ul>
 * <p>
 * The cut-off value has been introduced to handle the case of very small pivot elements
 * in the Simplex tableau, as these may lead to numerical instabilities and degeneracy.
 * Potential pivot elements smaller than this value will be treated as if they were zero
 * and are thus not considered by the pivot selection mechanism. The default value is safe
 * for many problems, but may need to be adjusted in case of very small coefficients
 * used in either the {@link Linear_Constraint} or {@link Linear_Objective_Function}.
 *
 */
class Simplex_Solver extends Linear_Optimizer 
{
    /** Default amount of error to accept in floating point comparisons (as ulps). */
    static const int DEFAULT_ULPS = 10;

    /** Default cut-off value. */
    static const double DEFAULT_CUT_OFF = 1e-10;

    /** Default amount of error to accept for algorithm convergence. */
    private static const double DEFAULT_EPSILON = 1.0e-6;

    /** Amount of error to accept for algorithm convergence. */
    private const double epsilon;

    /** Amount of error to accept in floating point comparisons (as ulps). */
    private const int max_ulps;

    /**
     * Cut-off value for entries in the tableau: values smaller than the cut-off
     * are treated as zero to improve numerical stability.
     */
    private const double cut_off;

    /** The pivot selection method to use. */
    private PivotSelection_rule pivot_selection;

    /**
     * The solution callback to access the best solution found so far in case
     * the optimizer fails to find an optimal solution within the iteration limits.
     */
    private Solution_callback solution_callback;

    /**
     * Builds a simplex solver with default settings.
     */
    public Simplex_Solver() 
    {
        this(DEFAULT_EPSILON, DEFAULT_ULPS, DEFAULT_CUT_OFF);
    }

    /**
     * Builds a simplex solver with a specified accepted amount of error.
     *
     * @param epsilon Amount of error to accept for algorithm convergence.
     */
    public Simplex_Solver(const double epsilon) 
    {
        this(epsilon, DEFAULT_ULPS, DEFAULT_CUT_OFF);
    }

    /**
     * Builds a simplex solver with a specified accepted amount of error.
     *
     * @param epsilon Amount of error to accept for algorithm convergence.
     * @param max_ulps Amount of error to accept in floating point comparisons.
     */
    public Simplex_Solver(const double epsilon, const int max_ulps) 
    {
        this(epsilon, max_ulps, DEFAULT_CUT_OFF);
    }

    /**
     * Builds a simplex solver with a specified accepted amount of error.
     *
     * @param epsilon Amount of error to accept for algorithm convergence.
     * @param max_ulps Amount of error to accept in floating point comparisons.
     * @param cut_off Values smaller than the cut_off are treated as zero.
     */
    public Simplex_Solver(const double epsilon, const int max_ulps, const double cut_off) 
    {
        this.epsilon = epsilon;
        this.max_ulps = max_ulps;
        this.cut_off = cut_off;
        this.pivot_selection = PivotSelection_rule.DANTZIG;
    }

    /**
     * {@inherit_doc}
     *
     * @param opt_data Optimization data. In addition to those documented in
     * {@link Linear_Optimizer#optimize(Optimization_data...)
     * Linear_Optimizer}, this method will register the following data:
     * <ul>
     *  <li>{@link Solution_callback}</li>
     *  <li>{@link PivotSelection_rule}</li>
     * </ul>
     *
     * @return {@inherit_doc}
     * @Math_Illegal_State_Exception if the maximal number of iterations is exceeded.
     * @org.hipparchus.exception. if the dimension
     * of the constraints does not match the dimension of the objective function
     */
    //override
    public Point_valuePair optimize(Optimization_data... opt_data)
        Math_Illegal_State_Exception 
        {
        // Set up base class and perform computation.
        return super.optimize(opt_data);
    }

    /**
     * {@inherit_doc}
     *
     * @param opt_data Optimization data.
     * In addition to those documented in
     * {@link Linear_Optimizer#parse_optimization_data(Optimization_data[])
     * Linear_Optimizer}, this method will register the following data:
     * <ul>
     *  <li>{@link Solution_callback}</li>
     *  <li>{@link PivotSelection_rule}</li>
     * </ul>
     */
    //override
    protected void parse_optimization_data(Optimization_data... opt_data) 
    {
        // Allow base class to register its own data.
        super.parse_optimization_data(opt_data);

        // reset the callback before parsing
        solution_callback = NULL;

        for (Optimization_data data : opt_data) 
        {
            if (data instanceof Solution_callback) 
            {
                solution_callback = (Solution_callback) data;
                continue;
            }
            if (data instanceof PivotSelection_rule) 
            {
                pivot_selection = (PivotSelection_rule) data;
                continue;
            }
        }
    }

    /**
     * Returns the column with the most negative coefficient in the objective function row.
     *
     * @param tableau Simple tableau for the problem.
     * @return the column with the most negative coefficient.
     */
    private Integer get_pivot_column(Simplex_Tableau tableau) 
    {
        double min_value = 0;
        Integer min_pos = NULL;
        for (int i = tableau.get_num_objective_functions(); i < tableau.get_width() - 1; i++) 
        {
            const double entry = tableau.get_entry(0, i);
            // check if the entry is strictly smaller than the current minimum
            // do not use a ulp/epsilon check
            if (entry < min_value) 
            {
                min_value = entry;
                min_pos = i;

                // Bland's rule: chose the entering column with the lowest index
                if (pivot_selection == PivotSelection_rule.BLAND && is_valid_pivot_column(tableau, i)) 
                {
                    break;
                }
            }
        }
        return min_pos;
    }

    /**
     * Checks whether the given column is valid pivot column, i.e. will result
     * in a valid pivot row.
     * <p>
     * When applying Bland's rule to select the pivot column, it may happen that
     * there is no corresponding pivot row. This method will check if the selected
     * pivot column will return a valid pivot row.
     *
     * @param tableau simplex tableau for the problem
     * @param col the column to test
     * @return {@code true} if the pivot column is valid, {@code false} otherwise
     */
    private bool is_valid_pivot_column(Simplex_Tableau tableau, int col) 
    {
        for (int i = tableau.get_num_objective_functions(); i < tableau.get_height(); i++) 
        {
            const double entry = tableau.get_entry(i, col);

            // do the same check as in get_pivot_row
            if (Precision.compare_to(entry, 0d, cut_off) > 0) 
            {
                return true;
            }
        }
        return false;
    }

    /**
     * Returns the row with the minimum ratio as given by the minimum ratio test (MRT).
     *
     * @param tableau Simplex tableau for the problem.
     * @param col Column to test the ratio of (see {@link #get_pivot_column(Simplex_Tableau)}).
     * @return the row with the minimum ratio.
     */
    private Integer get_pivot_row(Simplex_Tableau tableau, const int col) 
    {
        // create a list of all the rows that tie for the lowest score in the minimum ratio test
        List<Integer> min_ratio_positions = Array_list<>();
        double min_ratio = Double.MAX_VALUE;
        for (int i = tableau.get_num_objective_functions(); i < tableau.get_height(); i++) 
        {
            const double rhs = tableau.get_entry(i, tableau.get_width() - 1);
            const double entry = tableau.get_entry(i, col);

            // only consider pivot elements larger than the cut_off threshold
            // selecting others may lead to degeneracy or numerical instabilities
            if (Precision.compare_to(entry, 0d, cut_off) > 0) 
            {
                const double ratio = std::abs(rhs / entry);
                // check if the entry is strictly equal to the current min ratio
                // do not use a ulp/epsilon check
                const int cmp = Double.compare(ratio, min_ratio);
                if (cmp == 0) 
                {
                    min_ratio_positions.add(i);
                }
else if (cmp < 0) 
                {
                    min_ratio = ratio;
                    min_ratio_positions.clear();
                    min_ratio_positions.add(i);
                }
            }
        }

        if (min_ratio_positions.is_empty()) 
        {
            return NULL;
        }
else if (min_ratio_positions.size() > 1) 
        {
            // there's a degeneracy as indicated by a tie in the minimum ratio test

            // 1. check if there's an artificial variable that can be forced out of the basis
            if (tableau.get_num_artificial_variables() > 0) 
            {
                for (Integer row : min_ratio_positions) 
                {
                    for (int i{}; i < tableau.get_num_artificial_variables(); i++) 
                    {
                        int column = i + tableau.get_artificial_variable_offset();
                        const double entry = tableau.get_entry(row, column);
                        if (Precision.equals(entry, 1d, max_ulps) && row.equals(tableau.get_basic_row(column))) 
                        {
                            return row;
                        }
                    }
                }
            }

            // 2. apply Bland's rule to prevent cycling:
            //    take the row for which the corresponding basic variable has the smallest index
            //
            // see http://www.stanford.edu/class/msande310/blandrule.pdf
            // see http://en.wikipedia.org/wiki/Bland%27s_rule (not equivalent to the above paper)

            Integer min_row = NULL;
            int min_index = tableau.get_width();
            for (Integer row : min_ratio_positions) 
            {
                const int basic_var = tableau.get_basic_variable(row);
                if (basic_var < min_index) 
                {
                    min_index = basic_var;
                    min_row = row;
                }
            }
            return min_row;
        }
        return min_ratio_positions.get(0);
    }

    /**
     * Runs one iteration of the Simplex method on the given model.
     *
     * @param tableau Simple tableau for the problem.
     * @Math_Illegal_State_Exception if the allowed number of iterations has been exhausted.
     * @Math_Illegal_State_Exception if the model is found not to have a bounded solution.
     */
    protected void do_iteration(const Simplex_Tableau tableau)
        Math_Illegal_State_Exception 
        {

        increment_iteration_count();

        Integer pivot_col = get_pivot_column(tableau);
        Integer pivot_row = get_pivot_row(tableau, pivot_col);
        if (pivot_row == NULL) 
        {
            throw Math_Illegal_State_Exception(Localized_Optim_Formats.UNBOUNDED_SOLUTION);
        }

        tableau.perform_row_operations(pivot_col, pivot_row);
    }

    /**
     * Solves Phase 1 of the Simplex method.
     *
     * @param tableau Simple tableau for the problem.
     * @Math_Illegal_State_Exception if the allowed number of iterations has been exhausted, * or if the model is found not to have a bounded solution, or if there is no feasible solution
     */
    protected void solve_phase1(const Simplex_Tableau tableau)
        Math_Illegal_State_Exception 
        {

        // make sure we're in Phase 1
        if (tableau.get_num_artificial_variables() == 0) 
        {
            return;
        }

        while (!tableau.is_optimal()) 
        {
            do_iteration(tableau);
        }

        // if W is not zero then we have no feasible solution
        if (!Precision.equals(tableau.get_entry(0, tableau.get_rhs_offset()), 0d, epsilon)) 
        {
            throw Math_Illegal_State_Exception(Localized_Optim_Formats.NO_FEASIBLE_SOLUTION);
        }
    }

    /** {@inherit_doc} */
    //override
    public Point_valuePair do_optimize()
        Math_Illegal_State_Exception 
        {

        // reset the tableau to indicate a non-feasible solution in case
        // we do not pass phase 1 successfully
        if (solution_callback != NULL) 
        {
            solution_callback.set_tableau(null);
        }

        const Simplex_Tableau tableau =
            Simplex_Tableau(get_function(), get_constraints(), get_goal_type(), is_restricted_to_non_negative(), epsilon, max_ulps);

        solve_phase1(tableau);
        tableau.drop_phase1_objective();

        // after phase 1, we are sure to have a feasible solution
        if (solution_callback != NULL) 
        {
            solution_callback.set_tableau(tableau);
        }

        while (!tableau.is_optimal()) 
        {
            do_iteration(tableau);
        }

        // check that the solution respects the non_negative restriction in case
        // the epsilon/cut_off values are too large for the actual linear problem
        // (e.g. with very small constraint coefficients), the solver might actually
        // find a non-valid solution (with negative coefficients).
        const Point_valuePair solution = tableau.get_solution();
        if (is_restricted_to_non_negative()) 
        {
            const std::vector<double> coeff = solution.get_point();
            for (int i{}; i < coeff.size(); i++) 
            {
                if (Precision.compare_to(coeff[i], 0, epsilon) < 0) 
                {
                    throw Math_Illegal_State_Exception(Localized_Optim_Formats.NO_FEASIBLE_SOLUTION);
                }
            }
        }
        return solution;
    }
}


