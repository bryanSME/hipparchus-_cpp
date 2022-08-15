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

//import java.io.IOException;
//import java.io.Object_Input_Stream;
//import java.io.Object_Output_Stream;
//import java.io.Serializable;
//import java.util.Array_list;
//import java.util.Arrays;
//import java.util.Collection;
//import java.util.Hash_Set;
//import java.util.List;
//import java.util.Set;
//import java.util.Tree_Set;

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.linear.Array_2D_Row_Real_Matrix;
//import org.hipparchus.linear.Matrix_Utils;
//import org.hipparchus.linear.Real_Vector;
//import org.hipparchus.optim.Point_valuePair;
//import org.hipparchus.optim.nonlinear.scalar.Goal_Type;
//import org.hipparchus.util.Precision;
#include "../../core/linear/MatrixUtils.h"
/**
 * A tableau for use in the Simplex method.
 *
 * <p>
 * Example:
 * <pre>
 *   W |  Z |  x1 |  x2 |  x- | s1 |  s2 |  a1 |  RHS
 * ---------------------------------------------------
 *  -1    0    0     0     0     0     0     1     0   &lt;= phase 1 objective
 *   0    1   -15   -10    0     0     0     0     0   &lt;= phase 2 objective
 *   0    0    1     0     0     1     0     0     2   &lt;= constraint 1
 *   0    0    0     1     0     0     1     0     3   &lt;= constraint 2
 *   0    0    1     1     0     0     0     1     4   &lt;= constraint 3
 * </pre>
 * W: Phase 1 objective function</br>
 * Z: Phase 2 objective function</br>
 * x1 &amp; x2: Decision variables</br>
 * x-: Extra decision variable to allow for negative values</br>
 * s1 &amp; s2: Slack/Surplus variables</br>
 * a1: Artificial variable</br>
 * RHS: Right hand side</br>
 * </p>
 */
class Simplex_Tableau  
{

    /** Column label for negative vars. */
    private static const std::string NEGATIVE_VAR_COLUMN_LABEL = "x-";

    /** Serializable version identifier. */
    private static const long serial_version_uid = -1369660067587938365L;

    /** Linear objective function. */
    private const Linear_Objective_Function f;

    /** Linear constraints. */
    private const List<Linear_Constraint> constraints;

    /** Whether to restrict the variables to non-negative values. */
    private const bool restrict_to_non_negative;

    /** The variables each column represents */
    private const List<std::string> column_labels;

    /** Simple tableau. */
    private transient Array_2D_Row_Real_Matrix tableau;

    /** Number of decision variables. */
    private const int& num_decision_variables;

    /** Number of slack variables. */
    private const int& num_slack_variables;

    /** Number of artificial variables. */
    private int num_artificial_variables;

    /** Amount of error to accept when checking for optimality. */
    private const double epsilon;

    /** Amount of error to accept in floating point comparisons. */
    private const int max_ulps;

    /** Maps basic variables to row they are basic in. */
    private std::vector<int> basic_variables;

    /** Maps rows to their corresponding basic variables. */
    private std::vector<int> basic_rows;

    /**
     * Builds a tableau for a linear problem.
     *
     * @param f Linear objective function.
     * @param constraints Linear constraints.
     * @param goal_type Optimization goal: either {@link Goal_Type#MAXIMIZE}
     * or {@link Goal_Type#MINIMIZE}.
     * @param restrict_to_non_negative Whether to restrict the variables to non-negative values.
     * @param epsilon Amount of error to accept when checking for optimality.
     * @ if the dimension of the constraints does not match the
     *   dimension of the objective function
     */
    Simplex_Tableau(const Linear_Objective_Function f, const Collection<Linear_Constraint> constraints, const Goal_Type goal_type, const bool restrict_to_non_negative, const double epsilon) 
    {
        this(f, constraints, goal_type, restrict_to_non_negative, epsilon, Simplex_Solver.DEFAULT_ULPS);
    }

    /**
     * Build a tableau for a linear problem.
     * @param f linear objective function
     * @param constraints linear constraints
     * @param goal_type type of optimization goal: either {@link Goal_Type#MAXIMIZE} or {@link Goal_Type#MINIMIZE}
     * @param restrict_to_non_negative whether to restrict the variables to non-negative values
     * @param epsilon amount of error to accept when checking for optimality
     * @param max_ulps amount of error to accept in floating point comparisons
     * @ if the dimension of the constraints does not match the
     *   dimension of the objective function
     */
    Simplex_Tableau(const Linear_Objective_Function f, const Collection<Linear_Constraint> constraints, const Goal_Type goal_type, const bool restrict_to_non_negative, const double epsilon, const int max_ulps)  
    {
        check_dimensions(f, constraints);
        this.f                      = f;
        this.constraints            = normalize_constraints(constraints);
        this.restrict_to_non_negative  = restrict_to_non_negative;
        this.column_labels           = Array_list<>();
        this.epsilon                = epsilon;
        this.max_ulps                = max_ulps;
        this.num_decision_variables   = f.get_coefficients().get_dimension() + (restrict_to_non_negative ? 0 : 1);
        this.num_slack_variables      = get_constraint_type_counts(Relationship.LEQ) +
                                      get_constraint_type_counts(Relationship.GEQ);
        this.num_artificial_variables = get_constraint_type_counts(Relationship.EQ) +
                                      get_constraint_type_counts(Relationship.GEQ);
        this.tableau = create_tableau(goal_type == Goal_Type.MAXIMIZE);
        // initialize the basic variables for phase 1:
        //   we know that only slack or artificial variables can be basic
        initialize_basic_variables(get_slack_variable_offset());
        initialize_column_labels();
    }

    /**
     * Checks that the dimensions of the objective function and the constraints match.
     * @param objective_function the objective function
     * @param c the set of constraints
     * @ if the constraint dimensions do not match with the
     *   dimension of the objective function
     */
    private void check_dimensions(const Linear_Objective_Function objective_function, const Collection<Linear_Constraint> c) 
    {
        const int dimension = objective_function.get_coefficients().get_dimension();
        for (const Linear_Constraint constraint : c) 
        {
            const int constraint_dimension = constraint.get_coefficients().get_dimension();
            if (constraint_dimension != dimension) 
            {
                throw (Localized_Core_Formats.DIMENSIONS_MISMATCH, constraint_dimension, dimension);
            }
        }
    }
    /**
     * Initialize the labels for the columns.
     */
    protected void initialize_column_labels() 
    {
      if (get_num_objective_functions() == 2) 
      {
        column_labels.add("W");
      }
      column_labels.add("Z");
      for (int i{}; i < get_original_num_decision_variables(); i++) 
      {
        column_labels.add("x" + i);
      }
      if (!restrict_to_non_negative) 
      {
        column_labels.add(NEGATIVE_VAR_COLUMN_LABEL);
      }
      for (int i{}; i < get_num_slack_variables(); i++) 
      {
        column_labels.add("s" + i);
      }
      for (int i{}; i < get_num_artificial_variables(); i++) 
      {
        column_labels.add("a" + i);
      }
      column_labels.add("RHS");
    }

    /**
     * Create the tableau by itself.
     * @param maximize if true, goal is to maximize the objective function
     * @return created tableau
     */
    protected Array_2D_Row_Real_Matrix create_tableau(const bool maximize) 
    {

        // create a matrix of the correct size
        int width = num_decision_variables + num_slack_variables +
        num_artificial_variables + get_num_objective_functions() + 1; // + 1 is for RHS
        int height = constraints.size() + get_num_objective_functions();
        Array_2D_Row_Real_Matrix matrix = Array_2D_Row_Real_Matrix(height, width);

        // initialize the objective function rows
        if (get_num_objective_functions() == 2) 
        {
            matrix.set_entry(0, 0, -1);
        }

        int z_index = (get_num_objective_functions() == 1) ? 0 : 1;
        matrix.set_entry(z_index, z_index, maximize ? 1 : -1);
        Real_Vector objective_coefficients = maximize ? f.get_coefficients().map_multiply(-1) : f.get_coefficients();
        copy_array(objective_coefficients.to_array(), matrix.get_data_ref()[z_index]);
        matrix.set_entry(z_index, width - 1, maximize ? f.get_constant_term() : -1 * f.get_constant_term());

        if (!restrict_to_non_negative) 
        {
            matrix.set_entry(z_index, get_slack_variable_offset() - 1, get_inverted_coefficient_sum(objective_coefficients));
        }

        // initialize the constraint rows
        int slack_var = 0;
        int artificial_var = 0;
        for (int i{}; i < constraints.size(); i++) 
        {
            Linear_Constraint constraint = constraints.get(i);
            int row = get_num_objective_functions() + i;

            // decision variable coefficients
            copy_array(constraint.get_coefficients().to_array(), matrix.get_data_ref()[row]);

            // x-
            if (!restrict_to_non_negative) 
            {
                matrix.set_entry(row, get_slack_variable_offset() - 1, get_inverted_coefficient_sum(constraint.get_coefficients()));
            }

            // RHS
            matrix.set_entry(row, width - 1, constraint.get_value());

            // slack variables
            if (constraint.get_relationship() == Relationship.LEQ) 
            {
                matrix.set_entry(row, get_slack_variable_offset() + slack_var++, 1);  // slack
            }
else if (constraint.get_relationship() == Relationship.GEQ) 
            {
                matrix.set_entry(row, get_slack_variable_offset() + slack_var++, -1); // excess
            }

            // artificial variables
            if ((constraint.get_relationship() == Relationship.EQ) ||
                (constraint.get_relationship() == Relationship.GEQ)) 
                {
                matrix.set_entry(0, get_artificial_variable_offset() + artificial_var, 1);
                matrix.set_entry(row, get_artificial_variable_offset() + artificial_var++, 1);
                matrix.set_row_vector(0, matrix.get_row_vector(0).subtract(matrix.get_row_vector(row)));
            }
        }

        return matrix;
    }

    /**
     * Get versions of the constraints which have positive right hand sides.
     * @param original_constraints original (not normalized) constraints
     * @return versions of the constraints
     */
    public List<Linear_Constraint> normalize_constraints(Collection<Linear_Constraint> original_constraints) 
    {
        List<Linear_Constraint> normalized = Array_list<>(original_constraints.size());
        for (Linear_Constraint constraint : original_constraints) 
        {
            normalized.add(normalize(constraint));
        }
        return normalized;
    }

    /**
     * Get a equation equivalent to this one with a positive right hand side.
     * @param constraint reference constraint
     * @return equation
     */
    private Linear_Constraint normalize(const Linear_Constraint constraint) 
    {
        if (constraint.get_value() < 0) 
        {
            return Linear_Constraint(constraint.get_coefficients().map_multiply(-1), constraint.get_relationship().opposite_relationship(), -1 * constraint.get_value());
        }
        return Linear_Constraint(constraint.get_coefficients(), constraint.get_relationship(), constraint.get_value());
    }

    /**
     * Get the number of objective functions in this tableau.
     * @return 2 for Phase 1.  1 for Phase 2.
     */
    protected const int get_num_objective_functions() 
    {
        return this.num_artificial_variables > 0 ? 2 : 1;
    }

    /**
     * Get a count of constraints corresponding to a specified relationship.
     * @param relationship relationship to count
     * @return number of constraint with the specified relationship
     */
    private int get_constraint_type_counts(const Relationship relationship) 
    {
        int count{};
        for (const Linear_Constraint constraint : constraints) 
        {
            if (constraint.get_relationship() == relationship) 
            {
                ++count;
            }
        }
        return count;
    }

    /**
     * Get the -1 times the sum of all coefficients in the given array.
     * @param coefficients coefficients to sum
     * @return the -1 times the sum of all coefficients in the given array.
     */
    protected static double get_inverted_coefficient_sum(const Real_Vector coefficients) 
    {
        double sum{};
        for (double coefficient : coefficients.to_array()) 
        {
            sum -= coefficient;
        }
        return sum;
    }

    /**
     * Checks whether the given column is basic.
     * @param col index of the column to check
     * @return the row that the variable is basic in.  NULL if the column is not basic
     */
    protected Integer get_basic_row(const int col) 
    {
        const int row = basic_variables[col];
        return row == -1 ? NULL : row;
    }

    /**
     * Returns the variable that is basic in this row.
     * @param row the index of the row to check
     * @return the variable that is basic for this row.
     */
    protected int get_basic_variable(const int row) 
    {
        return basic_rows[row];
    }

    /**
     * Initializes the basic variable / row mapping.
     * @param start_column the column to start
     */
    private void initialize_basic_variables(const int& start_column) 
    {
        basic_variables = int[get_width() - 1];
        basic_rows = int[get_height()];

        Arrays.fill(basic_variables, -1);

        for (int i = start_column; i < get_width() - 1; i++) 
        {
            Integer row = find_basic_row(i);
            if (row != NULL) 
            {
                basic_variables[i] = row;
                basic_rows[row] = i;
            }
        }
    }

    /**
     * Returns the row in which the given column is basic.
     * @param col index of the column
     * @return the row that the variable is basic in, or {@code NULL} if the variable is not basic.
     */
    private Integer find_basic_row(const int col) 
    {
        Integer row = NULL;
        for (int i{}; i < get_height(); i++) 
        {
            const double entry = get_entry(i, col);
            if (Precision.equals(entry, 1d, max_ulps) && (row == NULL)) 
            {
                row = i;
            }
else if (!Precision.equals(entry, 0d, max_ulps)) 
            {
                return NULL;
            }
        }
        return row;
    }

    /**
     * Removes the phase 1 objective function, positive cost non-artificial variables, * and the non-basic artificial variables from this tableau.
     */
    protected void drop_phase1_objective() 
    {
        if (get_num_objective_functions() == 1) 
        {
            return;
        }

        const Set<Integer> columns_to_drop = Tree_Set<>();
        columns_to_drop.add(0);

        // positive cost non-artificial variables
        for (int i = get_num_objective_functions(); i < get_artificial_variable_offset(); i++) 
        {
            const double entry = get_entry(0, i);
            if (Precision.compare_to(entry, 0d, epsilon) > 0) 
            {
                columns_to_drop.add(i);
            }
        }

        // non-basic artificial variables
        for (int i{}; i < get_num_artificial_variables(); i++) 
        {
            int col = i + get_artificial_variable_offset();
            if (get_basic_row(col) == NULL) 
            {
                columns_to_drop.add(col);
            }
        }

        const std::vector<std::vector<double>> matrix = std::vector<double>(get_height() - 1][get_width() - columns_to_drop.size()];
        for (int i{ 1 }; i < get_height(); i++) 
        {
            int col = 0;
            for (int j{}; j < get_width(); j++) 
            {
                if (!columns_to_drop.contains(j)) 
                {
                    matrix[i - 1][col++] = get_entry(i, j);
                }
            }
        }

        // remove the columns in reverse order so the indices are correct
        Integer[] drop = columns_to_drop.to_array(new Integer[0]);
        for (int i = drop.size() - 1; i >= 0; i--) 
        {
            column_labels.remove(static_cast<int>( drop[i]);
        }

        this.tableau = Array_2D_Row_Real_Matrix(matrix);
        this.num_artificial_variables = 0;
        // need to update the basic variable mappings as row/columns have been dropped
        initialize_basic_variables(get_num_objective_functions());
    }

    /**
     * @param src the source array
     * @param dest the destination array
     */
    private void copy_array(const std::vector<double> src, const std::vector<double> dest) 
    {
        System.arraycopy(src, 0, dest, get_num_objective_functions(), src.size());
    }

    /**
     * Returns whether the problem is at an optimal state.
     * @return whether the model has been solved
     */
    bool is_optimal() 
    {
        const std::vector<double> objective_function_row = get_row(0);
        const int end = get_rhs_offset();
        for (int i = get_num_objective_functions(); i < end; i++) 
        {
            const double entry = objective_function_row[i];
            if (Precision.compare_to(entry, 0d, epsilon) < 0) 
            {
                return false;
            }
        }
        return true;
    }

    /**
     * Get the current solution.
     * @return current solution
     */
    protected Point_valuePair get_solution() 
    {
        int negative_var_column = column_labels.index_of(NEGATIVE_VAR_COLUMN_LABEL);
        Integer negative_var_basic_row = negative_var_column > 0 ? get_basic_row(negative_var_column) : NULL;
        double most_negative = negative_var_basic_row == NULL ? 0 : get_entry(negative_var_basic_row, get_rhs_offset());

        const Set<Integer> used_basic_rows = Hash_Set<>();
        const std::vector<double> coefficients = std::vector<double>(get_original_num_decision_variables()];
        for (int i{}; i < coefficients.size(); i++) 
        {
            int col_index = column_labels.index_of("x" + i);
            if (col_index < 0) 
            {
                coefficients[i] = 0;
                continue;
            }
            Integer basic_row = get_basic_row(col_index);
            if (basic_row != NULL && basic_row == 0) 
            {
                // if the basic row is found to be the objective function row
                // set the coefficient to 0 -> this case handles unconstrained
                // variables that are still part of the objective function
                coefficients[i] = 0;
            }
else if (used_basic_rows.contains(basic_row)) 
            {
                // if multiple variables can take a given value
                // then we choose the first and set the rest equal to 0
                coefficients[i] = 0 - (restrict_to_non_negative ? 0 : most_negative);
            }
else 
            {
                used_basic_rows.add(basic_row);
                coefficients[i] =
                    (basic_row == NULL ? 0 : get_entry(basic_row, get_rhs_offset())) -
                    (restrict_to_non_negative ? 0 : most_negative);
            }
        }
        return Point_valuePair(coefficients, f.value(coefficients));
    }

    /**
     * Perform the row operations of the simplex algorithm with the selected
     * pivot column and row.
     * @param pivot_col the pivot column
     * @param pivot_row the pivot row
     */
    protected void perform_row_operations(const int& pivot_col, int pivot_row) 
    {
        // set the pivot element to 1
        const double pivot_val = get_entry(pivot_row, pivot_col);
        divide_row(pivot_row, pivot_val);

        // set the rest of the pivot column to 0
        for (int i{}; i < get_height(); i++) 
        {
            if (i != pivot_row) 
            {
                const double multiplier = get_entry(i, pivot_col);
                if (multiplier != 0.0) 
                {
                    subtract_row(i, pivot_row, multiplier);
                }
            }
        }

        // update the basic variable mappings
        const int previous_basic_variable = get_basic_variable(pivot_row);
        basic_variables[previous_basic_variable] = -1;
        basic_variables[pivot_col] = pivot_row;
        basic_rows[pivot_row] = pivot_col;
    }

    /**
     * Divides one row by a given divisor.
     * <p>
     * After application of this operation, the following will hold:
     * <pre>dividend_row = dividend_row / divisor</pre>
     *
     * @param dividend_row_index index of the row
     * @param divisor value of the divisor
     */
    protected void divide_row(const int dividend_row_index, const double divisor) 
    {
        const std::vector<double> dividend_row = get_row(dividend_row_index);
        for (int j{}; j < get_width(); j++) 
        {
            dividend_row[j] /= divisor;
        }
    }

    /**
     * Subtracts a multiple of one row from another.
     * <p>
     * After application of this operation, the following will hold:
     * <pre>minuend_row = minuend_row - multiple * subtrahend_row</pre>
     *
     * @param minuend_row_index row index
     * @param subtrahend_row_index row index
     * @param multiplier multiplication factor
     */
    protected void subtract_row(const int minuend_row_index, const int subtrahend_row_index, const double multiplier) 
    {
        const std::vector<double> minuend_row = get_row(minuend_row_index);
        const std::vector<double> subtrahend_row = get_row(subtrahend_row_index);
        for (int i{}; i < get_width(); i++) 
        {
            minuend_row[i] -= subtrahend_row[i] * multiplier;
        }
    }

    /**
     * Get the width of the tableau.
     * @return width of the tableau
     */
    protected const int get_width() 
    {
        return tableau.get_column_dimension();
    }

    /**
     * Get the height of the tableau.
     * @return height of the tableau
     */
    protected const int get_height() 
    {
        return tableau.get_row_dimension();
    }

    /**
     * Get an entry of the tableau.
     * @param row row index
     * @param column column index
     * @return entry at (row, column)
     */
    protected const double get_entry(const int& row, const int column) 
    {
        return tableau.get_entry(row, column);
    }

    /**
     * Set an entry of the tableau.
     * @param row row index
     * @param column column index
     * @param value for the entry
     */
    protected const void set_entry(const int& row, const int& column, const double value) 
    {
        tableau.set_entry(row, column, value);
    }

    /**
     * Get the offset of the first slack variable.
     * @return offset of the first slack variable
     */
    protected const int get_slack_variable_offset() 
    {
        return get_num_objective_functions() + num_decision_variables;
    }

    /**
     * Get the offset of the first artificial variable.
     * @return offset of the first artificial variable
     */
    protected const int get_artificial_variable_offset() 
    {
        return get_num_objective_functions() + num_decision_variables + num_slack_variables;
    }

    /**
     * Get the offset of the right hand side.
     * @return offset of the right hand side
     */
    protected const int get_rhs_offset() 
    {
        return get_width() - 1;
    }

    /**
     * Get the number of decision variables.
     * <p>
     * If variables are not restricted to positive values, this will include 1 extra decision variable to represent
     * the absolute value of the most negative variable.
     *
     * @return number of decision variables
     * @see #get_original_num_decision_variables()
     */
    protected const int get_num_decision_variables() 
    {
        return num_decision_variables;
    }

    /**
     * Get the original number of decision variables.
     * @return original number of decision variables
     * @see #get_num_decision_variables()
     */
    protected const int get_original_num_decision_variables() 
    {
        return f.get_coefficients().get_dimension();
    }

    /**
     * Get the number of slack variables.
     * @return number of slack variables
     */
    protected const int get_num_slack_variables() 
    {
        return num_slack_variables;
    }

    /**
     * Get the number of artificial variables.
     * @return number of artificial variables
     */
    protected const int get_num_artificial_variables() 
    {
        return num_artificial_variables;
    }

    /**
     * Get the row from the tableau.
     * @param row the row index
     * @return the reference to the underlying row data
     */
    protected const std::vector<double> get_row(const int& row) 
    {
        return tableau.get_data_ref()[row];
    }

    /**
     * Get the tableau data.
     * @return tableau data
     */
    protected const std::vector<std::vector<double>> get_data() 
    {
        return tableau.get_data();
    }

    /** {@inherit_doc} */
    //override
    public bool equals(Object other) 
    {

      if (this == other) 
      {
        return true;
      }

      if (other instanceof Simplex_Tableau) 
      {
          Simplex_Tableau rhs = (Simplex_Tableau) other;
          return (restrict_to_non_negative  == rhs.restrict_to_non_negative) &&
                 (num_decision_variables   == rhs.num_decision_variables) &&
                 (num_slack_variables      == rhs.num_slack_variables) &&
                 (num_artificial_variables == rhs.num_artificial_variables) &&
                 (epsilon                == rhs.epsilon) &&
                 (max_ulps                == rhs.max_ulps) &&
                 f.equals(rhs.f) &&
                 constraints.equals(rhs.constraints) &&
                 tableau.equals(rhs.tableau);
      }
      return false;
    }

    /** {@inherit_doc} */
    //override
    public int hash_code() 
    {
        return Boolean.value_of(restrict_to_non_negative).hash_code() ^
               num_decision_variables ^
               num_slack_variables ^
               num_artificial_variables ^
               static_cast<double>(epsilon).hash_code() ^
               max_ulps ^
               f.hash_code() ^
               constraints.hash_code() ^
               tableau.hash_code();
    }

    /**
     * Serialize the instance.
     * @param oos stream where object should be written
     * @IOException if object cannot be written to stream
     */
    private void write_object(Object_Output_Stream oos)
        IOException 
        {
        oos.default_write_object();
        Matrix_Utils.serialize_real__matrix(tableau, oos);
    }

    /**
     * Deserialize the instance.
     * @param ois stream from which the object should be read
     * @Class_Not_Found_Exception if a class in the stream cannot be found
     * @IOException if object cannot be read from the stream
     */
    private void read_object(Object_Input_Stream ois)
      Class_Not_Found_Exception, IOException 
      {
        ois.default_read_object();
        Matrix_Utils.deserialize_real__matrix(this, "tableau", ois);
    }
}


