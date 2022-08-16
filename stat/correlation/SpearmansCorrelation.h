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

  //package org.hipparchus.stat.correlation;

  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.linear.Block_Real_Matrix;
  //import org.hipparchus.linear.Real_Matrix;
  //import org.hipparchus.stat.Localized_Stat_Formats;
  //import org.hipparchus.stat.ranking.NaN_Strategy;
  //import org.hipparchus.stat.ranking.Natural_Ranking;
  //import org.hipparchus.stat.ranking.Ranking_Algorithm;
  //import org.hipparchus.util.Math_Arrays;
#include <vector>
#include "../../core/linear/RealMatrix.h"
#include "../ranking/RankingAlgorithm.h"
#include "PearsonsCorrelation.h"
#include "../ranking/NaturalRanking.h"
#include "../../core/linear/BlockRealMatrix.h"
#include "../../core/util/MathArrays.h"

/**
 * Spearman's rank correlation. This implementation performs a rank
 * transformation on the input data and then computes {@link Pearsons_Correlation}
 * on the ranked data.
 * <p>
 * By default, ranks are computed using {@link Natural_Ranking} with default
 * strategies for handling NaNs and ties in the data (NaNs maximal, ties averaged).
 * The ranking algorithm can be set using a constructor argument.
 */
class Spearmans_Correlation
{
private:
	/** Input data */
	const Real_Matrix my_data;

	/** Ranking algorithm  */
	const Ranking_Algorithm my_ranking_algorithm;

	/** Rank correlation */
	Pearsons_Correlation my_rank_correlation;

	/**
	 * Applies rank transform to each of the columns of <code>matrix</code>
	 * using the current <code>ranking_algorithm</code>.
	 *
	 * @param matrix matrix to transform
	 * @return a rank-transformed matrix
	 */
	Real_Matrix rank_transform(const Real_Matrix& matrix)
	{
		Real_Matrix transformed = matrix;
		for (int i{}; i < transformed.get_column_dimension(); i++)
		{
			transformed.set_column(i, my_ranking_algorithm.rank(transformed.get_column(i)));
		}

		return transformed;
	}

public:
	/**
	 * Create a Spearmans_Correlation without data.
	 */
	Spearmans_Correlation()
	{
		Spearmans_Correlation(Natural_Ranking());
	}

	/**
	 * Create a Spearmans_Correlation with the given ranking algorithm.
	 *
	 * @param ranking_algorithm ranking algorithm
	 * @ if the provided {@link Ranking_Algorithm} is of
	 * type {@link Natural_Ranking} and uses a {@link NaN_Strategy#REMOVED} strategy
	 */
	Spearmans_Correlation(const Ranking_Algorithm& ranking_algorithm)
	{
		if (ranking_algorithm instanceof Natural_Ranking && NaN_Strategy::REMOVED == ((Natural_Ranking)ranking_algorithm).get_nan_strategy())
		{
			throw std::exception("not implemented");
			//throw (Localized_Stat_Formats::NOT_SUPPORTED_NAN_STRATEGY, NaN_Strategy::REMOVED);
		}

		my_data = NULL;
		my_ranking_algorithm = ranking_algorithm;
		my_rank_correlation = NULL;
	}

	/**
	 * Create a Spearmans_Correlation from the given data matrix.
	 *
	 * @param data_matrix matrix of data with columns representing
	 * variables to correlate
	 */
	Spearmans_Correlation(const Real_Matrix& data_matrix)
	{
		Spearmans_Correlation(data_matrix, Natural_Ranking());
	}

	/**
	 * Create a Spearmans_Correlation with the given input data matrix
	 * and ranking algorithm.
	 *
	 * @param data_matrix matrix of data with columns representing
	 * variables to correlate
	 * @param ranking_algorithm ranking algorithm
	 * @ if the provided {@link Ranking_Algorithm} is of
	 * type {@link Natural_Ranking} and uses a {@link NaN_Strategy#REMOVED} strategy
	 */
	Spearmans_Correlation(const Real_Matrix& data_matrix, const Ranking_Algorithm& ranking_algorithm)
	{
		if (ranking_algorithm instanceof Natural_Ranking && NaN_Strategy::REMOVED == ((Natural_Ranking)ranking_algorithm).get_nan_strategy())
		{
			throw std::exception("not implemented");
			//throw (Localized_Stat_Formats::NOT_SUPPORTED_NAN_STRATEGY, NaN_Strategy::REMOVED);
		}

		my_ranking_algorithm = ranking_algorithm;
		my_data = rank_transform(data_matrix);
		my_rank_correlation = Pearsons_Correlation(my_data);
	}

	/**
	 * Calculate the Spearman Rank Correlation Matrix.
	 *
	 * @return Spearman Rank Correlation Matrix
	 * @Null_Pointer_Exception if this instance was created with no data
	 */
	Real_Matrix get_correlation_matrix()
	{
		return my_rank_correlation.get_correlation_matrix();
	}

	/**
	 * Returns a {@link Pearsons_Correlation} instance constructed from the
	 * ranked input data. That is, * <code>new Spearmans_Correlation(matrix).get_rank_correlation()</code>
	 * is equivalent to
	 * <code>new Pearsons_Correlation(rank_transform(matrix))</code> where
	 * <code>rank_transform(matrix)</code> is the result of applying the
	 * configured <code>Ranking_Algorithm</code> to each of the columns of
	 * <code>matrix.</code>
	 *
	 * <p>Returns NULL if this instance was created with no data.</p>
	 *
	 * @return Pearsons_Correlation among ranked column data
	 */
	Pearsons_Correlation get_rank_correlation() const
	{
		return my_rank_correlation;
	}

	/**
	 * Computes the Spearman's rank correlation matrix for the columns of the
	 * input matrix.
	 *
	 * @param matrix matrix with columns representing variables to correlate
	 * @return correlation matrix
	 */
	Real_Matrix compute_correlation_matrix(const Real_Matrix& matrix)
	{
		const Real_Matrix matrix_copy = rank_transform(matrix);
		return Pearsons_Correlation().compute_correlation_matrix(matrix_copy);
	}

	/**
	 * Computes the Spearman's rank correlation matrix for the columns of the
	 * input rectangular array.  The columns of the array represent values
	 * of variables to be correlated.
	 *
	 * @param matrix matrix with columns representing variables to correlate
	 * @return correlation matrix
	 */
	Real_Matrix compute_correlation_matrix(const std::vector<std::vector<double>>& matrix)
	{
		return compute_correlation_matrix(Block_Real_Matrix(matrix));
	}

	/**
	 * Computes the Spearman's rank correlation coefficient between the two arrays.
	 *
	 * @param x_array first data array
	 * @param y_array second data array
	 * @return Returns Spearman's rank correlation coefficient for the two arrays
	 * @ if the arrays lengths do not match
	 * @ if the array length is less than 2
	 */
	double correlation(const std::vector<double>& x_array, const std::vector<double>& y_array)
	{
		Math_Arrays::check_equal_length(x_array, y_array);
		if (x_array.size() < 2)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::INSUFFICIENT_DIMENSION, x_array.size(), 2);
		}

		return Pearsons_Correlation().correlation(my_ranking_algorithm.rank(x_array), my_ranking_algorithm.rank(y_array));
	}
};