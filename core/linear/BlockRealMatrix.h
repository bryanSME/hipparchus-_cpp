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

  //package org.hipparchus.linear;

  //import java.io.Serializable;
  //import java.util.Arrays;
#include <vector>
#include <algorithm>

#include  "MatrixUtils.h"
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Utils;

/**
 * Cache-friendly implementation of Real_Matrix using a flat arrays to store
 * square blocks of the matrix.
 * <p>
 * This implementation is specially designed to be cache-friendly. Square blocks are
 * stored as small arrays and allow efficient traversal of data both in row major direction
 * and columns major direction, one block at a time. This greatly increases performances
 * for algorithms that use crossed directions loops like multiplication or transposition.
 * </p>
 * <p>
 * The size of square blocks is a static parameter. It may be tuned according to the cache
 * size of the target computer processor. As a rule of thumbs, it should be the largest
 * value that allows three blocks to be simultaneously cached (this is necessary for example
 * for matrix multiplication). The default value is to use 52x52 blocks which is well suited
 * for processors with 64k L1 cache (one block holds 2704 values or 21632 bytes). This value
 * could be lowered to 36x36 for processors with 32k L1 cache.
 * </p>
 * <p>
 * The regular blocks represent {@link #BLOCK_SIZE} x {@link #BLOCK_SIZE} squares. Blocks
 * at right hand side and bottom side may be smaller to fit matrix dimensions. The square
 * blocks are flattened in row major order in single dimension arrays which are therefore
 * {@link #BLOCK_SIZE}<sup>2</sup> elements long for regular blocks. The blocks are themselves
 * organized in row major order.
 * </p>
 * <p>
 * As an example, for a block size of 52x52, a 100x60 matrix would be stored in 4 blocks.
 * Block 0 would be a {@code double[2704]} array holding the upper left 52x52 square, block 1
 * would be a {@code double[416]} array holding the upper right 52x8 rectangle, block 2 would
 * be a {@code double[2496]} array holding the lower left 48x52 rectangle and block 3 would
 * be a {@code double[384]} array holding the lower right 48x8 rectangle.
 * </p>
 * <p>
 * The layout complexity overhead versus simple mapping of matrices to java
 * arrays is negligible for small matrices (about 1%). The gain from cache efficiency leads
 * to up to 3-fold improvements for matrices of moderate to large size.
 * </p>
 */
class Block_Real_Matrix : public Abstract_Real_Matrix
{
private:
	/** Blocks of matrix entries. */
	std::vector<std::vector<double>> my_blocks;
	/** Number of rows of the matrix. */
	const int my_rows;
	/** Number of columns of the matrix. */
	const int my_columns;
	/** Number of block rows of the matrix. */
	const int my_block_rows;
	/** Number of block columns of the matrix. */
	const int my_block_columns;

	/**
	 * Get the height of a block.
	 * @param block_row row index (in block sense) of the block
	 * @return height (number of rows) of the block
	 */
	int block_height(const int& block_row)
	{
		return (block_row == my_block_rows - 1)
			? my_rows - block_row * BLOCK_SIZE
			: BLOCK_SIZE;
	};

	/**
	 * Get the width of a block.
	 * @param block_column column index (in block sense) of the block
	 * @return width (number of columns) of the block
	 */
	int block_width(const int& block_column)
	{
		return (block_column == my_block_columns - 1)
			? my_columns - block_column * BLOCK_SIZE
			: BLOCK_SIZE;
	};

	/**
	 * Copy a part of a block into another one
	 * <p>This method can be called only when the specified part fits in both
	 * blocks, no verification is done here.</p>
	 * @param src_block source block
	 * @param src_width source block width ({@link #BLOCK_SIZE} or smaller)
	 * @param src_start_row start row in the source block
	 * @param src_end_row end row (exclusive) in the source block
	 * @param src_start_column start column in the source block
	 * @param src_end_column end column (exclusive) in the source block
	 * @param dst_block destination block
	 * @param dst_width destination block width ({@link #BLOCK_SIZE} or smaller)
	 * @param dst_start_row start row in the destination block
	 * @param dst_start_column start column in the destination block
	 */
	void copy_block_part(const std::vector<double>& src_block, const int& src_width, const int& src_start_row, const int& src_end_row, const int& src_start_column, const int& src_end_column, const std::vector<double>& dst_block, const int& dst_width, const int& dst_start_row, const int& dst_start_column)
	{
		const int length = src_end_column - src_start_column;
		int src_pos = src_start_row * src_width + src_start_column;
		int dst_pos = dst_start_row * dst_width + dst_start_column;
		for (int src_row = src_start_row; src_row < src_end_row; ++src_row)
		{
			System.arraycopy(src_block, src_pos, dst_block, dst_pos, length);
			src_pos += src_width;
			dst_pos += dst_width;
		}
	};

public:
	/** Block size. */
	static constexpr int BLOCK_SIZE = 52;
	/**
	 * Create a matrix with the supplied row and column dimensions.
	 *
	 * @param rows  the number of rows in the matrix
	 * @param columns  the number of columns in the matrix
	 * @ if row or column dimension is not
	 * positive.
	 */
	Block_Real_Matrix(const int rows, const int columns)
		:
		my_rows{ rows },
		my_columns{ columns },
		my_block_rows{ (rows + BLOCK_SIZE - 1) / BLOCK_SIZE },
		my_block_columns{ (columns + BLOCK_SIZE - 1) / BLOCK_SIZE }
	{
		super(my_rows, my_columns);

		// allocate storage blocks, taking care of smaller ones at right and bottom
		my_blocks = create_blocks_layout(my_rows, my_columns);
	};

	/**
	 * Create a dense matrix copying entries from raw layout data.
	 * <p>The input array <em>must</em> already be in raw layout.</p>
	 * <p>Calling this constructor is equivalent to call:
	 * <pre>matrix = Block_Real_Matrix(raw_data.size(), raw_data[0].size(), *                                   to_blocks_layout(raw_data), false);</pre>
	 * </p>
	 *
	 * @param raw_data data for matrix, in raw layout
	 * @ if the shape of {@code block_data} is
	 * inconsistent with block layout.
	 * @ if row or column dimension is not
	 * positive.
	 * @see #Block_Real_Matrix(int, int, std::vector<std::vector<double>>, bool)
	 */
	Block_Real_Matrix(const std::vector<std::vector<double>>& raw_data)
	{
		Block_Real_Matrix(raw_data.size(), raw_data[0].size(), to_blocks_layout(raw_data), false);
	}

	/**
	 * Create a dense matrix copying entries from block layout data.
	 * <p>The input array <em>must</em> already be in blocks layout.</p>
	 *
	 * @param rows Number of rows in the matrix.
	 * @param columns Number of columns in the matrix.
	 * @param block_data data for matrix
	 * @param copy_array Whether the input array will be copied or referenced.
	 * @ if the shape of {@code block_data} is
	 * inconsistent with block layout.
	 * @ if row or column dimension is not
	 * positive.
	 * @see #create_blocks_layout(int, int)
	 * @see #to_blocks_layout(std::vector<std::vector<double>>)
	 * @see #Block_Real_Matrix(std::vector<std::vector<double>>)
	 */
	Block_Real_Matrix(const int& rows, const int& columns, const std::vector<std::vector<double>>& block_data, const bool copy_array) // NOPMD - array copy is taken care of by parameter
		:
		// number of blocks
		my_block_rows{ (rows + BLOCK_SIZE - 1) / BLOCK_SIZE },
		my_block_columns{ (columns + BLOCK_SIZE - 1) / BLOCK_SIZE },
		my_rows{ rows },
		my_columns{ columns }
	{
		super(my_rows, my_columns);

		if (copy_array)
		{
			// allocate storage blocks, taking care of smaller ones at right and bottom
			my_blocks = std::vector<std::vector<double>>(my_block_rows * my_block_columns);
		}
		else
		{
			// reference existing array
			my_blocks = block_data;
		}

		int index{};
		for (int i_block{}; i_block < my_block_rows; ++i_block)
		{
			const int i_height = block_height(i_block);
			for (int j_block{}; j_block < my_block_columns; ++j_block, ++index)
			{
				if (block_data[index].size() != i_height * block_width(j_block))
				{
					throw std::exception("not implemented");
					//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, block_data[index].size(), i_height * block_width(j_block));
				}
				if (copy_array)
				{
					blocks[index] = block_data[index].clone();
				}
			}
		}
	}

	/**
	 * Convert a data array from raw layout to blocks layout.
	 * <p>
	 * Raw layout is the straightforward layout where element at row i and
	 * column j is in array element <code>raw_data[i][j]</code>. Blocks layout
	 * is the layout used in {@link Block_Real_Matrix} instances, where the matrix
	 * is split in square blocks (except at right and bottom side where blocks may
	 * be rectangular to fit matrix size) and each block is stored in a flattened
	 * one-dimensional array.
	 * </p>
	 * <p>
	 * This method creates an array in blocks layout from an input array in raw layout.
	 * It can be used to provide the array argument of the {@link
	 * #Block_Real_Matrix(int, int, std::vector<std::vector<double>>, bool)} constructor.
	 * </p>
	 * @param raw_data Data array in raw layout.
	 * @return a data array containing the same entries but in blocks layout.
	 * @ if {@code raw_data} is not rectangular.
	 * @see #create_blocks_layout(int, int)
	 * @see #Block_Real_Matrix(int, int, std::vector<std::vector<double>>, bool)
	 */
	static std::vector<std::vector<double>> to_blocks_layout(const std::vector<std::vector<double>> raw_data)

	{
		const int rows = raw_data.size();
		const int columns = raw_data[0].size();
		const int& my_block_rows = (rows + BLOCK_SIZE - 1) / BLOCK_SIZE;
		const int& my_block_columns = (columns + BLOCK_SIZE - 1) / BLOCK_SIZE;

		// safety checks
		for (int i{}; i < raw_data.size(); ++i)
		{
			const int length = raw_data[i].size();
			if (length != columns)
			{
				throw std::exception("not implemented");
				// throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, columns, length);
			}
		}

		// convert array
		auto blocks = std::vector<std::vector<double>>(my_block_rows * my_block_columns, std::vector<double>());
		int block_index{};
		for (int i_block{}; i_block < my_block_rows; ++i_block)
		{
			const int p_start{ i_block * BLOCK_SIZE };
			const int p_end{ std::min(p_start + BLOCK_SIZE, my_rows) };
			const int i_height = p_end - p_start;
			for (int j_block{}; j_block < my_block_columns; ++j_block)
			{
				const int q_start{ j_block * BLOCK_SIZE };
				const int q_end{ std::min(q_start + BLOCK_SIZE, columns) };
				const int j_width = q_end - q_start;

				// allocate block
				const auto block = std::vector<double>(i_height * j_width);
				blocks[block_index] = block;

				// copy data
				int index{};
				for (int p{ p_start }; p < p_end; ++p)
				{
					System.arraycopy(raw_data[p], q_start, block, index, j_width);
					index += j_width;
				}
				++block_index;
			}
		}

		return blocks;
	}

	/**
	 * Create a data array in blocks layout.
	 * <p>
	 * This method can be used to create the array argument of the {@link
	 * #Block_Real_Matrix(int, int, std::vector<std::vector<double>>, bool)} constructor.
	 * </p>
	 * @param rows Number of rows in the matrix.
	 * @param columns Number of columns in the matrix.
	 * @return a data array in blocks layout.
	 * @see #to_blocks_layout(std::vector<std::vector<double>>)
	 * @see #Block_Real_Matrix(int, int, std::vector<std::vector<double>>, bool)
	 */
	static std::vector<std::vector<double>> create_blocks_layout(const int rows, const int columns)
	{
		const int& my_block_rows = (rows + BLOCK_SIZE - 1) / BLOCK_SIZE;
		const int& my_block_columns = (columns + BLOCK_SIZE - 1) / BLOCK_SIZE;

		const std::vector<std::vector<double>> blocks = std::vector<double>(my_block_rows * my_block_columns][];
		int block_index{};
		for (int i_block{}; i_block < my_block_rows; ++i_block)
		{
			const int p_start{ i_block * BLOCK_SIZE };
			const int p_end{ std::min(p_start + BLOCK_SIZE, my_rows) };
			const int i_height = p_end - p_start;
			for (int j_block{}; j_block < my_block_columns; ++j_block)
			{
				const int q_start{ j_block * BLOCK_SIZE };
				const int q_end{ std::min(q_start + BLOCK_SIZE, columns) };
				const int j_width = q_end - q_start;
				blocks[block_index] = std::vector<double>(i_height * j_width];
				++block_index;
			}
		}

		return blocks;
	}

	/** {@inherit_doc} */
	//override
	Block_Real_Matrix create_matrix(const int& row_dimension, const int& column_dimension)

	{
		return Block_Real_Matrix(row_dimension, column_dimension);
	}

	/** {@inherit_doc} */
	//override
	Block_Real_Matrix copy()
	{
		// create an empty matrix
		Block_Real_Matrix copied = Block_Real_Matrix(my_rows, my_columns);

		// copy the blocks
		for (int i{}; i < blocks.size(); ++i)
		{
			System.arraycopy(blocks[i], 0, copied.blocks[i], 0, blocks[i].size());
		}

		return copied;
	}

	/** {@inherit_doc} */
	//override
	Block_Real_Matrix add(const Real_Matrix& m)

	{
		if (m instanceof Block_Real_Matrix)
		{
			return add((Block_Real_Matrix)m);
		}
		else
		{
			// safety check
			Matrix_Utils::check_addition_compatible(this, m);

			const auto out = Block_Real_Matrix(my_rows, my_columns);

			// perform addition block-wise, to ensure good cache behavior
			int block_index{};
			for (int i_block{}; i_block < out.my_block_rows; ++i_block)
			{
				const int p_start{ i_block * BLOCK_SIZE };
				const int p_end{ std::min(p_start + BLOCK_SIZE, my_rows) };
				for (int j_block{}; j_block < out.my_block_columns; ++j_block)
				{
					// perform addition on the current block
					auto out_block = out.blocks[block_index];
					const std::vector<double> t_block = my_blocks[block_index];
					const int q_start{ j_block * BLOCK_SIZE };
					const int q_end{ std::min(q_start + BLOCK_SIZE, columns) };
					int k{};
					for (int p{ p_start }; p < p_end; ++p)
					{
						for (int q{ q_start }; q < q_end; ++q)
						{
							out_block[k] = t_block[k] + m.get_entry(p, q);
							++k;
						}
					}
					// go to next block
					++block_index;
				}
			}

			return out;
		}
	}

	/**
	 * Compute the sum of this matrix and {@code m}.
	 *
	 * @param m Matrix to be added.
	 * @return {@code this} + m.
	 * @ if {@code m} is not the same
	 * size as this matrix.
	 */
	Block_Real_Matrix add(const Block_Real_Matrix m)

	{
		// safety check
		Matrix_Utils::check_addition_compatible(this, m);

		const auto out = Block_Real_Matrix(my_rows, my_columns);

		// perform addition block-wise, to ensure good cache behavior
		for (int block_index{}; block_index < out.blocks.size(); ++block_index)
		{
			std::vector<double> out_block = out.blocks[block_index];
			const auto t_block = my_blocks[block_index];
			const std::vector<double> m_block = m.blocks[block_index];
			for (int k{}; k < out_block.size(); ++k)
			{
				out_block[k] = t_block[k] + m_block[k];
			}
		}

		return out;
	}

	/** {@inherit_doc} */
	//override
	Block_Real_Matrix subtract(const Real_Matrix& m)

	{
		if (m instanceof Block_Real_Matrix)
		{
			return subtract((Block_Real_Matrix)m);
		}
		else
		{
			// safety check
			Matrix_Utils::check_subtraction_compatible(this, m);

			const auto out = Block_Real_Matrix(my_rows, my_columns);

			// perform subtraction block-wise, to ensure good cache behavior
			int block_index{};
			for (int i_block{}; i_block < out.my_block_rows; ++i_block)
			{
				const int p_start{ i_block * BLOCK_SIZE };
				const int p_end{ std::min(p_start + BLOCK_SIZE, my_rows) };
				for (int j_block{}; j_block < out.my_block_columns; ++j_block)
				{
					// perform subtraction on the current block
					auto out_block = out.blocks[block_index];
					const auto t_block = my_blocks[block_index];
					const int q_start{ j_block * BLOCK_SIZE };
					const int q_end{ std::min(q_start + BLOCK_SIZE, columns) };
					int k{};
					for (int p{ p_start }; p < p_end; ++p)
					{
						for (int q{ q_start }; q < q_end; ++q)
						{
							out_block[k] = t_block[k] - m.get_entry(p, q);
							++k;
						}
					}
					// go to next block
					++block_index;
				}
			}

			return out;
		}
	}

	/**
	 * Subtract {@code m} from this matrix.
	 *
	 * @param m Matrix to be subtracted.
	 * @return {@code this} - m.
	 * @ if {@code m} is not the
	 * same size as this matrix.
	 */
	Block_Real_Matrix subtract(const Block_Real_Matrix m)

	{
		// safety check
		Matrix_Utils::check_subtraction_compatible(this, m);

		const auto out = Block_Real_Matrix(my_rows, my_columns);

		// perform subtraction block-wise, to ensure good cache behavior
		for (int block_index{}; block_index < out.blocks.size(); ++block_index)
		{
			auto out_block = out.blocks[block_index];
			const auto t_block = my_blocks[block_index];
			const auto m_block = m.blocks[block_index];
			for (int k{}; k < out_block.size(); ++k)
			{
				out_block[k] = t_block[k] - m_block[k];
			}
		}

		return out;
	}

	/** {@inherit_doc} */
	//override
	Block_Real_Matrix scalar_add(const double d)
	{
		const auto out = Block_Real_Matrix(my_rows, my_columns);

		// perform subtraction block-wise, to ensure good cache behavior
		for (int block_index{}; block_index < out.blocks.size(); ++block_index)
		{
			auto out_block = out.blocks[block_index];
			const auto t_block = my_blocks[block_index];
			for (int k{}; k < out_block.size(); ++k)
			{
				out_block[k] = t_block[k] + d;
			}
		}

		return out;
	}

	/** {@inherit_doc} */
	//override
	Real_Matrix scalar_multiply(const double d)
	{
		const auto out = Block_Real_Matrix(my_rows, my_columns);

		// perform subtraction block-wise, to ensure good cache behavior
		for (int block_index{}; block_index < out.blocks.size(); ++block_index)
		{
			auto out_block = out.blocks[block_index];
			const auto t_block = my_blocks[block_index];
			for (int k{}; k < out_block.size(); ++k)
			{
				out_block[k] = t_block[k] * d;
			}
		}

		return out;
	}

	/** {@inherit_doc} */
	//override
	Block_Real_Matrix multiply(const Real_Matrix& m)

	{
		if (m instanceof Block_Real_Matrix)
		{
			return multiply((Block_Real_Matrix)m);
		}
		else
		{
			// safety check
			Matrix_Utils::check_multiplication_compatible(this, m);

			const auto out = Block_Real_Matrix(my_rows, m.get_column_dimension());

			// perform multiplication block-wise, to ensure good cache behavior
			int block_index{};
			for (int i_block{}; i_block < out.my_block_rows; ++i_block)
			{
				const int p_start{ i_block * BLOCK_SIZE };
				const int p_end{ std::min(p_start + BLOCK_SIZE, my_rows) };

				for (int j_block{}; j_block < out.my_block_columns; ++j_block)
				{
					const int q_start{ j_block * BLOCK_SIZE };
					const int q_end = std::min(q_start + BLOCK_SIZE, m.get_column_dimension());

					// select current block
					auto out_block = out.blocks[block_index];

					// perform multiplication on current block
					for (int k_block = 0; k_block < my_block_columns; ++k_block)
					{
						const int k_width = block_width(k_block);
						const auto t_block = my_blocks[i_block * my_block_columns + k_block];
						const int r_start = k_block * BLOCK_SIZE;
						int k{};
						for (int p{ p_start }; p < p_end; ++p)
						{
							const int l_start = (p - p_start) * k_width;
							const int l_end = l_start + k_width;
							for (int q{ q_start }; q < q_end; ++q)
							{
								double sum{};
								int r{ r_start };
								for (int l{ l_start }; l < l_end; ++l)
								{
									sum += t_block[l] * m.get_entry(r++, q);
								}
								out_block[k] += sum;
								++k;
							}
						}
					}
					// go to next block
					++block_index;
				}
			}

			return out;
		}
	}

	/**
	 * Returns the result of postmultiplying this by {@code m}.
	 *
	 * @param m Matrix to postmultiply by.
	 * @return {@code this} * m.
	 * @ if the matrices are not compatible.
	 */
	Block_Real_Matrix multiply(Block_Real_Matrix m)

	{
		// safety check
		Matrix_Utils::check_multiplication_compatible(this, m);

		const auto out = Block_Real_Matrix(my_rows, m.columns);

		// perform multiplication block-wise, to ensure good cache behavior
		int block_index{};
		for (int i_block{}; i_block < out.my_block_rows; ++i_block)
		{
			const int p_start{ i_block * BLOCK_SIZE };
			const int p_end{ std::min(p_start + BLOCK_SIZE, my_rows) };

			for (int j_block{}; j_block < out.my_block_columns; ++j_block)
			{
				const int j_width = out.block_width(j_block);
				const int j_width2 = j_width + j_width;
				const int j_width3 = j_width2 + j_width;
				const int j_width4 = j_width3 + j_width;

				// select current block
				auto out_block = out.blocks[block_index];

				// perform multiplication on current block
				for (int k_block = 0; k_block < my_block_columns; ++k_block)
				{
					const int k_width = block_width(k_block);
					const auto t_block = my_blocks[i_block * my_block_columns + k_block];
					const auto m_block = m.blocks[k_block * m.my_block_columns + j_block];
					int k{};
					for (int p{ p_start }; p < p_end; ++p)
					{
						const int l_start = (p - p_start) * k_width;
						const int l_end = l_start + k_width;
						for (int n_start = 0; n_start < j_width; ++n_start)
						{
							double sum{};
							int l = l_start;
							int n = n_start;
							while (l < l_end - 3)
							{
								sum += t_block[l] * m_block[n] +
									t_block[l + 1] * m_block[n + j_width] +
									t_block[l + 2] * m_block[n + j_width2] +
									t_block[l + 3] * m_block[n + j_width3];
								l += 4;
								n += j_width4;
							}
							while (l < l_end)
							{
								sum += t_block[l++] * m_block[n];
								n += j_width;
							}
							out_block[k] += sum;
							++k;
						}
					}
				}
				// go to next block
				++block_index;
			}
		}

		return out;
	}

	/**
	 * Returns the result of postmultiplying {@code this} by {@code m^T}.
	 * @param m matrix to first transpose and second postmultiply by
	 * @return {@code this * m^T}
	 * @ if
	 * {@code column_dimension(this) != column_dimension(m)}
	 * @since 1.3
	 */
	Block_Real_Matrix multiply_transposed(Block_Real_Matrix m)

	{
		// safety check
		Matrix_Utils::check_same_column_dimension(this, m);

		const auto out = Block_Real_Matrix(my_rows, m.rows);

		// perform multiplication block-wise, to ensure good cache behavior
		int block_index{};
		for (int i_block{}; i_block < out.my_block_rows; ++i_block)
		{
			const int p_start{ i_block * BLOCK_SIZE };
			const int p_end{ std::min(p_start + BLOCK_SIZE, my_rows) };

			for (int j_block{}; j_block < out.my_block_columns; ++j_block)
			{
				const int j_width = out.block_width(j_block);

				// select current block
				auto out_block = out.blocks[block_index];

				// perform multiplication on current block
				for (int k_block = 0; k_block < my_block_columns; ++k_block)
				{
					const int k_width = block_width(k_block);
					const auto t_block = my_blocks[i_block * my_block_columns + k_block];
					const std::vector<double> m_block = m.blocks[j_block * m.my_block_columns + k_block];
					int k{};
					for (int p{ p_start }; p < p_end; ++p)
					{
						const int l_start = (p - p_start) * k_width;
						const int l_end = l_start + k_width;
						for (int n_start = 0; n_start < j_width * k_width; n_start += k_width)
						{
							double sum{};
							int l = l_start;
							int n = n_start;
							while (l < l_end - 3)
							{
								sum += t_block[l] * m_block[n] +
									t_block[l + 1] * m_block[n + 1] +
									t_block[l + 2] * m_block[n + 2] +
									t_block[l + 3] * m_block[n + 3];
								l += 4;
								n += 4;
							}
							while (l < l_end)
							{
								sum += t_block[l++] * m_block[n++];
							}
							out_block[k] += sum;
							++k;
						}
					}
				}
				// go to next block
				++block_index;
			}
		}

		return out;
	}

	/** {@inherit_doc} */
	//override
	Block_Real_Matrix multiply_transposed(const Real_Matrix& m)

	{
		if (m instanceof Block_Real_Matrix)
		{
			return multiply_transposed((Block_Real_Matrix)m);
		}
		else
		{
			// safety check
			Matrix_Utils::check_same_column_dimension(this, m);

			const auto out = Block_Real_Matrix(my_rows, m.get_row_dimension());

			// perform multiplication block-wise, to ensure good cache behavior
			int block_index{};
			for (int i_block{}; i_block < out.my_block_rows; ++i_block)
			{
				const int p_start{ i_block * BLOCK_SIZE };
				const int p_end{ std::min(p_start + BLOCK_SIZE, my_rows) };

				for (int j_block{}; j_block < out.my_block_columns; ++j_block)
				{
					const int q_start{ j_block * BLOCK_SIZE };
					const int q_end = std::min(q_start + BLOCK_SIZE, m.get_row_dimension());

					// select current block
					auto out_block = out.blocks[block_index];

					// perform multiplication on current block
					for (int k_block = 0; k_block < my_block_columns; ++k_block)
					{
						const int k_width = block_width(k_block);
						const auto t_block = my_blocks[i_block * my_block_columns + k_block];
						const int r_start = k_block * BLOCK_SIZE;
						int k{};
						for (int p{ p_start }; p < p_end; ++p)
						{
							const int l_start{ (p - p_start) * k_width };
							const int l_end{ l_start + k_width };
							for (int q{ q_start }; q < q_end; ++q)
							{
								double sum{};
								int r{ r_start };
								for (int l{ l_start }; l < l_end; ++l)
								{
									sum += t_block[l] * m.get_entry(q, r++);
								}
								out_block[k] += sum;
								++k;
							}
						}
					}
					// go to next block
					++block_index;
				}
			}
			return out;
		}
	}

	/**
	 * Returns the result of postmultiplying {@code this^T} by {@code m}.
	 * @param m matrix to postmultiply by
	 * @return {@code this^T * m}
	 * @ if
	 * {@code column_dimension(this) != column_dimension(m)}
	 * @since 1.3
	 */
	Block_Real_Matrix transpose_multiply(const Block_Real_Matrix& m)

	{
		// safety check
		Matrix_Utils::check_same_row_dimension(this, m);

		const auto out = Block_Real_Matrix(my_columns, m.columns);

		// perform multiplication block-wise, to ensure good cache behavior
		int block_index{};
		for (int i_block{}; i_block < out.my_block_rows; ++i_block)
		{
			const int i_height = out.block_height(i_block);
			const int i_height2 = i_height + i_height;
			const int i_height3 = i_height2 + i_height;
			const int i_height4 = i_height3 + i_height;
			const int p_start{ i_block * BLOCK_SIZE };
			const int p_end = std::min(p_start + BLOCK_SIZE, columns);

			for (int j_block{}; j_block < out.my_block_columns; ++j_block)
			{
				const int j_width = out.block_width(j_block);
				const int j_width2 = j_width + j_width;
				const int j_width3 = j_width2 + j_width;
				const int j_width4 = j_width3 + j_width;

				// select current block
				auto out_block = out.blocks[block_index];

				// perform multiplication on current block
				for (int k_block = 0; k_block < my_block_rows; ++k_block)
				{
					const int k_height = block_height(k_block);
					const auto t_block = my_blocks[k_block * my_block_columns + i_block];
					const auto m_block = m.blocks[k_block * m.my_block_columns + j_block];
					int k{};
					for (int p{ p_start }; p < p_end; ++p)
					{
						const int l_start{ p - p_start };
						const int l_end{ l_start + i_height * k_height };
						for (int n_start{}; n_start < j_width; ++n_start)
						{
							double sum{};
							int l{ l_start };
							int n{ n_start };
							while (l < l_end - i_height3)
							{
								sum += t_block[l] * m_block[n] +
									t_block[l + i_height] * m_block[n + j_width] +
									t_block[l + i_height2] * m_block[n + j_width2] +
									t_block[l + i_height3] * m_block[n + j_width3];
								l += i_height4;
								n += j_width4;
							}
							while (l < l_end)
							{
								sum += t_block[l] * m_block[n];
								l += i_height;
								n += j_width;
							}
							out_block[k] += sum;
							++k;
						}
					}
				}
				// go to next block
				++block_index;
			}
		}

		return out;
	}

	/** {@inherit_doc} */
	//override
	Block_Real_Matrix transpose_multiply(const Real_Matrix& m)

	{
		if (m instanceof Block_Real_Matrix)
		{
			return transpose_multiply((Block_Real_Matrix)m);
		}
		else
		{
			// safety check
			Matrix_Utils::check_same_row_dimension(this, m);

			const auto out = Block_Real_Matrix(my_columns, m.get_column_dimension());

			// perform multiplication block-wise, to ensure good cache behavior
			int block_index{};
			for (int i_block{}; i_block < out.my_block_rows; ++i_block)
			{
				const int i_height = out.block_height(i_block);
				const int p_start{ i_block * BLOCK_SIZE };
				const int p_end = std::min(p_start + BLOCK_SIZE, my_columns);

				for (int j_block{}; j_block < out.my_block_columns; ++j_block)
				{
					const int q_start{ j_block * BLOCK_SIZE };
					const int q_end = std::min(q_start + BLOCK_SIZE, m.get_column_dimension());

					// select current block
					auto out_block = out.blocks[block_index];

					// perform multiplication on current block
					for (int k_block = 0; k_block < my_block_rows; ++k_block)
					{
						const int      k_height = block_height(k_block);
						const auto t_block{ my_blocks[k_block * my_block_columns + i_block] };
						const int      r_start{ k_block * BLOCK_SIZE };
						int k{};
						for (int p{ p_start }; p < p_end; ++p)
						{
							const int l_start{ p - p_start };
							const int l_end{ l_start + i_height * k_height };
							for (int q{ q_start }; q < q_end; ++q)
							{
								double sum{};
								int r{ r_start };
								for (int l{ l_start }; l < l_end; l += i_height)
								{
									sum += t_block[l] * m.get_entry(r++, q);
								}
								out_block[k] += sum;
								++k;
							}
						}
					}
					// go to next block
					++block_index;
				}
			}
			return out;
		}
	}

	/** {@inherit_doc} */
	//override
	std::vector<std::vector<double>> get_data()
	{
		const auto data = std::vector<std::vector<double>>(get_row_dimension(), std::vector<double>(get_column_dimension());
		const int last_columns = my_columns - (my_block_columns - 1) * BLOCK_SIZE;

		for (int i_block{}; i_block < my_block_rows; ++i_block)
		{
			const int p_start{ i_block * BLOCK_SIZE };
			const int p_end{ std::min(p_start + BLOCK_SIZE, my_rows) };
			int regular_pos{};
			int last_pos{};
			for (int p{ p_start }; p < p_end; ++p)
			{
				const auto data_p = data[p];
				int block_index{ i_block * my_block_columns };
				int data_pos{};
				for (int j_block{}; j_block < my_block_columns - 1; ++j_block)
				{
					System.arraycopy(blocks[block_index++], regular_pos, data_p, data_pos, BLOCK_SIZE);
					data_pos += BLOCK_SIZE;
				}
				System.arraycopy(blocks[block_index], last_pos, data_p, data_pos, last_columns);
				regular_pos += BLOCK_SIZE;
				last_pos += last_columns;
			}
		}

		return data;
	}

	/** {@inherit_doc} */
	//override
	double get_norm1()
	{
		const std::vector<double> col_sums = std::vector<double>(BLOCK_SIZE];
		double max_col_sum = 0;
		for (int j_block{}; j_block < my_block_columns; j_block++)
		{
			const int j_width = block_width(j_block);
			Arrays.fill(col_sums, 0, j_width, 0.0);
			for (int i_block{}; i_block < my_block_rows; ++i_block)
			{
				const int i_height = block_height(i_block);
				const auto block = my_blocks[i_block * my_block_columns + j_block];
				for (int j{}; j < j_width; ++j)
				{
					double sum{};
					for (int i{}; i < i_height; ++i)
					{
						sum += std::abs(block[i * j_width + j]);
					}
					col_sums[j] += sum;
				}
			}
			for (int j{}; j < j_width; ++j)
			{
				max_col_sum = std::max(max_col_sum, col_sums[j]);
			}
		}
		return max_col_sum;
	}

	/** {@inherit_doc} */
	//override
	double get_norm_infty()
	{
		const std::vector<double> row_sums = std::vector<double>(BLOCK_SIZE];
		double max_row_sum = 0;
		for (int i_block{}; i_block < my_block_rows; ++i_block)
		{
			const int i_height = block_height(i_block);
			Arrays.fill(row_sums, 0, i_height, 0.0);
			for (int j_block{}; j_block < my_block_columns; j_block++)
			{
				const int j_width = block_width(j_block);
				const auto block = my_blocks[i_block * my_block_columns + j_block];
				for (int i{}; i < i_height; ++i)
				{
					double sum{};
					for (int j{}; j < j_width; ++j)
					{
						sum += std::abs(block[i * j_width + j]);
					}
					row_sums[i] += sum;
				}
			}
			for (int i{}; i < i_height; ++i)
			{
				max_row_sum = std::max(max_row_sum, row_sums[i]);
			}
		}
		return max_row_sum;
	}

	/** {@inherit_doc} */
	//override
	double get_frobenius_norm()
	{
		double sum2 = 0;
		for (int block_index{}; block_index < blocks.size(); ++block_index)
		{
			for (const double entry : blocks[block_index])
			{
				sum2 += entry * entry;
			}
		}
		return std::sqrt(sum2);
	}

	/** {@inherit_doc} */
	//override
	Block_Real_Matrix get_sub_matrix(const int& start_row, const int& end_row, const int& start_column, const int& end_column)

	{
		// safety checks
		Matrix_Utils::check_sub_matrix_index(this, start_row, end_row, start_column, end_column);

		// create the output matrix
		const Block_Real_Matrix out =
			Block_Real_Matrix(end_row - start_row + 1, end_column - start_column + 1);

		// compute blocks shifts
		const int block_start_row = start_row / BLOCK_SIZE;
		const int rows_shift = start_row % BLOCK_SIZE;
		const int block_start_column = start_column / BLOCK_SIZE;
		const int columns_shift = start_column % BLOCK_SIZE;

		// perform extraction block-wise, to ensure good cache behavior
		int p_block = block_start_row;
		for (int i_block{}; i_block < out.my_block_rows; ++i_block)
		{
			const int i_height = out.block_height(i_block);
			int q_block = block_start_column;
			for (int j_block{}; j_block < out.my_block_columns; ++j_block)
			{
				const int j_width = out.block_width(j_block);

				// handle one block of the output matrix
				const int out_index = i_block * out.my_block_columns + j_block;
				auto out_block = out.blocks[out_index];
				const int index = p_block * my_block_columns + q_block;
				const int width = block_width(q_block);

				const int height_excess = i_height + rows_shift - BLOCK_SIZE;
				const int width_excess = j_width + columns_shift - BLOCK_SIZE;
				if (height_excess > 0)
				{
					// the submatrix block spans on two blocks rows from the original matrix
					if (width_excess > 0)
					{
						// the submatrix block spans on two blocks columns from the original matrix
						const int width2 = block_width(q_block + 1);
						copy_block_part(blocks[index], width, rows_shift, BLOCK_SIZE, columns_shift, BLOCK_SIZE, out_block, j_width, 0, 0);
						copy_block_part(blocks[index + 1], width2, rows_shift, BLOCK_SIZE, 0, width_excess, out_block, j_width, 0, j_width - width_excess);
						copy_block_part(blocks[index + my_block_columns], width, 0, height_excess, columns_shift, BLOCK_SIZE, out_block, j_width, i_height - height_excess, 0);
						copy_block_part(blocks[index + my_block_columns + 1], width2, 0, height_excess, 0, width_excess, out_block, j_width, i_height - height_excess, j_width - width_excess);
					}
					else
					{
						// the submatrix block spans on one block column from the original matrix
						copy_block_part(blocks[index], width, rows_shift, BLOCK_SIZE, columns_shift, j_width + columns_shift, out_block, j_width, 0, 0);
						copy_block_part(blocks[index + my_block_columns], width, 0, height_excess, columns_shift, j_width + columns_shift, out_block, j_width, i_height - height_excess, 0);
					}
				}
				else
				{
					// the submatrix block spans on one block row from the original matrix
					if (width_excess > 0)
					{
						// the submatrix block spans on two blocks columns from the original matrix
						const int width2 = block_width(q_block + 1);
						copy_block_part(blocks[index], width, rows_shift, i_height + rows_shift, columns_shift, BLOCK_SIZE, out_block, j_width, 0, 0);
						copy_block_part(blocks[index + 1], width2, rows_shift, i_height + rows_shift, 0, width_excess, out_block, j_width, 0, j_width - width_excess);
					}
					else
					{
						// the submatrix block spans on one block column from the original matrix
						copy_block_part(blocks[index], width, rows_shift, i_height + rows_shift, columns_shift, j_width + columns_shift, out_block, j_width, 0, 0);
					}
				}
				++q_block;
			}
			++p_block;
		}

		return out;
	}

	/** {@inherit_doc} */
	//override
	void set_sub_matrix(const std::vector<std::vector<double>> sub_matrix, const int& row, const int column)

	{
		// safety checks
		//Math_Utils::check_not_null(sub_matrix);
		const int ref_length = sub_matrix[0].size();
		if (ref_length == 0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::AT_LEAST_ONE_COLUMN);
		}
		const int end_row = row + sub_matrix.size() - 1;
		const int& end_column = column + ref_length - 1;
		Matrix_Utils::check_sub_matrix_index(this, row, end_row, column, end_column);
		for (const std::vector<double> sub_row : sub_matrix)
		{
			if (sub_row.size() != ref_length)
			{
				throw std::exception("not implemented");
				//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, ref_length, sub_row.size());
			}
		}

		// compute blocks bounds
		const int block_start_row = row / BLOCK_SIZE;
		const int block_end_row = (end_row + BLOCK_SIZE) / BLOCK_SIZE;
		const int block_start_column = column / BLOCK_SIZE;
		const int block_end_column = (end_column + BLOCK_SIZE) / BLOCK_SIZE;

		// perform copy block-wise, to ensure good cache behavior
		for (int i_block = block_start_row; i_block < block_end_row; ++i_block)
		{
			const int i_height = block_height(i_block);
			const int first_row{ i_block * BLOCK_SIZE };
			const int i_start = std::max(row, first_row);
			const int i_end = std::min(end_row + 1, first_row + i_height);

			for (int j_block = block_start_column; j_block < block_end_column; ++j_block)
			{
				const int j_width = block_width(j_block);
				const int first_column = j_block * BLOCK_SIZE;
				const int j_start = std::max(column, first_column);
				const int j_end = std::min(end_column + 1, first_column + j_width);
				const int j_length = j_end - j_start;

				// handle one block, row by row
				const auto block = my_blocks[i_block * my_block_columns + j_block];
				for (int i = i_start; i < i_end; ++i)
				{
					System.arraycopy(sub_matrix[i - row], j_start - column, block, (i - first_row) * j_width + (j_start - first_column), j_length);
				}
			}
		}
	}

	/** {@inherit_doc} */
	//override
	Block_Real_Matrix get_row_matrix(const int row)

	{
		Matrix_Utils::check_row_index(this, row);
		const auto out = Block_Real_Matrix(1, columns);

		// perform copy block-wise, to ensure good cache behavior
		const int i_block{ row / BLOCK_SIZE };
		const int i_row{ row - i_block * BLOCK_SIZE };
		int out_block_index = 0;
		int out_index{};
		std::vector<double> out_block = out.blocks[out_block_index];
		for (int j_block{}; j_block < my_block_columns; ++j_block)
		{
			const int j_width = block_width(j_block);
			const auto block = my_blocks[i_block * my_block_columns + j_block];
			const int available = out_block.size() - out_index;
			if (j_width > available)
			{
				System.arraycopy(block, i_row * j_width, out_block, out_index, available);
				out_block = out.blocks[++out_block_index];
				System.arraycopy(block, i_row * j_width, out_block, 0, j_width - available);
				out_index = j_width - available;
			}
			else
			{
				System.arraycopy(block, i_row * j_width, out_block, out_index, j_width);
				out_index += j_width;
			}
		}

		return out;
	}

	/** {@inherit_doc} */
	//override
	void set_row_matrix(const int& row, const Real_Matrix matrix)

	{
		if (matrix instanceof Block_Real_Matrix)
		{
			set_row_matrix(row, (Block_Real_Matrix)matrix);
		}
		else
		{
			super.set_row_matrix(row, matrix);
		}
	}

	/**
	 * Sets the entries in row number <code>row</code>
	 * as a row matrix.  Row indices start at 0.
	 *
	 * @param row the row to be set
	 * @param matrix row matrix (must have one row and the same number of columns
	 * as the instance)
	 * @ if the specified row index is invalid.
	 * @ if the matrix dimensions do
	 * not match one instance row.
	 */
	void set_row_matrix(const int& row, const Block_Real_Matrix& matrix)
	{
		Matrix_Utils::check_row_index(this, row);
		const int n_cols = get_column_dimension();
		if ((matrix.get_row_dimension() != 1) || (matrix.get_column_dimension() != n_cols))
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH_2x2, matrix.get_row_dimension(), matrix.get_column_dimension(), 1, n_cols);
		}

		// perform copy block-wise, to ensure good cache behavior
		const int i_block{ row / BLOCK_SIZE };
		const int i_row{ row - i_block * BLOCK_SIZE };
		int m_block_index{};
		int m_index{};
		std::vector<double> m_block = matrix.blocks[m_block_index];
		for (int j_block{}; j_block < my_block_columns; ++j_block)
		{
			const int j_width = block_width(j_block);
			const auto block = my_blocks[i_block * my_block_columns + j_block];
			const int available = m_block.size() - m_index;
			if (j_width > available)
			{
				System.arraycopy(m_block, m_index, block, i_row * j_width, available);
				m_block = matrix.blocks[++m_block_index];
				System.arraycopy(m_block, 0, block, i_row * j_width, j_width - available);
				m_index = j_width - available;
			}
			else
			{
				System.arraycopy(m_block, m_index, block, i_row * j_width, j_width);
				m_index += j_width;
			}
		}
	}

	/** {@inherit_doc} */
	//override
	Block_Real_Matrix get_column_matrix(const int& column)

	{
		Matrix_Utils::check_column_index(this, column);
		const auto out = Block_Real_Matrix(my_rows, 1);

		// perform copy block-wise, to ensure good cache behavior
		const int j_block = column / BLOCK_SIZE;
		const int j_column = column - j_block * BLOCK_SIZE;
		const int j_width = block_width(j_block);
		int out_block_index = 0;
		int out_index{};
		std::vector<double> out_block = out.blocks[out_block_index];
		for (int i_block{}; i_block < my_block_rows; ++i_block)
		{
			const int i_height = block_height(i_block);
			const auto block = my_blocks[i_block * my_block_columns + j_block];
			for (int i{}; i < i_height; ++i)
			{
				if (out_index >= out_block.size())
				{
					out_block = out.blocks[++out_block_index];
					out_index = 0;
				}
				out_block[out_index++] = block[i * j_width + j_column];
			}
		}

		return out;
	}

	/** {@inherit_doc} */
	//override
	void set_column_matrix(const int& column, const Real_Matrix matrix)

	{
		if (matrix instanceof Block_Real_Matrix)
		{
			set_column_matrix(column, (Block_Real_Matrix)matrix);
		}
		else
		{
			super.set_column_matrix(column, matrix);
		}
	}

	/**
	 * Sets the entries in column number <code>column</code>
	 * as a column matrix.  Column indices start at 0.
	 *
	 * @param column the column to be set
	 * @param matrix column matrix (must have one column and the same number of rows
	 * as the instance)
	 * @ if the specified column index is invalid.
	 * @ if the matrix dimensions do
	 * not match one instance column.
	 */
	void set_column_matrix(const int& column, const Block_Real_Matrix& matrix)
	{
		Matrix_Utils::check_column_index(this, column);
		const int n_rows = get_row_dimension();
		if ((matrix.get_row_dimension() != n_rows) || (matrix.get_column_dimension() != 1))
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH_2x2, matrix.get_row_dimension(), matrix.get_column_dimension(), n_rows, 1);
		}

		// perform copy block-wise, to ensure good cache behavior
		const int j_block{ column / BLOCK_SIZE };
		const int j_column{ column - j_block * BLOCK_SIZE };
		const int j_width{ block_width(j_block) };
		int m_block_index{};
		int m_index{};
		auto m_block = matrix.blocks[m_block_index];
		for (int i_block{}; i_block < my_block_rows; ++i_block)
		{
			const int i_height{ block_height(i_block) };
			const auto block = my_blocks[i_block * my_block_columns + j_block];
			for (int i{}; i < i_height; ++i)
			{
				if (m_index >= m_block.size())
				{
					m_block = matrix.blocks[++m_block_index];
					m_index = 0;
				}
				block[i * j_width + j_column] = m_block[m_index++];
			}
		}
	}

	/** {@inherit_doc} */
	//override
	Real_Vector get_row_vector(const int& row)
	{
		Matrix_Utils::check_row_index(this, row);
		const auto out_data = std::vector<double>(columns];

		// perform copy block-wise, to ensure good cache behavior
		const int i_block{ row / BLOCK_SIZE };
		const int i_row{ row - i_block * BLOCK_SIZE };
		int out_index{};
		for (int j_block{}; j_block < my_block_columns; ++j_block)
		{
			const int j_width = block_width(j_block);
			const auto block = my_blocks[i_block * my_block_columns + j_block];
			System.arraycopy(block, i_row * j_width, out_data, out_index, j_width);
			out_index += j_width;
		}

		return Array_Real_Vector(out_data, false);
	}

	/** {@inherit_doc} */
	//override
	void set_row_vector(const int& row, const Real_Vector vector)

	{
		if (vector instanceof Array_Real_Vector)
		{
			set_row(row, ((Array_Real_Vector)vector).get_data_ref());
		}
		else
		{
			super.set_row_vector(row, vector);
		}
	}

	/** {@inherit_doc} */
	//override
	Real_Vector get_column_vector(const int& column)

	{
		Matrix_Utils::check_column_index(this, column);
		const auto out_data = std::vector<double>(rows];

		// perform copy block-wise, to ensure good cache behavior
		const int j_block = column / BLOCK_SIZE;
		const int j_column = column - j_block * BLOCK_SIZE;
		const int j_width = block_width(j_block);
		int out_index{};
		for (int i_block{}; i_block < my_block_rows; ++i_block)
		{
			const int i_height = block_height(i_block);
			const auto block = my_blocks[i_block * my_block_columns + j_block];
			for (int i{}; i < i_height; ++i)
			{
				out_data[out_index++] = block[i * j_width + j_column];
			}
		}

		return Array_Real_Vector(out_data, false);
	}

	/** {@inherit_doc} */
	//override
	void set_column_vector(const int& column, const Real_Vector vector)

	{
		if (vector instanceof Array_Real_Vector)
		{
			set_column(column, ((Array_Real_Vector)vector).get_data_ref());
		}
		else
		{
			super.set_column_vector(column, vector);
		}
	}

	/** {@inherit_doc} */
	//override
	std::vector<double> get_row(const int row)
	{
		Matrix_Utils::check_row_index(this, row);
		const std::vector<double> out = std::vector<double>(columns];

		// perform copy block-wise, to ensure good cache behavior
		const int i_block{ row / BLOCK_SIZE };
		const int i_row{ row - i_block * BLOCK_SIZE };
		int out_index{};
		for (int j_block{}; j_block < my_block_columns; ++j_block)
		{
			const int j_width = block_width(j_block);
			const auto block = my_blocks[i_block * my_block_columns + j_block];
			System.arraycopy(block, i_row * j_width, out, out_index, j_width);
			out_index += j_width;
		}

		return out;
	}

	/** {@inherit_doc} */
	//override
	void set_row(const int& row, const std::vector<double> array)

	{
		Matrix_Utils::check_row_index(this, row);
		const int n_cols = get_column_dimension();
		if (array.size() != n_cols)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH_2x2, 1, array.size(), 1, n_cols);
		}

		// perform copy block-wise, to ensure good cache behavior
		const int i_block{ row / BLOCK_SIZE };
		const int i_row{ row - i_block * BLOCK_SIZE };
		int out_index{};
		for (int j_block{}; j_block < my_block_columns; ++j_block)
		{
			const int j_width = block_width(j_block);
			const auto block = my_blocks[i_block * my_block_columns + j_block];
			System.arraycopy(array, out_index, block, i_row * j_width, j_width);
			out_index += j_width;
		}
	}

	/** {@inherit_doc} */
	//override
	std::vector<double> get_column(const int& column)
	{
		Matrix_Utils::check_column_index(this, column);
		const std::vector<double> out = std::vector<double>(rows];

		// perform copy block-wise, to ensure good cache behavior
		const int j_block{ column / BLOCK_SIZE };
		const int j_column{ column - j_block * BLOCK_SIZE };
		const int j_width{ block_width(j_block) };
		int out_index{};
		for (int i_block{}; i_block < my_block_rows; ++i_block)
		{
			const int i_height = block_height(i_block);
			const auto block = my_blocks[i_block * my_block_columns + j_block];
			for (int i{}; i < i_height; ++i)
			{
				out[out_index++] = block[i * j_width + j_column];
			}
		}

		return out;
	}

	/** {@inherit_doc} */
	//override
	void set_column(const int& column, const std::vector<double>& arr)

	{
		Matrix_Utils::check_column_index(this, column);
		const int n_rows = get_row_dimension();
		if (arr.size() != n_rows)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH_2x2, arr.size(), 1, n_rows, 1);
		}

		// perform copy block-wise, to ensure good cache behavior
		const int j_block{ column / BLOCK_SIZE };
		const int j_column = column - j_block * BLOCK_SIZE;
		const int j_width = block_width(j_block);
		int out_index{};
		for (int i_block{}; i_block < my_block_rows; ++i_block)
		{
			const int i_height = block_height(i_block);
			const auto block = my_blocks[i_block * my_block_columns + j_block];
			for (int i{}; i < i_height; ++i)
			{
				block[i * j_width + j_column] = arr[out_index++];
			}
		}
	}

	/** {@inherit_doc} */
	//override
	double get_entry(const int& row, const int& column)
	{
		Matrix_Utils::check_matrix_index(this, row, column);
		const int i_block{ row / BLOCK_SIZE };
		const int j_block = column / BLOCK_SIZE;
		const int k{ (row - i_block * BLOCK_SIZE) * block_width(j_block) + (column - j_block * BLOCK_SIZE) };
		return blocks[i_block * my_block_columns + j_block][k];
	}

	/** {@inherit_doc} */
	//override
	void set_entry(const int& row, const int& column, const double& value)
	{
		Matrix_Utils::check_matrix_index(this, row, column);
		const int i_block{ row / BLOCK_SIZE };
		const int j_block = column / BLOCK_SIZE;
		const int k{ (row - i_block * BLOCK_SIZE) * block_width(j_block) + (column - j_block * BLOCK_SIZE) };
		blocks[i_block * my_block_columns + j_block][k] = value;
	}

	/** {@inherit_doc} */
	//override
	void add_to_entry(const int& row, const int& column, const double& increment)
	{
		Matrix_Utils::check_matrix_index(this, row, column);
		const int i_block{ row / BLOCK_SIZE };
		const int j_block{ column / BLOCK_SIZE };
		const int k{ (row - i_block * BLOCK_SIZE) * block_width(j_block) + (column - j_block * BLOCK_SIZE) };
		blocks[i_block * my_block_columns + j_block][k] += increment;
	}

	/** {@inherit_doc} */
	//override
	void multiply_entry(const int& row, const int& column, const double& factor)
	{
		Matrix_Utils::check_matrix_index(this, row, column);
		const int i_block{ row / BLOCK_SIZE };
		const int j_block = column / BLOCK_SIZE;
		const int k = (row - i_block * BLOCK_SIZE) * block_width(j_block) + (column - j_block * BLOCK_SIZE);
		blocks[i_block * my_block_columns + j_block][k] *= factor;
	}

	/** {@inherit_doc} */
	//override
	Block_Real_Matrix transpose()
	{
		const int n_rows = get_row_dimension();
		const int n_cols = get_column_dimension();
		const auto out = Block_Real_Matrix(n_cols, n_rows);

		// perform transpose block-wise, to ensure good cache behavior
		int block_index{};
		for (int i_block{}; i_block < my_block_columns; ++i_block)
		{
			const int p_start{ i_block * BLOCK_SIZE };
			const int p_end{ std::min(p_start + BLOCK_SIZE, columns) };
			for (int j_block{}; j_block < my_block_rows; ++j_block)
			{
				// transpose current block
				auto out_block = out.blocks[block_index];
				const auto t_block = my_blocks[j_block * my_block_columns + i_block];
				const int q_start{ j_block * BLOCK_SIZE };
				const int q_end = std::min(q_start + BLOCK_SIZE, my_rows);
				int k{};
				for (int p{ p_start }; p < p_end; ++p)
				{
					const int l_inc{ p_end - p_start };
					int l{ p - p_start };
					for (int q{ q_start }; q < q_end; ++q)
					{
						out_block[k] = t_block[l];
						++k;
						l += l_inc;
					}
				}
				// go to next block
				++block_index;
			}
		}

		return out;
	}

	/** {@inherit_doc} */
	//override
	int get_row_dimension() const
	{
		return my_rows;
	}

	/** {@inherit_doc} */
	//override
	int get_column_dimension() const
	{
		return my_columns;
	}

	/** {@inherit_doc} */
	//override
	std::vector<double> operate(const std::vector<double>& v)
	{
		if (v.size() != columns)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, v.size(), columns);
		}
		const auto out = std::vector<double>(rows];

		// perform multiplication block-wise, to ensure good cache behavior
		for (int i_block{}; i_block < my_block_rows; ++i_block)
		{
			const int p_start{ i_block * BLOCK_SIZE };
			const int p_end{ std::min(p_start + BLOCK_SIZE, my_rows) };
			for (int j_block{}; j_block < my_block_columns; ++j_block)
			{
				const auto block = my_blocks[i_block * my_block_columns + j_block];
				const int q_start{ j_block * BLOCK_SIZE };
				const int q_end{ std::min(q_start + BLOCK_SIZE, columns) };
				int k{};
				for (int p{ p_start }; p < p_end; ++p)
				{
					double sum{};
					int q = q_start;
					while (q < q_end - 3)
					{
						sum += block[k] * v[q] +
							block[k + 1] * v[q + 1] +
							block[k + 2] * v[q + 2] +
							block[k + 3] * v[q + 3];
						k += 4;
						q += 4;
					}
					while (q < q_end)
					{
						sum += block[k++] * v[q++];
					}
					out[p] += sum;
				}
			}
		}

		return out;
	}

	/** {@inherit_doc} */
	//override
	std::vector<double> pre_multiply(const std::vector<double>& v)
	{
		if (v.size() != rows)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, v.size(), rows);
		}
		const std::vector<double> out = std::vector<double>(columns];

		// perform multiplication block-wise, to ensure good cache behavior
		for (int j_block{}; j_block < my_block_columns; ++j_block)
		{
			const int j_width = block_width(j_block);
			const int j_width2 = j_width + j_width;
			const int j_width3 = j_width2 + j_width;
			const int j_width4 = j_width3 + j_width;
			const int q_start{ j_block * BLOCK_SIZE };
			const int q_end{ std::min(q_start + BLOCK_SIZE, columns) };
			for (int i_block{}; i_block < my_block_rows; ++i_block)
			{
				const auto block = my_blocks[i_block * my_block_columns + j_block];
				const int p_start{ i_block * BLOCK_SIZE };
				const int p_end{ std::min(p_start + BLOCK_SIZE, my_rows) };
				for (int q{ q_start }; q < q_end; ++q)
				{
					int k{ q - q_start };
					double sum{};
					int p{ p_start };
					while (p < p_end - 3)
					{
						sum += block[k] * v[p] +
							block[k + j_width] * v[p + 1] +
							block[k + j_width2] * v[p + 2] +
							block[k + j_width3] * v[p + 3];
						k += j_width4;
						p += 4;
					}
					while (p < p_end)
					{
						sum += block[k] * v[p++];
						k += j_width;
					}
					out[q] += sum;
				}
			}
		}

		return out;
	}

	/** {@inherit_doc} */
	//override
	double walk_in_row_order(const Real_Matrix_Changing_Visitor& visitor)
	{
		visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
		for (int i_block{}; i_block < my_block_rows; ++i_block)
		{
			const int p_start{ i_block * BLOCK_SIZE };
			const int p_end{ std::min(p_start + BLOCK_SIZE, my_rows) };
			for (int p{ p_start }; p < p_end; ++p)
			{
				for (int j_block{}; j_block < my_block_columns; ++j_block)
				{
					const int j_width{ block_width(j_block) };
					const int q_start{ j_block * BLOCK_SIZE };
					const int q_end{ std::min(q_start + BLOCK_SIZE, columns) };
					const auto block = my_blocks[i_block * my_block_columns + j_block];
					int k{ (p - p_start) * j_width };
					for (int q{ q_start }; q < q_end; ++q)
					{
						block[k] = visitor.visit(p, q, block[k]);
						++k;
					}
				}
			}
		}
		return visitor.end();
	}

	/** {@inherit_doc} */
	//override
	double walk_in_row_order(const Real_Matrix_Preserving_Visitor& visitor)
	{
		visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
		for (int i_block{}; i_block < my_block_rows; ++i_block)
		{
			const int p_start{ i_block * BLOCK_SIZE };
			const int p_end{ std::min(p_start + BLOCK_SIZE, my_rows) };
			for (int p{ p_start }; p < p_end; ++p)
			{
				for (int j_block{}; j_block < my_block_columns; ++j_block)
				{
					const int j_width = block_width(j_block);
					const int q_start{ j_block * BLOCK_SIZE };
					const int q_end{ std::min(q_start + BLOCK_SIZE, columns) };
					const auto block = my_blocks[i_block * my_block_columns + j_block];
					int k{ (p - p_start) * j_width };
					for (int q{ q_start }; q < q_end; ++q)
					{
						visitor.visit(p, q, block[k]);
						++k;
					}
				}
			}
		}
		return visitor.end();
	}

	/** {@inherit_doc} */
	//override
	double walk_in_row_order(const Real_Matrix_Changing_Visitor& visitor, const int& start_row, const int& end_row, const int& start_column, const int& end_column)
	{
		Matrix_Utils::check_sub_matrix_index(this, start_row, end_row, start_column, end_column);
		visitor.start(rows, columns, start_row, end_row, start_column, end_column);
		for (int i_block = start_row / BLOCK_SIZE; i_block < 1 + end_row / BLOCK_SIZE; ++i_block)
		{
			const int p0{ i_block * BLOCK_SIZE };
			const int p_start = std::max(start_row, p0);
			const int p_end = std::min((i_block + 1) * BLOCK_SIZE, 1 + end_row);
			for (int p{ p_start }; p < p_end; ++p)
			{
				for (int j_block = start_column / BLOCK_SIZE; j_block < 1 + end_column / BLOCK_SIZE; ++j_block)
				{
					const int j_width = block_width(j_block);
					const int q0{ j_block * BLOCK_SIZE };
					const int q_start{ std::max(start_column, q0) };
					const int q_end{ std::min((j_block + 1) * BLOCK_SIZE, 1 + end_column) };
					const auto block{ my_blocks[i_block * my_block_columns + j_block] };
					int k{ (p - p0) * j_width + q_start - q0 };
					for (int q{ q_start }; q < q_end; ++q)
					{
						block[k] = visitor.visit(p, q, block[k]);
						++k;
					}
				}
			}
		}
		return visitor.end();
	}

	/** {@inherit_doc} */
	//override
	double walk_in_row_order(const Real_Matrix_Preserving_Visitor& visitor, const int& start_row, const int& end_row, const int& start_column, const int& end_column)
	{
		Matrix_Utils::check_sub_matrix_index(this, start_row, end_row, start_column, end_column);
		visitor.start(rows, columns, start_row, end_row, start_column, end_column);
		for (int i_block{ start_row / BLOCK_SIZE }; i_block < 1 + end_row / BLOCK_SIZE; ++i_block)
		{
			const int p0{ i_block * BLOCK_SIZE };
			const int p_start = std::max(start_row, p0);
			const int p_end = std::min((i_block + 1) * BLOCK_SIZE, 1 + end_row);
			for (int p{ p_start }; p < p_end; ++p)
			{
				for (int j_block{ start_column / BLOCK_SIZE }; j_block < 1 + end_column / BLOCK_SIZE; ++j_block)
				{
					const int j_width = block_width(j_block);
					const int q0{ j_block * BLOCK_SIZE };
					const int q_start = std::max(start_column, q0);
					const int q_end = std::min((j_block + 1) * BLOCK_SIZE, 1 + end_column);
					const auto block = my_blocks[i_block * my_block_columns + j_block];
					int k = (p - p0) * j_width + q_start - q0;
					for (int q{ q_start }; q < q_end; ++q)
					{
						visitor.visit(p, q, block[k]);
						++k;
					}
				}
			}
		}
		return visitor.end();
	}

	/** {@inherit_doc} */
	//override
	double walk_in_optimized_order(const Real_Matrix_Changing_Visitor& visitor)
	{
		visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
		int block_index{};
		for (int i_block{}; i_block < my_block_rows; ++i_block)
		{
			const int p_start{ i_block * BLOCK_SIZE };
			const int p_end{ std::min(p_start + BLOCK_SIZE, my_rows) };
			for (int j_block{}; j_block < my_block_columns; ++j_block)
			{
				const int q_start{ j_block * BLOCK_SIZE };
				const int q_end{ std::min(q_start + BLOCK_SIZE, columns) };
				const auto block = my_blocks[block_index];
				int k{};
				for (int p{ p_start }; p < p_end; ++p)
				{
					for (int q{ q_start }; q < q_end; ++q)
					{
						block[k] = visitor.visit(p, q, block[k]);
						++k;
					}
				}
				++block_index;
			}
		}
		return visitor.end();
	}

	/** {@inherit_doc} */
	//override
	double walk_in_optimized_order(const Real_Matrix_Preserving_Visitor& visitor)
	{
		visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
		int block_index{};
		for (int i_block{}; i_block < my_block_rows; ++i_block)
		{
			const int p_start{ i_block * BLOCK_SIZE };
			const int p_end{ std::min(p_start + BLOCK_SIZE, my_rows) };
			for (int j_block{}; j_block < my_block_columns; ++j_block)
			{
				const int q_start{ j_block * BLOCK_SIZE };
				const int q_end{ std::min(q_start + BLOCK_SIZE, columns) };
				const auto block = my_blocks[block_index];
				int k{};
				for (int p{ p_start }; p < p_end; ++p)
				{
					for (int q{ q_start }; q < q_end; ++q)
					{
						visitor.visit(p, q, block[k]);
						++k;
					}
				}
				++block_index;
			}
		}
		return visitor.end();
	}

	/** {@inherit_doc} */
	//override
	double walk_in_optimized_order(const Real_Matrix_Changing_Visitor& visitor, const int& start_row, const int& end_row, const int& start_column, const int& end_column)
	{
		Matrix_Utils::check_sub_matrix_index(this, start_row, end_row, start_column, end_column);
		visitor.start(rows, columns, start_row, end_row, start_column, end_column);
		for (int i_block = start_row / BLOCK_SIZE; i_block < 1 + end_row / BLOCK_SIZE; ++i_block)
		{
			const int p0{ i_block * BLOCK_SIZE };
			const int p_start = std::max(start_row, p0);
			const int p_end = std::min((i_block + 1) * BLOCK_SIZE, 1 + end_row);
			for (int j_block = start_column / BLOCK_SIZE; j_block < 1 + end_column / BLOCK_SIZE; ++j_block)
			{
				const int j_width = block_width(j_block);
				const int q0{ j_block * BLOCK_SIZE };
				const int q_start = std::max(start_column, q0);
				const int q_end = std::min((j_block + 1) * BLOCK_SIZE, 1 + end_column);
				const auto block = my_blocks[i_block * my_block_columns + j_block];
				for (int p{ p_start }; p < p_end; ++p)
				{
					int k{ (p - p0) * j_width + q_start - q0 };
					for (int q{ q_start }; q < q_end; ++q)
					{
						block[k] = visitor.visit(p, q, block[k]);
						++k;
					}
				}
			}
		}
		return visitor.end();
	}

	/** {@inherit_doc} */
	//override
	double walk_in_optimized_order(const Real_Matrix_Preserving_Visitor& visitor, const int& start_row, const int& end_row, const int& start_column, const int& end_column)
	{
		Matrix_Utils::check_sub_matrix_index(this, start_row, end_row, start_column, end_column);
		visitor.start(rows, columns, start_row, end_row, start_column, end_column);
		for (int i_block = start_row / BLOCK_SIZE; i_block < 1 + end_row / BLOCK_SIZE; ++i_block)
		{
			const int p0{ i_block * BLOCK_SIZE };
			const int p_start = std::max(start_row, p0);
			const int p_end = std::min((i_block + 1) * BLOCK_SIZE, 1 + end_row);
			for (int j_block{ start_column / BLOCK_SIZE }; j_block < 1 + end_column / BLOCK_SIZE; ++j_block)
			{
				const int j_width = block_width(j_block);
				const int q0{ j_block * BLOCK_SIZE };
				const int q_start = std::max(start_column, q0);
				const int q_end = std::min((j_block + 1) * BLOCK_SIZE, 1 + end_column);
				const auto block = my_blocks[i_block * my_block_columns + j_block];
				for (int p{ p_start }; p < p_end; ++p)
				{
					int k{ (p - p0) * j_width + q_start - q0 };
					for (int q{ q_start }; q < q_end; ++q)
					{
						visitor.visit(p, q, block[k]);
						++k;
					}
				}
			}
		}
		return visitor.end();
	}
};