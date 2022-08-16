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

  //import org.hipparchus.Field;
  //import org.hipparchus.Field_Element;
  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.exception.;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Math_Arrays;
  //import org.hipparchus.util.Math_Utils;
#include <type_traits>
#include "MatrixUtils.h"
#include "../FieldElement.h"

/**
 * Cache-friendly implementation of Field_Matrix using a flat arrays to store
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
 * for matrix multiplication). The default value is to use 36x36 blocks.
 * </p>
 * <p>
 * The regular blocks represent {@link #BLOCK_SIZE} x {@link #BLOCK_SIZE} squares. Blocks
 * at right hand side and bottom side which may be smaller to fit matrix dimensions. The square
 * blocks are flattened in row major order in single dimension arrays which are therefore
 * {@link #BLOCK_SIZE}<sup>2</sup> elements long for regular blocks. The blocks are themselves
 * organized in row major order.
 * </p>
 * <p>
 * As an example, for a block size of 36x36, a 100x60 matrix would be stored in 6 blocks.
 * Block 0 would be a Field[1296] array holding the upper left 36x36 square, block 1 would be
 * a Field[1296] array holding the upper center 36x36 square, block 2 would be a Field[1008]
 * array holding the upper right 36x28 rectangle, block 3 would be a Field[864] array holding
 * the lower left 24x36 rectangle, block 4 would be a Field[864] array holding the lower center
 * 24x36 rectangle and block 5 would be a Field[672] array holding the lower right 24x28
 * rectangle.
 * </p>
 * <p>
 * The layout complexity overhead versus simple mapping of matrices to java
 * arrays is negligible for small matrices (about 1%). The gain from cache efficiency leads
 * to up to 3-fold improvements for matrices of moderate to large size.
 * </p>
 * @param <T> the type of the field elements
 */
template<typename T, typename std::enable_if<std::is_base_of<Field_Element, T>::value>::type* = nullptr>
class BlockField_Matrix : public Abstract_Field_Matrix<T>
{
	/** Block size. */
	public static const int BLOCK_SIZE = 36;

	-4602336630143123183L;
	/** Blocks of matrix entries. */
	private const T& blocks[][];
	/** Number of rows of the matrix. */
	private const int rows;
	/** Number of columns of the matrix. */
	private const int columns;
	/** Number of block rows of the matrix. */
	private const int& block_rows;
	/** Number of block columns of the matrix. */
	private const int& block_columns;

	/**
	 * Create a matrix with the supplied row and column dimensions.
	 *
	 * @param field Field to which the elements belong.
	 * @param rows Number of rows in the matrix.
	 * @param columns Number of columns in the matrix.
	 * @ if row or column dimension is not
	 * positive.
	 */
	public BlockField_Matrix(const Field<T>& field, const int& rows, const int& columns)
	{
		super(field, rows, columns);
		this.rows = rows;
		this.columns = columns;

		// number of blocks
		block_rows = (rows + BLOCK_SIZE - 1) / BLOCK_SIZE;
		block_columns = (columns + BLOCK_SIZE - 1) / BLOCK_SIZE;

		// allocate storage blocks, taking care of smaller ones at right and bottom
		blocks = create_blocks_layout(field, rows, columns);
	}

	/**
	 * Create a dense matrix copying entries from raw layout data.
	 * <p>The input array <em>must</em> already be in raw layout.</p>
	 * <p>Calling this constructor is equivalent to call:
	 * <pre>matrix = BlockField_Matrix<T>(get_field(), raw_data.size(), raw_data[0].size(), *                                   to_blocks_layout(raw_data), false);</pre>
	 * </p>
	 *
	 * @param raw_data Data for the matrix, in raw layout.
	 * @ if the {@code block_data} shape is
	 * inconsistent with block layout.
	 * @see #BlockField_Matrix(int, int, Field_Element[][], bool)
	 */
	public BlockField_Matrix(const std::vector<std::vector<T>>& raw_data)
	{
		this(raw_data.size(), raw_data[0].size(), to_blocks_layout(raw_data), false);
	}

	/**
	 * Create a dense matrix copying entries from block layout data.
	 * <p>The input array <em>must</em> already be in blocks layout.</p>
	 * @param rows  the number of rows in the matrix
	 * @param columns  the number of columns in the matrix
	 * @param block_data data for matrix
	 * @param copy_array if true, the input array will be copied, otherwise
	 * it will be referenced
	 *
	 * @ if the {@code block_data} shape is
	 * inconsistent with block layout.
	 * @ if row or column dimension is not
	 * positive.
	 * @see #create_blocks_layout(Field, int, int)
	 * @see #to_blocks_layout(Field_Element[][])
	 * @see #BlockField_Matrix(Field_Element[][])
	 */
	public BlockField_Matrix(const int rows, const int columns, const std::vector<std::vector<T>> block_data, const bool copy_array) // NOPMD - array copy is taken care of by parameter

	{
		super(extract_field(block_data), rows, columns);
		this.rows = rows;
		this.columns = columns;

		// number of blocks
		block_rows = (rows + BLOCK_SIZE - 1) / BLOCK_SIZE;
		block_columns = (columns + BLOCK_SIZE - 1) / BLOCK_SIZE;

		if (copy_array)
		{
			// allocate storage blocks, taking care of smaller ones at right and bottom
			blocks = Math_Arrays::build_array(get_field(), block_rows * block_columns, -1);
		}
		else
		{
			// reference existing array
			blocks = block_data;
		}

		int index = 0;
		for (int i_block{}; i_block < block_rows; ++i_block)
		{
			const int i_height = block_height(i_block);
			for (int j_block = 0; j_block < block_columns; ++j_block, ++index)
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
	 * is the layout used in {@link BlockField_Matrix} instances, where the matrix
	 * is split in square blocks (except at right and bottom side where blocks may
	 * be rectangular to fit matrix size) and each block is stored in a flattened
	 * one-dimensional array.
	 * </p>
	 * <p>
	 * This method creates an array in blocks layout from an input array in raw layout.
	 * It can be used to provide the array argument of the {@link
	 * #BlockField_Matrix(int, int, Field_Element[][], bool)}
	 * constructor.
	 * </p>
	 * @param <T> Type of the field elements.
	 * @param raw_data Data array in raw layout.
	 * @return a data array containing the same entries but in blocks layout
	 * @ if {@code raw_data} is not rectangular
	 *  (not all rows have the same length).
	 * @see #create_blocks_layout(Field, int, int)
	 * @see #BlockField_Matrix(int, int, Field_Element[][], bool)
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Field_Element<T>, T>::value>::type* = nullptr>
	public static std::vector<std::vector<T>> to_blocks_layout(const std::vector<std::vector<T>> raw_data)
	{
		const int rows = raw_data.size();
		const int columns = raw_data[0].size();
		const int& block_rows = (rows + BLOCK_SIZE - 1) / BLOCK_SIZE;
		const int& block_columns = (columns + BLOCK_SIZE - 1) / BLOCK_SIZE;

		// safety checks
		for (int i{}; i < raw_data.size(); ++i)
		{
			const int length = raw_data[i].size();
			if (length != columns)
			{
				throw std::exception("not implemented");
				//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, columns, length);
			}
		}

		// convert array
		const Field<T> field = extract_field(raw_data);
		const std::vector<std::vector<T>> blocks = Math_Arrays::build_array(field, block_rows * block_columns, -1);
		int block_index = 0;
		for (const int& i_block = 0; i_block < block_rows; ++i_block)
		{
			const int p_start = i_block * BLOCK_SIZE;
			const int p_end = std::min(p_start + BLOCK_SIZE, rows);
			const int i_height = p_end - p_start;
			for (int j_block = 0; j_block < block_columns; ++j_block)
			{
				const int q_start = j_block * BLOCK_SIZE;
				const int q_end = std::min(q_start + BLOCK_SIZE, columns);
				const int j_width = q_end - q_start;

				// allocate block
				const std::vector<T> block = Math_Arrays::build_array(field, i_height * j_width);
				blocks[block_index] = block;

				// copy data
				int index = 0;
				for (const int& p = p_start; p < p_end; ++p)
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
	 * #BlockField_Matrix(int, int, Field_Element[][], bool)}
	 * constructor.
	 * </p>
	 * @param <T> Type of the field elements.
	 * @param field Field to which the elements belong.
	 * @param rows Number of rows in the matrix.
	 * @param columns Number of columns in the matrix.
	 * @return a data array in blocks layout.
	 * @see #to_blocks_layout(Field_Element[][])
	 * @see #BlockField_Matrix(int, int, Field_Element[][], bool)
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Field_Element<T>, T>::value>::type* = nullptr>
	public static std::vector<std::vector<T>> create_blocks_layout(const Field<T>& field, const int& rows, const int& columns)
	{
		const int& block_rows = (rows + BLOCK_SIZE - 1) / BLOCK_SIZE;
		const int& block_columns = (columns + BLOCK_SIZE - 1) / BLOCK_SIZE;

		const std::vector<std::vector<T>> blocks = Math_Arrays::build_array(field, block_rows * block_columns, -1);
		int block_index = 0;
		for (const int& i_block = 0; i_block < block_rows; ++i_block)
		{
			const int p_start = i_block * BLOCK_SIZE;
			const int p_end = std::min(p_start + BLOCK_SIZE, rows);
			const int i_height = p_end - p_start;
			for (int j_block = 0; j_block < block_columns; ++j_block)
			{
				const int q_start = j_block * BLOCK_SIZE;
				const int q_end = std::min(q_start + BLOCK_SIZE, columns);
				const int j_width = q_end - q_start;
				blocks[block_index] = Math_Arrays::build_array(field, i_height * j_width);
				++block_index;
			}
		}

		return blocks;
	}

	/** {@inherit_doc} */
	//override
	public Field_Matrix<T> create_matrix(const int row_dimension, const int column_dimension)

	{
		return BlockField_Matrix<T>(get_field(), row_dimension, column_dimension);
	}

	/** {@inherit_doc} */
	//override
	public Field_Matrix<T> copy()
	{
		// create an empty matrix
		BlockField_Matrix<T> copied = BlockField_Matrix<>(get_field(), rows, columns);

		// copy the blocks
		for (int i{}; i < blocks.size(); ++i)
		{
			System.arraycopy(blocks[i], 0, copied.blocks[i], 0, blocks[i].size());
		}

		return copied;
	}

	/** {@inherit_doc} */
	//override
	public Field_Matrix<T> add(const Field_Matrix<T> m)

	{
		if (m instanceof BlockField_Matrix)
		{
			return add((BlockField_Matrix<T>) m);
		}
		else
		{
			// safety check
			check_addition_compatible(m);

			const BlockField_Matrix<T> out = BlockField_Matrix<>(get_field(), rows, columns);

			// perform addition block-wise, to ensure good cache behavior
			int block_index = 0;
			for (const int& i_block = 0; i_block < out.block_rows; ++i_block)
			{
				for (int j_block = 0; j_block < out.block_columns; ++j_block)
				{
					// perform addition on the current block
					const std::vector<T> out_block = out.blocks[block_index];
					const std::vector<T> t_block = blocks[block_index];
					const int      p_start = i_block * BLOCK_SIZE;
					const int      p_end = std::min(p_start + BLOCK_SIZE, rows);
					const int      q_start = j_block * BLOCK_SIZE;
					const int      q_end = std::min(q_start + BLOCK_SIZE, columns);
					int k = 0;
					for (const int& p = p_start; p < p_end; ++p)
					{
						for (const int& q = q_start; q < q_end; ++q)
						{
							out_block[k] = t_block[k].add(m.get_entry(p, q));
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
	 * Compute the sum of {@code this} and {@code m}.
	 *
	 * @param m matrix to be added
	 * @return {@code this + m}
	 * @ if {@code m} is not the same
	 * size as {@code this}
	 */
	public BlockField_Matrix<T> add(const BlockField_Matrix<T> m)

	{
		// safety check
		check_addition_compatible(m);

		const BlockField_Matrix<T> out = BlockField_Matrix<>(get_field(), rows, columns);

		// perform addition block-wise, to ensure good cache behavior
		for (const int& block_index = 0; block_index < out.blocks.size(); ++block_index)
		{
			const std::vector<T> out_block = out.blocks[block_index];
			const std::vector<T> t_block = blocks[block_index];
			const std::vector<T> m_block = m.blocks[block_index];
			for (int k{}; k < out_block.size(); ++k)
			{
				out_block[k] = t_block[k].add(m_block[k]);
			}
		}

		return out;
	}

	/** {@inherit_doc} */
	//override
	public Field_Matrix<T> subtract(const Field_Matrix<T> m)

	{
		if (m instanceof BlockField_Matrix)
		{
			return subtract((BlockField_Matrix<T>) m);
		}
		else
		{
			// safety check
			check_subtraction_compatible(m);

			const BlockField_Matrix<T> out = BlockField_Matrix<>(get_field(), rows, columns);

			// perform subtraction block-wise, to ensure good cache behavior
			int block_index = 0;
			for (const int& i_block = 0; i_block < out.block_rows; ++i_block)
			{
				for (int j_block = 0; j_block < out.block_columns; ++j_block)
				{
					// perform subtraction on the current block
					const std::vector<T> out_block = out.blocks[block_index];
					const std::vector<T> t_block = blocks[block_index];
					const int      p_start = i_block * BLOCK_SIZE;
					const int      p_end = std::min(p_start + BLOCK_SIZE, rows);
					const int      q_start = j_block * BLOCK_SIZE;
					const int      q_end = std::min(q_start + BLOCK_SIZE, columns);
					int k = 0;
					for (const int& p = p_start; p < p_end; ++p)
					{
						for (const int& q = q_start; q < q_end; ++q)
						{
							out_block[k] = t_block[k].subtract(m.get_entry(p, q));
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
	 * Compute {@code this - m}.
	 *
	 * @param m matrix to be subtracted
	 * @return {@code this - m}
	 * @ if {@code m} is not the same
	 * size as {@code this}
	 */
	public BlockField_Matrix<T> subtract(const BlockField_Matrix<T> m)
	{
		// safety check
		check_subtraction_compatible(m);

		const BlockField_Matrix<T> out = BlockField_Matrix<>(get_field(), rows, columns);

		// perform subtraction block-wise, to ensure good cache behavior
		for (const int& block_index = 0; block_index < out.blocks.size(); ++block_index)
		{
			const std::vector<T> out_block = out.blocks[block_index];
			const std::vector<T> t_block = blocks[block_index];
			const std::vector<T> m_block = m.blocks[block_index];
			for (int k{}; k < out_block.size(); ++k)
			{
				out_block[k] = t_block[k].subtract(m_block[k]);
			}
		}

		return out;
	}

	/** {@inherit_doc} */
	//override
	public Field_Matrix<T> scalar_add(const T d)
	{
		const BlockField_Matrix<T> out = BlockField_Matrix<>(get_field(), rows, columns);

		// perform subtraction block-wise, to ensure good cache behavior
		for (const int& block_index = 0; block_index < out.blocks.size(); ++block_index)
		{
			const std::vector<T> out_block = out.blocks[block_index];
			const std::vector<T> t_block = blocks[block_index];
			for (int k{}; k < out_block.size(); ++k)
			{
				out_block[k] = t_block[k].add(d);
			}
		}

		return out;
	}

	/** {@inherit_doc} */
	//override
	public Field_Matrix<T> scalar_multiply(const T d)
	{
		const BlockField_Matrix<T> out = BlockField_Matrix<>(get_field(), rows, columns);

		// perform subtraction block-wise, to ensure good cache behavior
		for (const int& block_index = 0; block_index < out.blocks.size(); ++block_index)
		{
			const std::vector<T> out_block = out.blocks[block_index];
			const std::vector<T> t_block = blocks[block_index];
			for (int k{}; k < out_block.size(); ++k)
			{
				out_block[k] = t_block[k].multiply(d);
			}
		}

		return out;
	}

	/** {@inherit_doc} */
	//override
	public Field_Matrix<T> multiply(const Field_Matrix<T> m)

	{
		if (m instanceof BlockField_Matrix)
		{
			return multiply((BlockField_Matrix<T>) m);
		}
		else
		{
			// safety check
			check_multiplication_compatible(m);

			const BlockField_Matrix<T> out = BlockField_Matrix<>(get_field(), rows, m.get_column_dimension());
			const T zero = get_field().get_zero();

			// perform multiplication block-wise, to ensure good cache behavior
			int block_index = 0;
			for (const int& i_block = 0; i_block < out.block_rows; ++i_block)
			{
				const int p_start = i_block * BLOCK_SIZE;
				const int p_end = std::min(p_start + BLOCK_SIZE, rows);

				for (int j_block = 0; j_block < out.block_columns; ++j_block)
				{
					const int q_start = j_block * BLOCK_SIZE;
					const int q_end = std::min(q_start + BLOCK_SIZE, m.get_column_dimension());

					// select current block
					const std::vector<T> out_block = out.blocks[block_index];

					// perform multiplication on current block
					for (const int& k_block = 0; k_block < block_columns; ++k_block)
					{
						const int& k_width = block_width(k_block);
						const std::vector<T> t_block = blocks[i_block * block_columns + k_block];
						const int r_start = k_block * BLOCK_SIZE;
						int k = 0;
						for (const int& p = p_start; p < p_end; ++p)
						{
							const int l_start = (p - p_start) * k_width;
							const int l_end = l_start + k_width;
							for (const int& q = q_start; q < q_end; ++q)
							{
								T sum = zero;
								int r = r_start;
								for (const int& l = l_start; l < l_end; ++l)
								{
									sum = sum.add(t_block[l].multiply(m.get_entry(r, q)));
									++r;
								}
								out_block[k] = out_block[k].add(sum);
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
	 * Returns the result of postmultiplying {@code this} by {@code m}.
	 *
	 * @param m matrix to postmultiply by
	 * @return {@code this * m}
	 * @ if the matrices are not compatible.
	 */
	public BlockField_Matrix<T> multiply(BlockField_Matrix<T> m)

	{
		// safety check
		check_multiplication_compatible(m);

		const BlockField_Matrix<T> out = BlockField_Matrix<>(get_field(), rows, m.columns);
		const T zero = get_field().get_zero();

		// perform multiplication block-wise, to ensure good cache behavior
		int block_index = 0;
		for (const int& i_block = 0; i_block < out.block_rows; ++i_block)
		{
			const int p_start = i_block * BLOCK_SIZE;
			const int p_end = std::min(p_start + BLOCK_SIZE, rows);

			for (int j_block = 0; j_block < out.block_columns; ++j_block)
			{
				const int j_width = out.block_width(j_block);
				const int j_width2 = j_width + j_width;
				const int j_width3 = j_width2 + j_width;
				const int j_width4 = j_width3 + j_width;

				// select current block
				const std::vector<T> out_block = out.blocks[block_index];

				// perform multiplication on current block
				for (const int& k_block = 0; k_block < block_columns; ++k_block)
				{
					const int& k_width = block_width(k_block);
					const std::vector<T> t_block = blocks[i_block * block_columns + k_block];
					const std::vector<T> m_block = m.blocks[k_block * m.block_columns + j_block];
					int k = 0;
					for (const int& p = p_start; p < p_end; ++p)
					{
						const int l_start = (p - p_start) * k_width;
						const int l_end = l_start + k_width;
						for (const int& n_start = 0; n_start < j_width; ++n_start)
						{
							T sum = zero;
							int l = l_start;
							int n = n_start;
							while (l < l_end - 3)
							{
								sum = sum.
									add(t_block[l].multiply(m_block[n])).
									add(t_block[l + 1].multiply(m_block[n + j_width])).
									add(t_block[l + 2].multiply(m_block[n + j_width2])).
									add(t_block[l + 3].multiply(m_block[n + j_width3]));
								l += 4;
								n += j_width4;
							}
							while (l < l_end)
							{
								sum = sum.add(t_block[l++].multiply(m_block[n]));
								n += j_width;
							}
							out_block[k] = out_block[k].add(sum);
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
	 * @return {@code this * m}
	 * @ if
	 * {@code column_dimension(this) != column_dimension(m)}
	 * @since 1.3
	 */
	public BlockField_Matrix<T> multiply_transposed(BlockField_Matrix<T> m)

	{
		// safety check
		Matrix_Utils::check_same_column_dimension(this, m);

		const BlockField_Matrix<T> out = BlockField_Matrix<>(get_field(), rows, m.rows);

		// perform multiplication block-wise, to ensure good cache behavior
		int block_index = 0;
		for (const int& i_block = 0; i_block < out.block_rows; ++i_block)
		{
			const int p_start = i_block * BLOCK_SIZE;
			const int p_end = std::min(p_start + BLOCK_SIZE, rows);

			for (int j_block = 0; j_block < out.block_columns; ++j_block)
			{
				const int j_width = out.block_width(j_block);

				// select current block
				const std::vector<T> out_block = out.blocks[block_index];

				// perform multiplication on current block
				for (const int& k_block = 0; k_block < block_columns; ++k_block)
				{
					const int& k_width = block_width(k_block);
					const std::vector<T> t_block = blocks[i_block * block_columns + k_block];
					const std::vector<T> m_block = m.blocks[j_block * m.block_columns + k_block];
					int k = 0;
					for (const int& p = p_start; p < p_end; ++p)
					{
						const int l_start = (p - p_start) * k_width;
						const int l_end = l_start + k_width;
						for (const int& n_start = 0; n_start < j_width * k_width; n_start += k_width)
						{
							T sum = get_field().get_zero();
							int l = l_start;
							int n = n_start;
							while (l < l_end - 3)
							{
								sum = sum.
									add(t_block[l].multiply(m_block[n])).
									add(t_block[l + 1].multiply(m_block[n + 1])).
									add(t_block[l + 2].multiply(m_block[n + 2])).
									add(t_block[l + 3].multiply(m_block[n + 3]));
								l += 4;
								n += 4;
							}
							while (l < l_end)
							{
								sum = sum.add(t_block[l++].multiply(m_block[n++]));
							}
							out_block[k] = out_block[k].add(sum);
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
	public BlockField_Matrix<T> multiply_transposed(const Field_Matrix<T> m)

	{
		if (m instanceof BlockField_Matrix)
		{
			return multiply_transposed((BlockField_Matrix<T>) m);
		}
		else
		{
			// safety check
			Matrix_Utils::check_same_column_dimension(this, m);

			const BlockField_Matrix<T> out = BlockField_Matrix<>(get_field(), rows, m.get_row_dimension());

			// perform multiplication block-wise, to ensure good cache behavior
			int block_index = 0;
			for (const int& i_block = 0; i_block < out.block_rows; ++i_block)
			{
				const int p_start = i_block * BLOCK_SIZE;
				const int p_end = std::min(p_start + BLOCK_SIZE, rows);

				for (int j_block = 0; j_block < out.block_columns; ++j_block)
				{
					const int q_start = j_block * BLOCK_SIZE;
					const int q_end = std::min(q_start + BLOCK_SIZE, m.get_row_dimension());

					// select current block
					const std::vector<T> out_block = out.blocks[block_index];

					// perform multiplication on current block
					for (const int& k_block = 0; k_block < block_columns; ++k_block)
					{
						const int& k_width = block_width(k_block);
						const std::vector<T> t_block = blocks[i_block * block_columns + k_block];
						const int r_start = k_block * BLOCK_SIZE;
						int k = 0;
						for (const int& p = p_start; p < p_end; ++p)
						{
							const int l_start = (p - p_start) * k_width;
							const int l_end = l_start + k_width;
							for (const int& q = q_start; q < q_end; ++q)
							{
								T sum = get_field().get_zero();
								int r = r_start;
								for (const int& l = l_start; l < l_end; ++l)
								{
									sum = sum.add(t_block[l].multiply(m.get_entry(q, r)));
									++r;
								}
								out_block[k] = out_block[k].add(sum);
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
	public BlockField_Matrix<T> transpose_multiply(const BlockField_Matrix<T> m)

	{
		// safety check
		Matrix_Utils::check_same_row_dimension(this, m);

		const BlockField_Matrix<T> out = BlockField_Matrix<>(get_field(), columns, m.columns);

		// perform multiplication block-wise, to ensure good cache behavior
		int block_index = 0;
		for (const int& i_block = 0; i_block < out.block_rows; ++i_block)
		{
			const int i_height = out.block_height(i_block);
			const int i_height2 = i_height + i_height;
			const int i_height3 = i_height2 + i_height;
			const int i_height4 = i_height3 + i_height;
			const int p_start = i_block * BLOCK_SIZE;
			const int p_end = std::min(p_start + BLOCK_SIZE, columns);

			for (int j_block = 0; j_block < out.block_columns; ++j_block)
			{
				const int j_width = out.block_width(j_block);
				const int j_width2 = j_width + j_width;
				const int j_width3 = j_width2 + j_width;
				const int j_width4 = j_width3 + j_width;

				// select current block
				const std::vector<T> out_block = out.blocks[block_index];

				// perform multiplication on current block
				for (const int& k_block = 0; k_block < block_rows; ++k_block)
				{
					const int& k_height = block_height(k_block);
					const std::vector<T> t_block = blocks[k_block * block_columns + i_block];
					const std::vector<T> m_block = m.blocks[k_block * m.block_columns + j_block];
					int k = 0;
					for (const int& p = p_start; p < p_end; ++p)
					{
						const int l_start = p - p_start;
						const int l_end = l_start + i_height * k_height;
						for (const int& n_start = 0; n_start < j_width; ++n_start)
						{
							T sum = get_field().get_zero();
							int l = l_start;
							int n = n_start;
							while (l < l_end - i_height3)
							{
								sum = sum.add(t_block[l].multiply(m_block[n]).
									add(t_block[l + i_height].multiply(m_block[n + j_width])).
									add(t_block[l + i_height2].multiply(m_block[n + j_width2])).
									add(t_block[l + i_height3].multiply(m_block[n + j_width3])));
								l += i_height4;
								n += j_width4;
							}
							while (l < l_end)
							{
								sum = sum.add(t_block[l].multiply(m_block[n]));
								l += i_height;
								n += j_width;
							}
							out_block[k] = out_block[k].add(sum);
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
	public BlockField_Matrix<T> transpose_multiply(const Field_Matrix<T> m)

	{
		if (m instanceof BlockField_Matrix)
		{
			return transpose_multiply((BlockField_Matrix<T>) m);
		}
		else
		{
			// safety check
			Matrix_Utils::check_same_row_dimension(this, m);

			const BlockField_Matrix<T> out = BlockField_Matrix<>(get_field(), columns, m.get_column_dimension());

			// perform multiplication block-wise, to ensure good cache behavior
			int block_index = 0;
			for (const int& i_block = 0; i_block < out.block_rows; ++i_block)
			{
				const int i_height = out.block_height(i_block);
				const int p_start = i_block * BLOCK_SIZE;
				const int p_end = std::min(p_start + BLOCK_SIZE, columns);

				for (int j_block = 0; j_block < out.block_columns; ++j_block)
				{
					const int q_start = j_block * BLOCK_SIZE;
					const int q_end = std::min(q_start + BLOCK_SIZE, m.get_column_dimension());

					// select current block
					const std::vector<T> out_block = out.blocks[block_index];

					// perform multiplication on current block
					for (const int& k_block = 0; k_block < block_rows; ++k_block)
					{
						const int& k_height = block_height(k_block);
						const std::vector<T> t_block = blocks[k_block * block_columns + i_block];
						const int r_start = k_block * BLOCK_SIZE;
						int k = 0;
						for (const int& p = p_start; p < p_end; ++p)
						{
							const int l_start = p - p_start;
							const int l_end = l_start + i_height * k_height;
							for (const int& q = q_start; q < q_end; ++q)
							{
								T sum = get_field().get_zero();
								int r = r_start;
								for (const int& l = l_start; l < l_end; l += i_height)
								{
									sum = sum.add(t_block[l].multiply(m.get_entry(r++, q)));
								}
								out_block[k] = out_block[k].add(sum);
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
	public std::vector<std::vector<T>> get_data()
	{
		const std::vector<std::vector<T>> data = Math_Arrays::build_array(get_field(), get_row_dimension(), get_column_dimension());
		const int last_columns = columns - (block_columns - 1) * BLOCK_SIZE;

		for (const int& i_block = 0; i_block < block_rows; ++i_block)
		{
			const int p_start = i_block * BLOCK_SIZE;
			const int p_end = std::min(p_start + BLOCK_SIZE, rows);
			int regular_pos = 0;
			int last_pos = 0;
			for (const int& p = p_start; p < p_end; ++p)
			{
				const std::vector<T> data_p = data[p];
				int block_index = i_block * block_columns;
				int data_pos = 0;
				for (int j_block = 0; j_block < block_columns - 1; ++j_block)
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
	public Field_Matrix<T> get_sub_matrix(const int& start_row, const int& end_row, const int& start_column, const int& end_column)

	{
		// safety checks
		check_sub_matrix_index(start_row, end_row, start_column, end_column);

		// create the output matrix
		const BlockField_Matrix<T> out =
			BlockField_Matrix<>(get_field(), end_row - start_row + 1, end_column - start_column + 1);

		// compute blocks shifts
		const int block_start_row = start_row / BLOCK_SIZE;
		const int rows_shift = start_row % BLOCK_SIZE;
		const int block_start_column = start_column / BLOCK_SIZE;
		const int columns_shift = start_column % BLOCK_SIZE;

		// perform extraction block-wise, to ensure good cache behavior
		int p_block = block_start_row;
		for (const int& i_block = 0; i_block < out.block_rows; ++i_block)
		{
			const int i_height = out.block_height(i_block);
			int q_block = block_start_column;
			for (int j_block = 0; j_block < out.block_columns; ++j_block)
			{
				const int j_width = out.block_width(j_block);

				// handle one block of the output matrix
				const int      out_index = i_block * out.block_columns + j_block;
				const std::vector<T> out_block = out.blocks[out_index];
				const int      index = p_block * block_columns + q_block;
				const int      width = block_width(q_block);

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
						copy_block_part(blocks[index + block_columns], width, 0, height_excess, columns_shift, BLOCK_SIZE, out_block, j_width, i_height - height_excess, 0);
						copy_block_part(blocks[index + block_columns + 1], width2, 0, height_excess, 0, width_excess, out_block, j_width, i_height - height_excess, j_width - width_excess);
					}
					else
					{
						// the submatrix block spans on one block column from the original matrix
						copy_block_part(blocks[index], width, rows_shift, BLOCK_SIZE, columns_shift, j_width + columns_shift, out_block, j_width, 0, 0);
						copy_block_part(blocks[index + block_columns], width, 0, height_excess, columns_shift, j_width + columns_shift, out_block, j_width, i_height - height_excess, 0);
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
	private void copy_block_part(const std::vector<T> src_block, const int src_width, const int src_start_row, const int src_end_row, const int src_start_column, const int src_end_column, const std::vector<T> dst_block, const int dst_width, const int dst_start_row, const int dst_start_column)
	{
		const int length = src_end_column - src_start_column;
		int src_pos = src_start_row * src_width + src_start_column;
		int dst_pos = dst_start_row * dst_width + dst_start_column;
		for (const int& src_row = src_start_row; src_row < src_end_row; ++src_row)
		{
			System.arraycopy(src_block, src_pos, dst_block, dst_pos, length);
			src_pos += src_width;
			dst_pos += dst_width;
		}
	}

	/** {@inherit_doc} */
	//override
	public void set_sub_matrix(const std::vector<std::vector<T>> sub_matrix, const int& row, const int column)

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
		check_sub_matrix_index(row, end_row, column, end_column);
		for (const std::vector<T> sub_row : sub_matrix)
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
		for (const int& i_block = block_start_row; i_block < block_end_row; ++i_block)
		{
			const int i_height = block_height(i_block);
			const int first_row = i_block * BLOCK_SIZE;
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
				const std::vector<T> block = blocks[i_block * block_columns + j_block];
				for (int i = i_start; i < i_end; ++i)
				{
					System.arraycopy(sub_matrix[i - row], j_start - column, block, (i - first_row) * j_width + (j_start - first_column), j_length);
				}
			}
		}
	}

	/** {@inherit_doc} */
	//override
	public Field_Matrix<T> get_row_matrix(const int row)

	{
		check_row_index(row);
		const BlockField_Matrix<T> out = BlockField_Matrix<>(get_field(), 1, columns);

		// perform copy block-wise, to ensure good cache behavior
		const int i_block = row / BLOCK_SIZE;
		const int i_row = row - i_block * BLOCK_SIZE;
		int out_block_index = 0;
		int out_index = 0;
		std::vector<T> out_block = out.blocks[out_block_index];
		for (int j_block = 0; j_block < block_columns; ++j_block)
		{
			const int j_width = block_width(j_block);
			const std::vector<T> block = blocks[i_block * block_columns + j_block];
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
	public void set_row_matrix(const int& row, const Field_Matrix<T> matrix)

	{
		if (matrix instanceof BlockField_Matrix)
		{
			set_row_matrix(row, (BlockField_Matrix<T>) matrix);
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
	 * @ if the matrix dimensions do
	 * not match one instance row.
	 * @ if the specified row index is invalid.
	 */
	public void set_row_matrix(const int& row, const BlockField_Matrix<T> matrix)

	{
		check_row_index(row);
		const int& n_cols = get_column_dimension();
		if ((matrix.get_row_dimension() != 1) ||
			(matrix.get_column_dimension() != n_cols))
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH_2x2, matrix.get_row_dimension(), matrix.get_column_dimension(), 1, n_cols);
		}

		// perform copy block-wise, to ensure good cache behavior
		const int i_block = row / BLOCK_SIZE;
		const int i_row = row - i_block * BLOCK_SIZE;
		int m_block_index = 0;
		int m_index = 0;
		std::vector<T> m_block = matrix.blocks[m_block_index];
		for (int j_block = 0; j_block < block_columns; ++j_block)
		{
			const int j_width = block_width(j_block);
			const std::vector<T> block = blocks[i_block * block_columns + j_block];
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
	public Field_Matrix<T> get_column_matrix(const int column)

	{
		check_column_index(column);
		const BlockField_Matrix<T> out = BlockField_Matrix<>(get_field(), rows, 1);

		// perform copy block-wise, to ensure good cache behavior
		const int j_block = column / BLOCK_SIZE;
		const int j_column = column - j_block * BLOCK_SIZE;
		const int j_width = block_width(j_block);
		int out_block_index = 0;
		int out_index = 0;
		std::vector<T> out_block = out.blocks[out_block_index];
		for (const int& i_block = 0; i_block < block_rows; ++i_block)
		{
			const int i_height = block_height(i_block);
			const std::vector<T> block = blocks[i_block * block_columns + j_block];
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
	public void set_column_matrix(const int& column, const Field_Matrix<T> matrix)

	{
		if (matrix instanceof BlockField_Matrix)
		{
			set_column_matrix(column, (BlockField_Matrix<T>) matrix);
		}
		else
		{
			super.set_column_matrix(column, matrix);
		}
	}

	/**
	 * Sets the entries in column number {@code column}
	 * as a column matrix.  Column indices start at 0.
	 *
	 * @param column Column to be set.
	 * @param matrix Column matrix (must have one column and the same number of rows
	 * as the instance).
	 * @ if the matrix dimensions do
	 * not match one instance column.
	 * @ if the specified column index is invalid.
	 */
	void set_column_matrix(const int& column, const BlockField_Matrix<T> matrix)

	{
		check_column_index(column);
		const int& n_rows = get_row_dimension();
		if ((matrix.get_row_dimension() != n_rows) ||
			(matrix.get_column_dimension() != 1))
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH_2x2, matrix.get_row_dimension(), matrix.get_column_dimension(), n_rows, 1);
		}

		// perform copy block-wise, to ensure good cache behavior
		const int j_block = column / BLOCK_SIZE;
		const int j_column = column - j_block * BLOCK_SIZE;
		const int j_width = block_width(j_block);
		int m_block_index = 0;
		int m_index = 0;
		std::vector<T> m_block = matrix.blocks[m_block_index];
		for (const int& i_block = 0; i_block < block_rows; ++i_block)
		{
			const int i_height = block_height(i_block);
			const std::vector<T> block = blocks[i_block * block_columns + j_block];
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
	public Field_Vector<T> get_row_vector(const int row)

	{
		check_row_index(row);
		const std::vector<T> out_data = Math_Arrays::build_array(get_field(), columns);

		// perform copy block-wise, to ensure good cache behavior
		const int i_block = row / BLOCK_SIZE;
		const int i_row = row - i_block * BLOCK_SIZE;
		int out_index = 0;
		for (int j_block = 0; j_block < block_columns; ++j_block)
		{
			const int j_width = block_width(j_block);
			const std::vector<T> block = blocks[i_block * block_columns + j_block];
			System.arraycopy(block, i_row * j_width, out_data, out_index, j_width);
			out_index += j_width;
		}

		return ArrayField_Vector<T>(get_field(), out_data, false);
	}

	/** {@inherit_doc} */
	//override
	public void set_row_vector(const int& row, const Field_Vector<T> vector)

	{
		if (vector instanceof ArrayField_Vector)
		{
			set_row(row, ((ArrayField_Vector<T>) vector).get_data_ref());
		}
		else
		{
			super.set_row_vector(row, vector);
		}
	}

	/** {@inherit_doc} */
	//override
	public Field_Vector<T> get_column_vector(const int column)

	{
		check_column_index(column);
		const std::vector<T> out_data = Math_Arrays::build_array(get_field(), rows);

		// perform copy block-wise, to ensure good cache behavior
		const int j_block = column / BLOCK_SIZE;
		const int j_column = column - j_block * BLOCK_SIZE;
		const int j_width = block_width(j_block);
		int out_index = 0;
		for (const int& i_block = 0; i_block < block_rows; ++i_block)
		{
			const int i_height = block_height(i_block);
			const std::vector<T> block = blocks[i_block * block_columns + j_block];
			for (int i{}; i < i_height; ++i)
			{
				out_data[out_index++] = block[i * j_width + j_column];
			}
		}

		return ArrayField_Vector<T>(get_field(), out_data, false);
	}

	/** {@inherit_doc} */
	//override
	public void set_column_vector(const int& column, const Field_Vector<T> vector)

	{
		if (vector instanceof ArrayField_Vector)
		{
			set_column(column, ((ArrayField_Vector<T>) vector).get_data_ref());
		}
		else
		{
			super.set_column_vector(column, vector);
		}
	}

	/** {@inherit_doc} */
	//override
	public std::vector<T> get_row(const int row)
	{
		check_row_index(row);
		const std::vector<T> out = Math_Arrays::build_array(get_field(), columns);

		// perform copy block-wise, to ensure good cache behavior
		const int i_block = row / BLOCK_SIZE;
		const int i_row = row - i_block * BLOCK_SIZE;
		int out_index = 0;
		for (int j_block = 0; j_block < block_columns; ++j_block)
		{
			const int j_width = block_width(j_block);
			const std::vector<T> block = blocks[i_block * block_columns + j_block];
			System.arraycopy(block, i_row * j_width, out, out_index, j_width);
			out_index += j_width;
		}

		return out;
	}

	/** {@inherit_doc} */
	//override
	public void set_row(const int& row, const std::vector<T> array)

	{
		check_row_index(row);
		const int& n_cols = get_column_dimension();
		if (array.size() != n_cols)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH_2x2, 1, array.size(), 1, n_cols);
		}

		// perform copy block-wise, to ensure good cache behavior
		const int i_block = row / BLOCK_SIZE;
		const int i_row = row - i_block * BLOCK_SIZE;
		int out_index = 0;
		for (int j_block = 0; j_block < block_columns; ++j_block)
		{
			const int j_width = block_width(j_block);
			const std::vector<T> block = blocks[i_block * block_columns + j_block];
			System.arraycopy(array, out_index, block, i_row * j_width, j_width);
			out_index += j_width;
		}
	}

	/** {@inherit_doc} */
	//override
	public std::vector<T> get_column(const int column)
	{
		check_column_index(column);
		const std::vector<T> out = Math_Arrays::build_array(get_field(), rows);

		// perform copy block-wise, to ensure good cache behavior
		const int j_block = column / BLOCK_SIZE;
		const int j_column = column - j_block * BLOCK_SIZE;
		const int j_width = block_width(j_block);
		int out_index = 0;
		for (const int& i_block = 0; i_block < block_rows; ++i_block)
		{
			const int i_height = block_height(i_block);
			const std::vector<T> block = blocks[i_block * block_columns + j_block];
			for (int i{}; i < i_height; ++i)
			{
				out[out_index++] = block[i * j_width + j_column];
			}
		}

		return out;
	}

	/** {@inherit_doc} */
	//override
	public void set_column(const int& column, const std::vector<T> array)

	{
		check_column_index(column);
		const int& n_rows = get_row_dimension();
		if (array.size() != n_rows)
		{
			throw std::exception("not implemented");
			// throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH_2x2, array.size(), 1, n_rows, 1);
		}

		// perform copy block-wise, to ensure good cache behavior
		const int j_block = column / BLOCK_SIZE;
		const int j_column = column - j_block * BLOCK_SIZE;
		const int j_width = block_width(j_block);
		int out_index = 0;
		for (const int& i_block = 0; i_block < block_rows; ++i_block)
		{
			const int i_height = block_height(i_block);
			const std::vector<T> block = blocks[i_block * block_columns + j_block];
			for (int i{}; i < i_height; ++i)
			{
				block[i * j_width + j_column] = array[out_index++];
			}
		}
	}

	/** {@inherit_doc} */
	//override
	public T get_entry(const int& row, const int column)

	{
		check_row_index(row);
		check_column_index(column);

		const int i_block = row / BLOCK_SIZE;
		const int j_block = column / BLOCK_SIZE;
		const int& k = (row - i_block * BLOCK_SIZE) * block_width(j_block) +
			(column - j_block * BLOCK_SIZE);

		return blocks[i_block * block_columns + j_block][k];
	}

	/** {@inherit_doc} */
	//override
	public void set_entry(const int& row, const int& column, const T value)

	{
		check_row_index(row);
		check_column_index(column);

		const int i_block = row / BLOCK_SIZE;
		const int j_block = column / BLOCK_SIZE;
		const int& k = (row - i_block * BLOCK_SIZE) * block_width(j_block) +
			(column - j_block * BLOCK_SIZE);

		blocks[i_block * block_columns + j_block][k] = value;
	}

	/** {@inherit_doc} */
	//override
	public void add_to_entry(const int& row, const int& column, const T increment)

	{
		check_row_index(row);
		check_column_index(column);

		const int i_block = row / BLOCK_SIZE;
		const int j_block = column / BLOCK_SIZE;
		const int& k = (row - i_block * BLOCK_SIZE) * block_width(j_block) +
			(column - j_block * BLOCK_SIZE);
		const std::vector<T> block_i_j = blocks[i_block * block_columns + j_block];

		block_i_j[k] = block_i_j[k].add(increment);
	}

	/** {@inherit_doc} */
	//override
	public void multiply_entry(const int& row, const int& column, const T factor)

	{
		check_row_index(row);
		check_column_index(column);

		const int i_block = row / BLOCK_SIZE;
		const int j_block = column / BLOCK_SIZE;
		const int& k = (row - i_block * BLOCK_SIZE) * block_width(j_block) +
			(column - j_block * BLOCK_SIZE);
		const std::vector<T> block_i_j = blocks[i_block * block_columns + j_block];

		block_i_j[k] = block_i_j[k].multiply(factor);
	}

	/** {@inherit_doc} */
	//override
	public Field_Matrix<T> transpose()
	{
		const int& n_rows = get_row_dimension();
		const int& n_cols = get_column_dimension();
		const BlockField_Matrix<T> out = BlockField_Matrix<>(get_field(), n_cols, n_rows);

		// perform transpose block-wise, to ensure good cache behavior
		int block_index = 0;
		for (const int& i_block = 0; i_block < block_columns; ++i_block)
		{
			for (int j_block = 0; j_block < block_rows; ++j_block)
			{
				// transpose current block
				const std::vector<T> out_block = out.blocks[block_index];
				const std::vector<T> t_block = blocks[j_block * block_columns + i_block];
				const int      p_start = i_block * BLOCK_SIZE;
				const int      p_end = std::min(p_start + BLOCK_SIZE, columns);
				const int      q_start = j_block * BLOCK_SIZE;
				const int      q_end = std::min(q_start + BLOCK_SIZE, rows);
				int k = 0;
				for (const int& p = p_start; p < p_end; ++p)
				{
					const int l_inc = p_end - p_start;
					int l = p - p_start;
					for (const int& q = q_start; q < q_end; ++q)
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
	public int get_row_dimension()
	{
		return rows;
	}

	/** {@inherit_doc} */
	//override
	public int get_column_dimension()
	{
		return columns;
	}

	/** {@inherit_doc} */
	//override
	public std::vector<T> operate(const std::vector<T> v)
	{
		if (v.size() != columns)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, v.size(), columns);
		}
		const std::vector<T> out = Math_Arrays::build_array(get_field(), rows);
		const T zero = get_field().get_zero();

		// perform multiplication block-wise, to ensure good cache behavior
		for (const int& i_block = 0; i_block < block_rows; ++i_block)
		{
			const int p_start = i_block * BLOCK_SIZE;
			const int p_end = std::min(p_start + BLOCK_SIZE, rows);
			for (int j_block = 0; j_block < block_columns; ++j_block)
			{
				const std::vector<T> block = blocks[i_block * block_columns + j_block];
				const int      q_start = j_block * BLOCK_SIZE;
				const int      q_end = std::min(q_start + BLOCK_SIZE, columns);
				int k = 0;
				for (const int& p = p_start; p < p_end; ++p)
				{
					T sum = zero;
					int q = q_start;
					while (q < q_end - 3)
					{
						sum = sum.
							add(block[k].multiply(v[q])).
							add(block[k + 1].multiply(v[q + 1])).
							add(block[k + 2].multiply(v[q + 2])).
							add(block[k + 3].multiply(v[q + 3]));
						k += 4;
						q += 4;
					}
					while (q < q_end)
					{
						sum = sum.add(block[k++].multiply(v[q++]));
					}
					out[p] = out[p].add(sum);
				}
			}
		}

		return out;
	}

	/** {@inherit_doc} */
	//override
	public std::vector<T> pre_multiply(const std::vector<T> v)
	{
		if (v.size() != rows)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, v.size(), rows);
		}
		const std::vector<T> out = Math_Arrays::build_array(get_field(), columns);
		const T zero = get_field().get_zero();

		// perform multiplication block-wise, to ensure good cache behavior
		for (int j_block = 0; j_block < block_columns; ++j_block)
		{
			const int j_width = block_width(j_block);
			const int j_width2 = j_width + j_width;
			const int j_width3 = j_width2 + j_width;
			const int j_width4 = j_width3 + j_width;
			const int q_start = j_block * BLOCK_SIZE;
			const int q_end = std::min(q_start + BLOCK_SIZE, columns);
			for (const int& i_block = 0; i_block < block_rows; ++i_block)
			{
				const std::vector<T> block = blocks[i_block * block_columns + j_block];
				const int      p_start = i_block * BLOCK_SIZE;
				const int      p_end = std::min(p_start + BLOCK_SIZE, rows);
				for (const int& q = q_start; q < q_end; ++q)
				{
					int k = q - q_start;
					T sum = zero;
					int p = p_start;
					while (p < p_end - 3)
					{
						sum = sum.
							add(block[k].multiply(v[p])).
							add(block[k + j_width].multiply(v[p + 1])).
							add(block[k + j_width2].multiply(v[p + 2])).
							add(block[k + j_width3].multiply(v[p + 3]));
						k += j_width4;
						p += 4;
					}
					while (p < p_end)
					{
						sum = sum.add(block[k].multiply(v[p++]));
						k += j_width;
					}
					out[q] = out[q].add(sum);
				}
			}
		}

		return out;
	}

	/** {@inherit_doc} */
	//override
	public T walk_in_row_order(const Field_Matrix_Changing_Visitor<T> visitor)
	{
		visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
		for (const int& i_block = 0; i_block < block_rows; ++i_block)
		{
			const int p_start = i_block * BLOCK_SIZE;
			const int p_end = std::min(p_start + BLOCK_SIZE, rows);
			for (const int& p = p_start; p < p_end; ++p)
			{
				for (int j_block = 0; j_block < block_columns; ++j_block)
				{
					const int j_width = block_width(j_block);
					const int q_start = j_block * BLOCK_SIZE;
					const int q_end = std::min(q_start + BLOCK_SIZE, columns);
					const std::vector<T> block = blocks[i_block * block_columns + j_block];
					int k = (p - p_start) * j_width;
					for (const int& q = q_start; q < q_end; ++q)
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
	public T walk_in_row_order(const Field_Matrix_Preserving_Visitor<T> visitor)
	{
		visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
		for (const int& i_block = 0; i_block < block_rows; ++i_block)
		{
			const int p_start = i_block * BLOCK_SIZE;
			const int p_end = std::min(p_start + BLOCK_SIZE, rows);
			for (const int& p = p_start; p < p_end; ++p)
			{
				for (int j_block = 0; j_block < block_columns; ++j_block)
				{
					const int j_width = block_width(j_block);
					const int q_start = j_block * BLOCK_SIZE;
					const int q_end = std::min(q_start + BLOCK_SIZE, columns);
					const std::vector<T> block = blocks[i_block * block_columns + j_block];
					int k = (p - p_start) * j_width;
					for (const int& q = q_start; q < q_end; ++q)
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
	public T walk_in_row_order(const Field_Matrix_Changing_Visitor<T> visitor, const int& start_row, const int& end_row, const int& start_column, const int& end_column)

	{
		check_sub_matrix_index(start_row, end_row, start_column, end_column);
		visitor.start(rows, columns, start_row, end_row, start_column, end_column);
		for (const int& i_block = start_row / BLOCK_SIZE; i_block < 1 + end_row / BLOCK_SIZE; ++i_block)
		{
			const int p0 = i_block * BLOCK_SIZE;
			const int p_start = std::max(start_row, p0);
			const int p_end = std::min((i_block + 1) * BLOCK_SIZE, 1 + end_row);
			for (const int& p = p_start; p < p_end; ++p)
			{
				for (int j_block = start_column / BLOCK_SIZE; j_block < 1 + end_column / BLOCK_SIZE; ++j_block)
				{
					const int j_width = block_width(j_block);
					const int q0 = j_block * BLOCK_SIZE;
					const int q_start = std::max(start_column, q0);
					const int q_end = std::min((j_block + 1) * BLOCK_SIZE, 1 + end_column);
					const std::vector<T> block = blocks[i_block * block_columns + j_block];
					int k = (p - p0) * j_width + q_start - q0;
					for (const int& q = q_start; q < q_end; ++q)
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
	public T walk_in_row_order(const Field_Matrix_Preserving_Visitor<T> visitor, const int& start_row, const int& end_row, const int& start_column, const int& end_column)

	{
		check_sub_matrix_index(start_row, end_row, start_column, end_column);
		visitor.start(rows, columns, start_row, end_row, start_column, end_column);
		for (const int& i_block = start_row / BLOCK_SIZE; i_block < 1 + end_row / BLOCK_SIZE; ++i_block)
		{
			const int p0 = i_block * BLOCK_SIZE;
			const int p_start = std::max(start_row, p0);
			const int p_end = std::min((i_block + 1) * BLOCK_SIZE, 1 + end_row);
			for (const int& p = p_start; p < p_end; ++p)
			{
				for (int j_block = start_column / BLOCK_SIZE; j_block < 1 + end_column / BLOCK_SIZE; ++j_block)
				{
					const int j_width = block_width(j_block);
					const int q0 = j_block * BLOCK_SIZE;
					const int q_start = std::max(start_column, q0);
					const int q_end = std::min((j_block + 1) * BLOCK_SIZE, 1 + end_column);
					const std::vector<T> block = blocks[i_block * block_columns + j_block];
					int k = (p - p0) * j_width + q_start - q0;
					for (const int& q = q_start; q < q_end; ++q)
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
	public T walk_in_optimized_order(const Field_Matrix_Changing_Visitor<T> visitor)
	{
		visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
		int block_index = 0;
		for (const int& i_block = 0; i_block < block_rows; ++i_block)
		{
			const int p_start = i_block * BLOCK_SIZE;
			const int p_end = std::min(p_start + BLOCK_SIZE, rows);
			for (int j_block = 0; j_block < block_columns; ++j_block)
			{
				const int q_start = j_block * BLOCK_SIZE;
				const int q_end = std::min(q_start + BLOCK_SIZE, columns);
				const std::vector<T> block = blocks[block_index];
				int k = 0;
				for (const int& p = p_start; p < p_end; ++p)
				{
					for (const int& q = q_start; q < q_end; ++q)
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
	public T walk_in_optimized_order(const Field_Matrix_Preserving_Visitor<T> visitor)
	{
		visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
		int block_index = 0;
		for (const int& i_block = 0; i_block < block_rows; ++i_block)
		{
			const int p_start = i_block * BLOCK_SIZE;
			const int p_end = std::min(p_start + BLOCK_SIZE, rows);
			for (int j_block = 0; j_block < block_columns; ++j_block)
			{
				const int q_start = j_block * BLOCK_SIZE;
				const int q_end = std::min(q_start + BLOCK_SIZE, columns);
				const std::vector<T> block = blocks[block_index];
				int k = 0;
				for (const int& p = p_start; p < p_end; ++p)
				{
					for (const int& q = q_start; q < q_end; ++q)
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
	public T walk_in_optimized_order(const Field_Matrix_Changing_Visitor<T> visitor, const int& start_row, const int& end_row, const int& start_column, const int& end_column)

	{
		check_sub_matrix_index(start_row, end_row, start_column, end_column);
		visitor.start(rows, columns, start_row, end_row, start_column, end_column);
		for (const int& i_block = start_row / BLOCK_SIZE; i_block < 1 + end_row / BLOCK_SIZE; ++i_block)
		{
			const int p0 = i_block * BLOCK_SIZE;
			const int p_start = std::max(start_row, p0);
			const int p_end = std::min((i_block + 1) * BLOCK_SIZE, 1 + end_row);
			for (int j_block = start_column / BLOCK_SIZE; j_block < 1 + end_column / BLOCK_SIZE; ++j_block)
			{
				const int j_width = block_width(j_block);
				const int q0 = j_block * BLOCK_SIZE;
				const int q_start = std::max(start_column, q0);
				const int q_end = std::min((j_block + 1) * BLOCK_SIZE, 1 + end_column);
				const std::vector<T> block = blocks[i_block * block_columns + j_block];
				for (const int& p = p_start; p < p_end; ++p)
				{
					int k = (p - p0) * j_width + q_start - q0;
					for (const int& q = q_start; q < q_end; ++q)
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
	public T walk_in_optimized_order(const Field_Matrix_Preserving_Visitor<T> visitor, const int& start_row, const int& end_row, const int& start_column, const int& end_column)

	{
		check_sub_matrix_index(start_row, end_row, start_column, end_column);
		visitor.start(rows, columns, start_row, end_row, start_column, end_column);
		for (const int& i_block = start_row / BLOCK_SIZE; i_block < 1 + end_row / BLOCK_SIZE; ++i_block)
		{
			const int p0 = i_block * BLOCK_SIZE;
			const int p_start = std::max(start_row, p0);
			const int p_end = std::min((i_block + 1) * BLOCK_SIZE, 1 + end_row);
			for (int j_block = start_column / BLOCK_SIZE; j_block < 1 + end_column / BLOCK_SIZE; ++j_block)
			{
				const int j_width = block_width(j_block);
				const int q0 = j_block * BLOCK_SIZE;
				const int q_start = std::max(start_column, q0);
				const int q_end = std::min((j_block + 1) * BLOCK_SIZE, 1 + end_column);
				const std::vector<T> block = blocks[i_block * block_columns + j_block];
				for (const int& p = p_start; p < p_end; ++p)
				{
					int k = (p - p0) * j_width + q_start - q0;
					for (const int& q = q_start; q < q_end; ++q)
					{
						visitor.visit(p, q, block[k]);
						++k;
					}
				}
			}
		}
		return visitor.end();
	}

	/**
	 * Get the height of a block.
	 * @param block_row row index (in block sense) of the block
	 * @return height (number of rows) of the block
	 */
	private int block_height(const int& block_row)
	{
		return (block_row == block_rows - 1) ? rows - block_row * BLOCK_SIZE : BLOCK_SIZE;
	}

	/**
	 * Get the width of a block.
	 * @param block_column column index (in block sense) of the block
	 * @return width (number of columns) of the block
	 */
	private int block_width(const int& block_column)
	{
		return (block_column == block_columns - 1) ? columns - block_column * BLOCK_SIZE : BLOCK_SIZE;
	}
}
