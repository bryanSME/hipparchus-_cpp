#pragma once

#include <vector>

/** Local iterator for sub-intervals. */
class Sub_Intervals_Iterator : std::vector<double>::Iterator
{
private:
	/** Current node. */
	BSP_Tree<Euclidean_1D> my_current;

	/** Sub-interval no yet returned. */
	std::vector<double> my_pending;

	/** 
	 * Walk the tree to select the pending sub-interval.
	 */
	void select_pending()
	{
		// look for the start of the interval
		BSP_Tree<Euclidean_1D> start = my_current;
		while (start != NULL && !is_interval_start(start))
		{
			start = next_internal_node(start);
		}

		if (start == NULL)
		{
			// we have exhausted the iterator
			my_current = NULL;
			my_pending = NULL;
			return;
		}

		// look for the end of the interval
		BSP_Tree<Euclidean_1D> end = start;
		while (end != NULL && !is_interval_end(end))
		{
			end = next_internal_node(end);
		}

		if (end != NULL)
		{
			// we have identified the interval
			my_pending = std::vector<double>
			{
				get_angle(start), get_angle(end)
			};

			// prepare search for next interval
			my_current = end;
			return;
		}
		
		// the const interval is open toward infinity
		my_pending = std::vector<double>
		{
			get_angle(start), std::numeric_limits<double>::infinity()
		};

		// there won't be any other intervals
		my_current = NULL;
	}

public:
	/** Simple constructor.
	 */
	Sub_Intervals_Iterator()
	{
		my_current = get_first_interval_boundary();

		if (my_current == NULL)
		{
			// all the leaf tree nodes share the same inside/outside status
			if ((Boolean)get_first_leaf(get_tree(false)).get_attribute())
			{
				// it is an inside node, it represents the full real line
				my_pending = std::vector<double>
				{
					-std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()
				};
			}
			else
			{
				pending = NULL;
			}
		}
		else if (is_interval_end(current))
		{
			// the first boundary is an interval end, // so the first interval starts at infinity
			my_pending = std::vector<double>
			{
				-std::numeric_limits<double>::infinity(), get_angle(my_current)
			};
		}
		else
		{
			select_pending();
		}
	}

	/** {@inherit_doc} */
	//override
	bool has_next() const
	{
		return my_pending != NULL;
	}

	/** {@inherit_doc} */
	//override
	std::vector<double> next()
	{
		if (my_pending == NULL)
		{
			throw std::exception("not implemented");
			//throw No_Such_Element_Exception();
		}
		const std::vector<double> next = my_pending;
		select_pending();
		return next;
	}

	/** {@inherit_doc} */
	//override
	void remove() const
	{
		throw std::exception("not implemented");
		//throw Unsupported_Operation_Exception();
	}
};