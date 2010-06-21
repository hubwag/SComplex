#ifndef CELLDIMINDEXER_HPP_INCLUDED
#define CELLDIMINDEXER_HPP_INCLUDED

#include <vector>

// No filtering for now, so it requires manuas skipping of cells of
// color != 1
template<typename CellT>
class CellDimIndexer
{
	std::vector<std::vector<CellT> > cellsByDim;
	std::vector<CellT> sortedCells;

	struct IsAvailable
	{
		bool operator()(CellT c) const
		{
			return c.getColor() == 1;
		}
	};

	public:

	typedef typename std::vector<CellT>::const_iterator iterator;

	CellDimIndexer()
	{
	}

	template<typename iterT>
	CellDimIndexer(iterT begin, iterT end, int maxDim)
	{
		cellsByDim.resize(maxDim + 1, std::vector<CellT>());

		for (; begin != end; ++begin)
		{
			cellsByDim[begin->getDim()].push_back(*begin);
		}

		for (size_t d = 0; d < cellsByDim.size(); d++)
		{
			const int n = cellsByDim[d].size();
			for (size_t i = 0; i < n; i++)
			{
				sortedCells.push_back(cellsByDim[d][i]);
			}
		}
	}

	iterator begin()
	{
		return sortedCells.begin();
	}

	iterator end()
	{
		return sortedCells.end();
	}
};

#endif
