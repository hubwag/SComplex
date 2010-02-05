#ifndef SIMPLEXCELL_HPP
#define SIMPLEXCELL_HPP

#include <cstdlib>

struct Simplex;

class SimplexCell
{
	Simplex *simp;
public:

	explicit SimplexCell(Simplex *simplex = 0) : simp(simplex)
	{
	}

	template<typename dummy>
	SimplexCell(dummy &) : simp(0)
	{
	}

	operator Simplex&()
	{
		return *simp;
	}

	operator const Simplex&() const
	{
		return *simp;
	}

	const Simplex& getImpl() const
	{
		return *simp;
	}

	Simplex& getImpl()
	{
		return *simp;
	}

	SimplexCell& operator=(Simplex &s);

	int getColor() const;

    void setColor(int col);

    template<int color>
    void setColor();

	size_t getDim();
};

#endif
