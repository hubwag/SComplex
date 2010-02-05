#include "SimplexCell.hpp"
#include "Simplex.hpp"

int SimplexCell::getColor() const
{
    return simp->getColor();
}

/*
template<>
void SimplexCell::setColor<2>()
{
    return simp->setColor<2>();
}
*/

void SimplexCell::setColor(int col)
{
    simp->setColor(col);
}

size_t SimplexCell::getDim()
{
    return simp->getDim();
}
