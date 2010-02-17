#include "SimplexCell.hpp"
#include "Simplex.hpp"

int SimplexCell::getColor() const
{
    return simp->getColor();
}

void SimplexCell::setColor(int col)
{
    simp->setColor(col);
}

size_t SimplexCell::getDim()
{
    return simp->getDim();
}
