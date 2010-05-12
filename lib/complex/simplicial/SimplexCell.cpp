#include <redHom/complex/simplicial/SimplexCell.hpp>
#include <redHom/complex/simplicial/Simplex.hpp>

int SimplexCell::getColor() const
{
    return simp->getColor();
}

void SimplexCell::setColor(int col)
{
    simp->setColor(col);
}

size_t SimplexCell::getDim() const
{
    return simp->getDim();
}

bool SimplexCell::operator<(const SimplexCell& b) const {
  return (this->simp->nrs) < (b.simp->nrs);
}
