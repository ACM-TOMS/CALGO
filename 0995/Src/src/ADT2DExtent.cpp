//This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
//To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/
//or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

#include "ADT2DExtent.h"

ADT2DExtent::ADT2DExtent(const Extent& domain) : space_transformer(domain) {}

void ADT2DExtent::store(int id, const Extent& extent) {
    Extent adt_extent;
    adt_extent.lo = space_transformer.toUnitSpace(extent.lo);
    adt_extent.hi = space_transformer.toUnitSpace(extent.hi);
    adt.store(id, adt_extent.lo.data());
}

std::vector<int> ADT2DExtent::retrieve(const Extent& domain) const {
    if(not space_transformer.getDomain().contains(domain))
        return std::vector<int>();

    Extent search_domain(space_transformer.toUnitSpace(domain.lo), space_transformer.toUnitSpace(domain.hi));
    return adt.retrieve(search_domain.lo.data());
}

void ADT2DExtent::removeFirst(int id, const Extent& extent) {
    if(not space_transformer.getDomain().contains(extent))
        return;

    Extent search_domain(space_transformer.toUnitSpace(extent.lo), space_transformer.toUnitSpace(extent.hi));
    adt.removeFirst(id, search_domain.lo.data());
}
