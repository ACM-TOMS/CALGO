//This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
//To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/
//or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

#include "ADT2D.h"
#include <tuple>

using namespace std;

ADT2D::ADT2D(ADT2D&& other) {
    root = other.root;
    other.root = nullptr;
}

ADT2D::~ADT2D() {
    if(root != nullptr)
        delete root;
}

void ADT2D::store(int id, const double* x) {
    if(root != nullptr)
        store(root, id, x);
    else root = new ADT2DElement{4, {{0.0, 0.0, 0.0, 0.0}}, {{1.0, 1.0, 1.0, 1.0}}, id, x};
}

std::vector<int> ADT2D::retrieve(const double* extent) const {
    vector<int> ids;
    if(root == nullptr)
        return ids;

    array<double, 4> a;
    array<double, 4> b;
    tie(a, b) = createHyperRectangle(extent);
    retrieve(root, ids, a, b);
    return ids;
}

void ADT2D::removeFirst(int id, const double* extent) {
    if(root == nullptr)
        return;

    array<double, 4> a;
    array<double, 4> b;
    tie(a, b) = createHyperRectangle(extent);
    removeFirst(root, nullptr, id, a, b);
}

void ADT2D::store(ADT2DElement* element, int id, const double* x) {
    if(determineChild(element, x) == 0) {
        if(element->left_child != nullptr) {
            store(element->left_child, id, x);
        } else {
            array<double, 4> x_max = element->x_max;
            int split_axis = element->level % 2;
            x_max[split_axis] = 0.5 * (element->x_min[split_axis] + element->x_max[split_axis]);
            element->left_child = new ADT2DElement{element->level + 1, element->x_min, x_max, id, x};
        }
    } else {
        if(element->right_child != nullptr) {
            store(element->right_child, id, x);
        } else {
            array<double, 4> x_min = element->x_min;
            int split_axis = element->level % 4;
            x_min[split_axis] = 0.5 * (element->x_min[split_axis] + element->x_max[split_axis]);
            element->right_child = new ADT2DElement{element->level + 1, x_min, element->x_max, id, x};
        }
    }
}

void ADT2D::retrieve(ADT2DElement* element, std::vector<int>& ids, const std::array<double, 4>& a,
                     const std::array<double, 4>& b) const {
    if(not element->containsHyperRectangle(a, b))
        return;
    if(element->hyperRectangleContainsObject(a, b))
        ids.push_back(element->id);

    if(element->left_child != nullptr)
        retrieve(element->left_child, ids, a, b);
    if(element->right_child != nullptr)
        retrieve(element->right_child, ids, a, b);
}

void ADT2D::removeFirst(ADT2DElement* element, ADT2DElement* parent, int id, const std::array<double, 4>& a,
                        const std::array<double, 4>& b) {
    if(not element->containsHyperRectangle(a, b))
        return;
    if(element->hyperRectangleContainsObject(a, b) and element->id == id) {
        replaceElementWithLeaf(element, parent);
        return;
    }

    if(element->left_child != nullptr)
        removeFirst(element->left_child, element, id, a, b);
    if(element->right_child != nullptr)
        removeFirst(element->right_child, element, id, a, b);
}

std::tuple<std::array<double, 4>, std::array<double, 4>> ADT2D::createHyperRectangle(const double* extent) const {
    array<double, 4> a;
    array<double, 4> b;
    double dx = extent[2] - extent[0];
    double dy = extent[3] - extent[1];

    a[0] = 0.0;
    a[1] = 0.0;
    a[2] = extent[2] - dx;
    a[3] = extent[3] - dy;

    b[0] = extent[0] + dx;
    b[1] = extent[1] + dy;
    b[2] = 1.0;
    b[3] = 1.0;

    return make_tuple(a, b);
}

int ADT2D::determineChild(ADT2DElement* element, double const* x) {
    int split_axis = element->level % 4;
    return x[split_axis] > 0.5 * (element->x_max[split_axis] + element->x_min[split_axis]);
}

void ADT2D::replaceElementWithLeaf(ADT2DElement* element, ADT2DElement* parent) {
    ADT2DElement* leaf;
    tie(leaf, parent) = findFirstLeaf(element, parent);
    element->id = leaf->id;
    element->object = leaf->object;

    if(parent == nullptr)
        root = nullptr;
    else if(parent->left_child == leaf)
        parent->left_child = nullptr;
    else parent->right_child = nullptr;

    delete leaf;
}

std::tuple<ADT2DElement*, ADT2DElement*> ADT2D::findFirstLeaf(ADT2DElement* element, ADT2DElement* parent) {
    if(isLeaf(element))
        return std::make_tuple(element, parent);
    if(element->left_child != nullptr)
        return findFirstLeaf(element->left_child, element);
    return findFirstLeaf(element->right_child, element);
}

bool ADT2D::isLeaf(const ADT2DElement* element) const {
    return element->left_child == nullptr and element->right_child == nullptr;
}
