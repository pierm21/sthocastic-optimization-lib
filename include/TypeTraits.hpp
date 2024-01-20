#pragma once

#include <Eigen/Dense>
#include <functional>

template <size_t dim>
using RealVector = Eigen::Matrix<double, dim, 1>;
template <size_t dim>
using RealFunction = std::function<double(const RealVector<dim>&)>;