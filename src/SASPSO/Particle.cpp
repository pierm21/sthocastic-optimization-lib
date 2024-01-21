#include "SASPSO/Particle.hpp"

using namespace type_traits;

template <size_t dim>
void Particle<dim>::initialize()
{
    // Create the uniform double generator in the range [lower_bound_i, upper_bound_i] for each dimension i
    for(size_t i = 0; i < dim; ++i)
    {
        std::uniform_real_distribution<double> distr(lower_bound_[i], upper_bound_[i]);
        // Initialize the position and velocity vectors
        position_[i] = distr(random_generator_);
        velocity_[i] = distr(random_generator_);
    }
    // Initialize the best position
    best_position_ = position_;
    // Initialize the best value
    best_value_ = problem_->get_fitness_function()(position_);
    // Initialize the constraint violation
    update_constraint_violation();
}

template <std::size_t dim>
void Particle<dim>::update(const RealVector<dim> &global_best_position, int iteration, int max_iter)
{
    double beta = (global_best_position - best_position_).norm();
    // Compute omega according to the current iteration
    double delta = (omega_s_ - omega_f_) / max_iter;
    double omega = (omega_s_ - omega_f_) * exp(-delta * iteration / beta) + omega_f_;
    // Compute phi1 and phi2 according to the current iteration
    delta = (phi1_s_ - phi1_f_) / max_iter;
    double phi1 = (phi1_s_ - phi1_f_) * exp(-delta * iteration / beta) + phi1_f_;
    // Compute phi2 according to the current iteration
    delta = (phi2_s_ - phi2_f_) / max_iter;
    double phi2 = (phi2_s_ - phi2_f_) * exp(-delta * iteration / beta) + phi2_f_;

    // Compute the isobarycenter of the particle
    RealVector<dim> G = position_;
    G += (phi1 * (best_position_ - position_) + phi2 * (global_best_position - position_))/3.0;

    // Generate a point in the hypersphere of radius ||G-position|| centered in G
    std::uniform_real_distribution<double> distr(0, 1);
    //  1. Generate a random direction
    RealVector<dim> random_direction = RealVector<dim>::NullaryExpr([&] (int) {return distr(random_generator_);});
    random_direction.normalize();
    //  2. Generate a random radius  with probability proportional to the surface area of a ball with a given radius
    double radius = pow(distr(random_generator_), 1.0/dim);




/*
    //  Create the uniform double generator in the range [lower_bound, upper_bound]
    std::uniform_real_distribution<double> distr(0, 1);
    // Initialize the position and velocity vectors
    std::generate(r1_.begin(), r1_.end(), [&]()
                  { return distr(random_generator_); });
    std::generate(r2_.begin(), r2_.end(), [&]()
                  { return distr(random_generator_); });

    // For each dimension
    for (std::size_t i = 0; i < velocity_.size(); ++i)
    {
        // Velocity update
        velocity_[i] = w * velocity_[i] + c1 * r1_[i] * (best_position_[i] - position_[i]) + c2 * r2_[i] * (global_best_position[i] - position_[i]);

        // Position update
        position_[i] += velocity_[i];
        // Check boundaries
        if (position_[i] < lower_bound_)
            position_[i] = lower_bound_;
        else if (position_[i] > upper_bound_)
            position_[i] = upper_bound_;
    }

    // Update best position if necessary
    double current_value = fitness_function_(position_);
    if (current_value < best_value_)
    {
        best_position_ = position_;
        best_value_ = current_value;
    }

    // Return the value of the fitness function in the best position
    return best_value_;
*/
}

template <std::size_t dim>
void Particle<dim>::print() const
{
    std::cout << "Position:\t(";
    for (auto &i : position_)
    {
        std::cout << i << ", ";
    }
    std::cout << "\b\b)" << std::endl;

    std::cout << "Velocity:\t(";
    for (auto &i : velocity_)
    {
        std::cout << i << ", ";
    }
    std::cout << "\b\b)" << std::endl;

    std::cout << "Best position:\t(";
    for (auto &i : best_position_)
    {
        std::cout << i << ", ";
    }
    std::cout << "\b\b)" << std::endl;

    std::cout << "Best value:\t" << best_value_ << std::endl;
    std::cout << std::endl;
}

template<size_t dim>
void Particle<dim>::update_constraint_violation()
{
    constraint_violation_ = 0.0;
    if(problem_->has_constraints())
    {
    for(RealFunction constraint : problem_->get_equality_constraints())
        constraint_violation += std::abs(constraint(position_));
    for(RealFunction constraint : problem_->get_inquality_constraints())
        constraint_violation += std::max(0, constraint(position_));
    }
}