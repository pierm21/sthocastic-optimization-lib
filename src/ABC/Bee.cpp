#pragma once

#include "ABC/Bee.hpp"

using namespace type_traits;

template <size_t dim>
void Bee<dim>::initialize()
{
    // Draw uniformly a value in the range [lower_bound_i, upper_bound_i] for each dimension i
    for (size_t i = 0; i < dim; ++i)
    {
        std::uniform_real_distribution<double> distr(this->problem_->get_lower_bound(i), this->problem_->get_upper_bound(i));
        // Initialize the position and velocity vectors
        this->position_[i] = distr(*this->random_generator_);
    }
    // Initialize the failure counter
    failure_counter_ = 0;
    // Initialize the value
    cost_value_ = this->problem_->get_fitness_function()(this->position_);
    // Initialize the constraint violation
    constraint_violation_ = compute_constraint_violation(this->position_);
}

template <size_t dim>
void Bee<dim>::update_position(const double MR, double violation_threshold, const std::vector<Bee<dim>>& colony_, double tol)
{
    int change_occurred = 0;
    int neighbour_index = -1;
    std::uniform_real_distribution<double> distr(0, 1.0);
    std::uniform_real_distribution<double> distr2(-1.0, 1.0);
    std::uniform_int_distribution<int> distr3(0, colony_.size() - 1);
    RealVector<dim> new_position;

    for (int i = 0; i < dim; ++i)
    {
        // Generate a random number in the range [0, 1]
        double r = distr(*this->random_generator_);
        if (r < MR)
        {
            change_occurred = 1;
            // Generate a random number in the range [-1, 1], used as a factor to the position update
            double phi = distr2(*this->random_generator_);
            
            // Choose a random neighbour different from the actual bee
            do{
                neighbour_index = distr3(*this->random_generator_);
            } while (r == i);

            new_position[i] = this->position_[i] + phi * (this->position_[i] - colony_[neighbour_index].get_position()[i]);
        }
    }

    if (change_occurred == 0)
    {
        // Choose a random neighbour
        neighbour_index = distr3(*this->random_generator_);

        // Choose a random dimension to be modified
        std::uniform_int_distribution<int> distr4(0, dim - 1);
        int i = distr4(*this->random_generator_); 

        // Generate a random number in the range [-1, 1], used as a factor to the position update
        double phi = distr2(*this->random_generator_);

        // Update position
        new_position[i] = this->position_[i] + phi * (this->position_[i] - colony_[neighbour_index].get_position()[i]);
    }

    // Compute the new position's behaviour
    double new_position_value_ = this->problem_->get_fitness_function()(new_position);
    double new_position_constraint_violation_ = compute_constraint_violation(new_position);

    // Update the position if the new position is better
    if(feasibility_rule(cost_value_, new_position_value_, constraint_violation_, new_position_constraint_violation_, violation_threshold, tol))
    {
        this->position_ = new_position;
        this->cost_value_ = new_position_value_;
        this->constraint_violation_ = new_position_constraint_violation_;
        failure_counter_ = 0;
    }
    else
    {
        failure_counter_++;
    }

}

template <size_t dim>
void Bee<dim>::compute_probability(const double total_fitness_value, const double total_constraint_violation, const double violation_threshold, const int colony_size)
{
    // If the solution is feasible, the probability is proportional to the fitness value
    if (constraint_violation_ <= violation_threshold)
    {
        fitness_probability_ = 0.5 + 0.5 * (fitness_value_ / total_fitness_value);
    }
    // If the solution is infeasible, the probability is inversly proportional to the constraint violation
    else
    {
        fitness_probability_ = 0.5 * (1 - (constraint_violation_ / total_constraint_violation));
    }
}


template <size_t dim>
double Bee<dim>::compute_constraint_violation(const RealVector<dim> &position) const
{
    double violation_ = 0.0;
    if (this->problem_->has_constraints())
    {
        for (RealFunction<dim> constraint : this->problem_->get_equality_constraints())
            violation_ += std::abs(constraint(position));
        for (RealFunction<dim> constraint : this->problem_->get_inequality_constraints())
            violation_ += std::max(0.0, constraint(position));
    }

    return violation_;
}

template <size_t dim>
bool Bee<dim>::feasibility_rule(double value1, double value2, double viol1, double viol2, double violation_threshold, double tol) const
{
    double ub = violation_threshold + tol;
    double lb = std::max(violation_threshold - tol, 0.0);

    // (a) a feasible solution is preferred over an infeasible solution
    if (viol1 <= lb && viol2 > ub)
        return true;
    else if (viol1 > ub && viol2 <= lb)
        return false;
    // (b) among two feasible solutions, the one with better objective function value is preferred
    else if (viol1 <= lb && viol2 <= lb)
        return value1 < value2;
    // (c) among two infeasible solutions, the one with smaller TAV is chosen
    else
        return viol1 < viol2;
}

template <std::size_t dim>
void Bee<dim>::print(std::ostream &out) const
{
    out << "Position: (";
    for (std::size_t i = 0; i < dim; ++i)
        out << this->position_[i] << ", ";
    out << "\b\b)" << std::endl;

    out << "Constraint violation:\t" << constraint_violation_ << std::endl;
}

template <std::size_t dim>
void Bee<dim>::compute_fitness_value()
{
        if (cost_value_ >= 0)
            fitness_value_ = 1.0 / (1.0 + cost_value_);
        else
            fitness_value_ = 1.0 + std::abs(cost_value_);
}

/*int main(){
    Bee<2> a;
    std::vector<Bee<2>> colony;
    colony.push_back(a);
    a.initialize();
}*/
