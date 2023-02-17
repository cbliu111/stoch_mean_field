//
// Created by cbl on 4/23/22.
//

#ifndef SIM_MODEL_H
#define SIM_MODEL_H

#include "base.h"
#include <fstream>


class Model {
public:
    int M = 0;
    int N = 0;
    int steps = 0;
    double end_time = 0;
    RandEngine rand;
    std::vector<Species> species;
    std::vector<double> populations;
    std::vector<Reaction> reactions; // get number of reactions by reactions.pos()
    std::vector<State> trajectory;
    std::vector<std::vector<double>> sm; // stoichiometric matrix

    explicit Model() = default;
    explicit Model(const char *filename);
    void add_species(const std::string &line);
    void add_reaction(const std::string &line);
    void add_condition(const std::string &line);
    int find_specie_id(const std::string &name);
    static void calculate_reaction_properties(Reaction &reaction);
    void build_reaction_reaction_dependency_graph();
    void build_species_reaction_dependency_graph();
    static bool is_linked(const Reaction &m, const Reaction &n);
    bool is_linked(const Species &m, const Reaction &n);
    void update_populations(double* populations, int index);
    void update_populations_leap(double* input_populations, int index, double n);
    void update_populations_sm(double* populations, int index);
    void record_trajectory(double* populations, double time, double weight);
    void output_trajectory_to_file(const std::string& filename, bool app=false);
    void print() const;
    void print_latex() const;
    void generate_cyclic_chain_model(int num_species);
    void generate_colloidal_aggregation_model(int num_species);
    void generate_random_connect_model(int num_species, int couple_degree);
    void clear();
    bool is_empty() const;
    ~Model() = default;
};

#endif //SIM_MODEL_H
