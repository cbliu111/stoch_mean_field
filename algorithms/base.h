//
// Created by cbl on 4/23/22.
//

#ifndef SIM_BASE_H
#define SIM_BASE_H

#include <cstdio>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_set>
#include "random_engine.h"

enum ReactionType {
    SOURCE,
    UNIMOLECULAR,
    BIMOLECULAR,
    COMPLEX
};

struct Reactant {
    int index = 0;
    double stoichiometric = 0;
};

struct Product {
    int index = 0;
    double stochiometric = 0;
};

struct Species {
    std::string name;
    bool is_constant;
    // linked reactions are reactions that has this species as one of its reactant
    std::vector<int> linked_reactions;
};

struct JumpVector {
    int index = 0;
    double step_size = 0;
};

struct State {
    double time = 0;
    std::vector<double> populations;
    double weight = 0;
};

class Reaction {
public:
    ReactionType reaction_type = SOURCE;
    int reaction_index = 0;
    double rate = 0;
    std::vector<int> species;
    std::vector<Reactant> reactants;
    std::vector<Product> products;
    std::vector<JumpVector> jump_vectors;
    std::vector<int> linked_reactions;
    bool operator==(const Reaction& other) const {
        return reaction_index == other.reaction_index;
    }
};

static inline double calculate_combinatorial(const double n, const double k) {
    /*
     * Boundary condition is controlled by the value of reaction val.
     * If number of molecules is not enough to fire a particular reaction,
     * the val is set to 0.
     * A 0 val will have no weight when using DM based algorithms,
     * therefore will not be selected.
     * A 0 val will also be assign infinite large reaction time,
     * therefore will also not be selected by FRM based algorithms.
     * */
    if (n < k) return 0;
    switch (static_cast<int>(k)) {
        case 1: return n;
        case 2: return n * (n-1) / 2;
        default: {
            unsigned int o = 1;
            for (int i = 1; i <= k; ++i) {
                o *= n - i + 1;
                o /= i;
            }
            return o;
        }
    }
}

static double calculate_propensity(const double* populations, Reaction& r) {
    switch (r.reaction_type) {
        case SOURCE: return r.rate;
        case UNIMOLECULAR: {
            double n = populations[r.reactants[0].index];
            double k = r.reactants[0].stoichiometric;
            return calculate_combinatorial(n, k) * r.rate;
        }
        case BIMOLECULAR: {
            double o = 1.0;
            double n0 = populations[r.reactants[0].index];
            double n1 = populations[r.reactants[1].index];
            double k0 = r.reactants[0].stoichiometric;
            double k1 = r.reactants[1].stoichiometric;
            o *= calculate_combinatorial(n0, k0);
            o *= calculate_combinatorial(n1, k1);
            return o * r.rate;
        }
        case COMPLEX: {
            double o = 1.0;
            for (auto &i : r.reactants) {
                double n = populations[i.index];
                double k = i.stoichiometric;
                o *= calculate_combinatorial(n, k);
            }
            return o * r.rate;
        }
        default: {
            return 0;
        }
    }
}

static double calculate_propensity_with_delta(const double* populations, Reaction& r, double delta) {
    switch (r.reaction_type) {
        case SOURCE: return r.rate;
        case UNIMOLECULAR: {
            double n = std::floor(populations[r.reactants[0].index] * (1 + delta));
            double k = r.reactants[0].stoichiometric;
            return calculate_combinatorial(n, k) * r.rate;
        }
        case BIMOLECULAR: {
            double o = 1.0;
            double n0 = std::floor(populations[r.reactants[0].index] * (1 + delta));
            double n1 = std::floor(populations[r.reactants[1].index] * (1 + delta));
            double k0 = r.reactants[0].stoichiometric;
            double k1 = r.reactants[1].stoichiometric;
            o *= calculate_combinatorial(n0, k0);
            o *= calculate_combinatorial(n1, k1);
            return o * r.rate;
        }
        case COMPLEX: {
            double o = 1.0;
            for (auto &i : r.reactants) {
                double n = std::floor(populations[i.index] * (1 + delta));
                double k = i.stoichiometric;
                o *= calculate_combinatorial(n, k);
            }
            return o * r.rate;
        }
        default: {
            return 0;
        }
    }
}

static void calculate_propensity_bound_factors(double delta, Reaction& r, double& upper_bound, double& lower_bound) {
    double n = std::ceil(1/delta);
    double mlb = std::ceil((1-delta)*n);
    double mub = std::floor((1+delta)*n);
    upper_bound = 1;
    lower_bound = 1;
    for (Reactant& react : r.reactants) {
        double k = react.stoichiometric;
        for (int i = 1; i < k+1; ++i) {
            lower_bound *= (mlb - (i - 1));
            lower_bound /= (n - (i - 1));
            upper_bound *= (mub - (i - 1));
            upper_bound /= (n - (i - 1));
        }
    }
}

/*
static void calculate_propensity_bound_factors(double bound_factors, Reaction& r, double& upper_bound, double& lower_bound) {
    unsigned int n = std::ceil(1/bound_factors);
    unsigned int mlb = std::ceil((1-bound_factors)*n);
    unsigned int mub = std::floor((1+bound_factors)*n);
    upper_bound = 1;
    lower_bound = 1;
    for (Reactant& i : r.reactants) {
        unsigned int k = i.stoichiometric;
        lower_bound *= static_cast<double>(calculate_combinatorial(mlb, k));
        lower_bound /= static_cast<double>(calculate_combinatorial(n, k));
        upper_bound *= static_cast<double>(calculate_combinatorial(mub, k));
        upper_bound /= static_cast<double>(calculate_combinatorial(n, k));
    }
}
*/

static double calculate_reaction_time(double rand, double propensity, double sys_time) {
    if (propensity == 0) {
        return std::numeric_limits<double>::infinity();
    }
    else {
        return std::log(1. / rand) / propensity + sys_time;
    }
}

static double update_reaction_time(RandEngine& rand, double propensity, double old_propensity, double old_time, double sys_time) {
    /*
     * if previous val is 0, then recalculate reaction time from val
     * if current val is 0, then assign end time
     * otherwise, update reaction time based on previous, current propensities and current time of the system
     * */
    if (propensity == 0) {
        return std::numeric_limits<double>::infinity();
    }
    if (old_propensity == 0) {
        return std::log(1. / rand.uniform01()) / propensity + sys_time;
    }
    return old_propensity / propensity * (old_time - sys_time) + sys_time;
}

template<class T>
class Matrix {
private:
    T* array; 
    int w;
    int h;
public: 
    Matrix() = default;
    Matrix(int width, int height) {
        array = new T [w * h];
        w = width;
        h = height;
    } 
    ~Matrix() {
        delete[] array;
    } 
    T at(int x, int y) const { 
        return array[y + w * x]; 
    } 
    void set(int x, int y, T value) {
        array[y + w * x] = value;
    }
    void assign(T value) {
        for (int i = 0; i < h; ++i) {
            for (size_t j = 0; j < w; ++j) {
                array[j + w * i] = value;
            }
        }
    }
    void assign_column(int x, T value) {
        for (int i = 0; i < h; ++i) {
            array[x + w * i] = value;
        }
    }
};

template<class T>
class Vector {
private:
    T* array;
    size_t s;
public:
    Vector() = default;
    Vector(size_t size) : s(size), array(new T [size]) {}
    ~Vector() {
        delete [] array;
    }
    T at(size_t x) const {
        return array[x];
    }
    void set(size_t x, T value) {
        array[x] = value;
    }
    void assign(T value) {
        for (size_t i = 0; i < s; ++i) {
            array[i] = value;
        }
    }
};


#endif //SIM_BASE_H
