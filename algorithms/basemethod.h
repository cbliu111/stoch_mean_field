//
// Created by cbl on 12/18/22.
//

#ifndef SIM_BASEMETHOD_H
#define SIM_BASEMETHOD_H

#include "model.h"

enum REACTANT_TYPE {

    HOR_3_COFF_3 = 1,
    HOR_3_COFF_2 = 2,
    HOR_3_COFF_1 = 3,

    HOR_2_COFF_2 = 4,
    HOR_2_COFF_1 = 5,

    HOR_1_COFF_1 = 6,

    NON_DEFINE = 7,
};

class BaseSimMethod {
public:
    Model model;
    RandEngine rand;
    size_t M;
    size_t N;
    double* populations;
    double weight;
    double time;
    double step;
    double end_time;
    int steps;
    double record_start_time;
    double time_check_point;
    double record_time_interval;
    int repeats;
    std::string file;
    long long rand_seed;
    explicit BaseSimMethod(Model& m) {
        model = m;
        M = model.reactions.size();
        N = model.species.size();
        populations = new double [N];
        for (int i = 0; i < N; ++i) {
            populations[i] = model.populations[i];
        }
        time = 0;
        // set initial weight equals to 1
        weight = 0;
        step = 0;
        record_start_time = 0;
        time_check_point = record_start_time;
        record_time_interval = std::numeric_limits<double>::infinity();
        repeats = 0;
        end_time = model.end_time;
        steps = model.steps;
        rand_seed = 0;
    }
    void set_seed(long long s) {
        // set the seed of first run
        // if multiple trajectories are needed, seed will be picked randomly
        rand_seed = s;
    }
    void set_end_time(double t) {
        end_time = t;
    }
    void set_total_steps(int s) {
        steps = s;
    }
    void set_number_data_points(int n) {
        if (n < 0) {
            return;
        }
        switch (n) {
            case 1:
                // if rst is equal to end_time, only record the final state
                record_start_time = end_time;
                break;
            case 0:
                // if rst is inf, no recording
                record_start_time = std::numeric_limits<double>::infinity();
                break;
            default:
                // record 0, ..., final state, there are (n-2) points between 0 and final time
                record_start_time = 0;
                record_time_interval = end_time / static_cast<double>(n - 1);
        }
    }
    void set_number_trajectory(int n) {
        // if n is 1, single trajectory simulation
        repeats = n;
    }
    void set_output_file(const std::string& filename) {
        file = filename;
    }
    /*
    set time = 0
    set step = 0
    clear trajectory
    init population from original model
    set weight = 0
    set time check point back to record start time
    */
    void reinitialize() {
        time = 0;
        step = 0;
        model.trajectory.clear();
        // assign populations from the original model to set simulation at the same initial state
        for (int i = 0; i < N; ++i) {
            populations[i] = model.populations[i];
        }
        time_check_point = record_start_time;
        // reset weight to 1
        weight = 1;
    }
    void set_random_seed_by_clock() {
        rand_seed = clock();
        rand.seed(rand_seed);
    }
    void record_trajectory() {
        // if the trajectory does not meet the check point time, its value wouldn't be saved
        if (time >= time_check_point) {
            model.record_trajectory(populations, time, weight);
            time_check_point += record_time_interval;
        }
    }
    void output_trajectory_to_file(int repeat) {
        // if there are no data in trajectory, continue to next repeat
        if (model.trajectory.empty()) {
            return;
        }
        // append all results in the same file, if multi-trajectories are simulated
        if (repeat == 0) model.output_trajectory_to_file(file, false);
        else model.output_trajectory_to_file(file, true);
    }
    virtual void sim() {}
    void run() {
        boost::timer::auto_cpu_timer timer;
        for (int r = 0; r < repeats; ++r) {
            reinitialize();
            sim();
            set_random_seed_by_clock();
            if (record_start_time > end_time) continue;
            output_trajectory_to_file(r);
        }
    }
    ~BaseSimMethod() {
        delete [] populations;
    }
};

/*
 * Base method for aleap, ileap and hleap
 * Implement compute_delta, compute_a0, set_probability_bound
 * generate_event_sequence, next_event, and check_out_range
 */
class BaseLeapMethod : public BaseSimMethod {
public:
    double event_count; // total event number
    int num_dm; // number of direct method steps
    double ua0; // sum of propensity upper bounds
    double a0; // total propensity
    double bw; // bound weight
    double epsilon; // error factor for tau leap
    double tau; // time interval for time line division
    double T; // time axis check point
    double default_leap_size; // preset leap length
    double leap_threshold; // criterion for leaping
    double* pop_lm; // leap mean population
    double* prop_mf; // mean field propensities
    double* prop; // propensities of current state
    double* delta; // species specific state bound factors
    std::vector<REACTANT_TYPE> reactant_types; // reactant types for compute tau
    std::vector<int> events; // list of event_vector
    std::vector<double> dn; // population difference
    std::vector<double> event_num; // number of events for channels
    std::vector<double> pop_copy; // backup population
    std::vector<double> max_reactant_change; // maximum population change given epsilon
    std::vector<double> species_net_change; // net species change given current state
    std::vector<double> species_net_change_squared; // net species change given current state
    explicit BaseLeapMethod(Model& m, 
        double bound_weight=-1, 
        double error_factor=-1, 
        double time_interval=-1) : BaseSimMethod(m) {
        // constructor for weight bound method
        event_count = 0;
        num_dm = 50;
        ua0 = 0;
        a0 = 0;
        bw = bound_weight;
        epsilon = error_factor;
        tau = time_interval;
        T = tau;
        event_count = 0;
        default_leap_size = 20;
        leap_threshold = 10;
        pop_lm = new double [N];
        prop_mf = new double [M];
        prop = new double [M];
        delta = new double [N];
        reactant_types.assign(N, NON_DEFINE);
        events.assign(M, 0);
        dn.assign(N, 0);
        event_num.assign(M, 0);
        pop_copy.assign(N, 0);
        max_reactant_change.assign(N, 0);
        species_net_change.assign(N, 0);
        species_net_change_squared.assign(N, 0);
        for (int i = 0; i < M; ++i) {
            compute_reaction_HOR(model.reactions[i]);
        }
        if (bw > 0) {
            // bound weight method, compute epsilon
            compute_epsilon();
        }
    }
    [[nodiscard]] double func(double x) const {
        return x * std::exp(1-x) - bw;
    }
    void compute_epsilon() {
        double a = 0;
        double b = 1;
        while ((b - a) >= 1e-3) {
            epsilon = (a + b) / 2;
            if (std::abs(func(epsilon)) < 1e-5) {
                break;
            }
            else if (func(epsilon) * func(a) < 0) {
                b = epsilon;
            }
            else {
                a = epsilon;
            }
        }
        epsilon = 1 - epsilon;
    }
    void compute_reaction_HOR(Reaction& r) {
        double hor = 0;
        for (Reactant& i : r.reactants) {
            hor += i.stoichiometric;
        }
        for (Reactant& i : r.reactants) {
            double coff = i.stoichiometric;
            if (hor == 3) {
                if (coff == 3) {
                    REACTANT_TYPE x = HOR_3_COFF_3;
                    REACTANT_TYPE y = reactant_types[i.index];
                    reactant_types[i.index] = (y < x) ? y : x;
                }
                else if (coff == 2) {
                    REACTANT_TYPE x = HOR_3_COFF_2;
                    REACTANT_TYPE y = reactant_types[i.index];
                    reactant_types[i.index] = (y < x) ? y : x;
                }
                else if (coff == 1) {
                    REACTANT_TYPE x = HOR_3_COFF_1;
                    REACTANT_TYPE y = reactant_types[i.index];
                    reactant_types[i.index] = (y < x) ? y : x;
                }
                else {
                    reactant_types[i.index] = NON_DEFINE;
                }
            }
            else if (hor == 2) {
                if (coff == 2) {
                    REACTANT_TYPE x = HOR_2_COFF_2;
                    REACTANT_TYPE y = reactant_types[i.index];
                    reactant_types[i.index] = (y < x) ? y : x;
                }
                else if (coff == 1) {
                    REACTANT_TYPE x = HOR_2_COFF_1;
                    REACTANT_TYPE y = reactant_types[i.index];
                    reactant_types[i.index] = (y < x) ? y : x;
                }
                else {
                    reactant_types[i.index] = NON_DEFINE;
                }
            }
            else if (hor == 1) {
                if (coff == 1) {
                    REACTANT_TYPE x = HOR_1_COFF_1;
                    REACTANT_TYPE y = reactant_types[i.index];
                    reactant_types[i.index] = (y < x) ? y : x;
                }
                else {
                    reactant_types[i.index] = NON_DEFINE;
                }
            }
            else {
                reactant_types[i.index] = NON_DEFINE;
            }
        }
    }
    void compute_max_change_by_reactant(Reaction& r) {
        for (Reactant& i : r.reactants) {
            int species = i.index;
            double pop = populations[species];
            double g;
            switch (reactant_types[i.index]) {
            case HOR_3_COFF_3:
                if (pop > 2) {
                    g = 3.0 + 1.0 / (pop -1) + 2.0 / (pop - 2);
                }
                else {
                    g = std::numeric_limits<double>::max();
                }
                break;
            case HOR_3_COFF_2:
                if (pop > 1) {
                    g = 3.0 + 3.0 / (2*(pop -1));
                }
                else {
                    g = std::numeric_limits<double>::max();
                }
                break;
            case HOR_3_COFF_1:
                g = 3.0;
                break;
            case HOR_2_COFF_2:
                if (pop > 1) {
                    g = 2.0 + 1.0 / (pop - 1);
                }
                else {
                    g = std::numeric_limits<double>::max();
                }
                break;
            case HOR_2_COFF_1:
                g = 2.0;
                break;
            default:
                g = 1.0;
                break;
            }
            double max_change = (epsilon / g) * pop;
            max_change = (max_change > 1) ? max_change : 1;
            double pre_change = max_reactant_change[i.index];
            max_reactant_change[i.index] = (pre_change > max_change) ? pre_change : max_change;
        }
    }
    void compute_propensities() {
        a0 = 0;
        for (int i = 0; i < M; ++i) {
            double p = calculate_propensity(populations, model.reactions[i]);
            prop[i] = p;
            a0 += p;
        }
    }
    void compute_net_change() {
        species_net_change.assign(N, 0);
        species_net_change_squared.assign(N, 0);
        for (int i = 0; i < M; ++i) {
            for (Reactant& r : model.reactions[i].reactants) {
                double net_change = r.stoichiometric * prop_mf[i];
                double net_change_squared = r.stoichiometric * r.stoichiometric * prop_mf[i];
                species_net_change[r.index] += net_change;
                species_net_change_squared[r.index] += net_change_squared;
            }
        }
    }
    void compute_tau() {
        tau = std::numeric_limits<double>::max();
        compute_mean_field_propensities();
        if (ua0 == 0) {
            tau = 0;
            return;
        }
        compute_net_change();
        for (int i = 0; i < M; ++i) {
            compute_max_change_by_reactant(model.reactions[i]);
        }
        for (int i = 0; i < N; ++i) {
            double c = max_reactant_change[i];
            double nc = species_net_change[i];
            double nc2 = species_net_change_squared[i];
            double possible_time = c / (2 * std::abs(nc));
            double possible_time_squared = c * c / (4 * std::abs(nc2));
            double min_time = (possible_time > possible_time_squared) ? possible_time_squared : possible_time;
            if (tau > min_time) {
                tau = min_time;
            }
        }
        if (tau == std::numeric_limits<double>::max()) {
            tau = default_leap_size / ua0;
        }
    }
    void backup_state() {
        for (int i = 0; i < N; ++i) {
            pop_copy[i] = populations[i];
        }
    }
    void restore_state() {
        for (int i = 0; i < N; ++i) {
            populations[i] = pop_copy[i];
        }
    }
    bool check_state_bound() {
        bool flag = true; // true is inside state boundaries
        for (int i = 0; i < N; ++i) {
            if (populations[i] < 0) {
                flag = false;
                break;
            }
        }
        return flag;
    }
    void do_ssa(int n) {
        for (int s = 0; s < n; ++s) {
            double sum = 0;
            for (int i = 0; i < M; ++i) {
                double p = calculate_propensity(populations, model.reactions[i]);
                sum += p;
                prop[i] = p;
            }
            if (sum == 0) break;
            time += std::log(1.0 / rand.uniform01()) / sum;
            double r = rand.uniform01() * sum;
            sum = 0;
            for (int i = 0; i < M; ++i) {
                sum += prop[i];
                if (r < sum) {
                    model.update_populations(populations, i);
                    break;
                }
            }
            step += 1;
            record_trajectory();
            if (time > end_time) break;
        }
    } 
    /*
    calculate bound_factors for each reaction channel based on bounded acceptance probability
     */
    void compute_delta() {
        auto* bound_factors = new double [M]; // bound factor for each reaction
        for (int i = 0; i < M; ++i) {
            // calculate sum of the stoichiometric number
            // and the maximum stoichiometric number
            double s = 0;
            double max_s = 0;
            for (Reactant& r : model.reactions[i].reactants) {
                s += r.stoichiometric;
                max_s = (r.stoichiometric > max_s) ? r.stoichiometric : max_s;
            }
            double pbs = std::pow(bw, s);
            bound_factors[i] = (1 - pbs) / ((1 - pbs) * max_s + 2 * pbs);
        }
        for (int i = 0; i < N; ++i) {
            double d = std::numeric_limits<double>::max();
            for (unsigned int j : model.species[i].linked_reactions) {
                d = (d < bound_factors[j]) ? d : bound_factors[j];
            }
            delta[i] = d;
        }
    }
    void compute_mean_field_propensities() {
        ua0 = 0;
        for (int i = 0; i < M; ++i) {
            double p = calculate_propensity(populations, model.reactions[i]);
            prop_mf[i] = p;
            ua0 += p;
        }
    }
    /*
     * events are generated by
     * Poisson: for single leap approximation
     * Gaussian: for Langevin approximation
     * Deterministic: for ODE approximation
     */
    void generate_events_poisson() {
        event_num.clear();
        event_count = 0;
        for (int i = 0; i < M; ++i) {
            double n = rand.poisson(prop_mf[i] * tau);
            event_num.push_back(n);
            event_count += n;
        }
    }
    void generate_events_langevin() {
        event_num.clear();
        event_count = 0;
        for (int i = 0; i < M; ++i) {
            double lambda = prop_mf[i] * tau;
            double n = rand.normal(lambda, std::sqrt(lambda));
            event_num.push_back(n);
            event_count += n;
        }
    }
    void generate_events_ode() {
        event_num.clear();
        event_count = 0;
        events.clear();
        for (int i = 0; i < M; ++i) {
            double n = prop_mf[i] * tau;
            event_num.push_back(n);
            event_count += n;
        }
    }
    /*
    turn event number into a vector of event index
    */
    void generate_event_sequence() {
        events.clear();
        for (int i = 0; i < M; ++i) {
            double n = event_num[i];
            for (int j = 0; j < n; ++j) {
                events.push_back(i);
            }
        }
    }
    /*
    choose next event by Fisher-Yates algorithm
    */
    int next_event(size_t i) {
        int e = rand.discrete_uniform(0, static_cast<int>(i)-1); // include i - 1
        // unsigned int e = rand.mt19937.operator()() % i;
        int id = events[e];
        events[e] = events[i-1];
        return id;
    }
    /*
    population differences are calculated for a leap
    if some population has negative value, then return true
    */
    bool check_out_range() {
        bool flag = false;
        dn.assign(N, 0);
        for (int i = 0; i < M; ++i) {
            for (JumpVector& v : model.reactions[i].jump_vectors) {
                if (!model.species[v.index].is_constant) {
                    dn[v.index] += event_num[i] * v.step_size;
                }
            }
        }
        for (int i = 0; i < N; ++i) {
            if (populations[i] + dn[i] < 0) {
                flag = true;
                break;
            }
        }
        return flag;
    }
    /*
    If out range, return true
    */
    bool check_reaction(Reaction& r) {
        bool flag = false;
        for (JumpVector& v : r.jump_vectors) {
            if (!model.species[v.index].is_constant) {
                if (populations[v.index] + v.step_size < 0) {
                    flag = true;
                    break;
                }
            }
        }
        return flag;
    }
    double compute_leap_weight() {
        double w = 1;
        double a_lm = 0;
        // leap to leap mean state
        // record population before leap
        for (int i = 0; i < N; ++i) {
            // pop_lm[i] = populations[i] + dn[i] / 2;
            pop_lm[i] = (pop_copy[i] + populations[i]) / 2;
        }
        // compute leap weight
        for (int i = 0; i < M; ++i) {
            double p = calculate_propensity(pop_lm, model.reactions[i]);
            a_lm += p;
            w *= std::pow(p / prop_mf[i], event_num[i]);
        }
        w *= std::exp(tau * (ua0 - a_lm));
        return w;
    }
    /*
    leap mean weight is computed by using the leap mean population
    and leap mean population is the mean state between start and end state
    mean weight is returned
    */
    double compute_leap_mean_weight() {
        double w = 1;
        double a_lm = 0;
        // leap to leap mean state
        // record population before leap
        for (int i = 0; i < N; ++i) {
            pop_lm[i] = populations[i] + dn[i];
        }
        // compute leap weight
        for (int i = 0; i < M; ++i) {
            double p = calculate_propensity(pop_lm, model.reactions[i]);
            a_lm += p;
            w *= std::pow(p / prop_mf[i], event_num[i]);
        }
        w *= std::exp(tau * (ua0 - a_lm));
        // return leap mean weight
        w = std::pow(w, 1./event_count);
        return w;
    }
    /*
    step weight is computed using the current state and the mean-field state
    */
    double compute_step_weight(unsigned int id, double dt) {
        double a = 0;
        for (int i = 0; i < M; ++i) {
            double p = calculate_propensity(populations, model.reactions[i]);
            a += p;
        }
        double p = calculate_propensity(populations, model.reactions[id]);
        double w = p / prop_mf[id] * std::exp(-(a-ua0)*dt);
        return w;
    }
    ~BaseLeapMethod() {
        delete [] pop_lm;
        delete [] prop_mf;
        delete [] prop;
        delete [] delta;
    }
};

#endif //SIM_BASEMETHOD_H
