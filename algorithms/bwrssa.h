//
// Created by cbl on 1/23/23.
//

#ifndef SIM_BWRSSA_H
#define SIM_BWRSSA_H

/*
if using only the bounded acceptance probability to select the step,
there will be a bias for the mean value
should using the complete weight to select the step
bound acceptance probability rssa algorithm
*/

#include "timer.h"
#include "basemethod.h"

class BWRSSA : public BaseSimMethod {
public:
    double* prop;
    double* m_lb;
    double* m_ub;
    double alpha;
    double a0;
    std::vector<double> delta;
    explicit BWRSSA(Model& m, double bound_weight)
            : BaseSimMethod(m) {
        alpha = bound_weight;
        a0 = 0;
        prop = new double [M];
        m_lb = new double [N];
        m_ub = new double [N];
        double t = std::numeric_limits<double>::max();
        delta.assign(N, t);
    }
    void compute_delta() {
        for (int i = 0; i < M; ++i) {
            size_t num_reactant = model.reactions[i].reactants.size();
            for (Reactant& j : model.reactions[i].reactants) {
                int idx = j.index;
                int coff = j.stoichiometric;
                if (num_reactant == 1) {
                    double v = (1-alpha) / (1+alpha);
                    delta[idx] = (delta[idx] < v) ? delta[idx] : v;
                }
                else if (num_reactant == 2) {
                    if (coff == 1) {
                        double v = (1 - std::sqrt(alpha)) / (1 + std::sqrt(alpha));
                        delta[idx] = (delta[idx] < v) ? delta[idx] : v;
                    }
                    else if (coff == 2) {
                        double v = (1 - std::sqrt(alpha)) / (1 + std::sqrt(alpha)) * (1 - 1 / populations[idx]);
                        delta[idx] = (delta[idx] < v) ? delta[idx] : v;
                    }
                }
            }
        }
        for (int i = 0; i < N; ++i) {
            // a more restrictive region should be used to reduce sample variance
            delta[i] = (1 - std::sqrt(alpha)) / (1 + std::sqrt(alpha));
        }
    }
    void compute_bounds() {
        for (int i = 0; i < N; ++i) {
            m_lb[i] = std::ceil((1 - delta[i]) * populations[i]);
            m_ub[i] = std::floor((1 + delta[i]) * populations[i]);
        }
    }
    bool check_range() {
        bool flag = true;
        for (int i = 0; i < N; ++i) {
            if (populations[i] >= m_ub[i] || populations[i] <= m_lb[i] || populations[i] <= 0) {
                flag = false;
                break;
            }
        }
        return flag;
    }
    void compute_propensities() {
        a0 = 0;
        for (int i = 0; i < M; ++i) {
            double p = calculate_propensity(populations, model.reactions[i]);
            prop[i] = p;
            a0 += p;
        }
    }
    void sim() override {
        while ((time < end_time) || (step < steps)) {
            compute_propensities();
            bool in_range = true;
            while (in_range) {
                double r = rand.uniform01() * a0;
                double sum = 0;
                int idx = 0;
                for (; idx < M; ++idx) {
                    sum += prop[idx];
                    if (r <= sum) break;
                }
                
                model.update_populations(populations, idx);
                time += std::log(1.0 / rand.uniform01()) / a0;
                step += 1;
                in_range = check_range();
            }
            record_trajectory();
        }
    }
    ~BWRSSA() = default;
};

void test_bwrssa(const std::string& model, double end_time, double accept_prob_bound) {
    const std::string x = "../models/" + model + ".txt";
    const char* model_file = x.c_str();
    const std::string traj_file = "../data/" + model + "-bwrssa-single_trajectory_out.csv";
    const std::string dist_file = "../data/" + model + "-bwrssa-population_dist.csv";
    Model model_t(model_file);
    auto engine = BWRSSA(model_t, accept_prob_bound);
    double et = end_time; // end time

    // single trajectory
    engine.set_end_time(et);
    engine.set_output_file(traj_file);
    engine.set_number_data_points(100);
    engine.set_number_trajectory(1);
    engine.set_seed(0);
    engine.run();
    // bench mark
    engine.set_end_time(et);
    engine.set_output_file("");
    engine.set_number_data_points(0);
    engine.set_number_trajectory(100);
    engine.set_seed(0);
    engine.run();
    // population distribution at end time
    engine.set_end_time(et);
    engine.set_output_file(dist_file);
    engine.set_number_data_points(1);
    engine.set_number_trajectory(100000);
    engine.set_seed(0);
    engine.run();
}
#endif //SIM_BWRSSA_H
