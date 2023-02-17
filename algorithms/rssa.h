//
// Created by cbl on 6/21/22.
//

#ifndef SIM_RSSA_H
#define SIM_RSSA_H

#include "timer.h"
#include "basemethod.h"

class RSSA : public BaseSimMethod {
public:
    // propensity low bounds, upper bounds and molecular number low and upper bounds
    double* prop_lb;
    double* prop_ub;
    double* m_lb;
    double* m_ub;
    double upper_factor;
    double lower_factor;
    explicit RSSA(Model& m, double delta) : BaseSimMethod(m) {
        prop_lb = new double [M];
        prop_ub = new double [M];
        m_lb = new double [N];
        m_ub = new double [N];
        upper_factor = 1+delta;
        lower_factor = 1-delta;
        init_bounds();
    }
    void init_bounds() {
        for (int i = 0; i < N; ++i) {
            m_lb[i] = std::ceil(lower_factor * populations[i]);
            m_ub[i] = std::floor(upper_factor * populations[i]);
        }
        for (int i = 0; i < M; ++i) {
            double lp = calculate_propensity(m_lb, model.reactions[i]);
            double up = calculate_propensity(m_ub, model.reactions[i]);
            prop_lb[i] = lp;
            prop_ub[i] = up;
        }
    }
    void sim() override {
        init_bounds();
        while ((time < end_time) || (step < steps)) {
            double ua0 = 0;
            for (int i = 0; i < M; ++i) {
                ua0 += prop_ub[i];
            }
            if (ua0 == 0) break;
            bool out_range = false;
            int idx;
            while ((!out_range) && (time < end_time)) {
                double u = 1;
                bool accept = false;
                double sum;
                while (!accept) {
                    double r1 = rand.uniform01() * ua0;
                    double r2 = rand.uniform01();
                    double r3 = rand.uniform01();
                    idx = 0;
                    sum = 0;
                    for (; idx < M; ++idx) {
                        sum += prop_ub[idx];
                        if (r1 <= sum) break;
                    }
                    if (r2 <= prop_lb[idx] / prop_ub[idx]) accept = true;
                    else {
                        double p = calculate_propensity(populations, model.reactions[idx]);
                        if (r2 <= p / prop_ub[idx]) accept = true;
                    }
                    u *= r3;
                }
                step += 1;
                time += std::log(1.0 / u) / ua0;
                model.update_populations(populations, idx);
                record_trajectory();
                
                for (int i : model.reactions[idx].species) {
                    if ((populations[i] > m_ub[i]) || (populations[i] < m_lb[i])) {
                        out_range = true;
                        m_lb[i] = std::ceil(lower_factor * populations[i]);
                        m_ub[i] = std::floor(upper_factor * populations[i]);
                    }
                }
            }
            for (int i : model.reactions[idx].species) {
                for (int j : model.species[i].linked_reactions) {
                    double lp = calculate_propensity(m_lb, model.reactions[j]);
                    double up = calculate_propensity(m_ub, model.reactions[j]);
                    prop_lb[j] = lp;
                    prop_ub[j] = up;
                }
            }
            if (time > end_time) break;
        }
    }
    ~RSSA() {
        delete [] prop_lb;
        delete [] prop_ub;
        delete [] m_lb;
        delete [] m_ub;
    }
};

void test_rssa(const std::string& model, double end_time, double delta) {
    const std::string x = "../models/" + model + ".txt";
    const char* model_file = x.c_str();
    const std::string traj_file = "../data/" + model + "-rssa-single_trajectory_out.csv";
    const std::string dist_file = "../data/" + model + "-rssa-population_dist.csv";
    Model model_t(model_file);
    auto engine = RSSA(model_t, delta);
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

void profile_rssa(const std::string& model, double end_time, double delta, int runs) {
    const std::string x = "../models/" + model + ".txt";
    const char* model_file = x.c_str();
    const std::string traj_file = "../data/" + model + "-rssa-single_trajectory_out.csv";
    const std::string dist_file = "../data/" + model + "-rssa-population_dist.csv";
    Model model_t(model_file);
    auto engine = RSSA(model_t, delta);
    double et = end_time; // end time

    // single trajectory
    engine.set_end_time(et);
    engine.set_output_file(traj_file);
    engine.set_number_data_points(1000);
    engine.set_number_trajectory(1);
    engine.set_seed(0);
    engine.run();
    // bench mark
    if (runs <= 0) return;
    engine.set_end_time(et);
    engine.set_output_file("");
    engine.set_number_data_points(0);
    engine.set_number_trajectory(runs);
    engine.set_seed(0);
    engine.run();
}

#endif //SIM_RSSA_H
