//
// Created by cbl on 6/14/22.
//

#ifndef SIM_RNRM_H
#define SIM_RNRM_H

#include "timer.h"
#include "basemethod.h"

class RNRM : public BaseSimMethod {
public:
    double* m_lb;
    double* m_ub;
    double* ubf;
    double* lbf;
    double* prop_ub;
    double* old_time;
    double delta;
    explicit RNRM(Model& m, double d) : BaseSimMethod(m) {
        m_lb = new double [N];
        m_ub = new double [N];
        ubf = new double[M];
        lbf = new double[M];
        prop_ub = new double[M];
        old_time = new double[M];
        delta = d;
        init_bounds_times();
    }
    void init_bounds_times() {
        for (int i = 0; i < N; ++i) {
            unsigned int x = model.populations[i];
            populations[i] = x;
            m_lb[i] = std::ceil((1-delta) * x);
            m_ub[i] = std::floor((1+delta) * x);
        }
        for (int i = 0; i < M; ++i) {
            calculate_propensity_bound_factors(delta, model.reactions[i], ubf[i], lbf[i]);
            double up = ubf[i] * calculate_propensity(populations, model.reactions[i]);
            prop_ub[i] = up;
            double r0 = rand.uniform01();
            double t = calculate_reaction_time(r0, up, time);
            old_time[i] = t;
        }
    }
    void sim() override {
        init_bounds_times();
        while ((time < end_time) || (step < steps)) {
            unsigned int idx;
            bool out_range = false;
            while ((!out_range) && (time < end_time)) {
                bool accept = false;
                while (!accept) {
                    double r1 = rand.uniform01();
                    double r2 = rand.uniform01();
                    idx = 0;
                    for (int i = 0; i < M; ++i) {
                        idx = (old_time[idx] > old_time[i]) ? i : idx;
                    }
                    time = old_time[idx];
                    double up = prop_ub[idx];
                    if (r1 <= lbf[idx]/ubf[idx]) accept = true;
                    else {
                        double p = calculate_propensity(populations, model.reactions[idx]);
                        if (r1 <= p / up) accept = true;
                    }
                    double t = calculate_reaction_time(r2, up, time);
                    old_time[idx] = t;
                }
                step += 1;
                model.update_populations(populations, idx);
                record_trajectory();
                for (int i : model.reactions[idx].species) {
                    double x = populations[i];
                    if ((x > m_ub[i]) || (x < m_lb[i])) {
                        out_range = true;
                        m_lb[i] = std::ceil((1-delta) * x);
                        m_ub[i] = std::floor((1+delta) * x);
                        for (int j : model.species[i].linked_reactions) {
                            double up = ubf[j] * calculate_propensity(populations, model.reactions[j]);
                            double t;
                            if (up == 0) {
                                t = std::numeric_limits<double>::infinity();
                            }
                            else if (prop_ub[j] == 0) {
                                double r0 = rand.uniform01();
                                t = std::log(1/r0) / up + time;
                            }
                            else {
                                t = prop_ub[j] / up * (old_time[j] - time) + time;
                            }
                            prop_ub[j] = up;
                            old_time[j] = t;
                        }
                    }
                }
            }
        }
    }
    ~RNRM() {
        delete [] m_lb;
        delete [] m_ub;
        delete [] ubf;
        delete [] lbf;
        delete [] prop_ub;
    }
};

void test_rnrm(const std::string& model, double end_time, double delta) {
    const std::string x = "../models/" + model + ".txt";
    const char* model_file = x.c_str();
    const std::string traj_file = "../data/" + model + "-rnrm-single_trajectory_out.csv";
    const std::string dist_file = "../data/" + model + "-rnrm-population_dist.csv";
    Model model_t(model_file);
    auto engine = RNRM(model_t, delta);
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

void profile_rnrm(const std::string& model, double end_time, double delta, int runs) {
    const std::string x = "../models/" + model + ".txt";
    const char* model_file = x.c_str();
    const std::string traj_file = "../data/" + model + "-rnrm-single_trajectory_out.csv";
    const std::string dist_file = "../data/" + model + "-rnrm-population_dist.csv";
    Model model_t(model_file);
    auto engine = RNRM(model_t, delta);
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
#endif //SIM_RNRM_H
