//
// Created by cbl on 12/20/22.
//

#ifndef SIM_BLWLEAP_H
#define SIM_BLWLEAP_H

#include "timer.h"
#include "basemethod.h"

/*
approximate leap weight method
instead of bound step weight, bound leap weight to increase sampling efficiency
*/

class BLWLEAP : public BaseLeapMethod {
public:
    explicit BLWLEAP(Model& m, double time_interval, double bound_weight) : BaseLeapMethod(m, bound_weight, -1, time_interval) {}
    void step_wise_update() {
    }
    double compute_mean_weight()  {
        for (int i = 0; i < N; ++i) {
            pop_lm[i] = (pop_copy[i] + populations[i]) / 2;
        }
        double a_lm = 0;
        double w = 1;
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
    void sim() override {
        bool retry = false;
        double dt = 0;
        while ((time < end_time) || (step < steps)) {
            // compute properties of mean-field state
            compute_mean_field_propensities();
            if (!retry) {
                compute_tau();
                dt = tau;
            }
            backup_state();
            for (int i = 0; i < M; ++i) {
                double n = rand.poisson(prop_mf[i] * dt);
                event_count += n;
                model.update_populations_leap(populations, i, n);
            }
            bool success = check_state_bound();
            if (success) {
                double w = compute_mean_weight();
                if (w > bw) {
                    retry = false;
                    step += event_count;
                    time += dt;
                }
                else {
                    restore_state();
                    dt /= 2;
                    retry = true;
                    continue;
                }
            }
            else {
                restore_state();
                dt /= 2;
                retry = true;
                continue;
            }
            record_trajectory();
            if (time >= end_time) break;
        }
    }
    ~BLWLEAP() = default;
};

void test_blwleap(const std::string& model, double end_time, double time_interval, double bound_weight) {
    const std::string x = "../models/" + model + ".txt";
    const char* model_file = x.c_str();
    const std::string traj_file = "../data/" + model + "-blwleap-bw_" + std::to_string(bound_weight) + "-single_trajectory_out.csv";
    const std::string dist_file = "../data/" + model + "-blwleap-bw_" + std::to_string(bound_weight) + "-population_dist.csv";
    Model model_t(model_file);
    auto engine = BLWLEAP(model_t, time_interval, bound_weight);
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
#endif //SIM_BLWLEAP_H
