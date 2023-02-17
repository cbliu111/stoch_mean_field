//
// Created by cbl on 1/23/23.
//

#ifndef SIM_ILEAP_H
#define SIM_ILEAP_H

/*
sequential importance sampling with mean-field randomization
rejection control is used to avoid small weight samples
this method is exact
*/

#include "timer.h"
#include "basemethod.h"

class ILEAP : public BaseLeapMethod {
public:
    explicit ILEAP(Model& m, double bound_weight) : BaseLeapMethod(m, bound_weight) {}
    void sim() override {
        while ((time < end_time) || (step < steps)) {
            compute_mean_field_propensities();
            tau = default_leap_size / ua0;
            generate_events_poisson();
            generate_event_sequence();
            for (size_t i = events.size(); i > 0; --i) {
                int id = next_event(i);
                double dt = std::log(1.0 / rand.uniform01()) / ua0;
                double w = compute_step_weight(id, dt);
                w *= weight; // limiting the trajectory weight
                double th = (1 < w / bw) ? 1 : w / bw;
                double r = rand.uniform01();
                if (r < th) {
                    weight = (w > bw) ? w : bw;
                    time += dt;
                    step += 1;
                    model.update_populations(populations, id);
                }
                else break;
                if (time >= end_time) break;
            }
            record_trajectory();
        }
    }
    ~ILEAP() = default;
};

void test_ileap(const std::string& model, double end_time, double bound_weight) {
    const std::string x = "../models/" + model + ".txt";
    const char* model_file = x.c_str();
    const std::string traj_file = "../data/" + model + "-ileap-single_trajectory_out.csv";
    const std::string dist_file = "../data/" + model + "-ileap-population_dist.csv";
    Model model_t(model_file);
    auto engine = ILEAP(model_t, bound_weight);
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
#endif //SIM_ILEAP_H
