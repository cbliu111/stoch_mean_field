//
// Created by cbl on 12/20/22.
//

#ifndef SIM_BTWLEAP_H
#define SIM_BTWLEAP_H

#include "timer.h"
#include "basemethod.h"

/*
approximate leap weight method
instead of bound step weight, bound leap weight to increase sampling efficiency
*/

class BTWLEAP : public BaseLeapMethod {
public:
    explicit BTWLEAP(Model& m, double time_interval, double bound_weight) : BaseLeapMethod(m, bound_weight, -1, time_interval) {}
    void sim() override {
        while ((time < end_time) || (step < steps)) {
            // compute properties of mean-field state
            compute_mean_field_propensities();
            generate_events_poisson();
            generate_event_sequence();
            for (size_t i = events.size(); i > 0; --i) {
                int id = next_event(i);
                double dt = std::log(1.0 / rand.uniform01()) / ua0;
                double w = compute_step_weight(id, dt);
                weight *= w; // accumulate trajectory weight
                double th = (1 < weight / bw) ? 1 : weight / bw;
                double r = rand.uniform01();
                if (r < th) {
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
    ~BTWLEAP() = default;
};

void test_btwleap(const std::string& model, double end_time, double time_interval, double bound_weight) {
    const std::string x = "../models/" + model + ".txt";
    const char* model_file = x.c_str();
    const std::string traj_file = "../data/" + model + "-btwleap-single_trajectory_out.csv";
    const std::string dist_file = "../data/" + model + "-btwleap-population_dist.csv";
    Model model_t(model_file);
    auto engine = BTWLEAP(model_t, time_interval, bound_weight);
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
#endif //SIM_BTWLEAP_H
