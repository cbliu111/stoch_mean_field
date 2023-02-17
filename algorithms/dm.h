//
// Created by cbl on 12/18/22.
//

#ifndef SIM_DM_H
#define SIM_DM_H

#include "timer.h"
#include "basemethod.h"

class DirectMethod : public BaseSimMethod {
public:
    double* propensities;
    explicit DirectMethod(Model& m) : BaseSimMethod(m) {
        propensities = new double [M];
    }
    void sim() override {
        while ((time < end_time) || (step < steps)) {
            double a0 = 0;
            for (int i = 0; i < M; ++i) {
                double temp = calculate_propensity(populations, model.reactions[i]);
                a0 += temp;
                propensities[i] = temp;
            }
            if (a0 == 0) break;
            // time should be calculated before setting sum to 0
            double r1 = rand.uniform01();
            double dt = std::log(1.0 / r1) / a0;
            time += dt;
            double r = rand.uniform01() * a0;
            double partial_prop = 0;
            for (int i = 0; i < M; ++i) {
                partial_prop += propensities[i];
                if (partial_prop >= r) {
                    model.update_populations(populations, i);
                    break;
                }
            }
            record_trajectory();
            step += 1;
            if (time > end_time) break;
        }
    }
    ~DirectMethod() {
        delete [] propensities;
    }
};

void test_dm(const std::string& model, double end_time) {
    const std::string x = "../models/" + model + ".txt";
    const char* model_file = x.c_str();
    const std::string traj_file = "../data/" + model + "-dm-single_trajectory_out.csv";
    const std::string dist_file = "../data/" + model + "-dm-population_dist.csv";
    Model model_t(model_file);
    auto engine = DirectMethod(model_t);
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

void profile_dm(const std::string& model, double end_time, int runs) {
    const std::string x = "../models/" + model + ".txt";
    const char* model_file = x.c_str();
    const std::string traj_file = "../data/" + model + "-dm-single_trajectory_out.csv";
    const std::string dist_file = "../data/" + model + "-dm-population_dist.csv";
    Model model_t(model_file);
    auto engine = DirectMethod(model_t);
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
#endif //SIM_DM_H
