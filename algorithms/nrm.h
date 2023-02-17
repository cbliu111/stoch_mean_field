//
// Created by cbl on 12/18/22.
//

#ifndef SIM_NRM_H
#define SIM_NRM_H

#include "timer.h"
#include "basemethod.h"

class NextReactionMethod : public BaseSimMethod {
public:
    double* old_time;
    double* old_prop;
    int id;
    explicit NextReactionMethod(Model& m) : BaseSimMethod(m) {
        old_time = new double [M];
        old_prop = new double [N];
        id = 0;
    }
    void sim() override {
        for (int i = 0; i < M; ++i) {
            double p = calculate_propensity(populations, model.reactions[i]);
            old_prop[i] = p;
            double r0 = rand.uniform01();
            old_time[i] = calculate_reaction_time(r0, p, time);
        }
        while ((time < end_time) || (step < steps)) {
            id = 0;
            for (int i = 0; i < M; ++i) {
                // find minimum time
                id = (old_time[id] < old_time[i]) ? id : i;
            }
            model.update_populations(populations, id);
            time = old_time[id];
            for (int i = 0; i < M; ++i) {
                double p = calculate_propensity(populations, model.reactions[i]);
                double new_time;
                if (i == id) {
                    double r1 = rand.uniform01();
                    new_time = calculate_reaction_time(r1, p, time);
                }
                else {
                    new_time = update_reaction_time(rand, p, old_prop[i], old_time[i], time);
                }
                old_prop[i] = p;
                old_time[i] = new_time;
            }
            step += 1;
        }
    }
    ~NextReactionMethod() {
        delete [] old_prop;
        delete [] old_time;
    }
};


void test_nrm(const std::string& model, double end_time) {
    const std::string x = "../models/" + model + ".txt";
    const char* model_file = x.c_str();
    const std::string traj_file = "../data/" + model + "-nrm-single_trajectory_out.csv";
    const std::string dist_file = "../data/" + model + "-nrm-population_dist.csv";
    Model model_t(model_file);
    auto engine = NextReactionMethod(model_t);
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

#endif //SIM_NRM_H
