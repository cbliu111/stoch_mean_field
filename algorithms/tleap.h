
//
// Created by cbl on 12/20/22.
//

#ifndef SIM_TLEAP_H
#define SIM_TLEAP_H

#include "timer.h"
#include "basemethod.h"

class TLEAP : public BaseLeapMethod {
public:
    explicit TLEAP(Model& m, double epsilon) : BaseLeapMethod(m, -1, epsilon, -1) {}
    void sim() override {
        while ((time < end_time) || (step < steps)) {
            compute_tau();
            if (ua0 == 0) break;
            backup_state();
            bool do_leap = true;
            while (do_leap) {
                if (tau > leap_threshold / ua0) {
                    event_count = 0;
                    for (int i = 0; i < M; ++i) {
                        double n = rand.poisson(prop_mf[i] * tau);
                        event_count += n;
                        model.update_populations_leap(populations, i, n);
                    }
                    bool success = check_state_bound();
                    if (success) {
                        time += tau;
                        step += event_count;
                        record_trajectory();
                        do_leap = false;
                    }
                    else {
                        restore_state();
                        tau /= 2;
                        continue;
                    }
                }
                else {
                    do_ssa(num_dm);
                    do_leap = false;
                }
            }
            if (time > end_time) break;
        }
    }
    ~TLEAP() = default;
};

void test_tleap(const std::string& model, double end_time, double epsilon) {
    const std::string x = "../models/" + model + ".txt";
    const char* model_file = x.c_str();
    const std::string traj_file = "../data/" + model + "-tleap-tau_" + std::to_string(epsilon) + "-single_trajectory_out.csv";
    const std::string dist_file = "../data/" + model + "-tleap-tau_" + std::to_string(epsilon) + "-population_dist.csv";
    Model model_t(model_file);
    auto engine = TLEAP(model_t, epsilon);
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

void profile_tleap(const std::string& model, double end_time, double epsilon, int runs) {
    const std::string x = "../models/" + model + ".txt";
    const char* model_file = x.c_str();
    const std::string traj_file = "../data/" + model + "-tleap-tau_" + std::to_string(epsilon) + "-single_trajectory_out.csv";
    const std::string dist_file = "../data/" + model + "-tleap-tau_" + std::to_string(epsilon) + "-population_dist.csv";
    Model model_t(model_file);
    auto engine = TLEAP(model_t, epsilon);
    double et = end_time; // end time

    // single trajectory
    engine.set_end_time(et);
    engine.set_output_file(traj_file);
    engine.set_number_data_points(1000);
    engine.set_number_trajectory(1);
    engine.set_seed(clock());
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
#endif //SIM_TLEAP_H
