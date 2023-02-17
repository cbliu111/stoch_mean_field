//
// Created by cbl on 12/20/22.
//

#ifndef SIM_SCALELEAP_H
#define SIM_SCALELEAP_H

#include "timer.h"
#include "basemethod.h"

/*
learn a proper proposal distribution for importance sampling 
time intervals of the history sampling are recorded and averaged for the next sampling
*/

class SCALELEAP : public BaseLeapMethod {
public:
    int n_leap; // number of leap during simulation
    int n_traj; // number of sampled trajectory
    int n_prerun; // number of prerun trajectory
    std::vector<double> dts; // time intervals of leaps of coarse grained process
    std::vector<double> ns; // number of parameter samples for coarse grained process
    explicit SCALELEAP(Model& m, double bound_weight) : BaseLeapMethod(m, bound_weight) {
        n_leap = 0;
        n_traj = 0;
        n_prerun = 100;
        std::cout << epsilon << std::endl;
    }
    void do_leap() {
        for (int i = 0; i < M; ++i) {
            double n = rand.poisson(prop_mf[i] * tau);
            event_count += n;
            model.update_populations_leap(populations, i, n);
        }
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
        while ((time < end_time) || (step < steps)) {
            if (!retry) {
                compute_tau();
            }
            if (ua0 == 0) break;
            event_count = 0;
            backup_state();
            do_leap();
            bool success = check_state_bound();
            if (success) {
                double w = compute_mean_weight();
                if (w > bw) {
                    retry = false;
                    step += event_count;
                    time += tau;
                    n_leap += 1;
                    if (time >= end_time) break;
                }
                else {
                    restore_state();
                    tau /= 2;
                    retry = true;
                    continue;
                    //do_ssa(num_dm);
                    //retry = false;
                }
            }
            else {
                restore_state();
                tau /= 2;
                retry = true;
                continue;
            }
            record_trajectory();
            if (time >= end_time) break;
        }
    }
    ~SCALELEAP() = default;
};

void test_scaleleap(const std::string& model, double end_time, double bound_weight) {
    const std::string x = "../models/" + model + ".txt";
    const char* model_file = x.c_str();
    const std::string traj_file = "../data/" + model + "-scaleleap-bw_" + std::to_string(bound_weight) + "-single_trajectory_out.csv";
    const std::string dist_file = "../data/" + model + "-scaleleap-bw_" + std::to_string(bound_weight) + "-population_dist.csv";
    Model model_t(model_file);
    auto engine = SCALELEAP(model_t, bound_weight);
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

void profile_scaleleap(const std::string& model, double end_time, double bound_weight, int runs) {
    const std::string x = "../models/" + model + ".txt";
    const char* model_file = x.c_str();
    const std::string traj_file = "../data/" + model + "-scaleleap-bw_" + std::to_string(bound_weight) + "-single_trajectory_out.csv";
    const std::string dist_file = "../data/" + model + "-scaleleap-bw_" + std::to_string(bound_weight) + "-population_dist.csv";
    Model model_t(model_file);
    auto engine = SCALELEAP(model_t, bound_weight);
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
#endif //SIM_SCALELEAP_H
