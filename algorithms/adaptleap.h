//
// Created by cbl on 12/20/22.
//

#ifndef SIM_ADAPTLEAP_H
#define SIM_ADAPTLEAP_H

#include "timer.h"
#include "basemethod.h"

/*
learn a proper proposal distribution for importance sampling 
time intervals of the history sampling are recorded and averaged for the next sampling
*/

class ADAPTLEAP : public BaseLeapMethod {
public:
    int n_leap; // number of leap during simulation
    int n_traj; // number of sampled trajectory
    int n_prerun; // number of prerun trajectory
    std::vector<double> dts; // time intervals of leaps of coarse grained process
    std::vector<double> ns; // number of parameter samples for coarse grained process
    explicit ADAPTLEAP(Model& m, double bound_weight) : BaseLeapMethod(m, bound_weight) {
        n_leap = 0;
        n_traj = 0;
        //n_prerun = 100;
        n_prerun = 1;
        std::cout << epsilon << std::endl;
    }
    void collect_parameters(double t) {
        if (n_leap >= ns.size()) {
            dts.push_back(t);
            ns.push_back(1);
        }
        else {
            dts[n_leap] += t;
            ns[n_leap] += 1;
        }
    }
    void sim() override {
        n_traj += 1;
        std::cout << "traj : " << n_traj << std::endl;
        n_leap = 0;
        // prerun steps 
        if (n_traj <= n_prerun) {
            while ((time < end_time) || (step < steps)) {
                compute_tau();
                if (ua0 == 0) break;
                generate_events_poisson();
                generate_event_sequence();
                double lt = 0;
                for (size_t i = events.size(); i > 0; --i) {
                    int id = next_event(i);
                    double dt = std::log(1.0 / rand.uniform01()) / ua0;
                    double w = compute_step_weight(id, dt);
                    if (w > bw) {
                        lt += dt;
                        step += 1;
                        model.update_populations(populations, id);
                    }
                    else break;
                    if (time >= end_time) break;
                }
                time += lt;
                collect_parameters(lt);
                n_leap += 1;
                record_trajectory();
                if (time >= end_time) break;
            }
            if (n_traj == n_prerun) {
                // calculate trajectory mean of propensities and dts
                for (int j = 0; j < ns.size(); ++j) {
                    dts[j] = dts[j] / ns[j];
                }
            }
        }
        else {
            while ((time < end_time) || (step < steps)) {
                compute_mean_field_propensities();
                if (ua0 == 0) break;
                if (n_leap >= dts.size()) {
                    compute_tau();
                    //tau = end_time / 100;
                }
                else {
                    tau = dts[n_leap];
                }
                backup_state();
                for (int i = 0; i < M; ++i) {
                    double p = prop_mf[i];
                    double t = tau;
                    double n = rand.poisson(p * t);
                    step += n;
                    model.update_populations_leap(populations, i, n);
                }
                bool success = check_state_bound();
                if (success) {
                    time += tau;
                    n_leap += 1;
                }
                else {
                    restore_state();
                    generate_event_sequence();
                    event_count = 0;
                    for (size_t i = events.size(); i > 0; --i) {
                        int id = next_event(i);
                        double dt = std::log(1.0 / rand.uniform01()) / ua0;
                        double w = compute_step_weight(id, dt);
                        if (w > bw) {
                            event_count += 1;
                            step += 1;
                            model.update_populations(populations, id);
                        }
                        else break;
                        if (time >= end_time) break;
                    }
                    time += rand.gamma(event_count, 1/ua0);
                    n_leap += 1;
                }
                record_trajectory();
                if (time >= end_time) break;
            }
        }
    }
    ~ADAPTLEAP() = default;
};

void test_adaptleap(const std::string& model, double end_time, double bound_weight) {
    const std::string x = "../models/" + model + ".txt";
    const char* model_file = x.c_str();
    const std::string traj_file = "../data/" + model + "-adaptleap-bw_" + std::to_string(bound_weight) + "-single_trajectory_out.csv";
    const std::string dist_file = "../data/" + model + "-adaptleap-bw_" + std::to_string(bound_weight) + "-population_dist.csv";
    Model model_t(model_file);
    auto engine = ADAPTLEAP(model_t, bound_weight);
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

void profile_adaptleap(const std::string& model, double end_time, double bound_weight, int runs) {
    const std::string x = "../models/" + model + ".txt";
    const char* model_file = x.c_str();
    const std::string traj_file = "../data/" + model + "-adaptleap-bw_" + std::to_string(bound_weight) + "-single_trajectory_out.csv";
    const std::string dist_file = "../data/" + model + "-adaptleap-bw_" + std::to_string(bound_weight) + "-population_dist.csv";
    Model model_t(model_file);
    auto engine = ADAPTLEAP(model_t, bound_weight);
    double et = end_time; // end time

    // single trajectory
    engine.set_end_time(et);
    engine.set_output_file(traj_file);
    engine.set_number_data_points(1000);
    engine.set_number_trajectory(1);
    engine.set_seed(0);
    engine.run();
    engine.set_number_trajectory(1);
    engine.set_seed(0);
    engine.run();
    // bench mark
    if (runs <= 0) return;
    engine.set_end_time(et);
    engine.set_output_file("");
    engine.set_number_data_points(0);
    engine.set_number_trajectory(100);
    engine.set_seed(0);
    engine.run();

    engine.set_number_trajectory(runs);
    engine.set_seed(0);
    engine.run();

}
#endif //SIM_ADAPTLEAP_H
