//
// Created by cbl on 6/27/22.
//

#ifndef SIM_RLEAP_H
#define SIM_RLEAP_H

#include "timer.h"
#include "basemethod.h"

class RLEAP : public BaseSimMethod {
public:
    double delta; // population coarse-graining factor
    double tau; // time interval for generating event sequence
    double ua0; // total mean-field propensity
    double T; // time check point
    double min_leap_size;
    double max_leap_size;
    double n_rssa;
    double* pop_ub; // population upper bound
    double* pop_lb; // population lower bound
    double* prop_ub; // propensity upper bound
    double* prop_lb; // propensity lower bound
    std::vector<int> events; // list of event_vector
    explicit RLEAP(Model& m, double d) : BaseSimMethod(m) {
        delta = d;
        tau = 0;
        ua0 = 0;
        T = 0;
        max_leap_size = 1e4;
        n_rssa = 100;
        min_leap_size = n_rssa * delta;
        pop_ub = new double [N];
        pop_lb = new double [N];
        prop_ub = new double [M];
        prop_lb = new double [M];
    }
    void compute_propensities() {
        ua0 = 0;
        for (int i = 0; i < M; ++i) {
            double up = calculate_propensity(pop_ub, model.reactions[i]);
            prop_lb[i] = calculate_propensity(pop_lb, model.reactions[i]);
            prop_ub[i] = up;
            ua0 += up;
        }
    }
    void compute_bounds() {
        for (int i = 0; i < N; ++i) {
            double x = populations[i];
            if (x < n_rssa) {
                pop_lb[i] = (x > min_leap_size) ? (x - min_leap_size) : 0;
                pop_ub[i] = x + min_leap_size;
            }
            else {
                pop_lb[i] = std::ceil(x * (1 - delta));
                pop_ub[i] = std::floor(x * (1 + delta));
            }
        }
    }
    void generate_event_sequence() {
        events.clear();
        for (int i = 0; i < M; ++i) {
            int n = static_cast<int>(rand.poisson(prop_ub[i] * (T-time)));
            for (int j = 0; j < n; ++j) {
                events.push_back(i);
            }
        }
    }
    int next_event(size_t i) {
        int e = rand.discrete_uniform(0, static_cast<int>(i)-1); // include i - 1
        int id = events[e];
        events[e] = events[i-1];
        return id;
    }
    bool check_accept(int id) {
        bool accept = false;
        double r1 = rand.uniform01();
        if (r1 <= prop_lb[id]/prop_ub[id]) accept = true;
        else {
            double p = calculate_propensity(populations, model.reactions[id]);
            if (r1 <= p / prop_ub[id]) accept = true;
        }
        return accept;
    }
    bool check_out_range(int id) {
        bool flag = false;
        for (int j : model.reactions[id].species) {
            double x = populations[j];
            if ((x > pop_ub[j]) || (x < pop_lb[j])) {
                flag = true;
                break;
            }
        }
        return flag;
    }

    void sim() override {
        T = 0;
        while ((time < end_time) || (step < steps)) {
            compute_bounds();
            compute_propensities();
            tau = max_leap_size / ua0;
            while (time >= T) {
                T += tau;
            }
            if (ua0 == 0) break;
            double event_count = 0;
            bool out_range = false;
            if ((T - time) < min_leap_size / ua0 || M < 5) {
                while (!out_range) {
                    double sum = 0;
                    int idx = 0;
                    double r1 = rand.uniform01() * ua0;
                    for (; idx < M; ++idx) {
                        sum += prop_ub[idx];
                        if (r1 < sum) break;
                    }
                    event_count += 1;
                    bool accept = check_accept(idx);
                    if (accept) {
                        step += 1;
                        model.update_populations(populations, idx);
                        out_range = check_out_range(idx);
                    }
                }
            }
            else {
                generate_event_sequence();
                if (events.empty()) {
                    // do not set time = T, not accurate
                    continue;
                }
                for (size_t i = events.size(); i > 0; --i) {
                    int id = next_event(i);
                    event_count += 1;
                    bool accept = check_accept(id);
                    if (accept) {
                        step += 1;
                        model.update_populations(populations, id);
                        out_range = check_out_range(id);
                        if (out_range) break;
                    }
                }
            }
            if (event_count == 0) continue;
            time += rand.gamma(event_count, 1/ua0);
            record_trajectory();
            if (time >= end_time) break;
        }
    }
    ~RLEAP() {
        delete [] pop_ub;
        delete [] pop_lb;
        delete [] prop_ub;
        delete [] prop_lb;
    }
};

void test_rleap(const std::string& model, double end_time, double delta) {
    const std::string x = "../models/" + model + ".txt";
    const char* model_file = x.c_str();
    const std::string traj_file = "../data/" + model + "-rleap-single_trajectory_out.csv";
    const std::string dist_file = "../data/" + model + "-rleap-population_dist.csv";
    Model model_t(model_file);
    auto engine = RLEAP(model_t, delta);
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

void profile_rleap(const std::string& model, double end_time, double delta, int runs) {
    const std::string x = "../models/" + model + ".txt";
    const char* model_file = x.c_str();
    const std::string traj_file = "../data/" + model + "-rleap-single_trajectory_out.csv";
    const std::string dist_file = "../data/" + model + "-rleap-population_dist.csv";
    Model model_t(model_file);
    auto engine = RLEAP(model_t, delta);
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
#endif //SIM_RLEAP_H
