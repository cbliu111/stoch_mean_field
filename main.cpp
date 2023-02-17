#include <iostream>
#include <fstream>
#include <random>
#include "algorithms/dm.h"
#include "algorithms/nrm.h"
#include "algorithms/rssa.h"
#include "algorithms/rnrm.h"
#include "algorithms/rleap.h"
#include "algorithms/tleap.h"
#include "algorithms/cme.h"
#include "algorithms/barssa.h"
#include "algorithms/ileap.h"
#include "algorithms/nbwleap.h"
#include "algorithms/bswleap.h"
#include "algorithms/blwleap.h"
#include "algorithms/btwleap.h"
#include "algorithms/adaptleap.h"
#include "algorithms/scaleleap.h"

class ParameterDict {
public:
    std::vector<std::string> methods;
    double end_time{};
    double leap_time_interval{};
    double state_delta{};
    double tau_epsilon{};
    double ba_alpha{};
    double bound_weight{};
    explicit ParameterDict() = default;
};

void test(const std::string& model, ParameterDict& param) {
    std::cout << "Model: " << model << std::endl;
    for (std::string& m : param.methods) {
        if (m == "DM") {
            std::cout << "----------" << "DM" << "----------" << std::endl;
            test_dm(model, param.end_time);
        }
        else if (m == "RNRM") {
            std::cout << "----------" << "RNRM" << "----------" << std::endl;
            test_rnrm(model, param.end_time, param.state_delta);
        }
        else if (m == "RSSA") {
            std::cout << "----------" << "RSSA" << "----------" << std::endl;
            test_rssa(model, param.end_time, param.state_delta);
        }
        else if (m == "RLEAP") {
            std::cout << "----------" << "RLEAP" << "----------" << std::endl;
            test_rleap(model, param.end_time, param.state_delta);
        }
        else if (m == "ILEAP") {
            std::cout << "----------" << "ILEAP" << "----------" << std::endl;
            test_ileap(model, param.end_time, param.bound_weight);
        }
        else if (m == "TLEAP") {
            std::cout << "----------" << "TLEAP" << "----------" << std::endl;
            test_tleap(model, param.end_time, param.tau_epsilon);
        }
        else if (m == "BARSSA") {
            std::cout << "----------" << "BARSSA" << "----------" << std::endl;
            test_barssa(model, param.end_time, param.ba_alpha);
        }
        else if (m == "NBWLEAP") {
            std::cout << "----------" << "NBWLEAP" << "----------" << std::endl;
            test_nbwleap(model, param.end_time, param.leap_time_interval);
        }
        else if (m == "BSWLEAP") {
            std::cout << "----------" << "BSWLEAP" << "----------" << std::endl;
            test_bswleap(model, param.end_time, param.leap_time_interval, param.bound_weight);
        }
        else if (m == "BLWLEAP") {
            std::cout << "----------" << "BLWLEAP" << "----------" << std::endl;
            test_blwleap(model, param.end_time, param.leap_time_interval, param.bound_weight);
        }
        else if (m == "BTWLEAP") {
            std::cout << "----------" << "BTWLEAP" << "----------" << std::endl;
            test_btwleap(model, param.end_time, param.leap_time_interval, param.bound_weight);
        }
        else if (m == "ADAPTLEAP") {
            std::cout << "----------" << "ADAPTLEAP" << "----------" << std::endl;
            test_adaptleap(model, param.end_time, param.bound_weight);
        }
    }
}

void accuracy_efficiency_nonstiff_dimer() {
    const std::string model = "dimer";
    std::vector<std::string> ms;
    std::vector<std::string> conditions;
    ms.emplace_back("DM");
    ms.emplace_back("RSSA");
    ms.emplace_back("RLEAP");
    ParameterDict p;
    p.methods = ms;
    p.end_time = 9.999;
    p.leap_time_interval = p.end_time / 100;
    p.state_delta = 0.01;
    p.tau_epsilon = 0.001;
    p.ba_alpha = 0.6;
    p.bound_weight = 0.93;
    // run exact methods
    //test(model, p); // completed
    // run tau leap
    ms.clear();
    ms.emplace_back("TLEAP");
    p.methods = ms;
    conditions.clear();
    conditions.assign({"0.001", "0.005", "0.01", "0.03", "0.04", "0.05"});
    for (std::string& x : conditions) {
        std::cout << "tau = " << x << std::endl;
        p.tau_epsilon = std::stod(x);
        //test(model, p);
    }
    // run barssa
    ms.clear();
    ms.emplace_back("BARSSA");
    p.methods = ms;
    conditions.clear();
    conditions.assign({"0.95", "0.9", "0.8", "0.7", "0.6"});
    for (std::string& x : conditions) {
        std::cout << "alpha = " << x << std::endl;
        p.ba_alpha = std::stod(x);
        test(model, p);
    }
    // run adaptleap
    ms.clear();
    ms.emplace_back("ADAPTLEAP");
    p.methods = ms;
    conditions.clear();
    conditions.assign({"0.99", "0.97", "0.95", "0.93", "0.9"});
    for (std::string& x : conditions) {
        std::cout << "bw = " << x << std::endl;
        p.bound_weight = std::stod(x);
        test(model, p);
    }
}

void accuracy_efficiency_schlogl() {
    const std::string model = "schlogl";
    std::vector<std::string> ms;
    std::vector<std::string> conditions;
    ms.emplace_back("DM");
    //ms.emplace_back("RSSA");
    //ms.emplace_back("RLEAP");
    ParameterDict p;
    p.methods = ms;
    p.end_time = 9.999;
    p.leap_time_interval = p.end_time / 100;
    p.state_delta = 0.05;
    p.tau_epsilon = 0.001;
    p.ba_alpha = 0.6;
    p.bound_weight = 0.93;
    // run exact methods
    test(model, p); //completed
    // run tau leap
    ms.clear();
    ms.emplace_back("TLEAP");
    p.methods = ms;
    conditions.clear();
    //conditions.assign({"0.001", "0.01", "0.1", "0.15", "0.18", "0.21", "0.24", "0.27", "0.3"});
    conditions.assign({"0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"});
    for (std::string& x : conditions) {
        std::cout << "tau = " << x << std::endl;
        p.tau_epsilon = std::stod(x);
        //test(model, p);
    }
    // run barssa
    ms.clear();
    ms.emplace_back("BARSSA");
    p.methods = ms;
    conditions.clear();
    conditions.assign({"0.95", "0.9", "0.8", "0.7", "0.6", "0.5"});
    for (std::string& x : conditions) {
        std::cout << "alpha = " << x << std::endl;
        p.ba_alpha = std::stod(x);
        //test(model, p);
    }
    // run adaptleap
    ms.clear();
    ms.emplace_back("ADAPTLEAP");
    p.methods = ms;
    conditions.clear();
    //conditions.assign({"0.99", "0.97", "0.95", "0.93", "0.9", "0.8", "0.7", "0.6", "0.5"});
    conditions.assign({"0.4", "0.3", "0.2", "0.1", "0.05"});
    for (std::string& x : conditions) {
        std::cout << "bw = " << x << std::endl;
        p.bound_weight = std::stod(x);
        //test(model, p);
    }
}

void scale_MAPK() {
    std::vector<std::string> sf {"", "x10", "x20", "x40", "x60", "x80", "x100"};
    for (size_t s = 0; s < sf.size(); ++s) {
        int runs = 1;
        const std::string model = "MAPK" + sf[s];
        double et = 0.05;
        runs = 0;
        std::cout << "scale factor : " << sf[s] << std::endl;
        std::cout << "----------DM----------" << std::endl;
        //profile_dm(model, et, runs);
        std::cout << "----------RSSA----------" << std::endl;
        //profile_rssa(model, et, 0.01, runs);
        std::cout << "----------RLEAP----------" << std::endl;
        //profile_rleap(model, et, 0.01, runs);
        std::cout << "----------TLEAP----------" << std::endl;
        //profile_tleap(model, et, 0.01, runs);
        std::cout << "----------BARSSA----------" << std::endl;
        //profile_barssa(model, et, 0.9, runs);
        std::cout << "----------ADAPTLEAP----------" << std::endl;
        profile_adaptleap(model, et, 0.99, runs);
        std::cout << "----------SCALELEAP----------" << std::endl;
        //profile_scaleleap(model, et, 0.99, runs);
    }
}

void accuracy_bsw_blw() {
    std::vector<std::string> conditions {"0.99", "0.95", "0.9", "0.8", "0.7", "0.6"};
    std::string model = "birth_death";
    std::cout << "Model: " << model << std::endl;
    double leap_time = 10;
    double et = 999.9;
    std::cout << "Method: No weight bound" << std::endl;
    //test_nbwleap(model,  et, leap_time);
    for (auto& s : conditions) {
        std::cout << "Method: bound step weight at : " << s << std::endl;
        test_bswleap(model, et, leap_time, std::stod(s));
        std::cout << "Method: bound mean leap weight at : " << s << std::endl;
        //test_blwleap(model, et, leap_time, std::stod(s));
    }
    model = "isomerization";
    leap_time = 10;
    et = 99.99;
    std::cout << "Model: " << model << std::endl;
    std::cout << "Method: No weight bound" << std::endl;
    //test_nbwleap(model,  et, leap_time);
    for (auto& s : conditions) {
        std::cout << "Method: bound step weight at : " << s << std::endl;
        test_bswleap(model, et, leap_time, std::stod(s));
        std::cout << "Method: bound mean leap weight at : " << s << std::endl;
        //test_blwleap(model, et, leap_time, std::stod(s));
    }
    model = "sequence";
    leap_time = 0.001;
    et = 0.0999;
    std::cout << "Model: " << model << std::endl;
    std::cout << "Method: No weight bound" << std::endl;
    //test_nbwleap(model,  et, leap_time);
    for (auto& s : conditions) {
        std::cout << "Method: bound step weight at : " << s << std::endl;
        test_bswleap(model, et, leap_time, std::stod(s));
        std::cout << "Method: bound mean leap weight at : " << s << std::endl;
        //test_blwleap(model, et, leap_time, std::stod(s));
    }
    model = "schlogl";
    leap_time = 9.999;
    et = 1;
    std::cout << "Model: " << model << std::endl;
    std::cout << "Method: No weight bound" << std::endl;
    //test_nbwleap(model,  et, leap_time);
    for (auto& s : conditions) {
        std::cout << "Method: bound step weight at : " << s << std::endl;
        test_bswleap(model, et, leap_time, std::stod(s));
        std::cout << "Method: bound mean leap weight at : " << s << std::endl;
        //test_blwleap(model, et, leap_time, std::stod(s));
    }
}

int main() {
    accuracy_efficiency_nonstiff_dimer();
    accuracy_efficiency_schlogl();
    scale_MAPK();
    accuracy_bsw_blw();
}

