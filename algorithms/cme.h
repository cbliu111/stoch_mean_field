//
// Created by cbl on 7/5/22.
//

#ifndef SIM_CME_H
#define SIM_CME_H

#include <fstream>

void CME_sequential_model() {
    double dt = 0.00001;
    double** pden = new double*[1100];
    for (int i = 0; i < 1100; ++i) {
        pden[i] = new double[1100];
    }
    for (int i = 0; i < 1100; ++i) {
        for (int j = 0; j < 1100; ++j) {
            pden[i][j] = 0;
        }
    }
    pden[1000][0] = 1;
    double t = 0;
    while (t <= 0.10001) {
        t += dt;
        for (int i = 0; i < 1001; ++i) {
            for (int j = 0; j < 1001; ++j) {
                double a = (i+9000+1) * pden[i+1][j-1];
                double b = (j+1) * pden[i][j+1];
                double c = (i+9000+j) * pden[i][j];
                pden[i][j] += (a + b - c) * dt;
            }
        }
    }
    double* phist = new double[1001];
    for (int i = 0; i < 1000; ++i) {
        for (int j = 0; j < 1000; ++j) {
            phist[1000-i-j] += pden[i][j];
        }
    }
    std::string filename = "sequence_CME.csv";
    std::ofstream fo;
    fo.precision(17);
    fo.open(filename, std::ios::out);
    for (int i = 0; i < 1001; ++i) {
        fo << phist[i] << std::endl;
    }
    fo << std::endl;
    fo.close();
    for (int i = 0; i < 1100; ++i) {
        delete [] pden[i];
    }
    delete [] pden;
    delete [] phist;
}

#endif //SIM_CME_H
