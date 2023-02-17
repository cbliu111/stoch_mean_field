//
// Created by cbl on 4/23/22.
//

#ifndef SIM_RANDOM_ENGINE_H
#define SIM_RANDOM_ENGINE_H

#include <iostream>
#include <limits>
#include <random>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/binomial_distribution.hpp>
#include <boost/random/gamma_distribution.hpp>

class RandEngine {
private:
public:
    boost::mt19937 mt19937;
    boost::uniform_01<> u01;
    explicit RandEngine() = default;
    explicit RandEngine(long long seed) {mt19937.seed(seed);}

    // change seed
    void seed(long long seed) {mt19937.seed(seed);}

    // discrete uniform random number
    double discreteUniform() {return mt19937();}

    // discrete uniform random number int [a, b]
    int discrete_uniform(int a, int b) {
        boost::uniform_int<> DiscreteUniform(a, b);
        boost::variate_generator<boost::mt19937&, boost::uniform_int<> > generatorDiscreteUniform(mt19937, DiscreteUniform);
        return (generatorDiscreteUniform());
    }

    // continuous uniform random number in (a,b)
    double continuous_open(double a, double b) {
        boost::uniform_01<> ContinuousZeroOneOpen;
        boost::variate_generator<boost::mt19937&, boost::uniform_01<> > generatorContinuousZeroOneOpen(mt19937, ContinuousZeroOneOpen);
        return (generatorContinuousZeroOneOpen()*(b-a) + a);
    }

    double uniform01() {
        // return a random value in range [0, 1)
        boost::variate_generator<boost::mt19937&, boost::uniform_01<> > g(mt19937, u01);
        return 1 - g();
    }

    // exponential random number
    double exponential(double mean) {
        if (mean < 0) {
            std::cerr << "Negative mean is given to the exponential random number generator" << std::endl;
            std::abort();
        }
        double lambda;
        if(mean==0) {
            lambda=std::numeric_limits<double>::max();
            std::cout<<"0 mean is given to the exponential random number generator"<<"\n";
        }
        else if(mean>=std::numeric_limits<double>::max()) {
            return mean;
        }
        else {
            lambda=1.0/mean;
        }
        boost::exponential_distribution<> ExponentialDistribution(lambda);
        boost::variate_generator<boost::mt19937&, boost::exponential_distribution<> > generatorExponential(mt19937, ExponentialDistribution);
        return generatorExponential();
    }

    // normal random number
    double normal(double mean, double var) {
        boost::normal_distribution<> NormalDistribution(mean, sqrt(var));
        boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > generatorNormal(mt19937, NormalDistribution);
        return generatorNormal();
    }

    // poisson random number
    double poisson(double mean) {
        if (mean == 0) {
            return 0;
        }
        else if(mean <= 700) {
            boost::poisson_distribution<> PoissonDistribution(mean);
            boost::variate_generator<boost::mt19937&, boost::poisson_distribution<> > generatorPoisson(mt19937, PoissonDistribution);
            return generatorPoisson();
        }
        else {
            double result, rounded_result;
            result = normal(mean, mean);
            rounded_result = floor(result);
            return (result - rounded_result) < 0.5 ? rounded_result : (rounded_result+1);
        }
    }

    // binomial random number
    double binomial(double n, double p) {
        double result, rounded_result;
        p=std::min(1.0, p);
        if(n<=20) {
            boost::binomial_distribution<> BinomialDistribution((int)n, p);
            boost::variate_generator<boost::mt19937&, boost::binomial_distribution<> > generatorBinomial(mt19937, BinomialDistribution);
            return generatorBinomial();
        }
        else {
            result=std::min(n, normal(n * p, n * p * (1 - p)));
            rounded_result=floor(result);
            return (result-rounded_result)<0.5?rounded_result:(rounded_result+1);
        }
    }

    double gamma(double a, double b) {
        boost::gamma_distribution<> GammaDistribution(a, b);
        boost::variate_generator<boost::mt19937&, boost::gamma_distribution<> > generatorGamma(mt19937, GammaDistribution);
        return generatorGamma();
    }
};

#endif //SIM_RANDOM_ENGINE_H
