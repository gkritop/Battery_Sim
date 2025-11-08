#pragma once
#include <cmath>
#include "params.hpp"
#include "grid.hpp"

namespace mbh {

inline double arrhenius_rate(double T, const Arrhenius& a){
    // Very simplified toy reaction rate. In real life you'd have R*T etc.
    // Make sure it's numerically tame near 300-500K.
    const double R = 8.314;
    double k = a.A * std::exp(-a.Ea/(R*std::max(1.0, T)));
    return k;
}

inline void source_arrhenius(const Params& p, const Field& T, Field& S){
    int nx=T.g.nx, ny=T.g.ny, nz=T.g.nz;
    for(int k=1;k<=nz;++k){
        for(int j=1;j<=ny;++j){
            for(int i=1;i<=nx;++i){
                double r = arrhenius_rate(T(i,j,k), p.arr);
                S(i,j,k) = (p.enable_source ? p.arr.H * r / (p.rho*p.cp) : 0.0);
            }
        }
    }
}

}
