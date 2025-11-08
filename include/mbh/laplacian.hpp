#pragma once
#include "grid.hpp"

namespace mbh {

inline void laplace7(const Field& u, Field& Lu){
    const int nx=u.g.nx, ny=u.g.ny, nz=u.g.nz;
    const double idx2 = 1.0/(u.g.dx*u.g.dx);
    const double idy2 = 1.0/(u.g.dy*u.g.dy);
    const double idz2 = 1.0/(u.g.dz*u.g.dz);

    for(int k=1;k<=nz;++k){
        for(int j=1;j<=ny;++j){
            for(int i=1;i<=nx;++i){

                double c = -2.0*(idx2+idy2+idz2)*u(i,j,k);

                double s = (u(i-1,j,k)+u(i+1,j,k))*idx2
                         + (u(i,j-1,k)+u(i,j+1,k))*idy2
                         + (u(i,j,k-1)+u(i,j,k+1))*idz2;
                         
                Lu(i,j,k) = c + s;
                
            }
        }
    }
}

}
