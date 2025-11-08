#pragma once
#include "laplacian.hpp"
#include "source.hpp"
#include "halo.hpp"

namespace mbh {

struct Coeff {
    double alpha; // k/(rho*cp)
};

inline Coeff make_coeff(const Params& p){
    return { p.k/(p.rho*p.cp) };
}

// Simple RK2 (Heun) for dT/dt = alpha*L(T) + S(T)
inline void step_rk2(const Params& p, const MPIBox& box, const Coeff& c, Field& T, Field& tmp, Field& Lu, Field& S){
    // halos needed for Laplacian
    exchange_halo(box, T);
    apply_dirichlet_bc(box, T, p.T_env);

    laplace7(T, Lu);
    source_arrhenius(p, T, S);

    int nx=T.g.nx, ny=T.g.ny, nz=T.g.nz;
    for(int k=1;k<=nz;++k) for(int j=1;j<=ny;++j) for(int i=1;i<=nx;++i){
        tmp(i,j,k) = T(i,j,k) + p.dt * (c.alpha*Lu(i,j,k) + S(i,j,k));
    }

    exchange_halo(box, tmp);
    apply_dirichlet_bc(box, tmp, p.T_env);

    laplace7(tmp, Lu);
    source_arrhenius(p, tmp, S);

    for(int k=1;k<=nz;++k) for(int j=1;j<=ny;++j) for(int i=1;i<=nx;++i){
        double f1 = c.alpha*Lu(i,j,k) + S(i,j,k);
        // re-use Lu,S as stage1 if we wanted. keeping it simple
        // Need f(T^n) from first stage -> recompute quickly (cheap)
    }
    // Recompute f(T^n)
    exchange_halo(box, T);
    apply_dirichlet_bc(box, T, p.T_env);
    laplace7(T, Lu);
    source_arrhenius(p, T, S);

    for(int k=1;k<=nz;++k) for(int j=1;j<=ny;++j) for(int i=1;i<=nx;++i){
        double f0 = c.alpha*Lu(i,j,k) + S(i,j,k);
        // f(tmp)
        // We already computed laplace(tmp) and source(tmp) into Lu,S above, but we overwrote Lu,S.
        // To keep code short, recompute (small cost).
    }
    // Recompute for tmp again (yes, redundant, but simple)
    exchange_halo(box, tmp);
    apply_dirichlet_bc(box, tmp, p.T_env);
    laplace7(tmp, Lu);
    source_arrhenius(p, tmp, S);

    for(int k=1;k<=nz;++k) for(int j=1;j<=ny;++j) for(int i=1;i<=nx;++i){
        double f0 = 0; // unused after rewrite below
        double ftmp = c.alpha*Lu(i,j,k) + S(i,j,k);

        // Heun: T^{n+1} = T^n + dt/2 * (f(T^n) + f(T^n + dt f(T^n)))
        // quick local Laplacian of original T was lost; do minimal store next time maybe.
        // To avoid more buffers, do a cheap approximation: average with first stage we stored in tmp - T
        // tmp = T + dt f0 -> so f0 ≈ (tmp - T)/dt
        // NOTE: this is exact for explicit Euler stage.

        double f0_est = (tmp(i,j,k) - T(i,j,k)) / p.dt;
        T(i,j,k) = T(i,j,k) + 0.5*p.dt*(f0_est + ftmp);
    }
}

// "IMEX": treat diffusion implicitly by Jacobi iterations on (I - theta*dt*alpha*L) T^{n+1} = T^n + dt*S(T^n)
inline void step_imex_jacobi(const Params& p, const MPIBox& box, const Coeff& c, Field& T, Field& rhs, Field& Lu, Field& S){
    // Build RHS = T^n + dt*S(T^n)
    exchange_halo(box, T);
    apply_dirichlet_bc(box, T, p.T_env);
    source_arrhenius(p, T, S);

    int nx=T.g.nx, ny=T.g.ny, nz=T.g.nz;

    for(int k=1;k<=nz;++k) for(int j=1;j<=ny;++j) for(int i=1;i<=nx;++i){
        rhs(i,j,k) = T(i,j,k) + p.dt * S(i,j,k);
    }

    // Jacobi sweeps
    Field Tnew = T;
    
    const double theta = p.imex_theta;
    const double ax = c.alpha*theta*p.dt/(T.g.dx*T.g.dx);
    const double ay = c.alpha*theta*p.dt/(T.g.dy*T.g.dy);
    const double az = c.alpha*theta*p.dt/(T.g.dz*T.g.dz);
    const double diag = 1.0 + 2.0*(ax+ay+az);

    for(int it=0; it<p.jacobi_iters; ++it){
        exchange_halo(box, T);
        apply_dirichlet_bc(box, T, p.T_env);

        for(int k=1;k<=nz;++k){
            for(int j=1;j<=ny;++j){
                for(int i=1;i<=nx;++i){
                    double nb = ax*(T(i-1,j,k)+T(i+1,j,k)) + ay*(T(i,j-1,k)+T(i,j+1,k)) + az*(T(i,j,k-1)+T(i,j,k+1));
                    Tnew(i,j,k) = (rhs(i,j,k) + nb) / diag;
                }
            }
        }
        
        std::swap(T.data, Tnew.data);
    }
}

}
