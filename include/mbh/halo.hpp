#pragma once
#include "grid.hpp"
#include <mpi.h>
#include <vector>
#include <algorithm>

namespace mbh {

inline void exchange_halo(const MPIBox& box, Field& T){
    const int nx=T.g.nx, ny=T.g.ny, nz=T.g.nz;

    // X faces: size = ny*nz
    std::vector<double> send_xm(ny*nz), send_xp(ny*nz);
    std::vector<double> recv_xm(ny*nz), recv_xp(ny*nz);

    for(int k=1;k<=nz;++k) for(int j=1;j<=ny;++j){
        send_xm[(k-1)*ny + (j-1)] = T(1, j, k);
        send_xp[(k-1)*ny + (j-1)] = T(nx, j, k);
    }

    MPI_Request reqs[12]; int rq=0;

    if(box.nbr_xm!=MPI_PROC_NULL){
        MPI_Irecv(recv_xm.data(), (int)recv_xm.size(), MPI_DOUBLE, box.nbr_xm, 10, box.comm, &reqs[rq++]);
        MPI_Isend(send_xm.data(), (int)send_xm.size(), MPI_DOUBLE, box.nbr_xm, 11, box.comm, &reqs[rq++]);
    }

    if(box.nbr_xp!=MPI_PROC_NULL){
        MPI_Irecv(recv_xp.data(), (int)recv_xp.size(), MPI_DOUBLE, box.nbr_xp, 11, box.comm, &reqs[rq++]);
        MPI_Isend(send_xp.data(), (int)send_xp.size(), MPI_DOUBLE, box.nbr_xp, 10, box.comm, &reqs[rq++]);
    }

    // Y faces: size = nx*nz
    std::vector<double> send_ym(nx*nz), send_yp(nx*nz);
    std::vector<double> recv_ym(nx*nz), recv_yp(nx*nz);

    for(int k=1;k<=nz;++k) for(int i=1;i<=nx;++i){
        send_ym[(k-1)*nx + (i-1)] = T(i, 1, k);
        send_yp[(k-1)*nx + (i-1)] = T(i, ny, k);
    }

    if(box.nbr_ym!=MPI_PROC_NULL){
        MPI_Irecv(recv_ym.data(), (int)recv_ym.size(), MPI_DOUBLE, box.nbr_ym, 20, box.comm, &reqs[rq++]);
        MPI_Isend(send_ym.data(), (int)send_ym.size(), MPI_DOUBLE, box.nbr_ym, 21, box.comm, &reqs[rq++]);
    }

    if(box.nbr_yp!=MPI_PROC_NULL){
        MPI_Irecv(recv_yp.data(), (int)recv_yp.size(), MPI_DOUBLE, box.nbr_yp, 21, box.comm, &reqs[rq++]);
        MPI_Isend(send_yp.data(), (int)send_yp.size(), MPI_DOUBLE, box.nbr_yp, 20, box.comm, &reqs[rq++]);
    }

    // Z faces: size = nx*ny
    std::vector<double> send_zm(nx*ny), send_zp(nx*ny);
    std::vector<double> recv_zm(nx*ny), recv_zp(nx*ny);

    for(int j=1;j<=ny;++j) for(int i=1;i<=nx;++i){
        send_zm[(j-1)*nx + (i-1)] = T(i, j, 1);
        send_zp[(j-1)*nx + (i-1)] = T(i, j, nz);
    }

    if(box.nbr_zm!=MPI_PROC_NULL){
        MPI_Irecv(recv_zm.data(), (int)recv_zm.size(), MPI_DOUBLE, box.nbr_zm, 30, box.comm, &reqs[rq++]);
        MPI_Isend(send_zm.data(), (int)send_zm.size(), MPI_DOUBLE, box.nbr_zm, 31, box.comm, &reqs[rq++]);
    }

    if(box.nbr_zp!=MPI_PROC_NULL){
        MPI_Irecv(recv_zp.data(), (int)recv_zp.size(), MPI_DOUBLE, box.nbr_zp, 31, box.comm, &reqs[rq++]);
        MPI_Isend(send_zp.data(), (int)send_zp.size(), MPI_DOUBLE, box.nbr_zp, 30, box.comm, &reqs[rq++]);
    }

    if(rq>0) MPI_Waitall(rq, reqs, MPI_STATUSES_IGNORE);

    // unpack into ghost cells
    if(box.nbr_xm!=MPI_PROC_NULL){
        for(int k=1;k<=nz;++k) for(int j=1;j<=ny;++j){
            T(0,j,k) = recv_xm[(k-1)*ny + (j-1)];
        }
    }
    if(box.nbr_xp!=MPI_PROC_NULL){
        for(int k=1;k<=nz;++k) for(int j=1;j<=ny;++j){
            T(nx+1,j,k) = recv_xp[(k-1)*ny + (j-1)];
        }
    }
    if(box.nbr_ym!=MPI_PROC_NULL){
        for(int k=1;k<=nz;++k) for(int i=1;i<=nx;++i){
            T(i,0,k) = recv_ym[(k-1)*nx + (i-1)];
        }
    }
    if(box.nbr_yp!=MPI_PROC_NULL){
        for(int k=1;k<=nz;++k) for(int i=1;i<=nx;++i){
            T(i,ny+1,k) = recv_yp[(k-1)*nx + (i-1)];
        }
    }
    if(box.nbr_zm!=MPI_PROC_NULL){
        for(int j=1;j<=ny;++j) for(int i=1;i<=nx;++i){
            T(i,j,0) = recv_zm[(j-1)*nx + (i-1)];
        }
    }
    if(box.nbr_zp!=MPI_PROC_NULL){
        for(int j=1;j<=ny;++j) for(int i=1;i<=nx;++i){
            T(i,j,nz+1) = recv_zp[(j-1)*nx + (i-1)];
        }
    }
}

inline void apply_dirichlet_bc(const MPIBox& box, Field& T, double Tenv){
    // If boundary (no neighbor), set ghost to Tenv (simple Dirichlet to env)
    int nx=T.g.nx, ny=T.g.ny, nz=T.g.nz;
    
    if(box.nbr_xm==MPI_PROC_NULL){ for(int k=1;k<=nz;++k) for(int j=1;j<=ny;++j) T(0,j,k)=Tenv; }
    if(box.nbr_xp==MPI_PROC_NULL){ for(int k=1;k<=nz;++k) for(int j=1;j<=ny;++j) T(nx+1,j,k)=Tenv; }
    if(box.nbr_ym==MPI_PROC_NULL){ for(int k=1;k<=nz;++k) for(int i=1;i<=nx;++i) T(i,0,k)=Tenv; }
    if(box.nbr_yp==MPI_PROC_NULL){ for(int k=1;k<=nz;++k) for(int i=1;i<=nx;++i) T(i,ny+1,k)=Tenv; }
    if(box.nbr_zm==MPI_PROC_NULL){ for(int j=1;j<=ny;++j) for(int i=1;i<=nx;++i) T(i,j,0)=Tenv; }
    if(box.nbr_zp==MPI_PROC_NULL){ for(int j=1;j<=ny;++j) for(int i=1;i<=nx;++i) T(i,j,nz+1)=Tenv; }
}

}