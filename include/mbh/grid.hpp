#pragma once
#include <vector>
#include <cmath>
#include <cassert>
#include <mpi.h>

namespace mbh {

struct Grid {
    // local sizes (without ghosts)
    int nx=0, ny=0, nz=0;
    // origin index (global start indices)
    int ox=0, oy=0, oz=0;
    // spacings
    double dx=0, dy=0, dz=0;
    // global sizes
    int NX=0, NY=0, NZ=0;
};

struct MPIBox {
    MPI_Comm comm = MPI_COMM_NULL;

    int rank=0, size=1;
    int dims[3]{0,0,0};
    int periods[3]{0,0,0};
    int coords[3]{0,0,0};

    int nbr_xm=-1, nbr_xp=-1, nbr_ym=-1, nbr_yp=-1, nbr_zm=-1, nbr_zp=-1;
};

inline void create_cartesian(int px,int py,int pz, MPIBox& box){
    MPI_Comm_size(MPI_COMM_WORLD, &box.size);

    if(px*py*pz==0){
        int dims[3] = {0,0,0};
        MPI_Dims_create(box.size, 3, dims);

        box.dims[0]=dims[0]; box.dims[1]=dims[1]; box.dims[2]=dims[2];
    }
    else{
        box.dims[0]=px; box.dims[1]=py; box.dims[2]=pz;
    }
    box.periods[0]=box.periods[1]=box.periods[2]=0;
    MPI_Cart_create(MPI_COMM_WORLD, 3, box.dims, box.periods, 0, &box.comm);
    MPI_Comm_rank(box.comm, &box.rank);
    MPI_Cart_coords(box.comm, box.rank, 3, box.coords);

    // neighbors
    MPI_Cart_shift(box.comm, 0, 1, &box.nbr_xm, &box.nbr_xp);
    MPI_Cart_shift(box.comm, 1, 1, &box.nbr_ym, &box.nbr_yp);
    MPI_Cart_shift(box.comm, 2, 1, &box.nbr_zm, &box.nbr_zp);
}

inline void split_sizes(int N, int P, int c, int& nloc, int& off){
    int base = N / P;
    int rem  = N % P;

    if(c < rem){ nloc = base+1; off = c*(base+1); }
    else { nloc = base; off = rem*(base+1) + (c-rem)*base; }
}

struct Field {
    Grid g;
    std::vector<double> data; // includes ghosts (1-layer)

    inline int idx(int i,int j,int k) const {
        int nxg = g.nx + 2;
        int nyg = g.ny + 2;

        return (k*nyg + j)*nxg + i;
    }

    void resize(const Grid& grid){
        g = grid;
        data.assign((g.nx+2)*(g.ny+2)*(g.nz+2), 0.0);
    }

    inline double& operator()(int i,int j,int k){ return data[idx(i,j,k)]; }
    
    inline double  operator()(int i,int j,int k) const { return data[idx(i,j,k)]; }
};

}
