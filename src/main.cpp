#include <mpi.h>
#include <cmath>
#include <iostream>
#include <sstream>
#include <filesystem>

#include "mbh/params.hpp"
#include "mbh/grid.hpp"
#include "mbh/halo.hpp"
#include "mbh/laplacian.hpp"
#include "mbh/source.hpp"
#include "mbh/io.hpp"
#include "mbh/time_integrators.hpp"

using namespace mbh;

int main(int argc, char** argv){
    MPI_Init(&argc,&argv);
    MPIBox box;

    MPI_Comm_rank(MPI_COMM_WORLD, &box.rank);
    MPI_Comm_size(MPI_COMM_WORLD, &box.size);

    Params p; p.parse_cli(argc, argv);

    if(box.rank==0){
        std::cout<<"battery_heat\n";
    }

    create_cartesian(p.px, p.py, p.pz, box);

    // global grid
    Grid g; g.NX=p.nx; g.NY=p.ny; g.NZ=p.nz;
    // split
    int nx,ny,nz,ox,oy,oz;

    split_sizes(g.NX, box.dims[0], box.coords[0], nx, ox);
    split_sizes(g.NY, box.dims[1], box.coords[1], ny, oy);
    split_sizes(g.NZ, box.dims[2], box.coords[2], nz, oz);

    g.nx=nx; g.ny=ny; g.nz=nz; g.ox=ox; g.oy=oy; g.oz=oz;
    g.dx = p.Lx / p.nx;
    g.dy = p.Ly / p.ny;
    g.dz = p.Lz / p.nz;

    Field T, tmp, Lu, S, rhs;
    T.resize(g); tmp.resize(g); Lu.resize(g); S.resize(g); rhs.resize(g);

    // init: uniform T0 with one spherical hotspot
    for(int k=0;k<=g.nz+1;++k){
        double z = (g.oz + (k-1))*g.dz;
        for(int j=0;j<=g.ny+1;++j){
            double y = (g.oy + (j-1))*g.dy;
            for(int i=0;i<=g.nx+1;++i){
                double x = (g.ox + (i-1))*g.dx;

                T(i,j,k) = p.T0;

                double cx = p.hotspot_x*p.Lx, cy = p.hotspot_y*p.Ly, cz = p.hotspot_z*p.Lz;
                double r2 = (x-cx)*(x-cx) + (y-cy)*(y-cy) + (z-cz)*(z-cz);
                
                if(r2 < p.hotspot_radius*p.hotspot_radius) T(i,j,k) = p.hotspot_T;
            }
        }
    }
    apply_dirichlet_bc(box, T, p.T_env);

    Coeff c = make_coeff(p);
    ensure_dir(p.outdir);

    double t=0.0;
    int step=0;
    double t0 = MPI_Wtime();

    p.print_short(box.rank);

    while(t < p.t_final + 1e-12){
        if(p.scheme=="rk2"){
            step_rk2(p, box, c, T, tmp, Lu, S);
        }else{
            step_imex_jacobi(p, box, c, T, rhs, Lu, S);
        }
        t += p.dt;
        ++step;

        if(step % p.output_every == 0 || t >= p.t_final){
            std::stringstream ss1, ss2;

            ss1<<p.outdir<<"/rank"<<box.rank<<"_step"<<step<<".bin";
            ss2<<p.outdir<<"/rank"<<box.rank<<"_xy_mid"<<step<<".csv";
            
            write_rankbin(ss1.str(), T);
            int kmid = T.g.nz/2;

            write_raw_slice_xy(ss2.str(), T, kmid);
            if(box.rank==0){
                std::cout<<"[io] wrote step "<<step<<" at t="<<t<<"\n";
            }
        }
    }

    double t1 = MPI_Wtime();
    double loc = t1 - t0;
    double max_t=0.0;
    
    MPI_Reduce(&loc, &max_t, 1, MPI_DOUBLE, MPI_MAX, 0, box.comm);
    
    if(box.rank==0){
        std::cout<<"[timing] wall max = "<<max_t<<" s\n";
        std::cout<<"Done.\n";
    }

    MPI_Finalize();
    return 0;
}
