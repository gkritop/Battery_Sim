#pragma once
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <filesystem>
#include <mpi.h>
#include "grid.hpp"

namespace mbh {

inline void write_raw_slice_xy(const std::string& fname, const Field& T, int kmid){
    std::ofstream os(fname);
    int nx=T.g.nx, ny=T.g.ny;

    for(int j=1;j<=ny;++j){
        for(int i=1;i<=nx;++i){
            os<< std::setprecision(6)<<T(i,j,kmid);
            
            if(i<nx) os<<",";
        }
        os<<"\n";
    }
}

inline void write_rankbin(const std::string& fname, const Field& T){
    std::ofstream os(fname, std::ios::binary);
    int nx=T.g.nx, ny=T.g.ny, nz=T.g.nz;

    os.write((char*)&nx, sizeof(int));
    os.write((char*)&ny, sizeof(int));
    os.write((char*)&nz, sizeof(int));

    for(int k=1;k<=nz;++k) for(int j=1;j<=ny;++j) for(int i=1;i<=nx;++i){
        double v=T(i,j,k);
        os.write((char*)&v, sizeof(double));
    }
}

inline void ensure_dir(const std::string& d){
    std::error_code ec; std::filesystem::create_directories(d, ec);
    (void)ec;
}

}