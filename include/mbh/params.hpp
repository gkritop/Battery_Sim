#pragma once
#include <string>
#include <vector>
#include <cstdlib>
#include <iostream>

namespace mbh {

struct Arrhenius {
    double A  = 1.0e7;     // pre-exponential (made smaller to be stable for explicit)
    double Ea = 8.0e4;     // activation energy [J/mol] (not scaled by R here, we keep simple)
    double H  = 3.0e5;     // heat release per reaction "unit"
};

struct Params {
    // domain
    int nx = 96, ny = 96, nz = 96;
    double Lx = 0.064, Ly = 0.064, Lz = 0.064; // 64 mm cube

    // material
    double rho = 2500.0;
    double cp  = 1000.0;
    double k   = 1.0; // W/mK

    // time
    double dt = 5e-5;
    double t_final = 0.01;
    std::string scheme = "rk2"; // or "imex"

    // initial BCs
    double T0 = 300.0;
    double T_env = 300.0;
    double hotspot_T = 380.0;
    double hotspot_radius = 0.004;
    double hotspot_x = 0.5, hotspot_y = 0.5, hotspot_z = 0.5;

    // source
    bool enable_source = true;
    Arrhenius arr;

    // output
    std::string outdir = "out";
    int output_every = 20;

    // MPI decomposition (0 means cartesian)

    int px = 0, py = 0, pz = 0;

    // imex "solver" (very light-weight)
    int jacobi_iters = 20;
    double imex_theta = 1.0; // backward Euler-ish

    // crude CLI parsing: --key=value (strings) or --key value
    void parse_cli(int argc, char** argv) {
        auto setKV = [&](const std::string& key, const std::string& val){
            if(key=="nx") nx = std::atoi(val.c_str());
            else if(key=="ny") ny = std::atoi(val.c_str());
            else if(key=="nz") nz = std::atoi(val.c_str());
            else if(key=="Lx") Lx = std::atof(val.c_str());
            else if(key=="Ly") Ly = std::atof(val.c_str());
            else if(key=="Lz") Lz = std::atof(val.c_str());
            else if(key=="rho") rho = std::atof(val.c_str());
            else if(key=="cp") cp = std::atof(val.c_str());
            else if(key=="k") k = std::atof(val.c_str());
            else if(key=="dt") dt = std::atof(val.c_str());
            else if(key=="t_final") t_final = std::atof(val.c_str());
            else if(key=="scheme") scheme = val;
            else if(key=="T0") T0 = std::atof(val.c_str());
            else if(key=="T_env") T_env = std::atof(val.c_str());
            else if(key=="hotspot_T") hotspot_T = std::atof(val.c_str());
            else if(key=="hotspot_radius") hotspot_radius = std::atof(val.c_str());
            else if(key=="hotspot_x") hotspot_x = std::atof(val.c_str());
            else if(key=="hotspot_y") hotspot_y = std::atof(val.c_str());
            else if(key=="hotspot_z") hotspot_z = std::atof(val.c_str());
            else if(key=="outdir") outdir = val;
            else if(key=="output_every") output_every = std::atoi(val.c_str());
            else if(key=="px") px = std::atoi(val.c_str());
            else if(key=="py") py = std::atoi(val.c_str());
            else if(key=="pz") pz = std::atoi(val.c_str());
            else if(key=="jacobi_iters") jacobi_iters = std::atoi(val.c_str());
            else if(key=="imex_theta") imex_theta = std::atof(val.c_str());
            else if(key=="A") arr.A = std::atof(val.c_str());
            else if(key=="Ea") arr.Ea = std::atof(val.c_str());
            else if(key=="H") arr.H = std::atof(val.c_str());
        };

        for(int i=1;i<argc;++i){
            std::string a = argv[i];
            if(a.rfind("--",0)==0){
                auto eq = a.find('=');
                if(eq!=std::string::npos){
                    setKV(a.substr(2,eq-2), a.substr(eq+1));
                }else{
                    std::string key = a.substr(2);
                    if(i+1<argc){
                        std::string val = argv[i+1];
                        if(val.rfind("--",0)!=0){
                            setKV(key, val); ++i;
                        }else{
                            setKV(key, "1"); // flag
                        }
                    }else{
                        setKV(key, "1");
                    }
                }
            }
        }
    }

    void print_short(int rank) const {
        if(rank==0){
            std::cout<<"[params] nx,ny,nz="<<nx<<","<<ny<<","<<nz
                     <<"L="<<Lx<<","<<Ly<<","<<Lz
                     <<"dt="<<dt<<"  t_final="<<t_final
                     <<"scheme="<<scheme
                     <<"outdir="<<outdir<<"\n";
        }
    }
};

}
