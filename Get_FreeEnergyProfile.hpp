//
//  Get_FreeEnergyProfile_header.hpp
//  Clusters
//
//  Created by Sarwar Hussain on 3/20/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#ifndef Get_FreeEnergyProfile_hpp
#define Get_FreeEnergyProfile   _hpp

#include <stdio.h>
#include "FreeEnergy.hpp"
#include "miscformulas.hpp"

class GetFreeEnergyProfile{
    
public:
    GetFreeEnergyProfile(FILE* input, FILE* output, int min_N, int max_N, double d_N, double Rho, double Mu, double Sigma, double T, double theta_cg, double theta_cb);
    ~GetFreeEnergyProfile();
    void ReadProfileSphericalCaps();
    void ReadProfileSpherocylinder();
    
protected:
    const char* format;
    FILE* input_file;
    FILE* output_file;
    double Rho;         //Moles per m3
    double Mu;
    double Sigma;       //Liquid-crystal surface tension
    double T;           //Kelvin
    double theta_cg, theta_cb;
    int Nmin, Nmax;
    double dN;
    int len_N;
    std::vector<double> NArray_Gmin;    //Vector of minimum G value for each integer N value
    std::vector<double> NArray_confs;   //Number of distinct points (i.e. clusters) for each integer value of N
    std::vector<std::vector<double> > NArray_quant; //[N, Gmin, Rg, Rb, Rg_secondary, db, dg_secondary, V, SA]. Quantities such as G, Rg, Rb etc. correponding to the minimum G cluster for each integer value of N
};

#endif /* Get_FreeEnergyProfile_hpp */
