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
#include <string>

using namespace std;
class GetFreeEnergyProfile{
    
public:
    GetFreeEnergyProfile(FILE* input, FILE* output, int min_N, int max_N, double d_N, double Rho, double Mu, double Sigma, double T, double theta_cg, double theta_cb);
    ~GetFreeEnergyProfile();
    
    virtual void read(){};
    virtual void printToFile(){};
    
    
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

    void add_to_format_string(int repeats, string s);
    void setup_NArray_quant(int);

    //Template functions related to scaning
    template <class T> int scan(T* output, string format)
    {
        int fscanf_result = fscanf(input_file, format.c_str(), output);
        return fscanf_result;
    }
    
    template <class T> int scan(std::vector<T>* output, string individual_format)
    {
        int total_matches = 0;
        int size = (int)output->size();
        int fscanf_result;
        for(int i=0 ; i<size; i++)
        {
            fscanf_result = fscanf(input_file, individual_format.c_str(), &((*output)[i]));
            if(fscanf_result == EOF)
            {
                return EOF;
            }
            else {total_matches += fscanf_result;}
        }
        return total_matches;
    }
    
    template <class T> int scan_end(std::vector<T>* output, string individual_format, string end_format)
    {
        int total_matches = 0;
        int size = (int)output->size();
        int fscanf_result;
        for(int i=0 ; i<size; i++)
        {
            if(i == size-1)
            {
                fscanf_result = fscanf(input_file, end_format.c_str(), &((*output)[i]));
            }
            else
            {
                fscanf_result = fscanf(input_file, individual_format.c_str(), &((*output)[i]));
            }
    
            if(fscanf_result == EOF)
            {
                return EOF;
            }
            else {total_matches += fscanf_result;}
        }
        return total_matches;
        
    }
};

//////////////////////////////////////// GetSphericalCapsFreeEnergy //////////////////////////////////////

class GetSphericalCapsFreeEnergy : public GetFreeEnergyProfile {
   
public:
    GetSphericalCapsFreeEnergy(FILE* input, FILE* output, int min_N, int max_N, double d_N, double Rho, double Mu, double Sigma, double T, double theta_cg, double theta_cb, int patches);
    ~GetSphericalCapsFreeEnergy();
    void read();
    void printToFile();
    
    void convert_quants_to_metre_cube();

private:
    double G, N, Volume, SA;
    std::vector<double> projected_SA ; //Size of number of caps
    std::vector<int> clstr_centre_location_modifier;
    std::vector<double> Radii;
    
    int n_patches, n_unique_patches ; 
    //Scan quantity from the file
    void GetMinimumPerN();
    bool isSingleCap(std::vector<double>);
    void addToNQuants(int);
    int scan_manager_spherical_caps();
    
    
    struct MinQuants
    {
        double G, N, Volume, SA;
        std::vector<double> projected_SA ; //Size of number of caps
        std::vector<int> clstr_centre_location_modifier;
        std::vector<double> Radii;
    }min_quants;
    void setup_min_quants();
    void update_min_quants_per_N(int);
    void print_min_quants_per_N(int);

};


//class GetSpherocylinderFreeEnergy : public GetFreeEnergyProfile {
//
//public:
//    GetSpherocylinderFreeEnergy();
//    ~GetSpherocylinderFreeEnergy();
//    void read();
//
//private:
//};

#endif /* Get_FreeEnergyProfile_hpp */

