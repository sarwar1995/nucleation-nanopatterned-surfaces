//
//  Get_FreeEnergyProfile.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 2/16/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#include <stdio.h>
#include "FreeEnergy.hpp"
#include "miscformulas.hpp"
#include "Get_FreeEnergyProfile.hpp"

GetFreeEnergyProfile::GetFreeEnergyProfile(FILE* input, FILE* output, int min_N, int max_N, double d_N, double rho, double mu, double sigma, double t, double theta_g, double theta_b)
: Nmin(min_N), Nmax (max_N), dN(d_N), input_file(input), output_file(output), Rho(rho), Mu(mu), Sigma(sigma), T(t), theta_cg(theta_g), theta_cb(theta_b)
{
    len_N  = (int) ((Nmax-Nmin)/dN) + 1;
    NArray_Gmin.resize(len_N, 0.0);    //Vector of minimum G value for each integer N value
    NArray_confs.resize (len_N , 0.0);   //Number of distinct points (i.e. clusters) for each integer value of N
    NArray_quant.resize (len_N); //[N, Gmin, Rg, Rb, Rg_secondary, db, dg_secondary, V, SA]. Quantities such as G, Rg, Rb etc. correponding to the minimum G cluster for each integer value of N
}

GetFreeEnergyProfile::~GetFreeEnergyProfile()
{
    
}


void GetFreeEnergyProfile::setup_NArray_quant(int length)
{
    
    for (int i=0 ; i<len_N; i++)
    {
        NArray_quant[i] = std::vector<double> (length,0.0);
        if(i==0)
        {
            NArray_quant[i][0] = Nmin;
        }
        else
        {
            NArray_quant[i][0]=NArray_quant[i-1][0]+dN ;
        }
    }
}

void GetFreeEnergyProfile::add_to_format_string(int repeats, string s)
{
    string format_str;
    for(int i=0; i<repeats; i++)
        format_str = format + s;
    format = format_str.c_str();
    
}

void GetFreeEnergyProfile::ReadProfileSphericalCaps()
{
    double N, Volume, SA;
    double projected_SA_single_cap ;
    int db, dg_secondary;
    
    std::vector<double> projected_SA(5, 0.0);   //This is obviously only for the 5 patch case with 3 pseudo-indipendent radii.
    std::vector<double> Radii (3, 0.0);
    /* projected_SA = [centre_good, bad_left, bad_right, good_left, good_right]
     Radii        = [centre_good, bad, secondary_good]*/
    double G;
    
    //NArray_quant will have 9 values only for the all spherical caps case
    setup_NArray_quant(9);
    
    format = "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf%d\t%d\n" ;
    if(input_file == NULL)
    {
        printf("Error opening input file\n");
        exit(1);
    }
    int fscanf_result = fscanf(input_file, format, &N, &Volume, &SA, &projected_SA[0], &projected_SA[1], &projected_SA[2], &projected_SA[3], &projected_SA[4], &Radii[0], &Radii[1], &Radii[2], &db, &dg_secondary);
    
    while(fscanf_result != EOF)
    {
        Volume = Volume * 1e-30;    //This converts the volume from Ang^3 to m^3
        SA = SA* 1e-20;
        for(size_t y=0; y<projected_SA.size(); y++)
        {
            projected_SA[y] = projected_SA[y]*1e-20 ;
        }
        
        if(Radii[1]==0.0 && Radii[2] == 0.0) //Just a single good cap case
        {
            projected_SA_single_cap = projected_SA[0];
            if(projected_SA[1] != 0.0 || projected_SA[2] != 0.0)
            {
                printf("projected SA of bad and good secondary are not zero for Rb=Rg_second=0.0, projected_SA[1]=%10.10f\t projected_SA[2]=%10.10f\n", projected_SA[1], projected_SA[2]);
                abort();
            }
            G = free_energy_singlecap(Rho, Mu, Sigma, Volume, SA, projected_SA_single_cap, theta_cg, T);
            add_to_N(N,  G, Radii[0], Radii[1], Radii[2], db, dg_secondary, Volume, SA, dN, Nmin, len_N, NArray_Gmin, NArray_confs,  NArray_quant);
        }
        
        else
        {
            G = free_energy(Rho, Mu, Sigma, Volume, SA, projected_SA, theta_cg, theta_cb, T);
            add_to_N(N,  G, Radii[0], Radii[1], Radii[2], db, dg_secondary, Volume, SA, dN, Nmin, len_N, NArray_Gmin, NArray_confs,  NArray_quant);
        }
        fscanf_result = fscanf(input_file, format, &N, &Volume, &SA, &projected_SA[0], &projected_SA[1], &projected_SA[2], &projected_SA[3], &projected_SA[4], &Radii[0], &Radii[1], &Radii[2], &db, &dg_secondary);
    }
    print_NGDataFile (NArray_quant, output_file);
    fclose(output_file);
    fclose(input_file);
}

void GetFreeEnergyProfile::ReadProfileSpherocylinder()
{
    //Input file should be in the format:
    // Rg, cyl_length, chord_length, Rb, dB, N, Volume, SA, ProjSA[0], ProjSA[1], ProjSA[2]
    
    format = "%lf\t%lf\t%lf\t%lf\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n";
    double Rg, cyl_length, chord_length, Rb, N, Volume, SA;
    std::vector<double> ProjSA (3,0.0);
    double G;
    int dB;
    //NArray_quant will have 11 values only for the spherocylinder case
    setup_NArray_quant(11);
    
    if(input_file == NULL)
    {
        printf("Error opening input file\n");
        exit(1);
    }
    
    
    int fscanf_result = fscanf(input_file, format, &Rg, &cyl_length, &chord_length, &Rb, &dB, &N, &Volume, &SA, &ProjSA[0], &ProjSA[1], &ProjSA[2]);
    
    
    while(fscanf_result != EOF)
    {
        Volume = Volume * 1e-30;    //This converts the volume from Ang^3 to m^3
        SA = SA * 1e-20;
        
        for(int j=0; j< (int)ProjSA.size(); j++){
            ProjSA[j] = ProjSA[j] * 1e-20;
            
        }
        G = free_energy_spherocylinder (Rho, Mu, Sigma, Volume, SA, ProjSA, theta_cg, theta_cb, T);
//        G = free_energy_singlecap(Rho, Mu, Sigma, Volume, SA, ProjSA, theta_cg, T);

        add_to_N_spherocylinder (N, G, Rg, cyl_length, chord_length, Rb, dB, Volume,  SA, ProjSA,  dN,  Nmin,  len_N, NArray_Gmin, NArray_confs, NArray_quant);
        
        fscanf_result = fscanf(input_file, format, &Rg, &cyl_length, &chord_length, &Rb, &dB, &N, &Volume, &SA, &ProjSA[0], &ProjSA[1], &ProjSA[2]);
    }
    print_NGDataFile_spherocylinder (NArray_quant, output_file);
    fclose(output_file);
    fclose(input_file);
}

//////////////////////////////////////// GetSphericalCapsFreeEnergy //////////////////////////////////////

GetSphericalCapsFreeEnergy::GetSphericalCapsFreeEnergy(FILE* input, FILE* output, int min_N, int max_N, double d_N, double Rho, double Mu, double Sigma, double T, double theta_cg, double theta_cb, int patches):GetFreeEnergyProfile( input,  output,  min_N,  max_N,  d_N,  Rho,  Mu,  Sigma,  T,  theta_cg,  theta_cb), n_patches(patches)
{
    if(n_patches < 3) {printf("Constructing GetSphericalCapsFreeEnergy, n_patches is less than 3, error\n"); abort();}
    n_unique_patches = ((n_patches-1)/2) + 1 ;
    
    int NArray_quant_size = 4 + n_unique_patches + n_patches + n_unique_patches;
    setup_NArray_quant(NArray_quant_size);
    projected_SA.resize(n_patches, 0.0);
    clstr_centre_location_modifier.resize(n_unique_patches, 0.0);
    Radii.resize(n_unique_patches, 0.0);
    setup_min_quants();
}

GetSphericalCapsFreeEnergy::~GetSphericalCapsFreeEnergy(){}

void GetSphericalCapsFreeEnergy::setup_min_quants()
{
    min_quants.projected_SA.resize(n_patches, 0.0);
    min_quants.clstr_centre_location_modifier.resize(n_unique_patches, 0.0);
    min_quants.Radii.resize(n_unique_patches, 0.0);
}


void GetSphericalCapsFreeEnergy::read()
{
    // G + N + V + SA = 4, [R + clstr_centre_location_modifier + proj_SA] = 3*n_unique_patches
    int fscanf_result = scan_manager_spherical_caps();
    while(fscanf_result != EOF)
    {
        convert_quants_to_metre_cube();
        if (isSingleCap(Radii))
        {
            G = free_energy_singlecap(Rho, Mu, Sigma, Volume, SA, projected_SA[0], theta_cg, T);
        }
        else
        {
            G = free_energy(Rho, Mu, Sigma, Volume, SA, projected_SA, theta_cg, theta_cb, T);
        }
        GetMinimumPerN();
        fscanf_result = scan_manager_spherical_caps();
    }
    printToFile();
    fclose(output_file);
    fclose(input_file);
}

//double* N, double* Volume, double* SA,  std::vector<double>* projected_SAs, std::vector<int>* clstr_centre_location_modifier, std::vector<double>* Radii
int GetSphericalCapsFreeEnergy::scan_manager_spherical_caps()
{
    if(input_file == NULL)
    {
        printf("Error opening input file\n");
        exit(1);
    }
    int scan_N, scan_V, scan_SA, scan_proj, scan_clstr_centre, scan_radii;
    string double_string ("%lf\t");
    string int_string ("%d\t");
    string double_end_string ("%lf\n");
    scan_N = scan<double>(&N, double_string);
    scan_V = scan<double>(&Volume, double_string);
    scan_SA = scan<double>(&SA, double_string);
    scan_proj = scan<double>(&projected_SA, double_string);
    scan_clstr_centre  = scan<int>(&clstr_centre_location_modifier, int_string);
    scan_radii = scan_end<double>(&Radii, double_string, double_end_string);
//    printf("scan_N=%d scan_V=%d scan_SA=%d scan_proj=%d scan_clstr_centre=%d scan_radii=%d\n", scan_N, scan_V, scan_SA, scan_proj, scan_clstr_centre, scan_radii);
//
//    printf("N=%10.5f V=%10.10f SA=%10.10f proj_SA=[%10.10f %10.10f %10.10f] clstr_centre=[%d %d] radii=[%10.10f %10.10f]\n",
//           N, Volume, SA, projected_SA[0], projected_SA[1], projected_SA[2], clstr_centre_location_modifier[0], clstr_centre_location_modifier[1], Radii[0], Radii[1]);
    
    
    if (scan_N == EOF || scan_V == EOF || scan_SA == EOF || scan_proj == EOF || scan_clstr_centre == EOF || scan_radii == EOF){return EOF;}
    else return (scan_N + scan_V + scan_SA + scan_proj + scan_clstr_centre + scan_radii);

}

void GetSphericalCapsFreeEnergy::convert_quants_to_metre_cube()
{
    Volume = Volume * 1e-30;    //This converts the volume from Ang^3 to m^3
    SA = SA* 1e-20;
    for(size_t y=0; y<projected_SA.size(); y++)
    {
        projected_SA[y] = projected_SA[y]*1e-20 ;
    }
}

bool GetSphericalCapsFreeEnergy::isSingleCap(std::vector<double> radii)
{
    if(radii.empty())
    {
        printf("radii vector is empty\n"); abort();
    }
    else
    {
        if(radii[0] != 0.0)
        {
            for(int i=1; i<radii.size(); i++)
            {
                if(radii[i] != 0.0) {return false;}
            }
            return true;
        }
        else
        {
            printf("No first cap as well\n");
            return true;
        }
    }
    
}

void GetSphericalCapsFreeEnergy::GetMinimumPerN ()
{
    double rNhere = round_nearN (N , dN);   //This rounds to nearest N.
    int indN = (int) ((rNhere-Nmin)/dN);
    if(indN >= len_N){printf("Number of particles are not enough rNhere = %10.5f\t N=%10.5f\n", rNhere, N); abort();}
    
    double minNi = NArray_Gmin[indN]; //This is 0 when no free energy has been added
    if(G <= minNi || NArray_confs[indN] == 0.0)
    {
        addToNQuants(indN);
    }
    else
    {
        NArray_quant[indN][1] = minNi;
    }
    NArray_confs[indN] = NArray_confs[indN] + 1;
}

void GetSphericalCapsFreeEnergy::addToNQuants(int indN)
{
    NArray_Gmin[indN] = G ;
    NArray_quant[indN][1] = G;
    NArray_quant[indN][2] = Volume ;
    NArray_quant[indN][3] = SA ;
    
    for(int i=0; i<n_patches; i++)
    {
        int ind = i+4;
        NArray_quant[indN][ind] = projected_SA[i];
    }
    for(int i=0; i<n_unique_patches; i++)
    {
        int ind = i+4+n_patches;
        NArray_quant[indN][ind] = (double)clstr_centre_location_modifier[i];
    }
    for(int i=0; i<n_unique_patches; i++)
    {
        int ind = i+4+n_patches+n_unique_patches;
        NArray_quant[indN][ind] = Radii[i];
    }
}

void GetSphericalCapsFreeEnergy::printToFile()
{
    for(int i = 0; i<len_N ; i++)
    {
        update_min_quants_per_N(i);
        print_min_quants_per_N(i);
    }
}

void GetSphericalCapsFreeEnergy::update_min_quants_per_N(int i)
{
    min_quants.N = NArray_quant[i][0];
    min_quants.G = NArray_quant[i][1];
    min_quants.Volume = NArray_quant[i][2];
    min_quants.SA = NArray_quant[i][3];
//    printf("min N, G, Volume, SA = %10.5f\t%10.10f\t%10.10f\t%10.10f\n", min_quants.N, min_quants.G, min_quants.Volume, min_quants.SA);
    for(int j=0; j<n_patches; j++)
    {
        int ind = j+4;
        min_quants.projected_SA[j] = NArray_quant[i][ind];
//        printf("ind = %dmin projected_SA[%d]=%10.10f\t",ind, j, min_quants.projected_SA[j]);
    }
    for(int j=0; j<n_unique_patches; j++)
    {
        int ind = j+4+n_patches;
        min_quants.clstr_centre_location_modifier[j] = (int)NArray_quant[i][ind];
//        printf("ind=%d min clstr_centre_location_modifier[%d]=%d\t",ind, j, min_quants.clstr_centre_location_modifier[j]);
    }
    for(int j=0; j<n_unique_patches; j++)
    {
        int ind = j+4+n_patches+n_unique_patches;
        min_quants.Radii[j] = NArray_quant[i][ind];
//        printf("ind=%d min Radii[%d]=%10.10f\t", ind, j, min_quants.Radii[j]);
    }
}

void GetSphericalCapsFreeEnergy::print_min_quants_per_N(int i)
{
    if(output_file == NULL)
    {
        printf("Error opening output file\n");
        exit(1);
    }
    //[N, Gmin (/kbT), V, SA, proj_SA, clstr_centre_modifier, Radii]
    fprintf(output_file,"%10.5f\t%10.10f\t%10.10f\t%10.10f\t" , min_quants.N , min_quants.G , min_quants.Volume * 1e30 , min_quants.SA * 1e20);
    
    for(int j=0; j<n_patches; j++)
    {
        fprintf(output_file,"%10.10f\t", min_quants.projected_SA[j] * 1e20);
    }
    for(int j=0; j<n_unique_patches; j++)
    {
        fprintf(output_file,"%d\t", min_quants.clstr_centre_location_modifier[j]);
    }
    for(int j=0; j<n_unique_patches; j++)
    {
        if(j == n_unique_patches-1)
        {
            fprintf(output_file,"%10.10f\n", min_quants.Radii[j]);
        }
        else
        {fprintf(output_file,"%10.10f\t", min_quants.Radii[j]);}
    }
}


int main(int argc, const char * argv[])
{
    /* Energy related variables*/
    double Rho;         //Moles per m3
    double Mu;
    double Sigma;       //Liquid-crystal surface tension
    double T;           //Kelvin
    double theta_cg, theta_cb;
    
    int Nmin, Nmax;
    double dN;
    int isSpherocylinderFile;
    FILE* input_file;
    FILE* output_file;
    
    int n_patches;
    
    Nmin        = atoi(argv[1]);
    dN          = atof(argv[2]);
    Nmax        = atoi(argv[3]);
    Rho         = atof(argv[4]);   //Moles per m3
    Mu          = atof(argv[5]);
    Sigma       = atof(argv[6]);   //Liquid-crystal surface tension
    T           = atof(argv[7]);  //Kelvin
    theta_cg    = atof(argv[8]);
    theta_cb    = atof(argv[9]);
    input_file  = fopen(argv[10], "r");
    output_file = fopen(argv[11], "w");
    isSpherocylinderFile = atof(argv[12]);
    n_patches   = atoi(argv[13]);
    
    GetFreeEnergyProfile get_free_energy_profile (input_file, output_file, Nmin, Nmax, dN, Rho, Mu, Sigma, T, theta_cg, theta_cb);
    if (isSpherocylinderFile==1)
    {
        get_free_energy_profile.ReadProfileSpherocylinder();
    }
    else
    {
        GetSphericalCapsFreeEnergy get_spherical_caps_free_energy_profile (input_file, output_file, Nmin, Nmax, dN, Rho, Mu, Sigma, T, theta_cg, theta_cb, n_patches);
        get_spherical_caps_free_energy_profile.read();
    }
    
    
    
    return 0;
}


//printf("%10.3f\t%10.20f\t%10.20f\t%10.20f\t%10.20f\t%10.20f\t%10.20f\t%10.20f\t%10.20f\t%10.20f\t%10.20f\t%d\t%d\n", N, Volume, SA, projected_SA[0], projected_SA[1], projected_SA[2], projected_SA[3], projected_SA[4], Radii[0], Radii[1], Radii[2], db, dg_secondary);
