//
//  main.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 7/1/20.
//  Copyright © 2020 Sarwar Hussain. All rights reserved.
//

#include <iostream>
#include <stdio.h>
#include <vector>
#include <string.h>
#include "Patch.hpp"
#include "Stripes.hpp"
#include "Spherical_cap.hpp"
#include "Composite_cluster.hpp"
#include "MC.hpp"
#include <chrono>
#include "FreeEnergy.hpp"
#include "VolumeSA_calculations.hpp"

//Define these but CONSIDER taking input
//double Rho = 5.4944e4 ; moles/m3
//double Mu = 641.4332101443232;
//double Sigma = 25.19e-03; // +- 2.95 Remains unchanged with Temperature
//double avogadro = 6.022140857e23;
//double kb = 1.38064852e-23;
//double T = 170.0;

using namespace std;

/* In testing of hourglass, it was found that delta = 0.05 works to give better surface area for a radii length scale of a a few units i.e. r=3 or 4. */

int main(int argc, const char * argv[]) {
    
    /// Inputs from command line ///
    
    double Rho ;
    /* Surface related variables */
    double z_surface = 0.0; //Unless stated otherwise
    double pG_width_x, pG_width_y, pB_width_x, pB_width_y, theta_cg, theta_cb;
    int num_patches;
    
    /* Cluster related variables */
    double d_Rg, d_Rb;
    double Rg_max, Rb_max;

    
    /* Monte-  Carlo related variables */
    int n_points;
    int aSeed[3];
    double delta;
    
    /* CellList related variables */
    std::vector<double> R_cutoff(3, 0.0);
    
    /* I/O related variables */
//    FILE* outputfile;
//    FILE* outputfile_spreadout;
//    FILE* outputfile_not_spreadout;
    FILE* V_SA_DataFile;   //The file containing the Volume and SA for each grid point
    
    
    /* Reading all variables from command line */
    pG_width_x  = atof(argv[1]);     //x width of the good patches (same for all good patches) (Angstroms)
    pG_width_y  = atof(argv[2]);     //y width of the good patches (same for all good patches) (Angstroms)
    pB_width_x  = atof(argv[3]);     //x width of the good patches (same for all good patches) (Angstroms)
    pB_width_y  = atof(argv[4]);     //y width of the good patches (same for all good patches) (Angstroms)
    theta_cg    = atof(argv[5]);       //Good patch contact angle
    theta_cb    = atof(argv[6]);       //Bad patch contact angle
    d_Rg        = atof(argv[7]);       //increments in good patch cap radius (Angstroms i.e. d_Rg=0.05 A)
    d_Rb        = atof(argv[8]);       //increments in bad patch cap radius  (Angstroms)
    Rg_max      = atof(argv[9]);       //maximum limit of good patch cap.    (Angstroms)
    Rb_max      = atof(argv[10]);       //maximum limit of bad patch cap      (Angstroms)
    n_points    = atoi(argv[11]);       //Total points in millions in the box for MC
    n_points    = n_points*1e06;
    aSeed[0]    = atoi(argv[12]);       //Seed in x-direction for MC
    aSeed[1]    = atoi(argv[13]);       //Seed in y-direction for MC
    aSeed[2]    = atoi(argv[14]);       //Seed in z-direction for MC
    delta       = atof(argv[15]);      //buffer region width=(2*\delta) for surface points (Angstroms)
    R_cutoff[0] = atof(argv[16]);   //Cell size in x for CellList       (Angstroms)
    R_cutoff[1] = atof(argv[17]);   //Cell size in y for CellList       (Angstroms)
    R_cutoff[2] = atof(argv[18]);   //Cell size in z for CellList       (Angstroms)
    Rho         = atof(argv[19]);   //Moles per m3
    num_patches = atoi(argv[20]);
    std::string tag(argv[21]);
    //double startingRg = atof(argv[25]);
    
    
    
//    Nmin        = atoi(argv[19]);
//    dN          = atof(argv[20]);
//    Nmax        = atoi(argv[21]);
//    Mu          = atof(argv[23]);
//    Sigma       = atof(argv[24]);   //Liquid-crystal surface tension
//    T           = atof(argv[25]);  //Kelvin

    //outputfile = fopen (argv[17],"w"); // Output file
    //pointsFile = fopen (argv[12],"w");
   
    
//    std::string outFileName_spreadout = tag + "-spreadout.txt";
//    std::string outFileName_not_spreadout = tag + "-not_spreadout.txt";
    std::string V_SA_DataFileName = tag + "_V_SA_data.txt";
    
//    std::cout<<"outFileName_spreadout = "<<outFileName_spreadout<<std::endl;
//    std::cout<<"outFileName_not_spreadout = "<<outFileName_not_spreadout<<std::endl;
//    outputfile_spreadout = fopen (outFileName_spreadout.c_str(),"w"); // Output file for the centres of bad patches that are on the bad patch
//    outputfile_not_spreadout = fopen (outFileName_not_spreadout.c_str(),"w"); // Output file for the centres of bad patches that are on the good patch
    std::cout<<"tag = "<<tag<<std::endl;
    std::cout<<"V_SA_DataFileName = "<<V_SA_DataFileName<<std::endl;
    V_SA_DataFile = fopen(V_SA_DataFileName.c_str(), "w");
    if(V_SA_DataFile == NULL)
    {
        printf("Error opening output file\n");
        exit(1);
    }
    
    printf("d_Rg, d_Rb = %10.5f %10.5f\n",d_Rg, d_Rb);
    printf("Rg_max, Rb_max = %10.5f %10.5f\n",Rg_max, Rb_max);
    printf("p_width_good = [%10.5f %10.5f] p_width_bad = [%10.5f %10.5f]\n", pG_width_x, pG_width_y, pB_width_x, pB_width_y);
    printf("delta=%10.5f\n",delta);
    int len_Rg = (int) ((Rg_max-0.0)/d_Rg) + 1; //This 1 is added because for loop is i<len_Rg
    int len_Rb; //Calculated seperately for each Rg. Rb's are symmetric around the central patch.
    int len_Rg_secondary; //Calculated seperately for each Rb. They also use the same Rg_max
    
    /* Settings for 5 Patch surface*/
    
    /* Setting the centreal good patch */
    std::vector<double> centre_good(3);   //Centre of the good patch. Only x-y
    centre_good[0] = 0.0; centre_good[1] = 0.0; centre_good[2] = z_surface;
    std::vector<double> dim_Good(2,0.0);
    std::vector<double> dim_Bad(2,0.0);
    
    dim_Good[0] = pG_width_x; dim_Good[1] = pG_width_y; //Common dimensions for all Good patches
    dim_Bad[0] = pB_width_x; dim_Good[1] = pB_width_y; //Common dimensions for all Bad patches
    Patch good_patch (theta_cg, centre_good, dim_Good);
    
    /* Setting the secondry good patches */
    std::vector<double> centre_good_left(3); //(-2a, 0, 0)
    std::vector<double> centre_good_right(3);//(2a, 0, 0)
    centre_good_left[0] = -(pG_width_x + pB_width_x); centre_good_left[1] = 0.0; centre_good_left[2] = z_surface;
    centre_good_right[0] = pG_width_x + pB_width_x; centre_good_right[1] = 0.0; centre_good_right[2] = z_surface;
    Patch good_patch_left (theta_cg, centre_good_left, dim_Good);
    Patch good_patch_right (theta_cg, centre_good_right, dim_Good);

    /* Setting the bad patches */
    std::vector<double> centre_bad_left(3); //(-a, 0, 0)
    std::vector<double> centre_bad_right(3);//(a, 0, 0)
    centre_bad_left[0] = -(pG_width_x + pB_width_x)/2.0; centre_bad_left[1] = 0.0; centre_bad_left[2] = z_surface;
    centre_bad_right[0] = (pG_width_x + pB_width_x)/2.0; centre_bad_right[1] = 0.0; centre_bad_right[2] = z_surface;
    Patch bad_patch_left (theta_cb, centre_bad_left, dim_Bad);
    Patch bad_patch_right (theta_cb, centre_bad_right, dim_Bad);
    
    /* Setting the surface */
    std::vector<Patch> list_of_patches {good_patch, bad_patch_left, bad_patch_right, good_patch_left, good_patch_right};
    std::vector<std::vector<double> > orientations {centre_good, centre_bad_left, centre_bad_right, centre_good_left, centre_good_right};
    
    Stripes stripes (list_of_patches, orientations, z_surface);
    
    /* Setting the height of the box to R_max */
    double box_height = (Rg_max >= Rb_max) ? 2.0*Rg_max : 2.0*Rb_max ;
    //(theta_cg < pi/2.0) ? Rg_max : 2.0*Rg_max ;
    
    /* Calculating box. Here box[i] = 2-dimensional and represents the boundaries of the box in each direction. */
    stripes.calc_box(box_height);
    for(size_t i =0; i<stripes.box.size(); i++)
    {
        for(size_t j =0; j<stripes.box[i].size(); j++)
        {
            printf("box[%d][%d] = %10.5f\t",(int)i,(int)j,stripes.box[i][j]);
        }
    }
    
    /* Setting the MC engine */
    MC mc_engine (n_points, stripes.box, aSeed);
    printf("box_volume = %10.5f\n", mc_engine.box_volume());
    printf("potency good=%10.10f\t potency bad=%10.10f\n",potency_factor(theta_cg), potency_factor(theta_cb));
    
    double density = ((double)n_points)/mc_engine.box_volume() ;
    printf("density = %10.5f\n", density);
    /* Generating the cellList */
    //CellList cell_list (mc_engine.get_points(), R_cutoff, mc_engine.get_box());
     /*Setting up vectors for N_arrays */
    
    /* Setting up required variables and output data structures*/
    double Volume;
    double SA;
    double projected_SA;
    double Volume_MC;
    double SA_MC;
    double N;   //Number of particles per point (i.e. per cluster)
    double Rg, Rb, Rb_min, Rg_secondary, Rg_secondary_min;
    double d_B, d_G1;
    std::vector<double> c_bad_left(3, 0.0), c_bad_right(3, 0.0); //Projected centres of bad spherical caps. x-values will depend on Rb values.
    std::vector<double> c_good_left(3,0.0), c_good_right(3, 0.0);
    
    c_bad_left[1] = 0.0; c_bad_right[1] = 0.0; c_good_left[1] = 0.0; c_good_right[1] = 0.0;
    
    c_bad_left[2] = z_surface; c_bad_right[2] = z_surface; c_good_left[2] = z_surface; c_good_right[2] = z_surface;
    
    
    std::vector<int> stripes_bounds; //array of (0,1) to check crossing of boundaries
    std::vector<int> stripes_box_breach; //array of (0,1) to check box surface breach
    
    std::vector<double> mc_volume_SA(2,0.0);
    
    printf("Before starting Rg loop\n");
    
    Shape* Cluster_shape_ptr;
    std::vector<double> Radii (3, 0.0); //This is 3 for Rg (centre), Rb and Rg (secondary)
    std::vector<double>compcluster_projected_SA (num_patches, 0.0); //This is specifically for the 5 patch case just for the central patch case
    
    
    /* Starting the loop for Rg.*/
    for(int i = 0 ; i< len_Rg; i++) // 1
    {
        printf("i=%d\t",i);
        Rg = 0.0 + i*d_Rg ;  //  startingRg + i*d_Rg            //sphere's radius
        Radii[0] = Rg;
        
        double projected_rg = Rg*sin(theta_cg); //projected radii of the circles
        printf("projected_rg = %10.10f\n",projected_rg);
        
        Spherical_cap GoodCap (centre_good, projected_rg, theta_cg, z_surface);
        
        Cluster_shape_ptr= &GoodCap;
        
        stripes_bounds = stripes.monitor_cluster_spread(Cluster_shape_ptr);
        
        if(stripes_bounds[0] == 0) //Cap on central cluster is within bounds
        {
            Radii[1] = 0.0; Radii[2] = 0.0; //i.e. no cluster on the bad or the secondary good patch
            int zero = 0;
//            auto t_begin = std::chrono::high_resolution_clock::now();
            //Here I am using the mc_engine to calculate spherical caps volume as a second line of checking that the mc code is working correctly with the given parameters
            //mc_volume_SA = mc_engine.calc_volume_SA(Cluster_shape_ptr, delta);
            Volume = GoodCap.getVolume(); //mc_volume_SA[0];
            SA = GoodCap.getSA(); //mc_volume_SA[1];
            mc_volume_SA[0]= Volume; mc_volume_SA[1] = SA;
            projected_SA = GoodCap.projected_SA();
            compcluster_projected_SA[0] = projected_SA;
            
            add_Volume_SA (Radii, mc_volume_SA, compcluster_projected_SA, Rho, zero, zero, V_SA_DataFile);
            
//            auto t_end = std::chrono::high_resolution_clock::now();
//            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_begin);
//
//        fprintf(outputfile_spreadout, "%10.10f\t%d\t%10.10f\t%10.10f\t%lld\t%d\t%d\n", Rg, zero, Volume*1e30, SA*1e20, duration.count(), stripes_bounds[0], zero);
//            fprintf(outputfile_not_spreadout, "%10.10f\t%d\t%10.10f\t%10.10f\t%lld\t%d\t%d\n", Rg, zero, Volume*1e30, SA*1e20, duration.count(), stripes_bounds[0], zero);
        }
        else
        {
            
            /* This is for projected radii i.e. rb and rg are projected radii not sphere radius.
             rb >= sqrt(rg^2 - (a/2)^2 ) */
            
            double projected_rb_min = sqrt(projected_rg*projected_rg - ((pG_width_x*pG_width_x)/4.0));
            Rb_min = projected_rb_min/(double)sin(theta_cb) ; //Be careful to not have a theta that is 0 or 180
            
            len_Rb = (int) ((Rb_max-Rb_min)/d_Rb) + 1;
            
            printf("projected_rb_min = %10.5f\t Rb_min=%10.5f\n", projected_rb_min, Rb_min);
            
            //Two cases for each Rb value. +1: spread out, -1: not spread out
            std::vector<int> dB_list = {+1, -1};
            for(int k=0; k<dB_list.size(); k++)
            {
                for(int j = 0 ; j<len_Rb; j++)
                {
//                    printf("j = %d\tdB_list = %d\n",j, dB_list[k]);
                    /* Have a file output the values of near surface points*/
                    Rb = Rb_min + j*d_Rb ;
                    Radii[1] = Rb; Radii[2] = 0.0;
                    
                    
                    if(Rb > Rb_max){throw std::invalid_argument("Invalid Rb > Rb_max");}
                    double projected_rb = Rb * (double)sin(theta_cb);
                
                    double inside_sq = (projected_rb*projected_rb) - (projected_rg*projected_rg) + ((pG_width_x*pG_width_x)/4.0);
                    if(abs(inside_sq)<1e-10 && inside_sq < 0.0) {inside_sq = 0.0;}
                    d_B = dB_list[k] * sqrt(inside_sq) ;
                    /* Symmetric caps on either side of central patch */
                    c_bad_left[0] = -(pG_width_x/2.0) - d_B;
                    c_bad_right[0] = (pG_width_x/2.0) + d_B;
                    
                    /*Checking that centres satisfy cb1 < 0 and cb2 > 0 and */
                    if(c_bad_left[0] >= 0.0 || c_bad_right[0] <= 0.0)
                    {
                        printf("cb1 = %10.10f\tcb2=%10.10f\n",c_bad_left[0],c_bad_right[0]);
                        break;
                    }
                    
                    /*Check for  rb1 + cb1 < rg < rb1 - cb1 which should be satisfied by the previous condition */
                    if(!(projected_rg < projected_rb - c_bad_left[0] && projected_rg > projected_rb + c_bad_left[0] && projected_rg < projected_rb + c_bad_right[0] && projected_rg > projected_rb - c_bad_right[0]))
                    {
                        printf("j=%d\t Project_Rg=%10.10f Project_Rb=%10.10f dB=%10.10f dB_list[k]=%d c_bad_left[0]=%10.10f c_bad_right[0]=%10.10f inside_sq=%10.20f\n",j, projected_rg, projected_rb, d_B, dB_list[k], c_bad_left[0], c_bad_right[0], inside_sq);
                        throw std::invalid_argument("Invalid Cb1 and Cb2");
                    }

                    Spherical_cap BadCap_left(c_bad_left, projected_rb, theta_cb, z_surface);
                    Spherical_cap BadCap_right(c_bad_right, projected_rb, theta_cb, z_surface);
                    std::vector<Spherical_cap> capsList = {GoodCap, BadCap_left, BadCap_right};
                    
                    Composite_cluster Comp_cluster (capsList, stripes);
                    Cluster_shape_ptr = &Comp_cluster;
                    stripes_bounds = stripes.monitor_cluster_spread(Cluster_shape_ptr);
                    stripes_box_breach = stripes.monitor_box_breach(Cluster_shape_ptr);
                    
                    for(size_t q=0; q<stripes_box_breach.size(); q++)
                    {
                        if(stripes_box_breach[q] == 1)
                        {
                            printf("Box was breached for Rg=%10.10f\tRb=%10.10f\t, db=%d c_bad_left=%10.10f\t c_bad_right = %10.10f\t and q=%d\n",Rg, Rb,dB_list[k],c_bad_left[0],c_bad_right[0], (int) q);
                            abort();
                        }
                    }
                    
                    if(stripes_bounds[1]==0 && stripes_bounds[2]==0) //i.e. the bad patches are not crossed
                    {
                        compcluster_projected_SA.clear();
                        compcluster_projected_SA = Comp_cluster.projected_SAs();
                        int diff_size = num_patches - (int)compcluster_projected_SA.size();
                        if(diff_size > 0)
                        {for(int x=0; x<diff_size;x++){compcluster_projected_SA.push_back(0.0);}}
                        else if (diff_size < 0){printf("compcluster proj SA size larger than num_patches\n"); abort();}
                        
//                        auto t1 = std::chrono::high_resolution_clock::now();
                        try
                        {
                                mc_volume_SA = mc_engine.calc_volume_SA(Cluster_shape_ptr, delta);
                                //mc_volume_SA = mc_engine.update_volume_SA(Cluster_shape_ptr, &cell_list, delta);

                        }
                        catch(const std::logic_error& e)
                        {
                            if(Rg==0.0 && Rb==0.0)
                            {
                                Volume_MC=0.0;
                                SA_MC = 0.0;
                            }
                            else
                            {
                                std::cerr << e.what();
                                abort();
                            }
                        }
//                        auto t2 = std::chrono::high_resolution_clock::now();
//                        auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
                        add_Volume_SA (Radii, mc_volume_SA, compcluster_projected_SA, Rho, dB_list[k], 0, V_SA_DataFile);
//                        printf("duration = %lld\n",duration_ms.count());
                        
//                        if(dB_list[k] == 1)
//                        {
//                            fprintf(outputfile_spreadout,"%10.10f\t%10.10f\t%10.10f\t%10.10f\t%lld\t%d\t%10.10f\n", Rg, Rb, Volume_MC * 1e30, SA_MC * 1e20, duration_ms.count(),stripes_bounds[0], d_B);
//                            
//                        }
//                        else
//                        {
//                            fprintf(outputfile_not_spreadout,"%10.10f\t%10.10f\t%10.10f\t%10.10f\t%lld\t%d\t%10.10f\n", Rg, Rb, Volume_MC * 1e30, SA_MC * 1e20, duration_ms.count(),stripes_bounds[0], d_B);
//                        }
                    }
                    else //i.e. the bad patches are crossed and so moving into the secondary good patch territory
                    {
                        if(stripes_bounds[1]==1 && stripes_bounds[2]==1)
                        {
                            double dG1 = pB_width_x - d_B ; /* b - dB1*/
                            double projected_rg_secondary_min = sqrt(projected_rb*projected_rb - dG1*dG1);
                            Rg_secondary_min = projected_rg_secondary_min/(double)sin(theta_cg) ; //Be careful to not have a theta that is 0 or 180
                            len_Rg_secondary = (int) ((Rg_max-Rg_secondary_min)/d_Rg) + 1;
                            
                            //Two cases for each Rg secondary value. +1: spread out, -1: not spread out
                            std::vector<int> dG_secondary_list = {+1, -1};
                            for(int m=0; m<dG_secondary_list.size(); m++)
                            {
                                for(int n = 0 ; n<len_Rg_secondary; n++)
                                {
//                                    printf("n = %d\n", n);
                                    Rg_secondary = Rg_secondary_min + n*d_Rg ;
                                    Radii[2] = Rg_secondary;
                                    
                                    double projected_rg_secondary = Rg_secondary * (double)sin(theta_cg);
                                    
                                    double inside_sq_Rg_second = (projected_rg_secondary*projected_rg_secondary) - (projected_rb*projected_rb) + (dG1*dG1);
                                    
                                    if(abs(inside_sq_Rg_second)<1e-10 && inside_sq_Rg_second < 0.0) {inside_sq_Rg_second = 0.0;}
                                    d_G1 = dG_secondary_list[m] * sqrt(inside_sq_Rg_second);
                                    
                                    c_good_left[0] = -(pB_width_x + (pG_width_x/2.0)) - d_G1;
                                    c_good_right[0] = (pB_width_x + (pG_width_x/2.0)) + d_G1;
                                    
                                    /*Checking that centres satisfy rb1 - (cb1 - cg1)< rg1 < rb1 + (cb1 - cg1)
                                     and rb2 - (cg2 - cb2)< rg2 < rb2 + (cg2 - cb2)*/
                                    bool good_left_and_bad_left = ((projected_rb - (c_bad_left[0] - c_good_left[0]) < projected_rg_secondary) && projected_rg_secondary < projected_rb + (c_bad_left[0] - c_good_left[0]));
                                    
                                    bool good_right_and_bad_right = (projected_rb - (c_good_right[0] - c_bad_right[0]) < projected_rg_secondary && projected_rg_secondary < projected_rb + (c_good_right[0] - c_bad_right[0]) ) ;
                                    
                                    if(!(good_left_and_bad_left && good_right_and_bad_right))
                                    {
                                        printf("good_left_bad_left =%d\t good_right_and_bad_right=%d\n",good_left_and_bad_left, good_right_and_bad_right);
                                        if(good_left_and_bad_left != good_right_and_bad_right)
                                        {
                                            printf("bad_left = %10.5f\t good_left = %10.5f\t bad_right = %10.5f\t good_right = %10.5f\tn", c_bad_left[0], c_good_left[0], c_bad_right[0], c_good_right[0]);
                                            printf("Symmetry is broken\n"); abort();
                                        }
                                        break;
                                    }
                                    
                                    Spherical_cap GoodCap_left(c_good_left, projected_rg_secondary, theta_cg, z_surface);
                                    Spherical_cap GoodCap_right(c_good_right, projected_rg_secondary, theta_cg, z_surface);
                                    std::vector<Spherical_cap> capsList_good_second = {GoodCap, BadCap_left, BadCap_right, GoodCap_left, GoodCap_right};
                                    
                                    Composite_cluster Comp_cluster_good_second (capsList_good_second, stripes);
                                    Cluster_shape_ptr = &Comp_cluster_good_second;
                                    stripes_bounds = stripes.monitor_cluster_spread(Cluster_shape_ptr);
                                    stripes_box_breach = stripes.monitor_box_breach(Cluster_shape_ptr);
                                    
                                    for(size_t q=0; q<stripes_box_breach.size(); q++)
                                    {
                                        if(stripes_box_breach[q] == 1 && !(stripes_bounds[3]==1 && stripes_bounds[4]==1))
                                        {
                                            printf("Box was breached for Rg=%10.10f\tRb=%10.10f\tRg_secondary=%10.10f\t, db=%d c_bad_left=%10.10f\t c_bad_right = %10.10f\t c_good_left=%10.10f\t c_good_right = %10.10f\t and q=%d\n",Rg, Rb, Rg_secondary, dB_list[k],c_bad_left[0],c_bad_right[0], c_good_left[0],c_good_right[0], (int) q);
                                            abort();
                                        }
                                    }
                                    if(stripes_bounds[3]==1 && stripes_bounds[4]==1) //i.e. the secondary good patches are  crossed
                                    {
                                        printf("Secondary good patch boundaries are crossed at c_good_left=%10.10f\t c_good_right=%10.10f\t and projected_rg_second=%10.10f\n", c_good_left[0], c_good_right[0], projected_rg_secondary);
                                        break;
                                    }
                                    else
                                    {
                                        std::vector<double> compcluster_projected_SA_good_second(capsList_good_second.size(),0.0);
                                        compcluster_projected_SA_good_second = Comp_cluster_good_second.projected_SAs();
                                        int diff_size = num_patches - (int)compcluster_projected_SA_good_second.size();
                                        if(diff_size > 0)
                                        {for(int x=0; x<diff_size;x++){compcluster_projected_SA_good_second.push_back(0.0);}}
                                        else if (diff_size < 0){printf("compcluster secondary good proj SA size larger than num_patches\n"); abort();}
                                        
                                        try
                                        {
                                            mc_volume_SA = mc_engine.calc_volume_SA(Cluster_shape_ptr, delta);
                                        }
                                        catch(const std::logic_error& e)
                                        {
                                            if(Rg==0.0 && Rb==0.0 && Rg_secondary == 0.0)
                                            {
                                                Volume_MC=0.0;
                                                SA_MC = 0.0;
                                            }
                                            else
                                            {
                                                std::cerr << e.what();
                                                abort();
                                            }
                                        }
                                        add_Volume_SA (Radii, mc_volume_SA, compcluster_projected_SA_good_second, Rho, dB_list[k], dG_secondary_list[m], V_SA_DataFile);
                                        
                                    }
                                }
                                
                            }
                            
                        }
                        else
                        {
                            printf("The bad patches are not crossed symmetrically \n");
                            abort();
                        }
                    }

                    
                }
            }
        }
        
    }
    printf("Loop finished\n");
    fclose(V_SA_DataFile);
    
//    fclose(outputfile_spreadout);
//    fclose(outputfile_not_spreadout);
    return 0;
}
