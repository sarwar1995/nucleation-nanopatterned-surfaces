//
//  main_parallel.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 2/7/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.


#include <iostream>
#include <stdio.h>
#include <vector>
#include <string.h>
//#include "Composite_cluster.hpp"
//#include "Surface/Surface.hpp"
#include "Patch.hpp"
#include "Stripes.hpp"
#include "Spherical_cap.hpp"
#include "Composite_cluster.hpp"
#include "MC.hpp"
#include "FreeEnergy.hpp"
#include <chrono>
#include <mpi.h>


//Define these but CONSIDER taking input
double Rho;         //Moles per m3
double Mu;
double Sigma;       //Liquid-crystal surface tension
double T;
#define avogadro 6.022140857e23;
#define kb 1.38064852e-23;


using namespace std;

/* In testing of hourglass, it was found that delta = 0.05 works to give better surface area for a radii length scale of a a few units i.e. r=3 or 4. */

int myRank, nProcs;
int main(int argc, char * argv[]) {
    
    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &myRank);
    MPI_Comm_size (MPI_COMM_WORLD, &nProcs);
    
    if(myRank == 0)
    {printf("nProcs = %d\n", nProcs);}
    if (nProcs > 2)
    {
        throw std::invalid_argument("Number of mpi processes must be 2");
    }
    /// Inputs from command line ///
    /* Surface related variables */
    double z_surface = 0.0; //Unless stated otherwise
    double p_width_x, p_width_y, theta_cg, theta_cb;
    
    /* Cluster related variables */
    double d_Rg, d_Rb;
    double Rg_max, Rb_max;
    //    int Nmin, Nmax;
    //    double dN;
    
    /* Monte-  Carlo related variables */
    int n_points;
    int aSeed[3];
    double delta;
    
    /* CellList related variables */
    std::vector<double> R_cutoff(3, 0.0);
    
    /* I/O related variables */
    //    FILE* outputfile;
    std::string tag;
    
    FILE* outputfile_spreadout;
    FILE* outputfile_not_spreadout;
    
    
    //    FILE* pointsFile_1;
    //    FILE* pointsFile_2;
    
    int len_Rg;
    int len_Rb; //Calculated seperately for each Rg. Rb's are symmetric around the central patch.
    double startingRg ;
    
    /* Reading all variables from command line */
    if(myRank == 0)
    {
        p_width_x  = atof(argv[1]);     //x width of the patches (same for all)
        p_width_y  = atof(argv[2]);     //y width of the patches (same for all)
        theta_cg = atof(argv[3]);       //Good patch contact angle
        theta_cb = atof(argv[4]);       //Bad patch contact angle
        d_Rg     = atof(argv[5]);       //increments in good patch cap radius
        d_Rb     = atof(argv[6]);       //increments in bad patch cap radius
        Rg_max   = atof(argv[7]);       //maximum limit of good patch cap
        Rb_max   = atof(argv[8]);       //maximum limit of bad patch cap
        n_points = atoi(argv[9]);       //Total points in millions in the box for MC
        n_points = n_points*1e06;
        aSeed[0] = atoi(argv[10]);       //Seed in x-direction for MC
        aSeed[1] = atoi(argv[11]);       //Seed in y-direction for MC
        aSeed[2] = atoi(argv[12]);       //Seed in z-direction for MC
        delta    = atof(argv[13]);      //buffer region width=(2*\delta) for surface points
        R_cutoff[0] = atof(argv[14]);   //Cell size in x for CellList
        R_cutoff[1] = atof(argv[15]);   //Cell size in y for CellList
        R_cutoff[2] = atof(argv[16]);   //Cell size in z for CellList
        //outputfile = fopen (argv[17],"w"); // Output file
        //pointsFile = fopen (argv[12],"w");
        tag = argv[17];
        startingRg = atof(argv[18]);
        
        printf("d_Rg, d_Rb = %10.5f %10.5f\n",d_Rg, d_Rb);
        printf("Rg_max, Rb_max = %10.5f %10.5f\n",Rg_max, Rb_max);
        printf("delta=%10.5f\n",delta);
    }
    
    //Broadcasting all user input to the parallel processes. ONLY 2 processes.
    MPI_Bcast (&p_width_x, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&p_width_y, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&theta_cg, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&theta_cb, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&d_Rg, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&d_Rb, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&Rg_max, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&Rb_max, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&n_points, 1, MPI_INT, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&aSeed[0], 3, MPI_INT, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&delta, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&R_cutoff[0], 3, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&startingRg, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    
    //This length of Rg has to be available to all processes
    len_Rg = (int) ((Rg_max-0.0)/d_Rg) + 1; //This 1 is added because for loop is i<len_Rg
    
    if(myRank == 0) //Rank 0 will be used for spreadout cluster
    {
        MPI_Send(tag.c_str(), (int)tag.size(), MPI_CHAR, 1, 0, MPI_COMM_WORLD);
        std::string outFileName_spreadout = tag + "-spreadout.txt";
        std::cout<<"outFileName_spreadout = "<<outFileName_spreadout<<std::endl;
        outputfile_spreadout = fopen (outFileName_spreadout.c_str(),"w"); // Output file for the centres of bad patches that are on the bad patch
    }
    if(myRank == 1) //Rank 1 will be used for not spreadout cluster
    {
        /* This is a way to dynamically recieve message without knowing the size of the incoming message. https://mpitutorial.com/tutorials/dynamic-receiving-with-mpi-probe-and-mpi-status/*/
        MPI_Status status;
        int count;
        MPI_Probe(0, 0, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_CHAR, &count);
        char* tag_cstr = new char [count];
        MPI_Recv(tag_cstr, count, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        tag = tag_cstr;
        std::cout<<"Tag in rank 1 = "<<tag<<std::endl ;
        delete [] tag_cstr;
        
        std::string outFileName_not_spreadout = tag + "-not_spreadout.txt";
        std::cout<<"outFileName_not_spreadout = "<<outFileName_not_spreadout<<std::endl;
        outputfile_not_spreadout = fopen (outFileName_not_spreadout.c_str(),"w"); // Output file for the centres of bad patches that are on the good patch
    }

    
    MPI_Barrier(MPI_COMM_WORLD);
    
    /* Setting the good patch */
    std::vector<double> centre_good(3);   //Centre of the good patch. Only x-y
    centre_good[0] = 0.0; centre_good[1] = 0.0; centre_good[2] = z_surface;
    std::vector<double> dim(2,0.0);
    dim[0] = p_width_x; dim[1] = p_width_y; //Common dimension for all patches
    Patch good_patch (theta_cg, centre_good, dim);
    
    
    /* Setting the bad patches */
    std::vector<double> centre_bad_left(3); //(-a, 0, 0)
    std::vector<double> centre_bad_right(3);//(a, 0, 0)
    centre_bad_left[0] = -p_width_x; centre_bad_left[1] = 0.0; centre_bad_left[2] = z_surface;
    centre_bad_right[0] = p_width_x; centre_bad_right[1] = 0.0; centre_bad_right[2] = z_surface;
    Patch bad_patch_left (theta_cb, centre_bad_left, dim);
    Patch bad_patch_right (theta_cb, centre_bad_right, dim);
    
    /* Setting the surface */
    std::vector<Patch> list_of_patches;
    std::vector<std::vector<double> > orientations;
    orientations.push_back(centre_good);        //Good Patch centre
    orientations.push_back(centre_bad_left);    //Bad left Patch centre
    orientations.push_back(centre_bad_right);   //Bad right Patch centre
    list_of_patches.push_back(good_patch);
    list_of_patches.push_back(bad_patch_left);
    list_of_patches.push_back(bad_patch_right);
    
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
    
    if(myRank == 0)
    {
        printf("box_volume = %10.5f\n", mc_engine.box_volume());
        printf("potency good=%10.10f\t potency bad=%10.10f\n",potency_factor(theta_cg), potency_factor(theta_cb));
    }
    /* Generating the cellList */
    CellList cell_list (mc_engine.get_points(), R_cutoff, mc_engine.get_box());
    
    /* Setting up vectors for N_arrays */
    //    std::vector<double> NArray_Gmin;
    //    std::vector<double> NArray_confs;
    //    std::vector<std::vector<double> > NArray_quant;
    
    /* Setting up required variables and output data structures*/
    double Volume;
    double SA;
    double Volume_MC;
    double SA_MC;
    double Rg, Rb, Rb_min;
    double d_B;
    std::vector<double> c_bad_left(3, 0.0), c_bad_right(3, 0.0); //Projected centres of bad spherical caps. x-values will depend on Rb values.
    c_bad_left[1] = 0.0; c_bad_right[1] = 0.0;
    c_bad_left[2] = z_surface; c_bad_right[2] = z_surface;
    
    
    std::vector<int> stripes_bounds; //array of (0,1) to check crossing of boundaries
    std::vector<int> stripes_box_breach; //array of (0,1) to check box surface breach
    
    std::vector<double> mc_volume_SA(2,0.0);
    
    if(myRank == 0)
    { printf("Before starting Rg loop\n"); }
    
    Shape* Cluster_shape_ptr;
    
    /* Starting the loop for Rg.*/
    for(int i = 0 ; i< 1; i++) // len_Rg
    {
        printf("i=%d\t",i);
        Rg = startingRg + i*d_Rg ;  //  0.0 + i*d_Rg            //sphere's radius
        
        double projected_rg = Rg*sin(theta_cg); //projected radii of the circles
        printf("projected_rg = %10.10f\n",projected_rg);
        
        Spherical_cap GoodCap (centre_good, projected_rg, theta_cg, z_surface);
        
        Cluster_shape_ptr= &GoodCap;
        
        stripes_bounds = stripes.monitor_cluster_spread(Cluster_shape_ptr);
        
        if(stripes_bounds[0] == 0) //Cap on central cluster is within bounds
        {
            
                auto t_begin = std::chrono::high_resolution_clock::now();
                Volume = GoodCap.getVolume();
                SA = GoodCap.getSA();
                auto t_end = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_begin);
                int zero = 0;
            if(myRank == 0)
            {
                fprintf(outputfile_spreadout, "%10.10f\t%d\t%10.10f\t%10.10f\t%lld\t%d\t%d\n", Rg, zero, Volume, SA, duration.count(), stripes_bounds[0], zero);
                
            }
            if(myRank == 1)
            {
                fprintf(outputfile_not_spreadout, "%10.10f\t%d\t%10.10f\t%10.10f\t%lld\t%d\t%d\n", Rg, zero, Volume, SA, duration.count(), stripes_bounds[0], zero);
            }
            
        }
        else
        {
            
            /* This is for projected radii i.e. rb and rg are projected radii not sphere radius.
             rb >= sqrt(rg^2 - (a/2)^2 ) */
            
            double projected_rb_min = sqrt(projected_rg*projected_rg - ((p_width_x*p_width_x)/4.0));
            Rb_min = projected_rb_min/(double)sin(theta_cb) ;
            
            len_Rb = (int) ((Rb_max-Rb_min)/d_Rb) + 1;
            
            //Two cases for each Rb value. +1: spread out, -1: not spread out
            
            int dB_list;
            if (myRank == 0)
            {
                dB_list = 1;
            }
            else if (myRank == 1)
            {
                dB_list = -1;
            }
            else {throw std::invalid_argument("Number of processes is more than two");}
            
//            std::vector<int> dB_list = {+1, -1};
            
                for(int j = 0 ; j<len_Rb; j++)
                {
                    printf("rank = %d\t j = %d\tdB_list = %d\n",myRank, j, dB_list);
                    /* Have a file output the values of near surface points*/
                    Rb = Rb_min + j*d_Rb ;
                    
                    if(Rb > Rb_max){throw std::invalid_argument("Invalid Rb > Rb_max");}
                    double projected_rb = Rb*sin(theta_cb);
                    
                    
                    d_B = dB_list * sqrt(projected_rb*projected_rb - projected_rg*projected_rg + ((p_width_x*p_width_x)/4.0)) ;
                    printf("dB = %10.10f\n",d_B);
                    /* Symmetric caps on either side of central patch */
                    c_bad_left[0] = -(p_width_x/2.0) - d_B;
                    c_bad_right[0] = (p_width_x/2.0) + d_B;
                    
                    /*Checking that centres satisfy cb1 < 0 and cb2 > 0 and */
                    if(c_bad_left[0] >= 0.0 || c_bad_right[0] <= 0.0)
                    {
                        printf("cb1 = %10.10f\tcb2=%10.10f\n",c_bad_left[0],c_bad_right[0]);
                        break;
                    }
                    
                    /*Check for  rb1 + cb1 < rg < rb1 - cb1 which should be satisfied by the previous condition */
                    if(!(projected_rg < projected_rb - c_bad_left[0] && projected_rg > projected_rb + c_bad_left[0]))
                    {
                        printf("j=%d\t Project_Rg=%10.10f Project_Rb=%10.10f dB=%10.10f dB_list=%d\n",j, projected_rg, projected_rb, d_B, dB_list);
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
                            printf("Box was breached for Rg=%10.10f\tRb=%10.10f\t and db=%d\n",Rg, Rb,dB_list);
                            abort();
                        }
                    }
                    
                    if(stripes_bounds[1]==1 && stripes_bounds[2]==1)
                    {
                        printf("Bad patch bounds are crossed at Rb=%10.10f\t and db=%d\n",Rb,dB_list);
                        break;
                    }
                    
                    
                    auto t1 = std::chrono::high_resolution_clock::now();
                    try
                    {
                        
                        if(j==0) //Beginning of growing bad patches
                        {
                            
                            mc_volume_SA = mc_engine.calc_volume_SA(Cluster_shape_ptr, delta);
                            
                        }
                        else
                        {
                            mc_volume_SA = mc_engine.calc_volume_SA(Cluster_shape_ptr, delta);
                            
                            //mc_volume_SA = mc_engine.update_volume_SA(Cluster_shape_ptr, &cell_list, delta);
                            
                            //                            //This printing part should be after the calc_volume_SA call
                            //                            if(dB_list == -1)
                            //                            {
                            //                                if(j == len_Rb-2)
                            //                                {
                            //                                    fprintf(pointsFile_1,"Rg=%10.10f\t Rb=%10.10f\t db=%d\n", Rg, Rb, dB_list);
                            //                                    mc_engine.print_surf_points(pointsFile_1);
                            //                                }
                            //                                if(j == len_Rb-3)
                            //                                {
                            //                                    fprintf(pointsFile_2,"Rg=%10.10f\t Rb=%10.10f\t db=%d\n", Rg, Rb, dB_list);
                            //                                    mc_engine.print_surf_points(pointsFile_2);
                            //                                }
                            //                            }
                            
                        }
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
                    auto t2 = std::chrono::high_resolution_clock::now();
                    auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
                    Volume_MC = mc_volume_SA[0];
                    SA_MC = mc_volume_SA[1];
                    if(myRank == 0)
                    {
                        fprintf(outputfile_spreadout,"%10.10f\t%10.10f\t%10.10f\t%10.10f\t%lld\t%d\t%10.10f\n", Rg, Rb, Volume_MC, SA_MC, duration_ms.count(),stripes_bounds[0], d_B);}
                    else if (myRank == 1)
                    {
                        fprintf(outputfile_not_spreadout,"%10.10f\t%10.10f\t%10.10f\t%10.10f\t%lld\t%d\t%10.10f\n", Rg, Rb, Volume_MC, SA_MC, duration_ms.count(),stripes_bounds[0], d_B);
                    }
                    //fprintf(outputfile, "%10.10f\t%10.10f\t%lld\t%d\n", Volume_MC, SA_MC, duration_ms.count(), stripes_bounds[0]);
                }
            }
        
        
    }
    
    if (myRank == 0)
    {
        fclose(outputfile_spreadout);
    }
    else if(myRank == 1)
    {
        fclose(outputfile_not_spreadout);
    }
    return 0;
    MPI_Finalize();
}
