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
#include "MC_parallel.hpp"
#include "VolumeSA_calculations.hpp"
#include "FreeEnergy.hpp"
#include <chrono>
#include <mpi.h>


using namespace std;

/* In testing of hourglass, it was found that delta = 0.05 works to give better surface area for a radii length scale of a a few units i.e. r=3 or 4. */

int myRank, nProcs;
int main(int argc, char * argv[]) {
    
    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &myRank);
    MPI_Comm_size (MPI_COMM_WORLD, &nProcs);
    
    if(myRank == 0)
    {printf("nProcs = %d\n", nProcs);}

    /* Creating a group of processes for each branch For 1e07 points, we want to have at least 10 processes per communicator group. This MPI_Comm_split uses color to assign a communicator to each process and their original rank as the key to assign a rank within the communicator */
    
    int BranchSize = 10;
    /* IMPORTANT: The way the code is setup now, ensure that there are only 4 branches in total i.e. total number of processors are BranchSize*4 */
    
    int color = myRank/BranchSize;
    MPI_Comm branch_comm;
    MPI_Comm_split(MPI_COMM_WORLD, color, myRank, &branch_comm);
    int branch_rank, branch_size;
    MPI_Comm_rank(branch_comm, &branch_rank);
    MPI_Comm_size(branch_comm, &branch_size);
    //printf("original rank=%d\t branch_rank=%d branch_size=%d\n", myRank, branch_rank, branch_size);
    
    /* Creating another communicator containing the roots of each communicator so that V, SA..etc. data can be transferred from them to the original root */
    int color_roots = ( myRank%BranchSize == 0) ? myRank%BranchSize : 1;
    MPI_Comm roots_comm;
    MPI_Comm_split(MPI_COMM_WORLD, color_roots, myRank, &roots_comm);
    int roots_rank, roots_size;
    MPI_Comm_rank(roots_comm, &roots_rank);
    MPI_Comm_size(roots_comm, &roots_size);
    //printf("original rank=%d\t color_roots=%d roots_rank=%d roots_size= %d\n", myRank,color_roots, roots_rank, roots_size);
    
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
    
    /* I/O related variable */
//    FILE* V_SA_spreadoutDataFile;
//    FILE* V_SA_not_spreadoutDataFile;
    FILE* V_SA_DataFile;
    std::string tag;
    

    
    
    /* Reading all variables from command line */
    if(myRank == 0)
    {
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
        tag.assign(argv[21]);
        printf("d_Rg, d_Rb = %10.5f %10.5f\n",d_Rg, d_Rb);
        printf("n_points = %d\n",n_points);
        printf("Rg_max, Rb_max = %10.5f %10.5f\n",Rg_max, Rb_max);
        printf("p_width_good = [%10.5f %10.5f] p_width_bad = [%10.5f %10.5f]\n", pG_width_x, pG_width_y, pB_width_x, pB_width_y);
        printf("delta=%10.5f\n",delta);
        
    }

    //Broadcasting all user input to the parallel processes.
    MPI_Bcast (&pG_width_x, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&pG_width_y, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&pB_width_x, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&pB_width_y, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
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
    MPI_Bcast (&Rho, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&num_patches, 1, MPI_INT, 0,  MPI_COMM_WORLD);

    //MPI_Bcast (&startingRg, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);

    //This length of Rg has to be available to all processes
    int len_Rg = (int) ((Rg_max-0.0)/d_Rg) + 1; //This 1 is added because for loop is i<len_Rg
    int len_Rb; //Calculated seperately for each Rg. Rb's are symmetric around the central patch.
    int len_Rg_secondary; //Calculated seperately for each Rb. They also use the same Rg_max

    /*The I/O part of printing the output is done be rank=0 only*/
    if(myRank == 0)
    {
        std::string V_SA_DataFileName = tag + "_V_SA_data.txt";
        std::cout<<"tag = "<<tag<<std::endl;
        std::cout<<"V_SA_DataFileName = "<<V_SA_DataFileName<<std::endl;
        V_SA_DataFile = fopen(V_SA_DataFileName.c_str(), "w");
        printf("ptr to file = %p\n", V_SA_DataFile);
        if(V_SA_DataFile == NULL)
        {
            printf("Error opening output file\n");
            exit(1);
        }

        //MPI_Send(tag.c_str(), (int)tag.size(), MPI_CHAR, 1, 0, MPI_COMM_WORLD);

    }
//    if(myRank == 1) //Rank 1 will be used for not spreadout cluster
//    {
//        /* This is a way to dynamically recieve message without knowing the size of the incoming message. https://mpitutorial.com/tutorials/dynamic-receiving-with-mpi-probe-and-mpi-status/*/
//        MPI_Status status;
//        int count;
//        MPI_Probe(0, 0, MPI_COMM_WORLD, &status);
//        MPI_Get_count(&status, MPI_CHAR, &count);
//        char* tag_cstr = new char [count];
//        MPI_Recv(tag_cstr, count, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//        tag = tag_cstr;
//        std::cout<<"Tag in rank 1 = "<<tag<<std::endl ;
//        delete [] tag_cstr;
//
//        std::string V_SA_not_spredoutDataFileName = tag + "_V_SA_data_notspreadout.txt";
//        FILE* V_SA_not_spredoutDataFile = fopen(V_SA_not_spredoutDataFileName.c_str(), "w");
//        if(V_SA_not_spredoutDataFile == NULL)
//        {
//            printf("Error opening not spreadout output file\n");
//            exit(1);
//        }
//    }


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
     /*Setting the MC engine */
    
    MC mc_engine (n_points, stripes.box, aSeed, branch_comm);
    double density = ((double)n_points)/mc_engine.box_volume() ;

    if(myRank == 0)
    {
        for(size_t i =0; i<stripes.box.size(); i++)
        {
            for(size_t j =0; j<stripes.box[i].size(); j++)
            {
                printf("box[%d][%d] = %10.5f\t",(int)i,(int)j,stripes.box[i][j]);
            }
        }
        printf("box_volume = %10.5f\n", mc_engine.box_volume());
        printf("potency good=%10.10f\t potency bad=%10.10f\n",potency_factor(theta_cg), potency_factor(theta_cb));
        printf("density = %10.5f\n", density);
    }


    /* Generating the cellList */
    //CellList cell_list (mc_engine.get_points(), R_cutoff, mc_engine.get_box());
    /*Setting up vectors for N_arrays */

    /* Setting up required variables and output data structures*/
    std::vector<int> dB_array, dG_array;
    std::vector<double> Number_particles, Volume_array, SA_array, projected_SA_array, Radii_array;



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

    Shape* Cluster_shape_ptr;

    /* Starting the loop for Rg.*/
    for(int i = 0 ; i< len_Rg; i++) // 1
    {

        Rg = 0.0 + i*d_Rg ;  //  startingRg + i*d_Rg            //sphere's radius

        /* This is because Rg will be same for all leader (root) processes of each branch*/

        double projected_rg = Rg*sin(theta_cg); //projected radii of the circles
        if(myRank == 0)
        {printf("i=%d\t",i); printf("projected_rg = %10.10f\n",projected_rg);}

        Spherical_cap GoodCap (centre_good, projected_rg, theta_cg, z_surface);

        Cluster_shape_ptr= &GoodCap;

        stripes_bounds = stripes.monitor_cluster_spread(Cluster_shape_ptr);

        if(stripes_bounds[0] == 0) //Cap on central cluster is within bounds
        {
            if(myRank == 0)
            {

                Radii_array.push_back(Rg);
                Radii_array.push_back(0.0); Radii_array.push_back(0.0); //i.e. no cluster on the bad or the secondary good patch
                int zero = 0;

                //Here I am using the mc_engine to calculate spherical caps volume as a second line of checking that the mc code is working correctly with the given parameters

                //mc_volume_SA = mc_engine.calc_volume_SA(Cluster_shape_ptr, delta);

                Volume = GoodCap.getVolume(); //mc_volume_SA[0];
                SA = GoodCap.getSA(); //mc_volume_SA[1];
                mc_volume_SA[0]= Volume; mc_volume_SA[1] = SA;
                projected_SA = GoodCap.projected_SA();

                std::vector<double> local_compcluster_projected_SA(5,0.0);
                local_compcluster_projected_SA[0] = projected_SA;

                Volume_array.push_back(Volume);
                SA_array.push_back(SA);

                addQuant(projected_SA_array, local_compcluster_projected_SA, num_patches);

                N = (Volume * 1e-30) * Rho * avogadro ;
                Number_particles.push_back(N);
                dB_array.push_back(0.0);
                dG_array.push_back(0.0);
            }

        }
        else
        {

            /* This is for projected radii i.e. rb and rg are projected radii not sphere radius.
             rb >= sqrt(rg^2 - (a/2)^2 ) */

            double projected_rb_min = sqrt(projected_rg*projected_rg - ((pG_width_x*pG_width_x)/4.0));
            Rb_min = projected_rb_min/(double)sin(theta_cb) ; //Be careful to not have a theta that is 0 or 180

            len_Rb = (int) ((Rb_max-Rb_min)/d_Rb) + 1;

            //Two cases for each Rb value. +1: spread out, -1: not spread out
            /* Even though both color==1 and color==3 will do the same calculation as color==0 and color==2 respectively for Rg==0, the quantities are only added to the arrays for color==0 and color==2 for Rg==0*/
            int dB_list;
            if (color == 0 || color == 1)
            {
                dB_list = 1;
            }
            else if (color == 2 || color == 3)
            {
                dB_list = -1;
            }
            else {throw std::invalid_argument("Number of processes is not 4 for 5 patches");}

                for(int j = 0 ; j<len_Rb; j++)
                {
                    if(myRank == 0)
                    {printf("j = %d\t dB_list = %d\n",j, dB_list);}

                    /* Have a file output the values of near surface points*/
                    Rb = Rb_min + j*d_Rb ;

                    if(Rb > Rb_max){throw std::invalid_argument("Invalid Rb > Rb_max");}
                    double projected_rb = Rb * (double)sin(theta_cb);

                    double inside_sq = (projected_rb*projected_rb) - (projected_rg*projected_rg) + ((pG_width_x*pG_width_x)/4.0);

                    if(abs(inside_sq)<1e-10 && inside_sq < 0.0) {inside_sq = 0.0;}
                    d_B = dB_list * sqrt(inside_sq) ;
                    /* Symmetric caps on either side of central patch */
                    c_bad_left[0] = -(pG_width_x/2.0) - d_B;
                    c_bad_right[0] = (pG_width_x/2.0) + d_B;

                    /*Checking that centres satisfy cb1 < 0 and cb2 > 0 and */
                    if(c_bad_left[0] >= 0.0 || c_bad_right[0] <= 0.0)
                    {
//                        printf("cb1 = %10.10f\tcb2=%10.10f\n",c_bad_left[0],c_bad_right[0]);
                        break;
                    }

                    /*Check for  rb1 + cb1 < rg < rb1 - cb1 which should be satisfied by the previous condition */
                    if(!(projected_rg < projected_rb - c_bad_left[0] && projected_rg > projected_rb + c_bad_left[0] && projected_rg < projected_rb + c_bad_right[0] && projected_rg > projected_rb - c_bad_right[0]))
                    {
                        printf("j=%d\t Project_Rg=%10.10f Project_Rb=%10.10f dB=%10.10f dB_list=%d c_bad_left[0]=%10.10f c_bad_right[0]=%10.10f inside_sq=%10.20f\n",j, projected_rg, projected_rb, d_B, dB_list, c_bad_left[0], c_bad_right[0], inside_sq);
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
                            printf("Box was breached for Rg=%10.10f\tRb=%10.10f\t, db=%d c_bad_left=%10.10f\t c_bad_right = %10.10f\t and q=%d\n",Rg, Rb,dB_list,c_bad_left[0],c_bad_right[0], (int) q);
                            abort();
                            /* Using abort here because for the bad cap the inputs should be such that this does not happen*/
                        }
                    }

                    /* BAD PATCHES CROSSING CHECK */
                    if(stripes_bounds[1]==0 && stripes_bounds[2]==0) //i.e. the bad patches are not crossed
                    {
                        if(color_roots == 0 && (color==0 || color==2))
                        {
                            /* Now adding zero for Rg_secondary because bad patches are not crossed */
                            Radii_array.push_back(Rg);
                            Radii_array.push_back(Rb);
                            Radii_array.push_back(0.0);
//                            printf("roots_rank=%d\t Radii_array size = %zu\n", roots_rank, Radii_array.size());

                            std::vector<double> local_compcluster_projected_SA = Comp_cluster.projected_SAs();
                            int diff_size = num_patches - (int)local_compcluster_projected_SA.size();
                            if(diff_size > 0)
                            {
                                for(int x=0; x<diff_size;x++){local_compcluster_projected_SA.push_back(0.0);}

                            }
                            else if (diff_size < 0)
                            {
                                printf("compcluster proj SA size larger than num_patches\n"); abort();
                            }
                            addQuant(projected_SA_array, local_compcluster_projected_SA, num_patches);
                        }

                        if(color==0 || color==2) //The calculation for the cluster without any secondary good patch is only done for the first two groups
                        {
//                            auto t1 = std::chrono::high_resolution_clock::now();
                            try
                            {
                                mc_volume_SA = mc_engine.calc_volume_SA(Cluster_shape_ptr, delta);
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
//                            auto t2 = std::chrono::high_resolution_clock::now();
//                            auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
                            if(color_roots == 0)
                            {
//                                printf("With bad Volume = %10.10f SA=%10.10f\n", mc_volume_SA[0], mc_volume_SA[1]);
                                Volume_array.push_back(mc_volume_SA[0]);
                                SA_array.push_back(mc_volume_SA[1]);

                                N = (mc_volume_SA[0] * 1e-30) * Rho * avogadro ;
                                Number_particles.push_back(N);
                                dB_array.push_back(dB_list);
                                dG_array.push_back(0.0);

                            }

                        }

                    }

                    /* SECONDARY GOOD PATCH TERRITORY */
                    else //i.e. the bad patches are crossed and so moving into the secondary good patch territory
                    {
                        if(stripes_bounds[1]==1 && stripes_bounds[2]==1)
                        {
                            double dG1 = pB_width_x - d_B ; /* b - dB1*/
                            double projected_rg_secondary_min = sqrt(projected_rb*projected_rb - dG1*dG1);
                            Rg_secondary_min = projected_rg_secondary_min/(double)sin(theta_cg) ; //Be careful to not have a theta that is 0 or 180
                            len_Rg_secondary = (int) ((Rg_max-Rg_secondary_min)/d_Rg) + 1;

                            int dG_secondary_list;
                            if (color == 0 || color == 2)
                            {
                                dG_secondary_list = 1;
                            }
                            else if (color == 1 || color == 3)
                            {
                                dG_secondary_list = -1;
                            }

                            for(int n = 0 ; n<len_Rg_secondary; n++)
                            {
                                if(myRank == 0){printf("n = %d\n", n);}

                                Rg_secondary = Rg_secondary_min + n*d_Rg ;

                                if(color_roots == 0 && (color==0 || color==2))
                                {
                                    Radii_array.push_back(Rg);
                                    Radii_array.push_back(Rb);
                                    Radii_array.push_back(Rg_secondary);
                                    /* Now pushing secondary radius in the two groups that also calculate Rg_secondary=0.0 and already have Rg and Rb pushed to Radii array*/
                                }
                                else if(color_roots == 0 && (color==1 || color==3))
                                {
                                    /* Pushing all radii to other two groups of processes that do not have Rg or Rb pushed from earlier*/
                                    Radii_array.push_back(Rg);
                                    Radii_array.push_back(Rb);
                                    Radii_array.push_back(Rg_secondary);
                                }

                                double projected_rg_secondary = Rg_secondary * (double)sin(theta_cg);

                                double inside_sq_Rg_second = (projected_rg_secondary*projected_rg_secondary) - (projected_rb*projected_rb) + (dG1*dG1);

                                if(abs(inside_sq_Rg_second)<1e-10 && inside_sq_Rg_second < 0.0) {inside_sq_Rg_second = 0.0;}
                                d_G1 = dG_secondary_list * sqrt(inside_sq_Rg_second);

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
                                        printf("Box was breached for Rg=%10.10f\tRb=%10.10f\tRg_secondary=%10.10f\t, db=%d c_bad_left=%10.10f\t c_bad_right = %10.10f\t c_good_left=%10.10f\t c_good_right = %10.10f\t and q=%d\n",Rg, Rb, Rg_secondary, dB_list,c_bad_left[0],c_bad_right[0], c_good_left[0],c_good_right[0], (int) q);
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
                                    if(color_roots == 0)
                                    {
                                        /* Now adding zero for Rg_secondary because bad patches are not crossed */

                                        std::vector<double> local_compcluster_projected_SA = Comp_cluster_good_second.projected_SAs();
                                        int diff_size = num_patches - (int)local_compcluster_projected_SA.size();
                                        if(diff_size > 0)
                                        {
                                            for(int x=0; x<diff_size;x++){local_compcluster_projected_SA.push_back(0.0);}

                                        }
                                        else if (diff_size < 0)
                                        {
                                            printf("compcluster proj SA size larger than num_patches\n"); abort();
                                        }
                                        addQuant(projected_SA_array, local_compcluster_projected_SA, num_patches);
                                    }

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
                                    if(color_roots == 0)
                                    {
//                                        printf("With Good secondary Volume = %10.10f SA=%10.10f\n", mc_volume_SA[0], mc_volume_SA[1]);
                                        Volume_array.push_back(mc_volume_SA[0]);
                                        SA_array.push_back(mc_volume_SA[1]);

                                        N = (mc_volume_SA[0] * 1e-30) * Rho * avogadro ;
                                        Number_particles.push_back(N);
                                        dB_array.push_back(dB_list);
                                        dG_array.push_back(dG_secondary_list);
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
    MPI_Barrier(MPI_COMM_WORLD);

    if(color_roots == 0)
    {
        /*These counts are true for single element quantities per grid point like N, Volume, SA, dB and dG*/
        int counts[roots_size];
        int projected_SA_counts[roots_size];
        int radii_counts[roots_size];

        int nelements = (int)Volume_array.size();
        int projected_SA_nelements = (int)projected_SA_array.size();
        int radii_nelements = (int)Radii_array.size();

        printf("process=%d\t root_rank=%d roots_size=%d Volume, SA, N, dB, dG, proj_SA and Radii elements = %d %d %d %d %d %d %d\n", myRank, roots_rank, roots_size, (int)Volume_array.size(), (int)SA_array.size(), (int)Number_particles.size(), (int)dB_array.size(), (int)dG_array.size(), (int)projected_SA_array.size(), (int)Radii_array.size());

        MPI_Gather(&nelements, 1, MPI_INT, &counts[0], 1, MPI_INT, 0, roots_comm);
        MPI_Gather(&projected_SA_nelements, 1, MPI_INT, &projected_SA_counts[0], 1, MPI_INT, 0, roots_comm);
        MPI_Gather(&radii_nelements, 1, MPI_INT, &radii_counts[0], 1, MPI_INT, 0, roots_comm);



        int disps[roots_size];
        int projected_SA_disps[roots_size];
        int radii_disps[roots_size];
        // Displacement for the first chunk of data - 0
        for (int i = 0; i < roots_size; i++)
        {
            //printf("i = %d\t counts = %d\t projected_SA_counts=%d\t radii_counts=%d\n", i, counts[i], projected_SA_counts[i], radii_counts[i]);
            disps[i] = (i > 0) ? (disps[i - 1] + counts[i - 1]) : 0;

            projected_SA_disps[i] = (i > 0) ? (projected_SA_disps[i - 1] + projected_SA_counts[i - 1]) : 0;

            radii_disps[i] = (i > 0) ? (radii_disps[i - 1] + radii_counts[i - 1]) : 0;
            //printf("disps = %d\t projected_SA_disps=%d\t radii_disps=%d\n", disps[i], projected_SA_disps[i], radii_disps[i]);
        }

        std::vector<int> dB_global_array, dG_global_array;
        std::vector<double> Number_particles_global, Volume_global_array, SA_global_array, projected_global_SA_array, Radii_global_array;
        if(myRank == 0)
        {
            int V_SA_size = disps[roots_size - 1] + counts[roots_size - 1];
            int proj_SA_size = projected_SA_disps[roots_size - 1] + projected_SA_counts[roots_size - 1];
            int radii_size = radii_disps[roots_size - 1] + radii_counts[roots_size - 1];

            dB_global_array.resize(V_SA_size);
            dG_global_array.resize(V_SA_size);
            Number_particles_global.resize(V_SA_size);
            Volume_global_array.resize(V_SA_size);
            SA_global_array.resize(V_SA_size);

            projected_global_SA_array.resize(proj_SA_size);

            Radii_global_array.resize(radii_size);

        }


        MPI_Gatherv(Volume_array.data(), nelements, MPI_DOUBLE, Volume_global_array.data(), &counts[0], &disps[0], MPI_DOUBLE, 0, roots_comm);

        MPI_Gatherv(SA_array.data(), nelements, MPI_DOUBLE, SA_global_array.data(), &counts[0], &disps[0], MPI_DOUBLE, 0, roots_comm);

        MPI_Gatherv(Number_particles.data(), nelements, MPI_DOUBLE, Number_particles_global.data(), &counts[0], &disps[0], MPI_DOUBLE, 0, roots_comm);

        MPI_Gatherv(dB_array.data(), nelements, MPI_INT, dB_global_array.data(), &counts[0], &disps[0], MPI_INT, 0, roots_comm);

        MPI_Gatherv(dG_array.data(), nelements, MPI_INT, dG_global_array.data(), &counts[0], &disps[0], MPI_INT, 0, roots_comm);

        MPI_Gatherv(projected_SA_array.data(), projected_SA_nelements, MPI_DOUBLE, projected_global_SA_array.data(), &projected_SA_counts[0], &projected_SA_disps[0], MPI_DOUBLE, 0, roots_comm);

        MPI_Gatherv(Radii_array.data(), radii_nelements, MPI_DOUBLE, Radii_global_array.data(), &radii_counts[0], &radii_disps[0], MPI_DOUBLE, 0, roots_comm);


        if(myRank == 0)
        {
            if(V_SA_DataFile != NULL)
            {
                add_Volume_SA_parallel(Number_particles_global, Radii_global_array, Volume_global_array, SA_global_array, projected_global_SA_array, dB_global_array, dG_global_array , V_SA_DataFile);
                printf("After adding quants\n");
            }
            else
            {printf("File pointer is invalid\n") ;}


        }
    }

    if(myRank == 0)
    {
        if(V_SA_DataFile != NULL)
        {
            fclose(V_SA_DataFile);
            printf("After closing file\n");

        }

    }
    
    MPI_Comm_free(&branch_comm);
    MPI_Finalize();
    return EXIT_SUCCESS;
}
