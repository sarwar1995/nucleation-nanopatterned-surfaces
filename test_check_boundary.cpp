//
//  test_check_boundary.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 7/1/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//


#include <iostream>
#include <stdio.h>
#include <vector>
#include <cstring>
#include <string.h>
#include <cmath>
//#include "Composite_cluster.hpp"
//#include "Surface/Surface.hpp"
#include "Patch.hpp"
#include "Stripes.hpp"
#include "Spherical_cap.hpp"
#include "Composite_cluster.hpp"
#include "MC_parallel.hpp"
#include "surface_setup.hpp"
#include "DynamicBox.hpp"
#include "CheckBoundary.hpp"
#include "FreeEnergy.hpp"

int myRank, nProcs;

void get_cap_centres(double projected_radius_g, double projected_radius_b, double dist_from_boundary, double prev_patch_boundary, double theta_bad,  int is_bad_patch, std::vector<double>&centre_left_cap, std::vector<double>&centre_right_cap)
{
    double d_ch = sqrt(projected_radius_g*projected_radius_g - dist_from_boundary*dist_from_boundary);
    printf("is_bad_patch = %d d_ch = %10.10f\t min Radius=%10.5f\n", is_bad_patch, d_ch, (d_ch/sin(theta_bad)));
    
    double d_B = sqrt(projected_radius_b*projected_radius_b - d_ch*d_ch);
    printf("d_B=%10.5f\n", d_B);
    int dB_mod = 1;
    if(is_bad_patch==1){dB_mod = -1;}
    centre_left_cap[0] = -(prev_patch_boundary) - dB_mod*d_B ;
    centre_right_cap[0] = (prev_patch_boundary) + dB_mod*d_B ;
    printf("prev_patch_boundary=%10.5f\t centre_left_cap=%10.5f\t centre_right_cap=%10.5f\n", prev_patch_boundary, centre_left_cap[0], centre_right_cap[0]);
    
}

void calc_vol_SA_normal(Shape* Cluster_shape_ptr, CheckBoundary* check_boundary, MC*  mc_engine, std::vector<double>& mc_volume_SA)
{
    printf("Inside normal V_SA calc \n");
    check_boundary->ManageBoxBreach(Cluster_shape_ptr);
    mc_volume_SA = mc_engine->calc_volume_SA(Cluster_shape_ptr);
}

void calc_vol_SA_virtual_points(Shape* Cluster_shape_ptr, CheckBoundary* check_boundary, MC*  mc_engine, std::vector<double>& mc_volume_SA)
{
    printf("Inside virtual V_SA calc \n");
    int n_inside_virtual = mc_engine->get_n_inside_virtual();
    int n_surf_virtual = mc_engine->get_n_near_surf_virtual();
    if(n_inside_virtual != 0 && n_surf_virtual!=0)
    {
        printf("n_inside_virtual=%d\t n_surf_virtual=%d\n are not zero before adding points\n", n_inside_virtual, n_surf_virtual);
        abort();
    }
    int isBreach = check_boundary->ManageBoxBreach_and_calc_Vol_SA(Cluster_shape_ptr);
    if (isBreach == 1)
    {
        mc_volume_SA = mc_engine->calc_volume_SA_with_virtual_points(Cluster_shape_ptr);
        
    }
    else
    {
        mc_volume_SA = mc_engine->calc_volume_SA(Cluster_shape_ptr);
    }
    check_boundary->reset_n_breach_cycles();
}

int main(int argc, char * argv[])
{
    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &myRank);
    MPI_Comm_size (MPI_COMM_WORLD, &nProcs);
    
    double Rg, Rb;
    double Rg_second, Rb_third;
    double Rg_max, Rb_max, d_Rg, d_Rb;
    int num_patches;
    double z_surface = 0.0;
    double good_patch_x_width, good_patch_y_width, bad_patch_x_width, bad_patch_y_width;
    double theta_good;
    double theta_bad;
    double starting_box_dim;
    double extension_length;
    /* Patch variables*/
    int last_patch_isinfinite;
    /* MC variables */
    int n_starting_points;
    double point_density;
    int mc_seed[3];
    double delta;
    int n_max_points, starting_patch_is_bad;
    /* physical variables */
    double Rho ;

    double Mu, Sigma, T;
    std::string tag;
    
    if(myRank==0)
    {
        Rg = atof(argv[1]);
        Rb = atof(argv[2]);
        Rg_second = atof(argv[3]);
        Rb_third = atof(argv[4]);
        int start_index = 4;
        good_patch_x_width  = atof(argv[start_index + 1]);     //x width of the good patches (same for all good patches) (Angstroms)
        good_patch_y_width  = atof(argv[start_index + 2]);     //y width of the good patches (same for all good patches) (Angstroms)
        bad_patch_x_width   = atof(argv[start_index + 3]);     //x width of the good patches (same for all good patches) (Angstroms)
        bad_patch_y_width   = atof(argv[start_index + 4]);     //y width of the good patches (same for all good patches) (Angstroms)
        theta_good          = atof(argv[start_index + 5]);       //Good patch contact angle
        theta_bad           = atof(argv[start_index + 6]);       //Bad patch contact angle
        d_Rg                = atof(argv[start_index + 7]);       //increments in good patch cap radius (Angstroms i.e. d_Rg=0.05 A)
        d_Rb                = atof(argv[start_index + 8]);       //increments in bad patch cap radius  (Angstroms)
        Rg_max              = atof(argv[start_index + 9]);       //maximum limit of good patch cap.    (Angstroms)
        Rb_max              = atof(argv[start_index + 10]);       //maximum limit of bad patch cap      (Angstroms)
        point_density       = atof(argv[start_index + 11]);       //Total density of points.
        mc_seed[0]          = atoi(argv[start_index + 12]);       //Seed in x-direction for MC
        mc_seed[1]          = atoi(argv[start_index + 13]);       //Seed in y-direction for MC
        mc_seed[2]          = atoi(argv[start_index + 14]);       //Seed in z-direction for MC
        delta               = atof(argv[start_index + 15]);      //buffer region width=(2*\delta) for surface points (Angstroms)
        Rho                 = atof(argv[start_index + 16]);   //Moles per m3
        num_patches         = atoi(argv[start_index + 17]);
        extension_length    = atof(argv[start_index + 18]);
        tag.assign(argv[start_index + 19]);
        starting_box_dim = atof(argv[start_index + 20]);
        last_patch_isinfinite = atoi(argv[start_index + 21]); //Either 0 or 1
        n_max_points = atoi(argv[start_index + 22]);
        starting_patch_is_bad = atoi(argv[start_index + 23]);
        Mu = atof(argv[start_index + 24]);
        Sigma = atof(argv[start_index + 25]);
        T = atof(argv[start_index + 26]);
        printf("Rg, Rb = %10.10f %10.10f\n",Rg, Rb);
        printf("Rg_second, Rb_third = %10.10f %10.10f\n",Rg_second, Rb_third);
        printf("d_Rg, d_Rb = %10.5f %10.5f\n",d_Rg, d_Rb);
        printf("Rg_max, Rb_max = %10.5f %10.5f\n",Rg_max, Rb_max);
        printf("p_width_good = [%10.5f %10.5f] p_width_bad = [%10.5f %10.5f]\n", good_patch_x_width, good_patch_y_width, bad_patch_x_width, bad_patch_y_width);
        printf("delta=%10.5f\n",delta);
        printf("density = %10.5f\n", point_density);
        std::cout<<"tag = "<<tag<<std::endl;
        printf("Rho=%10.10f\n", Rho);
        printf("Mu=%10.10f\n", Mu);
        printf("Sigma=%10.10f\n", Sigma);
        printf("T=%10.10f\n", T);
        
        
        
    }
    
    
    MPI_Bcast (&Rg, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&Rb, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&Rg_second, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&Rb_third, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&good_patch_x_width, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&good_patch_y_width, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&bad_patch_x_width, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&bad_patch_y_width, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&theta_good, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&theta_bad, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&d_Rg, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&d_Rb, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&Rg_max, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&Rb_max, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&point_density, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&mc_seed[0], 3, MPI_INT, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&delta, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&Rho, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&num_patches, 1, MPI_INT, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&extension_length, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&starting_box_dim, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&last_patch_isinfinite, 1, MPI_INT, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&n_max_points, 1, MPI_INT, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&starting_patch_is_bad, 1, MPI_INT, 0,  MPI_COMM_WORLD);
    
    MPI_Bcast (&Mu, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&Sigma, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&T, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    /*Member classes*/
    CheckBoundary check_boundary;
    MC mc_engine;
    DynamicBox dynamic_box;
    Stripes stripes;
    Surface* surface_ptr;
    SurfaceSetup surface_setup (num_patches, z_surface, good_patch_x_width, good_patch_y_width, bad_patch_x_width, bad_patch_y_width, theta_good, theta_bad, last_patch_isinfinite);
    stripes = surface_setup.create_new_stripes(); //This should work because the stripes object already exists in this class and the data returned form create_new_stripes is getting copied into the existing stripes object.
    surface_ptr = &stripes;
    stripes.initial_box(starting_box_dim);
    
    std::vector<std::vector<double> > here_box;
    here_box = (*surface_ptr).box;
    if(myRank==0)
    {printf("box volume via surface ptr inside = %10.5f\n", ((here_box[0][1] - here_box[0][0]) * (here_box[1][1] - here_box[1][0]) * (here_box[2][1] - here_box[2][0])));}
    
    
    n_starting_points = (int) (stripes.box_volume() * point_density);
    dynamic_box = DynamicBox (stripes.box, extension_length);
    if(myRank==0)
    {dynamic_box.print_box();}
        
    
    
    mc_engine =  MC (n_starting_points, stripes.box, mc_seed, MPI_COMM_WORLD, delta);
    check_boundary = CheckBoundary (surface_ptr, &mc_engine , &dynamic_box, point_density);
    
    
    std::vector<double> centre_good{0.0, 0.0, z_surface};
    std::vector<double> centre_left_cap{0.0, 0.0, z_surface};
    std::vector<double> centre_right_cap{0.0, 0.0, z_surface};
    std::vector<double> centre_left_second_cap{0.0, 0.0, z_surface};
    std::vector<double> centre_right_second_cap{0.0, 0.0, z_surface};
    std::vector<double> centre_left_third_cap{0.0, 0.0, z_surface};
    std::vector<double> centre_right_third_cap{0.0, 0.0, z_surface};

    double projected_radius_g, projected_radius_b;
    double projected_radius_g_second, projected_radius_b_third;
    projected_radius_g = Rg*sin(theta_good);
    projected_radius_b = Rb*sin(theta_bad);
    projected_radius_g_second = Rg_second*sin(theta_good);
    projected_radius_b_third = Rb_third*sin(theta_bad);
    
    if(myRank==0)
    {printf("projected_radius_g = %10.5f\t projected_radius_b =%10.5f\n", projected_radius_g, projected_radius_b);
        printf("projected_radius_g_second = %10.5f\t projected_radius_b_third =%10.5f\n", projected_radius_g_second, projected_radius_b_third);}
    
//    double d_ch = sqrt(projected_radius_g*projected_radius_g - (good_patch_x_width/2.0)*(good_patch_x_width/2.0));
//    printf("min Radius=%10.5f\n", (d_ch/sin(theta_bad)));
//
//    double d_B = sqrt(projected_radius_b*projected_radius_b - d_ch*d_ch);
//    printf("d_B=%10.5f\n", d_B);
//
//    centre_left_cap[0] = -(good_patch_x_width/2.0) - d_B ;
//    centre_right_cap[0] = (good_patch_x_width/2.0) + d_B ;
//    printf("(good_patch_x_width/2.0)=%10.5f\t centre_left_cap=%10.5f\t centre_right_cap=%10.5f\n", (good_patch_x_width/2.0), centre_left_cap[0], centre_right_cap[0]);
    
    get_cap_centres(projected_radius_g, projected_radius_b, good_patch_x_width/2.0, good_patch_x_width/2.0, theta_bad, 1, centre_left_cap, centre_right_cap);
    
    double mod_prev_patch_boundary = ((good_patch_x_width/2.0) + bad_patch_x_width);
    double dist_from_boundary = centre_left_cap[0] - (-1*mod_prev_patch_boundary) ;
    
    if(myRank==0)
    {printf("prev_patch_boundary=%10.10f\tdist_from_boundary=%10.10f\n", mod_prev_patch_boundary, dist_from_boundary);}
    get_cap_centres(projected_radius_b, projected_radius_g_second, dist_from_boundary ,  mod_prev_patch_boundary, theta_good, 0,  centre_left_second_cap, centre_right_second_cap);
    
    mod_prev_patch_boundary = ((good_patch_x_width/2.0) + bad_patch_x_width + good_patch_x_width);
    dist_from_boundary = centre_left_second_cap[0] - (-1*mod_prev_patch_boundary) ;
    
    get_cap_centres(projected_radius_g_second, projected_radius_b_third, dist_from_boundary ,  mod_prev_patch_boundary, theta_bad, 1,  centre_left_third_cap, centre_right_third_cap);
    
    
    Spherical_cap GoodCap (centre_good, projected_radius_g, theta_good, z_surface);
    Spherical_cap Cap_left (centre_left_cap, projected_radius_b, theta_bad, z_surface);
    Spherical_cap Cap_right (centre_right_cap, projected_radius_b, theta_bad, z_surface);
    
    Spherical_cap Cap_second_left (centre_left_second_cap, projected_radius_g_second, theta_good, z_surface);
    Spherical_cap Cap_second_right (centre_right_second_cap, projected_radius_g_second, theta_good, z_surface);
    
    Spherical_cap Cap_third_left (centre_left_third_cap, projected_radius_b_third, theta_bad, z_surface);
    Spherical_cap Cap_third_right (centre_right_third_cap, projected_radius_b_third, theta_bad, z_surface);
    
    
    std::vector<Spherical_cap> capsList {GoodCap, Cap_left, Cap_right, Cap_second_left, Cap_second_right, Cap_third_left, Cap_third_right};
    
    Composite_cluster Comp_cluster (capsList, stripes);
    Shape* Cluster_shape_ptr = &Comp_cluster;
    
    check_boundary.ManageBoxBreach(Cluster_shape_ptr);
    
    printf("After manage box breach\n");
    double Volume, SA;
    double N;
    std::vector<double> mc_volume_SA (2, 0.0);
    
    int n_points = mc_engine.get_num_points ();
    printf("n_points = %d n_max_points=%d\n", n_points, n_max_points);
    if(n_points < n_max_points)
    {
        calc_vol_SA_normal(Cluster_shape_ptr, &check_boundary, &mc_engine, mc_volume_SA);
    }
    else
    {
        calc_vol_SA_virtual_points(Cluster_shape_ptr, &check_boundary, &mc_engine, mc_volume_SA);
    }
    
    if(myRank==0)
    {
        Volume = mc_volume_SA[0];
        SA = mc_volume_SA[1];
        N = (Volume * 1e-30) * Rho * avogadro ;
        double dN = 1;
        double rNhere = round_nearN (N , dN);
        std::vector<double> local_projected_SA = Comp_cluster.projected_SAs();
        int diff_size = num_patches - (int)local_projected_SA.size();
        if(diff_size > 0)
        {
            for(int x=0; x<diff_size;x++){local_projected_SA.push_back(0.0);}
            
        }
        else if (diff_size < 0)
        {
            printf("compcluster proj SA size larger than num_patches\n"); abort();
        }
        
        printf("Volume=%10.10f\t SA=%10.10f\t N=%10.5f\t N_rounded=%10.5f\n", Volume, SA, N, rNhere);
        printf("local_projected_SA = [");
        for(size_t i=0; i<local_projected_SA.size(); i++)
        {
            printf("%10.10f\t", local_projected_SA[i]);
            if(isnan(local_projected_SA[i]))
            {
                local_projected_SA[i] = 0.0;
            }
            local_projected_SA[i] = local_projected_SA[i]*1e-20;
        }
            printf("]\n");
        
        double G = free_energy(Rho, Mu, Sigma, Volume*1e-30, SA*1e-20, local_projected_SA, theta_good, theta_bad, T);
        printf("free_energy=%10.10f\n", G);
        
        std::string surfFileName = tag + "_surf_points.txt";
        FILE* SurfpointsFile = fopen(surfFileName.c_str(), "w");
        mc_engine.print_surf_points(Cluster_shape_ptr, SurfpointsFile);
    }
    
    //print_output_variables();

    
    
    MPI_Finalize();
    return EXIT_SUCCESS;
    
}



