//
//  main_spherocylinder_parallel.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 3/28/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#include <stdio.h>
#include <vector>
#include "Spherocylinder_cap_composite.hpp"
#include "Patch.hpp"
#include "Stripes.hpp"
#include "DynamicBox.hpp"
#include "MC_parallel.hpp"
#include "CheckBoundary.hpp"
#include "FreeEnergy.hpp"
#include "EvolveCluster.hpp"
#include "VolumeSA_calculations.hpp"
#include <mpi.h>

using namespace std;

Patch setup_patch(double width_x, double width_y, double x_coordinate_centre, double theta, double z_surface)
{
    std::vector<double> centre(3);   //Centre of the good patch. Only x-y
    centre[0] = x_coordinate_centre;
    centre[1] = 0.0;
    centre[2] = z_surface;
    std::vector<double> dim(2,0.0);
    
    dim[0] = width_x; dim[1] = width_y; //Common dimensions for all Good patches
    
    Patch patch (theta, centre, dim);
    return patch ;
}

int myRank, nProcs;
int main(int argc, char * argv[]) {
    
    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &myRank);
    MPI_Comm_size (MPI_COMM_WORLD, &nProcs);
    
    if(myRank == 0)
        {printf("nProcs = %d\n", nProcs);}
    
    /* Creating a group of processes for each branch, for dB_sign = {+1, -1}. Each branch's leader is responsible for checking the boundaries and managing the box.
     BranchSize is taken as an input.
     This MPI_Comm_split uses color to assign a communicator to each process and their original rank as the key to assign a rank within the communicator */
    
    int BranchSize;
    if(myRank == 0)
    {
        /* Reading BranchSize from command line */
        BranchSize  = atof(argv[1]);
    }
    MPI_Bcast (&BranchSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    int color = myRank/BranchSize;
    printf("rank = %d\t color = %d\n", myRank, color);
    MPI_Comm branch_comm;
    MPI_Comm_split(MPI_COMM_WORLD, color, myRank, &branch_comm);
    int branch_rank, branch_size;
    MPI_Comm_rank(branch_comm, &branch_rank);
    MPI_Comm_size(branch_comm, &branch_size);
    printf("original rank=%d\t branch_rank=%d branch_size=%d\n", myRank, branch_rank, branch_size);
    
    /* Creating another communicator containing the roots of each communicator so that V, SA..etc, data can be transferred from them to the original root */
    int color_roots = ( myRank%BranchSize == 0) ? myRank%BranchSize : 1;
    MPI_Comm roots_comm;
    MPI_Comm_split(MPI_COMM_WORLD, color_roots, myRank, &roots_comm);
    int roots_rank, roots_size;
    MPI_Comm_rank(roots_comm, &roots_rank);
    MPI_Comm_size(roots_comm, &roots_size);

    //Density
    double Rho ;
    
    //Variables for spherocylinder
    double d_cyl_length, d_chord_length;
    double cyl_max; //Maximum length of the cylinder
    
    //Variables for caps
    double d_Rg, d_Rb;
    double Rg_max, Rb_max;
    
    //Variables for Patches
    double z_surface = 0.0; //Unless stated otherwise
    
    double pG_width_x, pG_width_y, pB_width_x, pB_width_y, theta_cg, theta_cb;
    int num_patches;       //This is the maximum number of patches that are allowed.
    
    // Monte-Carlo related variables
    double point_density; //Constant point density for all additional points added dynamically
    int n_points;       //Initial number of points in the first good patch box calcualted from the density
    int aSeed[3];
    double delta;
    
    //Dynamic box related variables
    double extension_length; //The fixed length for box extension in each direction
    std::string tag;
    
    FILE* outputfile;
    FILE* output_points_file;
    
    if(myRank == 0)
    {
        /* Reading all variables from command line */
        pG_width_x  = atof(argv[2]);     //x width of the good patches (same for all good patches) (Angstroms)
        pG_width_y  = atof(argv[3]);     //y width of the good patches (same for all good patches) (Angstroms)
        pB_width_x  = atof(argv[4]);     //x width of the good patches (same for all good patches) (Angstroms)
        pB_width_y  = atof(argv[5]);     //y width of the good patches (same for all good patches) (Angstroms)
        theta_cg    = atof(argv[6]);       //Good patch contact angle
        theta_cb    = atof(argv[7]);       //Bad patch contact angle
        d_Rg        = atof(argv[8]);       //increments in good patch cap radius (Angstroms i.e. d_Rg=0.05 A)
        d_Rb        = atof(argv[9]);       //increments in bad patch cap radius  (Angstroms)
        d_cyl_length= atof(argv[10]);
        d_chord_length= atof(argv[11]);
        Rg_max      = atof(argv[12]);       //maximum limit of good patch cap.    (Angstroms)
        Rb_max      = atof(argv[13]);
        cyl_max     = atof(argv[14]);        //maximum limit of bad patch cap      (Angstroms)
        point_density = atof(argv[15]);     // Constant point density for MC. (points/Angstrom^3)
        aSeed[0]    = atoi(argv[16]);       //Seed in x-direction for MC
        aSeed[1]    = atoi(argv[17]);       //Seed in y-direction for MC
        aSeed[2]    = atoi(argv[18]);       //Seed in z-direction for MC
        delta       = atof(argv[19]);      //buffer region width=(2*\delta) for surface points (Angstroms)
        Rho         = atof(argv[20]);   //Moles per m3
        num_patches = atoi(argv[21]);
        extension_length = atof(argv[22]);
        tag.assign(argv[23]);
        printf("d_Rg, d_Rb = %10.5f %10.5f\n",d_Rg, d_Rb);
        printf("point_density = %10.10f\n",point_density);
        printf("Rg_max, Rb_max = %10.5f %10.5f\n",Rg_max, Rb_max);
        printf("p_width_good = [%10.5f %10.5f] p_width_bad = [%10.5f %10.5f]\n", pG_width_x, pG_width_y, pB_width_x, pB_width_y);
        printf("delta=%10.5f\n",delta);
    }
    MPI_Bcast (&pG_width_x, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&pG_width_y, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&pB_width_x, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&pB_width_y, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&theta_cg, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&theta_cb, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&d_Rg, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&d_Rb, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&d_cyl_length, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&d_chord_length, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    
    MPI_Bcast (&Rg_max, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&Rb_max, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&cyl_max, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    
    MPI_Bcast (&point_density, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&aSeed[0], 3, MPI_INT, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&delta, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&Rho, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&num_patches, 1, MPI_INT, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&extension_length, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    
    /* I/O part handled by rank = 0 only. */
    if(myRank == 0)
    {
//        std::string output_points_file_name = tag + "_points.txt";
        std::string outputfile_name = tag + "_Volume_SA_data.txt" ;
//        output_points_file = fopen(output_points_file_name.c_str(), "w");
        outputfile = fopen(outputfile_name.c_str(), "w");
    }
    
    Patch good_patch = setup_patch(pG_width_x, pG_width_y, 0.0, theta_cg, z_surface);
    std::vector<double> centre_good ({0.0, 0.0, z_surface});
    std::vector<Patch> list_of_patches {good_patch};
    std::vector<std::vector<double> > orientations {centre_good};
    
    Surface* surface_ptr;
    Stripes stripes (list_of_patches, orientations, z_surface);
    surface_ptr = &stripes;
    
    //Here starting with a much smaller box to check if dynamic box works.
    double box_height = (Rg_max >= Rb_max) ? (double)(Rg_max/6.0) : (double)(Rb_max/6.0) ;
    // ? 2.0*Rg_max : 2.0*Rb_max ;
    
    stripes.calc_box(box_height);
    double initial_box_volume = stripes.box_volume();
    n_points = (int) (initial_box_volume * point_density);
    
    if(myRank == 0)
    {
        printf("n_points = %d\n", n_points);
        for(size_t i =0; i<stripes.box.size(); i++)
        {
            for(size_t j =0; j<stripes.box[i].size(); j++)
            {
                printf("box[%d][%d] = %10.5f\t",(int)i,(int)j,stripes.box[i][j]);
            }
        }
    }
    
    std::vector<double> mc_volume_SA(2,0.0);
    //Instantiating dynamic_box.
    // Two separate dynamic box instances are instantiated for the two branches. 
    DynamicBox dynamic_box (stripes.box, extension_length);
    
    MC mc_engine (n_points, stripes.box, aSeed, branch_comm);
    MC* mc_engine_ptr = &mc_engine;
    
    CheckBoundary check_boundary (surface_ptr, mc_engine_ptr , &dynamic_box, point_density);
    CheckBoundary* check_boundary_ptr = &check_boundary ;
    Shape* Cluster_shape_ptr;
    
    std::vector<double> maximum_limits ({Rg_max, Rb_max , cyl_max});
    std::vector<double> increments({d_Rg, d_Rb , d_cyl_length});
    std::vector<double> patch_widths({pG_width_x, pG_width_y, pB_width_x, pB_width_y});
    
    //EvolveCluster evolve_cluster (check_boundary_ptr,  mc_engine_ptr, surface_ptr, mc_volume_SA, maximum_limits , increments, theta_cg,  theta_cb, patch_widths , z_surface, delta, Rho);
    
    std::vector<int> dB_array;
    std::vector<double> Number_particles, Volume_array, SA_array, projected_SA_array;
    std::vector<double> Rg_array, Rb_array, cyl_length_array, chord_length_array;
    
    if (color_roots==0){printf("Volume_array size before anything = %d\n", (int)Volume_array.size());}
    EvolveCluster evolve_cluster(check_boundary_ptr,  mc_engine_ptr, surface_ptr, mc_volume_SA, maximum_limits , increments, theta_cg,  theta_cb, patch_widths , z_surface, delta, Rho, myRank, color, color_roots, dB_array, Number_particles, Volume_array, SA_array,  projected_SA_array,  Rg_array, Rb_array, cyl_length_array, chord_length_array);
    
    /* Setting up required variables and output data structures*/
    //Rg, cyl_length, chord_length, Rb_here, dB_sign_here, N, Volume, SA, projected_SA
    
 
    
    double Rg, Rb;
    double cyl_length, chord_length;
    double N, Volume;
    double SA;
    double projected_SA;
    
    std::vector<int> stripes_bounds; //array of (0,1) to check crossing of boundaries
    std::vector<int> box_breach; //array of (0,1) to check box surface breach
    
    int len_Rg = (int) ((Rg_max-0.0)/d_Rg) + 1;
    int len_cyl = (int) ((cyl_max-0.0)/d_cyl_length) + 1;
    for(int i = 0 ; i< len_Rg; i++) // 1
    {
    
        Rg = 0.0 + i*d_Rg ;

        if(myRank == 0)
        {
            printf("Rg = %10.10f\n", Rg);
            
        }
        
        cyl_length = 0.0;
        chord_length = 0.0;
        
        double projected_rg = Rg*sin(theta_cg);
        Spherical_cap GoodCap (centre_good, projected_rg, theta_cg, z_surface);
        Cluster_shape_ptr= &GoodCap;
        check_boundary.ManageBoxBreach(Cluster_shape_ptr);
        
        std::vector<double> shape_xy_spread = Cluster_shape_ptr->xy_spread();
        std::vector<double> cPatchBounds = good_patch.patch_boundaries();
        
        if(cPatchBounds[0] < shape_xy_spread[0] && cPatchBounds[1] > shape_xy_spread[1]) //Only spherical cap within the boundary
        {
            if((cPatchBounds[0] > shape_xy_spread[0] - d_Rg*sin(theta_cg)) && (cPatchBounds[0] <= shape_xy_spread[0]) && (cPatchBounds[1] < shape_xy_spread[1] + d_Rg*sin(theta_cg)) && (cPatchBounds[1] >= shape_xy_spread[1])  )
            {
                if(myRank == 0)
                {printf("Inside spherocylinder patch \n");}
                for(int j = 0 ; j< len_cyl; j++) //j=0 is just the spherical cap that touches the boundary
                {
                    cyl_length = 0.0 + j*d_cyl_length;
                    
//                    if(myRank == 0)
//                    {printf("cyl_length = %10.10f\n", cyl_length);}
                    std::vector<double> normal_to_the_circle({0.0 , 1.0, 0.0});
                    
                    SpheroCylinder sphero_cylinder (GoodCap, cyl_length, normal_to_the_circle, z_surface);
                    
//                    printf("Volume_array size before evolve cluster = %d\n",(int)Volume_array.size());
                    evolve_cluster.EvolveBadCapWithSpherocylinder_parallel (&sphero_cylinder, Cluster_shape_ptr, cyl_length, d_chord_length, Rg, outputfile, output_points_file);
                }
            }
            else
            {
                if(color == 0) //myRank == 0
                {
                //Only spherical cap
                    mc_volume_SA = mc_engine.calc_volume_SA(Cluster_shape_ptr, delta);
                    Volume = mc_volume_SA[0]; //GoodCap.getVolume();
                    SA = mc_volume_SA[1];
                    if(color_roots == 0)
                    {
                            projected_SA = GoodCap.projected_SA();
                            N = (Volume * 1e-30) * Rho * avogadro ;
                            double Rb_here = 0.0;
                            int dB_sign_here = 0;

                            Rg_array.push_back(Rg);
                        
                            cyl_length_array.push_back(cyl_length);
                        
                            chord_length_array.push_back(chord_length);
                        
                            Rb_array.push_back(Rb_here);
                        
                            dB_array.push_back(dB_sign_here);
                        
                            Number_particles.push_back(N);
                        
                            Volume_array.push_back(Volume);
                        
                            SA_array.push_back(SA);
                        
                            //printf("Volume_array size = %d SA_array = %d Rg_array size=%d only color 0\n", (int)Volume_array.size(), (int)SA_array.size(), (int) Rg_array.size());
                        
                            std::vector<double> projected_SA_array_here ({0.0, projected_SA, 0.0});
                            addQuant(projected_SA_array, projected_SA_array_here, (int)projected_SA_array_here.size());
                            
                            //fprintf(outputfile, "%10.10f\t%10.10f\t%10.10f\t%10.10f\t%d\t%10.10f\t%10.10f\t%10.10f\t%10.10f\n" , Rg, cyl_length, chord_length, Rb_here, dB_sign_here, N, Volume, SA, projected_SA);
                    }
                    
                }
                
            }
        }
        else if (cPatchBounds[0] > shape_xy_spread[0] && cPatchBounds[1] < shape_xy_spread[1])
        {
            //Not allowing rg to grow beyond the good patch boundary
            break;
        }
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    if(color_roots == 0)
    {
        int counts[roots_size];
        int projected_SA_counts[roots_size];
        int radii_counts[roots_size];
        
        int nelements = (int)Volume_array.size();
        int projected_SA_nelements = (int)projected_SA_array.size();
        
       // printf("process=%d\t root_rank=%d roots_size=%d Volume=%d, SA=%d, N=%d, dB=%d, proj_SA=%d and Rg=%d elements\n", myRank, roots_rank, roots_size, (int)Volume_array.size(), (int)SA_array.size(), (int)Number_particles.size(), (int)dB_array.size(),  (int)projected_SA_array.size(), (int)Rg_array.size());
        
        MPI_Gather(&nelements, 1, MPI_INT, &counts[0], 1, MPI_INT, 0, roots_comm);
        MPI_Gather(&projected_SA_nelements, 1, MPI_INT, &projected_SA_counts[0], 1, MPI_INT, 0, roots_comm);
        
        if(myRank==0){printf("After MPIGather\n");}
        
        int disps[roots_size];
        int projected_SA_disps[roots_size];
        // Displacement for the first chunk of data - 0
        for (int i = 0; i < roots_size; i++)
        {
            //printf("i = %d\t counts = %d\t projected_SA_counts=%d\t radii_counts=%d\n", i, counts[i], projected_SA_counts[i], radii_counts[i]);
            disps[i] = (i > 0) ? (disps[i - 1] + counts[i - 1]) : 0;

            projected_SA_disps[i] = (i > 0) ? (projected_SA_disps[i - 1] + projected_SA_counts[i - 1]) : 0;

            //printf("disps = %d\t projected_SA_disps=%d\n", disps[i], projected_SA_disps[i]);
        }
        std::vector<int> dB_global_array;
        std::vector<double> Number_particles_global, Volume_global_array, SA_global_array, projected_global_SA_array;
        std::vector<double> Rg_global_array, Rb_global_array, cyl_length_global_array, chord_length_global_array;


        if(myRank == 0)
        {
            int V_SA_size = disps[roots_size - 1] + counts[roots_size - 1];
            int proj_SA_size = projected_SA_disps[roots_size - 1] + projected_SA_counts[roots_size - 1];

            dB_global_array.resize(V_SA_size);
            Number_particles_global.resize(V_SA_size);
            Volume_global_array.resize(V_SA_size);
            SA_global_array.resize(V_SA_size);
            Rg_global_array.resize(V_SA_size);
            Rb_global_array.resize(V_SA_size);
            cyl_length_global_array.resize(V_SA_size);
            chord_length_global_array.resize(V_SA_size);

            projected_global_SA_array.resize(proj_SA_size);

        }

        MPI_Gatherv(Volume_array.data(), nelements, MPI_DOUBLE, Volume_global_array.data(), &counts[0], &disps[0], MPI_DOUBLE, 0, roots_comm);

        MPI_Gatherv(SA_array.data(), nelements, MPI_DOUBLE, SA_global_array.data(), &counts[0], &disps[0], MPI_DOUBLE, 0, roots_comm);

        MPI_Gatherv(Number_particles.data(), nelements, MPI_DOUBLE, Number_particles_global.data(), &counts[0], &disps[0], MPI_DOUBLE, 0, roots_comm);

        MPI_Gatherv(Rg_array.data(), nelements, MPI_DOUBLE, Rg_global_array.data(), &counts[0], &disps[0], MPI_DOUBLE, 0, roots_comm);

        MPI_Gatherv(Rb_array.data(), nelements, MPI_DOUBLE, Rb_global_array.data(), &counts[0], &disps[0], MPI_DOUBLE, 0, roots_comm);

        MPI_Gatherv(cyl_length_array.data(), nelements, MPI_DOUBLE, cyl_length_global_array.data(), &counts[0], &disps[0], MPI_DOUBLE, 0, roots_comm);

        MPI_Gatherv(chord_length_array.data(), nelements, MPI_DOUBLE, chord_length_global_array.data(), &counts[0], &disps[0], MPI_DOUBLE, 0, roots_comm);

        MPI_Gatherv(dB_array.data(), nelements, MPI_INT, dB_global_array.data(), &counts[0], &disps[0], MPI_INT, 0, roots_comm);

        MPI_Gatherv(projected_SA_array.data(), projected_SA_nelements, MPI_DOUBLE, projected_global_SA_array.data(), &projected_SA_counts[0], &projected_SA_disps[0], MPI_DOUBLE, 0, roots_comm);

        if(myRank == 0)
        {
            if(outputfile != NULL)
            {
                add_Volume_SA_spherocylinder_parallel(dB_global_array,  Number_particles_global,  Volume_global_array,  SA_global_array,  projected_global_SA_array, Rg_global_array, Rb_global_array, cyl_length_global_array,  chord_length_global_array, outputfile);
                printf("After adding quants\n");
            }
            else
            {printf("File pointer is invalid\n") ;}

        }
    }
    
    
    
//
//    if(myRank == 0)
//    {
//        if(outputfile != NULL)
//        {
//            fclose(outputfile);
//            printf("After closing file\n");
//
//        }
//
//    }
    
    MPI_Comm_free(&branch_comm);
    MPI_Finalize();
    return EXIT_SUCCESS;

}
