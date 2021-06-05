//
//  manager_spherical_caps.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 5/22/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#include "manager_spherical_caps.hpp"

ManagerSphericalCaps::ManagerSphericalCaps():z_surface(0)
{
    /*
     Currently we have have 3 levels of parallelization, where at each level the loop is broken into two.
     1. Rg loop
     2. dB sign i.e. whether bad cap centre is on good patch or bad patch
     3. Rb loop
     NOTE: Each of the loops divide the total available workers from the previous level into 2. So as such at each level n, where the n=1...3 there are 2^n groups of workers available such that the workers in any group of the last level works on the MC volume calculation part.
     */
    MPI_Comm_rank (MPI_COMM_WORLD, &myRank);
    MPI_Comm_size (MPI_COMM_WORLD, &nProcs);
   
}


ManagerSphericalCaps::~ManagerSphericalCaps()
{

    delete []counts;
}

//color_roots(color_roots), roots_size(roots_size), roots_comm(comm)
ManagerSphericalCaps::ManagerSphericalCaps(int branches_per_node, int Levels):
z_surface(0)
{
    printf("Inside ManagerSphericalCaps constructor\n ");
    MPI_Comm_rank (MPI_COMM_WORLD, &myRank);
    MPI_Comm_size (MPI_COMM_WORLD, &nProcs);
    parallel_process = ParallelProcess(branches_per_node, Levels);
}

void ManagerSphericalCaps::read_from_command_line(char* argv[], int start_index)
{
    /* Reading all variables from command line */
    if(parallel_process.is_parent())
    {
        /* Reading all variables from command line */
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
        printf("d_Rg, d_Rb = %10.5f %10.5f\n",d_Rg, d_Rb);
        printf("Rg_max, Rb_max = %10.5f %10.5f\n",Rg_max, Rb_max);
        printf("p_width_good = [%10.5f %10.5f] p_width_bad = [%10.5f %10.5f]\n", good_patch_x_width, good_patch_y_width, bad_patch_x_width, bad_patch_y_width);
        printf("delta=%10.5f\n",delta);
        printf("density = %10.5f\n", point_density);
        
    }
}

void ManagerSphericalCaps::broadcast_data_to_workers()
{
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

}


void ManagerSphericalCaps::setup(char* argv[], int start_index)
{
    read_from_command_line(argv, start_index);
    broadcast_data_to_workers();
    setup_surface();
    setup_box();
    setup_mc_and_boundary();
    setup_input_variables();
    setup_evolve_spherical_caps();
}


void ManagerSphericalCaps::setup_surface()
{
    SurfaceSetup surface_setup (num_patches, z_surface, good_patch_x_width, good_patch_y_width, bad_patch_x_width, bad_patch_y_width, theta_good, theta_bad);
    stripes = surface_setup.create_new_stripes(); //This should work because the stripes object already exists in this class and the data returned form create_new_stripes is getting copied into the existing stripes object.
    surface_ptr = &stripes; 
    stripes.initial_box(starting_box_dim);
    
    if(myRank == 0)
    {
        std::vector<std::vector<double> > here_box;
        here_box = (*surface_ptr).box;
        printf("box volume via surface ptr inside = %10.5f\n", ((here_box[0][1] - here_box[0][0]) * (here_box[1][1] - here_box[1][0]) * (here_box[2][1] - here_box[2][0])));
        
    }
    
}

void ManagerSphericalCaps::setup_box()
{
    //stripes.initial_box(starting_box_dim);
    n_points = (int) (stripes.box_volume() * point_density);
    dynamic_box = DynamicBox (stripes.box, extension_length);
    if(myRank == 0)
    {
        printf("printing box inside box setup\n");
        dynamic_box.print_box();
        
    }
}

void ManagerSphericalCaps::setup_mc_and_boundary()
{
    mc_engine =  MC (n_points, stripes.box, mc_seed, parallel_process.get_last_lvl_branch_comm());
    check_boundary = CheckBoundary (surface_ptr, &mc_engine , &dynamic_box, point_density);
    if(myRank == 0)
    {printf("mc n_points inside=%d\n", mc_engine.get_num_points());
        printf("check boundary point density inside=%10.10f\n", check_boundary.get_point_density());}
    
}

void ManagerSphericalCaps::setup_input_variables()
{
    input_variables.Rg_max = Rg_max;
    input_variables.Rb_max = Rb_max;
    input_variables.d_Rg = d_Rg;
    input_variables.d_Rb = d_Rb;
    input_variables.theta_good = theta_good;
    input_variables.theta_bad = theta_bad;
    input_variables.delta = delta;
    input_variables.Rho = Rho;
}

void ManagerSphericalCaps::setup_evolve_spherical_caps()
{
    
    evolve_spherical_cap = EvolveSphericalCap(parallel_process, stripes, &output_variables, &input_variables , check_boundary, mc_engine, dynamic_box, num_patches);
    evolve_spherical_cap.print_initial_variables();
}


//evolve
void ManagerSphericalCaps::evolve()
{
    evolve_spherical_cap.evolve();
}


// Printing functions
void ManagerSphericalCaps::print_surface_ptr()
{
    if(myRank == 0)
    {
        std::vector<std::vector<double> > here_box;
        here_box = (*surface_ptr).box;
        printf("box volume via surface ptr outside = %10.5f\n", ((here_box[0][1] - here_box[0][0]) * (here_box[1][1] - here_box[1][0]) * (here_box[2][1] - here_box[2][0])));
        
    }
}
void ManagerSphericalCaps::print_box()
{
    if(myRank == 0)
    {
    printf("printing box outside\n");
        dynamic_box.print_box();}
}
void ManagerSphericalCaps::print_mc_and_check_boundary()
{
    if(myRank == 0)
    {
    printf("mc n_points outside=%d\n", mc_engine.get_num_points());
        printf("check boundary point density outside=%10.10f\n", check_boundary.get_point_density());}
}


void ManagerSphericalCaps::gather()
{
    
//    if(color_roots == 0)
//    {
//        /*These counts are true for single element quantities per grid point like N, Volume, SA, dB and dG*/
//        counts = new int[roots_size];
//        MPI_Gather(&nelements, 1, MPI_INT, counts, 1, MPI_INT, 0, roots_comm);
//    }
}


void ManagerSphericalCaps::print_quants(int rank)
{
    if(myRank==rank)
    {
        std::cout<<" good_patch_x_width = "<<good_patch_x_width<<std::endl;
        std::cout<<" good_patch_y_width = "<<good_patch_y_width<<std::endl;
        std::cout<<" bad_patch_x_width = "<<bad_patch_x_width<<std::endl;
        std::cout<<" bad_patch_y_width = "<<bad_patch_y_width<<std::endl;
        std::cout<<" theta_good = "<<theta_good<<std::endl;
        std::cout<<" theta_bad = "<<theta_bad<<std::endl;
        std::cout<<" d_Rg = "<<d_Rg<<std::endl;
        std::cout<<" d_Rb = "<<d_Rb<<std::endl;
        std::cout<<" Rg_max = "<<Rg_max<<std::endl;
        std::cout<<" Rb_max = "<<Rb_max<<std::endl;
        std::cout<<" point_density = "<<point_density<<std::endl;
        std::cout<<" mc_seed[0] = "<<mc_seed[0]<<std::endl;
        std::cout<<" mc_seed[1] = "<<mc_seed[1]<<std::endl;
        std::cout<<" mc_seed[2] = "<<mc_seed[2]<<std::endl;
        std::cout<<" delta = "<<delta<<std::endl;
        std::cout<<" Rho = "<<Rho<<std::endl;
        std::cout<<" num_patches = "<<num_patches<<std::endl;
        std::cout<<" extension_length = "<<extension_length<<std::endl;
        std::cout<<" starting_box_dim = "<<starting_box_dim<<std::endl;

    }
}



//MPI_Init (&argc, &argv);
//MPI_Comm_rank (MPI_COMM_WORLD, &myRank);
//MPI_Comm_size (MPI_COMM_WORLD, &nProcs);
//
//if(myRank == 0)
//{printf("nProcs = %d\n", nProcs);}
//
//int BranchSize;
//int levels;
//if(myRank == 0)
//{
//    /* Reading BranchSize from command line */
//    BranchSize  = atoi(argv[1]);
//    levels = atoi(argv[2]);
//}
//MPI_Bcast (&BranchSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
//MPI_Bcast (&levels, 1, MPI_INT, 0, MPI_COMM_WORLD);
//
//std::vector<int> worker_colors_per_level (levels); //This will have the color of each worker at each level.
//std::vector<int> num_workers_per_groups_per_level (levels);
//for(int i=0; i < levels; i++)
//{
//    int group_size = pow(BranchSize, i+1);
//    num_workers_per_groups_per_level[i] = nProcs/group_size; //Ensure that they are multiples
//    worker_colors_per_level[i] = myRank/num_workers_per_groups_per_level[i];
//    printf("level =%d\t rank = %d\t worker_colors_per_level = %d\n", i+1, myRank, worker_colors_per_level[i]);
//}
//
//std::vector<MPI_Comm> branch_comm_levels(levels);
//std::vector<int> level_branch_rank(levels);
//std::vector<int> level_branch_size(levels);
//for(int i=0; i < levels; i++)
//{
//    MPI_Comm_split(MPI_COMM_WORLD, worker_colors_per_level[i], myRank, &branch_comm_levels[i]); //color
//
//    MPI_Comm_rank(branch_comm_levels[i], &level_branch_rank[i]);
//    MPI_Comm_size(branch_comm_levels[i], &level_branch_size[i]);
//    printf("level = %d\t original rank=%d\t branch_rank=%d branch_size=%d\n", i, myRank, level_branch_rank[i], level_branch_size[i]);
//}
//
////MPI_Comm_split(MPI_COMM_WORLD, worker_colors_per_level[levels-1], myRank, &branch_comm); //color
//
//
///* Creating another communicator containing the roots of each communicator so that V, SA..etc, data can be transferred from them to the original root */
////int color_roots = ( myRank%BranchSize == 0) ? myRank%BranchSize : 1;
//int color_roots = ( myRank%num_workers_per_groups_per_level[levels-1] == 0) ? 0 : 1 ;
//MPI_Comm roots_comm;
//MPI_Comm_split(MPI_COMM_WORLD, color_roots, myRank, &roots_comm);
//int roots_rank, roots_size;
//MPI_Comm_rank(roots_comm, &roots_rank);
//MPI_Comm_size(roots_comm, &roots_size);
//
//
//
//
///* This is the older parallelisation scheme where I was dividing only the dB loop into two by making two branches.
// */
//
//int color = myRank/BranchSize;
//MPI_Comm branch_comm;
//MPI_Comm_split(MPI_COMM_WORLD, color, myRank, &branch_comm);
//int branch_rank, branch_size;
//MPI_Comm_rank(branch_comm, &branch_rank);
//MPI_Comm_size(branch_comm, &branch_size);
////printf("original rank=%d\t branch_rank=%d branch_size=%d\n", myRank, branch_rank, branch_size);
//
///* Creating another communicator containing the roots of each communicator so that V, SA..etc. data can be transferred from them to the original root */
//int color_roots = ( myRank%BranchSize == 0) ? myRank%BranchSize : 1;
//MPI_Comm roots_comm;
//MPI_Comm_split(MPI_COMM_WORLD, color_roots, myRank, &roots_comm);
//int roots_rank, roots_size;
//MPI_Comm_rank(roots_comm, &roots_rank);
//MPI_Comm_size(roots_comm, &roots_size);
