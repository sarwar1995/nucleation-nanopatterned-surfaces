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
ManagerSphericalCaps::ManagerSphericalCaps(int branch_size, int branches_per_node, int Levels):
z_surface(0)
{
    parallel_process = ParallelProcess(branch_size, branches_per_node, Levels);
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
    init_cap_identifier();
    setup_surface();
    setup_box();
    setup_mc_and_boundary();
}


void ManagerSphericalCaps::setup_surface()
{
    SurfaceSetup surface_setup (num_patches, z_surface, good_patch_x_width, good_patch_y_width, bad_patch_x_width, bad_patch_y_width, theta_good, theta_bad);
    stripes = surface_setup.create_new_stripes(); //This should work because the stripes object already exists in this class and the data returned form create_new_stripes is getting copied into the existing stripes object.
    surface_ptr = &stripes;
}

void ManagerSphericalCaps::setup_box()
{
    stripes.initial_box(starting_box_dim);
    n_points = (int) (stripes.box_volume() * point_density);
    dynamic_box = DynamicBox (stripes.box, extension_length);
}

void ManagerSphericalCaps::setup_mc_and_boundary()
{
    mc_engine =  MC (n_points, stripes.box, mc_seed, parallel_process.get_last_lvl_branch_comm());
    check_boundary = CheckBoundary (surface_ptr, &mc_engine , &dynamic_box, point_density);
}


void ManagerSphericalCaps::setup_output_variables()
{
    int n_unique_patches = stripes.get_n_unique_patches();
    projected_SA.resize(n_unique_patches) ;
    clstr_centre_location_modifier_global_arrays.resize(n_unique_patches); //{+1, -1}
    
    radii_global_arrays.resize(n_unique_patches);
    
    centre_left_cap[2] = z_surface;
    centre_right_cap[2] = z_surface;
}

void ManagerSphericalCaps::evolve_three_caps()
{
    double Rg;
    if(num_patches != 3)
    {
        printf("Can't use evolve_three_caps with more than three patches\n");
        abort();
    }
    else
    {
        std::vector<int> Rg_loop_start_end = getLoopStartEnd (len_Rg, worker_colors_per_level[0]);
        printf("rank = %d Rg loop start end = %d %d %d\n", myRank, (int)Rg_loop_start_end.size(), Rg_loop_start_end[0], Rg_loop_start_end[1]);
        
        /* Starting the loop for Rg.*/
        for(int i = Rg_loop_start_end[0] ; i <= Rg_loop_start_end[1]; i++) // len_Rg
        {
            cap_identifier = 0;
            cap_type = CENTRE_GOOD;
            Rg = 0.0 + i*d_Rg ;  //  startingRg + i*d_Rg            //sphere's radius
            double projected_rg = Rg*sin(theta_good); //projected radii of the circles
            if(level_branch_rank[0] == 0)
            {
                printf("i=%d\t Rg = %10.10f\tprojected_rg = %10.10f\n",i, Rg, projected_rg);
            }
            
            Spherical_cap GoodCap (centre_good, projected_rg, theta_good, z_surface);
            Cluster_shape_ptr= &GoodCap;
            check_boundary_ptr->ManageBoxBreach(Cluster_shape_ptr);
            if(check_breaking_condition())
            {
                break;
            }
            else if(check_bounds())
            {
                evolve_bad_cap(projected_rg);
            }
            else
            {
                calc_volume_SA();
            }
        }
    }
}

//void ManagerSphericalCaps::evolve_bad_cap(double projected_rg)
//{
//    int Color_per_group = worker_colors_per_level[1] % 2;
//    double projected_rb_min = sqrt(projected_rg*projected_rg - ((good_patch_x_width*good_patch_x_width)/4.0));
//    double Rb_min = projected_rb_min/(double)sin(theta_bad) ; //Be careful to not have a theta that is 0 or 180
//    int len_Rb = (int) ((Rb_max-Rb_min)/d_Rb) + 1;
//    std::vector<int> Rb_loop_start_end = getLoopStartEnd (len_Rb, worker_colors_per_level[2]);
//    int dB_sign;
//    double d_B;
//    if (Color_per_group == 0)
//    {
//        dB_sign = -1;
//    }
//    else if (Color_per_group == 1)
//    {
//        dB_sign = 1;
//    }
//
//    for(int l=Rb_loop_start_end[0]; l <= Rb_loop_start_end[1]; l++) //len_Rb //Adding another level for Rb parallelization
//    {
//        Rb = Rb_min + l * d_Rb ;
//        double projected_rb = Rb * sin(theta_cb);
//        double inside_sq = (projected_rb*projected_rb) - (projected_rg*projected_rg) + ((good_patch_x_width*good_patch_x_width)/4.0);
//
//        if(abs(inside_sq)<1e-10 && inside_sq < 0.0) {inside_sq = 0.0;}
//        d_B = dB_list * sqrt(inside_sq) ;
//        /* Symmetric caps on either side of central patch */
//        c_bad_left[0] = -(good_patch_x_width/2.0) - d_B;
//        c_bad_right[0] = (good_patch_x_width/2.0) + d_B;
//    }
//
//}

void ManagerSphericalCaps::evolve_cap()
{
    double R, R_min, d_R;
    double projected_r;
    double theta;
    int len_R;
    std::vector<int> R_loop_start_end;
    cap_type = (cap_identifier%2 == 0) ? GOOD : BAD ;
    theta = (cap_type==GOOD) ? theta_good : theta_bad;
    if (cap_identifier == 0)
    {
        cap_type = CENTRE_GOOD;
        len_R = (int) ((Rg_max-0.0)/d_Rg) + 1;
        R_min = 0.0;
        d_R = d_Rg;
    }
    else
    {

    }

    R_loop_start_end = getLoopStartEnd (len_R, worker_colors_per_level[0]);
    printf("rank = %d R loop start end = %d %d %d\n", myRank, (int)R_loop_start_end.size(), R_loop_start_end[0], R_loop_start_end[1]);
    for(int i = R_loop_start_end[0] ; i <= R_loop_start_end[1]; i++) // len_Rg
    {

        R = R_min + i*d_R ;  //  startingRg + i*d_Rg            //sphere's radius
        projected_r = R * sin(theta);
        if(level_branch_rank[0] == 0)
        {
            printf("i=%d\t Rg = %10.10f\tprojected_rg = %10.10f\n",i, R, projected_r);
        }
        if (cap_identifier == 0)
        {
            //Is this valid. ????
            Spherical_cap GoodCap (centre_good, projected_rg, theta_good, z_surface);

            Cluster_shape_ptr= &GoodCap;
        }
        else
        {

        }

        check_boundary_ptr->ManageBoxBreach(Cluster_shape_ptr);

        check_bounds();
    }
}

bool ManagerSphericalCaps::check_bounds()
{
    //returning false indicates bounds are NOT crossed and true indicates they ARE.
    stripes_bounds = surface_ptr->monitor_cluster_spread(Cluster_shape_ptr);
    if(cap_type == CENTRE_GOOD)
    {
        if(stripes_bounds[0] == 0)
        {return false;}
        else
        {return true;}
    }
    else
    {
        int stripes_index = 2*cap_identifier - 1;
        if (stripes_bounds [stripes_index] != stripes_bounds [stripes_index+1])
        {
            printf("Patch with identifier %d not crossed symmetrically\n", cap_identifier);
            abort();
        }
        else if (stripes_bounds [stripes_index] == 0 && stripes_bounds[stripes_index+1] == 0)
        {
            return false;
        }
        else  if (stripes_bounds [stripes_index] == 1 && stripes_bounds[stripes_index+1] == 1)
        {
            return true;
        }
    }
    
}


bool ManagerSphericalCaps::check_breaking_condition()
{
    //returning false indicates breaking will NOT occur and true indicates that it WILL.
    if (surface_ptr->surface_bounds_breach (Cluster_shape_ptr)) //This is checking the final
    {
        return true;
    }
    else if(cap_identifier == 0)
    {
        return false;
    }
    else
    {
        std::vector<double> current_centre = cap_centres[cap_identifier];
        double current_proj_radii = cap_projected_radii[cap_identifier];
        std::vector<double> previous_centre = cap_centres[cap_identifier-1];
        double previous_proj_radii = cap_projected_radii[cap_identifier-1];
        
        //The abs() value of centre differences accounts for both left and right images of the current cap.
        bool check_less_than = current_proj_radii < previous_proj_radii + abs(current_centre[0] - previous_centre[0]) ;
        bool check_greater_than = current_proj_radii > previous_proj_radii - abs(current_centre[0] - previous_centre[0]) ;
    
        if(!(check_less_than && check_greater_than))
        {
            return true;
        }
        else {return false;}
    }
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


void ManagerSphericalCaps::print_nelements()
{
//    if(myRank==0)
//    {
//        for(int i=0; i<roots_size; ++i)
//        {
//            printf("counts[%d] = %d\n", i, counts[i]);
//        }
//    }
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
