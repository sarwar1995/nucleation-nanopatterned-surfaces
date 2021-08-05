//
//  evolve_spherical_caps.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 5/28/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#include "evolve_spherical_caps.hpp"

EvolveSphericalCap::EvolveSphericalCap()
{
    MPI_Comm_rank (MPI_COMM_WORLD, &myRank);
    MPI_Comm_size (MPI_COMM_WORLD, &nProcs);
}

EvolveSphericalCap::~EvolveSphericalCap()
{
    
}

EvolveSphericalCap::EvolveSphericalCap(ParallelProcess* process, Stripes stripes, SphericalCapOutput* output_vars, SphericalCapInput* input_vars, CheckBoundary* check_bounds, MC* mc, int n_patches, PeriodicIO* periodic, int starting_bad_patch):
parallel_process(process), stripes(stripes), check_boundary(check_bounds), mc_engine(mc), num_patches(n_patches), output_variables(output_vars), dummy_decider(0), periodic_io(periodic), starting_patch_is_bad(starting_bad_patch)
{
    MPI_Comm_rank (MPI_COMM_WORLD, &myRank);
    MPI_Comm_size (MPI_COMM_WORLD, &nProcs);
    
    int n_unique_patches = ((num_patches-1)/2) + 1; //adding 1 for central patch.
    current_radii_array.resize(n_unique_patches , 0.0);
    current_projected_radii_array.resize(n_unique_patches , 0.0);
    current_clstr_centre_modifier_array.resize(n_unique_patches , 0.0);
    current_cap_centres.resize(n_unique_patches);
    current_centre_to_left_patch_boundary_distance.resize(n_unique_patches , 0.0);
    
    
    init_cap_identifier();
    parallelisation_lvl = 0;
    centre_left_cap.resize(3,0.0);
    centre_right_cap.resize(3, 0.0);
    mc_volume_SA.resize(2,0.0);
    Rg_max = input_vars->Rg_max;
    Rb_max = input_vars->Rb_max;
    d_Rg = input_vars->d_Rg;
    d_Rb = input_vars->d_Rb;
    theta_good = input_vars->theta_good;
    theta_bad = input_vars->theta_bad;
    delta = input_vars->delta;
    Rho = input_vars->Rho;
    MAX_MC_POINTS = input_vars->n_max_points;
}


void EvolveSphericalCap::init_cap_identifier()
{
    int n_unique_patches = stripes.get_n_unique_patches();
    current_growing_cap_identifier.resize(n_unique_patches);
    for(int i=0; i<n_unique_patches; i++)
    {
        current_growing_cap_identifier[i] = 0;
    }
    growing_cap_index = 0;
}


void EvolveSphericalCap::calc_vol_SA_normal()
{
    check_boundary->ManageBoxBreach(Cluster_shape_ptr);
    mc_volume_SA = mc_engine->calc_volume_SA(Cluster_shape_ptr);
}

void EvolveSphericalCap::calc_vol_SA_virtual_points()
{
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

void EvolveSphericalCap::calc_volume_SA_central_good_cap(Spherical_cap& GoodCap)
{
    double Volume, SA;
    double N;
    int n_points = mc_engine->get_num_points ();
    if(n_points < MAX_MC_POINTS)
    {
        calc_vol_SA_normal();
    }
    else
    {
        calc_vol_SA_virtual_points();
    }
    
    Volume = mc_volume_SA[0];
    SA = mc_volume_SA[1];
    
    std::vector<double> local_projected_SA(num_patches, 0.0);
    local_projected_SA[0] = GoodCap.projected_SA();
    
    if(parallel_process->is_current_lvl_root(parallelisation_lvl))
    {
        
        output_variables->Volume_global_array.push_back(Volume);
        output_variables->SA_global_array.push_back(SA);
        
        addQuant(output_variables->projected_SA_global_array, local_projected_SA, num_patches);
        N = (Volume * 1e-30) * Rho * avogadro ;
        output_variables->Number_particles_global_array.push_back(N);
        //printf("Volume = %10.10f\t SA=%10.10f\t proj_SA = %10.10f\t N=%10.10f\n", Volume, SA, local_projected_SA[0], N);
        add_radii_and_centre_modifier();
    }
}

void EvolveSphericalCap::calc_volume_SA(Composite_cluster& Comp_cluster)
{
    double Volume, SA;
    double N;
    
    int n_points = mc_engine->get_num_points ();
    if(n_points < MAX_MC_POINTS)
    {
        calc_vol_SA_normal();
    }
    else
    {
        calc_vol_SA_virtual_points();
    }
    
    Volume = mc_volume_SA[0];
    SA = mc_volume_SA[1];
    
    
    if(parallel_process->is_current_lvl_root(parallelisation_lvl))
    {
//       printf("Volume = %10.10f\t SA=%10.10f\n", Volume, SA);
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
        
      
        output_variables->Volume_global_array.push_back(Volume);
        output_variables->SA_global_array.push_back(SA);
        
        addQuant(output_variables->projected_SA_global_array, local_projected_SA, num_patches);
        N = (Volume * 1e-30) * Rho * avogadro ;
        output_variables->Number_particles_global_array.push_back(N);
        add_radii_and_centre_modifier();
        //print_output_variables();

    }
}
void EvolveSphericalCap::add_radii_and_centre_modifier()
{
    for(int i=0; i<current_radii_array.size(); i++)
    {
        output_variables->radii_global_array.push_back(current_radii_array[i]);
        output_variables->clstr_centre_location_modifier_global_array.push_back(current_clstr_centre_modifier_array[i]);
    }
    
}

void EvolveSphericalCap::print_output_variables()
{
    if(parallel_process->is_current_lvl_root(parallelisation_lvl))
    {
        printf("myRank=%d\t", myRank);
        printf("output_variables radii = ");
        for(int i=0; i<output_variables->radii_global_array.size(); i++)
        {
            printf("%10.5f\t", output_variables->radii_global_array[i]);
        }
        printf("\n");
    }
}

void EvolveSphericalCap::add_to_grid_points(int n_points_in_box)
{
    vector<double> this_grid_point;
    for(int i=0; i<current_radii_array.size(); i++)
    {
        this_grid_point.push_back(current_radii_array[i]);
    }
    for(int i=0; i<current_radii_array.size(); i++)
    {
        this_grid_point.push_back((double)current_clstr_centre_modifier_array[i]);
    }
    grid_points.push_back(this_grid_point);
    cost_per_grid_point.push_back(n_points_in_box);
}

void EvolveSphericalCap::clear_grid_points()
{
    grid_points.clear();
    cost_per_grid_point.clear();
}



int EvolveSphericalCap::dummy_evolve(int n_loops)
{
    dummy_set_next_evolve();
    for(int i=0; i<n_loops; i++)
    {
        printf("loop = %d\n", i);
//        setup_next_evolve();
//        identify_growing_cap();
//        identify_cap_type();
//
//        print_initial_variables();
//        print_growing_cap_identifier();
    }
    if(dummy_decider != -1) {dummy_evolve(n_loops);}
    
//     printf("After resetting \n");
//    for(int i=0; i<n_loops; i++)
//    {
//        printf("loop = %d\n", i);
//        reset_evolve();
//        identify_growing_cap();
//        identify_cap_type();
//        print_initial_variables();
//        print_growing_cap_identifier();
//    }
    return 1;
}

void EvolveSphericalCap::dummy_set_next_evolve()
{
    dummy_decider = (dummy_decider == 0)? (dummy_decider+1): -1;
}

int EvolveSphericalCap::evolve ()
{
    /*
     1. check the current evolving cap.
     2. see if it is bad or good.
     3. check box breach boundary
     3. if(box not breached):
        1. start the for loop
        2. if (check patch boundary crossing)
            1. update growing cluster identifier
            2. call evolve
            3. After evolve returns reset the growing cluster identifier.
        3. return EXIT at the end of the for loop.
     4. else: return EXIT
     */
    int check_radius_min_within_radius_max;
    int evolve_output;
    if(parallel_process->is_parent()){printf("Inside evolve\n");}
    setup_next_evolve();
    check_radius_min_within_radius_max = identify_current_cap_growth_conditions();
    if(check_radius_min_within_radius_max == 0)
    {
        reset_evolve();
        return 0;
    }
    
    double z_surface = 0.0;
    double Radius, projected_radius ;

    int size_loop = (int)Radius_loop_start_end.size();
    for(int j = 0 ; j < size_loop ; j++)
    {
        
        int i= Radius_loop_start_end[j];
        Radius = current_radius_min + i * d_current_radius ;  //  startingRg + i*d_Rg            //sphere's radius
        projected_radius = Radius * sin(current_theta);
        
        if(parallel_process->is_current_lvl_root(parallelisation_lvl)){printf("parallelisation_lvl = %d i=%d Radius=%10.3f\n", parallelisation_lvl, i, Radius);}
        
        clstr_centre_location_modifier = (growing_cap_index ==0) ? 0 : clstr_centre_location_modifier_array_load_balancing[j];
        current_projected_radii_array[growing_cap_index] = projected_radius;
        current_radii_array[growing_cap_index] = Radius;
        current_clstr_centre_modifier_array[growing_cap_index] = clstr_centre_location_modifier;
        
        int should_write_to_file;
        should_write_to_file = (j == size_loop-1 || j == size_loop/2) ? 1 : 0;
//        if(size_loop > 100)
//        {
//            should_write_to_file = (j%20 == 0  ||  j == size_loop-1) ? 1 : 0;
//        }
//        else
//        {
//            should_write_to_file = (j%(size_loop/4) == 0  ||  j == size_loop-1) ? 1 : 0;
//        }
 
        if(growing_cap_index == 0)
        {
            std::vector<double> centre_good(3,0.0);
            current_cap_centres[growing_cap_index] = centre_good;
    
            Spherical_cap GoodCap (centre_good, projected_radius, current_theta, z_surface);
            
            Cluster_shape_ptr= &GoodCap;

            if(check_breaking_condition())
            {
                printf("check breaking condition true for central cap\n");
                break;
            }
            else
            {
                if (patch_bounds_crossed())
                {
                    if(growing_capsList.empty())
                    {
                        growing_capsList.push_back(GoodCap);
                    }
                    else
                    {
                        growing_capsList[0] = GoodCap;
                    }
                    calc_and_set_centre_to_left_patch_boundary_distance();
                    write_to_file_periodically();
                    evolve_output = evolve();
                    if(evolve_output == 0){break;}
                }
                else
                {
                    calc_volume_SA_central_good_cap(GoodCap);
                    if(should_write_to_file==1)
                    {
                        write_to_file_periodically();
                    }
                }
            }
            
        }
        else
        {
            double inside_sq = projected_radius*projected_radius - previous_cap_projected_radius*previous_cap_projected_radius + previous_centre_to_left_patch_boundary_distance * previous_centre_to_left_patch_boundary_distance ;
            
            if(abs(inside_sq)<1e-10 && inside_sq < 0.0) {inside_sq = 0.0;}
            else if(inside_sq<0.0)
            {
                printf("Inside sq for centre calculation cannot be zero. Ether radius being used is less than minimum or wrong patch boundary crossing is occuring\n");
                abort();
                
            }
            current_centre_to_previous_patch_boundary_distance = clstr_centre_location_modifier * sqrt(inside_sq) ;
            
            /* Symmetric caps on either side of central patch */
            centre_left_cap[0] = -(previous_patch_symmetric_boundary) - current_centre_to_previous_patch_boundary_distance;
            centre_right_cap[0] = (previous_patch_symmetric_boundary) + current_centre_to_previous_patch_boundary_distance;
            
            current_cap_centres[growing_cap_index] = centre_left_cap;
            
            Spherical_cap Cap_left(centre_left_cap, projected_radius, current_theta, z_surface);
            Spherical_cap Cap_right(centre_right_cap, projected_radius, current_theta, z_surface);
            
            std::vector<Spherical_cap> capsList = growing_capsList;
            capsList.push_back(Cap_left);
            capsList.push_back(Cap_right);
            
            Composite_cluster Comp_cluster (capsList, stripes);
            Cluster_shape_ptr = &Comp_cluster;
            if(check_breaking_condition())
            {
                printf("check breaking condition true for growing cap=%d\n", growing_cap_index);
                break;
            }
            else
            {
                if (patch_bounds_crossed())
                {
                    //updating the growing_capsList to contain the current boundary crossed cap. push_back will work here since the growing_capsList vector is resized after every reset of evolve.
                    growing_capsList.push_back(Cap_left);
                    growing_capsList.push_back(Cap_right);
                    calc_and_set_centre_to_left_patch_boundary_distance();
                    write_to_file_periodically();
                    evolve_output = evolve();
                    if(evolve_output == 0){break;}
                }
                else
                {
                    calc_volume_SA(Comp_cluster);
                    if(should_write_to_file==1)
                    {
                        write_to_file_periodically();
                        
                    }
                }
            }
        }
        
    }//for loop
    
    reset_evolve();
    return 1;
}


int EvolveSphericalCap::evolve_profiling()
{
    int check_radius_min_within_radius_max;
    int evolve_output;
    if(parallel_process->is_parent()){printf("Inside evolve\n");}
    setup_next_evolve();
    check_radius_min_within_radius_max = identify_current_cap_growth_conditions();
    if(check_radius_min_within_radius_max == 0)
    {
        reset_evolve();
        return 0;
    }
    int n_points_in_box;
    double z_surface = 0.0;
    double Radius, projected_radius ;

    int size_loop = (int)Radius_loop_start_end.size();
    for(int j = 0 ; j < size_loop; j++) // len_Rg
    {
        
        int i = Radius_loop_start_end[j];
        
        Radius = current_radius_min + i * d_current_radius ;  //  startingRg + i*d_Rg            //sphere's radius
        projected_radius = Radius * sin(current_theta);
        
        if(parallel_process->is_current_lvl_root(parallelisation_lvl)){printf("parallelisation_lvl = %d i=%d Radius=%10.3f\n", parallelisation_lvl, i, Radius);}
        
        clstr_centre_location_modifier = (growing_cap_index ==0) ? 0 : clstr_centre_location_modifier_array_load_balancing[j];
        current_projected_radii_array[growing_cap_index] = projected_radius;
        current_radii_array[growing_cap_index] = Radius;
        current_clstr_centre_modifier_array[growing_cap_index] = clstr_centre_location_modifier;
        int should_write_to_file;
        should_write_to_file = (j == size_loop-1 || j == size_loop/2) ? 1 : 0;
//        if(size_loop > 100)
//        {
//            should_write_to_file = (j%20 == 0  ||  j == size_loop-1) ? 1 : 0;
//        }
//        else
//        {
//            should_write_to_file = (j%(size_loop/4) == 0  ||  j == size_loop-1) ? 1 : 0;
//        }
        
        if(growing_cap_index == 0)
        {
            std::vector<double> centre_good(3,0.0);
            current_cap_centres[growing_cap_index] = centre_good;
            
            Spherical_cap GoodCap (centre_good, projected_radius, current_theta, z_surface);
            
            Cluster_shape_ptr= &GoodCap;
            
            check_boundary->ManageBoxBreach(Cluster_shape_ptr);
            n_points_in_box = mc_engine->get_num_points();
            
            if(check_breaking_condition())
            {
                printf("check breaking condition true for central cap\n");
                break;
            }
            else
            {
                if (patch_bounds_crossed())
                {
                    //                    printf("Patch bounds crossed for good patch at i=%d\n", i);
                    if(growing_capsList.empty())
                    {
                        growing_capsList.push_back(GoodCap);
                    }
                    else
                    {
                        growing_capsList[0] = GoodCap;
                    }
                    calc_and_set_centre_to_left_patch_boundary_distance();
                    write_profiling_periodically(); //Writing before every level up i.e. before evolve.
                    evolve_output = evolve_profiling();
                    
                    if(evolve_output == 0){break;}
                }
                else
                {
                    add_to_grid_points(n_points_in_box);
                    if(should_write_to_file==1) //And writing at the end of every loop
                    {
                        write_profiling_periodically();
                    }
                }
            }
            
        }
        else
        {
            double inside_sq = projected_radius*projected_radius - previous_cap_projected_radius*previous_cap_projected_radius + previous_centre_to_left_patch_boundary_distance * previous_centre_to_left_patch_boundary_distance ;
            
            if(abs(inside_sq)<1e-10 && inside_sq < 0.0) {inside_sq = 0.0;}
            else if(inside_sq<0.0)
            {
                printf("Inside sq for centre calculation cannot be zero. Ether radius being used is less than minimum or wrong patch boundary crossing is occuring\n");
                abort();
                
            }
            current_centre_to_previous_patch_boundary_distance = clstr_centre_location_modifier * sqrt(inside_sq) ;
            
            /* Symmetric caps on either side of central patch */
            centre_left_cap[0] = -(previous_patch_symmetric_boundary) - current_centre_to_previous_patch_boundary_distance;
            centre_right_cap[0] = (previous_patch_symmetric_boundary) + current_centre_to_previous_patch_boundary_distance;
            
            current_cap_centres[growing_cap_index] = centre_left_cap;
            
            Spherical_cap Cap_left(centre_left_cap, projected_radius, current_theta, z_surface);
            Spherical_cap Cap_right(centre_right_cap, projected_radius, current_theta, z_surface);
            
            std::vector<Spherical_cap> capsList = growing_capsList;
            capsList.push_back(Cap_left);
            capsList.push_back(Cap_right);
            
            Composite_cluster Comp_cluster (capsList, stripes);
            Cluster_shape_ptr = &Comp_cluster;
            
            check_boundary->ManageBoxBreach(Cluster_shape_ptr);
            n_points_in_box = mc_engine->get_num_points();
            
            if(check_breaking_condition())
            {
                printf("check breaking condition true for growing cap=%d\n", growing_cap_index);
                break;
            }
            else
            {
                if (patch_bounds_crossed())
                {
                    //updating the growing_capsList to contain the current boundary crossed cap. push_back will work here since the growing_capsList vector is resized after every reset of evolve.
                    growing_capsList.push_back(Cap_left);
                    growing_capsList.push_back(Cap_right);
                    calc_and_set_centre_to_left_patch_boundary_distance();
                    write_profiling_periodically();
                    evolve_output = evolve_profiling();
                    if(evolve_output == 0){break;}
                }
                else
                {
                    add_to_grid_points(n_points_in_box);
                    if(should_write_to_file==1) //And writing at the end of every loop
                    {
                        write_profiling_periodically();
                    }
                }
            }
        }
        
    }//for loop
    
    reset_evolve();
    return 1;
}

int EvolveSphericalCap::artificial_evolve(std::vector<double> input_Radii)
{
    int check_radius_min_within_radius_max;
    int evolve_output;
    if(parallel_process->is_parent()){printf("Inside evolve\n");}
    setup_next_evolve_constant_modifier();
    check_radius_min_within_radius_max = artificial_growth_conditions(input_Radii);
    if(check_radius_min_within_radius_max == 0)
    {
        artificial_reset_evolve(input_Radii);
        return 0;
    }
    
    double z_surface = 0.0;
    double Radius, projected_radius ;
    
    int size_loop = (int)Radius_loop_start_end.size();
    for(int j = 0 ; j < size_loop ; j++)
    {
        
        int i= Radius_loop_start_end[j];
        if(i!=0){printf("more than 1 Rs for artificial, growing_cap_index=%d\n", growing_cap_index); abort();}
        Radius = current_radius_min + i * d_current_radius ;  //  startingRg + i*d_Rg            //sphere's radius
        projected_radius = Radius * sin(current_theta);
        
        printf("growing_cap_index=%d\t current_radius_min=%10.10f\t Radius=%10.5f\t projected_radius=%10.10f\t current_theta=%10.10f\t sin(current_theta)=%10.10f\n", growing_cap_index, current_radius_min, Radius, projected_radius, current_theta, sin(current_theta));
        
        //if(parallel_process->is_current_lvl_root(parallelisation_lvl))
        if(myRank==0){printf("parallelisation_lvl = %d i=%d Radius=%10.3f\n", parallelisation_lvl, i, Radius);}
        
        
        current_projected_radii_array[growing_cap_index] = projected_radius;
        current_radii_array[growing_cap_index] = Radius;
        current_clstr_centre_modifier_array[growing_cap_index] = clstr_centre_location_modifier;
        
        int should_write_to_file;
        should_write_to_file = (j == size_loop-1 || j == size_loop/2) ? 1 : 0;
        
        if(growing_cap_index == 0)
        {
            std::vector<double> centre_good(3,0.0);
            current_cap_centres[growing_cap_index] = centre_good;
            
            Spherical_cap GoodCap (centre_good, projected_radius, current_theta, z_surface);
            
            Cluster_shape_ptr= &GoodCap;
            if(check_breaking_condition())
            {
                printf("check breaking condition true for central cap\n");
                break;
            }
            else
            {
                if (patch_bounds_crossed())
                {
                    if(growing_capsList.empty())
                    {
                        growing_capsList.push_back(GoodCap);
                    }
                    else
                    {
                        growing_capsList[0] = GoodCap;
                    }
                    calc_and_set_centre_to_left_patch_boundary_distance();
                    write_to_file_periodically();
                    evolve_output = artificial_evolve(input_Radii);
                    if(evolve_output == 0){break;}
                }
                else
                {
                    calc_volume_SA_central_good_cap(GoodCap);
                    if(should_write_to_file == 1)
                    {
                        write_to_file_periodically();
                    }
                    
                }
            }
            
        }
        else
        {
            double inside_sq = projected_radius*projected_radius - previous_cap_projected_radius*previous_cap_projected_radius + previous_centre_to_left_patch_boundary_distance * previous_centre_to_left_patch_boundary_distance ;
            
            if(abs(inside_sq)<1e-10 && inside_sq < 0.0) {inside_sq = 0.0;}
            else if(inside_sq<0.0)
            {
                printf("Inside sq for centre calculation cannot be zero. Ether radius being used is less than minimum or wrong patch boundary crossing is occuring\n");
                abort();
                
            }
            current_centre_to_previous_patch_boundary_distance = clstr_centre_location_modifier * sqrt(inside_sq) ;
            
            /* Symmetric caps on either side of central patch */
            centre_left_cap[0] = -(previous_patch_symmetric_boundary) - current_centre_to_previous_patch_boundary_distance;
            centre_right_cap[0] = (previous_patch_symmetric_boundary) + current_centre_to_previous_patch_boundary_distance;
            
            current_cap_centres[growing_cap_index] = centre_left_cap;
            
            Spherical_cap Cap_left(centre_left_cap, projected_radius, current_theta, z_surface);
            Spherical_cap Cap_right(centre_right_cap, projected_radius, current_theta, z_surface);
            
            std::vector<Spherical_cap> capsList = growing_capsList;
            capsList.push_back(Cap_left);
            capsList.push_back(Cap_right);
            
            Composite_cluster Comp_cluster (capsList, stripes);
            Cluster_shape_ptr = &Comp_cluster;
            
            if(check_breaking_condition())
            {
                printf("check breaking condition true for growing cap=%d\n", growing_cap_index);
                break;
            }
            else
            {
                if (patch_bounds_crossed())
                {
                    //updating the growing_capsList to contain the current boundary crossed cap. push_back will work here since the growing_capsList vector is resized after every reset of evolve.
                    growing_capsList.push_back(Cap_left);
                    growing_capsList.push_back(Cap_right);
                    calc_and_set_centre_to_left_patch_boundary_distance();
                    write_to_file_periodically();
                    evolve_output = artificial_evolve(input_Radii);
                    if(evolve_output == 0){break;}
                }
                else
                {
                    calc_volume_SA(Comp_cluster);
                    if(should_write_to_file == 1)
                    {
                        write_to_file_periodically();
                    }
                }
            }
        }
        
    }//for loop
    
    artificial_reset_evolve(input_Radii);
    return 1;
}

int EvolveSphericalCap::evolve_constant_modifiers ()
{
    /*
     1. check the current evolving cap.
     2. see if it is bad or good.
     3. check box breach boundary
     3. if(box not breached):
     1. start the for loop
     2. if (check patch boundary crossing)
     1. update growing cluster identifier
     2. call evolve
     3. After evolve returns reset the growing cluster identifier.
     3. return EXIT at the end of the for loop.
     4. else: return EXIT
     */
    int check_radius_min_within_radius_max;
    int evolve_output;
    if(parallel_process->is_parent()){printf("Inside evolve\n");}
    setup_next_evolve_constant_modifier();
    check_radius_min_within_radius_max = identify_current_cap_growth_conditions_constant_modifiers();
    if(check_radius_min_within_radius_max == 0)
    {
        reset_evolve_constant_modifier();
        return 0;
    }
    
    double z_surface = 0.0;
    double Radius, projected_radius ;
    
    int size_loop = (int)Radius_loop_start_end.size();
    for(int j = 0 ; j < size_loop ; j++)
    {
        
        int i= Radius_loop_start_end[j];
        Radius = current_radius_min + i * d_current_radius ;  //  startingRg + i*d_Rg            //sphere's radius
        projected_radius = Radius * sin(current_theta);
        
        //if(parallel_process->is_current_lvl_root(parallelisation_lvl))
        if(myRank==0){printf("parallelisation_lvl = %d i=%d Radius=%10.3f\n", parallelisation_lvl, i, Radius);}
        
        
        current_projected_radii_array[growing_cap_index] = projected_radius;
        current_radii_array[growing_cap_index] = Radius;
        current_clstr_centre_modifier_array[growing_cap_index] = clstr_centre_location_modifier;
        
        int should_write_to_file;
        should_write_to_file = (j == size_loop-1 || j == size_loop/2) ? 1 : 0;
//        if(size_loop > 100)
//        {
//            should_write_to_file = (j%20 == 0  ||  j == size_loop-1) ? 1 : 0;
//        }
//        else
//        {
//            should_write_to_file = (j%(size_loop/4) == 0  ||  j == size_loop-1) ? 1 : 0;
//        }
        
        if(growing_cap_index == 0)
        {
            std::vector<double> centre_good(3,0.0);
            current_cap_centres[growing_cap_index] = centre_good;
            
            Spherical_cap GoodCap (centre_good, projected_radius, current_theta, z_surface);
            
            Cluster_shape_ptr= &GoodCap;
            if(check_breaking_condition())
            {
                printf("check breaking condition true for central cap\n");
                break;
            }
            else
            {
                if (patch_bounds_crossed())
                {
                    if(growing_capsList.empty())
                    {
                        growing_capsList.push_back(GoodCap);
                    }
                    else
                    {
                        growing_capsList[0] = GoodCap;
                    }
                    calc_and_set_centre_to_left_patch_boundary_distance();
                    write_to_file_periodically();
                    evolve_output = evolve_constant_modifiers();
                    if(evolve_output == 0){break;}
                }
                else
                {
                    calc_volume_SA_central_good_cap(GoodCap);
                    if(should_write_to_file == 1)
                    {
                        write_to_file_periodically();
                    }
                        
                }
            }
            
        }
        else
        {
            double inside_sq = projected_radius*projected_radius - previous_cap_projected_radius*previous_cap_projected_radius + previous_centre_to_left_patch_boundary_distance * previous_centre_to_left_patch_boundary_distance ;
            
            if(abs(inside_sq)<1e-10 && inside_sq < 0.0) {inside_sq = 0.0;}
            else if(inside_sq<0.0)
            {
                printf("Inside sq for centre calculation cannot be zero. Ether radius being used is less than minimum or wrong patch boundary crossing is occuring\n");
                abort();
                
            }
            current_centre_to_previous_patch_boundary_distance = clstr_centre_location_modifier * sqrt(inside_sq) ;
            
            /* Symmetric caps on either side of central patch */
            centre_left_cap[0] = -(previous_patch_symmetric_boundary) - current_centre_to_previous_patch_boundary_distance;
            centre_right_cap[0] = (previous_patch_symmetric_boundary) + current_centre_to_previous_patch_boundary_distance;
            
            current_cap_centres[growing_cap_index] = centre_left_cap;
            
            Spherical_cap Cap_left(centre_left_cap, projected_radius, current_theta, z_surface);
            Spherical_cap Cap_right(centre_right_cap, projected_radius, current_theta, z_surface);
            
            std::vector<Spherical_cap> capsList = growing_capsList;
            capsList.push_back(Cap_left);
            capsList.push_back(Cap_right);
            
            Composite_cluster Comp_cluster (capsList, stripes);
            Cluster_shape_ptr = &Comp_cluster;
            
            if(check_breaking_condition())
            {
                printf("check breaking condition true for growing cap=%d\n", growing_cap_index);
                break;
            }
            else
            {
                if (patch_bounds_crossed())
                {
                    //updating the growing_capsList to contain the current boundary crossed cap. push_back will work here since the growing_capsList vector is resized after every reset of evolve.
                    growing_capsList.push_back(Cap_left);
                    growing_capsList.push_back(Cap_right);
                    calc_and_set_centre_to_left_patch_boundary_distance();
                    write_to_file_periodically();
                    evolve_output = evolve_constant_modifiers();
                    if(evolve_output == 0){break;}
                }
                else
                {
                    calc_volume_SA(Comp_cluster);
                    if(should_write_to_file == 1)
                    {
                        write_to_file_periodically();
                    }
                }
            }
        }
        
    }//for loop
    
    reset_evolve_constant_modifier();
    return 1;
}

int EvolveSphericalCap::evolve_profiling_constant_modifiers()
{
    int check_radius_min_within_radius_max;
    int evolve_output;
    setup_next_evolve_constant_modifier();
    check_radius_min_within_radius_max = identify_current_cap_growth_conditions_constant_modifiers();
    if(check_radius_min_within_radius_max == 0)
    {
        reset_evolve_constant_modifier();
        return 0;
    }
    int n_points_in_box;
    double z_surface = 0.0;
    double Radius, projected_radius ;
    
    
    int size_loop = (int)Radius_loop_start_end.size();

    for(int j = 0 ; j < size_loop; j++) // len_Rg
    {
        
        int i = Radius_loop_start_end[j];
        
        Radius = current_radius_min + i * d_current_radius ;  //  startingRg + i*d_Rg            //sphere's radius
        projected_radius = Radius * sin(current_theta);
        
        if(myRank==12){printf("parallelisation_lvl = %d i=%d Radius=%10.3f\t rank=%d\n", parallelisation_lvl, i, Radius, myRank);}
            
        
        
        current_projected_radii_array[growing_cap_index] = projected_radius;
        current_radii_array[growing_cap_index] = Radius;
        current_clstr_centre_modifier_array[growing_cap_index] = clstr_centre_location_modifier;
        int should_write_to_file;
        should_write_to_file = (j == size_loop-1 || j == size_loop/2) ? 1 : 0; //
        
        if(growing_cap_index == 0)
        {
            std::vector<double> centre_good(3,0.0);
            current_cap_centres[growing_cap_index] = centre_good;
            
            Spherical_cap GoodCap (centre_good, projected_radius, current_theta, z_surface);
            
            Cluster_shape_ptr= &GoodCap;
            
            check_boundary->ManageBoxBreach(Cluster_shape_ptr);
            n_points_in_box = mc_engine->get_num_points();
            
            if(check_breaking_condition())
            {
                printf("check breaking condition true for central cap\n");
                break;
            }
            else
            {
                if (patch_bounds_crossed())
                {
                    //                    printf("Patch bounds crossed for good patch at i=%d\n", i);
                    if(growing_capsList.empty())
                    {
                        growing_capsList.push_back(GoodCap);
                    }
                    else
                    {
                        growing_capsList[0] = GoodCap;
                    }
                    calc_and_set_centre_to_left_patch_boundary_distance();
                    //write_profiling_periodically(); //Writing before every level up i.e. before evolve.
                    evolve_output = evolve_profiling_constant_modifiers();
                    
                    if(evolve_output == 0){break;}
                }
                else
                {
                    add_to_grid_points(n_points_in_box);
                    if(should_write_to_file==1) //And writing at the end of every loop
                    {
                        write_profiling_periodically();
                    }
                }
            }
            
        }
        else
        {
            double inside_sq = projected_radius*projected_radius - previous_cap_projected_radius*previous_cap_projected_radius + previous_centre_to_left_patch_boundary_distance * previous_centre_to_left_patch_boundary_distance ;
            
            if(abs(inside_sq)<1e-10 && inside_sq < 0.0) {inside_sq = 0.0;}
            else if(inside_sq<0.0)
            {
                printf("Inside sq for centre calculation cannot be zero. Ether radius being used is less than minimum or wrong patch boundary crossing is occuring\n");
                abort();
                
            }
            current_centre_to_previous_patch_boundary_distance = clstr_centre_location_modifier * sqrt(inside_sq) ;
            /* Symmetric caps on either side of central patch */
            centre_left_cap[0] = -(previous_patch_symmetric_boundary) - current_centre_to_previous_patch_boundary_distance;
            centre_right_cap[0] = (previous_patch_symmetric_boundary) + current_centre_to_previous_patch_boundary_distance;
            
            current_cap_centres[growing_cap_index] = centre_left_cap;
            
            Spherical_cap Cap_left(centre_left_cap, projected_radius, current_theta, z_surface);
            Spherical_cap Cap_right(centre_right_cap, projected_radius, current_theta, z_surface);
            
            std::vector<Spherical_cap> capsList = growing_capsList;
            capsList.push_back(Cap_left);
            capsList.push_back(Cap_right);
            
            Composite_cluster Comp_cluster (capsList, stripes);
            Cluster_shape_ptr = &Comp_cluster;
            
            check_boundary->ManageBoxBreach(Cluster_shape_ptr);
            n_points_in_box = mc_engine->get_num_points();
            
            if(check_breaking_condition())
            {
                
                if(myRank==12)
                {
                    printf("check breaking condition true for growing cap=%d\n", growing_cap_index);
                    printf("parallelisation level=%d\t Radius=%10.3f\t size_loop=%d\n", parallelisation_lvl, Radius, size_loop);
                }
                reset_evolve_constant_modifier();
                if(myRank==12)
                {printf("After reset evolve before returning 1\n");}
                
                return 1;
                //break;
            }
            else
            {
                if (patch_bounds_crossed())
                {
                    //updating the growing_capsList to contain the current boundary crossed cap. push_back will work here since the growing_capsList vector is resized after every reset of evolve.
                    growing_capsList.push_back(Cap_left);
                    growing_capsList.push_back(Cap_right);
                    calc_and_set_centre_to_left_patch_boundary_distance();
                    //write_profiling_periodically();
                    if(myRank == 12)
                    {
                        printf("BEFORE next evolve parallelisation level=%d\t growing_cap=%d\tsize_loop=%d\n", parallelisation_lvl, growing_cap_index, size_loop);
                    }
                    evolve_output = evolve_profiling_constant_modifiers();
                    if(myRank == 12)
                    {
                        printf("Out of evolve parallelisation level=%d\t growing_cap=%d\tsize_loop=%d\n", parallelisation_lvl, growing_cap_index, size_loop);
                    }
                    if(evolve_output == 0){break;}
                }
                else
                {
                    add_to_grid_points(n_points_in_box);
                    if(should_write_to_file==1) //And writing at the end of every loop
                    {
                        write_profiling_periodically();
                    }
                }
            }
        }
        
    }//for loop
    
    reset_evolve_constant_modifier();
    return 1;
}



bool EvolveSphericalCap::has_no_caps()
{
    std::vector<int>::iterator it_one = std::find(current_growing_cap_identifier.begin(), current_growing_cap_identifier.end(), 1);
    return (it_one == current_growing_cap_identifier.end());
}

bool EvolveSphericalCap::has_one_cap()
{
    std::vector<int>::iterator it_one = std::find(current_growing_cap_identifier.begin(), current_growing_cap_identifier.end(), 1);
    return (it_one == current_growing_cap_identifier.begin());
}

bool EvolveSphericalCap::has_more_than_one_caps()
{
    std::vector<int> reversed_current_growing_cap_identifier = current_growing_cap_identifier;
    std::reverse(reversed_current_growing_cap_identifier.begin(),reversed_current_growing_cap_identifier.end());
    
    std::vector<int>::iterator it_one = std::find(reversed_current_growing_cap_identifier.begin(), reversed_current_growing_cap_identifier.end(), 1);
    
    return (it_one != (reversed_current_growing_cap_identifier.end()-1));
}

void EvolveSphericalCap::setup_next_evolve()
{
    if(has_no_caps())
    {
         parallelisation_lvl = 0;
    }
    else
    {parallelisation_lvl += 2;} //a jump of two because for each new set of cap evolutions on a new patch there is a parallelisation level before it that corresponds to the two clstr centre location modifiers {+1,-1} . Increase parallelisation level by 2 because the intermediary level will be for clstr_centre identifier loop
    update_growing_cap_identifier();
    set_MC_comm();
}

void EvolveSphericalCap::reset_evolve()
{
    if(has_more_than_one_caps())
    {
       parallelisation_lvl -= 2;
    }
    else
    {
      parallelisation_lvl = 0;
    }
    reset_current_growth_conditions_array();
    reset_growing_cap_identifier();
    
    reset_growing_capsList();
    
    identify_current_cap_growth_conditions(); //This is important to keep running the for loop of the previous cap.
    set_MC_comm();
    
}

void EvolveSphericalCap::setup_next_evolve_constant_modifier()
{
    if(has_no_caps())
    {
        parallelisation_lvl = 0;
    }
    else
    {parallelisation_lvl += 1;} //a jump of one because no longer an extra level for two clstr centre modifiers since we are using constant modifier (1 for good patches, -1 for bad patches, except the central good patch).
    update_growing_cap_identifier();
    set_MC_comm();
}

void EvolveSphericalCap::reset_evolve_constant_modifier()
{
//    if(myRank==12)
//    {
//    printf("Inside reset_evolve_constant_modifier myRank=%d\n", myRank);
//        printf("parallelisation level=%d\t growing_cap_index=%d\t Radius_loop_start_end size=%d\t myRank=%d\n", parallelisation_lvl, growing_cap_index, (int)Radius_loop_start_end.size(), myRank);
//
//    }
    if(has_more_than_one_caps())
    {
        parallelisation_lvl -= 1;
    }
    else
    {
        parallelisation_lvl = 0;
    }
    reset_current_growth_conditions_array();
    reset_growing_cap_identifier();
    
    reset_growing_capsList();
    
    identify_current_cap_growth_conditions_constant_modifiers(); //This is important to keep running the for loop of the previous cap.
    set_MC_comm();
//    if(myRank==12)
//    {
//        printf("End of reset evolve constant\n");
//        printf("parallelisation level=%d\t growing_cap_index=%d\t Radius_loop_start_end size=%d\t  myRank=%d\n", parallelisation_lvl, growing_cap_index, (int)Radius_loop_start_end.size(), myRank);
//
//    }
    
    
    // reset growing_capsList by 2 spaces
    // reset parallelization level by 2
    // update current_growing_cap_identifier
}

void EvolveSphericalCap::artificial_reset_evolve(std::vector<double> input_Radii)
{
    if(has_more_than_one_caps())
    {
        parallelisation_lvl -= 1;
    }
    else
    {
        parallelisation_lvl = 0;
    }
    reset_current_growth_conditions_array();
    reset_growing_cap_identifier();
    
    reset_growing_capsList();
    
    artificial_growth_conditions(input_Radii); //This is important to keep running the for loop of the previous cap.
    set_MC_comm();
}

void  EvolveSphericalCap::update_growing_cap_identifier()
{
    std::vector<int>::iterator it_zero = std::find(current_growing_cap_identifier.begin(), current_growing_cap_identifier.end(), 0);
    (*it_zero) = 1;
}

void  EvolveSphericalCap::reset_growing_cap_identifier()
{
    if(!has_no_caps())
    {
        std::vector<int> reversed_current_growing_cap_identifier = current_growing_cap_identifier;
        std::reverse(reversed_current_growing_cap_identifier.begin(),reversed_current_growing_cap_identifier.end());
        std::vector<int>::iterator it_one = std::find(reversed_current_growing_cap_identifier.begin(), reversed_current_growing_cap_identifier.end(), 1);
        (*it_one) = 0;
    std::reverse(reversed_current_growing_cap_identifier.begin(),reversed_current_growing_cap_identifier.end());
        current_growing_cap_identifier = reversed_current_growing_cap_identifier;
    }
}

void EvolveSphericalCap::reset_growing_capsList()
{
    int size = (int)growing_capsList.size();
    if(!growing_capsList.empty())
    {
        if(size == 1)
        {
            growing_capsList.clear();
        }
        else if (size % 2 == 0)
        {
            printf("growing capsList has an even size, which is not possible\n");
            abort();
        }
        else
        {
            growing_capsList.resize(size-2); //reducing the size by 2.
        }
    }
//    else
//    {
//        printf("growing capsList is empty when trying to reset\n");
//    }
}

void EvolveSphericalCap::reset_current_growth_conditions_array()
{
    current_radii_array[growing_cap_index] = 0.0;
    current_projected_radii_array[growing_cap_index] = 0.0;
    current_clstr_centre_modifier_array[growing_cap_index] = 0.0;
    current_cap_centres[growing_cap_index].clear();
    current_centre_to_left_patch_boundary_distance[growing_cap_index] = 0.0;
}



int EvolveSphericalCap::identify_current_cap_growth_conditions()
{
    identify_growing_cap();
    
    double projected_radius_min;
    if (growing_cap_index == 0)
    {
        current_radius_min = 0.0;
        current_radius_max = Rg_max;
        d_current_radius = d_Rg;
        length_current_radius = (int) ((current_radius_max-current_radius_min)/d_current_radius) + 1;
        if(length_current_radius < 0){length_current_radius = 0.0;}
        Radius_loop_start_end = parallel_process->getLoopStartEnd_balanced_binary_split (length_current_radius, parallel_process->get_worker_color_for_level (parallelisation_lvl));
        current_theta = theta_good;
        clstr_centre_location_modifier = 0;
        
        if(current_radius_max < current_radius_min)
        {
//            printf("myRank=%d parallel_lvl=%d length_current_radius = %d\t current_radius_max = %10.5f current_radius_min=%10.5f d_current_radius=%10.5f\n", myRank, parallelisation_lvl, length_current_radius, current_radius_max, current_radius_min, d_current_radius);
            return 0;
        }
    }
    else
    {
        if (growing_cap_index % 2 == 0) //good patch
        {
            d_current_radius = d_Rg;
            current_theta = theta_good;
            current_radius_max = Rg_max;
        }
        else
        {
            d_current_radius = d_Rb;
            current_theta = theta_bad;
            current_radius_max = Rb_max; //This is untrue when breaking happens!!!!! Needed for load balancing
            
        }
        previous_cap_projected_radius = current_projected_radii_array[growing_cap_index - 1];
        previous_centre_to_left_patch_boundary_distance = current_centre_to_left_patch_boundary_distance[growing_cap_index - 1];
        projected_radius_min = sqrt(previous_cap_projected_radius*previous_cap_projected_radius - (previous_centre_to_left_patch_boundary_distance * previous_centre_to_left_patch_boundary_distance));
        
        int true_previous_left_patch_index  = (growing_cap_index == 1) ? 0 : 2*(growing_cap_index-1) - 1;
        std::vector<double> previous_left_patch_bounds = stripes.Patches[true_previous_left_patch_index].patch_boundaries();
        previous_patch_symmetric_boundary = abs(previous_left_patch_bounds[0]);
        
        current_radius_min = projected_radius_min/(double)sin(current_theta) ;
        length_current_radius = (int) ((current_radius_max-current_radius_min)/d_current_radius) + 1;
        if(length_current_radius < 0){length_current_radius = 0.0;}
        Radius_loop_start_end = parallel_process->getLoopStartEnd_balanced_binary_split (length_current_radius, parallel_process->get_worker_color_for_level (parallelisation_lvl));
        int prev_parallelisation_lvl = (parallelisation_lvl == 0) ? 0 : (parallelisation_lvl - 1);
        
        
        get_balanced_modifier_array((int)Radius_loop_start_end.size(), parallel_process->is_current_lvl_even_branch(prev_parallelisation_lvl));
//        if (parallel_process->is_current_lvl_even_branch(prev_parallelisation_lvl)) //Here parallelisation level corresponds to Radius loop parallelization. cap position parallelisation is one level above that.
//        {
//            clstr_centre_location_modifier = 1;
//        }
//        else
//        {
//            clstr_centre_location_modifier = -1;
//        }
        
        if(current_radius_max < current_radius_min)
        {
            return 0;
        }
    }
    return 1;
}

int EvolveSphericalCap::artificial_growth_conditions(std::vector<double> input_Radii)
{
    identify_growing_cap();
    double projected_radius_min;
    if (growing_cap_index == 0)
    {
        current_radius_min = input_Radii[growing_cap_index];
        current_radius_max = current_radius_min;
        d_current_radius = d_Rg;
        length_current_radius = (int) ((current_radius_max-current_radius_min)/d_current_radius) + 1;
        if(length_current_radius < 0){length_current_radius = 0.0;}
        Radius_loop_start_end = parallel_process->getLoopStartEnd_balanced_binary_split (length_current_radius, parallel_process->get_worker_color_for_level (parallelisation_lvl));
        current_theta = theta_good;
        clstr_centre_location_modifier = 0;
        
        if(current_radius_max < current_radius_min)
        {
            return 0;
        }
    }
    else
    {
        previous_cap_projected_radius = current_projected_radii_array[growing_cap_index - 1];
        previous_centre_to_left_patch_boundary_distance = current_centre_to_left_patch_boundary_distance[growing_cap_index - 1];
        projected_radius_min = sqrt(previous_cap_projected_radius*previous_cap_projected_radius - (previous_centre_to_left_patch_boundary_distance * previous_centre_to_left_patch_boundary_distance));
        
        if (growing_cap_index % 2 == 0) //good patch
        {
            d_current_radius = d_Rg;
            current_theta = theta_good;
            clstr_centre_location_modifier = 1;
            current_radius_max = input_Radii[growing_cap_index];
        }
        else
        {
            d_current_radius = d_Rb;
            current_theta = theta_bad;
            clstr_centre_location_modifier = -1;
            current_radius_max = input_Radii[growing_cap_index];
        }
        int true_previous_left_patch_index  = (growing_cap_index == 1) ? 0 : 2*(growing_cap_index-1) - 1;
        std::vector<double> previous_left_patch_bounds = stripes.Patches[true_previous_left_patch_index].patch_boundaries();
        previous_patch_symmetric_boundary = abs(previous_left_patch_bounds[0]);
        
        current_radius_min = input_Radii[growing_cap_index];
        length_current_radius = (int) ((current_radius_max-current_radius_min)/d_current_radius) + 1;
        if(length_current_radius < 0){length_current_radius = 0.0;}
        Radius_loop_start_end = parallel_process->getLoopStartEnd_balanced_binary_split (length_current_radius, parallel_process->get_worker_color_for_level (parallelisation_lvl));
        
        if(current_radius_max < current_radius_min)
        {
            return 0;
        }
    }
    return 1;
}



int EvolveSphericalCap::identify_current_cap_growth_conditions_constant_modifiers()
{
    identify_growing_cap();
    
    double projected_radius_min;
    if (growing_cap_index == 0)
    {
        current_radius_min = 0.0;
        current_radius_max = Rg_max;
        d_current_radius = d_Rg;
        length_current_radius = (int) ((current_radius_max-current_radius_min)/d_current_radius) + 1;
        if(length_current_radius < 0){length_current_radius = 0.0;}
        Radius_loop_start_end = parallel_process->getLoopStartEnd_balanced_binary_split (length_current_radius, parallel_process->get_worker_color_for_level (parallelisation_lvl));
        current_theta = theta_good;
        clstr_centre_location_modifier = 0;
        
        if(current_radius_max < current_radius_min)
        {
            return 0;
        }
    }
    else
    {
        previous_cap_projected_radius = current_projected_radii_array[growing_cap_index - 1];
        previous_centre_to_left_patch_boundary_distance = current_centre_to_left_patch_boundary_distance[growing_cap_index - 1];
        double inside_sqrt = previous_cap_projected_radius*previous_cap_projected_radius - (previous_centre_to_left_patch_boundary_distance * previous_centre_to_left_patch_boundary_distance);
        
        if(inside_sqrt>=0)
        {projected_radius_min = sqrt(inside_sqrt);}
        else
        {
            printf("Calculation of minimum radius has negative sqrt for growing_cap_index=%d\n", growing_cap_index);
            abort();
        }
        
        if (growing_cap_index % 2 == 0) //good patch
        {
            d_current_radius = d_Rg;
            current_theta = theta_good;
            current_radius_max = Rg_max;
            if(starting_patch_is_bad == 1) //When this is true then the even patches are bad and all the input quantities are switch ed.
            {clstr_centre_location_modifier = -1;} //This is the constant modifier for good patches since more surface area on good patch is favorable
            else{clstr_centre_location_modifier = 1;}
            if(clstr_centre_location_modifier == -1) //This is for better load balancing
            {
                current_radius_max = previous_cap_projected_radius/sin(current_theta);
            }
            else
            {
                current_radius_max = Rg_max;
            }
//            current_radius_max = Rg_max;
        }
        else
        {
            d_current_radius = d_Rb;
            current_theta = theta_bad;
            
            if(starting_patch_is_bad == 1)
            {
                clstr_centre_location_modifier = 1;
            }
            else
            {
                clstr_centre_location_modifier = -1;
                
            } //This is the constant modifier for bad patches since less surface area on good patch is favorable
            if(clstr_centre_location_modifier == -1)  //This is for better load balancing
            {
                current_radius_max = previous_cap_projected_radius/sin(current_theta);
            }
            else
            {
                current_radius_max = Rb_max;
            }
//             current_radius_max = Rb_max;
            
        }
        int true_previous_left_patch_index  = (growing_cap_index == 1) ? 0 : 2*(growing_cap_index-1) - 1;
        std::vector<double> previous_left_patch_bounds = stripes.Patches[true_previous_left_patch_index].patch_boundaries();
        previous_patch_symmetric_boundary = abs(previous_left_patch_bounds[0]);
        
        current_radius_min = projected_radius_min/(double)sin(current_theta) ;
//        if(growing_cap_index==3)
//        {
//            printf("current_theta=%10.10f\t previous_cap_projected_radius = %10.10f\t previous_centre_to_left_patch_boundary_distance=%10.10f\t projected_radius_min=%10.10f\t current_radius_min=%10.10f\n", current_theta, previous_cap_projected_radius, previous_centre_to_left_patch_boundary_distance, projected_radius_min, current_radius_min)  ;
//        }
        
        
        length_current_radius = (int) ((current_radius_max-current_radius_min)/d_current_radius) + 1;
        if(length_current_radius < 0){length_current_radius = 0.0;}
        Radius_loop_start_end = parallel_process->getLoopStartEnd_balanced_binary_split (length_current_radius, parallel_process->get_worker_color_for_level (parallelisation_lvl));
        
        if(current_radius_max < current_radius_min)
        {
            return 0;
        }
    }
    return 1;
}


void EvolveSphericalCap::get_balanced_modifier_array(int size, bool isEvenBranch)
{
    clstr_centre_location_modifier_array_load_balancing.empty();
    clstr_centre_location_modifier_array_load_balancing.resize(size);
    for(int i=0; i<size; i++)
    {
        if(isEvenBranch)
        {
            clstr_centre_location_modifier_array_load_balancing[i] = pow(-1.0, i);
        }
        else
        {
            clstr_centre_location_modifier_array_load_balancing[i] = pow(-1.0, i+1);
            
        }
    }
}

void EvolveSphericalCap::calc_and_set_centre_to_left_patch_boundary_distance()
{
    std::vector<double> current_centre = current_cap_centres[growing_cap_index];
    int true_left_patch_index  = (growing_cap_index == 0) ? 0 : (2*(growing_cap_index) - 1);
    std::vector<double> previous_left_patch_bounds = stripes.Patches[true_left_patch_index].patch_boundaries();
    double left_bound = previous_left_patch_bounds[0];
    
    //One dimensional distance.
    current_centre_to_left_patch_boundary_distance[growing_cap_index] = (current_centre[0] - left_bound);
//    printf("current_centre_to_left_patch_boundary_distance[%d] = %10.5f\n", growing_cap_index, current_centre_to_left_patch_boundary_distance[growing_cap_index]);
}


bool EvolveSphericalCap::patch_bounds_crossed()
{
    /**/
    std::vector<int> stripes_bounds = stripes.monitor_cluster_spread(Cluster_shape_ptr);
    int left_patch_index, right_patch_index;
    if (growing_cap_index == 0)
    {
        left_patch_index = 0;
        right_patch_index = 0;
    }
    else
    {
        left_patch_index = 2*growing_cap_index - 1;
        right_patch_index = 2*growing_cap_index;
    }
    if (stripes_bounds[left_patch_index] != stripes_bounds[right_patch_index])
    { printf ("stripes bounds in check bounds are not symmetric\n"); abort();}
    else
    {
        if(stripes_bounds[left_patch_index] == 1 && stripes_bounds[right_patch_index] == 1)
        {
            
//            if(growing_cap_index==2)
//            {
//                printf("patch bounds crossed for growing_cap_index=%d projected_radius=%10.10f\t cap_centre_x=%10.10f\n", growing_cap_index, current_projected_radii_array[growing_cap_index], current_cap_centres[growing_cap_index][0]);
//                printf("current_projected_radii_array= ");
//                for(size_t t=0; t<current_projected_radii_array.size(); t++)
//                {
//                    printf("%10.10f\t", current_projected_radii_array[t]);
//                }
//                printf("\n");
//                printf("Cap centres size=%zu\t x-coord=\n", current_cap_centres.size());
//                for(size_t t=0; t<=growing_cap_index; t++)
//                {
//                    printf("%10.10f\t", current_cap_centres[t][0]);
//                }
//                printf("\n");
//                std::vector<double> shape_xy_spread = Cluster_shape_ptr->xy_spread();
//                printf("Cluster xy spread = [%10.5f %10.5f][%10.5f %10.5f]\n", shape_xy_spread[0], shape_xy_spread[1], shape_xy_spread[2], shape_xy_spread[3]);
//
//            }
            return true;
        }
        else
        {
            return false;
        }
    }
}

void EvolveSphericalCap::identify_growing_cap()
{
    
    //int size = (int)current_growing_cap_identifier.size();
    if(is_valid_growing_cap_identifier())
    {
         std::vector<int>::iterator it_zero = std::find(current_growing_cap_identifier.begin(), current_growing_cap_identifier.end(), 0);
        if(it_zero == current_growing_cap_identifier.begin())
        {
            growing_cap_index = 0;
        }
        else if (it_zero == current_growing_cap_identifier.end())
        {
           growing_cap_index = (int)current_growing_cap_identifier.size() - 1;
        }
        else
        {
            int dist = std::distance(current_growing_cap_identifier.begin(),it_zero);
            growing_cap_index = dist - 1;
            
        }
        
    }
    else
    {
        printf("current_growing_cap_identifier is not valid\n");
        abort();
    }
    
}



bool EvolveSphericalCap::is_valid_growing_cap_identifier ()
{
    if(current_growing_cap_identifier.empty())
    {
        printf("growing cap identified is empty\n");
        abort();
    }
    else
    {
        std::vector<int>::iterator it_zero = std::find(current_growing_cap_identifier.begin(), current_growing_cap_identifier.end(), 0);
        
        if(it_zero == current_growing_cap_identifier.begin())
        {
            return compare_vector_elements_to_value (it_zero, current_growing_cap_identifier.end(), 0);
        }
        else
        {
            bool all_zeroes_after_first_zero = compare_vector_elements_to_value (it_zero, current_growing_cap_identifier.end(), 0);
            bool all_ones_before_first_zero = compare_vector_elements_to_value (current_growing_cap_identifier.begin(), it_zero, 1);
            return (all_zeroes_after_first_zero && all_ones_before_first_zero);
        }
    }
}
   
void EvolveSphericalCap::identify_cap_type()
{
    if(growing_cap_index == 0)
    {
        cap_type = CENTRE_GOOD;
    }
    else if (growing_cap_index % 2 == 0)
    {
        cap_type =  GOOD;
    }
    else
    {
        cap_type = BAD;
    }
}



bool EvolveSphericalCap::check_breaking_condition()
{
    //returning false indicates breaking will NOT occur and true indicates that it WILL.
    int cap_identifier = growing_cap_index;
    if (stripes.surface_bounds_breach (Cluster_shape_ptr)) //This is checking the final
    {
        printf("total surface bounds have been breached\n");
        return true;
    }
    else if(cap_identifier == 0)
    {
        return false;
    }
    else
    {
        std::vector<double> current_centre = current_cap_centres[cap_identifier];
        double current_proj_radii = current_projected_radii_array[cap_identifier];
        std::vector<double> previous_centre = current_cap_centres[cap_identifier-1];
        double previous_proj_radii = current_projected_radii_array[cap_identifier-1];
        
        //The abs() value of centre differences accounts for both left and right images of the current cap.
        bool check_less_than = (current_proj_radii < previous_proj_radii + previous_centre[0] - current_centre[0]) ;
        bool check_greater_than = (current_proj_radii > previous_proj_radii - previous_centre[0] + current_centre[0]) ;

        if(!(check_less_than && check_greater_than))
        {
//            printf("check_less_than = %d check_greater_than=%d\n", check_less_than, check_greater_than);
//            printf("current_proj_radii=%10.10f\t previous_proj_radii=%10.10f\t previous_centre[0] =%10.10f\t current_centre[0]=%10.10f\n", current_proj_radii, previous_proj_radii, previous_centre[0], current_centre[0]);
            return true;
        }
        else {return false;}
    }
}

void EvolveSphericalCap::print_initial_variables ()
{
    printf("Initial vars\n");
    printf("parallelization lvl= %d\n", parallelisation_lvl);
    printf("growing_cap_index = %d\n", growing_cap_index);
    printf("cap_type = %d\n", cap_type);
    std::cout<<" theta_good = "<<theta_good<<std::endl;
    std::cout<<" theta_bad = "<<theta_bad<<std::endl;
    std::cout<<" d_Rg = "<<d_Rg<<std::endl;
    std::cout<<" d_Rb = "<<d_Rb<<std::endl;
    std::cout<<" Rg_max = "<<Rg_max<<std::endl;
    std::cout<<" Rb_max = "<<Rb_max<<std::endl;
    std::cout<<" delta = "<<delta<<std::endl;
    std::cout<<" Rho = "<<Rho<<std::endl;
    std::cout<<" num_patches = "<<num_patches<<std::endl;
}

void EvolveSphericalCap::print_growing_cap_identifier()
{
    printf("current_growing_cap_identifier = ");
    for(int i=0; i<(int)current_growing_cap_identifier.size(); i++)
    {
        printf("%d\t", current_growing_cap_identifier[i]);
    }
    printf("\n");
    
}

void EvolveSphericalCap::set_MC_comm()
{
    MPI_Comm MC_comm = parallel_process->get_this_lvl_branch_comm(parallelisation_lvl);
    mc_engine->set_branch_comm(MC_comm);
}

int EvolveSphericalCap::write_to_file_periodically()
{
    int written;
    if(parallel_process->is_current_lvl_root(parallelisation_lvl))
    {
        written = periodic_io->write_individually();
    }
    MPI_Comm this_lvl_branch_comm = parallel_process->get_this_lvl_branch_comm(parallelisation_lvl);
    MPI_Bcast (&written, 1, MPI_INT, 0,  this_lvl_branch_comm);
    if(written == 1)
    {
        clear_output_variables();
    }
    return 1; //periodic_io->gather_and_write());
}

int EvolveSphericalCap::write_profiling_periodically()
{
    int written;
    if(parallel_process->is_current_lvl_root(parallelisation_lvl))
    {
        printf("writing for rank=%d\t level=%d\n", myRank, parallelisation_lvl);
        written = periodic_io->write_profiling_individually(grid_points, cost_per_grid_point);
        if(written == 1)
        {
            clear_grid_points();
        }
    }
//    MPI_Comm this_lvl_branch_comm = parallel_process->get_this_lvl_branch_comm(parallelisation_lvl);
//    MPI_Bcast (&written, 1, MPI_INT, 0,  this_lvl_branch_comm);
//    if(written == 1)
//    {
//        clear_grid_points();
//    }
    return 1;
    //periodic_io->gather_and_write());
}

void EvolveSphericalCap::clear_output_variables()
{
    output_variables->clstr_centre_location_modifier_global_array.clear();
    output_variables->radii_global_array.clear();
    output_variables->Number_particles_global_array.clear();
    output_variables->Volume_global_array.clear();
    output_variables->SA_global_array.clear();
    output_variables->projected_SA_global_array.clear();
}



//Temporary Setters
void EvolveSphericalCap::set_cluster_shape_ptr(Shape* shape)
{
    Cluster_shape_ptr = shape;
}

void EvolveSphericalCap::set_growing_cap_index(int index)
{
    growing_cap_index = index;
}
