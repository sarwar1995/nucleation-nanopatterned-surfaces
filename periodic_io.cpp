//
//  periodic_io.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 6/23/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#include "periodic_io.hpp"

PeriodicIO::PeriodicIO()
{
    MPI_Comm_rank (MPI_COMM_WORLD, &myRank);
    MPI_Comm_size (MPI_COMM_WORLD, &nProcs);
    V_SA_DataFile_per_root_rank = NULL;
}

PeriodicIO::~PeriodicIO()
{
    if (V_SA_DataFile_per_root_rank!= NULL)
    {
        fclose(V_SA_DataFile_per_root_rank);
        
    }
    if (profilingFile_per_root_rank != NULL)
    {
        fclose(profilingFile_per_root_rank);
    }
    
}

PeriodicIO::PeriodicIO(ParallelProcess* parallel, SphericalCapOutput* output, std::string filetag, int n_patches):
parallel_process(parallel), output_variables(output), tag(filetag), num_patches(n_patches)
{
    MPI_Comm_rank (MPI_COMM_WORLD, &myRank);
    MPI_Comm_size (MPI_COMM_WORLD, &nProcs);
    get_fileName();
    get_fileName_per_root_rank();
    get_profiling_fileName_per_root_rank();
    V_SA_DataFile_per_root_rank = NULL;
    profilingFile_per_root_rank = NULL;
}

int PeriodicIO::gather_and_write()
{
    int a = open_output_file_for_appending();
    int b = gather();
    int c = print_to_file();
    return (a && b && c);
}

int PeriodicIO::write_individually()
{
    int a = conditinoal_open_per_root_rank_output_file_for_appending();
    int b = print_per_root_rank_to_file();
    return (a && b);
}

int PeriodicIO::write_profiling_individually(std::vector<std::vector<double>>& grid_points, std::vector<int>& cost_per_grid_point)
{
    int a = conditinoal_open_per_root_rank_profiling_file_for_appending();
    int b = print_profiling_per_root(grid_points, cost_per_grid_point);
    return (a && b);
}


int PeriodicIO::gather()
{
    //Here we are looking at the last level because all previous level roots are also going to be last level roots as long as the parallelisation branching is same at each node, meaning same number of branches are created at each node. So that each level root will be a multiple of the branch_per_node variable.
    if(parallel_process->is_last_lvl_root())
    {
        int roots_size = parallel_process->get_last_lvl_roots_size();
        /*These counts are true for single element quantities per grid point like N, Volume, SA, dB and dG*/
        int counts[roots_size];
        int projected_SA_counts[roots_size];
        int radii_counts[roots_size];
        int clstr_centre_modifier_counts[roots_size];
        
        int nelements = (int)output_variables->Volume_global_array.size();
        int projected_SA_nelements = (int)output_variables->projected_SA_global_array.size();
        int radii_nelements = (int)output_variables->radii_global_array.size();
        int clstr_centre_modifier_nelements = (int)output_variables->clstr_centre_location_modifier_global_array.size();
        //        printf("nelements=%d, projected_SA_nelements=%d, radii_nelements=%d, clstr_centre_modifier_nelements=%d\n", nelements, projected_SA_nelements, radii_nelements, clstr_centre_modifier_nelements);
        
        MPI_Comm roots_comm = parallel_process->get_last_lvl_roots_comm();
        
        MPI_Gather(&nelements, 1, MPI_INT, &counts[0], 1, MPI_INT, 0, roots_comm);
        MPI_Gather(&projected_SA_nelements, 1, MPI_INT, &projected_SA_counts[0], 1, MPI_INT, 0, roots_comm);
        MPI_Gather(&radii_nelements, 1, MPI_INT, &radii_counts[0], 1, MPI_INT, 0, roots_comm);
        MPI_Gather(&clstr_centre_modifier_nelements, 1, MPI_INT, &clstr_centre_modifier_counts[0], 1, MPI_INT, 0, roots_comm);
        
        int disps[roots_size];
        int projected_SA_disps[roots_size];
        int radii_disps[roots_size];
        int clstr_centre_modifier_disps[roots_size];
        // Displacement for the first chunk of data - 0
        for (int i = 0; i < roots_size; i++)
        {
            //printf("i = %d\t counts = %d\t projected_SA_counts=%d\t radii_counts=%d\n", i, counts[i], projected_SA_counts[i], radii_counts[i]);
            disps[i] = (i > 0) ? (disps[i - 1] + counts[i - 1]) : 0;
            projected_SA_disps[i] = (i > 0) ? (projected_SA_disps[i - 1] + projected_SA_counts[i - 1]) : 0;
            radii_disps[i] = (i > 0) ? (radii_disps[i - 1] + radii_counts[i - 1]) : 0;
            clstr_centre_modifier_disps[i] = (i > 0) ? (clstr_centre_modifier_disps[i - 1] + clstr_centre_modifier_counts[i - 1]) : 0;
        }
        
        if(myRank == 0)
        {
            V_SA_size = disps[roots_size - 1] + counts[roots_size - 1];
            proj_SA_size = projected_SA_disps[roots_size - 1] + projected_SA_counts[roots_size - 1];
            radii_size = radii_disps[roots_size - 1] + radii_counts[roots_size - 1];
            clstr_centre_modifier_size = clstr_centre_modifier_disps[roots_size - 1] + clstr_centre_modifier_counts[roots_size - 1];
            
            gathered_output_variables.clstr_centre_location_modifier_global_array.resize(clstr_centre_modifier_size);
            gathered_output_variables.Number_particles_global_array.resize(V_SA_size);
            gathered_output_variables.Volume_global_array.resize(V_SA_size);
            gathered_output_variables.SA_global_array.resize(V_SA_size);
            gathered_output_variables.projected_SA_global_array.resize(proj_SA_size);
            gathered_output_variables.radii_global_array.resize(radii_size);
        }
        
        
        MPI_Gatherv(output_variables->Volume_global_array.data(), nelements, MPI_DOUBLE, gathered_output_variables.Volume_global_array.data(), &counts[0], &disps[0], MPI_DOUBLE, 0, roots_comm);
        
        MPI_Gatherv(output_variables->SA_global_array.data(), nelements, MPI_DOUBLE, gathered_output_variables.SA_global_array.data(), &counts[0], &disps[0], MPI_DOUBLE, 0, roots_comm);
        
        MPI_Gatherv(output_variables->Number_particles_global_array.data(), nelements, MPI_DOUBLE, gathered_output_variables.Number_particles_global_array.data(), &counts[0], &disps[0], MPI_DOUBLE, 0, roots_comm);
        
        MPI_Gatherv(output_variables->clstr_centre_location_modifier_global_array.data(), clstr_centre_modifier_nelements, MPI_INT, gathered_output_variables.clstr_centre_location_modifier_global_array.data(), &clstr_centre_modifier_counts[0], &clstr_centre_modifier_disps[0], MPI_INT, 0, roots_comm);
        
        MPI_Gatherv(output_variables->projected_SA_global_array.data(), projected_SA_nelements, MPI_DOUBLE, gathered_output_variables.projected_SA_global_array.data(), &projected_SA_counts[0], &projected_SA_disps[0], MPI_DOUBLE, 0, roots_comm);
        
        MPI_Gatherv(output_variables->radii_global_array.data(), radii_nelements, MPI_DOUBLE, gathered_output_variables.radii_global_array.data(), &radii_counts[0], &radii_disps[0], MPI_DOUBLE, 0, roots_comm);
        
    }
    if (gathered_output_variables_isempty())
    {
        return 0;
    }
    else {return 1;}
}

bool PeriodicIO::gathered_output_variables_isempty()
{
    return (gathered_output_variables.Volume_global_array.empty() || gathered_output_variables.SA_global_array.empty() || gathered_output_variables.Number_particles_global_array.empty() || gathered_output_variables.clstr_centre_location_modifier_global_array.empty() || gathered_output_variables.radii_global_array.empty() || gathered_output_variables.projected_SA_global_array.empty());
}

int PeriodicIO::print_to_file()
{
    if(myRank == 0)
    {
        if(V_SA_DataFile != NULL)
        {
//            gathered_output_variables.clstr_centre_location_modifier_global_array.resize(clstr_centre_modifier_size);
//            gathered_output_variables.Number_particles_global_array.resize(V_SA_size);
//            gathered_output_variables.Volume_global_array.resize(V_SA_size);
//            gathered_output_variables.SA_global_array.resize(V_SA_size);
//            gathered_output_variables.projected_SA_global_array.resize(proj_SA_size);
//            gathered_output_variables.radii_global_array.resize(radii_size);
            
            
            add_Volume_SA_parallel (gathered_output_variables.Number_particles_global_array, gathered_output_variables.radii_global_array, gathered_output_variables.Volume_global_array, gathered_output_variables.SA_global_array, gathered_output_variables.projected_SA_global_array, gathered_output_variables.clstr_centre_location_modifier_global_array, num_patches, V_SA_DataFile );
            
            fclose(V_SA_DataFile);
            return 1;
        }
        else
        {
            printf("File pointer is invalid\n") ;
            return 0;
        }
    }
    return 1;
}

int PeriodicIO::print_per_root_rank_to_file()
{
    if(V_SA_DataFile_per_root_rank != NULL)
    {
        
        //printf("print rank=%d\n", myRank);
        add_Volume_SA_parallel (output_variables->Number_particles_global_array, output_variables->radii_global_array, output_variables->Volume_global_array, output_variables->SA_global_array, output_variables->projected_SA_global_array, output_variables->clstr_centre_location_modifier_global_array, num_patches, V_SA_DataFile_per_root_rank );
        
        return 1;
    }
    else
    {
        printf("per root File pointer is invalid, myRank=%d\n",  myRank) ;
        return 0;
    }
    return 1;
}

int PeriodicIO::print_profiling_per_root(std::vector<std::vector<double>>& grid_points, std::vector<int>& cost_per_grid_point)
{
    if(profilingFile_per_root_rank != NULL)
    {
        if(grid_points.size() == cost_per_grid_point.size())
        {
            //printf("print rank=%d\n", myRank);
            for(size_t i=0; i<grid_points.size(); i++)
            {
                for(size_t j=0; j<grid_points[i].size(); j++)
                {
                    fprintf(profilingFile_per_root_rank, "%10.5f\t",grid_points[i][j]);
                }
                fprintf(profilingFile_per_root_rank, "%d\n", cost_per_grid_point[i]);
            }
            return 1;
        }
        else
        {
            printf("sizes of grid points and cost per grid point are not same rank=%d\n",  myRank) ;
            return 0;
        }
    }
    else
    {
        printf("per root File pointer is invalid, myRank=%d\n",  myRank) ;
        return 0;
    }
    
    return 1;
}


void PeriodicIO::get_fileName()
{
    if(parallel_process->is_parent())
    {
        V_SA_DataFileName.assign(tag + "_V_SA_data.txt");
//        std::cout<<"tag = "<<tag<<std::endl;
//        std::cout<<"V_SA_DataFileName = "<<V_SA_DataFileName<<std::endl;
    }
}

void PeriodicIO::get_fileName_per_root_rank()
{
    char rank_char[256];
    sprintf(rank_char, "%d", myRank);
    V_SA_DataFileName_per_root_rank.assign(tag + "_V_SA_data_rank_" + rank_char + ".txt");
//    std::cout<<"tag = "<<tag<<std::endl;
//    std::cout<<"V_SA_DataFileName_per_root_rank = "<<V_SA_DataFileName_per_root_rank<<std::endl;
}

void PeriodicIO::get_profiling_fileName_per_root_rank()
{
    char rank_char[256];
    sprintf(rank_char, "%d", myRank);
    profiling_FileName_per_root_rank.assign(tag + "_profiling_rank_" + rank_char + ".txt");
//    std::cout<<"tag = "<<tag<<std::endl;
//    std::cout<<"profiling_FileName_per_root_rank = "<<profiling_FileName_per_root_rank<<std::endl;
    
}

int PeriodicIO::open_output_file_for_appending()
{
    if(parallel_process->is_parent())
    {
        V_SA_DataFile = fopen(V_SA_DataFileName.c_str(), "a"); //OPEN IT FOR APPENDING
        //printf("ptr to file = %p\n", V_SA_DataFile);
        if(V_SA_DataFile == NULL)
        {
            printf("Error opening output file\n");
            return 0;
            //exit(1);
        }
        else {return 1;}
    }
    return 1;
}


int PeriodicIO::conditinoal_open_per_root_rank_output_file_for_appending()
{
    V_SA_DataFile_per_root_rank = (V_SA_DataFile_per_root_rank == NULL) ? fopen(V_SA_DataFileName_per_root_rank.c_str(), "a") : V_SA_DataFile_per_root_rank; //OPEN IT FOR APPENDING
        //printf("ptr to file = %p\n", V_SA_DataFile_per_root_rank);
        if(V_SA_DataFile_per_root_rank == NULL)
        {
            printf("Error opening output file\n");
            return 0;
            //exit(1);
        }
        else {return 1;}
    return 1;
}

int PeriodicIO::conditinoal_open_per_root_rank_profiling_file_for_appending()
{
    
    profilingFile_per_root_rank = (profilingFile_per_root_rank == NULL) ? fopen(profiling_FileName_per_root_rank.c_str(), "a") : profilingFile_per_root_rank; //OPEN IT FOR APPENDING
    //printf("ptr to file = %p\n", V_SA_DataFile_per_root_rank);
    if(profilingFile_per_root_rank == NULL)
    {
        printf("Error opening profiling file\n");
        return 0;
        //exit(1);
    }
    else {return 1;}
    return 1;
}
