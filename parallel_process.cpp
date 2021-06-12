//
//  parallel_process.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 5/28/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#include "parallel_process.hpp"

ParallelProcess::ParallelProcess()
{
    
}

ParallelProcess::~ParallelProcess()
{
    
}

ParallelProcess::ParallelProcess(int branches_per_node, int Levels):
levels(Levels), n_branches_per_node(branches_per_node)
{
    MPI_Comm_rank (MPI_COMM_WORLD, &myRank);
    MPI_Comm_size (MPI_COMM_WORLD, &nProcs);
    last_level = levels-1; //must be declared before create_root_comms_per_level
    calc_worker_details_per_level();
    create_branch_comms_per_level();
    create_root_comms_per_level();
   
}


void ParallelProcess::calc_worker_details_per_level()
{
    worker_branch_colors_per_level.resize(levels); //This will have the color of each worker at each level.
    num_workers_per_groups_per_level.resize(levels);
    worker_root_colors_per_level.resize(levels);
    
    for(int i=0; i < levels; i++)
    {
        int group_size = pow(n_branches_per_node, i+1);
        num_workers_per_groups_per_level[i] = nProcs/group_size; //Ensure that they are multiples
        worker_branch_colors_per_level[i] = myRank/num_workers_per_groups_per_level[i];
        worker_root_colors_per_level[i] = ((myRank % num_workers_per_groups_per_level[i]) == 0) ? 0 : 1 ;
        
//        printf("num_workers_per_groups_per_level[%d]= %d\n", i, num_workers_per_groups_per_level[i]);
//        printf("level =%d\t rank = %d\t worker_branch_colors_per_level = %d\t worker_root_colors_per_level = %d\n", i, myRank, worker_branch_colors_per_level[i], worker_root_colors_per_level[i]);
    }
    
}

void ParallelProcess::create_branch_comms_per_level()
{
    branch_comm_levels.resize(levels);
    level_branch_rank.resize(levels);
    level_branch_size.resize(levels);
    
    for(int i=0; i < levels; i++)
    {
        MPI_Comm_split(MPI_COMM_WORLD, worker_branch_colors_per_level[i], myRank, &branch_comm_levels[i]); //color
        
        MPI_Comm_rank(branch_comm_levels[i], &level_branch_rank[i]);
        MPI_Comm_size(branch_comm_levels[i], &level_branch_size[i]);
    }
    
}

void ParallelProcess::create_root_comms_per_level()
{
    /* Creating another communicator containing the roots of each communicator so that V, SA..etc, data can be transferred from them to the original root */
    //Here last_lvl_root_color correspond to the colors of processes in each group of the last level.
    
    root_comm_levels.resize(levels);
    level_root_rank.resize(levels);
    level_root_size.resize(levels);
    
    for(int i=0; i < levels; i++)
    {
        MPI_Comm_split(MPI_COMM_WORLD, worker_root_colors_per_level[i], myRank, &root_comm_levels[i]);
        MPI_Comm_rank(root_comm_levels[i], &level_root_rank[i]);
        MPI_Comm_size(root_comm_levels[i], &level_root_size[i]);
    }
}

bool ParallelProcess::is_parent()
{
    return myRank==0 ;
}

bool ParallelProcess::is_last_lvl_root()
{
    return (worker_root_colors_per_level[last_level] == 0);
}


bool ParallelProcess::is_last_lvl_first_group()
{
    return (worker_branch_colors_per_level[last_level] == 0) ;
}

bool ParallelProcess::is_current_lvl_root(int level)
{
    return (worker_root_colors_per_level[level] == 0);
}

bool ParallelProcess::is_current_lvl_even_branch(int level)
{
    return (worker_branch_colors_per_level[level] % 2 == 0);
}


void ParallelProcess::free_MPI_comms()
{
//    MPI_Comm_free(&roots_comm);
    for(int i=0; i<levels; i++)
    {
        MPI_Comm_free (&branch_comm_levels[i]);
        MPI_Comm_free (&root_comm_levels[i]);
    }
}
