//
//  parallel_process.hpp
//  Clusters
//
//  Created by Sarwar Hussain on 5/28/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#ifndef parallel_process_hpp
#define parallel_process_hpp

#include <stdio.h>
#include <vector>
#include <string.h>
#include <cmath>
#include <mpi.h>


//Parallelisation structure :
/*
                    Loop for Good patch radius (Rg)
                                 |
                                 | (divide this into n_branches_per_node = 2)
                                /\
                               /  \
     cap_location modifier loop        cap_location modifier loop
        /\                                   /\
       /  \                                 /  \
 
 
 */

class ParallelProcess
{
public:
    
    ParallelProcess();
    ~ParallelProcess();
    ParallelProcess(int , int );
    bool is_parent();
    bool is_last_lvl_root();
    bool is_current_lvl_root(int);
    bool is_current_lvl_even_branch (int i);
    bool is_last_lvl_first_group();
    int get_worker_color_for_level (int level) {return worker_colors_per_level[level];}
    MPI_Comm get_last_lvl_branch_comm (){return branch_comm_levels[last_level];}
    
    
protected:
    //Parallelisation related variables.
    MPI_Comm roots_comm, branch_comm;
    std::vector<MPI_Comm> branch_comm_levels;
    int myRank, nProcs;
    int n_branches_per_node,  levels;
    int last_level;
    int last_lvl_root_color, roots_rank, roots_size;
    std::vector<int> worker_colors_per_level; //This will have the color of each worker at each level.
    std::vector<int> num_workers_per_groups_per_level;
    std::vector<int> level_branch_rank;
    std::vector<int> level_branch_size;
    
    //Create parallelisation
    void calc_worker_details_per_level();
    void create_branch_comms_per_level();
    void create_root_comms();
    
};
#endif /* parallel_process_hpp */
