//
//  MC_parallel.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 2/14/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#include "MC_parallel.hpp"

MC::MC()
{}

MC::MC(MPI_Comm branch_comm_input):
n_inside_virtual(0.0), n_virtual_points(0.0), n_near_surf_virtual(0.0)
{
    branch_comm = branch_comm_input;
    MPI_Comm_rank (branch_comm, &branch_rank);
    MPI_Comm_size (branch_comm, &branch_size);
    
    MPI_Comm_rank (MPI_COMM_WORLD, &myRank);
    MPI_Comm_size (MPI_COMM_WORLD, &nProcs);
    
    box.resize(3);
    for(size_t i=0; i<3; i++){
        box[i].resize(2);
        for(size_t j=0; j<2; j++){
            if(j==0){box[i][j] = 0;}
            else{box[i][j] = 1;}
        }
    }
    virtual_box = box;
    n_points = 1e6;
    x_points.resize(n_points);
    points_chunk_per_procs = (int)(n_points/branch_size) ;
    if(n_points % branch_size != 0){printf("Warning: number of points=%d are not perfectly divisible by number of processes=%d\n",n_points, branch_size); abort();}
    
}

MC::~MC()
{
    x_points.clear();
    box.clear();
    virtual_box.clear();
}

//MCbox has the start and end boundaries of the box in all three dimensions
MC::MC(int num_points, std::vector<std::vector<double> >& MCbox, int aSeed[3], MPI_Comm branch_comm_input, double del):
delta(del), n_inside_virtual(0.0), n_virtual_points(0.0), n_near_surf_virtual(0.0)
{
    
    branch_comm = branch_comm_input;
    MPI_Comm_rank (branch_comm, &branch_rank);
    MPI_Comm_size (branch_comm, &branch_size);
    
    MPI_Comm_rank (MPI_COMM_WORLD, &myRank);
    MPI_Comm_size (MPI_COMM_WORLD, &nProcs);
    
    n_points = num_points;
    box = MCbox;
    virtual_box = box;
    x_points.resize(n_points);
    seed[0] = aSeed[0];
    seed[1] = aSeed[1];
    seed[2] = aSeed[2];
    points_chunk_per_procs = (int)(n_points/branch_size) ;
//    printf("branch_rank = %d\t branch_size = %d\n", branch_rank, branch_size);
//    printf("points_chunk_per_procs = %d\n", points_chunk_per_procs);
      if(n_points % branch_size != 0){printf("Warning: Initial number of points=%d are not perfectly divisible by number of processes=%d\n",n_points, branch_size); abort();}
    generate_points();
}

MC::MC(int num_points, std::vector<std::vector<double> >& MCbox, int aSeed[3], int onLattice, MPI_Comm branch_comm_input, double del):
delta(del), n_inside_virtual(0.0), n_virtual_points(0.0), n_near_surf_virtual(0.0)
{
    if(onLattice == 0)
    {
        MC(num_points, MCbox, aSeed, branch_comm_input, del);
    }
    else //Right now only have functionality for cubic lattice
    {
        branch_comm = branch_comm_input;
        MPI_Comm_rank (branch_comm, &branch_rank);
        MPI_Comm_size (branch_comm, &branch_size);
        
        MPI_Comm_rank (MPI_COMM_WORLD, &myRank);
        MPI_Comm_size (MPI_COMM_WORLD, &nProcs);
        
        
        n_points = num_points;
        n_virtual_points = 0.0;
        box = MCbox;
        virtual_box = box;
        //Here the number of points will actually change for the cubic lattic, the way it is set up currently and so that will create problems.
        x_points.resize(n_points);
        seed[0] = aSeed[0];
        seed[1] = aSeed[1];
        seed[2] = aSeed[2];
        points_chunk_per_procs = (int)(n_points/branch_size) ;
        printf("branch_rank = %d\t branch_size = %d\n", branch_rank, branch_size);
        printf("points_chunk_per_procs = %d\n", points_chunk_per_procs);
        if(n_points % branch_size != 0){printf("Warning: Initial number of points=%d are not perfectly divisible by number of processes=%d\n",n_points, branch_size); abort();}
        generate_points_on_lattice();
    }
    
}

void MC::update_points_chunk_per_procs()
{
    points_chunk_per_procs = (int)(n_points/branch_size) ;
    //printf("updated points_chunk_per_procs = %d\n", points_chunk_per_procs);
    
}

void MC::update_box (std::vector<std::vector<double>> new_box)
{
    box = new_box;
    virtual_box = box;
}

void MC::update_virtual_box (std::vector<std::vector<double>> new_box)
{
    virtual_box = new_box;
}


void MC::add_points(std::vector<std::vector<double>> new_box, int direction_to_expand_in, int extra_points_per_added_box)
{
    //printf("Inside mc add points. new_box=[%10.5f\t %10.5f\t %10.5f\t %10.5f\t %10.5f\t %10.5f\t], direction=%d, extra_points_per_added_box=%d\n", new_box[0][0], new_box[0][1], new_box[1][0], new_box[1][1], new_box[2][0], new_box[2][1], direction_to_expand_in, extra_points_per_added_box);
    std::mt19937 mt_engine_x (seed[0]); //I am using the same seed here but maybe for generating more points we should a different seed than the original one.
    std::mt19937 mt_engine_y (seed[1]);
    std::mt19937 mt_engine_z (seed[2]);
    std::vector<std::vector<double>> added_box_bounds(2);
    added_box_bounds[0].resize(2);
    added_box_bounds[1].resize(2);
    added_box_bounds[0][0] = new_box[direction_to_expand_in][0];
    added_box_bounds[0][1] = box[direction_to_expand_in][0]; //The original lower bound in the expanded direction is now the upper bound
    added_box_bounds[1][0] = new_box[direction_to_expand_in][1];
    added_box_bounds[1][1] = box[direction_to_expand_in][1];
    
//    if(myRank==0){
//        printf("box=[%10.2f %10.2f][%10.2f %10.2f][%10.2f %10.2f]\n", box[0][0], box[0][1], box[1][0], box[1][1], box[2][0], box[2][1]);
//        printf("direction_to_expand_in = %d\t added_box_bounds = [%10.3f %10.3f][%10.3f %10.3f]\t extra_points_per_added_box=%d\n", direction_to_expand_in, added_box_bounds[0][0], added_box_bounds[0][1], added_box_bounds[1][0], added_box_bounds[1][1], extra_points_per_added_box);
//        printf("new_box=[%10.2f %10.2f][%10.2f %10.2f][%10.2f %10.2f]\n", new_box[0][0], new_box[0][1], new_box[1][0], new_box[1][1], new_box[2][0], new_box[2][1]);
//    }
    
    if(direction_to_expand_in == 2)
    {
        std::uniform_real_distribution<double> dist_x(box[0][0], box[0][1]); //a is included and b is not dist_x [a,b)
        std::uniform_real_distribution<double> dist_y(box[1][0], box[1][1]);
        std::uniform_real_distribution<double> dist_z(added_box_bounds[1][0], added_box_bounds[1][1]); //Since only upper box will be added, so just one dist_z
        int new_total_points = n_points + extra_points_per_added_box;
        for(size_t i=0; i<extra_points_per_added_box; i++)
        {
            std::vector<double> new_point(3,0.0);
            new_point[0] = dist_x(mt_engine_x);
            new_point[1] = dist_y(mt_engine_y);
            new_point[2] = dist_z(mt_engine_z);
            x_points.push_back(new_point);
        }
        n_points = new_total_points;
    }
    else
    {
        // when direction_to_expand_in != 2 then we will have two boxes added.
        std::vector<std::uniform_real_distribution<double>> distributions(3);
        
        for(int n=0; n<2; n++) //This n is the two symmetric boxes that are being added
        {
            for(int i=0; i<3; i++)
            {
                if(i != direction_to_expand_in)
                {
                    std::uniform_real_distribution<double> dist(box[i][0], box[i][1]);
                    distributions[i] = dist ;
                }
                else
                {
                    
                    std::uniform_real_distribution<double> dist(added_box_bounds[n][0], added_box_bounds[n][1]);
                    distributions[i] = dist ;
                }
                
            }
            for(size_t i=0; i<extra_points_per_added_box; i++)
            {
                std::vector<double> new_point(3,0.0);
                new_point[0] = distributions[0](mt_engine_x);
                new_point[1] = distributions[1](mt_engine_y);
                new_point[2] = distributions[2](mt_engine_z);
                x_points.push_back(new_point);
            }
            
        }
        
        int new_total_points = n_points + 2*extra_points_per_added_box;
        n_points = new_total_points;
    }
    
    update_points_chunk_per_procs();
    update_box(new_box);
    //This MPI_Barrier at the end of add_points ensures that all processes reach this point i.e. addition of extra points is complete for all processes before we proceed forward.
    //MPI_Barrier(branch_comm);
}

void MC::generate_points()
{
    //Using the built in mt19937 rnd engine
    std::mt19937 mt_engine_x (seed[0]); //deterministic seed. [0 2^64)
    std::mt19937 mt_engine_y (seed[1]);
    std::mt19937 mt_engine_z (seed[2]);
    std::uniform_real_distribution<double> dist_x(box[0][0], box[0][1]); //both range numbers are inclusive
    std::uniform_real_distribution<double> dist_y(box[1][0], box[1][1]);
    std::uniform_real_distribution<double> dist_z(box[2][0], box[2][1]);
    
    for(size_t i=0; i<n_points; i++)
    {
        x_points[i].resize(3);
        x_points[i][0] = dist_x(mt_engine_x);
        x_points[i][1] = dist_y(mt_engine_y);
        x_points[i][2] = dist_z(mt_engine_z);
    }
    
}

void MC::generate_points_on_lattice()
{
    std::vector<int> Seed(3,0);
    Seed[0] = seed[0];
    Seed[1] = seed[1];
    Seed[2] = seed[2];
    CubicLattice cubic_lattice (n_points, box, Seed);
    cubic_lattice.CalcTranslationVector();
    cubic_lattice.GenerateLattice();
    n_points = cubic_lattice.GetTrueTotalPoints();
    x_points = cubic_lattice.get_lattice_points();
}


std::vector<double> MC::getMeasures (int n_inside, int n_near_surf)
{
    std::vector<double> measures(2,0.0);
    //if(myRank==0){printf("Inside MC get measures n_points=%d\t n_inside = %d\t n_near_surf=%d\n",n_points, n_inside, n_near_surf);}
    double fraction_volume =  ((double)n_inside/(double)n_points);
    //printf("fraction_volume = %10.10f\n",fraction_volume);
    double fraction_SA =  ((double)n_near_surf/(double)n_points);
    double volume = fraction_volume * box_volume() ;
    double volume_near_surf = fraction_SA * box_volume() ;
    double SA = volume_near_surf/(2.0*delta) ;  // (Volume/thickness) = SA
    measures[0] = volume;
    measures[1] = SA;
    return  measures;
}

std::vector<double> MC::getMeasures_virtual (int n_inside, int n_near_surf)
{
    std::vector<double> measures(2,0.0);
    int n_total_points = n_points + n_virtual_points;
//    if(myRank==0){printf("n_total_points=%d\t n_points=%d\t n_virtual_points=%d\t n_inside=%d\t n_near_surf=%d\n", n_total_points, n_points, n_virtual_points, n_inside, n_near_surf);}
    //printf("n_inside = %d\t n_near_surf=%d\n",n_inside, n_near_surf);
    double fraction_volume =  ((double)n_inside/(double)n_total_points);
    //printf("fraction_volume = %10.10f\n",fraction_volume);
    double fraction_SA =  ((double)n_near_surf/(double)n_total_points);
    double volume = fraction_volume * virtual_box_volume() ; //box_volume is also updated in the add_and_loop_points function
    double volume_near_surf = fraction_SA * virtual_box_volume() ;
    double SA = volume_near_surf/(2.0*delta) ;  // (Volume/thickness) = SA
    measures[0] = volume;
    measures[1] = SA;
    return  measures;
}




void MC::addPoint(std::vector<double> &destination, std::vector<double> &point)
{
    destination.push_back(point[0]);
    destination.push_back(point[1]);
    destination.push_back(point[2]);
}

std::vector<int> MC::loopPoints (Shape* cluster)
{
    /* NO INSERTION IN INTERIOR OR SURFACE POINTS HERE. DO NOT USE UPDATE WITH THIS*/
    
    std::vector<int> int_surf_size(2,-1);
    std::vector <double> point (3, 0.0);

    int count=0;
    int interior_count = 0;
    int surface_count = 0;
    int start, end;
    
    if(n_points % branch_size != 0)
    {
        int extra_points = n_points % branch_size;
        //Adding the remaining extra points to the last branch rank
        if(branch_rank == branch_size - 1)
        {
            start = branch_rank*points_chunk_per_procs;
            end = (branch_rank+1)*points_chunk_per_procs - 1;
            end = end + extra_points ;
        }
        else
        {
            start = branch_rank*points_chunk_per_procs;
            end = (branch_rank+1)*points_chunk_per_procs - 1;
        }
    }
    else
    {
        start = branch_rank*points_chunk_per_procs;
        end = (branch_rank+1)*points_chunk_per_procs - 1;
        
    }
    
    for (int i=start ; i <= end; i++)
    {
        count++;
        point[0] = x_points[i][0];
        point[1] = x_points[i][1];
        point[2] = x_points[i][2];
        int isinside = cluster->isInside(point);
        if(isinside)
        {
            interior_count++;
            //addPoint(interior_points_process, point);
        }
        
        int nearSurf = cluster->nearSurface(point, delta);
        if(nearSurf)
        {
            surface_count++;
            //addPoint(surface_points_process, point);
        }
    }

    /* This code is taken from stakeoverflow (https://stackoverflow.com/questions/12080845/mpi-receive-gather-dynamic-vector-length)*/
    
    //MPI_Barrier(branch_comm);
    // Each process tells the root how many elements it holds
    int interior_counts[branch_size];
    int surface_counts[branch_size];
    MPI_Gather(&interior_count, 1, MPI_INT, &interior_counts[0], 1, MPI_INT, 0, branch_comm);
    MPI_Gather(&surface_count, 1, MPI_INT, &surface_counts[0], 1, MPI_INT, 0, branch_comm);
    
    if(branch_rank == 0)
    {
        int_surf_size[0] = 0;
        int_surf_size[1] = 0;
        for(int w=0; w<branch_size; w++)
        {
            int_surf_size[0] = int_surf_size[0] + interior_counts[w];
            int_surf_size[1] = int_surf_size[1] + surface_counts[w];
        }
        
        if(int_surf_size[0] > n_points || int_surf_size[1] > n_points)
        {
            printf("num of interior or surface points are larger than total points; interior count=%d\t surface count=%d\n", int_surf_size[0], int_surf_size[1]);
            abort();
        }
    }
    
    return int_surf_size;
    
}

std::vector<int> MC::gather_points(int& interior_count, int& surface_count)
{
    std::vector<int> int_surf_size(2, 0.0);
    MPI_Allreduce(&interior_count, &int_surf_size[0], 1, MPI_INT, MPI_SUM, branch_comm);
    MPI_Allreduce(&surface_count, &int_surf_size[1], 1, MPI_INT, MPI_SUM, branch_comm);
    return int_surf_size;
}


void MC::check_point(Shape* cluster, std::vector<double>& point, int& interior_count, int& surface_count)
{
    int isinside = cluster->isInside(point);
    if(isinside)
    {
        interior_count++;
    }
    
    int nearSurf = cluster->nearSurface(point, delta);
    if(nearSurf)
    {
        surface_count++;
    }
}

bool MC::within_box(std::vector<double>& point)
{
    bool within_x = point[0] <= box[0][1] && point[0] >= box[0][0];
    bool within_y = point[1] <= box[1][1] && point[1] >= box[1][0];
    bool within_z = point[2] <= box[2][1] && point[2] >= box[2][0];
    
    return (within_x && within_y && within_z);
}

void MC::add_and_loop_points(Shape* cluster, std::vector<std::vector<double>> new_box, int direction_to_expand_in, int extra_points_per_added_box)
{
   // printf("Inside add and loop points\n");
    std::vector<int> int_surf_size(2, 0.0);
    
    int count=0;
    int interior_count = 0;
    int surface_count = 0;
    int start, end;
    
    std::mt19937 mt_engine_x (seed[0] + myRank + branch_rank); //I am using the same seed here but maybe for generating more points we should a different seed than the original one.
    std::mt19937 mt_engine_y (seed[1] + myRank + branch_rank);
    std::mt19937 mt_engine_z (seed[2] + myRank + branch_rank);
    std::vector<std::vector<double>> added_box_bounds(2);
    added_box_bounds[0].resize(2);
    added_box_bounds[1].resize(2);
    added_box_bounds[0][0] = new_box[direction_to_expand_in][0];
    added_box_bounds[0][1] = virtual_box[direction_to_expand_in][0]; //The original lower bound in the expanded direction is now the upper bound
    added_box_bounds[1][0] = new_box[direction_to_expand_in][1]; //boundary going from new box to box
    added_box_bounds[1][1] = virtual_box[direction_to_expand_in][1];
    
    int n_extra_points_per_box_per_rank;
    if(extra_points_per_added_box % branch_size != 0)
    {
        n_extra_points_per_box_per_rank = (branch_rank == branch_size-1) ? (int)(extra_points_per_added_box/branch_size) + (extra_points_per_added_box % branch_size) : (int)(extra_points_per_added_box/branch_size);
        //printf("Imperfect division n_extra_points_per_box_per_rank = %d myRank=%d\n", n_extra_points_per_box_per_rank, myRank);
        
    }
    else
    {
        n_extra_points_per_box_per_rank = (int)(extra_points_per_added_box/branch_size);
        //printf("perfect division n_extra_points_per_box_per_rank = %d myRank=%d branch_size=%d\n", n_extra_points_per_box_per_rank, myRank, branch_size);
    }
//    std::vector<std::vector<double>> added_box_bounds_per_rank(2);
//    added_box_bounds_per_rank[0].resize(2);
//    added_box_bounds_per_rank[1].resize(2);
//
//
//    double signed_length_lower = added_box_bounds[0][1] - added_box_bounds[0][0];
//    double signed_length_upper = added_box_bounds[1][1] - added_box_bounds[1][0];
//    double d_lower_regions = (double)(signed_length_lower/(double)branch_size);
//    double d_upper_regions = (double)(signed_length_upper/(double)branch_size);
//
//
//    added_box_bounds_per_rank[0][0] = added_box_bounds[0][0] + branch_rank*d_lower_regions;
//    added_box_bounds_per_rank[0][1] = added_box_bounds[0][0] + (branch_rank+1)*d_lower_regions;
//    added_box_bounds_per_rank[1][0] = added_box_bounds[1][0] + branch_rank*d_upper_regions;
//    added_box_bounds_per_rank[1][1] = added_box_bounds[1][0] + (branch_rank+1)*d_upper_regions;
    
//    if(myRank==0){
//        printf("virtual_box=[%10.2f %10.2f][%10.2f %10.2f][%10.2f %10.2f]\n", virtual_box[0][0], virtual_box[0][1], virtual_box[1][0], virtual_box[1][1], virtual_box[2][0], virtual_box[2][1]);
//        printf("direction_to_expand_in = %d\t added_box_bounds = [%10.3f %10.3f][%10.3f %10.3f]\t n_extra_points_per_box_per_rank=%d extra_points_per_added_box=%d\n", direction_to_expand_in, added_box_bounds[0][0], added_box_bounds[0][1], added_box_bounds[1][0], added_box_bounds[1][1], n_extra_points_per_box_per_rank, extra_points_per_added_box);
//        printf("new_box=[%10.2f %10.2f][%10.2f %10.2f][%10.2f %10.2f]\n", new_box[0][0], new_box[0][1], new_box[1][0], new_box[1][1], new_box[2][0], new_box[2][1]);
//    }
    
//    printf("per rank added_box_bounds = [%10.3f %10.3f][%10.3f %10.3f], myRank=%d\n", added_box_bounds_per_rank[0][0], added_box_bounds_per_rank[0][1], added_box_bounds_per_rank[1][0], added_box_bounds_per_rank[1][1], myRank);
//    if(n_inside_virtual != 0 || n_near_surf_virtual != 0)
//    {
//        printf("n_inside_virtual=%d\t n_near_surf_virtual=%d\n are not zero before adding points\n", n_inside_virtual, n_near_surf_virtual);
//        abort();
//    }
 
    if(direction_to_expand_in == 2)
    {
        std::uniform_real_distribution<double> dist_x(virtual_box[0][0], virtual_box[0][1]); //a is included and b is not dist_x [a,b)
        std::uniform_real_distribution<double> dist_y(virtual_box[1][0], virtual_box[1][1]);
        std::uniform_real_distribution<double> dist_z(added_box_bounds[1][0], added_box_bounds[1][1]); //Since only upper box will be added, so just one dist_z
        for(size_t i=0; i<n_extra_points_per_box_per_rank ; i++) //extra_points_per_added_box
        {
            std::vector<double> new_point(3,0.0);
            new_point[0] = dist_x(mt_engine_x);
            new_point[1] = dist_y(mt_engine_y);
            new_point[2] = dist_z(mt_engine_z);
            if(within_box(new_point))
            {
                printf("new_point=[%10.2f %10.2f %10.2f] is within box\n",  new_point[0],  new_point[1],  new_point[2]);
                abort();
            }
            check_point(cluster, new_point, interior_count, surface_count);
        }
        n_virtual_points += extra_points_per_added_box; //points that have not actually been stored in the "x" matrix
        //MPI_Barrier(branch_comm);
        
        int_surf_size = gather_points(interior_count, surface_count);
        n_inside_virtual += int_surf_size[0];
        n_near_surf_virtual += int_surf_size[1];
        //printf("direction_to_expand_in = 2 n_points=%d n_virtual_points=%d interior_count = %d\t surface_count=%d int_surf_size=[%d %d] n_inside_virtual=%d n_near_surf_virtual=%d myRank=%d branch_rank=%d\n", n_points, n_virtual_points, interior_count, surface_count, int_surf_size[0], int_surf_size[1], n_inside_virtual, n_near_surf_virtual, myRank, branch_rank);
    }
    else
    {
        // when direction_to_expand_in != 2 then we will have two boxes added.
        std::vector<std::uniform_real_distribution<double>> distributions(3);
        
        for(int n=0; n<2; n++) //This n is the two symmetric boxes that are being added
        {
            for(int i=0; i<3; i++)
            {
                if(i != direction_to_expand_in)
                {
                    std::uniform_real_distribution<double> dist(virtual_box[i][0], virtual_box[i][1]);
                    distributions[i] = dist ;
                }
                else
                {
                    std::uniform_real_distribution<double> dist(added_box_bounds[n][0], added_box_bounds[n][1]);
                    distributions[i] = dist ;
                }
                
            }
            for(size_t i=0; i<n_extra_points_per_box_per_rank; i++)
            {
                std::vector<double> new_point(3,0.0);
                new_point[0] = distributions[0](mt_engine_x);
                new_point[1] = distributions[1](mt_engine_y);
                new_point[2] = distributions[2](mt_engine_z);
                if(within_box(new_point))
                {
                    printf("new_point=[%10.2f %10.2f %10.2f] is within box\n",  new_point[0],  new_point[1],  new_point[2]);
                    abort();
                }
                check_point(cluster, new_point, interior_count, surface_count);
            }
            n_virtual_points += extra_points_per_added_box; //points that have not actually been stored in the "x" matrix
            //printf("direction_to_expand_in =%d n_points=%d n_virtual_points = %d\t interior_count = %d\t surface_count=%d int_surf_size=[%d %d] n_inside_virtual=%d n_near_surf_virtual=%d myRank=%d branch_rank=%d\n", direction_to_expand_in, n_points, n_virtual_points, interior_count, surface_count, int_surf_size[0], int_surf_size[1], n_inside_virtual, n_near_surf_virtual, myRank, branch_rank);
        }
        //MPI_Barrier(branch_comm);
        int_surf_size = gather_points(interior_count, surface_count);
        n_inside_virtual += int_surf_size[0];
        n_near_surf_virtual += int_surf_size[1];
    }
    
    update_virtual_box(new_box);
    //This MPI_Barrier at the end of add_points ensures that all processes reach this point i.e. addition of extra points is complete for all processes before we proceed forward.
}

std::vector<double> MC::calc_volume_SA(Shape* cluster)
{
    
    //Calculating volume and SA for the first time. Currently doing this whenever Rg is increased.
    printf("Inside calc_volume_SA MC\n");
    std::vector<double> measures(2,-1.0);
    n_inside = 0;
    n_near_surf = 0;
    std::vector<int> sizes = loopPoints(cluster);
    
    if(branch_rank == 0)
    {
        n_inside = sizes[0];
        n_near_surf = sizes[1];
        
        /*NOT THROWING THIS ERROR FOR TESTING SPHERICAL CAP WITH MC*/
//        if(n_inside==0 || n_near_surf==0)
//        {
//            throw std::logic_error("n_inside or n_near_surf are zero");
//        }
        
        measures = getMeasures (n_inside, n_near_surf);
    }
    return measures;
}

std::vector<double> MC::calc_volume_SA_with_virtual_points(Shape* cluster)
{
    std::vector<double> measures(2,-1.0);
    n_inside = 0;
    n_near_surf = 0;
    std::vector<int> sizes = loopPoints(cluster);
    
    if(branch_rank == 0)
    {
//        if(n_inside_virtual==0 || n_near_surf_virtual==0)
//        {
//            printf("n_inside_virtual=%d n_near_surf_virtual=%d cannot be zero for branch rank=0 in calc_volume_SA_with_virtual_points, sizes=[%d %d] myRank=%d", n_inside_virtual, n_near_surf_virtual, sizes[0], sizes[1], myRank);
//            //abort();
//        }
        n_inside = sizes[0] + n_inside_virtual;
        n_near_surf = sizes[1] + n_near_surf_virtual;
        measures = getMeasures_virtual (n_inside, n_near_surf);
    }
//    if(myRank==0){printf("n_inside_virtual=%d n_near_surf_virtual=%d sizes=[%d %d] myRank=%d\n", n_inside_virtual, n_near_surf_virtual, sizes[0], sizes[1], myRank);}
    reset_virtual();
    return measures;
}

void MC::reset_virtual()
{
    n_inside_virtual = 0.0;
    n_near_surf_virtual = 0.0;
    n_virtual_points = 0.0;
    virtual_box = box;
}

void MC::set_branch_comm (MPI_Comm new_branch_comm)
{
    int compare_result;
    branch_comm = new_branch_comm;
    MPI_Comm_compare(new_branch_comm, branch_comm, &compare_result);
//    printf("result is identical = %d\n", compare_result==MPI_IDENT);
    
    
    MPI_Comm_rank (branch_comm, &branch_rank);
    MPI_Comm_size (branch_comm, &branch_size);
    update_points_chunk_per_procs();
}

//This requires making sure that the growth (delta r) is smaller than the cell size so that neigbs indeed have all the extra points that need to be added or deleted
//std::vector<double> MC::update_volume_SA(Shape* new_cluster, CellList* NgbCellList, double delta)
//{
//    printf("Initial surface and interior_points = %d\t%d\n", (int)surface_points.size(), (int)interior_points.size());
//    //number of interior points will be added
//    n_near_surf = 0;    //number of surface points will change
//    std::vector<double> point(3,0.0);
//    std::vector<double> measures(2,0.0);
//    
//    //Updating surface and interior points count and vectors
//    std::vector<std::vector<double> > ngb_surface_points ; //Ngbs of existing surface points, excluding the points themselves
//    std::set<std::vector<double> > new_surface_points; //Updated set of surface points
//    auto t1 = std::chrono::high_resolution_clock::now();
//    ngb_surface_points = NgbCellList->calc_neighbors (surface_points);
//    auto t2 = std::chrono::high_resolution_clock::now();
//    auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
//    printf("duration for finding neighbors of surface points = %lld\n",duration_ms.count());
//    
//    
//    printf("number of ngbs of surface_points = %d\n", (int)ngb_surface_points.size());
//    for(size_t i=0; i<ngb_surface_points.size(); i++)
//    {
//        point = ngb_surface_points[i];
//        int isinside = new_cluster->isInside(point);
//        //This check is only being made for ngbs of surface points because all existing surface points that are already inside are expected to be part of interior_points already
//        if(isinside)
//        {
//            interior_points.insert(point);
//        }
//        int nearSurf = new_cluster->nearSurface(point, delta);
//        if(nearSurf)
//        {
//            new_surface_points.insert(point);
//        }
//    }
//    //checking those from old surface points and adding them to new surface points
//    for(std::set<std::vector<double> >::iterator it=surface_points.begin(); it!=surface_points.end();it++)
//    {
//        point = *it; //surface_points[i];
//        int nearSurf = new_cluster->nearSurface(point, delta);
//        if(nearSurf)
//        {
//            new_surface_points.insert(point);
//            
//        }
//    }
//    
//    //clear and update surface_points vector
//    surface_points.clear();
//    surface_points = new_surface_points;
//    n_inside = (int) interior_points.size();
//    n_near_surf = (int) surface_points.size();
//    printf("New surface and interior_points = %d\t%d\n",n_near_surf, n_inside);
//    measures = getMeasures (n_inside, n_near_surf, delta);
//    return measures;
//}


void MC::print_points(FILE * pointsFile)
{
    for(size_t i=0; i<x_points.size(); i++)
    {
        fprintf(pointsFile,"%d\t%10.10f\t%10.10f\t%10.10f\n",(int) i,x_points[i][0],x_points[i][1],x_points[i][2]);
    }
}

void MC::print_surf_points(Shape* cluster, FILE * SurfpointsFile)
{
    if(branch_rank == 0)
    {
        std::set<std::vector<double> >::iterator it;
        std::vector <double> point (3, 0.0);
        for (int i=0 ; i < n_points; i++)
        {
            point[0] = x_points[i][0];
            point[1] = x_points[i][1];
            point[2] = x_points[i][2];
            int isinside = cluster->isInside(point);
            if(isinside)
            {
                //interior_count++;
                //addPoint(interior_points_process, point);
//                interior_points.insert(point);
            }

            int nearSurf = cluster->nearSurface(point, delta);
            if(nearSurf)
            {
                //surface_count++;
                //addPoint(surface_points_process, point);
                surface_points.insert(point);

            }
        }
        int count = 0;
        for(it=surface_points.begin(); it!=surface_points.end(); it++)
        {
            count++;
            std::vector<double> point = *it ;
            fprintf(SurfpointsFile,"%d\t%10.10f\t%10.10f\t%10.10f\n",count, point[0], point[1], point[2]);
        }
    }

}


//void MC::print_volume_points(FILE* VolumePointsFile)
//{
//    std::set<std::vector<double> >::iterator it;
//    int count=0;
//
//    for(it=interior_points.begin(); it!=interior_points.end(); it++)
//    {
//        count++;
//        std::vector<double> point = *it ;
//        fprintf(VolumePointsFile,"%d\t%10.10f\t%10.10f\t%10.10f\n",count, point[0], point[1], point[2]);
//    }
//}


/* REALLY NICE TO HAVE THIS. DO NOT DELETE */
// Place to hold the gathered data. This is not really needed when we are simply counting the number of interior and surface points
// Allocate at root only
//    double *interior_points_array = NULL;
//    double *surface_points_array = NULL;
//    if (myRank == 0)
//    {
//        interior_points_array = new double[interior_disps[nProcs-1] + interior_counts[nProcs-1]];
//        surface_points_array = new double[surface_disps[nProcs-1] + surface_counts[nProcs-1]];
//    }
//
//    // Collect everything into the root
//    MPI_Gatherv(&interior_points_process[0], interior_nelements, MPI_DOUBLE,
//                interior_points_array, &interior_counts[0], &interior_disps[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);
//
//    MPI_Gatherv(&surface_points_process[0], surface_nelements, MPI_DOUBLE,
//                surface_points_array, &surface_counts[0], &surface_disps[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);
//
//    if(myRank == 0)
//    {
//        printf("Original interior points size on root = %zu\n", interior_points_process.size());
//        printf("Original surface points size on root = %zu\n", surface_points_process.size());
//        printf("Total interior and surface points = %d\t %d\n", interior_size, surface_size);
//        interior_points.resize(interior_size);
//        surface_points.resize(surface_size);
//    }
//
//    delete interior_points_array;
//    delete surface_points_array
