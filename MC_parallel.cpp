//
//  MC_parallel.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 2/14/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#include "MC_parallel.hpp"

MC::MC(MPI_Comm branch_comm_input)
{
    branch_comm = branch_comm_input;
    MPI_Comm_rank (branch_comm, &branch_rank);
    MPI_Comm_size (branch_comm, &branch_size);
    
    box.resize(3);
    for(size_t i=0; i<3; i++){
        box[i].resize(2);
        for(size_t j=0; j<2; j++){
            if(j==0){box[i][j] = 0;}
            else{box[i][j] = 1;}
        }
    }
    n_points = 1e6;
    x_points.resize(n_points);
    points_chunk_per_procs = (int)(n_points/branch_size) ;
    if(n_points % branch_size != 0){printf("Warning: number of points=%d are not perfectly divisible by number of processes=%d\n",n_points, branch_size); abort();}
    
}

MC::~MC()
{
    x_points.clear();
    box.clear();
}

//MCbox has the start and end boundaries of the box in all three dimensions
MC::MC(int num_points, std::vector<std::vector<double> >& MCbox, int aSeed[3], MPI_Comm branch_comm_input){
    
    branch_comm = branch_comm_input;
    MPI_Comm_rank (branch_comm, &branch_rank);
    MPI_Comm_size (branch_comm, &branch_size);
    
    n_points = num_points;
    box = MCbox;
    BoxVolume = (box[0][1] - box[0][0]) * (box[1][1] - box[1][0]) * (box[2][1] - box[2][0]) ;
    x_points.resize(n_points);
    seed[0] = aSeed[0];
    seed[1] = aSeed[1];
    seed[2] = aSeed[2];
    points_chunk_per_procs = (int)(n_points/branch_size) ;
    printf("branch_rank = %d\t branch_size = %d\n", branch_rank, branch_size);
    printf("points_chunk_per_procs = %d\n", points_chunk_per_procs);
      if(n_points % branch_size != 0){printf("Warning: number of points=%d are not perfectly divisible by number of processes=%d\n",n_points, branch_size); abort();}
    generate_points();
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


std::vector<double> MC::getMeasures (int n_inside, int n_near_surf, double delta)
{
    std::vector<double> measures(2,0.0);
    //printf("n_inside = %d\t n_near_surf=%d\n",n_inside, n_near_surf);
    double fraction_volume =  ((double)n_inside/(double)n_points);
    //printf("fraction_volume = %10.10f\n",fraction_volume);
    double fraction_SA =  ((double)n_near_surf/(double)n_points);
    double volume = fraction_volume * BoxVolume ;
    double volume_near_surf = fraction_SA * BoxVolume ;
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

std::vector<int> MC::loopPoints (Shape* cluster, double delta)
{
    
    
    std::vector<int> int_surf_size(2,-1);
    std::vector <double> point (3, 0.0);

    int count=0;
    int interior_count = 0;
    int surface_count = 0;
    int start = branch_rank*points_chunk_per_procs;
    int end = (branch_rank+1)*points_chunk_per_procs - 1;
    
    
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
    
    int interior_counts[branch_size];
    int surface_counts[branch_size];

//    //Remember that interior_points_process and surface...process are unrolled linear contiguous vectors
//    int interior_nelements = (int)interior_points_process.size();
//    int surface_nelements = (int)surface_points_process.size();
    
    MPI_Barrier(branch_comm);
    // Each process tells the root how many elements it holds
    MPI_Gather(&interior_count, 1, MPI_INT, &interior_counts[0], 1, MPI_INT, 0, branch_comm);
    MPI_Gather(&surface_count, 1, MPI_INT, &surface_counts[0], 1, MPI_INT, 0, branch_comm);
    
//    int interior_disps[branch_size];
//    int surface_disps[branch_size];
//
//    // Displacement for the first chunk of data - 0
//    for (int i = 0; i < branch_size; i++)
//    {
//        interior_disps[i] = (i > 0) ? (interior_disps[i-1] +  interior_counts[i-1]) : 0;
//        surface_disps[i] = (i > 0) ? (surface_disps[i-1] +  surface_counts[i-1]) : 0;
//
//    }
    
    if(branch_rank == 0)
    {
        int_surf_size[0] = 0;
        int_surf_size[1] = 0;
        for(int w=0; w<branch_size; w++)
        {
            int_surf_size[0] = int_surf_size[0] + interior_counts[w];
            int_surf_size[1] = int_surf_size[1] + surface_counts[w];
        }
        
//        int interior_size = interior_disps[branch_size-1] + interior_counts[branch_size-1];
//        int surface_size = surface_disps[branch_size-1] + surface_counts[branch_size-1] ;
//        if (interior_size%3 != 0 || surface_size%3 != 0 || interior_size < 0 || surface_size < 0)
//        {
//            printf("Interior size = %d and surface size = %d are not divisible by 3\n", interior_size, surface_size);
//            printf("branch_rank = %d\n", branch_rank);
//            abort();
//        }
//        int_surf_size[0] = (interior_size/3);
//        int_surf_size[1] = (surface_size/3);
        if(int_surf_size[0] >= n_points || int_surf_size[1] >= n_points)
        {
            printf("num points are larger than total points interior count=%d\t surface count=%d\n", int_surf_size[0], int_surf_size[1]);
            abort();
        }
    }
    
    return int_surf_size;
    
   
}

//This function is to calculate initial volume and SA
//Add Exceptions here try and catch block in main to catch this exception
std::vector<double> MC::calc_volume_SA(Shape* cluster, double delta)
{
    
    //Calculating volume and SA for the first time. Currently doing this whenever Rg is increased.
    std::vector<double> measures(2,-1.0);
    n_inside = 0;
    n_near_surf = 0;
    std::vector<int> sizes = loopPoints(cluster, delta);
    
    if(branch_rank == 0)
    {
        n_inside = sizes[0];
        n_near_surf = sizes[1];
        if(n_inside==0 || n_near_surf==0)
        {
            throw std::logic_error("n_inside or n_near_surf are zero");
        }
        
        measures = getMeasures (n_inside, n_near_surf, delta);
    }
    return measures;
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

//void MC::print_surf_points(FILE * SurfpointsFile)
//{
//    std::set<std::vector<double> >::iterator it;
//    int count=0;
//    
//    for(it=surface_points.begin(); it!=surface_points.end(); it++)
//    {
//        count++;
//        std::vector<double> point = *it ;
//        fprintf(SurfpointsFile,"%d\t%10.10f\t%10.10f\t%10.10f\n",count, point[0], point[1], point[2]);
//    }
//}


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
