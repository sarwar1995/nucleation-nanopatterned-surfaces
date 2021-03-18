//
//  MC.hpp
//  Clusters
//
//  Created by Sarwar Hussain on 8/3/20.
//  Copyright © 2020 Sarwar Hussain. All rights reserved.
//

#ifndef MC_hpp
#define MC_hpp

#include <stdio.h>
#include <vector>
#include <set>
#include <random>
#include <math.h>
#include <iostream>
#include <cstdlib>
#include "CellList.hpp"
#include "Shape.hpp"
#include "miscformulas.hpp"
#include "CubicLattice.hpp"
#include <stdexcept>
#include <mpi.h>

class MC
{
public:
    MC(MPI_Comm);
    ~MC();
    MC(int, std::vector<std::vector<double> >&, int[3], MPI_Comm);
    MC(int, std::vector<std::vector<double> >&, int[3], int, MPI_Comm);
    void generate_points();
    void generate_points_on_lattice();
    inline double box_volume() {return BoxVolume;}
    
    //Sets of Interior and surface points: Using sets to improve uniqueness checking. For the continuos update of volume and surface area by just checking the surface points. Only works for growth of the same type of cluster.
    std::vector<std::vector<double> > interior_points;
    std::vector<std::vector<double> > surface_points;
    
    
    //Volume and SA calculation
    void addPoint(std::vector<double>&, std::vector<double>&);
    
    std::vector<int> loopPoints (Shape* cluster, double delta);
    std::vector<double> calc_volume_SA(Shape*, double); //Returns volume and surface area
    std::vector<double> update_volume_SA(Shape*, CellList*, double);
    
    //Miscellaneous I/O
    void print_points(FILE*);
    void print_surf_points(FILE*);
    void print_volume_points(FILE*);
    
    //Getters
    std::vector<std::vector<double> > get_points() {return x_points;}
    std::vector<std::vector<double> > get_box() {return box;}
    
    
private:
    // internal MPI variables.
    MPI_Comm branch_comm;
    int branch_rank, branch_size;
    int points_chunk_per_procs;
    
    std::vector<double> interior_points_process;
    std::vector<double> surface_points_process;
    
    int seed[3] ;
    int n_points;
    int n_inside, n_near_surf;
    std::vector<std::vector<double> > box;
    /* Check if double** is faster or not*/
    std::vector<std::vector<double> > x_points;
    double BoxVolume;
    
    //Volume and SA given interior and surface points
    std::vector<double> getMeasures (int, int, double);
    
};

#endif /* MC_hpp */
