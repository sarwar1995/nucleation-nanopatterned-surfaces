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
    MC();
    MC(MPI_Comm);
    ~MC();
    MC(int, std::vector<std::vector<double> >&, int[3], MPI_Comm, double);
    MC(int, std::vector<std::vector<double> >&, int[3], int, MPI_Comm, double);
    void generate_points();
    void generate_points_on_lattice();
    void add_points(std::vector<std::vector<double>> new_box, int direction_to_expand_in, int extra_points);
    
    void update_box (std::vector<std::vector<double>> new_box);
    void update_virtual_box (std::vector<std::vector<double>> new_box);
    
    void update_points_chunk_per_procs();
    inline double box_volume() {return ((box[0][1] - box[0][0]) * (box[1][1] - box[1][0]) * (box[2][1] - box[2][0]));}
    inline double virtual_box_volume() {return ((virtual_box[0][1] - virtual_box[0][0]) * (virtual_box[1][1] - virtual_box[1][0]) * (virtual_box[2][1] - virtual_box[2][0]));}
    
    //Sets of Interior and surface points: Using sets to improve uniqueness checking. For the continuos update of volume and surface area by just checking the surface points. Only works for growth of the same type of cluster.
    std::set<std::vector<double> > interior_points;
    std::set<std::vector<double> > surface_points;
    
    
    //Volume and SA calculation
    void addPoint(std::vector<double>&, std::vector<double>&);
    
    std::vector<int> loopPoints (Shape* cluster);
    std::vector<double> calc_volume_SA(Shape*); //Returns volume and surface area
    std::vector<double> update_volume_SA(Shape*, CellList*, double);
    std::vector<double> calc_volume_SA_with_virtual_points(Shape* cluster);
    
    //Combined addition of points and looping at the same time
    void add_and_loop_points(Shape* cluster, std::vector<std::vector<double>> new_box, int direction_to_expand_in, int extra_points);
    
    
    //Miscellaneous I/O
    void print_points(FILE*);
    void print_surf_points(Shape* cluster, FILE*);
    void print_volume_points(FILE*);
    
    //Getters
    std::vector<std::vector<double> > get_points() {return x_points;}
    std::vector<std::vector<double> > get_box() {return box;}
    int get_num_points (){return n_points;};
    
    //Setters related to setting branch communicator
    void set_branch_comm (MPI_Comm);
    
    int get_n_inside_virtual(){return n_inside_virtual;}
    int get_n_near_surf_virtual() {return n_near_surf_virtual;}
    
private:
    // internal MPI variables.
     int myRank, nProcs;
    MPI_Comm branch_comm;
    int branch_rank, branch_size;
    int points_chunk_per_procs;
    
    std::vector<double> interior_points_process;
    std::vector<double> surface_points_process;
    
    int seed[3] ;
    int n_points;
    int n_virtual_points; //For direct calculation of adding and checking points without storing them, due to space limitations.
    int n_inside, n_near_surf;
    int n_inside_virtual, n_near_surf_virtual;
    std::vector<std::vector<double> > box;
    std::vector<std::vector<double> > virtual_box; //virtual_box encloses box and only box contains concrete points whereas virtual_box does not contain any points
    double delta;
    /* Check if double** is faster or not*/
    std::vector<std::vector<double> > x_points;
    
    //Volume and SA given interior and surface points
    std::vector<double> getMeasures (int, int);
    std::vector<double> getMeasures_virtual(int, int);
    
    void check_point(Shape* cluster, std::vector<double>& point, int& interior_count, int& surface_count);
    std::vector<int> gather_points(int& , int&);
    void reset_virtual();
    bool within_box(std::vector<double>& point);
};

#endif /* MC_hpp */
