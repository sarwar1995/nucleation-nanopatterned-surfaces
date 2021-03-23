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
//#include <mpi.h>

class MC {
public:
    MC();
    ~MC();
    MC(int, std::vector<std::vector<double> >&, int[3]);
    MC(int, std::vector<std::vector<double> >&, int[3], int);
    void generate_points();
    void generate_points_on_lattice();
    void add_points(std::vector<std::vector<double>> new_box, int direction_to_expand_in, int extra_points);
    void update_box (std::vector<std::vector<double>> new_box);
    inline double box_volume() {return ((box[0][1] - box[0][0]) * (box[1][1] - box[1][0]) * (box[2][1] - box[2][0]));}
    
    //Sets of Interior and surface points: Using sets to improve uniqueness checking. For the continuos update of volume and surface area by just checking the surface points. Only works for growth of the same type of cluster.
    std::set<std::vector<double> > interior_points;
    std::set<std::vector<double> > surface_points;
    
    
    //Volume and SA calculation
    std::vector<double> calc_volume_SA(Shape*, double); //Returns volume and surface area
    std::vector<double> update_volume_SA(Shape*, CellList*, double);
    
    //Miscellaneous I/O
    void print_points(FILE*);
    void print_surf_points(FILE*);
    void print_volume_points(FILE*);
    
    //Getters
    std::vector<std::vector<double> > get_points() {return x_points;}
    std::vector<std::vector<double> > get_box() {return box;}
    int get_num_points (){return n_points;};
    
    
private:
    // internal MPI variables.
//    int myRank, nProcs;
//    int points_chunk_per_procs;
    
    int seed[3] ;
    int n_points;
    int n_inside, n_near_surf;
    std::vector<std::vector<double> > box;
    std::vector<std::vector<double> > x_points;
    
    //Volume and SA given interior and surface points
    std::vector<double> getMeasures (int, int, double);

};

#endif /* MC_hpp */
