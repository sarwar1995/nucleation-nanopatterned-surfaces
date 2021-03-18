//
//  CubicLattice.hpp
//  Clusters
//
//  Created by Sarwar Hussain on 3/6/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#ifndef CubicLattice_hpp
#define CubicLattice_hpp

#include <stdio.h>
#include <algorithm>
#include <vector>

class CubicLattice
{
public:
    CubicLattice();
    ~CubicLattice();
    CubicLattice(int points, std::vector<std::vector<double> > box, std::vector<int> Seed);
    
    void CalcTranslationVector ();
    std::vector<int> GetDirectionalPoints();
    int GetTotalPoints();
    int GetTrueTotalPoints();
    
    int GetPointsAlongAxis (int axis);
    void RefineCellSize(int);
    void RepeatLattice (int dimension, int axis);
    std::vector<double> Generate1DLattice(int);
    void GenerateLattice();
    
    //Getters and Setters
    std::vector<double> get_translation_vector(){return translation;}
    std::vector<std::vector<double> > get_lattice_points(){return lattice_points;}
    
    // I/O commands
    void PrintLattice (FILE*);
    
    
protected:
    int n_points;
    int n_lattice_points;
    double unit_cell_size;
    std::vector<double> translation;
    std::vector<std::vector<double> > box;
    std::vector<int> seed;
    double box_x, box_y, box_z;
    std::vector<std::vector<double> > lattice_points;
    
};

#endif /* CubicLattice_hpp */
