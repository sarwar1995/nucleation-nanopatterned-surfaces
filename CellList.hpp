//
//  CellList.hpp
//  Clusters
//
//  Created by Sarwar Hussain on 8/5/20.
//  Copyright © 2020 Sarwar Hussain. All rights reserved.
//

#ifndef CellList_hpp
#define CellList_hpp

#include <stdio.h>
#include <vector>
#include <set>
#include "miscformulas.hpp"

//struct CellList{
//    std::vector<double> head;
//    std::vector<double> list;
//};

class CellList{
public:
    CellList();
    ~CellList();
    CellList(std::vector<std::vector<double> >, std::vector<double>, std::vector<std::vector<double> >);
    
    //Cell index and neighbors
    int calc_cell_ind (std::vector<double>&);
    std::vector<std::vector<double> > cell_points (int);
    int cell_size (int);
    std::vector<int> cell_ngbs (int);
    
    //Calculate and return the coordinates of neighbors of a list of points
    std::vector<std::vector<double> > calc_neighbors (std::vector<std::vector<double> >);
    std::vector<std::vector<double> > calc_neighbors (std::set<std::vector<double> >);
    
    //total cells
    int get_total_cells(){return total_cells;}
    
    //Printing
    void print_cell (int, int);
    
private:
    void generate();
    std::vector<double> head;   //List of indices of first atoms of each cell. Size=total_cells
    std::vector<double> list;   //List of atom indices
    std::vector<int> m_cells;   //number of cells in each direction
    int total_cells;            //Total number of cells = M^3
    int n_points;               //Total number of points in the simulation box
    std::vector<std::vector<double> > box;  //2d array of [start, end] of each box dimension
    std::vector<std::vector<double> > x_points; //coordinates of all points
    std::vector<double> Rc;                     //Size of each cell
    
    
};

#endif /* CellList_hpp */
