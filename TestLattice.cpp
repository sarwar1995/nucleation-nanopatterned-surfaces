//
//  TestLattice.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 3/8/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#include <stdio.h>
#include "CubicLattice.hpp"


int main(int argc, char* argv[])
{
    
    int points;
    std::vector<int> Seed(3);
    std::vector<std::vector<double> > the_box(3);
    double box_x, box_y, box_z;
    FILE* outputfile;
    printf("number of arguments = %d\n", argc);
    
    points  = atoi(argv[1]);
    Seed[0] = atoi(argv[2]);
    Seed[1] = atoi(argv[3]);
    Seed[2] = atoi(argv[4]);
    box_x   = atof(argv[5]);
    box_y   = atof(argv[6]);
    box_z   = atof(argv[7]);
    //outputfile = fopen(argv[8],"w");
    
    
    
    the_box[0].resize(2);
    the_box[1].resize(2);
    the_box[2].resize(2);
    
    the_box[0][0] = -box_x/2.0;     the_box[0][1] = box_x/2.0;
    the_box[1][0] = -box_y/2.0;     the_box[1][1] = box_y/2.0;
    the_box[2][0] = -box_z/2.0;     the_box[2][1] = box_z/2.0;
    
    printf("box_x= %10.5f\t %10.5f\n", the_box[0][0], the_box[0][1]);
    
    CubicLattice cubic_lattice (points, the_box, Seed);
    cubic_lattice.CalcTranslationVector(Seed);
    std::vector<double> translation_vector = cubic_lattice.get_translation_vector();
    printf("translation = [%10.5f\t %10.5f\t %10.5f\n]", translation_vector[0], translation_vector[1], translation_vector[2]);
    cubic_lattice.GenerateLattice();
    //cubic_lattice.PrintLattice(outputfile);
}
