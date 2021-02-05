//
//  miscformulas.hpp
//  Clusters
//
//  Created by Sarwar Hussain on 7/1/20.
//  Copyright © 2020 Sarwar Hussain. All rights reserved.
//
#include <stdio.h>
#include <vector>
#include <math.h>

#ifndef miscformulas_hpp
#define miscformulas_hpp

double potency_factor(double);
double round_nearN (double, double);
double point_plane_dist(std::vector<double>&, std::vector<double>&, double);
double MagCrossProduct(double[3], double[3]);
double point_line_dist (double[3], std::vector<double>&, std::vector<double>&);
double point_line_dist (std::vector<double>&, std::vector<double>&, std::vector<double>&);

double dotprod (std::vector<double>& , std::vector<double>&);
void subtract_vectors(std::vector<double>&, std::vector<double>&, std::vector<double>&);
double vector_norm (std::vector<double>&);
void add_to_N (double, double, double, double, double, double, int, std::vector<double>& , std::vector<double>& , std::vector<std::vector<double> >& );


template <class T>
int InVector(std::vector<T> vec, T args)
{
    if(std::find(vec.begin(), vec.end(), args) != vec.end()) //If found before end then inVector
    {
        return 1;
    }
    else
    {return 0;}
}

template <class T>
void appendTovec (std::vector<T>& vec1, std::vector<T>& vec2)
{
    vec1.insert(vec1.end(),vec2.begin(),vec2.end());
}

void print_point(std::vector<double>);

#endif /* miscformulas_hpp */
