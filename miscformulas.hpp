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
std::vector<double> add_vectors (std::vector<double>&, std::vector<double>&);
double vector_norm (std::vector<double>&);
void add_to_N (double N, double G, double Rg, double Rb, double Rg_secondary, int db, int dg_secondary, double volume, double SA, double dN, double Nmin, int lenN, std::vector<double>& NArray_Gmin, std::vector<double>& NArray_confs, std::vector<std::vector<double> >& NArray_quant);
void print_NGDataFile (std::vector<std::vector<double> >&, FILE*);

void addQuant(std::vector<double> &destination, std::vector<double> &origin, int size_origin);

template <class T>
std::vector<T> scalar_mult_to_vector (T scalar, std::vector<T>& a)
{
    std::vector<double> result(a.size(),0.0);
    for(size_t i=0; i<a.size(); i++)
    {
        result[i] = scalar * a[i];
    }
    return result;
    
}

template <class T>
std::vector<T> add_vectors (std::vector<T>& a, std::vector<T>& b)
{
    std::vector<double> result(a.size(),0.0);
    if(a.size() != b.size())
    {
        throw std::logic_error ("vector size should be same for addition\n");
    }
    else
    {
        
        for(size_t i=0; i<a.size(); i++)
        {
            result[i] = a[i] + b[i];
        }
    }
    return result;
}

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

void print_points_tofile(std::vector<std::vector<double> >, FILE*);

#endif /* miscformulas_hpp */
