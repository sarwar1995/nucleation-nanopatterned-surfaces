//
//  CubicLattice.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 3/6/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#include "CubicLattice.hpp"
#include <stdio.h>
#include <vector>
#include <set>
#include <random>

CubicLattice::CubicLattice()
{
    
}

CubicLattice::~CubicLattice()
{
    
}

CubicLattice::CubicLattice(int points, std::vector<std::vector<double> > the_box, std::vector<int> aSeed)
{
    n_points = points;
    box=the_box;
    box_x = box[0][1] - box[0][0];
    box_y = box[1][1] - box[1][0];
    box_z = box[2][1] - box[2][0];
    seed = aSeed;
    printf("box_x= %10.10f\t box_y = %10.10f\t box_z = %10.10f\n", box_x, box_y, box_z);
    double box_volume = box_x * box_y * box_z;
    unit_cell_size = cbrt(box_volume/(double)n_points);
    printf("unit_cell_size = %10.10f\n", unit_cell_size);
    int current_points = GetTotalPoints();
    printf("current_points = %d\n", current_points);
    RefineCellSize(current_points);
    printf("unit_cell_size after refinement = %10.10f\n", unit_cell_size);
    n_lattice_points = GetTotalPoints();
    printf("number of lattice points after refinement = %d\n", n_lattice_points);
}

std::vector<int> CubicLattice::GetDirectionalPoints()
{
    //These do not take into account the traslation vector which can change the number of points
    std::vector<int> directional_points (3,0);
    directional_points[0] = (int) (box_x/unit_cell_size);
    directional_points[1] = (int) (box_y/unit_cell_size);
    directional_points[2] = (int) (box_z/unit_cell_size);
    return directional_points;
}
int CubicLattice::GetTotalPoints()
{
    //These do not take into account the traslation vector which can change the number of points
    std::vector<int> directional_points = GetDirectionalPoints();
    return directional_points[0] * directional_points[1] * directional_points[2];
}


void CubicLattice::RefineCellSize(int current_points)
{
    double refined_cell_size;
    int point_differece = (current_points < n_points) ? n_points - current_points : 0;
    std::vector<double> box_length_ratio (3,0.0);
    box_length_ratio[0] = box_x/box_z;
    box_length_ratio[1] = box_y/box_z;
    box_length_ratio[2] = box_z/box_z;
    
    
    std::vector<int> additional_points (3, 0);
    additional_points[2] = (int)(cbrt(point_differece/(box_length_ratio[0]*box_length_ratio[1])));
    additional_points[0] = (int)(box_length_ratio[0] * (double)additional_points[2]) ;
    additional_points[1] = (int)(box_length_ratio[1] * (double)additional_points[2]) ;
    
    std::vector<int> directional_points = GetDirectionalPoints();
    int new_points_x = directional_points[0] + additional_points[0];
    int new_points_y = directional_points[1] + additional_points[1];
    int new_points_z = directional_points[2] + additional_points[2];
    
    double refined_size_x = box_x/(double)new_points_x ;
    double refined_size_y = box_y/(double)new_points_y ;
    double refined_size_z = box_z/(double)new_points_z ;
    if((refined_size_x == refined_size_y) && (refined_size_y==refined_size_z))
    {
        refined_cell_size = refined_size_x;
    }
    else
    {
        refined_cell_size = std::min({refined_size_x, refined_size_y, refined_size_z});
    }
    unit_cell_size = refined_cell_size;
}

void CubicLattice::CalcTranslationVector ()
{
    //Here box is a 3d vector of 2d vectors having min and max boundaries
    std::vector<double> translation_vector (3,0.0);
    std::mt19937 mt_engine_x (seed[0]); //deterministic seed. [0 2^64)
    std::mt19937 mt_engine_y (seed[1]);
    std::mt19937 mt_engine_z (seed[2]);
    std::uniform_real_distribution<double> dist_x(0, unit_cell_size); //both range numbers are inclusive
    std::uniform_real_distribution<double> dist_y(0, unit_cell_size);
    std::uniform_real_distribution<double> dist_z(0, unit_cell_size);
    translation_vector[0] = dist_x(mt_engine_x);
    translation_vector[1] = dist_y(mt_engine_y);
    translation_vector[2] = dist_z(mt_engine_z);
    translation = translation_vector ;
}


int CubicLattice::GetPointsAlongAxis (int axis)
{
    //These take into account the translation vector and so reflect the true number of points along an axis
    std::vector<int> directional_points = GetDirectionalPoints();
    int n_points_along_axis;
    if (translation[axis] == 0 || translation[axis] == unit_cell_size)
    {
        n_points_along_axis = directional_points[axis] + 1;
    }
    else
    {
        n_points_along_axis = directional_points[axis];
    }
    return n_points_along_axis;
}

int CubicLattice::GetTrueTotalPoints()
{
    int x_axis_points = GetPointsAlongAxis(0);
    int y_axis_points = GetPointsAlongAxis(1);
    int z_axis_points = GetPointsAlongAxis(2);
    return (x_axis_points * y_axis_points * z_axis_points);
}

std::vector<double> CubicLattice::Generate1DLattice(int axis)
{
    int n_points_along_axis = GetPointsAlongAxis(axis);
    std::vector<double>points_along_axis(n_points_along_axis, 0.0);
    double starting_axis_point = box[axis][0];
    for(int i=0; i<n_points_along_axis; i++)
    {
        double translation_axis = (translation[axis] == 0 || translation[axis] == unit_cell_size) ? 0 : translation[axis] ;
        points_along_axis[i] = (i==0) ? (starting_axis_point + translation_axis) : (points_along_axis[i-1] + unit_cell_size);
        if(points_along_axis[i] > box[axis][1])
        {
            printf("point[%d] went outside box bound for translation ( %10.10f ) before total points along axis = %d\n", i, translation_axis, axis);
            abort();
        }
    }
    return points_along_axis;
}


void CubicLattice::GenerateLattice()
{
    std::vector<double> lattice_x = Generate1DLattice(0);
    printf("1D lattice generated size = %d\n", (int)lattice_x.size());
    std::vector<double> lattice_y = Generate1DLattice(1);
    std::vector<double> lattice_z = Generate1DLattice(2);
    
    n_lattice_points = (int)lattice_x.size() * (int)lattice_y.size() * (int)lattice_z.size() ;
    printf("n_points = %d\n", n_lattice_points);
    lattice_points.resize(n_lattice_points);
    
    int index;
    for(int i=0; i < lattice_x.size(); i++)
    {
        for(int j=0; j < lattice_y.size(); j++)
        {
            for(int k=0; k < lattice_z.size(); k++)
            {
                index = k + j * (int)lattice_z.size() + i * (int)lattice_y.size() * (int)lattice_z.size() ;
                if(index >= n_lattice_points)
                {
                    printf("index greater = %d\n", index);
                }
                lattice_points[index].push_back(lattice_x[i]);
                lattice_points[index].push_back(lattice_y[j]);
                lattice_points[index].push_back(lattice_z[k]);
            }
            
        }
    }
    printf("lattice generation done\n");
    
}

void CubicLattice::PrintLattice (FILE* outputfile)
{
    if(outputfile != NULL)
    {
        for(int i=0; i<(int)lattice_points.size(); i++)
        {
            fprintf(outputfile,"%10.10f\t%10.10f\t%10.10f\n" ,lattice_points[i][0], lattice_points[i][1], lattice_points[i][2]);
        }
        fclose(outputfile);
    }
    else
    {
        printf("File failed to open while printing cubic lattice \n");
        abort();
    }
    
}
