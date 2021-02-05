//
//  TestHourglass.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 2/4/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#include "Cylinder.hpp"

Cylinder::Cylinder()
{
    
}

Cylinder::~Cylinder()
{
    
}

Cylinder::Cylinder(double radius, std::vector<std::vector<double> >& circle_centres)
{
    rc = radius;
    if((int)circle_centres.size() != 2)
    {
        throw std::invalid_argument("Size of centres should be 2");
        abort();
    }
    else
    {centres = circle_centres;}
    
}

int Cylinder::isWithin(std::vector<double> & point)
{
    int yes = 0;
    std::vector<double> pointMinusP1 (3, 0.0);
    std::vector<double> pointMinusP2 (3, 0.0);
    std::vector<double> Axis (3, 0.0);
    for (size_t i=0; i<point.size(); i++)
    {
        pointMinusP1[i] = point[i] - centres[0][i];
        pointMinusP2[i] = point[i] - centres[1][i];
        Axis[i] = centres[1][i] - centres[0][i];
    }
    double dotProd_with_p1 = dotprod(pointMinusP1, Axis);
    double dotProd_with_p2 = dotprod(pointMinusP2, Axis);
    if(dotProd_with_p1 >= 0 && dotProd_with_p2 <= 0)
    {yes = 1;}
    else{yes = 0;}
    return yes;
}


int Cylinder::nearSurface(std::vector<double> & point,  double delta)
{
    std::vector<double> centre_near_wall = centres[0]; //p1
    std::vector<double> centre_near_free = centres[1]; //p2
    double point_axis_distance = point_line_dist(point, centre_near_wall, centre_near_free);
    int pointWithin = isWithin(point);
    if((pointWithin == 1) && (point_axis_distance <= rc + delta) && (rc - delta <= point_axis_distance))
    {
        return 1;
    }
    else{return 0;}    
}

int Cylinder::isInside(std::vector<double> & point)
{
    int result = 0;
    std::vector<double> centre_near_wall = centres[0]; //p1
    std::vector<double> centre_near_free = centres[1]; //p2
    
    double point_axis_distance = point_line_dist(point, centre_near_wall, centre_near_free);
    int pointWithin = isWithin(point);
    if(pointWithin == 1 && point_axis_distance < rc )
    {
        result = 1;
    }
    else{result = 0;}
    return result ;
    
}
