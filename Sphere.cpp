//
//  Sphere.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 7/1/20.
//  Copyright © 2020 Sarwar Hussain. All rights reserved.
//

#include "Sphere.hpp"

Sphere::Sphere(){
    radius = 1.0;
    centre[0]= 0.0;
    centre[1]= 0.0;
    centre[2]= 0.0;
    calc_volume();
    calc_SA();
}

Sphere::Sphere(double r, std::vector<double>& c)
{
    radius = r;
    centre[0] = c[0];
    centre[1] = c[1];
    centre[2] = c[2];
    calc_volume();
    calc_SA();
}

Sphere::~Sphere(){
    radius=0;
}

void Sphere::calc_volume()
{
    volume = (4.0/3.0)*pi*radius*radius*radius ;
}

void Sphere::calc_SA()
{
    SA = 4.0*pi*radius*radius ;
}

//Points on surface are not counted as being inside
int Sphere::isInside (std::vector<double>& point){
    
    if (centre_point_dist(point) < radius)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

int Sphere::nearSurface(std::vector<double>& point, double delta){
    double dist = centre_point_dist (point);
    
    if((dist <= radius + delta) && (radius - delta <= dist))
    {
        return 1;
    }
    else{return 0;}
    
}


double Sphere::centre_point_dist (std::vector<double>& point)
{
    double dist = sqrt((centre[0] - point[0]) * (centre[0] - point[0]) + (centre[1] - point[1]) * (centre[1] - point[1]) + (centre[2] - point[2])*(centre[2] - point[2])) ;
    return dist;
}

//void Sphere::getCentre(std::vector<double>& point)
//{
//    point[0] = centre[0];     point[1] = centre[1];     point[2] = centre[2];
//}
//
//void Sphere::getRadius(double & r)
//{
//    r = radius;
//}
