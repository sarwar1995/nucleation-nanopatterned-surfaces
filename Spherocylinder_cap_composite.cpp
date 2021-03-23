//
//  Spherocylinder_cap_composite.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 3/18/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#include "Spherocylinder_cap_composite.hpp"


Spherocylinder_cap_composite::Spherocylinder_cap_composite()
{
    
}


Spherocylinder_cap_composite::~Spherocylinder_cap_composite()
{
    
}

Spherocylinder_cap_composite::Spherocylinder_cap_composite(std::vector<Spherical_cap*> caps, SpheroCylinder* spherocylinder):
list_of_spherical_caps(caps),
sphero_cylinder(spherocylinder)
{
    Shape* sphero_cylinder_ptr = (Shape*) spherocylinder;
    list_of_shapes.push_back(sphero_cylinder_ptr);
    
    for (size_t i=0; i<caps.size(); i++)
    {
        Shape* cap_i_ptr = (Shape*) caps[i] ;
        list_of_shapes.push_back(cap_i_ptr);
    }
    
}

double Spherocylinder_cap_composite::projected_SA()
{
    return -1;
}



int Spherocylinder_cap_composite::isInside(std::vector<double>& point)
{
    //Write a test to check this i.e. if the point is inside any of the caps then it is inside the composite cluster
    for (size_t i=0; i<list_of_shapes.size(); i++)
    {
        if (list_of_shapes[i]->isInside(point) == 1)
        {
            return 1;
        }
        else{continue;}
    }
    return 0;
}

int Spherocylinder_cap_composite::nearSurface(std::vector<double> & point, double delta)
{
    for (size_t i=0; i<list_of_shapes.size(); i++)
    {
        if (list_of_shapes[i]->nearSurface(point, delta) == 1)
        {
            int outOfOthers = 1;
            for (size_t j=0; j<list_of_shapes.size(); j++)
            {
                if(i != j)
                {
                    if (list_of_shapes[j]->isInside(point) == 1)
                    {outOfOthers = 0;}
                }
            }
            if(outOfOthers==1)
            {return 1;}
            else{return 0;} //continue;
        }
    }
    return 0;  //If it gets out of the for loop without returning 1 within it then returns 0
}

std::vector<double> Spherocylinder_cap_composite::xy_spread()
{
    std::vector<double> thisSpread (4,0.0); //[Xmin, Xmax, Ymin, Ymax]
    thisSpread = list_of_shapes[0]->xy_spread();
    
    for(size_t i = 1; i<list_of_shapes.size(); i++)
    {
        std::vector<double> spread_ith_cap = list_of_shapes[i]->xy_spread();
        if(spread_ith_cap[0] <= thisSpread[0]){thisSpread[0] = spread_ith_cap[0];}
        if(spread_ith_cap[1] >= thisSpread[1]){thisSpread[1] = spread_ith_cap[1];}
        if(spread_ith_cap[2] <= thisSpread[2]){thisSpread[2] = spread_ith_cap[2];}
        if(spread_ith_cap[3] >= thisSpread[3]){thisSpread[3] = spread_ith_cap[3];}
    }
    return thisSpread;
}

std::vector<double> Spherocylinder_cap_composite::threeDim_spread()
{
    //To check the breach of 5 box planes, excluding the wall surface
    //[Xmin, Xmax, Ymin, Ymax, Zmax]
    //Start with central caps spread as the spread of cluster and then modify based on other caps
    std::vector<double> thisSpread = list_of_shapes[0]->threeDim_spread();
//    printf("Inside Spherocylinder_cap_composite three dim spread \n");
//    printf("list_of_shapes size = %d\n", (int)list_of_shapes.size());
    for(size_t i = 1; i<list_of_shapes.size(); i++)
    {
        std::vector<double> iCap3DSpread = list_of_shapes[i]->threeDim_spread();
        if(iCap3DSpread[0] <= thisSpread[0]) { thisSpread[0] = iCap3DSpread[0]; }
        if(iCap3DSpread[1] >= thisSpread[1]) { thisSpread[1] = iCap3DSpread[1]; }
        if(iCap3DSpread[2] <= thisSpread[2]) { thisSpread[2] = iCap3DSpread[2]; }
        if(iCap3DSpread[3] >= thisSpread[2]) { thisSpread[3] = iCap3DSpread[3]; }
        if(iCap3DSpread[4] >= thisSpread[4]) { thisSpread[4] = iCap3DSpread[4]; }
    }
    return thisSpread;
}
