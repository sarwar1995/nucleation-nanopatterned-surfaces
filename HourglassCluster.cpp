//
//  HourglassCluster.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 2/4/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#include "HourglassCluster.hpp"

Hourglass::Hourglass()
{
    
}

Hourglass::~Hourglass()
{
    
}

Hourglass::Hourglass(std::vector<Shape*>& list)
{
    Components = list;
}

int Hourglass::isInside(std::vector<double> & point)
{
    for(size_t i=0; i<Components.size(); i++)
    {
        if(Components[i]->isInside(point) == 1)
        {
            return 1;
        }
    }
    return 0;
    
}

int Hourglass::nearSurface(std::vector<double> & point, double delta)
{
    for (size_t i=0; i<Components.size(); i++)
    {
        if (Components[i]->nearSurface(point, delta) == 1)
        {
            int outOfOthers = 1;
            for (size_t j=0; j<Components.size(); j++)
            {
                if(i != j)
                {
                    if (Components[j]->isInside(point) == 1){outOfOthers = 0;}
                    else{continue;}
                }
                else{continue;}
            }
            if(outOfOthers==1)
            {return 1;}
            else{return 0;} //continue;
        }
        else{continue;}
    }
    return 0;  //If it gets out of the for loop without returning 1 within it then returns 0
}

