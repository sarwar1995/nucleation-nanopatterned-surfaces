//
//  CompositeShape.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 3/21/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#include "CompositeShape.hpp"
CompositeShape::CompositeShape(){
    
}
CompositeShape::~CompositeShape(){
    
}

CompositeShape::CompositeShape(std::vector<Shape*> shapes): list_of_shapes(shapes)
{
    
}

int CompositeShape::isInside(std::vector<double>& point)
{
    //Write a test to check this i.e. if the point is inside any of the caps then it is inside the composite cluster
    printf("Inside the isInside function of composite Shape\n");
    printf("list_of_shapes size = %d\n", (int)list_of_shapes.size());
    for (size_t i=0; i<list_of_shapes.size(); i++)
    {
        printf("list of shape i=%d\n", (int)i);
        
        if (list_of_shapes[i]->isInside(point) == 1)
        {
            return 1;
        }
        else{continue;}
    }
    return 0;
}

int CompositeShape::nearSurface(std::vector<double> & point, double delta)
{
    printf("Inside the nearSurface function of composite Shape\n");
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
    printf("After the nearSurface function of composite Shape\n");
    return 0;  //If it gets out of the for loop without returning 1 within it then returns 0
}

std::vector<double> CompositeShape::xy_spread()
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

std::vector<double> CompositeShape::threeDim_spread()
{
    //To check the breach of 5 box planes, excluding the wall surface
    //[Xmin, Xmax, Ymin, Ymax, Zmax]
    //Start with central caps spread as the spread of cluster and then modify based on other caps
    std::vector<double> thisSpread = list_of_shapes[0]->threeDim_spread();
    
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

