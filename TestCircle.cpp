//
//  TestCircle.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 2/11/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#include <iostream>
#include <stdio.h>
#include <vector>
#include <string.h>
#include "Circle.hpp"

using namespace std;


int main(int argc, const char * argv[])
{
    double R1, R2;
    double C1x, C2x;
    
    R1  = atof(argv[1]);
    R2  = atof(argv[2]);
    C1x = atof(argv[3]);
    C2x = atof(argv[4]);
    
    std::vector<double> centre_1 (3,0.0);
    centre_1[0] = C1x;
    std::vector<double> centre_2 (3,0.0);
    centre_2[0] = C2x;
    
    Circle circle_1 (centre_1 , R1);
    Circle circle_2 (centre_2 , R2);
    
    double intersect_area_from_1 = circle_1.intersection_area(circle_2);
    double intersect_area_from_2 = circle_2.intersection_area(circle_1);
    
    printf("area from 1 = %10.10f\n", intersect_area_from_1);
    printf("area from 2 = %10.10f\n", intersect_area_from_2);
    
    printf("area of self 1= %10.10f\n", circle_1.area());
    printf("area of self 2= %10.10f\n", circle_2.area());
    
    return 0;
}
