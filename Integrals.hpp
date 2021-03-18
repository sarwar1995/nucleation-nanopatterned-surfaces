//
//  Integrals.hpp
//  Clusters
//
//  Created by Sarwar Hussain on 7/6/20.
//  Copyright © 2020 Sarwar Hussain. All rights reserved.
//

#ifndef Integrals_hpp
#define Integrals_hpp

#include <stdio.h>
#include <math.h>

#endif /* Integrals_hpp */

double Volume_integral (double R, double z, double theta_1, double theta_2, double phi_1, double phi_2)
{
    double result;
    double R_cube = R*R*R ;
    double z_cube = z*z*z ;
    double term_1 = R_cube*(-cos(phi_2) + cos(phi_1));
    double sec_square_1 = (1/cos(phi_1)) * (1/cos(phi_1));
    double sec_square_2 = (1/cos(phi_2)) * (1/cos(phi_2));
    double term_2 = (-z_cube/2)*(sec_square_2 - sec_square_1);
    result = (1/3) * (theta_2 - theta_1) * (term_1 + term_2);
    return(result) ;
}

double Surface_integral (double R, double theta_1, double theta_2, double phi_1, double phi_2)
{
    double result;
    double R_square = R*R ;
    result = R_square * (theta_2 - theta_1) * -1*(cos(phi_2) - cos(phi_1)) ;
    return(result);
}
