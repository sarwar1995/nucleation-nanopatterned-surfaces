//
//  FreeEnergy.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 2/7/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#include "FreeEnergy.hpp"

double free_energy (double volume, double SA, std::vector<double>& projected_SA, double theta_good, double theta_bad)
{
    /* Here projected_SA has the projected surface area lying exclusively on the [0]: good, [1]:bad on left and [2]:bad on right patches */
    double G;
    double Vcomp = -1*volume* Rho * Mu;
    double Scomp = Sigma*SA;
    double SA_proj_good = projected_SA[0];
    double SA_proj_bad_1 = projected_SA[1];
    double SA_proj_bad_2 = projected_SA[2];
    double Sproj_comp = -Sigma*(cos(theta_good)*SA_proj_good + cos(theta_bad)*SA_proj_bad_1 + cos(theta_bad)*SA_proj_bad_2);
    G = Vcomp + Scomp + Sproj_comp;
    return G;
    
}
