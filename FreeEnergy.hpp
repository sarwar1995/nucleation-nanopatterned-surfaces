//
//  FreeEnergy.hpp
//  Clusters
//
//  Created by Sarwar Hussain on 2/7/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#ifndef FreeEnergy_hpp
#define FreeEnergy_hpp

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>

#define avogadro 6.022140857e23
#define kb 1.38064852e-23

// Each of these functions returns the free energy in (kb*T) and therefore need the volume and surface area inputs to be in SI units.
double free_energy (double, double, double, double , double , std::vector<double>& , double , double, double);
double free_energy_singlecap(double Rho, double Mu, double Sigma, double volume, double SA, double projected_SA, double theta, double T);
double free_energy_spherocylinder (double Rho, double Mu, double Sigma, double volume, double SA, std::vector<double>& projected_SA, double theta_good, double theta_bad, double T);
#endif /* FreeEnergy_hpp */
