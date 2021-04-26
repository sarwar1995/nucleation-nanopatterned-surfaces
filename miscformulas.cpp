//
//  miscformulas.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 7/1/20.
//  Copyright © 2020 Sarwar Hussain. All rights reserved.
//

#include "miscformulas.hpp"

double potency_factor(double theta)
{
    double oneMinusCosSquare = (1-cos(theta))*(1-cos(theta));
    return ((oneMinusCosSquare)*(2 + cos(theta)))/4.0 ;
}


double round_nearN (double x, double dN)
{
    double temp = x/dN ;
    double round_t = round(temp);
    double result = dN*round_t;
    return(result);
}

//std::vector<double> three_plane_intersection(std::vector<double>& normal_1, std::vector<double>& normal_2,  double constant_1, double constant_2) //std::vector<double>& normal_3, double constant_3
//{
//    //Third plane is y=0;
//    //    double B1, B2, B3;
//    //    double C1, C2, C3;
//    //    double D1, D2, D3;
//    double A1 = normal_1[0], B1 = normal_1[1], C1 = normal_1[2], D1 = constant_1;
//    double A2 = normal_2[0], B2 = normal_2[1], C2 = normal_2[2], D2 = constant_2;
//    double x = (C1*D2 - C2*D1)/(C2*A1 - A2*C1) ;
//    double z = (A2*D1 - A1*D2)/(C2*A1 - A2*C1) ;
//    std::vector<double> result(3,0.0);
//    result[0] = x; result[1] = 0; result[2] = z;
//    return(result);
//
//}

double point_plane_dist(std::vector<double>& normal, std::vector<double>& point, double constant){
    double A, B, C;
    double X, Y, Z;
    double D = constant;
    A = normal[0]; B = normal[1]; C = normal[2];
    X = point[0]; Y = point[1]; Z = point[2];
    double numerator =  A*X + B*Y + C*Z + D ;
    //numerator = (numerator >= 0) ? numerator : -1.0*numerator ;
    double denom = sqrt(A*A + B*B + C*C);
    double result = (double) (numerator/denom);
    return (result) ;
}


double MagCrossProduct(double V1[3], double V2[3])
{
    double cross[3];
    double v1x, v1y, v1z, v2x, v2y, v2z;
    v1x = V1[0]; v1y = V1[1]; v1z = V1[2];
    v2x = V2[0]; v2y = V2[1]; v2z = V2[2];
    cross[0] = v1y*v2z - v2y*v1z ;
    cross[1] = -1*(v1x*v2z - v2x*v1z);
    cross[2] = v1x*v2y - v2x*v1y ;
    double result = sqrt(cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2]);
    return result;
}


double point_line_dist (std::vector<double>& X0, std::vector<double> & X1, std::vector<double> & X2)
{
    double X0_array[3];
    for(size_t i=0; i<X0.size(); i++)
    {
        X0_array[i] = X0[i];
    }
    double result = point_line_dist (X0_array, X1, X2);
    return result;
}

double point_line_dist (double X0[3], std::vector<double> & X1, std::vector<double> & X2)
{
    //Line between X1= (x1,y1,z1) and X2 = (x1,y1,z1)
    double dist;
    double V1[3];
    double V2[3];
    double V3[3];
    
    for (int i =0; i<3; i++)
    {
        V1[i] = X0[i]-X1[i];
        V2[i] = X0[i]-X2[i];
        V3[i] = X2[i]-X1[i];
    }
    double num = MagCrossProduct(V1,V2) ;
    double denom = sqrt(V3[0]*V3[0] + V3[1]*V3[1] + V3[2]*V3[2]);
    dist = num/denom ;
    return dist;
}


double dotprod (std::vector<double>& V1 , std::vector<double>& V2)
{
    double result = 0.0;
    for (int i =0 ; i<3; i++)
    {
        result = result + V1[i]*V2[i];
    }
    return(result);
    
}

void subtract_vectors(std::vector<double>& a, std::vector<double>& b, std::vector<double>& c)
{
    if(a.size() == b.size())
    {
        c.resize(a.size());
        for(size_t i=0; i<a.size(); i++)
        {
            c[i] = b[i] - a[i];
        }
    }
    else{printf("Different sized vectors cannot be subtracted"); abort();}
    
}



double vector_norm (std::vector<double>& a)
{
    double norm = 0.0;
    for(size_t i=0; i<a.size(); i++)
    {
        norm = norm + (a[i]*a[i]);
    }
    return (sqrt(norm));
}



//double radius_of_intersection(double r1, double r2, double d) //d is the centre of r2. Of two spheres
//{
//    double dsquare = d*d ;
//    double r1square = r1*r1;
//    double r2square = r2*r2;
//    double dr1r2square = (dsquare - r2square + r1square)*(dsquare - r2square + r1square) ;
//    double inside = 4*dsquare*r1square - dr1r2square;
//    double twod = 2*d;
//    double result = sqrt(inside)/twod ;
//    return(result);
//
//}


// Coordinate of the plane of intersection of the two spheres if origin is at the centre of the first sphere
//double x_intersection (double r1, double r2, double d) //r2 second sphere
//{
//    double dsquare = d*d;
//    double result = (dsquare - r2*r2 + r1*r1)/(2*d);
//    return (result);
//}
//
//
//double circle_circle_int_radii(double R, double r, std::vector<double>& centre_1,  std::vector<double>& centre_2)
//{
//    std::vector<double> distance_c1_c2;
//    subtract_vectors(centre_1, centre_2, distance_c1_c2);
//    double d = sqrt(distance_c1_c2[0]*distance_c1_c2[0] + distance_c1_c2[1]*distance_c1_c2[1] + distance_c1_c2[2]*distance_c1_c2[2]) ;
//    double d_square = d*d;
//    double R_square = R*R;
//    double r_square = r*r;
//    double dRr_square = d_square - r_square + R_square;
//    double r_int = (sqrt(4 * d_square * R_square - dRr_square * dRr_square))/(2*d) ;
//    return(r_int);
//}
//
//double area_circular_segment (double R, double d)   //d = distance between centre and segment chord
//{
//    double dsquare = d*d;
//    double Rsquare = R*R;
//    double result = Rsquare*acos(d/R) - d*sqrt(Rsquare - dsquare);
//    return (result);
//}

//double volume_lens()
//{
//
//}
//
//double surface_area_lens()
//{
//
//}

void add_to_N (double N, double G, double Rg, double Rb, double Rg_secondary, int db, int dg_secondary, double volume, double SA, double dN, double Nmin, int lenN, std::vector<double>& NArray_Gmin, std::vector<double>& NArray_confs, std::vector<std::vector<double> >& NArray_quant)
{
    double rNhere = round_nearN (N , dN);   //This rounds to nearest N.
    int indN = (int) ((rNhere-Nmin)/dN);
    if(indN >= lenN){printf("Number of particles are not enough rNhere = %10.5f\n", rNhere); abort();}
    
    double minNi = NArray_Gmin[indN]; //This is 0 when no free energy has been added
    if(G <= minNi || NArray_confs[indN] == 0.0)
    {
        NArray_Gmin[indN] = G ;
        NArray_quant[indN][1] = G;
        NArray_quant[indN][2] = Rg;
        NArray_quant[indN][3] = Rb;
        NArray_quant[indN][4] = Rg_secondary;
        NArray_quant[indN][5] = (double) db;
        NArray_quant[indN][6] = (double) dg_secondary;
        NArray_quant[indN][7] = volume ;
        NArray_quant[indN][8] = SA ;
    }
    else
    {
        NArray_quant[indN][1] = minNi;
    }
    NArray_confs[indN] = NArray_confs[indN] + 1;
}



void add_to_N_spherocylinder (double N, double G, double Rg, double cyl_length, double chord_length, double Rb, int dB, double volume, double SA, std::vector<double> proj_SA, double dN, double Nmin, int lenN, std::vector<double>& NArray_Gmin, std::vector<double>& NArray_confs, std::vector<std::vector<double> >& NArray_quant)
{
    double rNhere = round_nearN (N , dN);   //This rounds to nearest N.
    int indN = (int) ((rNhere-Nmin)/dN);
    if(indN >= lenN){printf("Number of particles are not enough rNhere = %10.5f\n", rNhere); abort();}
    
    double minNi = NArray_Gmin[indN]; //This is 0 when no free energy has been added
    if(G <= minNi || NArray_confs[indN] == 0.0)
    {
        NArray_Gmin[indN] = G ;
        NArray_quant[indN][1] = G;
        NArray_quant[indN][2] = Rg;
        NArray_quant[indN][3] = cyl_length;
        NArray_quant[indN][4] = chord_length;
        NArray_quant[indN][5] = Rb;
        NArray_quant[indN][6] = (double) dB;
        NArray_quant[indN][7] = volume ;
        NArray_quant[indN][8] = SA ;
        NArray_quant[indN][9] = proj_SA[0] ;
        NArray_quant[indN][10] = proj_SA[1] ;
        NArray_quant[indN][11] = proj_SA[2] ;
    }
    else
    {
        NArray_quant[indN][1] = minNi;
    }
    NArray_confs[indN] = NArray_confs[indN] + 1;
}

void print_point(std::vector<double> point)
{
    printf("[%10.15f ,%10.15f ,%10.15f]\n", point[0],point[1],point[2]);
}

void print_points_tofile(std::vector<std::vector<double> > points, FILE* file)
{
    for(size_t i=0; i<points.size(); i++)
    {
        fprintf(file, "%10.10f\t%10.10f\t%10.10f\n", points[i][0], points[i][1] , points[i][2]);
    }
}

void print_NGDataFile (std::vector<std::vector<double> >& NArray_quant, FILE* NGDataFile)
{
    int lenN = (int) NArray_quant.size();
//    int N_quantities = (int) NArray_quant[0].size();
    for(int i = 0; i<lenN ; i++)
    {
        //[N, Gmin (/kbT), Rg, Rb, Rg_secondary, db, dg_secondary, V, SA]
        fprintf(NGDataFile,"%10.5f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.2f\t%10.2f\t%10.10f\t%10.10f\n" ,NArray_quant[i][0],NArray_quant[i][1], NArray_quant[i][2], NArray_quant[i][3], NArray_quant[i][4], NArray_quant[i][5], NArray_quant[i][6], NArray_quant[i][7] * 1e30, NArray_quant[i][8] * 1e20);
        //fprintf(NGDataFile,"%10.5f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.5f\n" ,NGmin[i][0],NGmin[i][1],NGmin[i][2],NGmin[i][3],NGmin[i][4],NGmin[i][5],NGmin[i][6]*1e30,NGmin[i][7]*1e20,NGmin[i][8]*1e20,NGmin[i][9]*1e20,NGmin[i][10]*1e30, (Nparticles[i]/1e+03));
       }
}

void print_NGDataFile_spherocylinder (std::vector<std::vector<double> >& NArray_quant, FILE* NGDataFile)
{
    int lenN = (int) NArray_quant.size();
    for(int i = 0; i<lenN ; i++)
    {
        //NArray_quant = [N, G, Rg, cyl_length, chord_length, Rb, dB, volume, SA, proj_SA [0, 1, 2]]
    
        fprintf(NGDataFile,"%10.5f\t %10.10f\t %10.10f\t %10.10f\t %10.10f\t %10.10f\t %10.2f\t %10.10f\t %10.10f\t %10.10f\t %10.10f\t %10.10f\n" ,NArray_quant[i][0],NArray_quant[i][1], NArray_quant[i][2], NArray_quant[i][3], NArray_quant[i][4], NArray_quant[i][5], NArray_quant[i][6], NArray_quant[i][7] * 1e30, NArray_quant[i][8] * 1e20, NArray_quant[i][9] * 1e20, NArray_quant[i][10] * 1e20, NArray_quant[i][11] * 1e20);
    }
}



void addQuant(std::vector<double> &destination, std::vector<double> &origin, int size_origin)
{
    for(size_t i=0; i<size_origin; i++)
    {
        destination.push_back(origin[i]);
    }
}


std::vector<double> add_double_vectors (std::vector<double>& a, std::vector<double>& b)
{
    std::vector<double> result(a.size(),0.0);
    if(a.size() != b.size())
    {
        printf("size of a = %d\t b=%d\n", (int)a.size(), (int)b.size());
        throw std::logic_error ("vector size should be same for addition\n");
    }
    else
    {
        
        for(size_t i=0; i<a.size(); i++)
        {
            result[i] = a[i] + b[i];
        }
    }
    return result;
}

std::vector<int> getLoopStartEnd (int length, int level_color)
{
    //Here I am using branch size as 2 because we are dividing into two branches at each level so making a binary division at each node.
    std::vector<int> loop_end_points (2);
    int start, end;
    int branch_rank;
    int branch_size = 2;
    int chunk_per_process = (int)(length/branch_size) ;
    
    if (level_color % 2 == 0)
    {
        branch_rank = 0;
    }
    else {branch_rank = 1;}
    
    if(length % branch_size != 0)
    {
        int extra_points = length % branch_size;
        //Adding the remaining extra points to the last branch rank
        if(branch_rank == branch_size - 1)
        {
            start = branch_rank*chunk_per_process;
            end = (branch_rank+1)*chunk_per_process - 1;
            end = end + extra_points ;
        }
        else
        {
            start = branch_rank*chunk_per_process;
            end = (branch_rank+1)*chunk_per_process - 1;
        }
    }
    else
    {
        start = branch_rank*chunk_per_process;
        end = (branch_rank+1)*chunk_per_process - 1;
        
    }
    loop_end_points[0] = start;
    loop_end_points[1] = end;
    return loop_end_points;
}
