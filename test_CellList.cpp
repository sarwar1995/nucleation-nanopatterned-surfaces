//
//  test_CellList.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 11/5/20.
//  Copyright © 2020 Sarwar Hussain. All rights reserved.
//

#include <stdio.h>
#include "test_CellList.hpp"
#include "MC.hpp"
#include "miscformulas.hpp"



int main(int argc, const char * argv[]){
    
    /* Setup */
    
    FILE* test_points_file = fopen(argv[1], "w");
    int num_points = 128;
    printf("It's in here = %d\n", num_points);
    std::vector<std::vector<double> > box;
    box.resize(3);
    printf("box size = %d\n", (int) box.size());
    int aSeed[3];
    aSeed[0] = 20201; aSeed[1] = 20202; aSeed[2] = 20203;
    double box_size = 2.0;
    std::vector<double> R_cutoff;
    R_cutoff.resize(3);
    R_cutoff[0] = 0.5; R_cutoff[1] = 0.5; R_cutoff[2] = 0.5;
    printf("R_cutoff[0] = %10.2f\tR_cutoff[1] = %10.2f\tR_cutoff[2] = %10.2f\n", R_cutoff[0],  R_cutoff[1],  R_cutoff[2]);
    
    for(size_t i=0; i<box.size(); i++){
        box[i].resize(2);
        box[i][0] = 0.0; box[i][1] = box_size;
        printf("box[%d][0] = %10.2f\t box[%d][1] = %10.2f\n", (int)i, box[i][0], (int)i, box[i][1]);
        //if(i!=2){box[i][0] = -box_size/2.0; box[i][1] = box_size/2.0;}
        //else{box[i][0] = 0.0; box[i][1] = box_size;}
    }
    
     printf("box done\n");
    
    MC myMC (num_points, box, aSeed);
    myMC.print_points(test_points_file);
    std::vector<std::vector<double> > x = myMC.get_points();
    printf("x.size()=%d\n",(int) x.size());
    CellList myCellList (x, R_cutoff, box);
    
    /****************************************************/
    printf("total_cells=%d\n",myCellList.get_total_cells());
    
    std::vector<double> aPoint = {1.0444962406  ,  0.9879912926  ,  0.9789679438};
    printf("aPoint = [%10.10f,\t%10.10f\t,%10.10f]\n",aPoint[0],aPoint[1],aPoint[2]);
//
    
    //calc_cell_ind
    int ind = myCellList.calc_cell_ind(aPoint);
    printf("ind=%d\n",ind);
    std::vector<std::vector<double> > cell_pts;
    
    //All points in the cell with index=ind cell_points
    cell_pts = myCellList.cell_points(ind);
    myCellList.print_cell(ind, 1);

    //cell_ngbs
    std::vector<int> ngbs = myCellList.cell_ngbs(ind); //51
    printf("ngbs.size() = %d\t",(int)ngbs.size());
    for(size_t k=0; k<ngbs.size(); k++)
    {
        printf("ngbs[%d] = %d\t",(int) k,(int)ngbs[k]);
    }
    
    //std::vector<double> aPoint2 = {0.3860448006 ,   1.6126372947 ,   1.8402668191};
    std::vector<std::vector<double> >points_coords = {cell_pts[4]};//{cell_pts[4], x[32]};
    //printf("x[32] =\t");print_point(x[32]);
    std::vector<std::vector<double> > ngb_points = myCellList.calc_neighbors (points_coords);
    printf("is in calc_neighbors\n");
    printf("%d\t%d\n", InVector(ngb_points, cell_pts[4]), InVector(ngb_points, x[32]));
    
    
    printf("ngb_points.size() = %d\t",(int)ngb_points.size());
    printf("calc_neighbors\n");
    for(size_t k=0; k<ngb_points.size(); k++)
    {
        print_point(ngb_points[k]);
    }
    
    printf("cells points\n");
    int total_points=0;
    for(size_t i=0; i<ngbs.size(); i++)
    {
        total_points += myCellList.cell_size((int)ngbs[i]);
        printf("i = %d\t cell_size=%d\n",(int)i,myCellList.cell_size((int)ngbs[i]));
        myCellList.print_cell((int)ngbs[i], 1);
    }
    printf("Total points = %d\n",total_points);
    
    /**********************************/
    for(size_t u=0; u<cell_pts.size();u++)
    {
        int isinvector = InVector(points_coords, cell_pts[u]);
        printf("isinvector = %d\n",isinvector);
        print_point(cell_pts[u]);
    }
    print_point(cell_pts[4]);
    print_point(points_coords[0]);
    //bool a = (cell_pts[4][0] == points_coords[0][0]);

}
