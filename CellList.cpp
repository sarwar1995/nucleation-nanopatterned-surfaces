//
//  CellList.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 8/5/20.
//  Copyright © 2020 Sarwar Hussain. All rights reserved.
//

#include "CellList.hpp"

CellList :: CellList(){
    head.resize(0);
    list.resize(0);
}

CellList::~CellList(){
    head.clear();
    list.clear();
}

CellList::CellList(std::vector<std::vector<double> > x, std::vector<double> R_cutoff, std::vector<std::vector<double> > MCbox)
{    
    x_points = x;
    Rc = R_cutoff;
    box = MCbox;
    m_cells.resize(Rc.size());
    n_points = (int) x_points.size();
    for (size_t i=0; i<Rc.size(); i++)
    {
        
        double box_i = MCbox[i][1] - MCbox[i][0] ; //length_x = box_end_x - box_start_x
        m_cells[i] = (int) (box_i/Rc[i]) ;
        
    }
    total_cells = m_cells[0]*m_cells[1]*m_cells[2];
    head.resize(total_cells);
    list.resize(n_points);
    generate();
    
}

int CellList::calc_cell_ind(std::vector<double> & X)
{
    // box[i][0] is not at origin, therefore \alpha_th component of index of a point is calculated by calculating the distance from box corner i.e. box[i][0]
    double x_i = X[0];
    double y_i = X[1];
    double z_i = X[2];
    int index_xi = (int) ((x_i - box[0][0])/Rc[0]);
    int index_yi = (int) ((y_i - box[1][0])/Rc[1]);
    int index_zi = (int) ((z_i - box[2][0])/Rc[2]);
    int c = index_zi + index_yi*m_cells[2] + index_xi*m_cells[2]*m_cells[1] ;
    return c;
}

void CellList::generate()
{
    for(int i=0; i<total_cells; i++)
    {
        head[i] = -1; //-1 indicates that the head list is empty
    }
    for (int i=0; i<n_points; i++)
    {
        int c = calc_cell_ind(x_points[i]);
        if(c<0 || c>=total_cells){printf("Incorrect c = %d value at x_points = [%10.10f, %10.10f, %10.10f]\n", c, x_points[i][0], x_points[i][1], x_points[i][2]); abort();}
        list[i] = head[c];
        head[c] = i;
    }
}

//Getting all points inside a cell
std::vector<std::vector<double> > CellList::cell_points(int c)
{
    std::vector<std::vector<double> > cell_points ;
    std::set<std::vector<double> > set_cell_points ;
    int list_i = head[c];
    while (list_i != -1)
    {
        cell_points.push_back(x_points[list_i]);
        list_i = list[list_i];
    }
    
    return cell_points;
    
}
//Getting the size of a cell
int CellList::cell_size(int c)
{
    int count=0;
    int list_i = head[c];
    while (list_i != -1)
    {
        count++;
        list_i = list[list_i];
    }
    return count;
}


//Gives the ngbs of a cell EXCLUDING itself
std::vector<int> CellList::cell_ngbs (int c)
{
    std::vector<int> CellNgbs;
    int Cx, Cy, Cz; //Vector cell indices
    Cx = (int) (c/(m_cells[2]*m_cells[1]));
    Cy = (int) ((c/m_cells[2])%m_cells[1]);
    Cz = (int) (c%m_cells[2]);
//    printf("Cx=%d\tCy=%d\tCz=%d\n",Cx,Cy,Cz);
//    printf("m_cells[0]=%d\tm_cells[1]=%d\tm_cells[2]=%d\n",m_cells[0],m_cells[1],m_cells[2]);
    int mc[3];  //vector indices of neighboring cells (including self)
    int c_ngb;
//    int count=0;
    for(mc[0] = Cx-1; mc[0]<=Cx+1; mc[0]++)
    {
        if(mc[0]<0 || mc[0] >= m_cells[0]){continue;}
        else{
            for(mc[1] = Cy-1; mc[1]<=Cy+1; mc[1]++)
            {
                if(mc[1]<0 || mc[1] >= m_cells[1]){continue;}
                else{
                    for(mc[2] = Cz-1; mc[2]<=Cz+1; mc[2]++)
                    {
                        if(mc[2]<0 || mc[2] >= m_cells[2]){continue;}
                        else{
                            c_ngb = mc[2] + mc[1]*m_cells[2] + mc[0]*m_cells[2]*m_cells[1] ;
                            if(c_ngb != c)
                            {CellNgbs.push_back(c_ngb);}
                        }
                    }
                }
            }
        }
        
    }
    return CellNgbs;
}

std::vector<std::vector<double> > CellList::get_cell_boundaries (int c, int yesPrint)
{
    std::vector<std::vector<double> > cell_bounds(3);
    cell_bounds[0].resize(2);
    cell_bounds[1].resize(2);
    cell_bounds[2].resize(2);
    int Cx, Cy, Cz; //Vector cell indices
    Cx = (int) (c/(m_cells[2]*m_cells[1]));
    Cy = (int) ((c/m_cells[2])%m_cells[1]);
    Cz = (int) (c%m_cells[2]);
    
    cell_bounds[0][0] = (Cx * Rc[0]) + box[0][0];
    cell_bounds[0][1] = cell_bounds[0][0] + Rc[0];
    
    
    cell_bounds[1][0] = (Cy * Rc[1]) + box[1][0];
    cell_bounds[1][1] = cell_bounds[1][0] + Rc[1];
    
    cell_bounds[2][0] = (Cz * Rc[2]) + box[2][0];
    cell_bounds[2][1] = cell_bounds[2][0] + Rc[2];
    
    if(yesPrint)
    {
        printf("cell_x = [%10.10f, %10.10f\n]", cell_bounds[0][0], cell_bounds[0][1]);
        printf("cell_y = [%10.10f, %10.10f\n]", cell_bounds[1][0], cell_bounds[1][1]);
        printf("cell_z = [%10.10f, %10.10f\n]", cell_bounds[2][0], cell_bounds[2][1]);
    }
    return cell_bounds;
}


std::vector<std::vector<double> > CellList::calc_neighbors (std::set<std::vector<double> >points_coords)
{
    std::vector<std::vector<double> >vector_points_coords(points_coords.begin(), points_coords.end());
    
    return calc_neighbors(vector_points_coords);
}


//Gives the coordinates of points in the neighboring cells and self cells besides the points already present in point_coords
std::vector<std::vector<double> > CellList::calc_neighbors (std::vector<std::vector<double> >points_coords)
{
    /*Create a set of point coords vector.*/
    std::set<std::vector<double> > set_points_coords (points_coords.begin(), points_coords.end());
    std::set<std::vector<double> > neigh_coords; //Coords of all points contained in the cells that are not in points_coords
    std::set<int> cells;         //All cells in which points are present
//    std::vector<int> new_ngb_cells; //New neighbor cells of points not already present in cells
    
    printf("points_coords size = %d\n",(int) points_coords.size());
    //Calculating all the distinct cells of points_coord
    for (size_t i=0; i<points_coords.size(); i++)
    {
        int c = calc_cell_ind(points_coords[i]);
        if(!cells.count(c))
        {
            cells.insert(c);
            std::vector<std::vector<double> > CellPoints = cell_points(c);
            for(size_t k=0; k<CellPoints.size(); k++)
            {
                if(!set_points_coords.count(CellPoints[k]))
                {
                    neigh_coords.insert(CellPoints[k]);
                }
            }
            
       }
    }
    
    //For these cells, calculating ngbs of each cell that are not already contained in the cell
    std::set<int>::iterator set_iterate;
    for(set_iterate = cells.begin(); set_iterate != cells.end(); set_iterate++)
    {
        std::vector<int> CellNgbs = cell_ngbs (*set_iterate);
        for(size_t j=0; j<CellNgbs.size(); j++)
        {
            if(!cells.count(CellNgbs[j]))
            {
                std::vector<std::vector<double> > CellNgb_points = cell_points(CellNgbs[j]);
                neigh_coords.insert(CellNgb_points.begin(),CellNgb_points.end());
            }
        }
        
//        std::vector<int> CellNgbs = cell_ngbs (cells[i]);
//        for(size_t j=0; j<CellNgbs.size(); j++)
//        {
//            if(!InVector(cells, CellNgbs[j]))
//            {
////                new_ngb_cells.push_back(CellNgbs[j]);
//                std::vector<std::vector<double> > CellNgb_points = cell_points(CellNgbs[j]);
//                appendTovec(neigh_coords, CellNgb_points); //Add all points in the new_ngb_cells to neigh_coords
//            }
//        }
    }
    printf("neigh_coords size = %d\n",(int) neigh_coords.size());
    std::vector<std::vector<double> > vector_neigh_coords (neigh_coords.begin(), neigh_coords.end());
    return vector_neigh_coords;
    
}

//Nicely print the indices of particles in cell c
void CellList::print_cell(int c, int yesPoints)
{
    int list_i = head[c];
    printf("Cell c = %d\n[",c);
    while (list_i != -1)
    {
        if(yesPoints==1)
        {printf("{%10.5f,\t%10.5f,\t%10.5f}\t",x_points[list_i][0],x_points[list_i][1],x_points[list_i][2]);}
        else{printf("%d,\t",list_i);}
        list_i = list[list_i];
    }
    printf("]\n");
}

//Print the points in the neighboring cells of a cell c
void CellList::print_neighboring_points(int c)
{
    
}
