////////////////////////////////////////////////////////////////////
//  Authors: Matthew Rames and Zachary Pierson
//          Computational Geometry Assignment.
////////////////////////////////////////////////////////////////////

#include <string>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
using namespace std;

//structure that holds x and y coordinates of a point
struct point
{
    double x;
    double y;
};

//this contains the number of interior and boundry points
struct IBpoints
{
    unsigned int bPoints = 0; //# of boundry points
    unsigned int iPoints = 0; //# of interior points
    vector<point> vecB;
    vector<point> vecI;
};

//contains min and max x and y values
struct minMax
{
    double minx = 999;
    double miny = 999;
    double maxx = -999;
    double maxy = -999;
};



////////////////////////////////////////////////////////////////////
//  Author: Dr. Corwin
//          NOT USED
////////////////////////////////////////////////////////////////////
double min(double a, double b)
{
    if (a < b)
        return a;
    return b;
}


////////////////////////////////////////////////////////////////////
//  Author: Dr. Corwin
//          NOT USED
////////////////////////////////////////////////////////////////////
double max(double a, double b)
{
    if (a > b)
        return a;
    return b;
}

////////////////////////////////////////////////////////////////////
//  Author: Dr. Corwin
//          Cross product of vectors a and b                      //
//          USED
////////////////////////////////////////////////////////////////////
double cross(point a, point b)
{
    return a.x * b.y - a.y * b.x;
}

////////////////////////////////////////////////////////////////////
//  Author: Dr. Corwin
//      Is c on the line segment from a to b?                     //
//      Assumes c is on the line from a to b                      //
//      NOT USED
////////////////////////////////////////////////////////////////////
bool on(point a, point b, point c)
{
    if (((min(a.x, b.x) <= c.x) && (c.x <= max(a.x, b.x))) &&
            ((min(a.y, b.y) <= c.y) && (c.y <= max(a.y, b.y))))
        return true;

    return false;
}

////////////////////////////////////////////////////////////////////
//  Author: Dr. Corwin
//       Compute cross product of vector from a to b and          //
//       vector from b to c                                       //
//       Set to zero if close to zero                             //
//       NOT USED
////////////////////////////////////////////////////////////////////
double direction(point a, point b, point c)
{
    point ab;
    point bc;
    double result;

    ab.x = b.x - a.x;
    ab.y = b.y - a.y;
    bc.x = c.x - b.x;
    bc.y = c.y - b.y;
    result =  cross(ab, bc);
    if (fabs(result) < 1.0e-6)
        result = 0.0;
    return result;
}

////////////////////////////////////////////////////////////////////
//  Author: Dr. Corwin
//    Is the turn from a to b to c clockwise?                     //
//    NOT USED
////////////////////////////////////////////////////////////////////
bool turn(point a, point b, point c)
{
    return (direction(a, b, c) > 0);
}

////////////////////////////////////////////////////////////////////
//  Author: Dr. Corwin
//    Is the turn from p[0] to p[1] to p[2] clockwise?            //
//    NOT USED
////////////////////////////////////////////////////////////////////
bool clockwise(point p[])
{
    return turn(p[0], p[1], p[2]);
}

////////////////////////////////////////////////////////////////////
//  Author: Dr. Corwin
//     From text, page 937 (second edition)                       //
//     Does line segment from p1 to p2 intersect                  //
//     line segment from p3 to p4?                                //
//     NOT USED
////////////////////////////////////////////////////////////////////
bool intersect(point p1, point p2, point p3, point p4)
{
    double d1;
    double d2;
    double d3;
    double d4;

    d1 = direction(p3, p4, p1);
    d2 = direction(p3, p4, p2);
    d3 = direction(p1, p2, p3);
    d4 = direction(p1, p2, p4);

    if ((((d1 > 0) && (d2 < 0)) || ((d1 < 0) && (d2 > 0))) &&
            (((d3 > 0) && (d4 < 0)) || ((d3 < 0) && (d4 > 0))))
        return true;

    if ((d1 == 0) && on(p3, p4, p1))
        return true;

    if ((d2 == 0) && on(p3, p4, p2))
        return true;

    if ((d3 == 0) && on(p1, p2, p3))
        return true;

    if ((d4 == 0) && on(p1, p2, p4))
        return true;

    return false;
}

////////////////////////////////////////////////////////////////////
//  Author: Dr. Corwin
//          Area of polygon by adding/subtracting trapezoids      //
//          NOT USED
////////////////////////////////////////////////////////////////////
double area(point p[], int n)
{
    int i;
    int j;
    double result;

    result = 0;
    for (i = 0; i < n; i++)
    {
        j = (i + 1) % n;
        result += p[i].x * p[j].y;
        result -= p[i].y * p[j].x;
    }
    return fabs(result / 2);
}


////////////////////////////////////////////////////////////////////
//  This function checks the number of arguments to main.
////////////////////////////////////////////////////////////////////
void check_args(int argc)
{
    if(argc != 3)
    {
        cout << "Error: invalid number of arguments\nUsage: geometry  <inputFile>  <outputFile>" << endl;
        exit(1);
    }
}

////////////////////////////////////////////////////////////////////
//  This function tries to open the files that are passed into 
//  the function. This function appends .in and .out to the input
//  and output files before opening
////////////////////////////////////////////////////////////////////
void open_files(ifstream &fin, ofstream &fout, char **argv)
{
    string infile = argv[1];
    string outfile = argv[2];
    

    //Open and check if input file exists
    fin.open(infile.append(".in").c_str(), ios::in);
    if(!fin)
    {
        cout << "File " << argv[1] << "could not be opened" << endl;
        exit (1);
    }




    outfile.append(".out");

    //Open and check if output file exists
    fout.open(outfile.c_str(), ios::out);
    if(!fout)
    {
        cout << "File " << argv[2] << "could not be opened" << endl;
        exit (1);
    }
}

////////////////////////////////////////////////////////////////////
//  This function takes in points from a file and stores them in
//  a vector of points. This also keeps track of the minimum and
//  maximum x and y values inputed.
////////////////////////////////////////////////////////////////////
void get_coords(ifstream &fin, vector<point> &input, minMax &boundry)
{
    point newPoint;
    double c;

    //Read points from input file
    while( fin >> c  )
    {
        newPoint.x = c;
        fin >> c;
        newPoint.y = c;

        input.push_back(newPoint);

        if(newPoint.x < boundry.minx)
            boundry.minx = newPoint.x;
        if(newPoint.y < boundry.miny)
            boundry.miny = newPoint.y;
        if(newPoint.x > boundry.maxx)
            boundry.maxx = newPoint.x;
        if(newPoint.y > boundry.maxy)
            boundry.maxy = newPoint.y;
    }
}

////////////////////////////////////////////////////////////////////
//  This function checks to see if two points are equal
////////////////////////////////////////////////////////////////////
bool isEqual(point a, point b)
{
    return ((a.x == b.x) && (a.y == b.y));
}

////////////////////////////////////////////////////////////////////
//  This function returns a vector c where c = a - b
////////////////////////////////////////////////////////////////////
point p_minus(point a, point b)
{
    point result;
    result.x = a.x - b.x;
    result.y = a.y - b.y;

    return result;
}

////////////////////////////////////////////////////////////////////
//  findIBpoints is the main function for finding the interior and
//  boundry points. It takes a vector of counterclockwise order 
//  data points and their associated minimum and maximum x and y 
//  values. By taking a point that is known to be outside the 
//  polygon, and count the number of intersections to the point in
//  question; if the number of intersections is odd, then the point
//  is inside the polygon else it is outside the polygon. However,
//  if the path to the point in question intersects a vertex, then 
//  we move our outer point in the posotive x direction and try 
//  again. 
//  This function returns an object that contains the number of 
//  interior and boundry points and vectors of corresponding points.
////////////////////////////////////////////////////////////////////
IBpoints findIBpoints(vector<point> &vert, minMax min_max)
{
    IBpoints points;    //contains info about interior and boundry points
    point p;            //the current point in question
    p.x = min_max.minx; //setting the current point to lowest right value
    p.y = min_max.miny;
    point q, r, s, next;    //values used to compute BI points
    double rxs, t, u;       //values used to compute BI points

    int intersections = 0;  //used to determine if point is an interior point
    point op;               //outter point known to be outside the polygon

    //set the outter point to be outside polygon
    op.x = min_max.maxx + 1;
    op.y = min_max.maxy + 1;
    
    //number of possible points to check
    int numBox = ((min_max.maxx+1) - min_max.minx) * ((min_max.maxy+1) - min_max.miny); 

    //set the number of boundry points to the number of verticies
    points.bPoints = vert.size();

    //loop though all questionable points
    for(int i = 0; i < numBox; i++)
    {
        intersections = 0;

        //loop through all lines in polygon
        //this is the heart of finding the BI points
        for(int j = 0; j < (int)vert.size(); j++)
        {
            q = vert[j]; //is the first point of the line segment
            
            r = p_minus(op, p); // r is a vector from p to the outer point

            //Checks to see if vector is at end if it is set the next point
            //to be the beginning point
            if(j == ((int)vert.size() - 1))
                next = vert[0];
            else
                next = vert[j+1];

            s = p_minus(next, q); // s is the vector from q to the next q


            //If the the point in question is a vertex push the point and break
            if(isEqual(p, q) || isEqual(p, next)) //point in question is vertex
            {
                points.vecB.push_back(p);
                intersections = 0;
                break;
            }

            //The main equations to dertmine intersection used are as follows
            //t = (q-p) x s / (r x s)
            //u = (p-q) x r / (s x r)
            //If r x s == 0 there is no intersection
            //If r x s != 0 && 0 <= t <= 1 && 0 <= u <= 1 then there is
            //an intersection.
            rxs = cross(r, s);
            if(rxs == 0)    //lines don't intersect
                continue;

            t = (cross(p_minus(q, p), s) / rxs);
            if( t < 0 || t > 1 )    //lines don't intersect
                continue;

            // *Note* s x r = -(r x s)
            u = (cross(p_minus(p, q), r) / (-1 * rxs));
            if( u < 0 || u > 1 )    //lines don't intersect
                continue;

            //the point at with intersection occurrs is (q + us) = (p + tr)
            //if u is 0 => intersection occurrs at q
            //if u is 1 => intersection occurrs at next q
            if(u == 1 || u == 0) //hit vertex
            {
                op.x++; //move outside point over in x direction by 1
                j--;    //check the same line again
                continue;
            }

            //the point at with intersection occurrs is (q + us) = (p + tr)
            //if t is 0, intersection occurrs at p which is an integer boundry point
            //if t is 1, intersection occurrs at op *this should never happen
            //because op is outside of the possible points*
            if(t==0)
            {
                points.vecB.push_back(p);
                points.bPoints++;
                intersections = 0;
                break;
            }
            intersections++;    //keep track of number of intersections
            
        }
        
        //if intersections is odd, then p is an interior point
        if(intersections & 1)
        {
            points.vecI.push_back(p);
            points.iPoints++;
        }

        //update current point
        if( p.y+1 > min_max.maxy )
        {
            p.y = min_max.miny;
            p.x++;
        }
        else
            p.y++;

    }
    
    return points;
}

////////////////////////////////////////////////////////////////////
//  This function outputs information regarding the area and the
//  number of interior and boundry points for picks theorem.
////////////////////////////////////////////////////////////////////
void outputfile(char *infile, ofstream &fout, IBpoints a)
{
    double area;
    area = a.iPoints + (a.bPoints / 2.0) - 1;

    fout << "Reading input from file " << infile << ".in" << endl;
    fout << "Area = " << area << endl;
    fout << "I(P) = " << a.iPoints << endl;
    fout << "B(P) = " << a.bPoints << endl;

    fout << "Interior points" << endl;
    for(unsigned int i = 0; i < a.vecI.size() ; i++)
    {
        fout << "\t" << a.vecI[i].x << " " << a.vecI[i].y << endl;
    }

    fout << "\nBoundary points" << endl;
    for(unsigned int i = 0; i < a.vecB.size() ; i++)
    {
        fout << "\t" << a.vecB[i].x << " " << a.vecB[i].y << endl;
    }
}

////////////////////////////////////////////////////////////////////
//  Main takes in arguments for an input and output file. The
//  first parameter is the name of the inputfile without the ".in"
//  extension. The second perameter is the name of the wanted output
//  file. This file will be appended with a ".out" extension.
//
//  This program will take a list of vertecies in counterclockwise
//  order that forms a simple polygon and compute the area using 
//  picks theorem.
////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
    //Check command line arguments
    check_args(argc);

    ifstream fin;
    ofstream fout;

    //Open and check input and output files
    open_files(fin, fout, argv);

    minMax boundry;

    //Fill vector with points from input file
    vector<point> input;
    get_coords(fin, input, boundry);

    IBpoints a = findIBpoints(input, boundry);

    outputfile(argv[1], fout, a);



    fin.close();
    fout.close();
    return 0;
}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

