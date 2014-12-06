////////////////////////////////////////////////////////////////////
//                                                                //
//            Functions for Computational Geometry assignments    //
//                                                                //
////////////////////////////////////////////////////////////////////

#include <string>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
using namespace std;

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

struct point
{
    double x;
    double y;
};

struct minMax
{
    double minx = 999;
    double miny = 999;
    double maxx = -999;
    double maxy = -999;
};



////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

double min(double a, double b)
{
    if (a < b)
        return a;
    return b;
}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

double max(double a, double b)
{
    if (a > b)
        return a;
    return b;
}

////////////////////////////////////////////////////////////////////
//          Cross product of vectors a and b                      //
////////////////////////////////////////////////////////////////////

double cross(point a, point b)
{
    return a.x * b.y - a.y * b.x;
}

////////////////////////////////////////////////////////////////////
//      Is c on the line segment from a to b?                     //
//      Assumes c is on the line from a to b                      //
////////////////////////////////////////////////////////////////////

bool on(point a, point b, point c)
{
    if (((min(a.x, b.x) <= c.x) && (c.x <= max(a.x, b.x))) &&
            ((min(a.y, b.y) <= c.y) && (c.y <= max(a.y, b.y))))
        return true;

    return false;
}

////////////////////////////////////////////////////////////////////
//       Compute cross product of vector from a to b and          //
//       vector from b to c                                       //
//       Set to zero if close to zero                             //
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
//    Is the turn from a to b to c clockwise?                     //
////////////////////////////////////////////////////////////////////

bool turn(point a, point b, point c)
{
    return (direction(a, b, c) > 0);
}

////////////////////////////////////////////////////////////////////
//    Is the turn from p[0] to p[1] to p[2] clockwise?            //
////////////////////////////////////////////////////////////////////

bool clockwise(point p[])
{
    return turn(p[0], p[1], p[2]);
}

////////////////////////////////////////////////////////////////////
//     From text, page 937 (second edition)                       //
//     Does line segment from p1 to p2 intersect                  //
//     line segment from p3 to p4?                                //
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
//          Area of polygon by adding/subtracting trapezoids      //
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


void check_args(int argc)
{
    if(argc != 3)
    {
        cout << "Error: invalid number of arguments\nUsage: geometry  <inputFile>  <outputFile>" << endl;
        exit(1);
    }
}

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


    for(unsigned int i = 0; i < input.size() ; i++)
        cout << input[i].x << " " << input[i].y << endl;

    cout << "boundry.maxx = " << boundry.maxx << endl;
    cout << "boundry.maxy = " << boundry.maxy << endl;
    cout << "boundry.minx = " << boundry.minx << endl;
    cout << "boundry.miny = " << boundry.miny << endl;





    fin.close();
    fout.close();
    return 0;
}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

