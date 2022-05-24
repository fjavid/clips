#ifndef CLIP_H
#define CLIP_H

#include <iostream>
#include <Eigen/Core>
#include <cmath>

class Body
{
    public:
    Body(){};
    Body(double _density, double _friction, double _restitution) : den(_density), fric(_friction), rest(_restitution){};
    virtual ~Body() = default;
    std::string vis_type = "none";
    double envelope = 0.001;
    double margin = 0.0001;
    double den;                                 // Density of the clip material
    double fric;                                // Friction of the clip material
    double rest;                             // Restitution of the clip mateial
    double mass;
    double volume;
    double inertia_xx;
    double inertia_yy;
    double inertia_zz;
    std::vector<std::vector<double>> vertices;
    std::vector<std::vector<int>> elems;
};

class Clip : public Body
{
    public:
    Clip(){}; 
    Clip(double _h, double _w, double _rod_r, double _gap_h, double _fil_r, double _density, double _friction, double _restitution)
    : h(_h), w(_w), r(_rod_r), g(_gap_h), fil_r(_fil_r), Body(_density, _friction, _restitution)
    {
        // Body(_density, _friction, _restitution);
        vis_type = "wire";
        cs_area = M_PI * r * r;
        double volume_h = cs_area * h;
        double volume_w = cs_area * w;
        double volume_g = cs_area * g;
        volume = 2.0*volume_h + 2.0*volume_w - volume_g;
        double mass_h = den * volume_h;
        double mass_w = den * volume_w;
        double mass_g = den * volume_g;
        mass = den * volume;
        inertia_xx = 2.0*mass_h/12.0 * h * h
                            + 2.0*mass_w/2.0 * r * r + 2.0*mass_w * (0.5*h) * (0.5*h)
                            - mass_g/12.0 * g * g;
        inertia_yy = 2.0*mass_h/2.0 * r * r + 2.0*mass_h * (0.5*w) * (0.5*w)
                            + 2.0*mass_w/12.0 * w * w
                            - mass_g/2.0 * r * r - mass_g * (0.5*w) * (0.5*w);
        inertia_zz = 2.0*mass_h/12.0 * h * h + 2.0*mass_h * (0.5*w) * (0.5*w)
                            + 2.0*mass_w/12.0 * w * w + 2.0*mass_w * (0.5*h) * (0.5*h)
                            - mass_g/12.0 * g * g - mass_g * (0.5*w) * (0.5*w);

        // to match the way clips are modelled in chrono we add a length of 'r' (radius) to the end of each section. The edges of each section will not coincide then...
        // top edge
        vertices.push_back(std::vector<double>({0.5*w+r, 0.5*h, 0.0}));
        vertices.push_back(std::vector<double>({-0.5*w-r, 0.5*h, 0.0}));
        // left edge
        vertices.push_back(std::vector<double>({-0.5*w, 0.5*h+r, 0.0}));
        vertices.push_back(std::vector<double>({-0.5*w, -0.5*h-r, 0.0}));
        // bottom edge
        vertices.push_back(std::vector<double>({-0.5*w-r, -0.5*h, 0.0}));
        vertices.push_back(std::vector<double>({0.5*w+r, -0.5*h, 0.0}));
        // right edges
        if (g > 2.0*r)
        {
            // lower part
            vertices.push_back(std::vector<double>({0.5*w, -0.5*h-r, 0.0}));
            vertices.push_back(std::vector<double>({0.5*w, -0.5*g, 0.0}));
            // upper part
            vertices.insert(vertices.begin(), std::vector<double>({0.5*w, 0.5*h+r, 0.0}));
            vertices.insert(vertices.begin(), std::vector<double>({0.5*w, 0.5*g, 0.0}));
        }
        else
        {
            g = 0.0;
            vertices.push_back(std::vector<double>({0.5*w, -0.5*h-r, 0.0}));
            vertices.push_back(std::vector<double>({0.5*w, 0.5*h+r, 0.0}));
        }
        for (int ed=0; ed<0.5*vertices.size(); ed++)
            elems.push_back(std::vector<int>({2*ed, 2*ed+1}));
        
        std::cout << "mass : " << mass << " inertia XX : " << inertia_xx << " inertia YY : " << inertia_yy << " inertia ZZ:  " << inertia_zz << std::endl;
    };
    double r;                                   // Wire radius
    double g;                                   // Clips opening 
    double h;                                       // X-direction span
    double w;                                       // Y-direction span
    double fil_r;                                   // fillet radius
    double cs_area;
};

class Box : public Body
{
    public:
    Box(){};
    Box(double _int_x, double _int_y, double _int_z, double _thick, double _density, double _friction, double _restitution)
    : int_x(_int_x), int_y(_int_y), int_z(_int_z), thick(_thick), Body(_density, _friction, _restitution)
    {
        // Body(_density, _friction, _restitution);
        vis_type = "face";
        if (thick <= 0.0)
        {
            thick = 0.1*std::min(std::min(int_x, int_y), int_z);
        }
        
        double volume_out = (int_x+2.*thick) * (int_y+2.*thick) * (int_z+2.*thick);
        double volume_in = int_x * int_y * int_z;
        double volume = volume_out - volume_in;
        double mass = den * volume;
        double mass_out = den * volume_out;
        double mass_in = den * volume_in;
        inertia_xx = mass_out/12.0 * ((int_y+2.*thick)*(int_y+2.*thick) + (int_z+2.*thick)*(int_z+2.*thick))
                            - mass_in/12.0 * (int_y*int_y + int_z*int_z);
        inertia_yy = mass_out/12.0 * ((int_x+2.*thick)*(int_x+2.*thick) + (int_z+2.*thick)*(int_z+2.*thick))
                            - mass_in/12.0 * (int_x*int_x + int_z*int_z);
        inertia_zz = mass_out/12.0 * ((int_y+2.*thick)*(int_y+2.*thick) + (int_x+2.*thick)*(int_x+2.*thick))
                            - mass_in/12.0 * (int_y*int_y + int_x*int_x);
        
        vertices.push_back(std::vector<double>({-0.5*int_x, -0.5*int_y, -0.5*int_z}));
        vertices.push_back(std::vector<double>({-0.5*int_x, 0.5*int_y, -0.5*int_z}));
        vertices.push_back(std::vector<double>({-0.5*int_x, 0.5*int_y, 0.5*int_z}));
        vertices.push_back(std::vector<double>({-0.5*int_x, -0.5*int_y, 0.5*int_z}));

        vertices.push_back(std::vector<double>({0.5*int_x, -0.5*int_y, -0.5*int_z}));
        vertices.push_back(std::vector<double>({0.5*int_x, 0.5*int_y, -0.5*int_z}));
        vertices.push_back(std::vector<double>({0.5*int_x, 0.5*int_y, 0.5*int_z}));
        vertices.push_back(std::vector<double>({0.5*int_x, -0.5*int_y, 0.5*int_z}));

        elems.push_back(std::vector<int>({0, 1, 2}));
        elems.push_back(std::vector<int>({0, 2, 3}));
        elems.push_back(std::vector<int>({0, 1, 5}));
        elems.push_back(std::vector<int>({0, 3, 7}));
        elems.push_back(std::vector<int>({0, 4, 5}));
        elems.push_back(std::vector<int>({0, 4, 7}));
        elems.push_back(std::vector<int>({4, 5, 6}));
        elems.push_back(std::vector<int>({4, 7, 6}));
        elems.push_back(std::vector<int>({5, 6, 1}));
        elems.push_back(std::vector<int>({6, 1, 2}));
        elems.push_back(std::vector<int>({2, 3, 6}));
        elems.push_back(std::vector<int>({3, 7, 6}));
        
    }
    double int_x;
    double int_y;
    double int_z;
    double thick;

};

class Cylinder : public Body
{
    public:
    Cylinder(){};
    Cylinder(double _int_r, double _int_h, double _thick, bool _cap, double _density, double _friction, double _restitution)
    : int_r(_int_r), int_h(_int_h), thick(_thick), cap(_cap), Body(_density, _friction, _restitution)
    {
        // Body(_density, _friction, _restitution);
        vis_type = "face";
        if (thick <= 0.0)
        {
            thick = 0.05*int_r;
        }
        
        double volume_hollow_cyl = M_PI * ((int_r + thick) * (int_r + thick) - int_r * int_r) * int_h;
        double volume_lcap = M_PI * (int_r + thick) * (int_r + thick) * thick;
        double volume_rcap =  cap ? volume_lcap : 0.0;
        double volume = volume_hollow_cyl + volume_lcap + volume_rcap;
        double mass = den * volume;
        double mass_hollow_cyl = den * volume_hollow_cyl;
        double mass_lcap = den * volume_lcap;
        double mass_rcap = den * volume_rcap;
        inertia_xx = mass_hollow_cyl/12.0 * (3.0*(int_r*int_r + (int_r+thick)*(int_r+thick)) + int_h*int_h)+
                     0.25 * mass_lcap * (int_r+thick)*(int_r+thick) + mass_lcap * (0.5*(int_h + thick)) * (0.5*(int_h + thick))+
                     0.25 * mass_rcap * (int_r+thick)*(int_r+thick) + mass_rcap * (0.5*(int_h + thick)) * (0.5*(int_h + thick));
        inertia_yy = 0.5 * mass_hollow_cyl * (int_r*int_r + (int_r+thick)*(int_r+thick))+
                     0.5 * mass_lcap * (int_r+thick)*(int_r+thick)+
                     0.5 * mass_rcap * (int_r+thick)*(int_r+thick);
        inertia_zz = mass_hollow_cyl/12.0 * (3.0*(int_r*int_r + (int_r+thick)*(int_r+thick)) + int_h*int_h)+
                     0.25 * mass_lcap * (int_r+thick)*(int_r+thick) + mass_lcap * (0.5*(int_h + thick)) * (0.5*(int_h + thick))+
                     0.25 * mass_rcap * (int_r+thick)*(int_r+thick) + mass_rcap * (0.5*(int_h + thick)) * (0.5*(int_h + thick));
        
        vertices.push_back(std::vector<double>({0.0, -0.5*int_h, 0.0}));
        for (int ni=0; ni<p_res; ni++)
        {
            double theta = 2.0*M_PI*ni/p_res;
            vertices.push_back(std::vector<double>({int_r*cos(theta), -0.5*int_h, int_r*sin(theta)}));
            vertices.push_back(std::vector<double>({int_r*cos(theta), 0.5*int_h, int_r*sin(theta)}));
        }
        if (cap)
            vertices.push_back(std::vector<double>({0.0, 0.5*int_h, 0.0}));
        
        for (int ne=0; ne<p_res; ne++)
        {
            elems.push_back(std::vector<int>({0, 2*ne+1, ne==p_res-1 ? 1 : 2*ne+3}));
            elems.push_back(std::vector<int>({2*ne+1, 2*ne+2, ne==p_res-1 ? 1 : 2*ne+3}));
            elems.push_back(std::vector<int>({2*ne+2, ne==p_res-1 ? 1 : 2*ne+3, ne==p_res-1 ? 2 : 2*ne+4}));
            if (cap)
                elems.push_back(std::vector<int>({2*p_res+1, 2*ne+2, ne==p_res-1 ? 2 : 2*ne+4}));
        }
        
    }
    double int_r;
    double int_h;
    double thick;
    bool cap;
    int p_res = 32;
};

#endif