#ifndef SYSTEM_H
#define SYSTEM_H

#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/ChSystemSMC.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/solver/ChSolverPSOR.h"
#include "chrono/assets/ChTexture.h"
#include "chrono/physics/ChLinkMate.h"
#include "chrono_irrlicht/ChIrrApp.h"
#include "utils/ChUtilsCreators.h"
#include "utils/ChUtilsInputOutput.h"
#include "utils/ChUtilsGenerators.h"
#include "chrono_multicore/physics/ChSystemMulticore.h"
// #include "chrono_multicore/collision/ChCollisionSystemChronoMulticore.h"
#include "chrono/physics/ChLinkMotorRotationAngle.h"

#include "body.h"
#include <random>
#include <chrono>
#include <Eigen/Core>
#include <cmath>

using namespace chrono;
using namespace chrono::irrlicht;
using namespace postprocess;

// Use the main namespaces of Irrlicht
using namespace irr;
// using namespace irr::core;
using namespace irr::scene;
using namespace irr::video;
using namespace irr::io;
using namespace irr::gui;

class ClipSystem
{
    public:
    ClipSystem(ChSystemMulticoreNSC &_sys) : sys(&_sys){};
    virtual ~ClipSystem() = default;
    ChSystemMulticoreNSC * sys;
    std::set<int> gids;
    int n_body = 0;
    double envelope = 0.001;
    double margin = 0.0001;

    // std::map<std::string, Clip> clips;
    // std::map<std::string, Box> boxes;
    std::map<std::string, std::shared_ptr<Body>> bodies;
    
    std::string name_it(std::string body_type, int gid) {return body_type+'_'+std::to_string(gid); }

    int make_clip(double heigth, double width, double rad, double gap, double f_r, double density, double friction, double restitution, std::string name)
    {
        if (bodies.find(name) != bodies.end())
            std::cout << "clip name " << name << " already exists. It will be modified to the new properties provided here! \n";
        // Clip clp = Clip(heigth, width, rad, gap, f_r, density, friction, restitution);
        bodies.insert(std::make_pair(name, std::make_shared<Clip>(heigth, width, rad, gap, f_r, density, friction, restitution)) );
        return 1;
    }

    int make_cylinder(double int_r, double int_h, double thick, bool cap, double density, double friction, double restitution, std::string name)
    {
        if (bodies.find(name) != bodies.end())
            std::cout << "cylinder name " << name << " already exists. It will be modified to the new properties provided here! \n";
        // Clip clp = Clip(heigth, width, rad, gap, f_r, density, friction, restitution);
        bodies.insert(std::make_pair(name, std::make_shared<Cylinder>(int_r, int_h, thick, cap, density, friction, restitution)) );
        return 1;
    }

    int make_box(double int_x, double int_y, double int_z, double thick, double density, double friction, double restitution, std::string name)
    {
        if (bodies.find(name) != bodies.end())
            std::cout << "box name " << name << " already exists. It will be modified to the new properties provided here! \n";
        // Box bx = Box(int_x, int_y, int_z, thick, density, friction, restitution);
        bodies.insert(std::make_pair(name, std::make_shared<Box>(int_x, int_y, int_z, thick, density, friction, restitution)) );
        // bodies.insert(std::make_pair(name, bx) );
        return 1;
    }

    int add_floor(ChVector<> dim, ChVector<> center, int gid)
    {
        double density = 3000.0;
        auto floor = chrono_types::make_shared<ChBodyEasyBox>(dim.x(), dim.y(), dim.z(), density, false, false);
        if (gids.find(gid) != gids.end())
        {
            int new_gid = *gids.end() + 1;
            std::cout << "the body id " << gid << " is already exist. The code assigns " << new_gid << " to the new body!\n";
            gid = new_gid;
        }
        gids.insert(gid);
        floor->SetGid(gid);
        floor->SetPos(center);
        floor->SetCollide(false);
        floor->SetBodyFixed(true);
        std::string bname = name_it("floor", gid);
        floor->SetName(bname.c_str());
        sys->Add(floor);
        n_body++;
        return 1;
    }

    int add_box(std::string name, int gid, ChVector<>& center, const ChQuaternion<double>& rot)
    {
        // if (wall_t==0.0)
        // {
        //     wall_t = 0.1*std::min(std::min(dim.x(), dim.y()), dim.z());
        // }
        auto box_mat = chrono_types::make_shared<ChMaterialSurfaceNSC>();
        std::shared_ptr<Box> bx = std::dynamic_pointer_cast<Box>(bodies[name]);
        box_mat->SetFriction(bx->fric);
        box_mat->SetRestitution(bx->rest);

        auto box = std::shared_ptr<ChBody>(sys->NewBody());
        box->GetCollisionModel()->ClearModel();
        box->SetCollide(true);
        box->GetCollisionModel()->SetDefaultSuggestedEnvelope(bx->envelope);
        box->GetCollisionModel()->SetDefaultSuggestedMargin(bx->margin);
        // std::cout << dim << std::endl;
        // std::cout << "thickness : " <<  wall_t  << "\n dims : " << dim.x() << ", " << dim.y() << ", " << dim.z() <<  std::endl;
        double wall_t = bx->thick;
        double int_x = bx->int_x;
        double int_y = bx->int_y;
        double int_z = bx->int_z;
        // x-normal
        box->GetCollisionModel()->AddBox(box_mat, 0.5*wall_t, 0.5*int_y+wall_t, 0.5*int_z+wall_t, ChVector<>(0.5*(int_x+wall_t), 0.0, 0.0), QUNIT);
        box->GetCollisionModel()->AddBox(box_mat, 0.5*wall_t, 0.5*int_y+wall_t, 0.5*int_z+wall_t, ChVector<>(-0.5*(int_x+wall_t), 0.0, 0.0), QUNIT);
        // y-normal
        box->GetCollisionModel()->AddBox(box_mat, 0.5*int_x+wall_t, 0.5*wall_t, 0.5*int_z+wall_t, ChVector<>(0.0, 0.5*(int_y+wall_t), 0.0), QUNIT);
        box->GetCollisionModel()->AddBox(box_mat, 0.5*int_x+wall_t, 0.5*wall_t, 0.5*int_z+wall_t, ChVector<>(0.0, -0.5*(int_y+wall_t), 0.0), QUNIT);
        // z-normal
        box->GetCollisionModel()->AddBox(box_mat, 0.5*int_x+wall_t, 0.5*int_y+wall_t, 0.5*wall_t, ChVector<>(0.0, 0.0, 0.5*(int_z+wall_t)), QUNIT);
        box->GetCollisionModel()->AddBox(box_mat, 0.5*int_x+wall_t, 0.5*int_y+wall_t, 0.5*wall_t, ChVector<>(0.0, 0.0, -0.5*(int_z+wall_t)), QUNIT);
        
        // double volume_out = (int_x+2.*wall_t) * (int_y+2.*wall_t) * (int_z+2.*wall_t);
        // double volume_in = int_x * int_y * int_z;
        // double volume = volume_out - volume_in;
        // double mass = density * volume;
        // double mass_out = density * volume_out;
        // double mass_in = density * volume_in;
        // double inertia_xx = mass_out/12.0 * ((dim.y()+2.*wall_t)*(dim.y()+2.*wall_t) + (dim.z()+2.*wall_t)*(dim.z()+2.*wall_t))
        //                     - mass_in/12.0 * (dim.y()*dim.y() + dim.z()*dim.z());
        // double inertia_yy = mass_out/12.0 * ((int_x+2.*wall_t)*(int_x+2.*wall_t) + (dim.z()+2.*wall_t)*(dim.z()+2.*wall_t))
        //                     - mass_in/12.0 * (int_x*int_x + dim.z()*dim.z());
        // double inertia_zz = mass_out/12.0 * ((dim.y()+2.*wall_t)*(dim.y()+2.*wall_t) + (int_x+2.*wall_t)*(int_x+2.*wall_t))
        //                     - mass_in/12.0 * (dim.y()*dim.y() + dim.x()*dim.x());
        box->SetMass(bx->mass);
        box->SetInertia(ChMatrix33<>(ChVector<>(bx->inertia_xx, bx->inertia_yy, bx->inertia_zz)));
        box->SetPos(center);
        box->SetRot(rot);
        if (gids.find(gid) != gids.end())
        {
            int new_gid = *gids.rbegin() + 1;
            std::cout << "the body gid " << gid << " is already exist. The code assigns " << new_gid << " to the new body!\n";
            gid = new_gid;
        }
        gids.insert(gid);
        std::string bname = name_it("box", gid);
        box->SetName(bname.c_str());
        box->SetGid(gid);
        // box->SetBodyFixed(true);
        box->GetCollisionModel()->BuildModel();
        sys->Add(box);
        n_body++;
        return 1;
    }

    int add_cylinder(std::string name, int gid, ChVector<>& center, const ChQuaternion<double>& rot)
    {
        auto cylin_mat = chrono_types::make_shared<ChMaterialSurfaceNSC>();
        std::shared_ptr<Cylinder> cln = std::dynamic_pointer_cast<Cylinder>(bodies[name]);
        cylin_mat->SetFriction(cln->fric);
        cylin_mat->SetRestitution(cln->rest);

        auto cylin = std::shared_ptr<ChBody>(sys->NewBody());
        cylin->GetCollisionModel()->ClearModel();
        cylin->SetCollide(true);
        cylin->GetCollisionModel()->SetDefaultSuggestedEnvelope(cln->envelope);
        cylin->GetCollisionModel()->SetDefaultSuggestedMargin(cln->margin);
        // std::cout << dim << std::endl;
        // std::cout << "thickness : " <<  wall_t  << "\n dims : " << dim.x() << ", " << dim.y() << ", " << dim.z() <<  std::endl;
        double thick = cln->thick;
        double int_r = cln->int_r;
        double int_h = cln->int_h;
        
        // cylin->GetCollisionModel()->AddCylindricalShell(cylin_mat, int_r, 0.5*int_h, center, QUNIT);
        double sec_z = 2.0*M_PI * int_r / cln->p_res + thick;
        for (int ni=0; ni<cln->p_res; ni++)
        {
            double theta = 2.0*M_PI * (ni+0.5) / cln->p_res;
            ChVector<> sec_cen((int_r+0.5*thick)*cos(theta), 0.0, (int_r+0.5*thick)*sin(theta));
            ChQuaternion<> sec_q;
            sec_q.Q_from_AngAxis(-theta, ChVector<>(0.0, 1.0, 0.0));
            cylin->GetCollisionModel()->AddBox(cylin_mat, 0.5*thick, 0.5*int_h, 0.5*sec_z, sec_cen, sec_q);
        }
        cylin->GetCollisionModel()->AddCylinder(cylin_mat, int_r, int_r, 0.5*thick, ChVector<>(0.0, -0.5*(int_h+thick), 0.0), QUNIT);
        if (cln->cap)
            cylin->GetCollisionModel()->AddCylinder(cylin_mat, int_r, int_r, 0.5*thick, ChVector<>(0.0, +0.5*(int_h+thick), 0.0), QUNIT);
        
        cylin->SetMass(cln->mass);
        cylin->SetInertia(ChMatrix33<>(ChVector<>(cln->inertia_xx, cln->inertia_yy, cln->inertia_zz)));
        cylin->SetPos(center);
        cylin->SetRot(rot);
        
        if (gids.find(gid) != gids.end())
        {
            int new_gid = *gids.rbegin() + 1;
            std::cout << "the body gid " << gid << " is already exist. The code assigns " << new_gid << " to the new body!\n";
            gid = new_gid;
        }
        gids.insert(gid);
        std::string bname = name_it("cylin", gid);
        cylin->SetName(bname.c_str());
        cylin->SetGid(gid);
        // box->SetBodyFixed(true);
        cylin->GetCollisionModel()->BuildModel();
        sys->Add(cylin);
        n_body++;
        return 1;
    }

    int add_clip(std::string name, int gid, ChVector<>& center, const ChQuaternion<double>& rot)
    {
        auto clip_mat = chrono_types::make_shared<ChMaterialSurfaceNSC>();
        std::shared_ptr<Clip> clp = std::dynamic_pointer_cast<Clip>(bodies[name]);
        clip_mat->SetFriction(clp->fric);
        clip_mat->SetRestitution(clp->rest);

        auto clip = std::shared_ptr<ChBody>(sys->NewBody());
        
        clip->SetCollide(true);
        clip->GetCollisionModel()->ClearModel();
        clip->GetCollisionModel()->SetEnvelope(0.001);
        clip->GetCollisionModel()->SetDefaultSuggestedMargin(0.0001);
        // std::cout << "num of edges si : " << bodies[name].elems.size() << std::endl;
        for (std::vector<int> edge : clp->elems)
        {
            ChVector<> n1(clp->vertices[edge[0]].at(0), clp->vertices[edge[0]].at(1), clp->vertices[edge[0]].at(2));
            ChVector<> n2(clp->vertices[edge[1]].at(0), clp->vertices[edge[1]].at(1), clp->vertices[edge[1]].at(2));
            double h = 0.5*(n2 - n1).Length();
            ChVector<> d = 0.5 * (n1 + n2);
            ChVector<> axis_angle;
            axis_angle.Cross((n2-n1).GetNormalized(), ChVector<>(0.0, 1.0, 0.0));
            // std::cout << "axis angle : " << axis_angle << std::endl;
            ChQuaternion<> q_rot;
            q_rot.Q_from_Rotv(axis_angle);
            // std::cout << "rotation is : " << q_rot.Q_to_Rotv() << std::endl;
            utils::AddCylinderGeometry(clip.get(), clip_mat, clp->r, h, d, q_rot, true);
        }

        ChVector<double> inertia = rot.Rotate(ChVector<>(clp->inertia_xx, clp->inertia_yy, clp->inertia_zz));

        clip->GetCollisionModel()->BuildModel();
        // std::cout << "clip mass is : " << clips[name].mass << std::endl;
        clip->SetMass((double)clp->mass);
        clip->SetInertiaXX(inertia);
        clip->SetPos(center);
        clip->SetRot(rot);
        // std::cout << "rotation is : " << rot << std::endl;
        if (gids.find(gid) != gids.end())
        {
            int new_gid = *gids.rbegin() + 1;
            std::cout << "the body gid " << gid << " is already exist. The code assigns " << new_gid << " to the new body!\n";
            gid = new_gid;
        }
        // std::cout << "gid is : " << gid << " and find " << *gids.find(gid) << std::endl;
        gids.insert(gid);
        clip->SetGid(gid);
        std::string bname = name_it(name, gid);
        clip->SetName(bname.c_str());
        // std::cout << "clip rot : " << clip->GetRot() << ", mass : " << clip->GetMass() << ", inertia : " << clip->GetInertia() << std::endl;
        sys->AddBody(clip);
        n_body++;
        return 1;
    }

    int add_hook(double h_rad, double rad, double gap, double density, int gid, 
                ChVector<>& center, const ChQuaternion<double>& rot=QUNIT, double friction=0.1, double restitution=0.2)
    {
        double nseg = 32;
        auto hook_mat = chrono_types::make_shared<ChMaterialSurfaceNSC>();
        hook_mat->SetFriction(friction);
        hook_mat->SetRestitution(restitution);

        auto hook = std::shared_ptr<ChBody>(sys->NewBody());
        
        hook->SetCollide(true);
        hook->GetCollisionModel()->ClearModel();
        hook->GetCollisionModel()->SetEnvelope(envelope);
        hook->GetCollisionModel()->SetDefaultSuggestedMargin(margin);
        ChQuaternion<double> gap_rot = Q_from_AngAxis(0.5*gap, ChVector<double>(0, 0, 1.0));
        utils::AddTorusGeometry(hook.get(), hook_mat, h_rad, rad, nseg, 360.0 - gap, ChVector<>(0, 0, 0), gap_rot,true);

        double volume = (2.0*M_PI-gap) * h_rad * M_PI * rad * rad;
        double mass = density * volume;
        // values should be updated...
        double inertia_d = mass/8.0 * (4.*h_rad*h_rad + 5.*rad*rad);
        double inertia_zz = mass/4.0 * (4.*h_rad*h_rad + 3.*rad*rad);
        ChVector<> inertia(inertia_d, inertia_d, inertia_zz);

        hook->GetCollisionModel()->BuildModel();
        hook->SetMass(mass);
        hook->SetInertiaXX(inertia);
        hook->SetPos(center);
        hook->SetRot(rot);
        if (gids.find(gid) != gids.end())
        {
            int new_gid = *gids.rbegin() + 1;
            std::cout << "the body gid " << gid << " is already exist. The code assigns " << new_gid << " to the new body!\n";
            gid = new_gid;
        }
        gids.insert(gid);
        hook->SetGid(gid);
        std::string bname = name_it("hook", gid);
        hook->SetName(bname.c_str());
        sys->AddBody(hook);
        n_body++;
        return 1;
    }

    int clip_in_box(ChVector<> dim, ChVector<> center, int num_clip, double min_offset, std::string name)
    {
        ChVector<> offset(min_offset, min_offset, min_offset);
        std::shared_ptr<Clip> clp = std::dynamic_pointer_cast<Clip>(bodies[name]);
        ChVector<> clip_boounding_box(clp->w, clp->h, 2.0*clp->r);
        ChVector<> clip_box = clip_boounding_box + offset;
        int max_nz;
        int max_ny;
        int max_nx;
        // std::cout << (dim.z() - 2.0*clips[name].r - wall_offset.z()) / clip_box.z() << std::endl;
        // std::cout << (dim.y() - clips[name].w - wall_offset.y()) / clip_box.y() << std::endl;
        // std::cout << (dim.x() - clips[name].h - wall_offset.x()) / clip_box.x() << std::endl;
        max_nz = (int) std::floor(dim.z() / clip_box.z()) - 2;
        max_ny = (int) std::floor(dim.y() / clip_box.y()) - 2;
        max_nx = (int) std::floor(dim.x() / clip_box.x()) - 2;
        if (num_clip > max_nz)
        {
            if(num_clip > max_nz * max_ny)
            {
                max_nx = num_clip % (max_nz*max_ny) == 0 ? (int) (num_clip / (max_nz*max_ny)) : (int) (num_clip / (max_nz*max_ny) + 1);
            }
            else
            {
                max_ny = num_clip % max_nz == 0 ? (int) (num_clip / max_nz) : (int) (num_clip / max_nz + 1);
                max_nx = 1;
            }
        }
        else
        {
            max_nz = num_clip;
            max_ny = 1;
            max_nx = 1;
        }
        
        std::cout << "num_clip: " << num_clip << "\n N_z : " << max_nz << "\n N_y : " << max_ny << "\n N_x : " << max_nx << std::endl;
        ChVector<> blb_corner = center - 0.5 * dim + offset + clip_box;
        int gid = *gids.rbegin() + 1;
        for (int nc=0; nc< num_clip; nc++)
        {
            int nx = nc / (max_nz * max_ny);
            int ny = (nc - nx * max_nz * max_ny) / max_nz;
            int nz = nc - nx * max_nz * max_ny - ny * max_nz;
            ChVector<> clip_center = blb_corner + ChVector<>(nx*clip_box.x(), ny*clip_box.y(), nz*clip_box.z());
            add_clip(name, gid+nc, clip_center, QUNIT);
            // std::cout << clip_center << std::endl;
        }

        return 1;
    }

    std::vector<std::vector<double>> get_updated_location(std::string body_name, int gid)
    {
        std::vector<std::vector<double>> body_verts;
        // body_verts.push_back({0.0, 0.0, 0.0});

        // for (auto body : sys->Get_bodylist())
        // {
        //     std::string bname = body->GetName();
        //     int gid = body->GetGid();
        //     int iid = body->GetId();
        //     // std::cout << "name : " << bname << ", GID : " << gid << ", ID : " << iid << std::endl;
        //     ChVector<> translation = body->GetPos();
        //     ChQuaternion<> rotation = body->GetRot();
        //     // std::cout << "translation : " << translation << std::endl;
        //     // std::cout << "rotation : " << rotation << std::endl;
        // }



        std::string bname = name_it(body_name, gid);
        auto body = sys->SearchBody(bname.c_str());
        // std::cout << "name : " << body->GetName() << ", GID : " << body->GetGid() << ", ID : " << body->GetId() << std::endl;
        // // ChMatrix33<> rotationZX;
        // // rotationZX << 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0;
        // std::cout << "body idi is " << id << std::endl;
        // auto body = sys->SearchBodyID(id);
        // std::cout << "body name : " << body->GetGid() << std::endl;
        ChVector<> translation = body->GetPos();
        ChMatrix33<> rotation = body->GetRot();
        for (auto stdv : bodies[body_name]->vertices)
        {
            ChVector<> gvert = rotation * ChVector<>(stdv[0], stdv[1], stdv[2]) + translation;
            body_verts.push_back({gvert.x(), gvert.y(), gvert.z()});
            // std::cout << "clip verts : " << stdv[0] << ", " << stdv[1] << ", " << stdv[2] << std::endl;
            // std::cout << "clip verts : " << body_verts.back()[0] << ", " << body_verts.back()[1] << ", " << body_verts.back()[2] << std::endl;
        }
        // std::cout << "translation : " << translation << std::endl;
        // std::cout << "rotation : " << rotation << std::endl;
        // std::cout << "3D vecrtor : " << clips[body_name].vertices.at(3).at(2) << std::endl;

        return body_verts;
    }

    void add_motion(std::string bname_m, std::string bname_f, std::string type, ChVector<> phase, ChVector<> freq, ChVector<> amp)
    {
        auto moving = sys->SearchBody(bname_m.c_str());
        auto fixed = sys->SearchBody(bname_f.c_str());

        

        if (type=="vibration")
        {
            auto motion = chrono_types::make_shared<ChLinkLockLock>();
            motion->Initialize(moving , fixed, ChCoordsys<>(ChVector<>(0, 0, 0)));
            auto mmotion_x = chrono_types::make_shared<ChFunction_Sine>(phase.x(), freq.x(), amp.x());  // phase freq ampl
            motion->SetMotion_X(mmotion_x);
            auto mmotion_y = chrono_types::make_shared<ChFunction_Sine>(phase.y(), freq.y(), amp.y());  // phase freq ampl
            motion->SetMotion_Y(mmotion_y);
            auto mmotion_z = chrono_types::make_shared<ChFunction_Sine>(phase.z(), freq.z(), amp.z());  // phase freq ampl
            motion->SetMotion_Z(mmotion_z);
            sys->Add(motion);
        }
        else if (type=="rotation")
        {
            // Create a motor between the two bodies, constrained to rotate at 90 deg/s
            auto motor = chrono_types::make_shared<ChLinkMotorRotationAngle>();
            motor->Initialize(moving, fixed, ChFrame<>(ChVector<>(0, 0, 0), Q_ROTATE_X_TO_Z));
            motor->SetAngleFunction(chrono_types::make_shared<ChFunction_Ramp>(0, 10.0));
            sys->AddLink(motor);
        }
    }
    
};
#endif