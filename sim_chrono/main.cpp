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

#include "include/visualize.h"

// #define PICOPARAVIEW_IMPLEMENTATION
// #include "include/picoparaview.h"
#include "include/system.h"

#include <random>
#include <chrono>
#include <Eigen/Core>

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



void create_box(ChSystemMulticoreNSC *sys, std::shared_ptr<ChBody> floor, std::shared_ptr<ChMaterialSurface> mat,
                const int id, ChVector<>& pos, const ChQuaternion<double>& rot,
                const double box_x, const double box_y, const double box_z, const double box_t,
                const double density)
{
    // auto texture = chrono_types::make_shared<ChTexture>();
    // texture->SetTextureFilename(GetChronoDataFile("textures/concrete.jpg"));
    // auto box = std::make_shared<ChBody>();
    auto box = std::shared_ptr<ChBody>(sys->NewBody());
	
    
    box->GetCollisionModel()->ClearModel();
    box->SetCollide(true);
    box->GetCollisionModel()->SetDefaultSuggestedEnvelope(0.001);
    box->GetCollisionModel()->SetDefaultSuggestedMargin(0.001);
    double alpha = 1.1;
    box->GetCollisionModel()->AddBox(mat, 0.5*box_x+box_t, 0.5*box_t, 0.5*box_z+box_t, ChVector<>(0.0, 0.5*(box_y+box_t), 0), QUNIT);
    box->GetCollisionModel()->AddBox(mat, 0.5*box_x+box_t, 0.5*box_t, 0.5*box_z+box_t, ChVector<>(0.0, -0.5*(box_y+box_t), 0), QUNIT);
    box->GetCollisionModel()->AddBox(mat, 0.5*box_x+box_t, 0.5*box_y+box_t, 0.5*box_t, ChVector<>(0.0, 0.0, 0.5*(box_z+box_t)), QUNIT);
    box->GetCollisionModel()->AddBox(mat, 0.5*box_x+box_t, 0.5*box_y+box_t, 0.5*box_t, ChVector<>(0.0, 0.0, -0.5*(box_z+box_t)), QUNIT);
    box->GetCollisionModel()->AddBox(mat, 0.5*box_t, 0.5*box_y+box_t, 0.5*box_z+box_t, ChVector<>(0.5*(box_x+box_t), 0.0, 0.0), QUNIT);
    box->GetCollisionModel()->AddBox(mat, 0.5*box_t, 0.5*box_y+box_t, 0.5*box_z+box_t, ChVector<>(-0.5*(box_x+box_t), 0.0, 0.0), QUNIT);
    
    double volume_bottop = (box_x+2.*box_t) * (box_z+2.*box_t) * box_t;
    double volume_lefright = (box_y+2.*box_t) * (box_z+2.*box_t) * box_t;
    double volume_fronback = (box_x+2.*box_t) * (box_y+2.*box_t) * box_t;
    double volume = 2. * (volume_fronback + volume_lefright + volume_bottop);
    double mass_bottop = density * volume_bottop;
    double mass_lefright = density * volume_lefright;
    double mass_fronback = density * volume_fronback;
    double mass = 2. * (mass_fronback + mass_lefright + mass_bottop);
    double iner_xyz = mass/12.0 * ((box_y+2.*box_t)*(box_y+2.*box_t) + (box_z+2.*box_t)*(box_z+2.*box_t));
    // double iner_y = mass/12.0 * (box_x*box_x + box_z*box_z);
    // double iner_z = mass/12.0 * (box_x*box_x + box_y*box_y);
    auto box_iner = ChMatrix33<>(ChVector<>(iner_xyz, iner_xyz, iner_xyz));
    box->SetMass(mass);
    box->SetInertia(box_iner);
    // box->SetMaterialSurface(mat);
    // std::cout << "box materila is : " << box->
    box->SetPos(pos);
    box->SetRot(rot);
    box->GetCollisionModel()->BuildModel();
    box->SetName("box");
    sys->Add(box);
    box->SetGid(id);
    
    // box->AddAsset(texture);



    // auto bot_viz = chrono_types::make_shared<ChBoxShape>();
    // bot_viz->GetBoxGeometry().SetLengths(ChVector<>(alpha*box_x, box_t, alpha*box_z));
    // bot_viz->Pos = ChVector<>(0.0, -0.5*(box_t+box_y), 0);
    // box->GetAssets().push_back(bot_viz);

    // auto top_viz = chrono_types::make_shared<ChBoxShape>();
    // top_viz->GetBoxGeometry().SetLengths(ChVector<>(alpha*box_x, box_t, alpha*box_z));
    // top_viz->Pos = ChVector<>(0.0, 0.5*(box_t+box_y), 0);
    // box->GetAssets().push_back(top_viz);

    // auto left_viz = chrono_types::make_shared<ChBoxShape>();
    // left_viz->GetBoxGeometry().SetLengths(ChVector<>(box_t, alpha*box_y, alpha*box_z));
    // left_viz->Pos = ChVector<>(-0.5*(box_t+box_x), 0.0, 0.0);
    // box->GetAssets().push_back(left_viz);

    // auto right_viz = chrono_types::make_shared<ChBoxShape>();
    // right_viz->GetBoxGeometry().SetLengths(ChVector<>(box_t, alpha*box_y, alpha*box_z));
    // right_viz->Pos = ChVector<>(0.5*(box_t+box_x), 0.0, 0.0);
    // box->GetAssets().push_back(right_viz);

    // auto front_viz = chrono_types::make_shared<ChBoxShape>();
    // front_viz->GetBoxGeometry().SetLengths(ChVector<>(alpha*box_x, alpha*box_y, box_t));
    // front_viz->Pos = ChVector<>(0.0, 0.0, 0.5*(box_t+box_z));
    // box->GetAssets().push_back(front_viz);

    // auto back_viz = chrono_types::make_shared<ChBoxShape>();
    // back_viz->GetBoxGeometry().SetLengths(ChVector<>(alpha*box_x, alpha*box_y, box_t));
    // back_viz->Pos = ChVector<>(0.0, 0.0, -0.5*(box_t+box_z));
    // box->GetAssets().push_back(back_viz);
    



    auto box_motion = chrono_types::make_shared<ChLinkLockLock>();
    box_motion->Initialize(box , floor, ChCoordsys<>(ChVector<>(0, 0, 0)));

    auto mmotion_x = chrono_types::make_shared<ChFunction_Sine>(0, 4.0, 0.005);  // phase freq ampl
    box_motion->SetMotion_X(mmotion_x);
    auto mmotion_y = chrono_types::make_shared<ChFunction_Sine>(0, 4.0, 0.01);  // phase freq ampl
    box_motion->SetMotion_Y(mmotion_y);
    auto mmotion_z = chrono_types::make_shared<ChFunction_Sine>(0, 4.0, 0.005);  // phase freq ampl
    box_motion->SetMotion_Z(mmotion_z);
    sys->Add(box_motion);

}

void create_clip(ChSystemMulticoreNSC *sys, std::shared_ptr<ChMaterialSurface> mat,
                const int id, ChVector<>& pos, const ChQuaternion<double>& rot,
                const double clip_w, const double clip_h, const double clip_r, const double clip_g,
                const double mass, const ChVector<> inertia)
{
    // std::random_device rd;
    // std::default_random_engine eng(rd());
    // std::uniform_real_distribution<float> distr(0.0, 1.0);

    // // auto clip_color = chrono_types::make_shared<ChColorAsset>(ChColor(0.8f, 0.1f, 0.1f));
    // // auto init_ang_velo = (10.0*M_PI) * ChVector<>(distr(eng), distr(eng), distr(eng));
    // // auto init_velo = (0.2) * ChVector<>(distr(eng), distr(eng), distr(eng));
    // int velo_scale = id>0 ? 1.0 : id<0 ? -1.0 : 0;
    // // std::cout << "velo scale : " << velo_scale << std::endl;
    // auto init_ang_velo = (5.0*M_PI*velo_scale) * ChVector<>(1.0, 0.0, 0.0);
    // auto init_velo = (0.05*velo_scale) * ChVector<>(-1.0, -1.0, 0.0);
    // // if (id % 2 != 0)
    // // {
    // //     init_ang_velo = -init_ang_velo;
    // //     init_velo = -init_velo;
    // // }
    // // auto comp_pos = ChVector<>(-0.5*clip_h, 0, 0);
    // // ChQuaternion<>(1, 0, 0, 0)
    // // if (id % 2 != 0)
    // // {
    // //     // init_ang_velo = ChVector<>((20.0*M_PI), 0.0, 10.0*M_PI);
    // //     // inertia = chVector<> (inertia(0), inertia(2), inertia(1));
    // //     clip_color = chrono_types::make_shared<ChColorAsset>(ChColor(0.1f, 0.1f, 0.8f));
    // // }
    
    // // auto ball = std::shared_ptr<ChBody>(sys->NewBody());
    // auto clip = std::shared_ptr<ChBody>(sys->NewBody());
	// clip->SetGid(id);
    // clip->SetCollide(true);
    // clip->GetCollisionModel()->ClearModel();
    // // clip->GetCollisionModel()->SetDefaultSuggestedEnvelope(0.001);
    // clip->GetCollisionModel()->SetEnvelope(0.03);
    // clip->GetCollisionModel()->SetDefaultSuggestedMargin(0.0001);
    
    
    // // clip->SetMaterialSurface(mat);

    // double added_l = 1.0*clip_r;

    // utils::AddCylinderGeometry(clip.get(),
    //                             mat,
    //                             clip_r,
    //                             0.5*clip_w,
    //                             ChVector<>(-0.5*clip_h, 0, 0),
    //                             ChQuaternion<>(1, 0, 0, 0),
    //                             true);
    // utils::AddCylinderGeometry(clip.get(),
    //                             mat,
    //                             clip_r,
    //                             0.5*clip_w,
    //                             ChVector<>(0.5*clip_h, 0, 0),
    //                             ChQuaternion<>(1, 0, 0, 0),
    //                             true);
    // utils::AddCylinderGeometry(clip.get(),
    //                             mat,
    //                             clip_r,
    //                             0.5*clip_h,
    //                             ChVector<>(0.0, 0.5*clip_w, 0),
    //                             Q_ROTATE_Y_TO_X,
    //                             true);
    // utils::AddCylinderGeometry(clip.get(),
    //                             mat,
    //                             clip_r,
    //                             0.25*(clip_h-clip_g),
    //                             ChVector<>(-0.25*(clip_h+clip_g), -0.5*clip_w, 0),
    //                             Q_ROTATE_Y_TO_X,
    //                             true);
    // utils::AddCylinderGeometry(clip.get(),
    //                             mat,
    //                             clip_r,
    //                             0.25*(clip_h-clip_g),
    //                             ChVector<>(0.25*(clip_h+clip_g), -0.5*clip_w, 0),
    //                             Q_ROTATE_Y_TO_X,
    //                             true);
    
    // clip->GetCollisionModel()->BuildModel();

    // clip->SetMass(mass);
    // clip->SetInertiaXX(inertia);
    // clip->SetWvel_par(init_ang_velo);
    // clip->SetPos_dt(init_velo);
    
    // clip->SetPos(pos);
    // clip->SetRot(rot);
    // sys->AddBody(clip);
    // // clip->AddAsset(clip_color);
    // // if (id % 2 == 0)
    // //     clip->AddAsset(chrono_types::make_shared<ChColorAsset>(ChColor(0.8f, 0.1f, 0.1f)));
    // // else
    // //     clip->AddAsset(chrono_types::make_shared<ChColorAsset>(ChColor(0.1f, 0.1f, 0.8f)));
}


void create_closedclip(ChSystemMulticoreNSC *sys, std::shared_ptr<ChMaterialSurface> mat,
                const int id, ChVector<>& pos, const ChQuaternion<double>& rot,
                const double clip_w, const double clip_h, const double clip_r,
                const double mass, const ChVector<> inertia)
{
    std::random_device rd;
    std::default_random_engine eng(rd());
    std::uniform_real_distribution<float> distr(0.0, 1.0);

    // auto clip_color = chrono_types::make_shared<ChColorAsset>(ChColor(0.8f, 0.1f, 0.1f));
    // int velo_scale = id>0 ? 1.0 : id<0 ? -1.0 : 0;
    // auto init_ang_velo = (5.0*M_PI*velo_scale) * ChVector<>(1.0, 0.0, 0.0);
    // auto init_velo = (0.05*velo_scale) * ChVector<>(-1.0, -1.0, 0.0);
    // auto init_ang_velo = (5.0*M_PI) * ChVector<>(distr(eng), distr(eng), distr(eng));
    // auto init_velo = (0.05) * ChVector<>(distr(eng), distr(eng), distr(eng));
    // if (id % 2 != 0)
    // {
    //     init_ang_velo = -init_ang_velo;
    //     init_velo = -init_velo;
    // }
    // auto comp_pos = ChVector<>(-0.5*clip_h, 0, 0);
    // ChQuaternion<>(1, 0, 0, 0)
    // if (id % 2 != 0)
    // {
    //     // init_ang_velo = ChVector<>((20.0*M_PI), 0.0, 10.0*M_PI);
    //     // inertia = chVector<> (inertia(0), inertia(2), inertia(1));
    //     clip_color = chrono_types::make_shared<ChColorAsset>(ChColor(0.1f, 0.1f, 0.8f));
    // }
    
    // auto ball = std::shared_ptr<ChBody>(sys->NewBody());
    auto clip = std::shared_ptr<ChBody>(sys->NewBody());
	
    clip->SetCollide(true);
    clip->GetCollisionModel()->ClearModel();
    // clip->GetCollisionModel()->SetDefaultSuggestedEnvelope(0.001);
    clip->GetCollisionModel()->SetEnvelope(0.01);
    clip->GetCollisionModel()->SetDefaultSuggestedMargin(0.0001);
    
    
    // clip->SetMaterialSurface(mat);

    double added_l = 1.078*clip_r;

    utils::AddCylinderGeometry(clip.get(),
                                mat,
                                clip_r,
                                0.5*clip_w+added_l,
                                ChVector<>(-0.5*clip_h, 0, 0),
                                ChQuaternion<>(1, 0, 0, 0),
                                true);
    utils::AddCylinderGeometry(clip.get(),
                                mat,
                                clip_r,
                                0.5*clip_w+added_l,
                                ChVector<>(0.5*clip_h, 0, 0),
                                ChQuaternion<>(1, 0, 0, 0),
                                true);
    utils::AddCylinderGeometry(clip.get(),
                                mat,
                                clip_r,
                                0.5*clip_h+added_l,
                                ChVector<>(0.0, 0.5*clip_w, 0),
                                Q_ROTATE_Y_TO_X,
                                true);
    utils::AddCylinderGeometry(clip.get(),
                                mat,
                                clip_r,
                                0.5*clip_h+added_l,
                                ChVector<>(0.0, -0.5*clip_w, 0),
                                Q_ROTATE_Y_TO_X,
                                true);
    
    clip->GetCollisionModel()->BuildModel();

    clip->SetMass(mass);
    clip->SetInertiaXX(inertia);
    clip->SetPos(pos);
    clip->SetRot(rot);
    // ChQuaternion<> Q_ang_vel;
    // Q_ang_vel.Q_from_Rotv(init_ang_velo);
    // clip->SetRot_dt(ChQuaternion<>(0.0, 2.5, 1.25, 0.0) );
    // clip->SetWvel_loc(ChVector<>(init_ang_velo));
    // clip->SetPos_dt(init_velo);
    clip->SetGid(id);
    std::string bname = "clip" + std::to_string(id);
    clip->SetName(bname.c_str());
    // if (id==0)
    // {
    //     clip->SetName("clip0");
    //     // clip->SetBodyFixed(true);
    // }
    sys->AddBody(clip);
    
    // clip->AddAsset(clip_color);
    // if (id % 2 == 0)
    //     clip->AddAsset(chrono_types::make_shared<ChColorAsset>(ChColor(0.8f, 0.1f, 0.1f)));
    // else
    //     clip->AddAsset(chrono_types::make_shared<ChColorAsset>(ChColor(0.1f, 0.1f, 0.8f)));
    
}

void create_torus(ChSystemMulticoreNSC& sys, std::shared_ptr<ChMaterialSurface> mat,
                const int id, ChVector<>& pos, const ChQuaternion<double>& rot,
                const double torus_r, const double torus_t, const int nseg, const double angle,
                const double mass, const ChVector<> inertia)
{


    // auto torus_color = chrono_types::make_shared<ChColorAsset>(ChColor(0.8f, 0.1f, 0.1f));
    auto init_ang_velo = ChVector<>(0.0, (20.0*M_PI), 0.0);
    if (id % 2 != 0)
    {
        init_ang_velo = ChVector<>((20.0*M_PI), 0.0, 0.0*M_PI);
        // torus_color = chrono_types::make_shared<ChColorAsset>(ChColor(0.1f, 0.1f, 0.8f));
    }
    

    auto torus = std::make_shared<ChBody>();
	torus->SetGid(id);
    // clip->GetCollisionModel()->SetDefaultSuggestedEnvelope(0.001);
    torus->GetCollisionModel()->SetEnvelope(0.001);
    torus->GetCollisionModel()->SetDefaultSuggestedMargin(0.0001);
    torus->SetCollide(true);
    torus->GetCollisionModel()->ClearModel();

     utils::AddTorusGeometry(torus.get(),
                        mat,
                        torus_r,
                        torus_t,
                        nseg, 
                        angle,
                        ChVector<>(0, 0, 0),
                        ChQuaternion<>(1, 0, 0, 0),
                        true);
    
    torus->GetCollisionModel()->BuildModel();

    torus->SetMass(mass);
    torus->SetInertiaXX(inertia);
    
    // torus->SetWvel_par(init_ang_velo);
    torus->SetGid(id);
    torus->SetPos(pos);
    torus->SetRot(rot);
    std::string bname = "clip" + std::to_string(id);
    torus->SetName(bname.c_str());
    sys.Add(torus);
    // torus->AddAsset(torus_color);
}

void create_chain_torus(ChSystemMulticoreNSC& mphysicalSystem, int n_clips = 10, std::string clip_type = "closed") 
{
    ChVector<> gravity(0, -9.81, 0);
    mphysicalSystem.Set_G_acc(gravity);

    auto floorBody = chrono_types::make_shared<ChBodyEasyBox>(1., 0.2, 1.0, 7800, false, false);
    floorBody->SetGid(10001);
    floorBody->SetPos(ChVector<>(0, -0.3, 0));
    floorBody->SetBodyFixed(true);
    mphysicalSystem.Add(floorBody);

    auto clip_mat = chrono_types::make_shared<ChMaterialSurfaceNSC>();
    clip_mat->SetFriction(0.05f);
    // clip_mat->SetRollingFriction(0.001);
    // clip_mat->SetSpinningFriction(0.01);
    // clip_mat->SetCompliance(0.0001);
    // clip_mat->SetComplianceT(0.02);
    // clip_mat->SetComplianceRolling(0.08);
    // clip_mat->SetComplianceSpinning(0.08);
    // clip_mat->SetDampingF(0.001);
    clip_mat->SetRestitution(0.8f);

    // double clip_w = 0.015;
    // double clip_h = 0.04;
    // double clip_r = 0.001;
    // double clip_g = 0.004;
    // double closed_mass = 0.002695;
    // double open_mass = 0.0026;
    // auto cclip_iner_odd = ChVector<>(0.125176E-6, 0.556291236E-6, 0.680161092E-6);
    // auto cclip_iner_even = ChVector<>(0.125176E-6, 0.680161092E-6, 0.556291236E-6);
    // auto rot_odd = Q_ROTATE_X_TO_Y; //Q_ROTATE_X_TO_Y; //ChQuaternion<double>(0.7565, 0, 0.6, 0.8);
    // auto rot_even = Q_ROTATE_X_TO_Y * Q_ROTATE_Y_TO_Z; //Q_ROTATE_X_TO_Y; //ChQuaternion<double>(0.7565, 0, 0.6, 0.8);


    // Adding torus///////////////////////////
    double torus_r = 0.01;
    double torus_t = 0.001;
    int torus_nseg = 20;
    double torus_angle = 360.;
    double torus_mass = 0.003;
    auto torus_iner = ChVector<>(0.12E-6, 0.7E-6, 0.7E-6);
    auto rot_odd = Q_ROTATE_X_TO_Y; //Q_ROTATE_X_TO_Y; //ChQuaternion<double>(0.7565, 0, 0.6, 0.8);
    auto rot_even = Q_ROTATE_X_TO_Y * Q_ROTATE_Y_TO_Z; //Q_ROTATE_X_TO_Y; //ChQuaternion<double>(0.7565, 0, 0.6, 0.8);

    // auto pos = ChVector<>(0, 0, 0);
    // auto rot_yz = Q_ROTATE_Y_TO_Z;

    // create_torus(mphysicalSystem, clip_mat, 1, 
    //         pos, rot_odd,
    //         torus_r, torus_t, torus_nseg, torus_angle,
    //         torus_mass, torus_iner);




    double height = n_clips * (torus_r-torus_t);
    ChVector<> pos = ChVector<>(0.0, height, 0.0);
    for (int nc=0; nc<n_clips; ++nc)
    {
        // std::cout << "making clip : " << nc << " \n";
        // double mass;
        ChQuaternion<double> rot;
        // ChVector<> inertia;
        if (nc % 2 == 0)
        {
            // inertia = cclip_iner_even;
            rot = rot_even;
        }
        else
        {
            // inertia = cclip_iner_odd;
            rot = rot_odd;
        }
        create_torus(mphysicalSystem, clip_mat, nc, 
            pos, rot,
            torus_r, torus_t, torus_nseg, torus_angle,
            torus_mass, torus_iner);
        // if (clip_type=="closed")
        // {
        //     // mass = closed_mass;
        //     create_closedclip(&mphysicalSystem, clip_mat, nc, 
        //         pos, rot, clip_w, clip_h, clip_r,
        //         mass, inertia);
        // }
        // else if (clip_type=="open")
        // {
        //     mass = open_mass;
        //     // std::cout << "right before calling clip \n";
        //     create_clip(&mphysicalSystem, clip_mat, nc, 
        //         pos, rot, clip_w, clip_h, clip_r, clip_g,
        //         mass, inertia);
        // }
        // else
        //     std::cout << "Wrong clip type, it should be either open or closed!.....\n";
        pos -= ChVector<>(0.0, (torus_r-torus_t), 0.0);
    }

    auto first_clip = mphysicalSystem.SearchBody("clip0");
    auto clip_motion = chrono_types::make_shared<ChLinkLockLock>();
    clip_motion->Initialize(first_clip , floorBody, ChCoordsys<>(ChVector<>(0, 0, 0)));
    // auto mmotion_x = chrono_types::make_shared<ChFunction_Sine>(0, 10.0, 0.005);  // phase freq ampl
    // clip_motion->SetMotion_X(mmotion_x);
    // auto mmotion_y = chrono_types::make_shared<ChFunction_Sine>(0, 8.0, 0.01);  // phase freq ampl
    // clip_motion->SetMotion_Y(mmotion_y);
    // auto mmotion_z = chrono_types::make_shared<ChFunction_Sine>(0, 5.0, 0.005);  // phase freq ampl
    // clip_motion->SetMotion_Z(mmotion_z);
    mphysicalSystem.Add(clip_motion);
}

void create_chain(ChSystemMulticoreNSC& mphysicalSystem, int n_clips = 10, std::string clip_type = "closed") {
    ChVector<> gravity(0, -9.81, 0);
    mphysicalSystem.Set_G_acc(gravity);

    auto floorBody = chrono_types::make_shared<ChBodyEasyBox>(1., 0.2, 1.0, 7800, false, false);
    floorBody->SetGid(10001);
    floorBody->SetPos(ChVector<>(0, -0.3, 0));
    floorBody->SetBodyFixed(true);
    mphysicalSystem.Add(floorBody);

    auto clip_mat = chrono_types::make_shared<ChMaterialSurfaceNSC>();
    clip_mat->SetFriction(0.05f);
    // clip_mat->SetRollingFriction(0.001);
    // clip_mat->SetSpinningFriction(0.01);
    // clip_mat->SetCompliance(0.0001);
    // clip_mat->SetComplianceT(0.02);
    // clip_mat->SetComplianceRolling(0.08);
    // clip_mat->SetComplianceSpinning(0.08);
    // clip_mat->SetDampingF(0.001);
    clip_mat->SetRestitution(0.8f);

    double clip_w = 0.015;
    double clip_h = 0.04;
    double clip_r = 0.001;
    double clip_g = 0.004;
    double closed_mass = 0.002695;
    double open_mass = 0.0026;
    auto cclip_iner_odd = ChVector<>(0.125176E-6, 0.556291236E-6, 0.680161092E-6);
    auto cclip_iner_even = ChVector<>(0.125176E-6, 0.680161092E-6, 0.556291236E-6);
    auto rot_odd = Q_ROTATE_X_TO_Y; //Q_ROTATE_X_TO_Y; //ChQuaternion<double>(0.7565, 0, 0.6, 0.8);
    auto rot_even = Q_ROTATE_X_TO_Y * Q_ROTATE_Y_TO_Z; //Q_ROTATE_X_TO_Y; //ChQuaternion<double>(0.7565, 0, 0.6, 0.8);

    double height = n_clips * (clip_h-3.0*clip_r);
    ChVector<> pos = ChVector<>(0.0, height, 0.0);
    std::random_device rd;
    std::default_random_engine eng(rd());
    std::uniform_real_distribution<float> distr(0.0, 1.0);
    for (int nc=0; nc<n_clips; ++nc)
    {
        // std::cout << "making clip : " << nc << " \n";
        double mass;
        ChQuaternion<double> rot;
        ChVector<> inertia;

    int velo_scale = nc>0 ? 1.0 : nc<0 ? -1.0 : 0;
    // auto init_ang_velo = (5.0*M_PI*velo_scale) * ChVector<>(1.0, 0.0, 0.0);
    // auto init_velo = (0.05*velo_scale) * ChVector<>(-1.0, -1.0, 0.0);
    auto init_ang_velo = (5.0*M_PI) * ChVector<>(distr(eng), distr(eng), distr(eng));
    auto init_velo = (0.05) * ChVector<>(distr(eng), distr(eng), distr(eng));
    // if (id % 2 != 0)
    // {
    //     init_ang_velo = -init_ang_velo;
    //     init_velo = -init_velo;
    // }
    // auto comp_pos = ChVector<>(-0.5*clip_h, 0, 0);
    // ChQuaternion<>(1, 0, 0, 0)
    // if (id % 2 != 0)
    // {
    //     // init_ang_velo = ChVector<>((20.0*M_PI), 0.0, 10.0*M_PI);
    //     // inertia = chVector<> (inertia(0), inertia(2), inertia(1));
    //     clip_color = chrono_types::make_shared<ChColorAsset>(ChColor(0.1f, 0.1f, 0.8f));
    // }



        if (nc % 2 == 0)
        {
            inertia = cclip_iner_even;
            rot = rot_even;
            init_ang_velo = -init_ang_velo;
            init_velo = -init_velo;
        }
        else
        {
            inertia = cclip_iner_odd;
            rot = rot_odd;
        }
        if (clip_type=="closed")
        {
            mass = closed_mass;
            create_closedclip(&mphysicalSystem, clip_mat, nc, 
                pos, rot, clip_w, clip_h, clip_r,
                mass, inertia);
        }
        else if (clip_type=="open")
        {
            mass = open_mass;
            // std::cout << "right before calling clip \n";
            create_clip(&mphysicalSystem, clip_mat, nc, 
                pos, rot, clip_w, clip_h, clip_r, clip_g,
                mass, inertia);
        }
        else
            std::cout << "Wrong clip type, it should be either open or closed!.....\n";
        pos -= ChVector<>(0.0, (clip_h-3.0*clip_r), 0.0);
        std::string bname = "clip" + std::to_string(nc);
        auto clip = mphysicalSystem.SearchBody(bname.c_str());
        // std::cout << "body is : " << clip->GetId() << std::end;
        if (nc == 0)
            clip->SetBodyFixed(true);
        else
        {
            clip->SetWvel_loc(ChVector<>(init_ang_velo));
            clip->SetPos_dt(init_velo);
        }
    }

    // auto first_clip = mphysicalSystem.SearchBody("clip0");
    // auto clip_motion = chrono_types::make_shared<ChLinkLockLock>();
    // clip_motion->Initialize(first_clip , floorBody, ChCoordsys<>(ChVector<>(0, 0, 0)));
    // auto mmotion_x = chrono_types::make_shared<ChFunction_Sine>(0, 15.0, 0.008);  // phase freq ampl
    // clip_motion->SetMotion_X(mmotion_x);
    // auto mmotion_y = chrono_types::make_shared<ChFunction_Sine>(0, 15.0, 0.008);  // phase freq ampl
    // clip_motion->SetMotion_Y(mmotion_y);
    // auto mmotion_z = chrono_types::make_shared<ChFunction_Sine>(0, 15.0, 0.008);  // phase freq ampl
    // clip_motion->SetMotion_Z(mmotion_z);
    // mphysicalSystem.Add(clip_motion);
}

void create_fountain(ChSystemMulticoreNSC& mphysicalSystem, int n_clips = 20, std::string clip_type = "closed") {
    ChVector<> gravity(0, -9.81, 0);
    mphysicalSystem.Set_G_acc(gravity);

    auto floor_mat = chrono_types::make_shared<ChMaterialSurfaceNSC>();
    floor_mat->SetFriction(0.01);
    // floor_mat->SetCompliance(0.000001);
    floor_mat->SetRestitution(0.2);
    auto floorBody = std::shared_ptr<ChBody>(mphysicalSystem.NewBody());
    floorBody->GetCollisionModel()->ClearModel();
    floorBody->SetCollide(true);
    floorBody->GetCollisionModel()->SetDefaultSuggestedEnvelope(0.001);
    floorBody->GetCollisionModel()->SetDefaultSuggestedMargin(0.0001);
    floorBody->GetCollisionModel()->AddBox(floor_mat, 0.09, 0.005, 0.09, ChVector<>(0.0, -0.0025, 0), QUNIT);


    // auto floorBody = chrono_types::make_shared<ChBodyEasyBox>(1., 0.005, 1.0, 7800, true, false);
    floorBody->SetGid(10001);
    // floorBody->SetPos(ChVector<>(0, -0.0025, 0));
    floorBody->SetBodyFixed(true);
    mphysicalSystem.Add(floorBody);

    auto clip_mat = chrono_types::make_shared<ChMaterialSurfaceNSC>();
    clip_mat->SetFriction(0.05f);
    // clip_mat->SetRollingFriction(0.001);
    // clip_mat->SetSpinningFriction(0.01);
    // clip_mat->SetCompliance(0.0001);
    // clip_mat->SetComplianceT(0.02);
    // clip_mat->SetComplianceRolling(0.08);
    // clip_mat->SetComplianceSpinning(0.08);
    // clip_mat->SetDampingF(0.001);
    clip_mat->SetRestitution(0.8f);

    double clip_w = 0.015;
    double clip_h = 0.04;
    double clip_r = 0.001;
    double clip_g = 0.004;
    double closed_mass = 0.002695;
    double open_mass = 0.0026;
    auto cclip_iner_odd = ChVector<>(0.125176E-6, 0.556291236E-6, 0.680161092E-6);
    auto cclip_iner_even = ChVector<>(0.125176E-6, 0.680161092E-6, 0.556291236E-6);
    auto rot_odd = Q_ROTATE_Y_TO_Z * Q_ROTATE_X_TO_Y; //Q_ROTATE_X_TO_Y; //ChQuaternion<double>(0.7565, 0, 0.6, 0.8);
    auto rot_even = Q_ROTATE_X_TO_Z; //Q_ROTATE_X_TO_Y; //ChQuaternion<double>(0.7565, 0, 0.6, 0.8);

    int n_start = (int) n_clips/2;
    double pitch = 0.025;
    double h_r = 2. * clip_h;
    double inc_theta = std::asin(0.75*clip_h / h_r);
    ChVector<> pos;
    for (int nc=0; nc<n_clips; ++nc)
    {
        // int nloc = nc - n_start;
        double mass;
        ChQuaternion<double> rot;
        ChVector<> inertia;
        ChQuaternion<double> inc_rot;
        ChQuaternion<double> inc_rot_h;
        if (nc % 2 == 0)
        {
            inertia = cclip_iner_even;
            rot = rot_even;
            inc_rot.Q_from_AngY(nc*inc_theta);
            inc_rot_h.Q_from_AngZ(-pitch);
        }
        else
        {
            inertia = cclip_iner_odd;
            rot = rot_odd;
            inc_rot.Q_from_AngZ(-nc*inc_theta);
            inc_rot_h.Q_from_AngY(-pitch);
        }
        rot = rot * inc_rot * inc_rot_h;
        pos = ChVector<>(h_r*std::cos(nc*inc_theta), 0.5*clip_h + 0.5*pitch*nc*inc_theta/M_PI, -h_r*std::sin(nc*inc_theta));
        if (clip_type=="closed")
        {
            mass = closed_mass;
            create_closedclip(&mphysicalSystem, clip_mat, nc, 
                pos, rot, clip_w, clip_h, clip_r,
                mass, inertia);
        }
        else if (clip_type=="open")
        {
            mass = open_mass;
            // std::cout << "right before calling clip \n";
            create_clip(&mphysicalSystem, clip_mat, nc, 
                pos, rot, clip_w, clip_h, clip_r, clip_g,
                mass, inertia);
        }
        else
            std::cout << "Wrong clip type, it should be either open or closed!.....\n";

    }
}



void create_model(ChSystemMulticoreNSC& mphysicalSystem, int n_clips = 10, std::string clip_type = "closed") {
    ChVector<> gravity(0, -9.81, 0);
    mphysicalSystem.Set_G_acc(gravity);
    // SetChronoDataPath("/home/res_data/");
    // SetChronoDataPath("/Users/farhad/work/cs_master/chrono/data/");
    // make the shaking box
    auto floorBody = chrono_types::make_shared<ChBodyEasyBox>(1., 0.2, 1.0, 7800, false, false);
    floorBody->SetGid(10001);
    floorBody->SetPos(ChVector<>(0, -0.3, 0));
    floorBody->SetBodyFixed(true);
    mphysicalSystem.Add(floorBody);

    double box_x = 0.125;
    double box_y = 0.125;
    double box_z = 0.15;
    double box_t = 0.01;
    int box_main_id = 10002;
    auto box_mat = chrono_types::make_shared<ChMaterialSurfaceNSC>();
    box_mat->SetFriction(0.2);
    // box_mat->SetCompliance(0.000001);
    box_mat->SetRestitution(0.6);
    auto box_pos = ChVector<>(0, 0, 0);
    create_box(&mphysicalSystem, floorBody, box_mat, box_main_id,
                box_pos, QUNIT, box_x, box_y, box_z, box_t, 7800);

    // model_box(mphysicalSystem, floorBody, box_mat, 
    //     box_main_id, box_pos, QUNIT,
    //     box_x, box_y, box_z, box_t,
    //     7800);

    auto clip_mat = chrono_types::make_shared<ChMaterialSurfaceNSC>();
    clip_mat->SetFriction(0.2f);
    // clip_mat->SetRollingFriction(0.001);
    // clip_mat->SetSpinningFriction(0.01);
    // clip_mat->SetCompliance(0.0001);
    // clip_mat->SetComplianceT(0.02);
    // clip_mat->SetComplianceRolling(0.08);
    // clip_mat->SetComplianceSpinning(0.08);
    // clip_mat->SetDampingF(0.001);
    clip_mat->SetRestitution(0.6f);

    double clip_w = 0.015;
    double clip_h = 0.04;
    double clip_r = 0.001;
    double clip_g = 0.004;
    double closed_mass = 0.002695;
    double open_mass = 0.0026;
    auto cclip_iner_odd = ChVector<>(0.125176E-6, 0.556291236E-6, 0.680161092E-6);
    auto cclip_iner_even = ChVector<>(0.125176E-6, 0.680161092E-6, 0.556291236E-6);
    auto rot_odd = QUNIT; //Q_ROTATE_X_TO_Y; //ChQuaternion<double>(0.7565, 0, 0.6, 0.8);
    auto rot_even = Q_ROTATE_Y_TO_Z; //Q_ROTATE_X_TO_Y; //ChQuaternion<double>(0.7565, 0, 0.6, 0.8);

    // int n_clips = 10;
    int n_start = (int) n_clips/2;
    ChVector<> pos = ChVector<>(-0.5*box_x+clip_h, 0, -0.5*box_z+clip_h);
    for (int nc=0; nc<n_clips; ++nc)
    {
        // int nloc = nc - n_start;
        double mass;
        ChQuaternion<double> rot;
        ChVector<> inertia;
        if (nc % 2 == 0)
        {
            inertia = cclip_iner_even;
            rot = rot_even;
        }
        else
        {
            inertia = cclip_iner_odd;
            rot = rot_odd;
        }
        if (clip_type=="closed")
        {
            mass = closed_mass;
            create_closedclip(&mphysicalSystem, clip_mat, nc, 
                pos, rot, clip_w, clip_h, clip_r,
                mass, inertia);
        }
        else if (clip_type=="open")
        {
            mass = open_mass;
            // std::cout << "right before calling clip \n";
            create_clip(&mphysicalSystem, clip_mat, nc, 
                pos, rot, clip_w, clip_h, clip_r, clip_g,
                mass, inertia);
        }
        else
            std::cout << "Wrong clip type, it should be either open or closed!.....\n";
        // if (nc % 2 == 0)
        //     create_closedclip(&mphysicalSystem, clip_mat, nc, 
        //         pos, rot_even,
        //         clip_w, clip_h, clip_r,
        //         cclip_mass, cclip_iner_even);
        // else
        //     create_closedclip(&mphysicalSystem, clip_mat, nc, 
        //         pos, rot_odd,
        //         clip_w, clip_h, clip_r,
        //         cclip_mass, cclip_iner_odd);
        if (pos.x()+clip_h>0.5*box_x)
        {
            if (pos.z()+clip_h>0.5*box_z)
            {
                pos = ChVector<>(-0.5*box_x+clip_h, pos.y()+0.55*clip_h, -0.5*box_z+clip_h);
            }
            else
            {
                pos = ChVector<>(-0.5*box_x+clip_h, pos.y(), pos.z()+0.55*clip_h);
            }
        }
        else
        {
            pos += ChVector<>(0.55*clip_h, 0, 0);
        }


        // if (pos.z()>0.5*box_z-1.1*clip_h)
        // {
        //     // std::cout << "IF Befroe : pos is : " << pos << std::endl;
        //     pos = ChVector<>(-0.5*box_x+clip_h, pos.y()+0.55*clip_h, -0.5*box_z+clip_h);
        //     continue;
        //     // std::cout << "pos is : " << pos << std::endl;
        // }
        // if (pos.x()>0.5*box_x-1.1*clip_h)
        // {
        //     // std::cout << "ELSEIF Befroe : pos is : " << pos << std::endl;
        //     pos = ChVector<>(-0.5*box_x+clip_h, pos.y(), pos.z()+0.55*clip_h);
        //     continue;
        //     // std::cout << "pos is : " << pos << std::endl;
        // }
        // if ((pos.x()<0.5*box_x-1.1*clip_h) && (pos.z()<0.5*box_z-1.1*clip_h))
        // {
        //     // std::cout << "ELSE Befroe : pos is : " << pos << std::endl;
        //     pos += ChVector<>(0.55*clip_h, 0, 0);
        //     continue;
        //     // std::cout << "pos is : " << pos << std::endl;
        // }
        // // std::cout << "pos is : " << pos << std::endl;
    }
    // std::cout << "after for loop ... \n";
    // Adding torus///////////////////////////
    // double torus_r = 0.01;
    // double torus_t = 0.001;
    // int torus_nseg = 20;
    // double torus_angle = 360.;
    // double torus_mass = 0.003;
    // auto torus_iner = ChVector<>(0.12E-6, 0.7E-6, 0.7E-6);

    // auto pos = ChVector<>(0, 0, 0);
    // auto rot_yz = Q_ROTATE_Y_TO_Z;

    // create_torus(mphysicalSystem, clip_mat, 1, 
    //         pos, rot_odd,
    //         torus_r, torus_t, torus_nseg, torus_angle,
    //         torus_mass, torus_iner);

    // pos = ChVector<>(1.5*torus_r, 0, 0);
    // create_torus(mphysicalSystem, clip_mat, 1, 
    //         pos, rot_yz,
    //         torus_r, torus_t, torus_nseg, torus_angle,
    //         torus_mass, torus_iner);
    
    // create_closedclip(mphysicalSystem, clip_mat, 1, 
    //         pos, rot,
    //         clip_w, clip_h, clip_r,
    //         cclip_mass, cclip_iner);
    // pos = ChVector<>(0.5*clip_h, 0, 0);
    // create_closedclip(mphysicalSystem, clip_mat, 2, 
    //         pos, rot_yz,
    //         clip_w, clip_h, clip_r,
    //         cclip_mass, cclip_iner_yz);
}

void create_box_clip(ChSystemMulticoreNSC& mphysicalSystem, int n_clips = 10, std::string clip_type = "closed") {
    ChVector<> gravity(0, -9.81, 0);
    mphysicalSystem.Set_G_acc(gravity);
    // SetChronoDataPath("/home/res_data/");
    // SetChronoDataPath("/Users/farhad/work/cs_master/chrono/data/");
    // make the shaking box
    auto floorBody = chrono_types::make_shared<ChBodyEasyBox>(1., 0.2, 1.0, 7800, false, false);
    floorBody->SetGid(10001);
    floorBody->SetName("floor");
    floorBody->SetPos(ChVector<>(0, -0.3, 0));
    floorBody->SetBodyFixed(true);
    mphysicalSystem.Add(floorBody);

    double box_x = 0.125;
    double box_y = 0.125;
    double box_z = 0.125;
    double box_t = 0.01;
    int box_main_id = 10002;
    auto box_mat = chrono_types::make_shared<ChMaterialSurfaceNSC>();
    box_mat->SetFriction(0.2);
    // box_mat->SetCompliance(0.000001);
    box_mat->SetRestitution(0.6);
    auto box_pos = ChVector<>(0, 0, 0);
    create_box(&mphysicalSystem, floorBody, box_mat, box_main_id,
                box_pos, QUNIT, box_x, box_y, box_z, box_t, 7800);

    auto clip_mat = chrono_types::make_shared<ChMaterialSurfaceNSC>();
    clip_mat->SetFriction(0.2f);
    // clip_mat->SetRollingFriction(0.001);
    // clip_mat->SetSpinningFriction(0.01);
    // clip_mat->SetCompliance(0.0001);
    // clip_mat->SetComplianceT(0.02);
    // clip_mat->SetComplianceRolling(0.08);
    // clip_mat->SetComplianceSpinning(0.08);
    // clip_mat->SetDampingF(0.001);
    clip_mat->SetRestitution(0.6f);

    double clip_w = 0.015;
    double clip_h = 0.04;
    double clip_r = 0.001;
    double clip_g = 0.004;
    double closed_mass = 0.002695;
    double open_mass = 0.0026;
    auto cclip_iner_odd = ChVector<>(0.125176E-6, 0.556291236E-6, 0.680161092E-6);
    auto cclip_iner_even = ChVector<>(0.125176E-6, 0.680161092E-6, 0.556291236E-6);
    auto rot_odd = QUNIT; //Q_ROTATE_X_TO_Y; //ChQuaternion<double>(0.7565, 0, 0.6, 0.8);
    auto rot_even = Q_ROTATE_Y_TO_Z; //Q_ROTATE_X_TO_Y; //ChQuaternion<double>(0.7565, 0, 0.6, 0.8);

    // int n_clips = 10;
    int n_start = (int) n_clips/2;
    double clip_d = 0.55*clip_h;
    ChVector<> pos = ChVector<>(-0.5*box_x+clip_d, -0.5*box_y+clip_d, -0.5*box_z+clip_d);
    for (int nc=0; nc<n_clips; ++nc)
    {
        // int nloc = nc - n_start;
        double mass;
        ChQuaternion<double> rot;
        ChVector<> inertia;
        if (nc % 2 == 0)
        {
            inertia = cclip_iner_even;
            rot = rot_even;
        }
        else
        {
            inertia = cclip_iner_odd;
            rot = rot_odd;
        }
        if (clip_type=="closed")
        {
            mass = closed_mass;
            create_closedclip(&mphysicalSystem, clip_mat, nc, 
                pos, rot, clip_w, clip_h, clip_r,
                mass, inertia);
        }
        else if (clip_type=="open")
        {
            mass = open_mass;
            // std::cout << "right before calling clip \n";
            create_clip(&mphysicalSystem, clip_mat, nc, 
                pos, rot, clip_w, clip_h, clip_r, clip_g,
                mass, inertia);
        }
        else
            std::cout << "Wrong clip type, it should be either open or closed!.....\n";

        if (pos.x()+2.0*clip_d>0.5*box_x)
        {
            if (pos.z()+2.0*clip_d>0.5*box_z)
            {
                pos = ChVector<>(-0.5*box_x+clip_d, pos.y()+clip_d, -0.5*box_z+clip_d);
            }
            else
            {
                pos = ChVector<>(-0.5*box_x+clip_d, pos.y(), pos.z()+clip_d);
            }
        }
        else
        {
            pos += ChVector<>(clip_d, 0, 0);
        }
    }
}

void writeVTU(ChSystemMulticoreNSC &sys, int nframe, std::string clip_type="closed")
{
    double clip_w = 0.015;
    double clip_h = 0.04;
    double clip_r = 0.001;
    double clip_g = 0.004;

    std::vector< std::vector<int> > clip_edges;
    std::vector<std::vector<double> > clip_vertices;
    std::vector<std::vector<double> > all_vertices;
    std::vector< std::vector<int> > all_edges;

    int nv = 6;
    int ne = 5;
    if (clip_type=="open")
    {
        clip_vertices.resize(nv);
        clip_vertices[0] = std::vector<double>({0.5*clip_g, -0.5*clip_w, 0.0});
        clip_vertices[1] = std::vector<double>({0.5*clip_h, -0.5*clip_w, 0.0});
        clip_vertices[2] = std::vector<double>({0.5*clip_h, 0.5*clip_w, 0.0});
        clip_vertices[3] = std::vector<double>({-0.5*clip_h, 0.5*clip_w, 0.0});
        clip_vertices[4] = std::vector<double>({-0.5*clip_h, -0.5*clip_w, 0.0});
        clip_vertices[5] = std::vector<double>({-0.5*clip_g, -0.5*clip_w, 0.0});
        clip_edges.resize(ne);
        clip_edges[0] = std::vector<int>({0, 1});
        clip_edges[1] = std::vector<int>({1, 2});
        clip_edges[2] = std::vector<int>({2, 3});
        clip_edges[3] = std::vector<int>({3, 4});
        clip_edges[4] = std::vector<int>({4, 5});
    }
    else if (clip_type=="closed")
    {
        nv = 4;
        ne = 4;
        clip_vertices.resize(nv);
        clip_vertices[0] = std::vector<double>({0.5*clip_h, -0.5*clip_w, 0.0});
        clip_vertices[1] = std::vector<double>({0.5*clip_h, 0.5*clip_w, 0.0});
        clip_vertices[2] = std::vector<double>({-0.5*clip_h, 0.5*clip_w, 0.0});
        clip_vertices[3] = std::vector<double>({-0.5*clip_h, -0.5*clip_w, 0.0});
        clip_edges.resize(ne);
        clip_edges[0] = std::vector<int>({0, 1});
        clip_edges[1] = std::vector<int>({1, 2});
        clip_edges[2] = std::vector<int>({2, 3});
        clip_edges[3] = std::vector<int>({3, 0});
    }
    else
    {
        std::cout << "Wrong value for clip_type, it should be either open or closed! ....\n";
    }

    ppv::fmt::e_VtkType elemType = ppv::fmt::LINE;
    std::ofstream vtu("/home/multicore_clip/data/results/cclip_" + std::to_string(nframe) + ".vtu");
    // std::cout << "num of clips : " << (sys.Get_bodylist().size()-2) << std::endl;
    int nclips = 0;//sys.Get_bodylist().size()-2;
    
    // std::cout << "Body list size  \n" << sys.Get_bodylist().size() << std::endl;
    ChMatrix33<> rotationZX;
    rotationZX << 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0;
    int body_count = 0;
    for (auto body : sys.Get_bodylist())
    {
        std::string bname = body->GetName(); 
        if (bname.substr(0, 4) == "clip")
        {
            nclips++;
            // std::cout << "body ID : " << body->GetPos() << std::endl;
            // if (body_count == 0)
            // {
            //     body_count++;
            //     continue;
            // }
            ChVector<> translation = body->GetPos();
            ChMatrix33<> rotation = body->GetRot();
            // if (body_count % 2 != 0)
            //         rotation = rotationZX * rotation; 
            ChVector<double> clip_vrtx;
            ChVector<double> global_vrtx;
            for (auto vrtx : clip_vertices)
            {
                clip_vrtx.x() = vrtx[0];
                clip_vrtx.y() = vrtx[1];
                clip_vrtx.z() = vrtx[2];
                global_vrtx = rotation * clip_vrtx + translation; 
                // std::cout << "test vrtx : " << vrtx[0] << ", " << vrtx[1] << ", " << vrtx[2] << std::endl;
                all_vertices.push_back({global_vrtx.x(), global_vrtx.y(), global_vrtx.z()});
            }

            for (auto edg : clip_edges)
            {
                // std::cout << "test edge : " << edg[0]+6*body_count << ", " << edg[1]+6*body_count << std::endl;
                all_edges.push_back({edg[0]+nv*body_count, edg[1]+nv*body_count});
            }
            body_count++;
        }
        else if (bname.substr(0, 3) == "box")
        {
            double box_x = 0.125;
            double box_y = 0.125;
            double box_z = 0.125;
            double box_t = 0.01;
            ChVector<> translation = body->GetPos();
            ChMatrix33<> rotation = body->GetRot();
            std::vector<std::vector<double> > box_vertices;
            std::vector< std::vector<int> > box_faces;
            box_vertices.resize(8);
            box_vertices[0] = std::vector<double>({-0.5*box_x, -0.5*box_y, -0.5*box_z});
            box_vertices[1] = std::vector<double>({0.5*box_x, -0.5*box_y, -0.5*box_z});
            box_vertices[2] = std::vector<double>({0.5*box_x, 0.5*box_y, -0.5*box_z});
            box_vertices[3] = std::vector<double>({-0.5*box_x, 0.5*box_y, -0.5*box_z});
            box_vertices[4] = std::vector<double>({-0.5*box_x, -0.5*box_y, 0.5*box_z});
            box_vertices[5] = std::vector<double>({0.5*box_x, -0.5*box_y, 0.5*box_z});
            box_vertices[6] = std::vector<double>({0.5*box_x, 0.5*box_y, 0.5*box_z});
            box_vertices[7] = std::vector<double>({-0.5*box_x, 0.5*box_y, 0.5*box_z});
            std::vector<std::vector<double> > box_vrtx_all;
            ChVector<double> box_vrtx;
            ChVector<double> global_vrtx_box;
            for (auto vrtx : box_vertices)
            {
                box_vrtx.x() = vrtx[0];
                box_vrtx.y() = vrtx[1];
                box_vrtx.z() = vrtx[2];
                global_vrtx_box = rotation * box_vrtx + translation; 
                box_vrtx_all.push_back({global_vrtx_box.x(), global_vrtx_box.y(), global_vrtx_box.z()});
            }
            box_faces.resize(6);
            box_faces[0] = std::vector<int>({0, 1, 2, 3});
            box_faces[1] = std::vector<int>({4, 7, 6, 5});
            box_faces[2] = std::vector<int>({0, 1, 5, 4});
            box_faces[3] = std::vector<int>({3, 2, 6, 7});
            box_faces[4] = std::vector<int>({0, 3, 7, 4});
            box_faces[5] = std::vector<int>({1, 5, 6, 2});
            ppv::fmt::e_VtkType elemType = ppv::fmt::QUAD;
            std::ofstream vtubox("/home/multicore_clip/data/results/box_" + std::to_string(nframe) + ".vtu");
            ppv::fmt::OpenVTU(vtubox, 8, 6);
            ppv::WriteGeometryVTU(vtubox, box_vrtx_all, box_faces, ppv::fmt::QUAD);
            ppv::fmt::CloseVTU(vtubox);
            vtubox.close();
        }
    }
    // std::cout << "all vertices size : \n" << all_vertices.size() << std::endl;
    ppv::fmt::OpenVTU(vtu, nv*nclips, ne*nclips);
    ppv::WriteGeometryVTU(vtu, all_vertices, all_edges, elemType);
    ppv::fmt::CloseVTU(vtu);
    vtu.close();
}

float ComputeTotalKE(ChSystemMulticoreNSC &sys)
{
    // A very simple simulation loop..
    float ke = 0.0;
    for (auto body : sys.Get_bodylist())
    {
        std::string bname = body->GetName(); 
        if (bname.substr(0, 4) == "clip")
        {
            ChVector<> velo = body->GetPos_dt();
            ChQuaternion<> qrot_dt = body->GetRot_dt();
            ChQuaternion<> rot = body->GetRot();
            // ChMatrix33<> rot_mat = ChMatrix33<>(rot);
            // ChQuaternion<> ang_vel_Q = rot.GetInverse() * (qrot_dt*(2.0));
            ChVector<> angVrel;
            qrot_dt.Qdt_to_Wrel(angVrel, rot); //ang_vel_Q.Q_to_Rotv(); //ChVector<>(ang_vel_Q.e1(), ang_vel_Q.e2(), ang_vel_Q.e3());
            // ang_vel *= std::sqrt(2.0);
            float mass = body->GetMass();
            ChMatrix33<> inertia_loc = body->GetInertia();
            // ChMatrix33<> inertia_global = (rot_mat)*(body->GetInertia()*(rot_mat.transpose()));
            // float inertia = 0.006 - (ang_vel.GetNormalized()).Dot( inertia_global * (ang_vel.GetNormalized()));
            ke += 0.5 * ( mass * velo.Length2() + angVrel.Dot(inertia_loc * angVrel));
            // std::cout << "\n body " << body->GetId() << ", angular velocity : \n" << ang_vel << std::endl;
            // std::cout << "\n body " << body->GetId() << ", inertia is : \n" << rot_mat << std::endl;
            // std::cout << "\n body rotation is : \n" << qrot_dt;
        }
    }
    return ke;
}


int main(int argc, char* argv[]) 
{
    GetLog() << "Copyright (c) 2017 projectchrono.org\nChrono version: " << CHRONO_VERSION << "\n\n";
    int num_of_threads = 1;
    int num_of_clips = 6;
    // std::string type_of_clips = "closed";

    ChSystemMulticoreNSC mphysicalSystem;
    mphysicalSystem.SetNumThreads(num_of_threads);
    ChVector<> gravity(0, -9.81, 0);
    mphysicalSystem.Set_G_acc(gravity);

    ClipSystem csys(mphysicalSystem);
    double clip_h = 0.04;
    double clip_w = 0.02;
    double clip_r = 0.001;
    double clip_g = 0.003;
    double clip_density = 4000.0;
    double clip_fric = 0.1;
    double clip_rest = 0.5;
    std::string clip_name = "open_clip";
    ChVector<> clip_center(0.0, 0.0, 0.0);
    
    csys.make_clip(clip_h, clip_w, clip_r, clip_g, 0.0, clip_density, clip_fric, clip_rest, clip_name);
    ChVector<> dim(0.25, 0.15, 0.1);
    ChVector<> center(0.0, 0.0, 0.0);
    double b_den = 7800.0;
    double b_fric = 0.2;
    double b_rest = 0.4;
    std::string box_name = "box";
    csys.make_box(dim.x(), dim.y(), dim.z(), 0.02, b_den, b_fric, b_rest, "box");

    int num_clip = 5;
    double box_wall_t = 0.02;
    csys.add_floor(ChVector<>(0.1, 0.001, 0.1), ChVector<>(0.0, 0.0, 0.0), 1);
    csys.add_box(box_name, 2, center, QUNIT);
    csys.add_motion("box_2", "floor_1", "vibration", ChVector<>(0, 0, 0), ChVector<>(2, 3, 4), ChVector<>(0.005, 0.008, 0.01));
    csys.add_clip(clip_name, 3, clip_center, QUNIT);
    
    csys.clip_in_box(dim, center, num_clip, 0.01, clip_name);
    std::cout << "num of bodies : " << csys.n_body << std::endl;

    int max_iteration = 200;
    double tolerance = 1e-7;
    csys.sys->GetSettings()->solver.solver_type = SolverType::APGD;
    // csys.sys->GetSettings()->solver.solver_mode = SolverMode::SLIDING;
    // mphysicalSystem.GetSettings()->solver.max_iteration_normal = max_iteration / 3;
    // mphysicalSystem.GetSettings()->solver.max_iteration_sliding = max_iteration / 3;
    // mphysicalSystem.GetSettings()->solver.max_iteration_spinning = 0;
    // mphysicalSystem.GetSettings()->solver.max_iteration_bilateral = max_iteration / 3;
    // csys.sys->GetSettings()->solver.max_iteration = max_iteration;
    // csys.sys->GetSettings()->solver.tolerance = tolerance;
    // mphysicalSystem.GetSettings()->solver.alpha = 0;
    // mphysicalSystem.GetSettings()->solver.contact_recovery_speed = 10.0;
    // mphysicalSystem.ChangeSolverType(SolverType::APGDREF);


    double dT = 0.005;
    double endT = 2.0;
    int vtu_interval = 10; //1./dT / 20.0 < 1 ? 1 : (int) 1./dT / 100.0;
    int txt_interval = 1./dT / 100.0 < 1 ? 1 : (int) 1./dT / 100.0;
    // GetLog() << "  List of bodies: " << mphysicalSystem.Get_bodylist()[4]->GetCoord().pos.x() << "\n";
    csys.sys->SetMaxPenetrationRecoverySpeed(10.0);
    csys.sys->SetMinBounceSpeed(0.0001);

    double chronoTime = 0.0;
    auto str_chrono = std::chrono::system_clock::now();
    // writeVTU(mphysicalSystem, 0, type_of_clips);
    int nframe = 0;
    int nvtu = 0;
    // std::ofstream kin_eng("/home/project/clips/sim_chrono/outputs/kin_eng.txt");
    double start = std::clock();
    // double ke = ComputeTotalKE(mphysicalSystem);
    // kin_eng << chronoTime << ", " << ke << std::endl;
    ParaVisual pvis("/home/project/clips/sim_chrono/outputs", "sys1", dT);
    while (chronoTime < endT) 
    {
        // if (nframe % txt_interval == 0)
        // {
        //     ke = ComputeTotalKE(mphysicalSystem);
        //     kin_eng << chronoTime << ", " << ke << std::endl;
        // }
        // std::vector<std::vector<double>> new_verts = csys.get_updated_location("open_clip", 3);
        // std::cout << "frame is : " << new_verts.at(2).at(1) << std::endl;
        if (nframe % vtu_interval == 0)
        {
            nvtu++;
            for (auto body : csys.sys->Get_bodylist())
            {
                std::string bname = body->GetName(); 
                if (bname.substr(0, clip_name.size()) == clip_name)
                {
                    int gid = body->GetGid();
                    std::vector<std::vector<double>> verts = csys.get_updated_location(clip_name, gid);
                    pvis.write_wire_obj(verts, csys.bodies[clip_name]->elems, "clip", gid, nframe);
                }
                if (bname.substr(0, box_name.size()) == box_name)
                {
                    int gid = body->GetGid();
                    std::vector<std::vector<double>> verts = csys.get_updated_location(box_name, gid);
                    // std::cout << "sample vert : " << verts[0].at(0) << " , " << verts[0].at(1) << " , " << verts[0].at(2) << std::endl;
                    pvis.write_face_obj(verts, csys.bodies[box_name]->elems, "box", gid, nframe);
                }
                
                    
            }
            // pvis.write_wire_obj(verts, csys.clips["open_clip"].edges, "clip", id, ts);
        }
        chronoTime += dT;
        csys.sys->DoStepDynamics(dT);
        nframe++;
            
    }
    // std::string pvdn = "PVD" + clip_name;
    // pvis.write_pvd(pvdn);

    // ke = ComputeTotalKE(mphysicalSystem);
    // kin_eng << chronoTime << ", " << ke << std::endl;
    // kin_eng.close();

    auto end_chorno = std::chrono::system_clock::now();
    double duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end_chorno - str_chrono);
    std::cout << "simulation run time is : " << elapsed.count() << " mili-seconds" << '\n';
    std::cout << "simulation run time is : " << (double) duration/num_of_threads << " seconds" << std::endl;
    return 0;
}