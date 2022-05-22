// #include "chrono/physics/ChSystemNSC.h"
// #include "chrono/physics/ChBodyEasy.h"
// #include "chrono/solver/ChSolverPSOR.h"
// #include "chrono/assets/ChTexture.h"
// #include "chrono/physics/ChLinkMate.h"
// #include "chrono_irrlicht/ChIrrApp.h"
// #include "utils/ChUtilsCreators.h"
// #include "utils/ChUtilsInputOutput.h"
// #include "utils/ChUtilsGenerators.h"
// #include "chrono_multicore/physics/ChSystemMulticore.h"

// #include <random>

// using namespace chrono;
// using namespace chrono::irrlicht;

// // Use the main namespaces of Irrlicht
// using namespace irr;
// // using namespace irr::core;
// using namespace irr::scene;
// using namespace irr::video;
// using namespace irr::io;
// using namespace irr::gui;

// void model_box(ChSystemNSC& sys, std::shared_ptr<ChBody> floor, std::shared_ptr<ChMaterialSurface> mat,
//                 const int id, ChVector<>& pos, const ChQuaternion<double>& rot,
//                 const double box_x, const double box_y, const double box_z, const double box_t,
//                 const double density)
// {
//     auto texture = chrono_types::make_shared<ChTexture>();
//     texture->SetTextureFilename(GetChronoDataFile("textures/concrete.jpg"));
//     auto box = std::make_shared<ChBody>();
// 	box->SetIdentifier(id);
    
//     box->GetCollisionModel()->ClearModel();
//     box->SetCollide(true);
//     box->GetCollisionModel()->SetDefaultSuggestedEnvelope(0.01);
//     box->GetCollisionModel()->SetDefaultSuggestedMargin(0.0001);
//     double alpha = 1.1;
//     box->GetCollisionModel()->AddBox(mat, box_x+2.*box_t, box_t, box_z+2.*box_t, ChVector<>(0.0, 0.5*(box_t+box_y), 0), ChQuaternion<>(1, 0, 0, 0));
//     box->GetCollisionModel()->AddBox(mat, box_x+2.*box_t, box_t, box_z+2.*box_t, ChVector<>(0.0, -0.5*(box_t+box_y), 0), ChQuaternion<>(1, 0, 0, 0));
//     box->GetCollisionModel()->AddBox(mat, box_x+2.*box_t, box_y+2.*box_t, box_t, ChVector<>(0.0, 0.0, 0.5*(box_t+box_z)), ChQuaternion<>(1, 0, 0, 0));
//     box->GetCollisionModel()->AddBox(mat, box_x+2.*box_t, box_y+2.*box_t, box_t, ChVector<>(0.0, 0.0, -0.5*(box_t+box_z)), ChQuaternion<>(1, 0, 0, 0));
//     box->GetCollisionModel()->AddBox(mat, box_t, box_y+2.*box_t, box_z+2.*box_t, ChVector<>(0.5*(box_t+box_x), 0.0, 0.0), ChQuaternion<>(1, 0, 0, 0));
//     box->GetCollisionModel()->AddBox(mat, box_t, box_y+2.*box_t, box_z+2.*box_t, ChVector<>(-0.5*(box_t+box_x), 0.0, 0.0), ChQuaternion<>(1, 0, 0, 0));
    
//     double volume_bottop = (box_x+2.*box_t) * (box_z+2.*box_t) * box_t;
//     double volume_lefright = (box_y+2.*box_t) * (box_z+2.*box_t) * box_t;
//     double volume_fronback = (box_x+2.*box_t) * (box_y+2.*box_t) * box_t;
//     double volume = 2. * (volume_fronback + volume_lefright + volume_bottop);
//     double mass_bottop = density * volume_bottop;
//     double mass_lefright = density * volume_lefright;
//     double mass_fronback = density * volume_fronback;
//     double mass = 2. * (mass_fronback + mass_lefright + mass_bottop);
//     double iner_xyz = mass/12.0 * ((box_y+2.*box_t)*(box_y+2.*box_t) + (box_z+2.*box_t)*(box_z+2.*box_t));
//     // double iner_y = mass/12.0 * (box_x*box_x + box_z*box_z);
//     // double iner_z = mass/12.0 * (box_x*box_x + box_y*box_y);
//     auto box_iner = ChMatrix33<>(ChVector<>(iner_xyz, iner_xyz, iner_xyz));
//     box->SetMass(mass);
//     box->SetInertia(box_iner);
//     // box->SetMaterialSurface(mat);
//     // std::cout << "box materila is : " << box->
//     box->GetCollisionModel()->BuildModel();
//     sys.Add(box);
//     // clip->SetPos(pos);
//     // clip->SetRot(rot);
//     box->AddAsset(texture);

//     auto bot_viz = chrono_types::make_shared<ChBoxShape>();
//     bot_viz->GetBoxGeometry().SetLengths(ChVector<>(alpha*box_x, box_t, alpha*box_z));
//     bot_viz->Pos = ChVector<>(0.0, -0.5*(box_t+box_y), 0);
//     box->GetAssets().push_back(bot_viz);

//     auto top_viz = chrono_types::make_shared<ChBoxShape>();
//     top_viz->GetBoxGeometry().SetLengths(ChVector<>(alpha*box_x, box_t, alpha*box_z));
//     top_viz->Pos = ChVector<>(0.0, 0.5*(box_t+box_y), 0);
//     box->GetAssets().push_back(top_viz);

//     auto left_viz = chrono_types::make_shared<ChBoxShape>();
//     left_viz->GetBoxGeometry().SetLengths(ChVector<>(box_t, alpha*box_y, alpha*box_z));
//     left_viz->Pos = ChVector<>(-0.5*(box_t+box_x), 0.0, 0.0);
//     box->GetAssets().push_back(left_viz);

//     auto right_viz = chrono_types::make_shared<ChBoxShape>();
//     right_viz->GetBoxGeometry().SetLengths(ChVector<>(box_t, alpha*box_y, alpha*box_z));
//     right_viz->Pos = ChVector<>(0.5*(box_t+box_x), 0.0, 0.0);
//     box->GetAssets().push_back(right_viz);

//     // auto front_viz = chrono_types::make_shared<ChBoxShape>();
//     // front_viz->GetBoxGeometry().SetLengths(ChVector<>(alpha*box_x, alpha*box_y, box_t));
//     // front_viz->Pos = ChVector<>(0.0, 0.0, 0.5*(box_t+box_z));
//     // box->GetAssets().push_back(front_viz);

//     auto back_viz = chrono_types::make_shared<ChBoxShape>();
//     back_viz->GetBoxGeometry().SetLengths(ChVector<>(alpha*box_x, alpha*box_y, box_t));
//     back_viz->Pos = ChVector<>(0.0, 0.0, -0.5*(box_t+box_z));
//     box->GetAssets().push_back(back_viz);
    
//     auto box_motion = chrono_types::make_shared<ChLinkLockLock>();
//     box_motion->Initialize(box , floor, ChCoordsys<>(ChVector<>(0, 0, 0)));

//     auto mmotion_x = chrono_types::make_shared<ChFunction_Sine>(0, 2.0, 0.01);  // phase freq ampl
//     box_motion->SetMotion_X(mmotion_x);
//     auto mmotion_y = chrono_types::make_shared<ChFunction_Sine>(0, 2.0, 0.01);  // phase freq ampl
//     box_motion->SetMotion_Y(mmotion_y);
//     sys.Add(box_motion);

// }

// // void create_clip(ChSystemNSC& sys, std::shared_ptr<ChMaterialSurface> mat,
// //                 const int id, ChVector<>& pos, const ChQuaternion<double>& rot,
// //                 const double clip_w, const double clip_h, const double clip_r, const double clip_g,
// //                 const double mass, const ChMatrix33<> inertia)
// // {
// //     auto texture = chrono_types::make_shared<ChTexture>();
// //     texture->SetTextureFilename(GetChronoDataFile("textures/cubetexture_borders.png"));
// //     auto clip = std::make_shared<ChBody>();
// // 	clip->SetIdentifier(id);
// //     clip->GetCollisionModel()->ClearModel();
// //     clip->SetCollide(true);
// //     clip->GetCollisionModel()->SetDefaultSuggestedEnvelope(0.025);
// //     clip->GetCollisionModel()->SetDefaultSuggestedMargin(0.0005);
// //     utils::AddCylinderGeometry(clip.get(),
// //                                 mat,
// //                                 clip_r,
// //                                 0.5*clip_w,
// //                                 ChVector<>(-0.5*clip_h, 0, 0),
// //                                 ChQuaternion<>(1, 0, 0, 0),
// //                                 true);
// //     utils::AddCylinderGeometry(clip.get(),
// //                                 mat,
// //                                 clip_r,
// //                                 0.5*clip_w,
// //                                 ChVector<>(0.5*clip_h, 0, 0),
// //                                 ChQuaternion<>(1, 0, 0, 0),
// //                                 true);
// //     utils::AddCylinderGeometry(clip.get(),
// //                                 mat,
// //                                 clip_r,
// //                                 0.5*clip_h,
// //                                 ChVector<>(0.0, 0.5*clip_w, 0),
// //                                 Q_ROTATE_Y_TO_X,
// //                                 true);
// //     utils::AddCylinderGeometry(clip.get(),
// //                                 mat,
// //                                 clip_r,
// //                                 0.25*(clip_h-clip_g),
// //                                 ChVector<>(-0.25*(clip_h+clip_g), -0.5*clip_w, 0),
// //                                 Q_ROTATE_Y_TO_X,
// //                                 true);
// //     utils::AddCylinderGeometry(clip.get(),
// //                                 mat,
// //                                 clip_r,
// //                                 0.25*(clip_h-clip_g),
// //                                 ChVector<>(0.25*(clip_h+clip_g), -0.5*clip_w, 0),
// //                                 Q_ROTATE_Y_TO_X,
// //                                 true);
    
// //     clip->SetMass(mass);
// //     clip->SetInertia(inertia);
// //     clip->GetCollisionModel()->BuildModel();
// //     sys.Add(clip);
// //     clip->SetPos(pos);
// //     // clip->SetRot(rot);
// //     clip->AddAsset(texture);
// // }


// void create_closedclip(ChSystemNSC& sys, std::shared_ptr<ChMaterialSurface> mat,
//                 const int id, ChVector<>& pos, const ChQuaternion<double>& rot,
//                 const double clip_w, const double clip_h, const double clip_r,
//                 const double mass, const ChVector<> inertia)
// {
//     // auto texture = chrono_types::make_shared<ChTexture>();
//     // if (id % 2 == 0)
//     //     texture->SetTextureFilename(GetChronoDataFile("textures/pink.png"));
//     // else
//     //     texture->SetTextureFilename(GetChronoDataFile("textures/blue.png"));

//     std::random_device rd;
//     std::default_random_engine eng(rd());
//     std::uniform_real_distribution<float> distr(0.0, 1.0);

//     auto clip_color = chrono_types::make_shared<ChColorAsset>(ChColor(0.8f, 0.1f, 0.1f));
//     auto init_ang_velo = (5.0*M_PI) * ChVector<>(distr(eng), distr(eng), distr(eng));
//     auto init_velo = (0.05) * ChVector<>(distr(eng), distr(eng), distr(eng));
//     // auto comp_pos = ChVector<>(-0.5*clip_h, 0, 0);
//     // ChQuaternion<>(1, 0, 0, 0)
//     if (id % 2 != 0)
//     {
//         // init_ang_velo = ChVector<>((20.0*M_PI), 0.0, 10.0*M_PI);
//         // inertia = chVector<> (inertia(0), inertia(2), inertia(1));
//         clip_color = chrono_types::make_shared<ChColorAsset>(ChColor(0.1f, 0.1f, 0.8f));
//     }
    

//     auto clip = std::make_shared<ChBody>();
// 	clip->SetIdentifier(id);
//     // clip->GetCollisionModel()->SetDefaultSuggestedEnvelope(0.001);
//     clip->GetCollisionModel()->SetEnvelope(0.03);
//     clip->GetCollisionModel()->SetDefaultSuggestedMargin(0.0001);
//     clip->SetCollide(true);
//     clip->GetCollisionModel()->ClearModel();
//     // clip->SetMaterialSurface(mat);

//     double added_l = 1.0*clip_r;

//     utils::AddCylinderGeometry(clip.get(),
//                                 mat,
//                                 clip_r,
//                                 0.5*clip_w+added_l,
//                                 ChVector<>(-0.5*clip_h, 0, 0),
//                                 ChQuaternion<>(1, 0, 0, 0),
//                                 true);
//     utils::AddCylinderGeometry(clip.get(),
//                                 mat,
//                                 clip_r,
//                                 0.5*clip_w+added_l,
//                                 ChVector<>(0.5*clip_h, 0, 0),
//                                 ChQuaternion<>(1, 0, 0, 0),
//                                 true);
//     utils::AddCylinderGeometry(clip.get(),
//                                 mat,
//                                 clip_r,
//                                 0.5*clip_h+added_l,
//                                 ChVector<>(0.0, 0.5*clip_w, 0),
//                                 Q_ROTATE_Y_TO_X,
//                                 true);
//     utils::AddCylinderGeometry(clip.get(),
//                                 mat,
//                                 clip_r,
//                                 0.5*clip_h+added_l,
//                                 ChVector<>(0.0, -0.5*clip_w, 0),
//                                 Q_ROTATE_Y_TO_X,
//                                 true);
    
//     clip->GetCollisionModel()->BuildModel();

//     clip->SetMass(mass);
//     clip->SetInertiaXX(inertia);
//     clip->SetWvel_par(init_ang_velo);
//     clip->SetPos_dt(init_velo);
//     sys.Add(clip);
//     clip->SetPos(pos);
//     clip->SetRot(rot);
//     clip->AddAsset(clip_color);
//     // if (id % 2 == 0)
//     //     clip->AddAsset(chrono_types::make_shared<ChColorAsset>(ChColor(0.8f, 0.1f, 0.1f)));
//     // else
//     //     clip->AddAsset(chrono_types::make_shared<ChColorAsset>(ChColor(0.1f, 0.1f, 0.8f)));
    
// }

// void create_torus(ChSystemNSC& sys, std::shared_ptr<ChMaterialSurface> mat,
//                 const int id, ChVector<>& pos, const ChQuaternion<double>& rot,
//                 const double torus_r, const double torus_t, const int nseg, const double angle,
//                 const double mass, const ChVector<> inertia)
// {


//     auto torus_color = chrono_types::make_shared<ChColorAsset>(ChColor(0.8f, 0.1f, 0.1f));
//     auto init_ang_velo = ChVector<>(0.0, (20.0*M_PI), 0.0);
//     if (id % 2 != 0)
//     {
//         init_ang_velo = ChVector<>((20.0*M_PI), 0.0, 0.0*M_PI);
//         torus_color = chrono_types::make_shared<ChColorAsset>(ChColor(0.1f, 0.1f, 0.8f));
//     }
    

//     auto torus = std::make_shared<ChBody>();
// 	torus->SetIdentifier(id);
//     // clip->GetCollisionModel()->SetDefaultSuggestedEnvelope(0.001);
//     torus->GetCollisionModel()->SetEnvelope(0.015);
//     torus->GetCollisionModel()->SetDefaultSuggestedMargin(0.0001);
//     torus->SetCollide(true);
//     torus->GetCollisionModel()->ClearModel();

//      utils::AddTorusGeometry(torus.get(),
//                         mat,
//                         torus_r,
//                         torus_t,
//                         nseg, 
//                         angle,
//                         ChVector<>(0, 0, 0),
//                         ChQuaternion<>(1, 0, 0, 0),
//                         true);
    
//     torus->GetCollisionModel()->BuildModel();

//     torus->SetMass(mass);
//     torus->SetInertiaXX(inertia);
    
//     torus->SetWvel_par(init_ang_velo);
//     sys.Add(torus);
//     torus->SetPos(pos);
//     torus->SetRot(rot);
//     torus->AddAsset(torus_color);
// }

// void create_model(ChSystemNSC& mphysicalSystem) {
//     ChVector<> gravity(0, 0, 0);
//     mphysicalSystem.Set_G_acc(gravity);
//     SetChronoDataPath("/Users/farhad/work/cs_master/clip_entanglement/data/");
//     // SetChronoDataPath("/Users/farhad/work/cs_master/chrono/data/");
//     // make the shaking box
//     auto floorBody = chrono_types::make_shared<ChBodyEasyBox>(1., 0.2, 1.0, 7800, false, false);
//     floorBody->SetPos(ChVector<>(0, -0.3, 0));
//     floorBody->SetBodyFixed(true);
//     mphysicalSystem.Add(floorBody);

//     double box_x = 0.25;
//     double box_y = 0.25;
//     double box_z = 0.25;
//     double box_t = 0.01;
//     int box_main_id = 9999;
//     auto box_mat = chrono_types::make_shared<ChMaterialSurfaceNSC>();
//     box_mat->SetFriction(0.2);
//     // box_mat->SetCompliance(0.000001);
//     box_mat->SetRestitution(0.1);
//     auto box_pos = ChVector<>(0, 0, 0);
//     // create_box(mphysicalSystem, floorBody, box_mat, box_main_id,
//     //             box_pos, QUNIT, box_x, box_y, box_z, box_t, 7800);
//     // model_box(mphysicalSystem, floorBody, box_mat, 
//     //     box_main_id, box_pos, QUNIT,
//     //     box_x, box_y, box_z, box_t,
//     //     7800);

//     auto clip_mat = chrono_types::make_shared<ChMaterialSurfaceNSC>();
//     clip_mat->SetFriction(0.6f);
//     // clip_mat->SetRollingFriction(0.001);
//     // clip_mat->SetSpinningFriction(0.01);
//     // clip_mat->SetCompliance(0.0001);
//     // clip_mat->SetComplianceT(0.02);
//     // clip_mat->SetComplianceRolling(0.08);
//     // clip_mat->SetComplianceSpinning(0.08);
//     // clip_mat->SetDampingF(0.001);
//     clip_mat->SetRestitution(0.5f);

//     double clip_w = 0.015;
//     double clip_h = 0.04;
//     double clip_r = 0.001;
//     double clip_g = 0.004;
//     double cclip_mass = 0.002695;
//     auto cclip_iner_odd = ChVector<>(0.125176E-6, 0.556291236E-6, 0.680161092E-6);
//     auto cclip_iner_even = ChVector<>(0.125176E-6, 0.680161092E-6, 0.556291236E-6);
//     auto rot_odd = QUNIT; //Q_ROTATE_X_TO_Y; //ChQuaternion<double>(0.7565, 0, 0.6, 0.8);
//     auto rot_even = Q_ROTATE_Y_TO_Z; //Q_ROTATE_X_TO_Y; //ChQuaternion<double>(0.7565, 0, 0.6, 0.8);

//     int n_clips = 10;
//     for (int nc=0; nc<n_clips; ++nc)
//     {
//         auto pos = ChVector<>(0.75*nc*clip_h, 0, 0);
//         if (nc % 2 == 0)
//             create_closedclip(mphysicalSystem, clip_mat, nc, 
//                 pos, rot_even,
//                 clip_w, clip_h, clip_r,
//                 cclip_mass, cclip_iner_even);
//         else
//             create_closedclip(mphysicalSystem, clip_mat, nc, 
//                 pos, rot_odd,
//                 clip_w, clip_h, clip_r,
//                 cclip_mass, cclip_iner_odd);
//     }
//     // double torus_r = 0.01;
//     // double torus_t = 0.001;
//     // int torus_nseg = 20;
//     // double torus_angle = 360.;
//     // double torus_mass = 0.003;
//     // auto torus_iner = ChVector<>(0.12E-6, 0.7E-6, 0.7E-6);

//     // auto pos = ChVector<>(0, 0, 0);
//     // auto rot_yz = Q_ROTATE_Y_TO_Z;

//     // create_torus(mphysicalSystem, clip_mat, 1, 
//     //         pos, rot,
//     //         torus_r, torus_t, torus_nseg, torus_angle,
//     //         torus_mass, torus_iner);

//     // pos = ChVector<>(1.5*torus_r, 0, 0);
//     // create_torus(mphysicalSystem, clip_mat, 1, 
//     //         pos, rot_yz,
//     //         torus_r, torus_t, torus_nseg, torus_angle,
//     //         torus_mass, torus_iner);
    
//     // create_closedclip(mphysicalSystem, clip_mat, 1, 
//     //         pos, rot,
//     //         clip_w, clip_h, clip_r,
//     //         cclip_mass, cclip_iner);
//     // pos = ChVector<>(0.5*clip_h, 0, 0);
//     // create_closedclip(mphysicalSystem, clip_mat, 2, 
//     //         pos, rot_yz,
//     //         clip_w, clip_h, clip_r,
//     //         cclip_mass, cclip_iner_yz);
// }

// int main(int argc, char* argv[]) {
//     GetLog() << "Copyright (c) 2017 projectchrono.org\nChrono version: " << CHRONO_VERSION << "\n\n";

//     // Create a ChronoENGINE physical system
//     ChSystemNSC mphysicalSystem;

//     // Create the Irrlicht visualization (open the Irrlicht device,
//     // bind a simple user interface, etc. etc.)
//     ChIrrApp application(&mphysicalSystem, L"Clip contact", core::dimension2d<u32>(800, 600));

//     application.AddTypicalLights(core::vector3df(0.f, 0.f, 10.f), core::vector3df(30.f, 80.f, 60.f), 200, 10);
//     application.AddTypicalCamera(core::vector3df(0.0, .0, 0.1), core::vector3df(0.75*5*0.04, 0.0, 0));

//     create_model(mphysicalSystem);

//     application.AssetBindAll();
//     application.AssetUpdateAll();


//     auto solver = chrono_types::make_shared<ChSolverPSOR>();
//     solver->SetMaxIterations(200);
//     // solver->SetTolerance(0.0001);
//     solver->EnableWarmStart(true);
//     mphysicalSystem.SetSolver(solver);

//     double dT = 0.001;
//     double endT = 10.0;
//     int n_frames = 200;
//     int f_interval = (int) (endT/dT/n_frames);
//     std::vector<double> time_array;
//     std::vector<double> ke_array;
//     application.SetTimestep(dT);
//     // application.SetTryRealtime(true);
//     // application.SetPOVraySave(true);
//     // application.SetPOVraySaveInterval(5);

//     application.SetVideoframeSave(true);
//     application.SetVideoframeSaveInterval(10);

//     // application.GetSystem()->SetChTime(1.0);
//     // std::cout << "maxi recovery spped is: " << application.GetSystem()->GetMaxPenetrationRecoverySpeed() << std::endl;
//     application.GetSystem()->SetMaxPenetrationRecoverySpeed(0.03);
//     application.GetSystem()->SetMinBounceSpeed(0.1);
//     double start = std::clock();
//     while (application.GetDevice()->run()) {
//         // std::cout << "system stepsize is : " << application.GetSystem()->GetStep() << std::endl;
//         application.BeginScene(true, true, SColor(255, 140, 161, 192));
//         application.DrawAll();

//         application.DoStep();

//         application.EndScene();
//         time_array.push_back(application.GetSystem()->GetChTime());
//         // ke_array.push_back(application.GetSystem()->ComputeTotalKE());
//         // if (application.GetSystem()->GetChTime() > endT)
//         //     break;
//     }
//     double duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
//     std::cout << "simulation run time is : " << duration << std::endl;
//     return 0;
// }















// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Radu Serban
// =============================================================================
//
// Chrono::Multicore test program using NSC method for frictional contact.
//
// The model simulated here consists of a number of spherical objects falling
// in a fixed container.
//
// The global reference frame has Z up.
//
// If available, OpenGL is used for run-time rendering. Otherwise, the
// simulation is carried out for a pre-defined duration and output files are
// generated for post-processing with POV-Ray.
// =============================================================================

#include <cstdio>
#include <vector>
#include <cmath>

#include "chrono_multicore/physics/ChSystemMulticore.h"

#include "chrono/ChConfig.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/utils/ChUtilsInputOutput.h"

#ifdef CHRONO_OPENGL
#include "chrono_opengl/ChOpenGLWindow.h"
#endif

using namespace chrono;
using namespace chrono::collision;

// Tilt angle (about global Y axis) of the container.
double tilt_angle = 1 * CH_C_PI / 20;

// Number of balls: (2 * count_X + 1) * (2 * count_Y + 1)
int count_X = 2;
int count_Y = 2;

// -----------------------------------------------------------------------------
// Create a bin consisting of five boxes attached to the ground.
// -----------------------------------------------------------------------------
void AddContainer(ChSystemMulticoreNSC* sys) {
    // IDs for the two bodies
    int binId = -200;

    // Create a common material
    auto mat = chrono_types::make_shared<ChMaterialSurfaceNSC>();
    mat->SetFriction(0.4f);

    // Create the containing bin (4 x 4 x 1)
    auto bin = std::shared_ptr<ChBody>(sys->NewBody());
    bin->SetIdentifier(binId);
    bin->SetMass(1);
    bin->SetPos(ChVector<>(0, 0, 0));
    bin->SetRot(Q_from_AngY(tilt_angle));
    bin->SetCollide(true);
    bin->SetBodyFixed(true);

    ChVector<> hdim(2, 2, 0.5);
    double hthick = 0.1;

    bin->GetCollisionModel()->ClearModel();
    utils::AddBoxGeometry(bin.get(), mat, ChVector<>(hdim.x(), hdim.y(), hthick), ChVector<>(0, 0, -hthick));
    utils::AddBoxGeometry(bin.get(), mat, ChVector<>(hthick, hdim.y(), hdim.z()),
                          ChVector<>(-hdim.x() - hthick, 0, hdim.z()));
    utils::AddBoxGeometry(bin.get(), mat, ChVector<>(hthick, hdim.y(), hdim.z()),
                          ChVector<>(hdim.x() + hthick, 0, hdim.z()));
    utils::AddBoxGeometry(bin.get(), mat, ChVector<>(hdim.x(), hthick, hdim.z()),
                          ChVector<>(0, -hdim.y() - hthick, hdim.z()));
    utils::AddBoxGeometry(bin.get(), mat, ChVector<>(hdim.x(), hthick, hdim.z()),
                          ChVector<>(0, hdim.y() + hthick, hdim.z()));
    bin->GetCollisionModel()->BuildModel();

    sys->AddBody(bin);
}

// -----------------------------------------------------------------------------
// Create the falling spherical objects in a uniform rectangular grid.
// -----------------------------------------------------------------------------
void AddFallingBalls(ChSystemMulticore* sys) {
    // Common material
    auto ballMat = chrono_types::make_shared<ChMaterialSurfaceNSC>();
    ballMat->SetFriction(0.4f);

    // Create the falling balls
    int ballId = 0;
    double mass = 1;
    double radius = 0.15;
    ChVector<> inertia = (2.0 / 5.0) * mass * radius * radius * ChVector<>(1, 1, 1);

    for (int ix = -count_X; ix <= count_X; ix++) {
        for (int iy = -count_Y; iy <= count_Y; iy++) {
            ChVector<> pos(0.4 * ix, 0.4 * iy, 1);

            auto ball = std::shared_ptr<ChBody>(sys->NewBody());

            ball->SetIdentifier(ballId++);
            ball->SetMass(mass);
            ball->SetInertiaXX(inertia);
            ball->SetPos(pos);
            ball->SetRot(ChQuaternion<>(1, 0, 0, 0));
            ball->SetBodyFixed(false);
            ball->SetCollide(true);

            ball->GetCollisionModel()->ClearModel();
            utils::AddSphereGeometry(ball.get(), ballMat, radius);
            ball->GetCollisionModel()->BuildModel();

            sys->AddBody(ball);
        }
    }
}

// -----------------------------------------------------------------------------
// Create the system, specify simulation parameters, and run simulation loop.
// -----------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    GetLog() << "Copyright (c) 2017 projectchrono.org\nChrono version: " << CHRONO_VERSION << "\n\n";

    // Simulation parameters
    // ---------------------

    double gravity = 9.81;
    double time_step = 1e-3;

    uint max_iteration = 300;
    real tolerance = 1e-3;

    // Create system
    // -------------

    ChSystemMulticoreNSC msystem;

    // Set number of threads
    msystem.SetNumThreads(8);

    // Set gravitational acceleration
    msystem.Set_G_acc(ChVector<>(0, 0, -gravity));

    // Set solver parameters
    msystem.GetSettings()->solver.solver_mode = SolverMode::SLIDING;
    msystem.GetSettings()->solver.max_iteration_normal = max_iteration / 3;
    msystem.GetSettings()->solver.max_iteration_sliding = max_iteration / 3;
    msystem.GetSettings()->solver.max_iteration_spinning = 0;
    msystem.GetSettings()->solver.max_iteration_bilateral = max_iteration / 3;
    msystem.GetSettings()->solver.tolerance = tolerance;
    msystem.GetSettings()->solver.alpha = 0;
    msystem.GetSettings()->solver.contact_recovery_speed = 10000;
    msystem.ChangeSolverType(SolverType::APGDREF);
    // msystem.GetSettings()->collision.narrowphase_algorithm = ChNarrowphase::Algorithm::HYBRID;

    msystem.GetSettings()->collision.collision_envelope = 0.01;
    msystem.GetSettings()->collision.bins_per_axis = vec3(10, 10, 10);

    // Create the fixed and moving bodies
    // ----------------------------------

    AddContainer(&msystem);
    AddFallingBalls(&msystem);

// Perform the simulation
// ----------------------

// #ifdef CHRONO_OPENGL
//     opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();
//     gl_window.Initialize(1280, 720, "ballsNSC", &msystem);
//     gl_window.SetCamera(ChVector<>(0, -6, 0), ChVector<>(0, 0, 0), ChVector<>(0, 0, 1));
//     gl_window.SetRenderMode(opengl::WIREFRAME);

//     // Uncomment the following two lines for the OpenGL manager to automatically
//     // run the simulation in an infinite loop.
//     // gl_window.StartDrawLoop(time_step);
//     // return 0;

//     while (true) {
//         if (gl_window.Active()) {
//             gl_window.DoStepDynamics(time_step);
//             gl_window.Render();
//             ////if (gl_window.Running()) {
//             ////  msystem.CalculateContactForces();
//             ////  real3 frc = msystem.GetBodyContactForce(0);
//             ////  std::cout << frc.x << "  " << frc.y << "  " << frc.z << std::endl;
//             ////}
//         } else {
//             break;
//         }
//     }
// #else
    // Run simulation for specified time
    double time_end = 2;
    int num_steps = (int)std::ceil(time_end / time_step);
    double time = 0;

    for (int i = 0; i < num_steps; i++) {
        msystem.DoStepDynamics(time_step);
        time += time_step;
    }
#endif

    return 0;
}
