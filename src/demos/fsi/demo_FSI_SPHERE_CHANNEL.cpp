// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All right reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Author: Milad Rakhsha
// =============================================================================

// General Includes
#include <assert.h>
#include <limits.h>
#include <stdlib.h>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "chrono/physics/ChSystemSMC.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/utils/ChUtilsGenerators.h"
#include "chrono/utils/ChUtilsGeometry.h"
#include "chrono/utils/ChUtilsInputOutput.h"
#include "chrono_thirdparty/filesystem/path.h"
#include "chrono_thirdparty/filesystem/resolver.h"

// Chrono fsi includes
#include "chrono_fsi/ChDeviceUtils.cuh"
#include "chrono_fsi/ChFsiTypeConvert.h"
#include "chrono_fsi/ChSystemFsi.h"
#include "chrono_fsi/utils/ChUtilsGeneratorFsi.h"
#include "chrono_fsi/utils/ChUtilsJsonInput.h"
#include "chrono_fsi/utils/ChUtilsPrintSph.cuh"

// Chrono namespaces
using namespace chrono;
using namespace collision;
using namespace filesystem;

using std::cout;
using std::endl;
std::ofstream simParams;
typedef fsi::Real Real;

//----------------------------
// output directories and settings
//----------------------------
const std::string out_dir = GetChronoOutputPath() + "FSI_SPHERE_CHANNEL/";
std::string demo_dir;

bool save_output = true;

Real bxDim = 4.5;
Real byDim = 0.06;
Real bzDim = 0.8;

Real fxDim = bxDim;
Real fyDim = byDim;
Real fzDim = bzDim;
double cyl_length = 0.07;
double cyl_radius = .15;
ChVector<> cyl_pos = ChVector<>(0.0);

void ShowUsage() {
    cout << "usage: ./demo_FSI_SPHERE_CHANNEL <json_file>" << endl;
}
/// Forward declaration of a helper function
void SaveParaViewFilesMBD(fsi::ChSystemFsi& myFsiSystem,
                          ChSystemSMC& mphysicalSystem,
                          fsi::SimParams* paramsH,
                          int next_frame,
                          double mTime);
//------------------------------------------------------------------
// Create the objects of the MBD system. Rigid bodies, and if fsi, their
// bce representation are created and added to the systems
//------------------------------------------------------------------
void CreateSolidPhase(ChSystemSMC& mphysicalSystem, fsi::ChSystemFsi& myFsiSystem, fsi::SimParams* paramsH) {
    std::shared_ptr<ChMaterialSurfaceNSC> mat_g(new ChMaterialSurfaceNSC);
    auto mysurfmaterial = chrono_types::make_shared<ChMaterialSurfaceSMC>();

    // Set common material Properties
    mysurfmaterial->SetYoungModulus(6e4);
    mysurfmaterial->SetFriction(0.3f);
    mysurfmaterial->SetRestitution(0.2f);
    mysurfmaterial->SetAdhesion(0);
    // Ground body
    auto ground = chrono_types::make_shared<ChBody>();
    ground->SetIdentifier(-1);
    ground->SetBodyFixed(true);
    ground->SetCollide(true);
    ground->SetMaterialSurface(mysurfmaterial);
    ground->GetCollisionModel()->ClearModel();

    // Bottom wall
    ChVector<> sizeBottom(bxDim / 2, byDim / 2 + 2. * paramsH->HSML, 2 * paramsH->HSML); // walls should be taller?
    ChVector<> posBottom(0, 0, -2 * paramsH->HSML);
    ChVector<> posTop(0, 0, bzDim + 2 * paramsH->HSML);

    // left and right Wall
    ChVector<> size_YZ(2 * paramsH->HSML, byDim / 2 + 3 * paramsH->HSML, bzDim / 2);
    ChVector<> pos_xp(bxDim / 2 + paramsH->HSML, 0.0, bzDim / 2 + 1 * paramsH->HSML);
    ChVector<> pos_xn(-bxDim / 2 - 3 * paramsH->HSML, 0.0, bzDim / 2 + 1 * paramsH->HSML);

    // Front and back Wall
    ChVector<> size_XZ(bxDim / 2, 2 * paramsH->HSML, bzDim / 2);
    ChVector<> pos_yp(0, byDim / 2 + paramsH->HSML, bzDim / 2 + 1 * paramsH->HSML);
    ChVector<> pos_yn(0, -byDim / 2 - 3 * paramsH->HSML, bzDim / 2 + 1 * paramsH->HSML);

     utils::AddBoxGeometry(ground.get(), sizeBottom, posBottom, QUNIT, true);
    //    utils::AddBoxGeometry(ground.get(), size_YZ, pos_xp, QUNIT, true);
    //    utils::AddBoxGeometry(ground.get(), size_YZ, pos_xn, QUNIT, true);
   //     utils::AddBoxGeometry(ground.get(), size_XZ, pos_yp, QUNIT, true);
   //    utils::AddBoxGeometry(ground.get(), size_XZ, pos_yn, QUNIT, true);
    ground->GetCollisionModel()->BuildModel();
    mphysicalSystem.AddBody(ground);

    // Add the boundaries to the FSI system
    fsi::utils::AddBoxBce(myFsiSystem.GetDataManager(), paramsH, ground, posBottom, QUNIT, sizeBottom);
    fsi::utils::AddBoxBce(myFsiSystem.GetDataManager(), paramsH, ground, posTop, QUNIT, sizeBottom);
    //    fsi::utils::AddBoxBce(myFsiSystem.GetDataManager(), paramsH, ground, pos_yp, QUNIT, size_XZ,
    //    13); 
    //	fsi::utils::AddBoxBce(myFsiSystem.GetDataManager(), paramsH, ground, pos_yn, QUNIT,
    //    size_XZ, 13);
    std::cout << "Inserted side boundaries" << std::endl;	
	// Now add the cylinder!!!!!!!!!!!!!!
    ChQuaternion<> cyl_rot = chrono::QUNIT;

    int ii;
    for (ii=0; ii<5; ++ii) {
    ChVector<> this_cyl = ChVector<>(0.0);
    this_cyl = cyl_pos + ChVector<> (0.7*ii, 0.0, 0.0);
    auto body = std::make_shared<ChBody>();
    // body->SetIdentifier(-1);
    body->SetBodyFixed(true);
    body->SetCollide(true);
    body->SetMaterialSurface(mysurfmaterial);

    body->SetPos(this_cyl + ChVector<>(0, 0.5 * paramsH->HSML, 0.0));

    double sphereRad = 0.3;
    double volume = utils::CalcSphereVolume(sphereRad);
    // ChVector<> gyration = utils::CalcSphereGyration(sphereRad).Get_Diag();
    double density = paramsH->rho0;
    double mass = density * volume;
    body->SetMass(mass);
    // body->SetInertiaXX(mass * gyration);
    body->GetCollisionModel()->ClearModel();
    utils::AddSphereGeometry(body.get(), sphereRad);
    body->GetCollisionModel()->BuildModel();
    int numRigidObjects = mphysicalSystem.Get_bodylist().size();
    mphysicalSystem.AddBody(body);

    chrono::fsi::utils::AddCylinderBce(myFsiSystem.GetDataManager(), paramsH, body, ChVector<>(0, 0, 0),
                                       ChQuaternion<>(1, 0, 0, 0), cyl_radius, cyl_length, paramsH->HSML, false);
    myFsiSystem.AddFsiBody(body);
    }
	
}

// =============================================================================
int main(int argc, char* argv[]) {
    ChSystemSMC mphysicalSystem;
    fsi::ChSystemFsi myFsiSystem(&mphysicalSystem, true, fsi::ChFluidDynamics::Integrator::IISPH);
    fsi::SimParams* paramsH = myFsiSystem.GetSimParams();

    // Use the default input file or you may enter your input parameters as a command line argument
    std::string inputJson = GetChronoDataFile("fsi/input_json/demo_FSI_SPHERE_CHANNEL.json");
    if (argc < 2 || chrono::fsi::utils::ParseJSON(argv[1], paramsH, chrono::fsi::mR3(bxDim, byDim, bzDim)) == false) {
        ShowUsage();
        return 1;
    }
    myFsiSystem.SetFluidSystemLinearSolver(paramsH->LinearSolver);

    // This is an important stage for neighbor search stage
    // If there are periodic boundary conditions in the problem, one should be careful to make sure
    // that the placements of the particles with respect to the c_min and c_max of the domain does not create a gap
    // nor it creates overlap between the particles on the period patches.
    // Remember that particles close to the one side of the period patch interact with the particles close to the other
    // side.
    paramsH->cMin = fsi::mR3(-bxDim / 2 - paramsH->HSML / 2, -byDim / 2 - paramsH->HSML, -bzDim * 2);
    paramsH->cMax = fsi::mR3(bxDim / 2 + paramsH->HSML / 2, byDim / 2 + paramsH->HSML, bzDim * 2);
    fsi::utils::FinalizeDomainCreating(paramsH);

    if (strcmp(paramsH->out_name, "Undefined") == 0)
        demo_dir = out_dir + "Paraview/";
    else
        demo_dir = out_dir + paramsH->out_name + "/";

    if (!filesystem::create_directory(filesystem::path(out_dir))) {
        std::cout << "Error creating directory " << out_dir << std::endl;
        return 1;
    }
    if (!filesystem::create_directory(filesystem::path(demo_dir))) {
        std::cout << "Error creating directory " << out_dir << std::endl;
        return 1;
    }
    const std::string removeFiles = (std::string("rm ") + demo_dir + "/* ");
    system(removeFiles.c_str());
    if (argc == 2) {
        const std::string copyInputs = (std::string("cp ") + argv[1] + " " + demo_dir);
        system(copyInputs.c_str());
    }

    // ******************************* Create Fluid region ****************************************
    Real initSpace0 = paramsH->MULT_INITSPACE * paramsH->HSML;
    Real initSpace1 = initSpace0;
    utils::GridSampler<> sampler(initSpace0);
    fsi::Real3 boxCenter = fsi::mR3(-bxDim / 2 + fxDim / 2, 0 * paramsH->HSML, fzDim / 2 + 1 * paramsH->HSML);
    fsi::Real3 boxHalfDim = fsi::mR3(fxDim / 2, fyDim / 2, fzDim / 2);
	
	cyl_pos = ChVector<>(boxCenter.x - (1.25/4.5) * fxDim, boxCenter.y, boxCenter.z);
    //    mbody->SetPos(cyl_pos);
    //
    //    cyl_pos = ChVector<>(-initSpace1 / 2, 0, 0.2 - initSpace1 / 2);
    //    mbody->SetPos(cyl_pos);

    int numhelperMarkers = myFsiSystem.GetDataManager()->sphMarkersH.posRadH.size();

    // Ghost particles for Variable Resolution SPH
    int numGhostMarkers = 0;
	
    utils::Generator::PointVector points = sampler.SampleBox(fsi::ChFsiTypeConvert::Real3ToChVector(boxCenter),
                                                             fsi::ChFsiTypeConvert::Real3ToChVector(boxHalfDim));
    int numPart = points.size();
	
    std::vector<fsi::Real3> particles_position; // used to mark particles to remove
    if (argc > 2) {
        fsi::Real3 this_particle;

        std::fstream fin(out_dir + "BCE_Rigid.csv");
        if (!fin.good())
            throw ChException("ERROR opening Mesh file: BCE_Rigid.csv \n");

        std::string line;
        getline(fin, line);

        while (getline(fin, line)) {
            int ntoken = 0;
            std::string token;
            std::istringstream ss(line);
            while (std::getline(ss, token, ' ') && ntoken < 4) {
                std::istringstream stoken(token);
                if (ntoken == 0)
                    stoken >> this_particle.x;
                if (ntoken == 1)
                    stoken >> this_particle.y;
                if (ntoken == 2)
                    stoken >> this_particle.z;
                ++ntoken;
            }
            particles_position.push_back(this_particle);
        }
    }

    int BCE_RIGID_SIZE = particles_position.size();
    for (int i = 0; i < numPart; i++) {
        particles_position.push_back(fsi::ChFsiTypeConvert::ChVectorToReal3(points[i]));
    }

    std::cout << "Set Removal points from: BCE_Rigid.csv" << std::endl;
    std::cout << "Searching among " << particles_position.size() << "BCE points" << std::endl;

    int numremove = 0;
    for (int i = 0; i < numPart; i++) {
        bool removeThis = false;
        Real h = paramsH->HSML;
        fsi::Real4 p = fsi::mR4(points[i].x(), points[i].y(), points[i].z(), h);
        for (int remove = 0; remove < BCE_RIGID_SIZE; remove++) {
            double dist = length(particles_position[remove] - mR3(p));
            if (dist < 0.9 * initSpace1) {
                removeThis = true;
                break;
            }
        }
        if (!removeThis) {
			myFsiSystem.GetDataManager()->AddSphMarker(fsi::mR4(points[i].x(), points[i].y(), points[i].z(), paramsH->HSML),
													   fsi::mR3(1e-10),
													   fsi::mR4(paramsH->rho0, paramsH->BASEPRES, paramsH->mu0, -1));
        } else
            numremove++;
    }
	std::cout << "Removed " << numremove << " Fluid particles from the simulation" << std::endl;

    particles_position.clear();

    int numPhases = myFsiSystem.GetDataManager()->fsiGeneralData.referenceArray.size();
	int numFluidPart = myFsiSystem.GetDataManager()->sphMarkersH.posRadH.size();
    std::cout << "Final fluid size:" << numFluidPart << std::endl;

    if (numPhases != 0) {
        std::cout << "Error! numPhases is wrong, thrown from main\n" << std::endl;
        std::cin.get();
        return -1;
    } else {
	std::cout << "Now start to push_back " << std::endl;
        //myFsiSystem.GetDataManager()->fsiGeneralData.referenceArray.push_back(mI4(0, numPart, -1, -1));
        //myFsiSystem.GetDataManager()->fsiGeneralData.referenceArray.push_back(mI4(numPart, numPart, 0, 0));
        myFsiSystem.GetDataManager()->fsiGeneralData.referenceArray.push_back(
        mI4(numhelperMarkers + numGhostMarkers, numFluidPart, -1, -1));
	}
    std::cout << "Push_back done" << std::endl;
	
	printf("Ref size=%d\n", myFsiSystem.GetDataManager()->fsiGeneralData.referenceArray.size());
    ChVector<> gravity = ChVector<>(paramsH->gravity.x, paramsH->gravity.y, paramsH->gravity.z);

    /// Create the objects of the Chsystem (Solid phases) using something like this function
    CreateSolidPhase(mphysicalSystem, myFsiSystem, paramsH);

    /// You must call the finalize before you start the simulation
    myFsiSystem.Finalize();
	
	printf("\n\n sphMarkersH end%d, numAllMarkers is %d \n\n\n",
           myFsiSystem.GetDataManager()->sphMarkersH.posRadH.size(),
           myFsiSystem.GetDataManager()->numObjects.numAllMarkers);

    if (myFsiSystem.GetDataManager()->sphMarkersH.posRadH.size() !=
        myFsiSystem.GetDataManager()->numObjects.numAllMarkers) {
        printf("\n\n\n\n Error! (2) numObjects is not set correctly \n %d, %d \n\n\n",
               myFsiSystem.GetDataManager()->sphMarkersH.posRadH.size(),
               myFsiSystem.GetDataManager()->numObjects.numAllMarkers);
        return -1;
    }
	

	// Time to solve
    double mTime = 0;
    int stepEnd = int(paramsH->tFinal / paramsH->dT);
    stepEnd = 1000000;
    std::vector<std::vector<double>> vCoor;
    std::vector<std::vector<int>> faces;

    SaveParaViewFilesMBD(myFsiSystem, mphysicalSystem, paramsH, 0, mTime);
	
	const std::string copyInitials =
        (std::string("cp ") + demo_dir + std::string("BCE_Rigid0.csv ") + out_dir + std::string("BCE_Rigid.csv "));
    cout << copyInitials;
    system(copyInitials.c_str());
    if (argc <= 2) {
        printf("now please run with an input argument\n");
        return 0;
    }
	
	
    Real time = 0;
    Real Global_max_dT = paramsH->dT_Max;
    for (int tStep = 0; tStep < stepEnd + 1; tStep++) {
        printf("\nstep : %d, time= : %f (s) \n", tStep, time);
        double frame_time = 1.0 / paramsH->out_fps;
        int next_frame = std::floor((time + 1e-6) / frame_time) + 1;
        double next_frame_time = next_frame * frame_time;
        double max_allowable_dt = next_frame_time - time;
        if (max_allowable_dt > 1e-6)
            paramsH->dT_Max = std::min(Global_max_dT, max_allowable_dt);
        else
            paramsH->dT_Max = Global_max_dT;

        printf("next_frame is:%d,  max dt is set to %f\n", next_frame, paramsH->dT_Max);
        myFsiSystem.DoStepDynamics_FSI();
        time += paramsH->dT;
        SaveParaViewFilesMBD(myFsiSystem, mphysicalSystem, paramsH, next_frame, time);
	if (time>paramsH->tFinal)
	   break;
    }

    return 0;
}

//------------------------------------------------------------------
// Function to save the povray files of the MBD
//------------------------------------------------------------------
void SaveParaViewFilesMBD(fsi::ChSystemFsi& myFsiSystem,
                          ChSystemSMC& mphysicalSystem,
                          fsi::SimParams* paramsH,
                          int next_frame,
                          double mTime) {
    static double exec_time;

    int out_steps = std::ceil((1.0 / paramsH->dT) / paramsH->out_fps);
    static int out_frame = 0;
    double frame_time = 1.0 / paramsH->out_fps;
    if (save_output && std::abs(mTime - (next_frame)*frame_time) < 1e-7) {
        fsi::utils::PrintToFile(
            myFsiSystem.GetDataManager()->sphMarkersD2.posRadD, myFsiSystem.GetDataManager()->sphMarkersD2.velMasD,
            myFsiSystem.GetDataManager()->sphMarkersD2.rhoPresMuD,
            myFsiSystem.GetDataManager()->fsiGeneralData.referenceArray, thrust::host_vector<int4>(), demo_dir, true);
        cout << "-------------------------------------\n" << endl;
        cout << "             Output frame:   " << next_frame << endl;
        cout << "             Time:           " << mTime << endl;
        cout << "-------------------------------------\n" << endl;
        out_frame++;
    }
}
