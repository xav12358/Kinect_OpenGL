#include "eigen_func.h"
#include "include/mapViewer.h"
#include "include/system.h"

#include <stdlib.h>
#include <iostream>
using namespace std;

#include <TooN/TooN.h>
#include <TooN/se3.h>
#include <TooN/Lapack_Cholesky.h>
using namespace TooN;


//////////////////////////////////////////
/// \brief System::System
/// \param GL
///
System::System(GLWindow2 &GL):
    GLwin(GL)
{

////    SE3<> Myposition = SE3<>::exp(makeVector(1,1.5,2.5,1,2,3));
//    Eigen::VectorXf mu(6);
//    mu(0) = 1;
//    mu(1) = 1.5;
//    mu(2) = 2.5;
//    mu(3) = 1;
//    mu(4) = 2;
//    mu(5) = 3;
//    Eigen::Matrix4f Myposition = exp(mu);
    std::cout << "System0 " << std::endl;

    SourceKinect = new kinect();
    std::cout << "System1 " << std::endl;
    mpMapViewer = new MapViewer(GLwin);
    std::cout << "System2 " << std::endl;
    mbDone = false;
}


//////////////////////////////////
/// \brief System::Run
///
void System::Run()
{
    static bool bFirstFrame = true;


    while(!mbDone)
    {
//        std::cout << "Run " << std::endl;
        glutMainLoopEvent();
        // We use two versions of each video frame:
        // One black and white (for processing by the tracker etc)
        // and one RGB, for drawing.

        // Grab new video frame...
        GetAndFillFrameBWandRGB();

        if(bFirstFrame)
        {
            bFirstFrame = false;
        }

        GLwin.SetupViewport();
        GLwin.SetupVideoOrtho();

        Eigen::VectorXf mu1(6);
        mu1(0) = 0;
        mu1(1) = 0;
        mu1(2) = 0;
        mu1(3) = 0;
        mu1(4) = 0;
        mu1(5) = 0;

        Eigen::VectorXf mu2(6);
        mu2(0) = 0;
        mu2(1) = 0;
        mu2(2) = 0;
        mu2(3) = 3.14/2;
        mu2(4) = 0;
        mu2(5) = 0;

        SourceKinect->depthCameraPlaneToDepthWorld();
        SourceKinect->depthWorldToRGBWorld();

        Eigen::Matrix4f MyPosition = exp(mu1)*exp(mu2);
        mpMapViewer->DrawMap(MyPosition,SourceKinect);
        if(GLwin.mbUserPressedSpacebar)
        {
//            GLwin.mbUserPressedSpacebar = false;
        }
    }
}


////////////////////////////////////////
/// \brief System::GetAndFillFrameBWandRGB
///
void System::GetAndFillFrameBWandRGB()
{

#if (ISSource ==1)
    // Grab the BW frame

#else

#endif

}





