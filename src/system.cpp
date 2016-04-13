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
    mbDone = false;

//    SE3<> Myposition = SE3<>::exp(makeVector(1,1.5,2.5,1,2,3));

    Eigen::VectorXf mu(6);
    mu(0) = 1;
    mu(1) = 1.5;
    mu(2) = 2.5;
    mu(3) = 1;
    mu(4) = 2;
    mu(5) = 3;

    Eigen::Matrix4f Myposition = exp(mu);

    SourceKinect = new kinect();
    std::cout << "Myposition " << std::endl<< Myposition<< std::endl;


}


//////////////////////////////////
/// \brief System::Run
///
void System::Run()
{
    static bool bFirstFrame = true;


    while(!mbDone)
    {
        glutMainLoopEvent();
        // We use two versions of each video frame:
        // One black and white (for processing by the tracker etc)
        // and one RGB, for drawing.

        // Grab new video frame...
        GetAndFillFrameBWandRGB();

        if(bFirstFrame)
        {
            //	  mpARDriver->Init();
            bFirstFrame = false;
        }

        GLwin.SetupViewport();
        GLwin.SetupVideoOrtho();


//        SE3<> Myposition = SE3<>::exp(makeVector(0,0,0,1,2,3))*SE3<>::exp(makeVector(0,0,-1,0,0,0));

        Eigen::VectorXf mu1(6);
        mu1(0) = 0;
        mu1(1) = 0;
        mu1(2) = 0;
        mu1(3) = 1;
        mu1(4) = 2;
        mu1(5) = 3;

        Eigen::VectorXf mu2(6);
        mu2(0) = 0;
        mu2(1) = 0;
        mu2(2) = -1;
        mu2(3) = 0;
        mu2(4) = 0;
        mu2(5) = 0;

        Eigen::Matrix4f MyPosition = exp(mu1)*exp(mu2);
        mpMapViewer->DrawMap(MyPosition,SourceKinect);

        if(GLwin.mbUserPressedSpacebar)
        {
            GLwin.mbUserPressedSpacebar = false;
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





