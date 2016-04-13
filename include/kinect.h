#ifndef KINECT_H
#define KINECT_H


//#include <XnOpenNI.h>
//#include <XnLog.h>
//#include <XnCppWrapper.h>
//#include <XnFPSCalculator.h>
//using namespace xn;


#include <opencv/cv.h>
#include <opencv/highgui.h>
using namespace cv;

#include <include/Eigen/Dense>

class kinect
{
public:

//    Context context;
//    ScriptNode scriptNode;
//    EnumerationErrors errors;
//    DepthGenerator depth;
//    DepthMetaData depthMD;

    /////////////////////////////
    // RGB Intrinsic Parameters
    float fx_rgb,fy_rgb,cx_rgb,cy_rgb;
    // RGB Distortion Parameters
    float k1_rgb,k2_rgb,p1_rgb,p2_rgb,k3_rgb;

    /////////////////////////////
    // Depth Intrinsic Parameters
    float fx_d,fy_d,cx_d,cy_d;
    // Depth Distortion Parameters
    float k1_d,k2_d,p1_d,p2_d,k3_d;


    // Position of the Plane Points of the depth camera
    Eigen::Vector2f * ptDepthPlanePointsMap;

    Eigen::Vector3f * pt3DPoints_DepthCamera;
    Eigen::Vector3f * pt3DPoints_RGBCamera;

    uint8_t *u8_RGBImgage;
    u_int16_t *u16_DImgage;

public:
    kinect();
    ~kinect();

    void depthCameraPlaneToDepthWorld(void);
    void depthWorldToRGBWorld(void);
    void initDepthMap(void);
    void projectToRGBCamera(void);

};

#endif // KINECT_H
