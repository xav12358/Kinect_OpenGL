#include "include/kinect.h"

//---------------------------------------------------------------------------
// Defines
//---------------------------------------------------------------------------
#define SAMPLE_XML_PATH "/home/xavier/kinect/OpenNI/Samples/Config/SamplesConfig.xml"

//---------------------------------------------------------------------------
// Macros
//---------------------------------------------------------------------------
#define CHECK_RC(rc, what) \
    if (rc != XN_STATUS_OK)\
{\
    printf("%s failed: %s\n", what, xnGetStatusString(rc));\
    }
//                            return rc;  \
//                             }

#include "iostream"
using namespace std;

kinect::kinect()
{

    // RGB Intrinsic Parameters
    fx_rgb = 5.1885790117450188e+02;
    fy_rgb = 5.1946961112127485e+02;
    cx_rgb = 3.2558244941119034e+02;
    cy_rgb = 2.5373616633400465e+02;

    // RGB Distortion Parameters
    k1_rgb =  2.0796615318809061e-01;
    k2_rgb = -5.8613825163911781e-01;
    p1_rgb = 7.2231363135888329e-04;
    p2_rgb = 1.0479627195765181e-03;
    k3_rgb = 4.9856986684705107e-01;

    // Depth Intrinsic Parameters
    fx_d = 5.8262448167737955e+02;
    fy_d = 5.8269103270988637e+02;
    cx_d = 3.1304475870804731e+02;
    cy_d = 2.3844389626620386e+02;

    // Depth Distortion Parameters
    k1_d = -9.9897236553084481e-02;
    k2_d = 3.9065324602765344e-01;
    p1_d = 1.9290592870229277e-03;
    p2_d = -1.9422022475975055e-03;
    k3_d = -5.1031725053400578e-01;

    //    XnStatus nRetVal = XN_STATUS_OK;
    //    nRetVal = context.InitFromXmlFile(SAMPLE_XML_PATH, scriptNode, &errors);
    //    nRetVal = context.FindExistingNode(XN_NODE_TYPE_DEPTH, depth);
    //    CHECK_RC(nRetVal, "Find depth generator");

    //    XnFPSData xnFPS;
    //    nRetVal = xnFPSInit(&xnFPS, 180);
    //    CHECK_RC(nRetVal, "FPS Init");

    ptDepthPlanePointsMap     = (Eigen::Vector2f *)malloc(640*480*sizeof(Eigen::Vector2f));
    pt3DPoints_DepthCamera    = (Eigen::Vector3f *)malloc(640*480*sizeof(Eigen::Vector3f));
    pt3DPoints_RGBCamera      = (Eigen::Vector3f *)malloc(640*480*sizeof(Eigen::Vector3f));

    //Create the map for the depth camera
    initDepthMap();

    cv::Mat depthCam =  imread("/home/xavier/Bureau/basements/basement_0001b/d-1316653651.126315-1259733252.pgm", 1);   // Read the file s
    cv::Mat depthCam2 =  imread("/home/xavier/Bureau/basements/basement_0001b/d-1316653651.126315-1259733252.pgm", 0);   // Read the file s

    std::cout << "image " << depthCam.cols << " " <<depthCam.rows << std::endl;
    std::cout << "image " << depthCam.channels() << " " <<  (int)(CV_16U) << " " << std::endl;

    cv::imshow("rrr",depthCam);
    cv::imshow("rrr2",depthCam2);

    cv::waitKey(-1);
}

kinect::~kinect()
{
    //    depth.Release();
    //    scriptNode.Release();
    //    context.Release();
}

/////////////////////////////////////////
/// \brief kinect::depthCameraPlaneTo3DCamera
///        Calculate the 3D position of the points in depth camera coordinate
void kinect::depthCameraPlaneToDepthWorld()
{
    for(int y=0; y < 480 ;y++)
    {
        for(int x = 0; x < 640 ;x++)
        {
            float Z_depth = u16_DImgage[x+y*640]/100.0;
            if(Z_depth !=0)
            {
                (pt3DPoints_DepthCamera[x+y*640])(0,0) = (ptDepthPlanePointsMap[x+y*640])(0,0)*Z_depth;
                (pt3DPoints_DepthCamera[x+y*640])(1,0) = (ptDepthPlanePointsMap[x+y*640])(1,0)*Z_depth;
                (pt3DPoints_DepthCamera[x+y*640])(2,0) = Z_depth;
            }
        }
    }
}


///////////////////////////////////////
/// \brief kinect::depthWorldToRGBWorld
///         Transform the 3D depth camera coordinate to RGB coordinate
void kinect::depthWorldToRGBWorld(void)
{
    for(int y=0; y < 480 ;y++)
    {
        for(int x = 0; x < 640 ;x++)
        {
            float Z_depth = (pt3DPoints_DepthCamera[x+y*640])(2,0);
            if(Z_depth !=0)
            {
                (pt3DPoints_DepthCamera[x+y*640])(0,0);
                (pt3DPoints_DepthCamera[x+y*640])(1,0);
                (pt3DPoints_DepthCamera[x+y*640])(2,0);
            }
        }
    }
}


////////////////////////////////////////////////////////////////////////////
/// \brief kinect::initDepthMap : Generate point coordonate in camera plane
///
void kinect::initDepthMap(void)
{
    float ifx = 1./fx_d;
    float ify = 1./fy_d;

    float xp,yp,x0,y0;
    double icdist,r2,deltaX,deltaY;
    int iiters = 10;
    for(int y=0; y < 480 ;y++)
    {
        for(int x = 0; x < 640 ;x++)
        {
            float xp = (x - cx_d)*ifx;
            float yp = (y - cy_d)*ify;

            float x0 = xp;
            float y0 = yp;

            // Compensate distortion iteratively
            for(int  j = 0; j < iiters; j++ )
            {
                r2 = xp*xp + yp*yp;
                icdist = (1 + ((k3_d*r2 + k2_d)*r2 + k1_d)*r2);///(1 + ((k[4]*r2 + k[1])*r2 + k[0])*r2);
                deltaX = 2*p1_d*xp*yp; + p2_d*(r2 + 2*xp*xp);//+ k[8]*r2+k[9]*r2*r2;
                deltaY = p1_d*(r2 + 2*yp*yp) + 2*p2_d*xp*yp;//+ k[10]*r2+k[11]*r2*r2;
                xp = (x0 - deltaX)*icdist;
                yp = (y0 - deltaY)*icdist;
            }
            (ptDepthPlanePointsMap[x+y*640])(0) = xp;
            (ptDepthPlanePointsMap[x+y*640])(1) = yp;
        }
    }
}

