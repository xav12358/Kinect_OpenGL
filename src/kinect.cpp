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
    


    Rotation(0,0) = -9.9997798940829263e-01;
    Rotation(0,1) = -5.0518419386157446e-03;
    Rotation(0,2) = -4.3011152014118693e-03;
    Rotation(1,0) = 5.0359919480810989e-03;
    Rotation(1,1) = -9.9998051861143999e-01;
    Rotation(1,2) = 3.6879781309514218e-03;
    Rotation(2,0) = 4.3196624923060242e-03;
    Rotation(2,1) = -3.6662365748484798e-03;
    Rotation(2,2) = -9.9998394948385538e-01;


    std::cout << "rotation " << Rotation << std::endl;
    std::cout << "Translation " << Translation << std::endl;

    Rotation.inverse();

    Translation(0) = 2.5031875059141302e-02;
    Translation(1) = 6.6238747008330102e-04;
    Translation(2) = -2.9342312935846411e-04;

    std::cout << "rotation " << Rotation << std::endl;
    std::cout << "Translation " << Translation << std::endl;

    se3DepthCamToRGBCam = Eigen::Matrix4f::Identity();
    se3DepthCamToRGBCam.block(0,0,3,3) = Rotation.inverse();
    se3DepthCamToRGBCam.block(0,3,3,1) = Translation;
    std::cout << "se3DepthCamToRGBCam " << se3DepthCamToRGBCam << std::endl;


    //    XnStatus nRetVal = XN_STATUS_OK;
    //    nRetVal = context.InitFromXmlFile(SAMPLE_XML_PATH, scriptNode, &errors);
    //    nRetVal = context.FindExistingNode(XN_NODE_TYPE_DEPTH, depth);
    //    CHECK_RC(nRetVal, "Find depth generator");
    
    //    XnFPSData xnFPS;
    //    nRetVal = xnFPSInit(&xnFPS, 180);
    //    CHECK_RC(nRetVal, "FPS Init");
    
    ptDepthPlanePointsMap     = (Eigen::Vector2f *)malloc(640*480*sizeof(Eigen::Vector2f));
    ptRGBPlanePointsMap       = (Eigen::Vector2f *)malloc(640*480*sizeof(Eigen::Vector2f));
    pt3DPoints_DepthCamera    = (Eigen::Vector3f *)malloc(640*480*sizeof(Eigen::Vector3f));
    pt3DPoints_RGBCamera      = (Eigen::Vector3f *)malloc(640*480*sizeof(Eigen::Vector3f));
    
    //Create the map for the depth camera
    initDepthMap();
    
    //    cv::Mat depthCam =  imread("/home/xavier/Bureau/basements/basement_0001b/d-1316653651.126315-1259733252.pgm", 1);   // Read the file s
    //    cv::Mat depthCam2 =  imread("/home/xavier/Bureau/basements/basement_0001b/d-1316653651.126315-1259733252.pgm", 0);   // Read the file s
    
    //    std::cout << "image " << depthCam.cols << " " <<depthCam.rows << std::endl;
    //    std::cout << "image " << depthCam.channels() << " " <<  (int)(CV_16U) << " " << std::endl;
    
    //    cv::imshow("rrr",depthCam);
    //    cv::imshow("rrr2",depthCam2);
    
    //    cv::waitKey(-1);
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
    
    //    cv::Mat depthCam2 =  imread("/home/xavier/Bureau/basements/basement_0001b/d-1316653651.126315-1259733252.pgm", -1);   // Read the file s
    cv::Mat depthCam2 =  imread("/home/xavier/Bureau/RGBD\ Dataset/living_room_traj0_frei_png/depth/0.png",-1);   // Read the file s
    u16_DImgage = (uint16_t*)(depthCam2.data);
    
    min= 99999; max =0;
    for(int y=0; y < 480 ;y++)
    {
        for(int x = 0; x < 640 ;x++)
        {
            float Z_depth = u16_DImgage[x+y*640]/10000.0;
            
            if(min>Z_depth)
                min=Z_depth;
            
            if(max<Z_depth)
                max=Z_depth;
            
            if(Z_depth !=0)
            {
//                std::cout << "depth " << u16_DImgage[x+y*640] << std::endl;

                (pt3DPoints_DepthCamera[x+y*640])(0) = (ptDepthPlanePointsMap[x+y*640])(0)*Z_depth;
                (pt3DPoints_DepthCamera[x+y*640])(1) = (ptDepthPlanePointsMap[x+y*640])(1)*Z_depth;
                (pt3DPoints_DepthCamera[x+y*640])(2) = Z_depth;
            }
        }
    }
}


///////////////////////////////////////
/// \brief kinect::depthWorldToRGBWorld
///         Transform the 3D depth camera coordinate to RGB coordinate
void kinect::depthWorldToRGBWorld(void)
{

    cv::Mat rgbCam      =  imread("/home/xavier/Bureau/RGBD\ Dataset/living_room_traj0_frei_png/rgb/0.png",0);   // Read the file s
//    cv::Mat depthCam    =  imread("/home/xavier/Bureau/RGBD\ Dataset/living_room_traj0_frei_png/depth/0.png",0);   // Read the file s

    Eigen::Vector4f VectorXYZ_RGB;
    float xp,yp;
    float xpp,ypp;
    float rp2;
    float u,v;
    for(int y=0; y < 480 ;y++)
    {
        for(int x = 0; x < 640 ;x++)
        {
            float Z_depth = (pt3DPoints_DepthCamera[x+y*640])(2);
            if(Z_depth !=0)
            {
                Eigen::Vector4f XYZW;
                XYZW(0) = pt3DPoints_DepthCamera[x+y*640](0);
                XYZW(1) = pt3DPoints_DepthCamera[x+y*640](1);
                XYZW(2) = pt3DPoints_DepthCamera[x+y*640](2);
                XYZW(3) = 1;
                VectorXYZ_RGB = se3DepthCamToRGBCam * XYZW;
//                VectorXYZ_RGB = VectorXYZ_RGB/VectorXYZ_RGB(3);
                xp = VectorXYZ_RGB(0)/VectorXYZ_RGB(2);
                yp = VectorXYZ_RGB(1)/VectorXYZ_RGB(2);
                rp2 = (xp*xp + yp*yp);

                float xp_yp = xp*yp;
                xpp = xp * (1 + (k1_rgb + (k2_rgb +k3_rgb *rp2)*rp2) *rp2) + 2*p1_rgb * xp_yp + p2_rgb *(rp2 +2*xp*xp);
                ypp = yp * (1 + (k1_rgb + (k2_rgb +k3_rgb *rp2)*rp2) *rp2) + p1_rgb * (rp2 +2*yp*yp) + 2*p2_rgb *xp_yp;

                u = (int)(fx_rgb*xpp + cx_rgb);
                v = (int)(fy_rgb*ypp + cy_rgb);

                if(u>0 && u<640 && v>0 && v<480)
                {
//                    std::cout << "u " << u << " v " << v  << " val  " << (int)rgbCam.data[(int)(u+v*640)] << std::endl;
//                    std::cout << "xpp " << xpp << " xp " << xp  << " val  " << (int)rgbCam.data[(int)(u+v*640)] << std::endl;
//                    std::cout << " pt3DPoints_DepthCamera[x+y*640] " << pt3DPoints_DepthCamera[x+y*640] << std::endl;
//                    std::cout << "XYZW " << XYZW << std::endl;
                     ptRGBPlanePointsMap[x+y*640](0) = (float)(rgbCam.data[(int)(u+v*640)]);//u;
//                    ptRGBPlanePointsMap[x+y*640](1) = v;

//                    cv::circle(rgbCam,cv::Point(u,v),2,cv::Scalar(255,255,255),1);
//                    cv::circle(depthCam,cv::Point(x,y),2,cv::Scalar(255,255,255),1);
//                    cv::imshow("rgbCam",rgbCam);
//                    cv::imshow("depthCam",depthCam);
//                    cv::waitKey(-1);
                }
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

