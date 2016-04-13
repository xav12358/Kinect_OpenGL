
#ifndef __GL_WINDOW_2_H
#define __GL_WINDOW_2_H

//#include <cvd/glwindow.h>
//#include <TooN/TooN.h>

//class GLWindowMenu;

#include <opencv/cv.h>
#include <opencv/highgui.h>

#include <GL/gl.h>
#include <GL/glut.h>
#include <GL/freeglut.h>


#include <include/Eigen/Dense>

class GLWindow2
{

    bool stateLeft,stateRight,stateMiddle,stateMiddleUP,stateMiddleDOWN;
public:
    bool mbUserPressedSpacebar;

    GLWindow2(Eigen::Vector2f irSize);

    // Signal handlers
    void GLMouseHandler( int button, int state, int x, int y );
    void GLMouseMotion(int x,int y);
    void GLkeyboardHandler( unsigned char key, int x, int y );
    void GLReshape(int w, int h);
    void GLDisplay();



    //  // Some OpenGL helpers:
    void SetupViewport();
    void SetupVideoOrtho();
    void SetupUnitOrtho();
    void SetupWindowOrtho();
    void SetupVideoRasterPosAndZoom();

    // Map viewer mouse interface:
    std::pair<Eigen::VectorXf, Eigen::VectorXf > GetMousePoseUpdate();


protected:
    // User interface menus:
    Eigen::Vector2f mirVideoSize;   // The size of the source video material.

    Eigen::Vector2f mirLastMousePos;

    // Storage for map viewer updates:
    Eigen::VectorXf mvMCPoseUpdate;
    Eigen::VectorXf mvLeftPoseUpdate;


};








#endif
