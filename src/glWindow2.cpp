#include "include/glWindow2.h"

#include <stdlib.h>
#include <iostream>
using namespace std;

GLWindow2 gGLWin(Eigen::Vector2f(640,480));


void gMouseHandler( int button, int state, int x, int y );
void gOnMotion(int x, int y);
void gkeyboardHandler( unsigned char key, int x, int y );
void gDisplay();
void gReshape(int w, int h);
void gIdle();


// definition the handlers

void gMouseHandler( int button, int state, int x, int y )
{
    gGLWin.GLMouseHandler(button, state, x,y);
}

void gOnMotion(int x, int y)
{
    gGLWin.GLMouseMotion(x,y);
}

void gkeyboardHandler( unsigned char key, int x, int y )
{
    gGLWin.GLkeyboardHandler( key, x, y );
}

void gDisplay( )
{
    gGLWin.GLDisplay();
    cout << "gDisplay" << endl;
}

void gReshape(int w, int h)
{
    gGLWin.GLReshape(w,h);
}



/////////////////////////////////


GLWindow2::GLWindow2(Eigen::Vector2f irSize):
    mvMCPoseUpdate(6),
    mvLeftPoseUpdate(6)
{

    stateLeft   =   0;
    stateRight  =   0;
    stateMiddle =   0;
    mirVideoSize= irSize;

    //  GUI.RegisterCommand("GLWindow.AddMenu", GUICommandCallBack, this);
    //  glSetFont("sans");
    mvMCPoseUpdate  = Eigen::VectorXf::Zero(6);
    mvLeftPoseUpdate= Eigen::VectorXf::Zero(6);

    //  // initialize GLUT
    int argc = 1;
    char *argv[1] = {(char*)"Something"};
    glutInit(&argc, argv);
    glutInitDisplayMode( GLUT_RGBA | GLUT_DOUBLE );
    glutInitWindowSize( irSize[0], irSize[1] );
    glutCreateWindow( "PTAM" );

    //  // set up GUI callback functions
    //    glutDisplayFunc( gDisplay );
    glutReshapeFunc( gReshape );
    glutMouseFunc( gMouseHandler );
    glutKeyboardFunc( gkeyboardHandler );
    glutMotionFunc(gOnMotion);
    mbUserPressedSpacebar = false;

}


void GLWindow2::GLMouseHandler( int button, int state, int x, int y )
{



    if((button & GLUT_RIGHT_BUTTON) == GLUT_RIGHT_BUTTON)
        stateRight = !state;
    else  if((button & GLUT_MIDDLE_BUTTON) == GLUT_MIDDLE_BUTTON)
        stateMiddle = !state;
    else if((button & GLUT_LEFT_BUTTON) ==  GLUT_LEFT_BUTTON)
        stateLeft = !state;


}
void GLWindow2::GLMouseMotion(  int x, int y )
{
    Eigen::Vector2f where;
    where[0] = x;
    where[1] = y;
    Eigen::Vector2f irMotion = where - mirLastMousePos;

    mirLastMousePos = where;

    double dSensitivity = 0.01;

    if((stateLeft) )
    {

        mvMCPoseUpdate[3] -= irMotion[1] * dSensitivity;
        mvMCPoseUpdate[4] += irMotion[0] * dSensitivity;
    }
    else if((stateRight))
    {

        mvLeftPoseUpdate[4] -= irMotion[0] * dSensitivity;
        mvLeftPoseUpdate[3] += irMotion[1] * dSensitivity;
    }
    else if(stateMiddle )
    {

        mvLeftPoseUpdate[5] -= irMotion[0] * dSensitivity;
        mvLeftPoseUpdate[2] += irMotion[1] * dSensitivity;
    }
}




void GLWindow2::GLkeyboardHandler( unsigned char key, int x, int y )
{

    switch ( key )
    {
    case 'q':
        // quit when q is pressed
        exit(0);
        break;

    case 32:
        mbUserPressedSpacebar = !mbUserPressedSpacebar;
        break;

    case 'o':
        mvMCPoseUpdate[2] += 0.2;
        break;

    case 'l':
        mvMCPoseUpdate[2] -= 0.2;
        break;

    default:
        break;
    }
}

void GLWindow2::GLDisplay()
{

}
void GLWindow2::GLReshape(int w,int h)
{

}


void GLWindow2::SetupUnitOrtho()
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0,1,1,0,0,1);
}

void GLWindow2::SetupWindowOrtho()
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    //  glOrtho(0,size());
    glOrtho(-(double)mirVideoSize[0]/2,(double)mirVideoSize[0]/2, -(double) mirVideoSize[1]/2 , (double) mirVideoSize[1]/2, -1.0, 1.0);

}

void GLWindow2::SetupVideoOrtho()
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    //  glOrtho(-0.5,(double)mirVideoSize.x - 0.5, (double) mirVideoSize.y - 0.5, -0.5, -1.0, 1.0);
    //    glOrtho(-(double)mirVideoSize.x/2,(double)mirVideoSize.x/2, -(double) mirVideoSize.y/2 , (double) mirVideoSize.y/2, -1.0, 1.0);
    glOrtho(0,(double)mirVideoSize[0], (double) mirVideoSize[1] , 0, -1.0, 1.0);

}

void GLWindow2::SetupVideoRasterPosAndZoom()
{
    //  glRasterPos2d(-0.5,-0.5);
    //  double adZoom[2];
    //  adZoom[0] = (double) size()[0] / (double) mirVideoSize[0];
    //  adZoom[1] = (double) size()[1] / (double) mirVideoSize[1];
    //  glPixelZoom(adZoom[0], -adZoom[1]);
}

void GLWindow2::SetupViewport()
{
    glViewport(0, 0, mirVideoSize[0], mirVideoSize[1]);
}



std::pair<Eigen::VectorXf, Eigen::VectorXf> GLWindow2::GetMousePoseUpdate()
{
    pair<Eigen::VectorXf, Eigen::VectorXf > result = make_pair(mvLeftPoseUpdate, mvMCPoseUpdate);
//    mvLeftPoseUpdate = Eigen::VectorXf::Zero(6);
    mvMCPoseUpdate = Eigen::VectorXf::Zero(6);
    return result;
}

