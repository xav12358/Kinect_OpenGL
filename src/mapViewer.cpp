#include "include/mapViewer.h"
#include <iomanip>


#include <GL/gl.h>
#include <GL/glut.h>
#include <GL/freeglut.h>

#include <math.h>
#include <iostream>
#include <string>

using namespace std;

#include <opencv/cv.h>
#include <opencv/highgui.h>


GLfloat ColorMW[4][4] = {{1.0f, 0.0f, 0.0f},
                         {1.0f, 1.0f, 0.0f},
                         {0.0f, 1.0f, 0.0f},
                         {0.0f, 0.0f, 0.7f} };

MapViewer::MapViewer(GLWindow2 &glw):
    mGLWindow(glw)
{
//    mse3ViewerFromWorld =
//            SE3<>::exp(makeVector(0,0,2,0,0,0)) * SE3<>::exp(makeVector(0,0,0,0.8 * M_PI,0,0));
}


void MapViewer::DrawGrid()
{
    SetupFrustum();
    SetupModelView(Eigen::Matrix4f::Identity());
    glLineWidth(1);

    glBegin(GL_LINES);

    // Draw a larger grid around the outside..
    double dGridInterval = 0.1;

    double dMin = -100.0 * dGridInterval;
    double dMax =  100.0 * dGridInterval;

    for(int x=-10;x<=10;x+=1)
    {
        if(x==0)
            glColor3f(1,1,1);
        else
            glColor3f(0.3,0.3,0.3);
        glVertex3d((double)x * 10 * dGridInterval, dMin, 0.0);
        glVertex3d((double)x * 10 * dGridInterval, dMax, 0.0);
    }
    for(int y=-10;y<=10;y+=1)
    {
        if(y==0)
            glColor3f(1,1,1);
        else
            glColor3f(0.3,0.3,0.3);
        glVertex3d(dMin, (double)y * 10 *  dGridInterval, 0.0);
        glVertex3d(dMax, (double)y * 10 * dGridInterval, 0.0);
    }

    glEnd();

    glBegin(GL_LINES);
    dMin = -10.0 * dGridInterval;
    dMax =  10.0 * dGridInterval;

    for(int x=-10;x<=10;x++)
    {
        if(x==0)
            glColor3f(1,1,1);
        else
            glColor3f(0.5,0.5,0.5);

        glVertex3d((double)x * dGridInterval, dMin, 0.0);
        glVertex3d((double)x * dGridInterval, dMax, 0.0);
    }
    for(int y=-10;y<=10;y++)
    {
        if(y==0)
            glColor3f(1,1,1);
        else
            glColor3f(0.5,0.5,0.5);
        glVertex3d(dMin, (double)y * dGridInterval, 0.0);
        glVertex3d(dMax, (double)y * dGridInterval, 0.0);
    }

    glColor3f(1,0,0);
    glVertex3d(0,0,0);
    glVertex3d(1,0,0);
    glColor3f(0,1,0);
    glVertex3d(0,0,0);
    glVertex3d(0,1,0);
    glColor3f(1,1,1);
    glVertex3d(0,0,0);
    glVertex3d(0,0,1);
    glEnd();

}

void MapViewer::DrawMap(Eigen::Matrix4f se3CamFromWorld, kinect *SourceK)
{

    // Update viewer position according to mouse input:
    pair<Eigen::VectorXf, Eigen::VectorXf > pv6 = mGLWindow.GetMousePoseUpdate();
    Eigen::Matrix3f Matrix4f;
//    se3CamFromMC.get_translation() = mse3ViewerFromWorld * mv3MassCenter;
//    mse3ViewerFromWorld = SE3<>::exp(pv6.first) *
//            se3CamFromMC * SE3<>::exp(pv6.second) * se3CamFromMC.inverse() * mse3ViewerFromWorld;

    //     affichage matrice
    //    for(int i=0;i<3;i++)
    //    {
    //        for(int j=0;j<3;j++)
    //            cout <<se3CamFromWorld.get_rotation().get_matrix()[i][j] << " ";

    //        cout << endl;
    //    }

    //    for(int j=0;j<3;j++)
    //        cout << se3CamFromWorld.get_translation()[j] << " ";
    //    cout << endl;

    glDisable(GL_DEPTH_TEST);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();


    mGLWindow.SetupViewport();
    glClearColor(0,0,0,0);
    glClearDepth(1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_LINE_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glColorMask(1,1,1,1);
    glEnable(GL_DEPTH_TEST);


    DrawGrid();
    DrawCamera(se3CamFromWorld, true);
    cout << "DrawCloud" << endl;
    DrawCloud(se3CamFromWorld,SourceK);

    // show the rendering on the screen
    glutSwapBuffers();
    // post the next redisplay
    glutPostRedisplay();

}


void MapViewer::DrawCamera(Eigen::Matrix4f se3CfromW, bool bSmall)
{

    SetupModelView(se3CfromW.inverse());
    SetupFrustum();

    if(bSmall)
        glLineWidth(1);
    else
        glLineWidth(3);

    glBegin(GL_LINES);
    glColor3f(1,0,0);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(0.1f, 0.0f, 0.0f);
    glColor3f(0,1,0);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 0.1f, 0.0f);
    glColor3f(1,1,1);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 0.0f, 0.1f);
    glEnd();


    //    TooN::SE3<> Se3Tranformtmp = SE3<>::exp(makeVector(0,0,-0.1,-M_PI/2,0,0));
    //    //     affichage matrice
    //    cout << "Transform1" << endl;
    //    for(int i=0;i<3;i++)
    //    {
    //        for(int j=0;j<3;j++)
    //            cout <<Se3Tranformtmp.get_rotation().get_matrix()[i][j] << " ";

    //        cout << endl;
    //    }

    //    for(int j=0;j<3;j++)
    //        cout << Se3Tranformtmp.get_translation()[j] << " ";
    //    cout << endl;


    //    cout << "Transform2" << endl;

    //    Se3Tranformtmp = SE3<>::exp(makeVector(0,0,-0.1,0,0,0))*SE3<>::exp(makeVector(0,0,0,-M_PI/2,0,0));

    //    for(int i=0;i<3;i++)
    //    {
    //        for(int j=0;j<3;j++)
    //            cout <<Se3Tranformtmp.get_rotation().get_matrix()[i][j] << " ";

    //        cout << endl;
    //    }

    //    for(int j=0;j<3;j++)
    //        cout << Se3Tranformtmp.get_translation()[j] << " ";
    //    cout << endl;


    //    Se3Tranformtmp = SE3<>::exp(makeVector(0,0,0,-M_PI/2,0,0))*SE3<>::exp(makeVector(0,0,-0.1,0,0,0));
    //    //     affichage matrice
    //    cout << "Transform3" << endl;
    //    for(int i=0;i<3;i++)
    //    {
    //        for(int j=0;j<3;j++)
    //            cout <<Se3Tranformtmp.get_rotation().get_matrix()[i][j] << " ";

    //        cout << endl;
    //    }

    //    for(int j=0;j<3;j++)
    //        cout << Se3Tranformtmp.get_translation()[j] << " ";
    //    cout << endl;


    /////////////////////////////////
//    glColor3f(1,1,1);
//    TooN::SE3<> Se3Tranform1 = SE3<>::exp(makeVector(0,0,-0.1,0,0,0));
//    SetupModelView((Se3Tranform1*se3CfromW).inverse());
//    SetupFrustum();
//    DrawCone(20,0.01);

//    glColor3f(1,0,0);
//    Se3Tranform1 = SE3<>::exp(makeVector(0,0,0,0,-M_PI/2,0))*SE3<>::exp(makeVector(-0.1,0,0,0,0,0));
//    SetupModelView((Se3Tranform1*se3CfromW).inverse());
//    SetupFrustum();
//    DrawCone(20,0.01);

//    glColor3f(0,1,0);
//    Se3Tranform1 = SE3<>::exp(makeVector(0,0,0,M_PI/2,0,0))*SE3<>::exp(makeVector(0,-0.1,0,0,0,0));
//    SetupModelView((Se3Tranform1*se3CfromW).inverse());
//    SetupFrustum();
//    DrawCone(20,0.01);


    if(!bSmall)
    {
        double dcoef = 0.1;
        glEnable(GL_ALPHA_TEST);
        glAlphaFunc(GL_GREATER, 0.0f);

        glBegin(GL_QUADS);

        glColor4f(1,1,0.6,0.25);
        // BACK
        glVertex3f(-0.5*dcoef, -0.5*dcoef, 0.0f);
        glVertex3f(0.5*dcoef, -0.5*dcoef, 0.0f);
        glVertex3f(0.5*dcoef, 0.5*dcoef, 0.0f);
        glVertex3f(-0.5*dcoef, 0.5*dcoef, 0.0f);


        //FRONTGL_ALPHA_TEST
        glVertex3f(-0.5*dcoef, -0.5*dcoef, 0.5*dcoef);
        glVertex3f(0.5*dcoef, -0.5*dcoef, 0.5*dcoef);
        glVertex3f(0.5*dcoef, 0.5*dcoef, 0.5*dcoef);
        glVertex3f(-0.5*dcoef, 0.5*dcoef, 0.5*dcoef);

        // LEFT
        glVertex3f(-0.5*dcoef, -0.5*dcoef, 0.0f);
        glVertex3f(-0.5*dcoef, 0.5*dcoef, 0.0f);
        glVertex3f(-0.5*dcoef, 0.5*dcoef, 0.5*dcoef);
        glVertex3f(-0.5*dcoef, -0.5*dcoef, 0.5*dcoef);

        // RIGHT
        glVertex3f(0.5*dcoef, -0.5*dcoef, 0.0f);
        glVertex3f(0.5*dcoef, 0.5*dcoef, 0.0f);
        glVertex3f(0.5*dcoef, 0.5*dcoef, 0.5*dcoef);
        glVertex3f(0.5*dcoef, -0.5*dcoef, 0.5*dcoef);

        // TOP
        glVertex3f(0.5*dcoef, -0.5*dcoef, 0.0f);
        glVertex3f(0.5*dcoef, -0.5*dcoef, 0.5*dcoef);
        glVertex3f(-0.5*dcoef, -0.5*dcoef, 0.5*dcoef);
        glVertex3f(-0.5*dcoef, -0.5*dcoef, 0.0f);

        // BUTTOM
        glVertex3f(0.5*dcoef, 0.5*dcoef, 0.0f);
        glVertex3f(0.5*dcoef, 0.5*dcoef, 0.5*dcoef);
        glVertex3f(-0.5*dcoef, 0.5*dcoef, 0.5*dcoef);
        glVertex3f(-0.5*dcoef, 0.5*dcoef, 0.0f);

        glEnd();
        glDisable(GL_ALPHA_TEST);

        glLineWidth(1);
        glColor3f(0.5,0.5,0.5);
//        SetupModelView();


//        Eigen::Vector2f v2CamPosXY = se3CfromW.inverse().get_translation().slice<0,2>();

//        glBegin(GL_LINES);
//        glColor3f(1,0,1);
//        glVertex2d(v2CamPosXY[0] - 0.04, v2CamPosXY[1] + 0.04);
//        glVertex2d(v2CamPosXY[0] + 0.04, v2CamPosXY[1] - 0.04);
//        glVertex2d(v2CamPosXY[0] - 0.04, v2CamPosXY[1] - 0.04);
//        glVertex2d(v2CamPosXY[0] + 0.04, v2CamPosXY[1] + 0.04);
//        glEnd();
    }

}


void MapViewer::DrawCone(int NbPoints, float size)
{
    //    glColor3f(0.5,0.5,0.5);

    glBegin(GL_TRIANGLES);

    float fstep = 2*M_PI /NbPoints;
    for(int i=0;i<NbPoints;i++)
    {
        glVertex3d(0, 0, 0);
        glVertex3d(size*cos(i*fstep), size*sin(i*fstep), 0);
        glVertex3d(size*cos((i+1)*fstep), size*sin((i+1)*fstep), 0);

        glVertex3d(0, 0, 2.5*size);
        glVertex3d(size*cos(i*fstep), size*sin(i*fstep), 0);
        glVertex3d(size*cos((i+1)*fstep), size*sin((i+1)*fstep), 0);

    }

    glEnd();

}


void MapViewer::DrawCloud(Eigen::Matrix4f se3CfromW,kinect *SourceK)
{
//    SetupModelView(se3CfromW.inverse());
//    SetupFrustum();

//    glPointSize(0.1);
//    glBegin(GL_POINTS);
//    glColor3f(1,1,0);

//    for(int i=0;i<SourceK->mvSize[0]*SourceK->mvSize[1];i++)
//    {
//        TooN::Vector<3> V3DPoint = SourceK->VTab3DImage[i];
//        if(V3DPoint[2]!=0)
//            glVertex3f(V3DPoint[0],V3DPoint[1],V3DPoint[2]);

//    }
//    glEnd();

}



void MapViewer::SetupFrustum()
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    double zNear = 0.03;
    glFrustum(-zNear, zNear, 0.75*zNear,-0.75*zNear,zNear,50);
    glScalef(1,1,-1);
    return;
}

void MapViewer::SetupModelView(Eigen::Matrix4f se3WorldFromCurrent)
{
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    Eigen::Matrix4f tmp = mse3ViewerFromWorld * se3WorldFromCurrent;

    GLfloat matrix[4][4];
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            matrix[j][i] = tmp(i,j);

    for(int i=0;i<3;i++)
        matrix[3][i] = tmp(i,3);

    for(int i=0;i<4;i++)
        matrix[i][3] = (i==3);

    glMultMatrixf(matrix[0]);

    return;
}
