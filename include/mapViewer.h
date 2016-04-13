// -*- c++ -*-
// Copyright 2008 Isis Innovation Limited
//
// MapViewer.h
//
// Defines the MapViewer class
//
// This defines a simple map viewer widget, which can draw the 
// current map and the camera/keyframe poses within it.
//
#ifndef __MAP_VIEWER_H
#define __MAP_VIEWER_H

//#include <TooN/TooN.h>
//using namespace TooN;
//#include <TooN/se3.h>
//#include <TooN/Lapack_Cholesky.h>


#include <sstream>

#include <include/glWindow2.h>
#include <include/kinect.h>
#include <include/Eigen/Dense>


class MapViewer
{
public:
  MapViewer(GLWindow2 &glw);
  void DrawMap(Eigen::Matrix4f se3CamFromWorld, kinect *SourceK);
  std::string GetMessageForUser();
  
protected:
  GLWindow2 &mGLWindow;
  std::vector<std::pair<GLuint,int > *>textures;
  
  void DrawGrid();
  void DrawCamera(Eigen::Matrix4f se3, bool bSmall=false);
  void DrawCone(int NbPoints,float size);
  void DrawCloud(Eigen::Matrix4f se3CfromW, kinect *SourceK);

  void SetupFrustum();
  void SetupModelView(Eigen::Matrix4f se3WorldFromCurrent);

  
  Eigen::Vector3f mv3MassCenter;
  Eigen::Matrix4f mse3ViewerFromWorld;
  Eigen::Matrix4f CovPos;

};

#endif
