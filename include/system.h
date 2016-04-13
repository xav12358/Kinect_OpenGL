// -*- c++ -*-
// Copyright 2008 Isis Innovation Limited
//
// System.h
//
// Defines the System class
//
// This stores the main functional classes of the system, like the
// mapmaker, map, tracker etc, and spawns the working threads.
//
#ifndef __SYSTEM_H
#define __SYSTEM_H

#include <opencv/cv.h>
#include <opencv/highgui.h>

#include <include/glWindow2.h>
#include <include/kinect.h>
using namespace cv;

class MapViewer;

class System
{
public:
  System(GLWindow2 &GL);
  void Run();
  void GetAndFillFrameBWandRGB();
  
private:
//  cv::Mat mimFrameRGB;
//  cv::Mat mimFrameBW;
  kinect *SourceKinect;
  
  GLWindow2 &GLwin;
  MapViewer *mpMapViewer;

  bool mbDone;
};



#endif
