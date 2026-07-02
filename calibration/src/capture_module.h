/* ------------------------------------------------------------------------- *\

    Eye tracking calibration: Capture eye images from eye-tracking camera
 
    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
    SPDX-License-Identifier: GPL-3.0-only
    All rights reserved.

    Author: Dmitry Shmunk

\* ------------------------------------------------------------------------- */

#pragma once

#include <opencv2/opencv.hpp>

#include "common.h"
#include "EyeTracking.h"


#define SOMNIUM_VR1

#ifdef SOMNIUM_VR1
#define CAM_FPS       50
#define CAM_WIDTH     640
#define CAM_HEIGHT    400
#define FRAME_WIDTH   1280
#define FRAME_HEIGHT  400
#define DEWARP_FRAME_WIDTH   540
#define DEWARP_FRAME_HEIGHT  400
#define CALIB_FPS     25
#endif

#define MAX_CAL_FRAMES  2000    // use at most that many frames for calibration (normally 25-30fps * 5 positions * 2 sec/pos = 300)

enum CaptureProgramState
{
    NoCapture,
    CollectCalibration,
    EstimateEyePose
};

class InteractionModule;

class CaptureModule {
public:
    std::mutex m_mutex;

    eye_cfg* ei = nullptr;
    float gaze_vector[3];

    uint8_t* frames[MAX_CAL_FRAMES] = { nullptr }; // array of pointers to frames captured for calibration
    int calib_frames;                              // currently captured frames for calibration

    typedef void (InteractionModule::* GazeReferenceFunction)(float*, int*);
    void setGazeReferenceSource(std::shared_ptr<InteractionModule> obj, GazeReferenceFunction grefFunc);

    int startCapture(void);
    void stopCapture(void);
    void setMode(int _state);
    void gazeData(float* gaze_vec);

private:
    cv::VideoCapture cap;

    bool captureRun;
    bool captureStopped;
    int state; // capture run stage

    uint8_t* dewarped_frame = nullptr;

    std::shared_ptr<InteractionModule> m_interacionModuleObj = nullptr;
    GazeReferenceFunction m_grefFunction = nullptr;

    void runCapture();
    void cleanResources();
};
