/* ------------------------------------------------------------------------- *\

    Eye tracking calibration: Points animation

    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
    SPDX-License-Identifier: GPL-3.0-only
    All rights reserved.

    Author: Bogdan Nudnenko

\* ------------------------------------------------------------------------- */

#pragma once

#include "common.h"
#include <xr_linear.h>

#define POINT_SWITCH_DURATION 1.6f

enum InteractionProgramState
{
    ProgramNotInitialized,
    ProgramStarted,
    ProgramStoped,
    CalibrationCaptured
};

struct AnglePointInDegrees {
    float horizontal = 0.0f;
    float vertical = 0.0f;
};

struct NormalizedScreenPoint {
    float x; 
    float y;
};

class CaptureModule;

class InteractionModule {
public:
    typedef void (CaptureModule::* GazeFunction)(float* gaze_vector);
    void setGazeSource(std::shared_ptr<CaptureModule> obj, GazeFunction gazeFunc);

    void startInteractionProgram(std::vector<XrFovf> fov);
    InteractionProgramState getInteractionState() const;
    void interactionData(float* coordinates, float* radius, int* state);
    void gazeReferenceData(float* coordinates, int *idx);
    void cleanResources();
    ~InteractionModule();

private:
    XrFovf mergedFov = {};

    std::vector<NormalizedScreenPoint> m_normalizedPoints {};

    InteractionProgramState m_currentInteractionState = InteractionProgramState::ProgramNotInitialized;

    std::shared_ptr<CaptureModule> m_gazeModuleObj = nullptr;
    GazeFunction m_gazeFunction = nullptr;

    int m_currentIndex = 0;
    float m_currentRadius = 0.0001f;
    bool m_running = false;

    const float m_minRadiusValue = 0.002f;
    const float m_maxRadiusValue = 0.01f;
    const float m_pi = 3.141592653589793f;

    float smooth_gaze[2] = { 0,0 };

    std::thread m_timerThread;
    std::mutex m_mutex;

    NormalizedScreenPoint anglesToNormalized(const AnglePointInDegrees& point) const;
    void startInteractionCycle();
};

extern std::vector<AnglePointInDegrees> calibrationPoints;
