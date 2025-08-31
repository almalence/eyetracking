/* ------------------------------------------------------------------------- *\

    Eye tracking calibration: Points animation

    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
    SPDX-License-Identifier: GPL-3.0-only
    All rights reserved.

    Authors: Bogdan Nudnenko
             Dmitry Shmunk

\* ------------------------------------------------------------------------- */

#include "interaction_module.h"

// 5 points for calibration
std::vector<AnglePointInDegrees> calibrationPoints =
{
    {-15.0f,  0.0f},    // LEFT
    {15.0f,   0.0f},    // RIGHT
    {0.0f,  15.5f},     // UP
    {0.0f,  -15.5f},    // DOWN
    {0.0f,  0.0f},      // CENTER
};


#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))


NormalizedScreenPoint InteractionModule::anglesToNormalized(const AnglePointInDegrees& point) const
{
    NormalizedScreenPoint normScreenPoint{};

    normScreenPoint.x = (point.horizontal * (m_pi / 180.0f) - mergedFov.angleLeft) / (mergedFov.angleRight - mergedFov.angleLeft);
    normScreenPoint.y = (point.vertical   * (m_pi / 180.0f) - mergedFov.angleDown) / (mergedFov.angleUp - mergedFov.angleDown);

    return normScreenPoint;
}

void InteractionModule::startInteractionCycle()
{
    while (m_running) {
        auto startTime = std::chrono::high_resolution_clock::now();

        while (true) {
            auto currentTime = std::chrono::high_resolution_clock::now();
            float elapsedTime = std::chrono::duration<float>(currentTime - startTime).count();
            elapsedTime = fminf(POINT_SWITCH_DURATION, elapsedTime);

            {
                std::lock_guard<std::mutex> lock(m_mutex);
                m_currentRadius = m_maxRadiusValue - (m_maxRadiusValue - m_minRadiusValue) * (elapsedTime / POINT_SWITCH_DURATION);
            }

            if (elapsedTime >= POINT_SWITCH_DURATION)
            {
                {
                    if (!m_normalizedPoints.empty())
                    {
                        std::lock_guard<std::mutex> lock(m_mutex);

                        if (m_currentIndex < m_normalizedPoints.size())
                            ++m_currentIndex;
                        if (m_currentIndex == m_normalizedPoints.size())
                            //m_currentInteractionState = InteractionProgramState::ProgramStoped;
                            m_currentInteractionState = InteractionProgramState::CalibrationCaptured;
                    }
                }
                break;
            }
            std::this_thread::sleep_for(std::chrono::milliseconds(16));
        }

        if (m_currentInteractionState == InteractionProgramState::ProgramStoped) {
            break;
        }
    }
}

void InteractionModule::cleanResources()
{
    m_running = false;
    if (m_timerThread.joinable()) {
        m_timerThread.join();
    }
}

void InteractionModule::startInteractionProgram(std::vector<XrFovf> fov)
{
    mergedFov.angleLeft = MIN(fov[0].angleLeft, fov[1].angleLeft);
    mergedFov.angleRight = MAX(fov[0].angleRight, fov[1].angleRight);
    mergedFov.angleUp = MAX(fov[0].angleUp, fov[1].angleUp);
    mergedFov.angleDown = MIN(fov[0].angleDown, fov[1].angleDown);

    if (mergedFov.angleLeft == 0 || mergedFov.angleRight == 0 ||
        mergedFov.angleUp == 0 || mergedFov.angleDown == 0) {
        std::cout << "Error: Wrong input FOV data. One or more values are 0" << std::endl;
        return;
    }

    for (auto point : calibrationPoints)
        m_normalizedPoints.push_back(anglesToNormalized(point));

    // Start the program on the first request
    if (!m_running) {
        m_running = true;
        m_timerThread = std::thread(&InteractionModule::startInteractionCycle, this);
    }

    m_currentInteractionState = InteractionProgramState::ProgramStarted;
}

InteractionProgramState InteractionModule::getInteractionState() const {
    return m_currentInteractionState;
}

void InteractionModule::setGazeSource(std::shared_ptr<CaptureModule> obj, GazeFunction gazeFunc)
{
    m_gazeModuleObj = obj;
    m_gazeFunction = gazeFunc;
}

void InteractionModule::interactionData(float* coordinates, float* radius, int* state)
{
    if ((m_currentInteractionState == InteractionProgramState::ProgramStarted) ||
        (m_currentInteractionState == InteractionProgramState::CalibrationCaptured))
    {

        std::lock_guard<std::mutex> lock(m_mutex);

        if (m_currentInteractionState == InteractionProgramState::CalibrationCaptured)
            *state = 1;
        else
            *state = 0;

        if ((!m_normalizedPoints.empty()) && (m_currentIndex < m_normalizedPoints.size()))
        {
            coordinates[0] = m_normalizedPoints[m_currentIndex].x;
            coordinates[1] = m_normalizedPoints[m_currentIndex].y;
            *radius = m_currentRadius;
        }
        else // after calibration completed - use estimated gaze coordinate
        {
            float gaze_vector[3] = { 0,0,1 };
            if (m_gazeFunction)
            {
                (m_gazeModuleObj.get()->*m_gazeFunction)(gaze_vector);
            }

            AnglePointInDegrees gaze_point;
            gaze_point.horizontal = 180 / (float)m_pi * asinf(-gaze_vector[0]);
            gaze_point.vertical = 180 / (float)m_pi * asinf(gaze_vector[1]);

            auto screenPoint = anglesToNormalized(gaze_point);

            coordinates[0] = screenPoint.x;
            coordinates[1] = screenPoint.y;
            *radius = m_maxRadiusValue;
        }
    }
}

void InteractionModule::gazeReferenceData(float* coordinates, int *idx)
{
    std::lock_guard<std::mutex> lock(m_mutex);

    if (m_currentIndex < calibrationPoints.size())
    {
        coordinates[0] = calibrationPoints[m_currentIndex].horizontal;
        coordinates[1] = calibrationPoints[m_currentIndex].vertical;
        *idx = m_currentIndex;
    }
}

InteractionModule::~InteractionModule()
{
    cleanResources();
}
