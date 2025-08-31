/* ------------------------------------------------------------------------- *\

    Eye tracking calibration: apps main()

    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
    SPDX-License-Identifier: GPL-3.0-only
    All rights reserved.

    Authors: Bogdan Nudnenko
             Dmitry Shmunk

\* ------------------------------------------------------------------------- */

#include "common.h"
#include "openxr_module.h"
#include "graphics_module.h"
#include "interaction_module.h"
#include "capture_module.h"

#ifdef XR_USE_PLATFORM_ANDROID

#include <android_native_app_glue.h>

struct AndroidAppState {
    ANativeWindow* NativeWindow = nullptr;
    bool Resumed = false;
};

struct PlatformData {
    void* applicationVM;
    void* applicationActivity;
};

static void app_handle_cmd(struct android_app* app, int32_t cmd) {

    AndroidAppState* appState = (AndroidAppState*)app->userData;

    switch (cmd) {
        // There is no APP_CMD_CREATE. The ANativeActivity creates the
        // application thread from onCreate(). The application thread
        // then calls android_main().
        case APP_CMD_START: {
            break;
        }
        case APP_CMD_RESUME: {
            appState->Resumed = true;
            break;
        }
        case APP_CMD_PAUSE: {
            appState->Resumed = false;
            break;
        }
        case APP_CMD_STOP: {
            break;
        }
        case APP_CMD_DESTROY: {
            appState->NativeWindow = NULL;
            break;
        }
        case APP_CMD_INIT_WINDOW: {
            appState->NativeWindow = app->window;
            break;
        }
        case APP_CMD_TERM_WINDOW: {
            appState->NativeWindow = NULL;
            break;
        }
    }
}

/**
 * This is the main entry point of a native application that is using
 * android_native_app_glue.  It runs in its own thread, with its own
 * event loop for receiving input events and doing other things.
 */
void android_main(struct android_app* app) {

    auto openXrModule = std::make_shared<OpenXRModule>();
    auto openGlModule = std::make_shared<GraphicsModule>();
    auto interactionModule = std::make_shared <InteractionModule>();
    auto captureModule = std::make_shared <CaptureModule>();

    AndroidAppState appState = {};
    app->userData = &appState;
    app->onAppCmd = app_handle_cmd;

    std::shared_ptr<PlatformData> data = std::make_shared<PlatformData>();
    data->applicationVM = app->activity->vm;
    data->applicationActivity = app->activity->clazz;


    PFN_xrInitializeLoaderKHR initializeLoader = nullptr;
    if (XR_SUCCEEDED(xrGetInstanceProcAddr(XR_NULL_HANDLE, "xrInitializeLoaderKHR", (PFN_xrVoidFunction*)(&initializeLoader)))) {
        XrLoaderInitInfoAndroidKHR loaderInitInfoAndroid = {XR_TYPE_LOADER_INIT_INFO_ANDROID_KHR};
        loaderInitInfoAndroid.applicationVM = app->activity->vm;
        loaderInitInfoAndroid.applicationContext = app->activity->clazz;
        initializeLoader((const XrLoaderInitInfoBaseHeaderKHR*)&loaderInitInfoAndroid);
    }


    openXrModule->createInstance(data->applicationVM, data->applicationActivity);
    openXrModule->initializeSystem();
    openXrModule->initOpenGlExtension();

    openGlModule->initializeDevice();
    openGlModule->initializeResources();

    openXrModule->completeGraphicsBinding(openGlModule->getWindow());

    openXrModule->createSession();
    openXrModule->configureReferenceSpace();

    openXrModule->createSwapchains(openGlModule->SupportedColorSwapchainFormats);

    openXrModule->setDrawFunction(openGlModule, &GraphicsModule::draw);
    openGlModule->setInteractionProgram(interactionModule, &InteractionModule::interactionData);

    captureModule->startCapture();

    bool exitRenderLoop = false;
    while (app->destroyRequested == 0) {

        for (;;) {
            int events;
            struct android_poll_source* source;
            // If the timeout is zero, returns immediately without blocking.
            // If the timeout is negative, waits indefinitely until an event appears.
            const int timeoutMilliseconds =
                    (!appState.Resumed && (interactionModule->getInteractionState() != InteractionProgramState::ProgramStarted) && app->destroyRequested == 0) ? -1 : 0;
            if (ALooper_pollOnce(timeoutMilliseconds, nullptr, &events, (void**)&source) < 0) {
                break;
            }

            // Process this event.
            if (source != nullptr) {
                source->process(app, source);
            }
        }

        openXrModule->pollEvents(&exitRenderLoop);
        if (exitRenderLoop) {
            break;
        }

        if (interactionModule->getInteractionState() == InteractionProgramState::ProgramStoped)
        {
            break;
        }

        openXrModule->renderFrame();

        if (interactionModule->getInteractionState() == InteractionProgramState::ProgramNotInitialized)
        {
            interactionModule->startInteractionProgram(openXrModule->getFov()); // FOV is available only after the first frame.
        }
    }

    app->activity->vm->DetachCurrentThread();

    captureModule->stopCapture();   // cleanResources called from within

    openXrModule->cleanResources();
    interactionModule->cleanResources();

    ANativeActivity_finish(app->activity);
}
#else
int main(int argc, char* argv[]) {

    auto openXrModule = std::make_shared<OpenXRModule>();
    auto openGlModule = std::make_shared<GraphicsModule>();
    auto interactionModule = std::make_shared <InteractionModule>();
    auto captureModule = std::make_shared <CaptureModule>();


    openXrModule->createInstance();
    openXrModule->initializeSystem();
    openXrModule->initOpenGlExtension();

    openGlModule->initializeDevice();
    openGlModule->initializeResources();

    openXrModule->completeGraphicsBinding(openGlModule->getWindow());

    openXrModule->createSession();
    openXrModule->configureReferenceSpace();

    openXrModule->createSwapchains(openGlModule->SupportedColorSwapchainFormats);

    openXrModule->setDrawFunction(openGlModule, &GraphicsModule::draw);
    openGlModule->setInteractionProgram(interactionModule, &InteractionModule::interactionData);

    interactionModule->setGazeSource(captureModule, &CaptureModule::gazeData);

    captureModule->startCapture();
    captureModule->setMode(CaptureProgramState::CollectCalibration);
    captureModule->setGazeReferenceSource(interactionModule, &InteractionModule::gazeReferenceData);

    bool calibrated = false;
    std::cout << "Press ESC key to exit..." << std::endl;
    while (!(GetAsyncKeyState(VK_ESCAPE) & 0x8000))
    {
    // while (!quitKeyPressed) {
        bool exitRenderLoop = false;
        openXrModule->pollEvents(&exitRenderLoop);
        if (exitRenderLoop) break;

        if (interactionModule->getInteractionState() == InteractionProgramState::ProgramStoped) break;

        // once calibration frames captured - calibrate and switch to testing
        if ((interactionModule->getInteractionState() == InteractionProgramState::CalibrationCaptured) && (!calibrated))
        {
            captureModule->setMode(CaptureProgramState::NoCapture);

            // perform calibration with the frames captured
            float anc[5][2];
            for (int i=0; i<5; ++i)
            {
                anc[i][0] = calibrationPoints[i].horizontal;
                anc[i][1] = calibrationPoints[i].vertical;
            }
            int err;
            eye_calib* ec = calibrateEye(captureModule->ei, captureModule->frames, captureModule->calib_frames, anc, &err);

            if (err != EYECAL_ALL_OK)
            {
                char calibrationErrorText[EYECAL_MAX_ERR][256] = { "",
                    "Not enough memory",
                    "Not enough frames with visible pupil",
                    "Eye rotational axes are out of range"
                };
                printf("Calibration failed: %s\n", calibrationErrorText[err]);
            }
            else
            {
                // print-out calibration results
                char calib_str[1024];
                snprintf(calib_str, 1024, "eye_calib ec = { "
                    ".iris2hrot = %5.2f, "
                    ".iris2vrot = %5.2f, "
                    ".iris2pupil = %4.2f, "
                    ".pupil_decenter = {%5.2f,%5.2f}, "
                    ".eyeball_ctr = {%4.0f,%4.0f}, "
                    ".pca_angle = {%7.2f,%7.2f}, "
                    ".pupil_refr = %5.3f, "
                    ".pupil_angle = %8.2f, "
                    ".pupil_stretch = %5.3f, "
                    ".gaze_alpha = {%7.2f,%7.2f} }\n",
                    ec->iris2hrot, ec->iris2vrot, ec->iris2pupil,
                    ec->pupil_decenter[0], ec->pupil_decenter[1], ec->eyeball_ctr[0], ec->eyeball_ctr[1], ec->pca_angle[0], ec->pca_angle[1],
                    ec->pupil_refr, ec->pupil_angle, ec->pupil_stretch, ec->gaze_alpha[0], ec->gaze_alpha[1]);
                std::cout << calib_str;

                // dump calibration string into file
                std::ofstream cal_file("calibration.txt");
                cal_file << calib_str;
                cal_file.close();
            }

            calibrated = true;

            applyCalibration(captureModule->ei, ec);
            releaseCalibration(ec);
            captureModule->setMode(CaptureProgramState::EstimateEyePose);
        }

        openXrModule->renderFrame();

        if (interactionModule->getInteractionState() == InteractionProgramState::ProgramNotInitialized)
            interactionModule->startInteractionProgram(openXrModule->getFov()); // FOV is available only after the first frame.

        std::this_thread::sleep_for(std::chrono::milliseconds(16));
    }

    captureModule->stopCapture();   // cleanResources() called from within

    openXrModule->cleanResources();
    openGlModule->cleanResources();
    interactionModule->cleanResources();
}
#endif