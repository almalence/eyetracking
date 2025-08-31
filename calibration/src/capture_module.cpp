/* ------------------------------------------------------------------------- *\

    Eye tracking calibration: Capture eye images from eye-tracking camera
 
    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
    SPDX-License-Identifier: GPL-3.0-only
    All rights reserved.

    Author: Dmitry Shmunk

\* ------------------------------------------------------------------------- */

#include "capture_module.h"
#include "interaction_module.h"


int CaptureModule::startCapture()
{
    captureRun = false;
    captureStopped = false;

    // -------- Initialize Video Capture

    // A stupid way of finding ET camera:
    // For OpenCV under windows there is no way to choose camera to be open by UUID or name, only by index
    // ET camera happens to be at arbitrary index, depending on presence of other USB cameras in the system
    // Therefore: trying to open camera, checking if default frame parameters looks proper,
    // and if so - decide that this is the ET camera we are looking for
    for (int cidx=0; cidx<32; ++cidx)
    {
        cap.open(cidx, cv::CAP_MSMF);

        if (cap.isOpened())
        {
            int width = static_cast<int>( cap.get(cv::CAP_PROP_FRAME_WIDTH) );
            int height = static_cast<int>( cap.get(cv::CAP_PROP_FRAME_HEIGHT) );
            int fps = static_cast<int>( cap.get(cv::CAP_PROP_FPS) );

            if ((width==640) && (height==400) && (fps==50)) break; // ET camera found
        }
        cap.release();
    }    

    cap.set(cv::CAP_PROP_CONVERT_RGB, false);
    cap.set(cv::CAP_PROP_FRAME_WIDTH, CAM_WIDTH);

    // Check if the camera opened successfully
    if (!cap.isOpened()) {
        std::cerr << "OpenCV Error: Could not open camera." << std::endl;
        return -1;
    }
    
    // To display the resulting frame (will not work as GL is used for VR)
    // cv::namedWindow("Video Frame", cv::WINDOW_AUTOSIZE);

    // -------- Initialize Eye Tracking algorithm

    dewarped_frame = (uint8_t*)malloc(DEWARP_FRAME_WIDTH * DEWARP_FRAME_HEIGHT);
    if (dewarped_frame == NULL) {
        std::cerr << "Frame allocation failed." << std::endl;
        return -1;
    }

    // initializing with values reasonable for off-axis camera, 640x480 camera frame
    int valid_area[6] = { 280, 200,  230, 200,  280,  128 };
    float cam2hmd[6] = { 0, 22, 0,  25, 0, 48 };
    // Negative sign in camera HFOV indicates chief ray converging from camera entrance pupil towards eye
    // This may happen in cases where camera is behind HMD lens and image is formed through some sort of concave mirror surface
    ei = initEyeConfig(DEWARP_FRAME_WIDTH, DEWARP_FRAME_HEIGHT, 65, 120, valid_area, -1, -30, 0, CAM_FPS, cam2hmd);

    if (ei == NULL) {
        std::cerr << "Initialization failed." << std::endl;
        return -1;
    }

    // start capture with gathering frames for calibration
    captureRun = true;
    state = CaptureProgramState::NoCapture;

    auto captureThread = std::thread{ [this] {
                runCapture();
                cleanResources();
                captureStopped = true;
            } };
    captureThread.detach();

    return 0;
}


void CaptureModule::setMode(int _state)
{
    state = _state;
}

void CaptureModule::gazeData(float* gaze_vec)
{
    std::lock_guard<std::mutex> lock(m_mutex);
    memcpy(gaze_vec, gaze_vector, 3*sizeof(float));
}

// supplementary function to rectify through-the-lens camera view of Somnium VR1
void somniumCamDewarp(uint8_t* in, uint8_t* out, int eye)
{
    // center of the output frame in input frame
    const int xc = 360;
    const int yc = 200;

    // width / height of frame from the camera (both eyes combined horizontally)
    const int w = 1280;
    const int h = 400;

    // width / height of the dewarped frame
    const int wd = 540;
    const int hd = 400;

    for (int yo = 0; yo < hd; ++yo)
    {
        int y = yo - hd / 2 + 1;

        for (int xo = 0; xo < wd; ++xo)
        {
            int x = xo - wd / 2 + 1;

            float r = sqrtf((float)x * x + (float)y * y);

            // pincussion and perspective (due to off-axis camera placement behind the lens) corrections
            float xd = x - 1 + xc - 0.16f * x + 0.0012f * x * x + 0.0007f * y * y + 0.00055f * x * r;
            float yd = y - 1 + yc - 0.10f * y + 0.0007f * x * y + 0.00055f * y * r;

            if (eye == 0) xd = w - 1 - xd;

            // nearest-neighbor interpolation
            // int xi = (int)(xd+0.5f);
            // int yi = (int)(yd+0.5f);
            // yi = yi<0 ? 0:yi;  yi = yi>h-1 ? h-1:yi;
            // xi = xi<0 ? 0:xi;  xi = xi>w-1 ? w-1:xi;
            // uint8_t px = in[xi+yi*w];

            // bilinear interpolation
            int xi = (int)xd;
            int yi = (int)yd;
            int xf = (int)((xd - xi) * 256);
            int yf = (int)((yd - yi) * 256);
            yi = yi < 0 ? 0 : yi;  yi = yi > h - 2 ? h - 2 : yi;
            xi = xi < 0 ? 0 : xi;  xi = xi > w - 2 ? w - 2 : xi;
            uint8_t px = ((in[xi + yi * w] * (256 - xf) + in[xi + 1 + yi * w] * xf) * (256 - yf) + (in[xi + (yi + 1) * w] * (256 - xf) + in[xi + 1 + (yi + 1) * w] * xf) * yf) / 65536;

            if (eye == 0) out[wd - 1 - xo + yo * wd] = px;
            else out[xo + yo * wd] = px;
        }

        // correcting light fall-off on the side farther away from camera
        for (int xo = 240; xo < wd; ++xo)
        {
            int x = xo;
            if (eye == 0) x = wd - 1 - xo;
            uint32_t px = out[x + yo * wd] * 5 * (205 + (xo - 240)) / 1024;
            out[x + yo * wd] = (uint8_t)(px > 255 ? 255 : px);
        }
    }
}

void CaptureModule::setGazeReferenceSource(std::shared_ptr<InteractionModule> obj, GazeReferenceFunction grefFunc)
{
    m_interacionModuleObj = obj;
    m_grefFunction = grefFunc;
}

// Loop to continuously capture frames
void CaptureModule::runCapture()
{
    int frame_idx = 0;
    calib_frames = 0;

    while (captureRun)
    {
        // Capture frame-by-frame
        cv::Mat rawframe;
        cap >> rawframe;

        // If the frame is captured correctly
        if (!rawframe.empty())
        {
            cv::Mat frame = rawframe.reshape(0, { FRAME_HEIGHT, FRAME_WIDTH });

#ifdef SOMNIUM_VR1
            // In Somnium VR1:
            // - right eye is on the left, left eye is on the right
            // - left eye image is rotated 180 degree, correcting for that
            cv::Mat lefteye;
            cv::rotate(frame(cv::Range(0, FRAME_HEIGHT), cv::Range(FRAME_WIDTH / 2, FRAME_WIDTH)), lefteye, cv::ROTATE_180);
            lefteye.copyTo(frame(cv::Range(0, FRAME_HEIGHT), cv::Range(FRAME_WIDTH / 2, FRAME_WIDTH)));
#endif

            // Display the resulting frame (will not work as GL is used for VR)
            //cv::imshow("Video Frame", frame);

            if (state == CaptureProgramState::CollectCalibration)
            {
                // dewarp frame and store for calibration, use 25-30fps for calibration
                if ((frame_idx % (CAM_FPS / 25) == 0) && (calib_frames < MAX_CAL_FRAMES))
                {
                    frames[calib_frames] = (uint8_t*)malloc(DEWARP_FRAME_WIDTH * DEWARP_FRAME_HEIGHT);
                    if (frames[calib_frames])
                    {
                        // dump frame
                        //std::stringstream frame_fname;
                        //frame_fname << std::setfill('0')<< std::setw(5) << calib_frames;
                        //frame_fname << ".bmp";
                        //cv::imwrite(frame_fname.str(), frame);
                        //// update dataset csv
                        //char csvName[] = "01.csv";
                        //std::ifstream csvFile(csvName);
                        //std::ofstream csv;
                        //csv.open(csvName, std::ios_base::app);
                        //float angle[2];
                        //int idx;
                        //if (m_grefFunction)
                        //    (m_interacionModuleObj.get()->*m_grefFunction)(angle, &idx);
                        //if (!csvFile.good())    // write header if new file
                        //    csv << "imagefile, eye, vid, gaze_x, gaze_y" << std::endl;
                        //csv << frame_fname.str() << ",LR," << idx << "," << angle[0] << "," << angle[1] << std::endl;
                        //csv.close();

                        somniumCamDewarp(frame.data, frames[calib_frames], 0);

                        // dump dewarped frame
                        //cv::Mat f_dewarp(DEWARP_FRAME_HEIGHT, DEWARP_FRAME_WIDTH, CV_8UC1, frames[calib_frames]);
                        //std::stringstream frame_fname;
                        //frame_fname << std::setfill('0') << std::setw(5) << calib_frames;
                        //frame_fname << ".bmp";
                        //cv::imwrite(frame_fname.str(), f_dewarp);

                        ++calib_frames;
                    }
                }
            }

            if (state == CaptureProgramState::EstimateEyePose)
            {
                // call ET processing
                somniumCamDewarp(frame.data, dewarped_frame, 0);
                getEyePose(ei, dewarped_frame, NULL, 0);

                {
                    std::lock_guard<std::mutex> lock(m_mutex);
                    memcpy(gaze_vector, ei->gaze_vector_smooth, 3 * sizeof(float));
                }

                // dump eye pose estimation if eye sufficiently open
                if (ei->blink < 0.5f)
                {
                    printf("Eye ball position: [%f, %f, %f]   ", ei->eyeball3d[0], ei->eyeball3d[1], ei->eyeball3d[2]);
                    printf("Iris position: [%f, %f, %f]   ", ei->iris3d[0], ei->iris3d[1], ei->iris3d[2]);
                    printf("Gaze vector: [%f, %f, %f]                \r", ei->gaze_vector[0], ei->gaze_vector[1], ei->gaze_vector[2]);
                }
            }

            ++frame_idx;
        }
        else
        {
            std::cerr << "OpenCV Error: Could not capture frame." << std::endl;
            break;
        }
    }
}


void CaptureModule::stopCapture()
{
    captureRun = false;

    // wait for capture thread to stop
    while (!captureStopped)
        std::this_thread::yield();
}


// Clean up
void CaptureModule::cleanResources()
{
    releaseEyeConfig(ei);

    if (dewarped_frame)
    {
        free(dewarped_frame);
        dewarped_frame = nullptr;
    }

    for (int i=0; i<calib_frames; ++i)
        if (frames[i])
        {
            free(frames[i]);
            frames[i] = nullptr;
        }

    // Release the capture
    cap.release();
    cv::destroyAllWindows();
}
