# eyetracking
Fast, Iris-based eye tracking, implemented in MATLAB (reference code) and as an optimized native C implementation.


https://github.com/user-attachments/assets/391336a3-5567-4fe7-afb3-873ab0674474



Check [publication](https://ieeexplore.ieee.org/document/11189046) for the method details.

Performance:
- On a 3GHz 12th generation Intel CPU the single 640x480 frame processing takes **350µs** on a **single core**.
- On XR-specific Snapdragon XR2, 640x480 camera frame: **465µs** on a **single core**. This means that **stereo eye tracking at a 120Hz** refresh rate utilizes under **2% of the total CPU at less than 0.5ms latency**.


# Build

## Prerequisites (for eye calibration application only)

- ### MSVC 2022 (MSVC 2019 should work as well, but never tested)

- ### OpenCV
Install OpenCV library as described here:
https://docs.opencv.org/4.x/d3/d52/tutorial_windows_install.html
Follow the sections `Installation by Using the Pre-built Libraries` and `Set the OpenCV environment variable and add it to the systems path`
Make sure OpenCV_DIR points to where the `OpenCVConfig.cmake` is!

*Note:* OpenCV library is only needed to demonstrate access to the camera frames. The eye tracking algorithm itself is self-contained and does not depend on any 3rd-party libraries.

- ### OpenXR
OpenXR SDK is automatically fetched and built during the cmake configure step.

- ### Eye tracking camera

Code was tested with the Somnium VR1 eye tracking camera.

For camera accessible via OpenCV - adjust the default camera parameters in `src/capture_module.cpp/CaptureModule::startCapture()`.
For other cases - modify `CaptureModule::startCapture()` and `CaptureModule::runCapture()` methods in `src/capture_module.cpp`.

Also, modify frame parameters in `src/capture_module.h` and camera parameters in call to `initEyeConfig()` in `src/main.cpp`.

## Command-line test/calibration tools

The easiest way is to open the folder `SDK` with VS Code and build from there with 'Cmake Tools' extension.
Cmake build from command line should work as well, but never tested.

Use gcc 12.0 or newer, older versions and other compilers are not guaranteed to perform auto-vectorization effectively.

To build for Android (you will need android NDK):
```
cd SDK
ndk-build.cmd
```
This will produce command-line test executable: `libs\arm64-v8a\eyetracking`

## Eye-tracking library

The eye tracking library is built together with the command-line tools.

To use newly built library in calibration application, copy  

`SDK/out/build/release/libEyeTrack.a` into `calibration/et-lib/`

and rename it to `libEyeTrack.lib`  
Also update `EyeTracking.h` in the same place, if needed.

## Calibration application, on Windows
```
cmake -S . -B ./CMAKE_BUILD -D CMAKE_BUILD_TYPE=Release
```
Open solution in MSVC: `CMAKE_BUILD/EyeTrackerCalibrationApp.sln`, select 'Release' solution configuration, build ALL_BUILD project


## Calibration application, Android
```
cmake -S . -B ./CMAKE_BUILD -G Ninja -DANDROID_ABI=arm64-v8a -DANDROID_PLATFORM=android-24 -DANDROID_NDK=[ndk folder] -DCMAKE_TOOLCHAIN_FILE=[ndk folder]/build/cmake/android.toolchain.cmake -DCMAKE_MAKE_PROGRAM=Android\Sdk\cmake\3.22.1\bin\ninja.exe
```
Note: Camera capture was not implemented in Android, source code modifications needed (see section [Eye tracking camera](#eye-tracking-camera) above)


# Run

## Running MATLAB code

Before running MATLAB code, combine zip volumes with captured frames into single zip file:
```
cd MATLAB\testing\dataset\SomniumVR1
copy /b 01.zip.001+01.zip.002+01.zip.003 01.zip
```
### Run calibration
in Matlab:
```
> cd testing
> calibration
```
### Run eye pose estimation test
in Matlab:
```
> cd testing
> test
```
Note: test run will show intermediate processing images requiring you to hit a key to proceed to the next frame for the first 20 frames.


## Running command-line C code
```
cd SDK\out
build\release\EyeCalibration.exe calib\00000.bmp 199 540 400
build\release\EyeTracking.exe 01\00000.bmp 1000 540 400
```

## Running calibration app on HMD

### Windows (Somnium VR1)

To turn on eye-tracking IR-leds on Somnium VR1:
- Start Somnium VR Tool
- Go to Settings / Eyes
- Perform the Eye Tracking calibration (has nothing to do with the eye tracking of this project, only to enable eye-tracking switch in settings).
- Start Windows 'Camera' application and switch to Somnium Eye Tracking camera (frames will look green).
- Enable eye-tracking switch in the Somnium VR Tool settings.
- Exit windows 'Camera' app.

The above sequence of steps will keep IR-leds turned on, while allowing 3rd party apps to access the eye-tracking camera at the same time.

Run: CMAKE_BUILD\Release\EyeTrackerCalibrationApp.exe 


### Android HMD

Use Android Studio to build and run the calibration app.
Note: only graphic (calibration dots) is displayed in Android version, camera capture was never tested, and no calibration code implemented (should be an easy copy/paste from Windows version).


# Use in your own projects


Calling sequence to perform eye calibration:  
```
    eye_cfg *ei = initEyeConfig(...);

    // ... show calibration points on screen,
    // capture eye tracking camera frames ...

    eye_calib *ec = calibrateEye(...);
    // ec holds calibrated eye parameters here

    releaseEyeConfig(ei);
```

Calling sequence to perform eye tracking:  

```
    eye_cfg *ei = initEyeConfig(...);
    applyCalibration(ei, &ec);

    while( frames )
    {
        // process captured frame to get eye pose
        getEyePose(ei, frame);
        // eye pose is in ei->iris3d, ei->eyeball3d, ei->gaze_vector
    }

    releaseEyeConfig(ei);
```


For more detailed information, take a look at:
- [API documentation](docs/api-documentation-for-eye-tracking.md)
- [SDK/EyeTracking.h](SDK/EyeTracking.h)
- [Command-line calibration code](SDK/calibration.c)
- [Command-line eye pose estimation test](SDK/test.c)
