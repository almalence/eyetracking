### Eye Tracking API Documentation

---

## Overview
The eye tracking algorithm from Almalence provides a framework for estimating 3d eye pose.
It's designed for applications requiring precise eye tracking with Near-Eye Display devices, such as virtual reality HMDs.

### Dependencies
None. The Eye Tracking algorithm is self-contained and does not rely on any 3rd-party libraries/code.

---

## API Reference

### Error codes
- `EYECAL_ALL_OK`: No error.
- `EYECAL_NO_MEMORY`: Not enough memory.
- `EYECAL_NO_PUPIL`: Not enough frames with visible pupil.
- `EYECAL_AXES_RANGE`: Eye rotational axes are out of range.

### Data Structures

#### `eye_cfg`
- **Description**: Structure for storing eye configuration and state data.
This API document provides only the data fields set during initialization and those that are relevant for obtaining the eye pose.
For more information about data fields in `eye_cfg` used during internal processing please see `EyeTracking.h` header file.
Fields in `eye_cfg` structure are initialized with `initEyeConfig()` call.

- **Fields**:
  - **Camera-related constants**:
    - `int w, h`: Dimensions of the camera frame (width, height).
    - `float camera_rot`: Rotation angle (degrees clockwise) to align the eye horizontally.
    - `float cam2hmd[3][4]`: Transformation matrix from camera to HMD coordinates.
    - `int   valid_area[6]`: valid_area: Array describing valid area to constrain search for pupil in a form of:<br>
       [ctr_x ctr_y horz_dist vert_dist diag_dist corner_brighten], where:
       - ctr_x ctr_y - area center
       - horz_dist vert_dist - distances from area center
       - specify diag_dist:
         - =horz_dist+vert_dist to have rectangular crop area
         - =horz_dist=vert_dist to have diamond-shaped crop area
         - = 0.7 * (horz_dist+vert_dist) to have nearly round crop area
       - set corner_brighten to:
         - 0 to assume equally-bright image from the camera
         - 128 to bright-up image corners by factor 1.5
         - 256 to bright-up image corners by factor 2.0,
         - etc (factor = 1+corner_brighten/256).
    - `int   eyeball_ctr[2]`: 'ideal' (designed) location of the eyeball center in camera frame, normally near the frame center.


  - **Iris-related fields**:
    - `int min_iris_r, max_iris_r`: Minimum and maximum possible iris radius in pixels.
    - `float iris_r`: Detected iris radius in pixels.
    - `float pupil_decenter_mm[2]`: Detected horizontal decenter in mm (>0 - to the right)

  - **Eye pose data (result of eye pose estimation process)**:
    - `float iris3d[3]`: 3D position of the iris center in HMD space (mm).
    - `float eyeball3d[3]`: 3D position of the eyeball center in HMD space (mm).
    - `float gaze_vector[3]`: Normalized gaze vector.
    - `float gaze_vector_smooth[3]`: Gaze vector with jitter smoothed at fixation positions

#### `eye_calib`
- **Description**: Structure for storing the calibrated eye parameters.

- **Fields**:
  - `float iris2hrot`: Distance from iris to horizontal eye rotation axis.
  - `float iris2vrot`: Distance from iris to vertical eye rotation axis.
  - `float iris2pupil`: How much pupil is closer to the eyeball center than the iris plane.
  - `float pupil_decenter[2]`: Pupil decentration (horizontal and vertical), in mm.
  - `float eyeball_ctr[2]`: 2D eyeball center on the imaging plane, detected during calibration.
  - `float pca_angle[2]`: Pupillary circular axis angles (horizontal, vertical).
  - `float pupil_refr`: Corneal bulge refraction correction coefficient (currently - fixed at 1.121).
  - `float pupil_angle`: Detected pupil elongation angle (clockwise from 12'o-clock).
  - `float pupil_stretch`: Pupil elongation (1.0 - no elongation).
  - `float gaze_alpha[2]`: Horizontal/vertical angle between optical and visual axes.

---

### Functions

#### `eye_cfg *initEyeConfig(...)`
```c
eye_cfg *initEyeConfig(
    int w,
    int h,
    int min_iris_r,
    int max_iris_r,
    int valid_area[6],
    int pup_thr,
    float camera_hfov,
    float camera_rot,
    int fps,
    float cam2hmd[6]
);
```
- **Description**: Initializes the eye pose estimator and allocates necessary resources.
- **Parameters**:
  - `w, h` (int): Camera frame dimensions.
  - `min_iris_r, max_iris_r` (int): Minimum and maximum iris radius (pixels).
  - `valid_area[6]` (int array): Defines a valid area for pupil detection. See description in `eye_cfg` section for more details.
  - `pup_thr` (int): Pupil brightness threshold. Pass -1 for auto-detection (recommended).
  - `camera_hfov` (float): Camera horizontal field of view in degrees.
  - `camera_rot` (float): How much to rotate (degree clockwise) camera frame to make eye nearly-horizontal.
  - `fps` (int): Eye tracking camera framerate.
  - `cam2hmd[6]` (float array): Camera-to-HMD coordinate space (X-Y-Z Euler rotations used, default in e.g. Blender, applied after camera_rot):
    - x/y/z rotations (degree, around horizontal axis, around vertical axis, around axis perpendicular to the sensor)
    - followed by x/y/z displacements (mm).
- **Returns**: Pointer to an initialized `eye_cfg` structure, or `NULL` if memory allocation fails.

---

#### `void applyCalibration(eye_cfg *ei, eye_calib *ec)`
```c
void applyCalibration(eye_cfg *ei, eye_calib *ec);
```
- **Description**: Apply calibrated eye parameters.
- **Parameters**:
  - `ei` (eye_cfg*): Pointer to the eye configuration.
  - `ec` (eye_calib*): Pointer to the calibration structure holding calibrated eye parameters.
- **Returns**: None.

---

#### `void releaseEyeConfig(eye_cfg *ei)`
```c
void releaseEyeConfig(eye_cfg *ei);
```
- **Description**: Frees all resources associated with the `eye_cfg` instance.
- **Parameters**:
  - `ei` (eye_cfg*): Pointer to the configuration structure to be released.
- **Returns**: None.

---

#### `void getEyePose(eye_cfg *ei, uint8_t *im, int *pup_param, int16_t *pup_edge)`
```c
void getEyePose(eye_cfg *ei, uint8_t *im, int *pup_param, int16_t *pup_edge);
```
- **Description**: Estimates the eye pose based on the input frame.
- **Parameters**:
  - `ei` (eye_cfg*): Pointer to the configuration structure containing state data.
  - `im` (uint8_t*): Pointer to the camera frame data (single channel, usually IR).
  - `pup_edge` (int): If non-zero - use externally detected list of pupil edge points (pairs of x,y coordinates)
  - `pup_param` (int16_t): Array with externally detected pupil parameters (default: NULL): [npoints ctr_x ctr_y area]:
    - npoints - number of edge points in `pup_edge` array (at least 8 should be provided).
    - ctr_x,ctr_y - approximate coordinates of the ellipse center (optional, 0 can be passed if unknown).
    - area - pupil area (optional, 0 can be passed if unknown).
- **Results**: Updates the following fields in `ei`:
  - `iris3d`: Iris center location in 3D space (mm).
  - `eyeball3d`: Eyeball center location in 3D space (mm).
  - `gaze_vector`: Normalized gaze vector.
  - `gaze_vector_smooth`: Normalized gaze vector, with jitter smoothed.

---

#### `eye_calib * calibrateEye(eye_cfg *ei, uint8_t **frames, int M, float anc[5][2], int *err)`
```c
eye_calib * calibrateEye(
    eye_cfg   *ei,
    uint8_t  **frames,
    int        M,
    float      anc[5][2],
    int       *err
    );
```
- **Description**: Perform eye calibration.
M input frames should contain eye views equally distributed between the anchor locations (approximately M/5 frames for each)
Anchor locations should be in the following order: Left of center, Right of center, Above center, Below center, Center
- **Parameters**:
  - `ei` (eye_cfg*): Pointer to the eye instance, pre-initialized with `initEyeConfig()`.
  - `frames` (uint6_t **): Array of pointers to frames captured from the eye-tracking camera.
  - `M` (int): Total number of frames.
  - `anc` (float [5][2]): Locations of the anchor points used for calibration, horizontal and vertical position in degree.
- **Output**:
  - `err`       - `EYECAL_ALL_OK` if no errors, or error code.
- **Returns**: Structure with calibrated eye parameters, to be used with `applyCalibration()`.

---

#### `void releaseCalibration(eye_calib * ec)`
```c
void releaseCalibration(eye_calib * ec);
```
- **Description**: Free memory allocated for calibration parameter structure in `calibrateEye()`.
- **Parameters**: `ec` (eye_calib*): Pointer to the calibration structure.
- **Returns**: None.

---

## Usage Examples

### Initialization and Cleanup
```c
#include "EyeTracking.h"

int main() {
    int valid_area[6] = {320, 240, 300, 200, 360, 128};
    int eyeball_ctr[2] = {320, 280};
    float cam2hmd[3][4] = { /* Transformation matrix. Usually rotation in 3d space in first 3 rows/cols, last column = x/y/z distances in mm */ };
    float gaze_alpha[2] = {-5, -4};
    float pupil_decenter[2] = {0.0, 0.0};

    eye_cfg *ei = initEyeConfig(640, 480, 80, 210, valid_area, eyeball_ctr, -1, 60.0, 0.0, cam2hmd, gaze_alpha, pupil_decenter);

    if (ei == NULL) {
        printf("Initialization failed.\n");
        return -1;
    }

    // Clean up
    releaseEyeConfig(ei);
    return 0;
}
```

### Eye Pose Estimation
```c
#include "EyeTracking.h"

void processFrame(uint8_t *frame) {
    eye_cfg *ei = /* Initialized instance */;
    getEyePose(ei, frame);

    printf("Iris 3D Position: [%f, %f, %f]\n", ei->iris3d[0], ei->iris3d[1], ei->iris3d[2]);
    printf("Gaze Vector: [%f, %f, %f]\n", ei->gaze_vector[0], ei->gaze_vector[1], ei->gaze_vector[2]);
}
```

### Notes
- Ensure the camera frame dimensions (`w`, `h`) match the actual input image.
- Use the correct transformation matrix (`cam2hmd`) for accurate HMD space mapping.
