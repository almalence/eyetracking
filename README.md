# eyetracking
Fast, Iris-based eye tracking, implemented in MATLAB and native C

# Running MATLAB code

Before running MATLAB code, combine zip volumes into single zip file:
```
cd MATLAB\testing\dataset\SomniumVR1
copy /b 01.zip.001+01.zip.002+01.zip.003 01.zip
```
## Run calibration:
in Matlab:
```
> cd testing
> calibration
```
## Run testing:
in Matlab:
```
> cd testing
> test
```
Note: test run will show intermediate processing images requiring you to hit a key to proceed to the next frame for the first 20 frames.

