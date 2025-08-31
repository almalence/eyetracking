#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

#include "EyeTracking.h"

// to compare matlab 'ground truth' results with optimized C version
#include "mlab_eyepos.inc"
#define mlab_gnd mlab_01


void print_usage()
{
	printf (
		"Usage: test in001.bmp N Width Height\n"
		"Will process N frames from IR eye-tracking camera: in001.bmp in002.bmp .. in[N].bmp\n"
		"Width and Height are frame width (stride) and height in pixels\n"
		"BMPs should be 8bit, monochrome\n"
		"File names should have sequential numbering just before extension\n"
		);
}


void print_stats(float *iris3d, float *eyeball3d, float *gaze_vector, int N)
{
	float iris_diff = 0;
	float eye_diff = 0;
	float gaze_diff = 0;
	float iris_maxd = 0;
	float eye_maxd = 0;
	float gaze_maxd = 0;
	float mean_h=0, mean_v=0;
	int Nd=0;
	for (int i=0; i<N; ++i)
	{
		printf ("%03d : ", i);
		printf ("iris: [%6.2f %6.2f %6.2f]  eyeball: [%6.2f %6.2f %6.2f]  gazevec: [%5.2f %5.2f %5.2f]\n",
		    iris3d[3*i+0], iris3d[3*i+1], iris3d[3*i+2],
		    eyeball3d[3*i+0], eyeball3d[3*i+1], eyeball3d[3*i+2],
		    gaze_vector[3*i+0], gaze_vector[3*i+1], gaze_vector[3*i+2] );

		// excluding gaze switching periods
		if ( (i>=8) && (mlab_gnd[(i-8)*11+9] == mlab_gnd[i*11+9]) && (mlab_gnd[(i-8)*11+10] == mlab_gnd[i*11+10]) )
		{
			// compare to matlab results
			iris_diff += sqrtf( (iris3d[3*i+0]-mlab_gnd[i*11+0])*(iris3d[3*i+0]-mlab_gnd[i*11+0]) + (iris3d[3*i+1]-mlab_gnd[i*11+1])*(iris3d[3*i+1]-mlab_gnd[i*11+1]) + (iris3d[3*i+2]-mlab_gnd[i*11+2])*(iris3d[3*i+2]-mlab_gnd[i*11+2]) );
			eye_diff += sqrtf( (eyeball3d[3*i+0]-mlab_gnd[i*11+3])*(eyeball3d[3*i+0]-mlab_gnd[i*11+3]) + (eyeball3d[3*i+1]-mlab_gnd[i*11+4])*(eyeball3d[3*i+1]-mlab_gnd[i*11+4]) + (eyeball3d[3*i+2]-mlab_gnd[i*11+5])*(eyeball3d[3*i+2]-mlab_gnd[i*11+5]) );
			gaze_diff += sqrtf( (gaze_vector[3*i+0]-mlab_gnd[i*11+6])*(gaze_vector[3*i+0]-mlab_gnd[i*11+6]) + (gaze_vector[3*i+1]-mlab_gnd[i*11+7])*(gaze_vector[3*i+1]-mlab_gnd[i*11+7]) + (gaze_vector[3*i+2]-mlab_gnd[i*11+8])*(gaze_vector[3*i+2]-mlab_gnd[i*11+8]) );

			iris_maxd = fmaxf(iris_maxd, fabsf(iris3d[3*i+0]-mlab_gnd[i*11+0])); iris_maxd = fmaxf(iris_maxd, fabsf(iris3d[3*i+1]-mlab_gnd[i*11+1])); iris_maxd = fmaxf(iris_maxd, fabsf(iris3d[3*i+2]-mlab_gnd[i*11+2]));
			eye_maxd = fmaxf(eye_maxd, fabsf(eyeball3d[3*i+0]-mlab_gnd[i*11+3])); eye_maxd = fmaxf(eye_maxd, fabsf(eyeball3d[3*i+1]-mlab_gnd[i*11+4])); eye_maxd = fmaxf(eye_maxd, fabsf(eyeball3d[3*i+2]-mlab_gnd[i*11+5]));
			gaze_maxd = fmaxf(gaze_maxd, fabsf(gaze_vector[3*i+0]-mlab_gnd[i*11+6])); gaze_maxd = fmaxf(gaze_maxd, fabsf(gaze_vector[3*i+1]-mlab_gnd[i*11+7])); gaze_maxd = fmaxf(gaze_maxd, fabsf(gaze_vector[3*i+2]-mlab_gnd[i*11+8]));

			// compare to ground-truth gaze angles
			float gaze_deg_h = 180/(float)M_PI * asinf(-gaze_vector[3*i+0]);
			float gaze_deg_v = 180/(float)M_PI * asinf(-gaze_vector[3*i+1]);

			mean_h += gaze_deg_h - mlab_gnd[i*11+9];
			mean_v += gaze_deg_v - mlab_gnd[i*11+10];

			++Nd;
		}
	}
	mean_h /= Nd;
	mean_v /= Nd;

	float deg_err_h=0, deg_err_v=0;
	for (int i=0; i<N; ++i)
	{
		// excluding gaze switching periods
		if ( (i>=8) && (mlab_gnd[(i-8)*11+9] == mlab_gnd[i*11+9]) && (mlab_gnd[(i-8)*11+10] == mlab_gnd[i*11+10]) )
		{
			float gaze_deg_h = 180/(float)M_PI * asinf(-gaze_vector[3*i+0]);
			float gaze_deg_v = 180/(float)M_PI * asinf(-gaze_vector[3*i+1]);

			deg_err_h += (gaze_deg_h-mean_h - mlab_gnd[i*11+9])  * (gaze_deg_h-mean_h - mlab_gnd[i*11+9]);
			deg_err_v += (gaze_deg_v-mean_v - mlab_gnd[i*11+10]) * (gaze_deg_v-mean_v - mlab_gnd[i*11+10]);
		}
	}
	deg_err_h = sqrtf(deg_err_h/Nd);
	deg_err_v = sqrtf(deg_err_v/Nd);

	printf ("MSE to matlab version, iris: %7.3f  eyeball: %7.3f  gazevec: %7.3f\n", iris_diff/Nd, eye_diff/Nd, gaze_diff/Nd);
	printf ("Max deviation to matlab version, iris: %7.3f  eyeball: %7.3f  gazevec: %7.3f\n", iris_maxd, eye_maxd, gaze_maxd);
	printf("Gaze mean degree diff: [%5.1f %5.1f]\n", mean_h, mean_v);
	printf("Gaze MSE error excluding mean degree: [%6.2f %6.2f]\n", deg_err_h, deg_err_v);
}


int main (int ac, char *av[])
{
	if (ac!=5) { print_usage(); return 0; }

	char *ext = strchr(av[1],'.');
	if (ext==NULL) { print_usage(); return 1; }

	// detect numerical field in filename
	int extpos = (int)(ext-av[1]);
	int i, numfld;
	for (numfld=extpos-1; numfld>=0; --numfld)
		if ((av[1][numfld]<'0') || (av[1][numfld]>'9')) break;
	++numfld;
	int numlen = extpos-numfld;
	// set up filename sequencing
	int st = atoi(av[1]+numfld);
	int N = atoi(av[2]);
	int w = atoi(av[3]);
	int h = atoi(av[4]);
	char fnstr[256], lenfmt[256];
	sprintf(lenfmt, "%%0%dd", numlen);

	// initializing with values reasonable for off-axis camera, 540x400 camera frame
    int valid_area[6] = { 280, 200,  230, 200,  280,  128 };
    float cam2hmd[6] = { 0, 22, 0,  25, 0, 48 };
    // Negative sign in camera HFOV indicates chief ray converging from camera entrance pupil towards eye
    // This may happen in cases where camera is behind HMD lens and image is formed through some sort of concave mirror surface
	eye_cfg *ei = initEyeConfig(w, h, 65, 120, valid_area, -1, -30, 0, 50, cam2hmd);

	// Calibrated data
	eye_calib ec = { .iris2hrot = 12.19, .iris2vrot =  10.15, .iris2pupil = 0.30, .pupil_decenter = { 0.72,-0.21}, .eyeball_ctr = {-102,  44}, .pca_angle = {   0.00,   4.75}, .pupil_refr = 1.121, .pupil_angle =   157.50, .pupil_stretch = 0.945, .gaze_alpha = {  -0.53,  -2.26} };
	applyCalibration(ei, &ec);

	if (ei==NULL) { printf("Unable to initialize eye instance, not enough memory\n"); return 1; }

	uint8_t **cam=(uint8_t**)malloc(N*sizeof(uint8_t*));
	if (cam==NULL) { printf("Not enough memory to hold camera frame\n"); return 2; }
	for (i=0; i<N; ++i)
	{
		cam[i]=(uint8_t*)malloc(w*h);
		if (cam[i]==NULL) { printf("Not enough memory to hold camera frame\n"); return 2; }
	}

	// read camera frames
	for (i=0; i<N; ++i)
	{
		strncpy(fnstr, av[1], numfld);
		sprintf(fnstr+numfld, lenfmt, i+st);
		strncpy(fnstr+extpos, av[1]+extpos, 256-extpos);

		FILE *f = fopen(fnstr, "rb");
		size_t n;
		if (f==NULL) { printf("Can't open input frame\n"); return 3; }
		fseek(f, -w*h, SEEK_END);	// header skip, assuming nothing after pixel values
		for (int y=h-1; y>=0; --y)	// bmp data is bottom-up
			n = fread(cam[i]+y*w, w, 1, f);
		fclose(f);
	}

	float *iris3d = malloc(3*N*sizeof(float));
	float *eyeball3d = malloc(3*N*sizeof(float));
	float *gaze_vector = malloc(3*N*sizeof(float));
	uint8_t *tmp = malloc(w*h);

	printf ("Input frames read complete...\n");

	// run eye-tracking processing 10 times for time-measurement purposes
	clock_t begin = clock();
	for (int t=0; t<10; ++t)
	{
		for (i=0; i<N; ++i)
		{
			// using tmp as frame holder, as getEyePose will detroy input (filter and brighten)
			memcpy(tmp, cam[i], w*h);
			getEyePose(ei, tmp, NULL, 0);

			// keep data to dump later
			if (t==0)
			{
				iris3d[3*i]      = ei->iris3d[0];      iris3d[3*i+1]      = ei->iris3d[1];      iris3d[3*i+2]      = ei->iris3d[2];
				eyeball3d[3*i]   = ei->eyeball3d[0];   eyeball3d[3*i+1]   = ei->eyeball3d[1];   eyeball3d[3*i+2]   = ei->eyeball3d[2];
				gaze_vector[3*i] = ei->gaze_vector[0]; gaze_vector[3*i+1] = ei->gaze_vector[1]; gaze_vector[3*i+2] = ei->gaze_vector[2];
			}
		}
	}
	clock_t end = clock();

	// print-out results
	print_stats(iris3d, eyeball3d, gaze_vector, N);

	printf ("\n\nDone...   Time: %5.2f\n", (float)(end-begin)/CLOCKS_PER_SEC);

	releaseEyeConfig(ei);

	free(iris3d);
	free(eyeball3d);
	free(gaze_vector);
	free(tmp);
	for (i=0; i<N; ++i)
		free(cam[i]);
	free(cam);

	return 0;
}
