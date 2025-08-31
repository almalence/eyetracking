#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

#include "EyeTracking.h"

#define IMG_DEBUG

#ifdef IMG_DEBUG
	size_t dbgimlen;
	uint8_t *dbgim;
#endif

void print_usage()
{
	printf (
		"Usage: calibration in001.bmp N Width Height\n"
		"Will obtain eye calibration parameters from N frames from IR eye-tracking camera: in001.bmp in002.bmp .. in[N].bmp\n"
		"Width and Height are frame width (stride) and height in pixels\n"
		"BMPs should be 8bit, monochrome\n"
		"File names should have sequential numbering just before extension\n"
		);
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


    // initializing with values reasonable for off-axis camera, 640x480 camera frame
    int valid_area[6] = { 280, 200,  230, 200,  280,  128 };
    float cam2hmd[6] = { 0, 22, 0,  25, 0, 48 };
    // Negative sign in camera HFOV indicates chief ray converging from camera entrance pupil towards eye
    // This may happen in cases where camera is behind HMD lens and image is formed through some sort of concave mirror surface
	eye_cfg *ei = initEyeConfig(w, h, 65, 120, valid_area, -1, -30, 0, 30, cam2hmd);

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
		int32_t wi, hi;
		fseek(f, 0x12, SEEK_SET);	// get width/height from the header
		fread(&wi, 4, 1, f);
		fread(&hi, 4, 1, f);
		fseek(f, -wi*hi, SEEK_END);	// header skip, assuming nothing after pixel values
		for (int y=hi-1; y>=0; --y)	// bmp data is bottom-up
			n = fread(cam[i]+y*wi, wi, 1, f);

#ifdef IMG_DEBUG
		if (i==0)
		{
			fseek(f, -wi*hi, SEEK_END);
			dbgimlen = ftell(f)+w*h;
			dbgim = (uint8_t*)malloc(dbgimlen);
			fseek(f, 0, SEEK_SET);
			n = fread(dbgim, dbgimlen, 1, f);
			*(uint32_t*)(dbgim+0x12) = w;
			*(uint32_t*)(dbgim+0x16) = h;
		}
#endif

		fclose(f);
	}

	printf ("Input frames read complete...\n");

    // ToDo: use tmp as frame holder within calibration, as getEyePose will detroy input (filter and brighten)
    float anc[5][2]={{-15,0}, {15,0}, {0,-15}, {0,15}, {0,0}};
	int err;
    eye_calib *ec = calibrateEye(ei, cam, N, anc, &err);

	if (err != EYECAL_ALL_OK)
	{
		char calibrationErrorText[EYECAL_MAX_ERR][256] = {"",
			"Not enough memory",
			"Not enough frames with visible pupil",
			"Eye rotational axes are out of range"
		};
		printf("%s\n", calibrationErrorText[err]);
	}
	else
	{
		// print-out results
		printf( "eye_calib ec = { "
				".iris2hrot = %5.2f, "
				".iris2vrot = %5.2f, "
				".iris2pupil = %4.2f, "
				".pupil_decenter = {%5.2f,%5.2f}, "
				".eyeball_ctr = {%4.0f,%4.0f}, "
				".pca_angle = {%7.2f,%7.2f}, "
				".pupil_refr = %5.3f, "
				".pupil_angle = %8.2f, "
				".pupil_stretch = %5.3f, "
				".gaze_alpha = {%7.2f,%7.2f} }",
				ec->iris2hrot, ec->iris2vrot, ec->iris2pupil,
				ec->pupil_decenter[0], ec->pupil_decenter[1], ec->eyeball_ctr[0], ec->eyeball_ctr[1], ec->pca_angle[0], ec->pca_angle[1],
				ec->pupil_refr, ec->pupil_angle, ec->pupil_stretch, ec->gaze_alpha[0], ec->gaze_alpha[1] );
	}

    releaseEyeConfig(ei);
	releaseCalibration(ec);
	
	// Note: cam[], tmp here are freed by the system upon exit

	return 0;
}
