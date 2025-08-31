
LOCAL_PATH := $(call my-dir)

include $(CLEAR_VARS)

LOCAL_MODULE := eyetracking

LOCAL_CFLAGS := -DHAVE_NEON=1 -Ofast -ftree-slp-vectorize

LOCAL_SRC_FILES := ../test.c ../EyeTracking.c ../ellipse.c ../eigen.c

LOCAL_LDLIBS := -llog

include $(BUILD_EXECUTABLE)

