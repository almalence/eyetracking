/* ------------------------------------------------------------------------- *\

    Eye tracking calibration: Graphics handling

    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
    SPDX-License-Identifier: GPL-3.0-only
    All rights reserved.

    Author: Bogdan Nudnenko

\* ------------------------------------------------------------------------- */

#pragma once

#include "common.h"
#include "gfxwrapper_opengl.h"

class InteractionModule;

class GraphicsModule {
public:

	const std::vector<int64_t> SupportedColorSwapchainFormats = {
#ifndef XR_USE_PLATFORM_ANDROID
		GL_RGB10_A2,
		GL_RGBA16F,
#endif
		GL_RGBA8,
		GL_RGBA8_SNORM,
	};

	void initializeDevice();
	void initializeResources();
	ksGpuWindow* getWindow();
	void draw(uint32_t colorTexture, int offset_X, int offset_Y, int32_t width, int32_t height);

	typedef void (InteractionModule::* InteractionFunction)(float*, float*, int *);
	void setInteractionProgram(std::shared_ptr<InteractionModule> obj, InteractionFunction interactionFunc);
	void cleanResources();

private:
	ksGpuWindow m_window{};

	GLuint m_quad_VBO{ 0 };
	GLuint m_quad_EBO{ 0 };
	GLuint m_quad_VAO{ 0 };

	GLuint m_swapchainFramebuffer{ 0 };
	GLuint m_circleProgram{ 0 };

	GLint m_u_circleCenterLocation { 0 };
	GLint m_u_radius { 0 };
	GLint m_u_resolution { 0 };
	GLint m_u_state { 0 };
	
	std::map<uint32_t, uint32_t> m_colorToDepthMap;

	std::shared_ptr<InteractionModule> m_interacionModuleObj = nullptr;
	InteractionFunction m_interactionFunction = nullptr;

	void compileShaderProgram();
	uint32_t getDepthTexture(uint32_t colorTexture);
};