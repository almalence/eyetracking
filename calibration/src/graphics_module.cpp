/* ------------------------------------------------------------------------- *\

    Eye tracking calibration: Graphics handling

    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
    SPDX-License-Identifier: GPL-3.0-only
    All rights reserved.

    Author: Bogdan Nudnenko

\* ------------------------------------------------------------------------- */

#include "graphics_module.h"
#include "shaders.h"

constexpr float BackgroundGray[] = { 0.04f, 0.06f, 0.06f, 1.0f };

// Vertex data for a full-screen quad
GLfloat vertices[] = {
     1.f,  1.f, 0.0f,
     1.f, -1.f, 0.0f,
    -1.f, -1.f, 0.0f,
    -1.f,  1.f, 0.0f
};

// Indices for the quad (two triangles)
unsigned int indices[] = {
    0, 1, 3,
    1, 2, 3
};

void CheckShader(GLuint shader) {
    GLint r = 0;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &r);
    if (r == GL_FALSE) {
        GLchar msg[4096] = {};
        GLsizei length;
        glGetShaderInfoLog(shader, sizeof(msg), &length, msg);
        THROW(Fmt("Compile shader failed: %s", msg));
    }
}

void CheckProgram(GLuint prog) {
    GLint r = 0;
    glGetProgramiv(prog, GL_LINK_STATUS, &r);
    if (r == GL_FALSE) {
        GLchar msg[4096] = {};
        GLsizei length;
        glGetProgramInfoLog(prog, sizeof(msg), &length, msg);
        THROW(Fmt("Link program failed: %s", msg));
    }
}

void GraphicsModule::initializeDevice()
{
    // Initialize the gl extensions. Note we have to open a window.
    ksDriverInstance driverInstance{};
    ksGpuQueueInfo queueInfo{};
    ksGpuSurfaceColorFormat colorFormat{ KS_GPU_SURFACE_COLOR_FORMAT_B8G8R8A8 };
    ksGpuSurfaceDepthFormat depthFormat{ KS_GPU_SURFACE_DEPTH_FORMAT_D24 };
    ksGpuSampleCount sampleCount{ KS_GPU_SAMPLE_COUNT_1 };
    if (!ksGpuWindow_Create(&m_window, &driverInstance, &queueInfo, 0, colorFormat, depthFormat, sampleCount, 640, 480, false)) {
        THROW("Unable to create GL context");
    }
    ShowWindow(m_window.hWnd, SW_SHOWMINIMIZED); // minimize window out of the way (only created for GL context)

    GLint major = 0;
    GLint minor = 0;
    glGetIntegerv(GL_MAJOR_VERSION, &major);
    glGetIntegerv(GL_MINOR_VERSION, &minor);

    glEnable(GL_DEBUG_OUTPUT);
}

void GraphicsModule::compileShaderProgram() {
    // Compile Vertex Shader
    unsigned int vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexShaderSource, nullptr);
    glCompileShader(vertexShader);
    CheckShader(vertexShader);

    // Compile Fragment Shader
    unsigned int fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentShaderSource, nullptr);
    glCompileShader(fragmentShader);
    CheckShader(fragmentShader);

    // Link Shaders
    m_circleProgram = glCreateProgram();
    glAttachShader(m_circleProgram, vertexShader);
    glAttachShader(m_circleProgram, fragmentShader);
    glLinkProgram(m_circleProgram);
    CheckProgram(m_circleProgram);

    // Set uniforms
    m_u_circleCenterLocation = glGetUniformLocation(m_circleProgram, "u_circleCenter");
    m_u_radius = glGetUniformLocation(m_circleProgram, "u_radius");
    m_u_resolution = glGetUniformLocation(m_circleProgram, "u_resolution");
    m_u_state = glGetUniformLocation(m_circleProgram, "u_state");


    // Clean up shaders as they're no longer needed after linking
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);
}

void GraphicsModule::initializeResources()
{
    glGenFramebuffers(1, &m_swapchainFramebuffer);

    compileShaderProgram();

    glGenVertexArrays(1, &m_quad_VAO);

    glGenBuffers(1, &m_quad_VBO);
    glGenBuffers(1, &m_quad_EBO);

    glBindVertexArray(m_quad_VAO);

    glBindBuffer(GL_ARRAY_BUFFER, m_quad_VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_quad_EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid*)0);
    glEnableVertexAttribArray(0);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
}

uint32_t GraphicsModule::getDepthTexture(uint32_t colorTexture)
{
    auto depthBufferIt = m_colorToDepthMap.find(colorTexture);
    if (depthBufferIt != m_colorToDepthMap.end()) {
        return depthBufferIt->second;
    }

    GLint width;
    GLint height;
    glBindTexture(GL_TEXTURE_2D, colorTexture);
    glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_WIDTH, &width);
    glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_HEIGHT, &height);

    uint32_t depthTexture;
    glGenTextures(1, &depthTexture);
    glBindTexture(GL_TEXTURE_2D, depthTexture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT32, width, height, 0, GL_DEPTH_COMPONENT, GL_FLOAT, nullptr);

    m_colorToDepthMap.insert(std::make_pair(colorTexture, depthTexture));

    return depthTexture;
}



void GraphicsModule::draw(uint32_t colorTexture, int offset_X, int offset_Y, int32_t width, int32_t height)
{

    glBindFramebuffer(GL_FRAMEBUFFER, m_swapchainFramebuffer);

    glViewport(static_cast<GLint>(offset_X),
               static_cast<GLint>(offset_Y),
               static_cast<GLsizei>(width),
               static_cast<GLsizei>(height));

    glFrontFace(GL_CW);
    glCullFace(GL_BACK);

    glDisable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);

    const uint32_t depthTexture = getDepthTexture(colorTexture);

    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, colorTexture, 0);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, depthTexture, 0);

    // Clear swapchain and depth buffer.
    glClearColor(BackgroundGray[0], BackgroundGray[1], BackgroundGray[2], BackgroundGray[3]);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

    // Set shaders and uniform variables.
    glUseProgram(m_circleProgram);

    float resolution[] = { (float)width, (float)height };
    float circleCenter[2] = {0.5f,0.5f};
    float radius = 0.005f;
    int state = 0;
    if (m_interactionFunction)
    {
        (m_interacionModuleObj.get()->*m_interactionFunction)(circleCenter, &radius, &state);
    }
    
    glUniform2fv(m_u_circleCenterLocation, 1, circleCenter);
    glUniform1f(m_u_radius, radius);
    glUniform2fv(m_u_resolution, 1, resolution);
    glUniform1i(m_u_state, state);

    glBindVertexArray(m_quad_VAO);
    glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);

    GLenum err;
    while ((err = glGetError()) != GL_NO_ERROR) {
        std::cerr << "OpenGL error: " << err << std::endl;
    }

    glUseProgram(0);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void GraphicsModule::setInteractionProgram(std::shared_ptr<InteractionModule> obj, InteractionFunction interactionFunc)
{
    m_interacionModuleObj = obj;
    m_interactionFunction = interactionFunc;
}

void GraphicsModule::cleanResources()
{
    if (m_quad_VBO != 0) {
        glDeleteBuffers(1, &m_quad_VBO);
        m_quad_VBO = 0;
    }
    if (m_quad_EBO != 0) {
        glDeleteBuffers(1, &m_quad_EBO);
        m_quad_EBO = 0;
    }
    if (m_quad_VAO != 0) {
        glDeleteVertexArrays(1, &m_quad_VAO);
        m_quad_VAO = 0;
    }
    if (m_swapchainFramebuffer != 0) {
        glDeleteFramebuffers(1, &m_swapchainFramebuffer);
        m_swapchainFramebuffer = 0;
    }
    if (m_circleProgram != 0) {
        glDeleteProgram(m_circleProgram);
        m_circleProgram = 0;
    }

    m_u_circleCenterLocation = -1;
    m_u_radius = -1;
    m_u_resolution = -1;
    m_u_state = -1;

    ksGpuWindow_Destroy(&m_window);
    m_window = {};
}

ksGpuWindow* GraphicsModule::getWindow()
{
    return &m_window;
}
