/* ------------------------------------------------------------------------- *\

    Eye tracking calibration: OpenXR code

    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
    SPDX-License-Identifier: GPL-3.0-only
    All rights reserved.

    Author: Bogdan Nudnenko

\* ------------------------------------------------------------------------- */

#pragma once

#include "common.h"
#include "openxr.h"
#include "openxr_platform.h"
#include "xr_generated_dispatch_table.h"
#include <gfxwrapper_opengl.h>
#include <functional>

class GraphicsModule;

class OpenXrSwapchainEnveropment
{
private:
    XrSwapchain m_swapchain;
    XrSwapchainCreateInfo m_createInfo;
    std::vector<XrSwapchainImageBaseHeader*> m_swapchainImages;
#ifdef XR_USE_PLATFORM_ANDROID
    std::list<std::vector<XrSwapchainImageOpenGLESKHR>> m_swapchainImageBuffers;
#else
    std::list<std::vector<XrSwapchainImageOpenGLKHR>> m_swapchainImageBuffers;
#endif
    void allocateSwapchainImages();
public:
    OpenXrSwapchainEnveropment(XrSession xrSession, XrSwapchainCreateInfo createInfo);
    XrSwapchain swapchain() const;
    XrSwapchainCreateInfo* createInfo();
    XrSwapchainImageBaseHeader* image(int index) const;
};

class OpenXRModule 
{
public:
#ifdef XR_USE_PLATFORM_ANDROID
    void createInstance(void*, void*);
#else
    void createInstance();
#endif
    void initializeSystem();
    void createSession();
    void configureReferenceSpace();
    void initOpenGlExtension();
    void createSwapchains(const std::vector<int64_t>& supportedColorSwapchainFormats);
    void completeGraphicsBinding(ksGpuWindow* m_window);
    void pollEvents(bool* exitRenderLoop);
    void renderFrame();
    void cleanResources();

    std::vector<XrFovf> getFov() const;

    typedef void (GraphicsModule::*DrawFunction)(uint32_t, int, int, int32_t, int32_t);
    void setDrawFunction(std::shared_ptr<GraphicsModule> obj, DrawFunction drawFunc);

private:
	XrInstance m_instance{XR_NULL_HANDLE};
	XrSession m_session{XR_NULL_HANDLE};
	XrSpace m_appSpace{XR_NULL_HANDLE};
	XrSystemId m_systemId{XR_NULL_SYSTEM_ID};
    XrSessionState m_sessionState{ XR_SESSION_STATE_UNKNOWN };
    XrEventDataBuffer m_eventDataBuffer{ XR_TYPE_EVENT_DATA_BUFFER };

    bool debugUtilsAvailable = false;
    PFN_xrCreateDebugUtilsMessengerEXT pfnCreateDebugUtilsMessengerEXT = nullptr;
    PFN_xrDestroyDebugUtilsMessengerEXT pfnDestroyDebugUtilsMessengerEXT = nullptr;
    XrDebugUtilsMessengerEXT debugMessenger = XR_NULL_HANDLE;

    std::vector<XrViewConfigurationView> m_configViews;
    std::vector<XrView> m_views;

    std::vector<std::shared_ptr<OpenXrSwapchainEnveropment>> m_swapchains;

#ifdef XR_USE_PLATFORM_WIN32
    XrGraphicsBindingOpenGLWin32KHR m_graphicsBinding{ XR_TYPE_GRAPHICS_BINDING_OPENGL_WIN32_KHR };
#elif XR_USE_PLATFORM_ANDROID
    XrGraphicsBindingOpenGLESAndroidKHR m_graphicsBinding{XR_TYPE_GRAPHICS_BINDING_OPENGL_ES_ANDROID_KHR};
#else
#error UNKNOWN PLATFORM
#endif

    XrVector3f m_worldObjectScale{ 10.f, 10.f, 1.f };

    std::shared_ptr<GraphicsModule> m_graphicsModuleObj = nullptr;
    DrawFunction m_drawFunction = nullptr;

    bool m_sessionRunning = false;

    XrEventDataBaseHeader* tryReadNextEvent();
    void OpenXRHandleSessionStateChangedEvent(const XrEventDataSessionStateChanged& stateChangedEvent, bool* exitRenderLoop);
};