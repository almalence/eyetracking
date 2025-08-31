/* ------------------------------------------------------------------------- *\

    Eye tracking calibration: OpenXR code

    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
    SPDX-License-Identifier: GPL-3.0-only
    All rights reserved.

    Author: Bogdan Nudnenko

\* ------------------------------------------------------------------------- */

#include "openxr_module.h"
#include <xr_linear.h>

std::string XrResultToString(XrResult result) {
    switch (result) {
    case XR_SUCCESS: return "XR_SUCCESS";
    case XR_ERROR_INITIALIZATION_FAILED: return "XR_ERROR_INITIALIZATION_FAILED";
    case XR_ERROR_RUNTIME_FAILURE: return "XR_ERROR_RUNTIME_FAILURE";
    case XR_ERROR_OUT_OF_MEMORY: return "XR_ERROR_OUT_OF_MEMORY";
    case XR_ERROR_POSE_INVALID: return "XR_ERROR_POSE_INVALID";
    default: return "UNKNOWN_ERROR";
    }
}

bool CheckXrResult(XrResult result, const char* funcName, const char* fileAndLine) {
    if (XR_SUCCEEDED(result)) {
        return true;
    }
    else {
        std::cerr << "OpenXR Error: " << XrResultToString(result)
            << " AT " << funcName << " IN " << fileAndLine << std::endl;
        return false;
    }
}

#define CHECK_XR_RESULT(x) \
    CheckXrResult((x), #x, FILE_AND_LINE);

XrPosef PoseIdentity() {
    XrPosef t{};
    t.orientation.w = 1;
    return t;
}

XrPosef PoseTranslation(const XrVector3f& translation) {
    XrPosef t = PoseIdentity();
    t.position = translation;
    return t;
}


void OpenXrSwapchainEnveropment::allocateSwapchainImages()
{
    uint32_t imageCount;
    CHECK_XR_RESULT(xrEnumerateSwapchainImages(m_swapchain, 0, &imageCount, nullptr));

#ifdef XR_USE_PLATFORM_ANDROID
    std::vector<XrSwapchainImageOpenGLESKHR> swapchainImageBuffer(imageCount, { XR_TYPE_SWAPCHAIN_IMAGE_OPENGL_ES_KHR });
    for (XrSwapchainImageOpenGLESKHR& image : swapchainImageBuffer) {
        image.type = XR_TYPE_SWAPCHAIN_IMAGE_OPENGL_ES_KHR;
        m_swapchainImages.push_back(reinterpret_cast<XrSwapchainImageBaseHeader*>(&image));
    }
#else
    std::vector<XrSwapchainImageOpenGLKHR> swapchainImageBuffer(imageCount, { XR_TYPE_SWAPCHAIN_IMAGE_OPENGL_KHR });
    for (XrSwapchainImageOpenGLKHR& image : swapchainImageBuffer) {
        image.type = XR_TYPE_SWAPCHAIN_IMAGE_OPENGL_KHR;
        m_swapchainImages.push_back(reinterpret_cast<XrSwapchainImageBaseHeader*>(&image));
    }
#endif
    // Keep the buffer alive by moving it into the list of buffers.
    m_swapchainImageBuffers.push_back(std::move(swapchainImageBuffer));

    CHECK_XR_RESULT(xrEnumerateSwapchainImages(m_swapchain, imageCount, &imageCount, m_swapchainImages[0]));
}

OpenXrSwapchainEnveropment::OpenXrSwapchainEnveropment(XrSession xrSession, XrSwapchainCreateInfo createInfo)
    : m_swapchain(XR_NULL_HANDLE), m_createInfo(createInfo)
{
    CHECK_XR_RESULT(xrCreateSwapchain(xrSession, &createInfo, &m_swapchain));
    allocateSwapchainImages();
    m_createInfo = createInfo;
}

XrSwapchain OpenXrSwapchainEnveropment::swapchain() const {
    return m_swapchain;
}

XrSwapchainCreateInfo* OpenXrSwapchainEnveropment::createInfo()
{
    return &m_createInfo;
}

XrSwapchainImageBaseHeader* OpenXrSwapchainEnveropment::image(int index) const {
    return m_swapchainImages[index];
}


#ifdef XR_USE_PLATFORM_ANDROID
void OpenXRModule::createInstance(void* applicationVM, void* applicationActivity)
{
    std::vector<const char*> extensions = { XR_KHR_OPENGL_ES_ENABLE_EXTENSION_NAME, XR_EXT_DEBUG_UTILS_EXTENSION_NAME, XR_KHR_ANDROID_CREATE_INSTANCE_EXTENSION_NAME };
#else
void OpenXRModule::createInstance()
{
    std::vector<const char*> extensions = { XR_KHR_OPENGL_ENABLE_EXTENSION_NAME, XR_EXT_DEBUG_UTILS_EXTENSION_NAME };
#endif

#ifdef XR_USE_PLATFORM_ANDROID
    XrInstanceCreateInfoAndroidKHR instanceCreateInfoAndroid;
    instanceCreateInfoAndroid = {XR_TYPE_INSTANCE_CREATE_INFO_ANDROID_KHR};
    instanceCreateInfoAndroid.applicationVM = applicationVM;
    instanceCreateInfoAndroid.applicationActivity = applicationActivity;

    XrInstanceCreateInfo createInfo{ XR_TYPE_INSTANCE_CREATE_INFO };
    createInfo.next = &instanceCreateInfoAndroid;
    createInfo.enabledExtensionCount = (uint32_t)extensions.size();
    createInfo.enabledExtensionNames = extensions.data();

    strcpy(createInfo.applicationInfo.applicationName, "Almalence-EyeTracking-Example");
    createInfo.applicationInfo.apiVersion = XR_API_VERSION_1_0;

    CHECK_XR_RESULT(xrCreateInstance(&createInfo, &m_instance))
#else
    XrInstanceCreateInfo createInfo{ XR_TYPE_INSTANCE_CREATE_INFO };
    createInfo.next = nullptr;  // Needs to be set on Android.
    createInfo.enabledExtensionCount = (uint32_t)extensions.size();
    createInfo.enabledExtensionNames = extensions.data();

    strcpy(createInfo.applicationInfo.applicationName, "Almalence-EyeTracking-Example");
    createInfo.applicationInfo.apiVersion = XR_API_VERSION_1_0;

    CHECK_XR_RESULT(xrCreateInstance(&createInfo, &m_instance));
#endif
}

void OpenXRModule::initializeSystem()
{
    XrSystemGetInfo systemInfo{ XR_TYPE_SYSTEM_GET_INFO };
    systemInfo.formFactor = XR_FORM_FACTOR_HEAD_MOUNTED_DISPLAY;
    CHECK_XR_RESULT(xrGetSystem(m_instance, &systemInfo, &m_systemId));
}

void OpenXRModule::createSession()
{
    XrSessionCreateInfo createInfo{ XR_TYPE_SESSION_CREATE_INFO };
    createInfo.next = reinterpret_cast<const XrBaseInStructure*>(&m_graphicsBinding);
    createInfo.systemId = m_systemId;
    CHECK_XR_RESULT(xrCreateSession(m_instance, &createInfo, &m_session));
}

void OpenXRModule::configureReferenceSpace()
{
    XrReferenceSpaceCreateInfo referenceSpaceCreateInfo{ XR_TYPE_REFERENCE_SPACE_CREATE_INFO };
    referenceSpaceCreateInfo.poseInReferenceSpace = PoseIdentity();
    referenceSpaceCreateInfo.referenceSpaceType = XR_REFERENCE_SPACE_TYPE_VIEW;
    CHECK_XR_RESULT(xrCreateReferenceSpace(m_session, &referenceSpaceCreateInfo, &m_appSpace));
}

void OpenXRModule::initOpenGlExtension()
{
#ifdef XR_USE_PLATFORM_ANDROID
    PFN_xrGetOpenGLESGraphicsRequirementsKHR pfnGetOpenGLESGraphicsRequirementsKHR = nullptr;
	CHECK_XR_RESULT(xrGetInstanceProcAddr(m_instance, "xrGetOpenGLESGraphicsRequirementsKHR",
		reinterpret_cast<PFN_xrVoidFunction*>(&pfnGetOpenGLESGraphicsRequirementsKHR)));

    XrGraphicsRequirementsOpenGLESKHR graphicsRequirements{ XR_TYPE_GRAPHICS_REQUIREMENTS_OPENGL_ES_KHR };
	CHECK_XR_RESULT(pfnGetOpenGLESGraphicsRequirementsKHR(m_instance, m_systemId, &graphicsRequirements));
#else
    PFN_xrGetOpenGLGraphicsRequirementsKHR pfnGetOpenGLGraphicsRequirementsKHR = nullptr;
	CHECK_XR_RESULT(xrGetInstanceProcAddr(m_instance, "xrGetOpenGLGraphicsRequirementsKHR",
		reinterpret_cast<PFN_xrVoidFunction*>(&pfnGetOpenGLGraphicsRequirementsKHR)));

	XrGraphicsRequirementsOpenGLKHR graphicsRequirements{ XR_TYPE_GRAPHICS_REQUIREMENTS_OPENGL_KHR };
	CHECK_XR_RESULT(pfnGetOpenGLGraphicsRequirementsKHR(m_instance, m_systemId, &graphicsRequirements));
#endif
}

void OpenXRModule::createSwapchains(const std::vector<int64_t>& supportedColorSwapchainFormats)
{
    CHECK(m_session != XR_NULL_HANDLE);
    CHECK(m_swapchains.empty());
    CHECK(m_configViews.empty());

    XrSystemProperties systemProperties{ XR_TYPE_SYSTEM_PROPERTIES };
    CHECK_XR_RESULT(xrGetSystemProperties(m_instance, m_systemId, &systemProperties));

    uint32_t viewCount;
    CHECK_XR_RESULT(xrEnumerateViewConfigurationViews(m_instance, m_systemId, XR_VIEW_CONFIGURATION_TYPE_PRIMARY_STEREO, 0, &viewCount, nullptr));
    m_configViews.resize(viewCount, { XR_TYPE_VIEW_CONFIGURATION_VIEW });
    CHECK_XR_RESULT(xrEnumerateViewConfigurationViews(m_instance, m_systemId, XR_VIEW_CONFIGURATION_TYPE_PRIMARY_STEREO, viewCount,
        &viewCount, m_configViews.data()));

    uint32_t swapchainFormatCount;
    CHECK_XR_RESULT(xrEnumerateSwapchainFormats(m_session, 0, &swapchainFormatCount, nullptr));
    std::vector<int64_t> swapchainFormats(swapchainFormatCount);
    CHECK_XR_RESULT(xrEnumerateSwapchainFormats(m_session, (uint32_t)swapchainFormats.size(), &swapchainFormatCount,
        swapchainFormats.data()));

    CHECK(swapchainFormatCount == swapchainFormats.size());

    auto swapchainFormatIt = std::find_first_of(swapchainFormats.begin(), swapchainFormats.end(), std::begin(supportedColorSwapchainFormats),
            std::end(supportedColorSwapchainFormats));
    if (swapchainFormatIt == swapchainFormats.end()) {
        THROW("No runtime swapchain format supported for color swapchain");
    }

    m_views.resize(viewCount, { XR_TYPE_VIEW });

    if (viewCount > 0) {
        for (uint32_t i = 0; i < viewCount; i++) {
            const XrViewConfigurationView& vp = m_configViews[i];

            XrSwapchainCreateInfo swapchainCreateInfo{ XR_TYPE_SWAPCHAIN_CREATE_INFO };
            swapchainCreateInfo.arraySize = 1;
            swapchainCreateInfo.format = { *swapchainFormatIt };
            swapchainCreateInfo.width = vp.recommendedImageRectWidth;
            swapchainCreateInfo.height = vp.recommendedImageRectHeight;
            swapchainCreateInfo.mipCount = 1;
            swapchainCreateInfo.faceCount = 1;
            swapchainCreateInfo.sampleCount = 1;
            swapchainCreateInfo.usageFlags = XR_SWAPCHAIN_USAGE_SAMPLED_BIT | XR_SWAPCHAIN_USAGE_COLOR_ATTACHMENT_BIT;
            auto swapchain = std::make_shared<OpenXrSwapchainEnveropment>(m_session, swapchainCreateInfo);
            m_swapchains.push_back(swapchain);
        }
    }
}

void OpenXRModule::OpenXRHandleSessionStateChangedEvent(const XrEventDataSessionStateChanged& stateChangedEvent, bool* exitRenderLoop)
{
    m_sessionState = stateChangedEvent.state;

    if ((stateChangedEvent.session != XR_NULL_HANDLE) && (stateChangedEvent.session != m_session)) {
        std::cerr << "XrEventDataSessionStateChanged for unknown session" << std::endl;
        return;
    }

    switch (m_sessionState) {
        case XR_SESSION_STATE_READY: {
            CHECK(m_session != XR_NULL_HANDLE);
            XrSessionBeginInfo sessionBeginInfo{ XR_TYPE_SESSION_BEGIN_INFO };
            sessionBeginInfo.primaryViewConfigurationType = { XR_VIEW_CONFIGURATION_TYPE_PRIMARY_STEREO };
             CHECK_XR_RESULT(xrBeginSession(m_session, &sessionBeginInfo));
            m_sessionRunning = true;
            break;
        }
        case XR_SESSION_STATE_STOPPING: {
            CHECK(m_session != XR_NULL_HANDLE);
            m_sessionRunning = false;
            CHECK_XR_RESULT(xrEndSession(m_session));
            break;
        }
        case XR_SESSION_STATE_EXITING: {
            *exitRenderLoop = true;
            break;
        }
        case XR_SESSION_STATE_LOSS_PENDING: {
            *exitRenderLoop = true;
            break;
        }
        default:
            break;
        }
}

XrEventDataBaseHeader* OpenXRModule::tryReadNextEvent()
{
    auto* baseHeader = reinterpret_cast<XrEventDataBaseHeader*>(&m_eventDataBuffer);
    *baseHeader = { XR_TYPE_EVENT_DATA_BUFFER };
    const XrResult xr = xrPollEvent(m_instance, &m_eventDataBuffer);
    if (xr == XR_SUCCESS) {
        if (baseHeader->type == XR_TYPE_EVENT_DATA_EVENTS_LOST) {
            const auto* const eventsLost = reinterpret_cast<const XrEventDataEventsLost*>(baseHeader);
            std::cerr << Fmt("%d events lost", eventsLost) << std::endl;
        }

        return baseHeader;
    }
    CHECK_XR_RESULT(xr);
    return nullptr;
}

void OpenXRModule::pollEvents(bool* exitRenderLoop)
{
    *exitRenderLoop = false;

    // Process all pending messages.
    while (const XrEventDataBaseHeader* event = tryReadNextEvent()) {
        switch (event->type) {
            case XR_TYPE_EVENT_DATA_INSTANCE_LOSS_PENDING: {
                const auto& instanceLossPending = *reinterpret_cast<const XrEventDataInstanceLossPending*>(event);
                std::cerr << Fmt("XrEventDataInstanceLossPending by %lld", instanceLossPending.lossTime) << std::endl;
                *exitRenderLoop = true;
                return;
            }
            case XR_TYPE_EVENT_DATA_SESSION_STATE_CHANGED: {
                auto sessionStateChangedEvent = *reinterpret_cast<const XrEventDataSessionStateChanged*>(event);
                OpenXRHandleSessionStateChangedEvent(sessionStateChangedEvent, exitRenderLoop);
                break;
            }
            case XR_TYPE_EVENT_DATA_INTERACTION_PROFILE_CHANGED:
                break;
            case XR_TYPE_EVENT_DATA_REFERENCE_SPACE_CHANGE_PENDING:
            default: {
                break;
            }
        }
    }
}

void OpenXRModule::renderFrame()
{
    CHECK(m_session != XR_NULL_HANDLE);

    XrFrameWaitInfo frameWaitInfo{ XR_TYPE_FRAME_WAIT_INFO };
    XrFrameState frameState{ XR_TYPE_FRAME_STATE };
    CHECK_XR_RESULT(xrWaitFrame(m_session, &frameWaitInfo, &frameState));

    XrFrameBeginInfo frameBeginInfo{ XR_TYPE_FRAME_BEGIN_INFO };
    CHECK_XR_RESULT(xrBeginFrame(m_session, &frameBeginInfo));

    std::vector<XrCompositionLayerBaseHeader*> layers;

    XrViewState viewState{ XR_TYPE_VIEW_STATE };
    auto viewCapacityInput = (uint32_t)m_views.size();
    uint32_t viewCountOutput;

    if (frameState.shouldRender == XR_TRUE) {

        XrViewLocateInfo viewLocateInfo{ XR_TYPE_VIEW_LOCATE_INFO };
        viewLocateInfo.viewConfigurationType = XR_VIEW_CONFIGURATION_TYPE_PRIMARY_STEREO;
        viewLocateInfo.displayTime = frameState.predictedDisplayTime;
        viewLocateInfo.space = m_appSpace;

        CHECK_XR_RESULT(xrLocateViews(m_session, &viewLocateInfo, &viewState, viewCapacityInput, &viewCountOutput, m_views.data()));

        CHECK(viewCountOutput == viewCapacityInput);
        CHECK(viewCountOutput == m_configViews.size());
        CHECK(viewCountOutput == m_swapchains.size());

        std::vector<XrCompositionLayerQuad> quadLayerViews(viewCountOutput);

        for (uint32_t viewCountIterator = 0; viewCountIterator < viewCountOutput; viewCountIterator++) {

            XrSwapchainImageAcquireInfo acquireInfo{ XR_TYPE_SWAPCHAIN_IMAGE_ACQUIRE_INFO };

            uint32_t swapchainImageIndex;
            CHECK_XR_RESULT(xrAcquireSwapchainImage(m_swapchains[viewCountIterator]->swapchain(), &acquireInfo, &swapchainImageIndex));

            XrSwapchainImageWaitInfo waitInfo{ XR_TYPE_SWAPCHAIN_IMAGE_WAIT_INFO };
            waitInfo.timeout = XR_INFINITE_DURATION;
            CHECK_XR_RESULT(xrWaitSwapchainImage(m_swapchains[viewCountIterator]->swapchain(), &waitInfo));

            auto width = (int32_t)m_swapchains[viewCountIterator]->createInfo()->width;
            auto height = (int32_t)m_swapchains[viewCountIterator]->createInfo()->height;

            if (m_drawFunction)
            {
                const XrSwapchainImageBaseHeader* const swapchainImage = m_swapchains[viewCountIterator]->image(swapchainImageIndex);
#ifdef XR_USE_PLATFORM_ANDROID
                const uint32_t colorTexture = reinterpret_cast<const XrSwapchainImageOpenGLESKHR*>(swapchainImage)->image;
#else
                const uint32_t colorTexture = reinterpret_cast<const XrSwapchainImageOpenGLKHR*>(swapchainImage)->image;
#endif
                (m_graphicsModuleObj.get()->*m_drawFunction)(colorTexture, 0, 0, width, height);
            }

            XrSwapchainImageReleaseInfo releaseInfo{ XR_TYPE_SWAPCHAIN_IMAGE_RELEASE_INFO };
            CHECK_XR_RESULT(xrReleaseSwapchainImage(m_swapchains[viewCountIterator]->swapchain(), &releaseInfo));

            XrCompositionLayerQuad& layer = quadLayerViews[viewCountIterator];
            layer.type = { XR_TYPE_COMPOSITION_LAYER_QUAD };
            layer.next = nullptr;
            layer.subImage.swapchain = m_swapchains[viewCountIterator]->swapchain();
            layer.subImage.imageRect.offset = { 0 ,0 };
            layer.subImage.imageRect.extent = { width ,height };
            layer.pose = XrPosef{ { 0, 0, 0, 1 }, { 0.0f, 0.0f, -1.0f } };
            layer.space = m_appSpace;
            layer.eyeVisibility = XR_EYE_VISIBILITY_BOTH;

            XrFovf fov = m_views[0].fov; // Using the first eye's FOV

            float fovHorizontal = fov.angleRight - fov.angleLeft;
            float fovVertical = fov.angleUp - fov.angleDown;

            float distance = 1.f; // Distance in meters
            float swidth = 2.0f * distance * tanf(fovHorizontal / 2.0f);
            float sheight = 2.0f * distance * tanf(fovVertical / 2.0f);

            layer.size.width = swidth;
            layer.size.height = sheight;

            layers.push_back(reinterpret_cast<XrCompositionLayerBaseHeader*>(&layer));

        }

        XrFrameEndInfo frameEndInfo{ XR_TYPE_FRAME_END_INFO };
        frameEndInfo.displayTime = frameState.predictedDisplayTime;
        frameEndInfo.environmentBlendMode = XR_ENVIRONMENT_BLEND_MODE_OPAQUE;
        frameEndInfo.layerCount = (uint32_t)layers.size();
        frameEndInfo.layers = layers.data();
        CHECK_XR_RESULT(xrEndFrame(m_session, &frameEndInfo));
    }
}

void OpenXRModule::setDrawFunction(std::shared_ptr<GraphicsModule> obj, DrawFunction drawFunc)
{
    m_graphicsModuleObj = obj;
    m_drawFunction = drawFunc;
}

void OpenXRModule::completeGraphicsBinding(ksGpuWindow* window)
{
#ifdef XR_USE_PLATFORM_ANDROID
    m_graphicsBinding.display = window->display;
    m_graphicsBinding.config = (EGLConfig)0;
    m_graphicsBinding.context = window->context.context;
#else
    m_graphicsBinding.hDC = window->context.hDC;
    m_graphicsBinding.hGLRC = window->context.hGLRC;
#endif
}

void OpenXRModule::cleanResources()
{
    for (const auto& swapchain : m_swapchains) {
        xrDestroySwapchain(swapchain->swapchain());
    }

    if (m_appSpace != XR_NULL_HANDLE) {
        xrDestroySpace(m_appSpace);
    }

    if (m_session != XR_NULL_HANDLE) {
        xrDestroySession(m_session);
    }

    if (m_instance != XR_NULL_HANDLE) {
        xrDestroyInstance(m_instance);
    }
}

std::vector<XrFovf> OpenXRModule::getFov() const {
    return std::vector<XrFovf> { m_views[0].fov , m_views[1].fov };
}
