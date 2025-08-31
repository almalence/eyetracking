/* ------------------------------------------------------------------------- *\

    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
    SPDX-License-Identifier: GPL-3.0-only
    All rights reserved.

    Author: Bogdan Nudnenko

\* ------------------------------------------------------------------------- */

#pragma once

#ifdef XR_USE_PLATFORM_WIN32
#include <windows.h>
#elif XR_USE_PLATFORM_ANDROID
#include <jni.h>
#endif
#ifdef XR_USE_GRAPHICS_API_OPENGL_ES
#include <EGL/egl.h>
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <locale>
#include <memory>
#include <stddef.h>
#include <stdarg.h>
#include <vector>
#include <map>
#include <list>
#include <algorithm>
#include <chrono>
#include <thread>
#include <mutex>

inline std::string Fmt(const char* fmt, ...) {
    va_list vl;
    va_start(vl, fmt);
    int size = std::vsnprintf(nullptr, 0, fmt, vl);
    va_end(vl);

    if (size != -1) {
        std::unique_ptr<char[]> buffer(new char[size + 1]);

        va_start(vl, fmt);
        size = std::vsnprintf(buffer.get(), size + 1, fmt, vl);
        va_end(vl);
        if (size != -1) {
            return std::string(buffer.get(), size);
        }
    }

    throw std::runtime_error("Unexpected vsnprintf failure");
}


[[noreturn]] inline void Throw(std::string failureMessage, const char* originator = nullptr, const char* sourceLocation = nullptr) {
    if (originator != nullptr) {
        failureMessage += Fmt("\n    Origin: %s", originator);
    }
    if (sourceLocation != nullptr) {
        failureMessage += Fmt("\n    Source: %s", sourceLocation);
    }

    throw std::logic_error(failureMessage);
}

#define ALM_STRINGIFY(x) #x
#define TOSTRING(x) ALM_STRINGIFY(x)
#define FILE_AND_LINE __FILE__ ":" TOSTRING(__LINE__)

#define THROW(msg) Throw(msg, nullptr, FILE_AND_LINE);

#define CHECK(exp)                                      \
    {                                                   \
        if (!(exp)) {                                   \
            Throw("Check failed", #exp, FILE_AND_LINE); \
        }                                               \
    }

#define CHECK_MSG(exp, msg)                  \
    {                                        \
        if (!(exp)) {                        \
            Throw(msg, #exp, FILE_AND_LINE); \
        }                                    \
    }