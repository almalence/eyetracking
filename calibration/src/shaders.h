/* ------------------------------------------------------------------------- *\

    Eye tracking calibration: Graphics handling, shaders

    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
    SPDX-License-Identifier: GPL-3.0-only
    All rights reserved.

    Authors: Bogdan Nudnenko
             Dmitry Shmunk

\* ------------------------------------------------------------------------- */

#pragma once

#ifdef XR_USE_PLATFORM_ANDROID
const char* glslVersion = "#version 320 es\nprecision highp float;\n";
#else
 const char* glslVersion = "#version 330 core\nprecision mediump float;\n";
#endif

// Shared vertex shader code
inline const char* vertexShaderBody = R"(
    layout(location = 0) in vec3 position;

    void main()
    {
        gl_Position = vec4(position.x, position.y, position.z, 1.0);
    }
)";

// Shared fragment shader code
inline const char* fragmentShaderBody = R"(
    out vec4 FragColor;

    uniform vec2 u_resolution;
    uniform vec2 u_circleCenter;
    uniform float u_radius;
    uniform int u_state;

    void main()
    {
        vec2 uv = gl_FragCoord.xy / u_resolution;
        float dist = distance(uv, u_circleCenter);

        if (u_state == 0)  // calibration stage
        {
            if (dist < u_radius)
                FragColor = vec4(0.6f, 0.3f, 0.1f, 1.0f);
            else
                discard;
        }
        else                // testing stage
        {
            #define A2N(x) ((x)/100.0f+0.5f)    // angles to normalized coordinates, assuming a typical 100 degree FOV of VR HMD
            vec2 testPoints[6];
            testPoints[0] = vec2(A2N( -5.0f), A2N(-8.0f));
            testPoints[1] = vec2(A2N(-10.0f), A2N( 0.0f));
            testPoints[2] = vec2(A2N( -5.0f), A2N( 8.0f));
            testPoints[3] = vec2(A2N(  5.0f), A2N(-8.0f));
            testPoints[4] = vec2(A2N( 10.0f), A2N( 0.0f));
            testPoints[5] = vec2(A2N(  5.0f), A2N( 8.0f));
            // distance to closest test point
            vec2 refCenter;
            float dist_ref = 100.0f;
            for(int i=0; i<6; ++i)
            {
                float d = distance(uv, testPoints[i]);
                if (d<dist_ref)
                {
                    refCenter = testPoints[i];
                    dist_ref = d;
                }
            }

            if ((dist < u_radius*5) || (dist_ref < u_radius))
            {
                vec4 clr = vec4(0.04f, 0.06f, 0.06f, 1.0f);

                if (dist_ref < u_radius)
                {
                    float circle2ref = distance(u_circleCenter, refCenter);
                    if (circle2ref < u_radius*5)
                        clr = clr + vec4(0.6f, 0.3f, 0.1f, 1.0f);
                    else
                        clr = vec4(0.0f, 0.0f, 0.0f, 1.0f);
                }

                FragColor = clr;
            }
            else
                discard;
        }
    }
)";

// Combine the version and shader body
inline std::string vertexShaderSourceStr = std::string(glslVersion) + vertexShaderBody;
inline std::string fragmentShaderSourceStr = std::string(glslVersion) + fragmentShaderBody;
inline const char* fragmentShaderSource = fragmentShaderSourceStr.c_str();
inline const char* vertexShaderSource = vertexShaderSourceStr.c_str();
