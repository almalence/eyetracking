/* ------------------------------------------------------------------------- *\

    Ellipse fitting functions
 
    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
    SPDX-License-Identifier: GPL-3.0-only
    All rights reserved.

    Author: Dmitry Shmunk

\* ------------------------------------------------------------------------- */

#pragma once

int EllipseDirectFit(float *xy, int xylen, int skip_st, int skip_en, float A[6]);

int EllipseGeomFromAlgebraic(float A[6], float *cx, float *cy, float *semx, float *semy, float *phi);
