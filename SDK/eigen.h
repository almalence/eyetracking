/* ------------------------------------------------------------------------- *\

    Eigenvalues and eigenvectors for 3x3 matrix
 
    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
    SPDX-License-Identifier: GPL-3.0-only
    All rights reserved.

    Author: Dmitry Shmunk

\* ------------------------------------------------------------------------- */

#pragma once

// Eigenvalues and eigenvectors of a matrix
// Eigenvalues computed with QR algorithm
// Eigenvectors computed with inverse iteration
void eigen_solve(float M[3][3], float tol, int maxiter, float eigenvalues[3], float eigenvectors[3][3]);
