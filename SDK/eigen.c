/* ------------------------------------------------------------------------- *\

    Eigenvalues and eigenvectors for 3x3 matrix
 
    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
    SPDX-License-Identifier: GPL-3.0-only
    All rights reserved.

    Author: Dmitry Shmunk

\* ------------------------------------------------------------------------- */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "eigen.h"

#define dot(v1,v2) ((v1)[0]*(v2)[0] + (v1)[1]*(v2)[1] + (v1)[2]*(v2)[2])

void qr_decompose(float M[3][3], float q[3][3], float r[3][3])
{
    memset(r,0,3*3*sizeof(float));

    float col[3], uv[3];
 
    for(int i = 0; i < 3; i++)
    {
        col[0] = M[0][i]; col[1] = M[1][i]; col[2] = M[2][i];
        for(int j = 0; j < i; j++)
        {
            uv[0] = q[0][j]; uv[1] = q[1][j]; uv[2] = q[2][j];
            float dp = dot(uv, col);
            col[0] -= uv[0]*dp; col[1] -= uv[1]*dp; col[2] -= uv[2]*dp;
            r[j][i] = dp;
        }
        float norm = sqrtf(dot(col, col));

        // ToDo: Fix possible Div0
        r[i][i] = norm;
        q[0][i] = col[0] / norm; q[1][i] = col[1] / norm; q[2][i] = col[2] / norm;
    }
}

void solve_qr(float M[3][3], float v[3], float sol[3])
{
    float q[3][3], r[3][3];
    qr_decompose(M, q, r);

    float rhs[3];
    rhs[0] = rhs[1] = rhs[2] = 0;
    for(int j = 0; j < 3; j++)
        for(int i = 0; i < 3; i++)
            rhs[i] += q[j][i] * v[j];

    for(int i = 2; i >= 0; i--)
    {
        float bk_sub = 0;
        for(int j = i+1; j <= 2; j++)
            bk_sub += sol[j] * r[i][j];
        sol[i] = (rhs[i] - bk_sub) / r[i][i];
    }
}

void backsolve(float M[3][3], float ev, float tol, int maxiter, float cur[3])
{
    cur[0] = cur[1] = cur[2] = 1;
    float prv[3];
    // Preturb eigenvalue: prevent from becoming singular
    float delta = ev + 0.000001f;
    
    float Md[3][3];
    memcpy (Md, M, 3*3*sizeof(float));
    Md[0][0] -= delta; Md[1][1] -= delta; Md[2][2] -= delta;

    for (int n=0; (n<maxiter) ; ++n)
    {
        prv[0] = cur[0]; prv[1] = cur[1]; prv[2] = cur[2];
        solve_qr(Md, prv, cur);
        // prevent oscillations due to sign change of the same vector
        if(cur[0] < 0)
        {
            cur[0] = -cur[0]; cur[1] = -cur[1]; cur[2] = -cur[2];
        }

        float norm = sqrtf(dot(cur, cur));
        cur[0] /= norm; cur[1] /= norm; cur[2] /= norm;

        if (fabs(cur[0] - prv[0]) <= tol &&
            fabs(cur[1] - prv[1]) <= tol &&
            fabs(cur[2] - prv[2]) <= tol ) break;
    }
}

void solve_eigenvectors(float M[3][3], float eigenvalues[3], float tol, int maxiter, float eigenvectors[3][3])
{
    for(int i = 0; i < 3; i++)
    {
        float evec[3];
        backsolve(M, eigenvalues[i], tol, maxiter, evec);
        eigenvectors[0][i] = evec[0];
        eigenvectors[1][i] = evec[1];
        eigenvectors[2][i] = evec[2];
    }
}

void solve_eigenvalues(float M[3][3], float tol, int maxiter, float evs[3])
{
    float X[3][3];
    memcpy (X, M, 3*3*sizeof(float));
    int upper_tri = 0;
    // QR algorithm iterations.
    for (int n=0; (n<maxiter) && !upper_tri; ++n)
    {
        float q[3][3], r[3][3];
        qr_decompose(X, q, r);

        memset(X, 0, 3*3*sizeof(float));
        for (int i = 0; i < 3; ++i)
            for (int k = 0; k < 3; ++k)
                for (int j = 0; j < 3; ++j)
                    X[i][j] += r[i][k] * q[k][j];

        upper_tri = 1;
        for (int i=0; (i<3) && upper_tri; ++i)
            for(int j=0; j<i; ++j)
                if (fabs(X[i][j]) > tol)
                    {
                        upper_tri = 0;
                        break;
                    }
    }
    evs[0] = X[0][0]; evs[1] = X[1][1]; evs[2] = X[2][2];
}

void eigen_solve(float M[3][3], float tol, int maxiter, float eigenvalues[3], float eigenvectors[3][3])
{
    solve_eigenvalues(M, tol, maxiter, eigenvalues);
    solve_eigenvectors(M, eigenvalues, tol, maxiter, eigenvectors);
}
