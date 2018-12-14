/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

//
// C++ Implementation: benchtimer
//
// Description:
//
//
// Author: Matthias Hess
// <matze@u-154-czaggh30.gpi.geowissenschaften.uni-tuebingen.de>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include <iostream>

#include <cstdlib>
#include <cstring>

#include "benchtimer.h"

BenchTimer::BenchTimer()
{
#ifndef _WIN32
    memset((void*)&start_val, 0, sizeof(struct timeval));
    memset((void*)&end_val, 0, sizeof(struct timeval));
#endif
    is_running = false;
}

BenchTimer::~BenchTimer() {}

double BenchTimer::precision_s()
{
    prec_s = 0.0;

    while (prec_s == 0.0)
    {
        start();
        stop();
        prec_s = time_s();
    }

    return prec_s;
}

double BenchTimer::precision_ms()
{
    return 1000.0 * precision_s();
}

void BenchTimer::start()
{
#ifndef _WIN32
    gettimeofday(&start_val, NULL);
#endif
    is_running = true;
}

void BenchTimer::stop()
{
#ifndef _WIN32
    gettimeofday(&end_val, NULL);
#endif
    is_running = false;
}

bool BenchTimer::running()
{
    return is_running;
}

double BenchTimer::time_s()
{
    double start_s, end_s;

#ifndef _WIN32
    start_s = (double)start_val.tv_sec + (double)start_val.tv_usec / 1000000.0;
    end_s = (double)end_val.tv_sec + (double)end_val.tv_usec / 1000000.0;
#else
    start_s = 0.0;
    end_s = 0.0;
#endif

    return end_s - start_s;
}

double BenchTimer::time_ms()
{
    double start_ms, end_ms;

#ifndef _WIN32
    start_ms =
        (double)start_val.tv_sec * 1000.0 + (double)start_val.tv_usec / 1000.0;
    end_ms = (double)end_val.tv_sec * 1000.0 + (double)end_val.tv_usec / 1000.0;
#else
    start_ms = 0.0;
    end_ms = 0.0;
#endif

    return end_ms - start_ms;
}
