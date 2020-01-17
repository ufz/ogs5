/**
 * \copyright
 * Copyright (c) 2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <iostream>
#include <string>

#include "tokenbuf.h"

TokenBuf::TokenBuf(std::istream& stream, int internal_buf_size)
    : input_stream(stream),
      buf_max_size(internal_buf_size),
      buf_read_pos(0),
      buf_num_bytes(0),
      last_token_valid(false)
{
    buf = new char[buf_max_size];
}

TokenBuf::~TokenBuf()
{
    delete[] buf;
}

bool TokenBuf::done()
{
    if ((buf_read_pos == buf_num_bytes && !last_token_valid &&
         input_stream.eof()) ||
        input_stream.bad())
        return true;
    else
        return false;
}

void TokenBuf::readBuf()
{
    buf_read_pos = 0;
    input_stream.read(buf, buf_max_size);
    buf_num_bytes = input_stream.gcount();
}

void TokenBuf::getNextToken(bool ignore_line_break)
{
    int i;
    bool got_it = false;
    int status = 0;

    last_token = "";

    // Skip non token characters

    while (status == 0)
    {
        for (i = buf_read_pos; i < buf_num_bytes; i++, buf_read_pos++)
            if ((buf[i] != ' ' && buf[i] != '\t'))
            {
                if (ignore_line_break)
                {
                    if (buf[i] != '\r' && buf[i] != '\n')
                    {
                        status = 1;
                        break;
                    }
                }
                else
                {
                    status = 1;
                    break;
                }
            }

        if (status == 0)
        {
            if (input_stream.eof())
                status = 4;
            readBuf();
            if (input_stream.bad())
                status = 255;
        }
    }

    // And copy token characters to 'last_token'

    if (!ignore_line_break &&
        (buf[buf_read_pos] == '\r' || buf[buf_read_pos] == '\n'))
    {
        last_token += buf[buf_read_pos];
        buf_read_pos++;
        status = 2;
    }

    while (status == 1)
    {
        for (i = buf_read_pos; i < buf_num_bytes; i++, buf_read_pos++)
        {
            if (buf[i] != ' ' && buf[i] != '\t' && buf[i] != '\r' &&
                buf[i] != '\n')
            {
                printf("%d ", buf[i]);
                last_token += buf[i];
            }
            else
            {
                status = 2;
                break;
            }
        }

        if (status == 1)
        {
            if (input_stream.eof())
                status = 4;
            readBuf();
            if (input_stream.bad())
                status = 255;
        }
    }

    if (status == 2)
        last_token_valid = true;
}

TokenBuf::Types TokenBuf::type()
{
    Types rv;
    int i;
    int status;

    return NONE;
}

std::string TokenBuf::peek(bool ignore_line_break)
{
    int i;

    if (!last_token_valid)
    {
        getNextToken(ignore_line_break);
        last_token_valid = true;
    }

    return last_token;
}

std::string TokenBuf::consume(bool ignore_line_break)
{
    std::string rv;

    if (last_token_valid)
        rv = last_token;
    else
        rv = peek(ignore_line_break);

    last_token_valid = false;

    return rv;
}

void TokenBuf::getline(char* linebuf, int max_num)
{
    int i, pos;
    bool reading = true;

    pos = 0;

    if (last_token_valid)
    {
        if (last_token[0] != '\n' && last_token[0] != '\r')
        {
            for (i = 0; i < last_token.size() && i < max_num - 1; i++)
                linebuf[i] = last_token[i];
            pos = i;
        }
        last_token_valid = false;
    }

    // If there is a rest of input buffer handle that first...

    do
    {
        for (i = buf_read_pos; i < buf_num_bytes && pos < max_num - 1; i++)
        {
            if (buf[i] != '\r' && buf[i] != '\n')
            {
                linebuf[pos] = buf[i];
                pos++;
            }
            else
            {
                reading = false;
                break;
            }
        }

        buf_read_pos = i + 1;

        if (reading)
        {
            readBuf();
            reading = !done();
        }
    } while (reading);

    linebuf[pos] = '\0';
}

bool TokenBuf::get_non_empty_line(char* linebuf, int max_num)
{
    bool valid = true;

    if (!done())
    {
        do
        {
            getline(linebuf, max_num);
        } while (!done() && linebuf[0] == '\0');
        valid = (!input_stream.bad()) && linebuf[0] != '\0';
    }
    else
        valid = false;

    return valid;
}
