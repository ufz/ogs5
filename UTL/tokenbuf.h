/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// -*- C++ -*-

#ifndef _TOKENBUF_H
#define _TOKENBUF_H

#include <iostream>
#include <string>

class TokenBuf
{
public:
    enum Types
    {
        NONE,
        INTEGER,
        REAL,
        CHARACTER,
        STRING
    };

    explicit TokenBuf(std::istream& stream, int internal_buf_size = 1024);
    ~TokenBuf();

    bool done();

    Types type();
    std::string peek(bool ignore_line_break = true);
    std::string consume(bool ignore_line_break = true);

    void getline(char* linebuf, int max_num);
    bool get_non_empty_line(char* linbuf, int max_num);

protected:
    void readBuf();
    void getNextToken(bool ignore_line_break);

    std::istream& input_stream;
    std::string last_token;
    bool last_token_valid;

    int buf_max_size;
    int buf_read_pos;
    int buf_num_bytes;
    char* buf;
};

#endif
