#pragma once

#include <string>
#include <unistd.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "timings.hpp"
#include "sdsl/int_vector.hpp"

struct gapped_pattern {
    std::string raw_regexp;
    sdsl::int_vector<0> sdsl_regexp;
    gapped_pattern(const std::string& p) : raw_regexp(p)
    {
        sdsl_regexp.resize(raw_regexp.size());
        for (size_t i=0; i<raw_regexp.size(); i++) {
            sdsl_regexp[i] = raw_regexp[i];
        }
    };
};

struct gapped_search_result {
    std::vector<uint64_t> positions;
    gapped_search_result() = default;
    gapped_search_result(size_t n)
    {
        positions.resize(n);
    }
};

namespace utils
{

std::vector<gapped_pattern>
parse_pattern_file(std::string file_name)
{
    std::vector<gapped_pattern> patterns;
    std::ifstream in(file_name);
    if (in) {
        std::string line;
        while (std::getline(in,line)) {
            patterns.emplace_back(gapped_pattern(line));
        }
    } else {
        LOG(FATAL) << "Cannot open pattern file '" << file_name << "'";
    }
    return patterns;
}

bool directory_exists(std::string dir)
{
    struct stat sb;
    const char* pathname = dir.c_str();
    if (stat(pathname, &sb) == 0 && (S_IFDIR & sb.st_mode)) {
        return true;
    }
    return false;
}

bool file_exists(std::string file_name)
{
    std::ifstream in(file_name);
    if (in) {
        in.close();
        return true;
    }
    return false;
}

void create_directory(std::string dir)
{
    if (!directory_exists(dir)) {
        if (mkdir(dir.c_str(), 0777) == -1) {
            LOG(FATAL) << "could not create directory";
        }
    }
}


}