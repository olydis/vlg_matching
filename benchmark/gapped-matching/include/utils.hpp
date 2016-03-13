#pragma once

#include <string>
#include <unistd.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <regex>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "timings.hpp"
#include "sdsl/int_vector.hpp"

typedef std::vector<uint64_t> string_type;

struct gapped_pattern {
    std::string raw_regexp;
    sdsl::int_vector<0> sdsl_regexp;
    std::vector<string_type> subpatterns;
    gapped_pattern(const std::string& p, bool string_patterns) : raw_regexp(p)
    {
        auto parse_subpattern = [&](std::string s) {
            if (string_patterns) {
                return string_type(s.begin(), s.end());
            } else {
                std::istringstream symbol_stream(s);
                string_type subpattern;
                uint64_t symbol;
                while (symbol_stream >> symbol)
                    subpattern.push_back(symbol);
                return subpattern;
            }
        };

        sdsl_regexp.resize(raw_regexp.size());
        for (size_t i=0; i<raw_regexp.size(); i++) {
            sdsl_regexp[i] = raw_regexp[i];
        }
        int64_t last_gap_end = 0;
        size_t gap_pos;
        while ((gap_pos = raw_regexp.find(".*", last_gap_end)) != std::string::npos) {
            auto subpattern = raw_regexp.substr(last_gap_end,gap_pos - last_gap_end);
            subpatterns.push_back(parse_subpattern(subpattern));
            last_gap_end = gap_pos + 2;
        }
        auto last_subpattern = raw_regexp.substr(last_gap_end);
        subpatterns.push_back(parse_subpattern(last_subpattern));
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
parse_pattern_file(std::string file_name, bool string_patterns)
{
    gm_timer tm("READ PATTERNS",true);
    std::vector<gapped_pattern> patterns;
    std::ifstream in(file_name);
    if (in) {
        std::string line;
        while (std::getline(in,line)) {
            try {
                gapped_pattern pat(line, string_patterns);
                patterns.push_back(pat);
            } catch (...) {
                LOG(ERROR) << "Could not parse pattern '" << line << "'. Skipped";
            }
        }
    } else {
        LOG(FATAL) << "Cannot open pattern file '" << file_name << "'";
    }
    LOG(INFO) << "read " << patterns.size() << " patterns";
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
