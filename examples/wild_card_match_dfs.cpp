#include <sdsl/suffix_arrays.hpp>
#include <vector>
#include <regex>
#include <boost/regex.hpp>
#include <iostream>
#include "../include/sdsl/matching.hpp"

// BOOST

using namespace sdsl;
using namespace std;

using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;

//#define VERBOSE

// TODO: escape
// TODO: time regex construction (execute after timing on small string)
size_t ctor_total;

template <typename type_matching_index>
size_t match_ref(string& text, type_matching_index& index, string s1, string s2, size_t min_gap, size_t max_gap, std::function<void(size_t a, size_t b)> callback)
{
    (void)index;
    
    auto start = timer::now();
    std::ostringstream regex_stream;
    regex_stream << "(" << s1 << ".{" << min_gap << "," << max_gap << "})(" << s2 << ")";
    basic_regex<char> regex(regex_stream.str());
    auto stop = timer::now();

#ifdef VERBOSE
    cout << "regex: " << regex_stream.str() << endl;
#endif

    match_results<std::string::iterator> match;
    
    ctor_total += duration_cast<milliseconds>(stop-start).count();

    size_t result = 0;
    size_t offset = 0;
    while (regex_search(text.begin() + offset, text.end(), match, regex))
    {
        ++result;
        callback(offset + match.position(), offset + match.position() + match.length() - s2.size());
        offset += match.position() + match.length();
    }
    return result;
}

template <typename type_matching_index>
size_t match_ref_boost(string& text, type_matching_index& index, string s1, string s2, size_t min_gap, size_t max_gap, std::function<void(size_t a, size_t b)> callback)
{
    (void)index;
    
    auto start = timer::now();
    std::ostringstream regex_stream;
    regex_stream << "(" << s1 << ".{" << min_gap << "," << max_gap << "})(" << s2 << ")";
    boost::basic_regex<char> regex(regex_stream.str());
    auto stop = timer::now();

#ifdef VERBOSE
    cout << "regex: " << regex_stream.str() << endl;
#endif

    boost::match_results<std::string::iterator> match;
    
    ctor_total += duration_cast<milliseconds>(stop-start).count();

    size_t result = 0;
    size_t offset = 0;
    while (boost::regex_search(text.begin() + offset, text.end(), match, regex))
    {
        ++result;
        callback(offset + match.position(), offset + match.position() + match.length() - s2.size());
        offset += match.position() + match.length();
    }
    return result;
}

template <typename type_matching_index>
size_t match_dfs(string& text, type_matching_index& index, string s1, string s2, size_t min_gap, size_t max_gap, std::function<void(typename type_matching_index::size_type a, typename type_matching_index::size_type b)> callback)
{
    (void)text;
    
    size_t cnt = 0;
    for (auto res : index.match2(s1, incremental_wildcard_pattern(s2, 
        s1.size() + min_gap, 
        s1.size() + max_gap)))
    {
        callback(res.first, res.second);
        ++cnt;
    }
    return cnt;
}

template <typename type_matching_index>
size_t match_dfs_ex(type_matching_index& index, 
    string s, 
    incremental_wildcard_pattern p1)
{
    size_t cnt = 0;
    for (auto res : index.match2(s, p1))
    {
        cout << res.first << " " << res.second << endl;
        ++cnt;
    }
    return cnt;
}

int main(int argc, char* argv[])
{
    if ( argc < 2 ){
        cout << "Usage: ./" << argv[0] << " file_to_index" << endl;
        cout << "    Input pat1 pat2 min_gap max_gap" << endl;
        return 1;
    }

    cout << "loading as string" << endl;
    ifstream tmp(argv[1]);
    string text((istreambuf_iterator<char>(tmp)),
                 istreambuf_iterator<char>());

    cout << "loading as matching_index" << endl;

    //string index_file = string(argv[1]) + ".sdsl";

    cache_config cc(false,".","WILD_CARD_MATCH_TMP");
    matching_index<> index;
    construct(index, argv[1], cc, 1);

    // CALLBACKS
    using size_type = matching_index<>::size_type;
    std::function<void(size_type a, size_type b)> callback_cout = [&](size_type a, size_type b) 
    { 
        if (b - a > 100)
            cout << "\t" << a << " " << b << endl; 
        else
            cout << "\t" << a << " " << b << ": " << text.substr(a, b - a + 1) << endl;
    };
    std::function<void(size_type a, size_type b)> callback_dbg = [&](size_type a, size_type b) 
    {
        cout << "\t" << a << " " << b << endl; 
    };
    std::function<void(size_type a, size_type b)> callback_nop = [](size_type a, size_type b) { (void)a; (void)b; };

#ifdef VERBOSE
    auto callback_auto = callback_cout;
#else
    auto callback_auto = callback_nop;
#endif

    auto compare = [&](string s1, string s2, size_t min_gap, size_t max_gap) 
    { 
#ifdef VERBOSE
        cout << "Comparing with: s1=" << s1 << ", s2=" << s2 << ", min_gap=" << min_gap << ", max_gap=" << max_gap << endl;
#endif

#ifdef VERBOSE
        cout << "[own]" << endl;
#endif

        size_t matches = match_dfs(text, index, s1, s2, min_gap, max_gap, callback_auto);

#ifdef VERBOSE
        cout << "[reference]" << endl;
#endif

        size_t matches_ref = match_ref(text, index, s1, s2, min_gap, max_gap, callback_auto);
            
        bool success = matches == matches_ref;
        if (!success)
        {
            cout << "Error" << endl;
            cout << "MATCHES:     " << matches << endl;
            cout << "MATCHES_ref: " << matches_ref << endl;
        }
        else
        {
#ifdef VERBOSE
            cout << matches << " found" << endl;
#endif
        }
        
        return success;
    };

    bool test = false;
    bool bench = true;
    vector<pair<string, string>> tcs;
    tcs.emplace_back("include", "class");
    tcs.emplace_back("unsigned", "double");
    tcs.emplace_back("return", "print");
    tcs.emplace_back("insert", "remove");
    tcs.emplace_back("static", "const");
    tcs.emplace_back("cout", "endl");
    tcs.emplace_back("while", "switch");
    tcs.emplace_back("malloc", "free");
    tcs.emplace_back("open", "close");
    tcs.emplace_back("struct", "struct");
    if (test)
    {
        // TESTS
        int max = 100;
        for (unsigned int tci = 0; tci < tcs.size(); tci++)
        {
            auto s1 = tcs[tci].first;
            auto s2 = tcs[tci].second;
            cout << "TEST CASE: " << s1 << " " << s2 << endl;
            for (int i = 0; i <= max; i += 17)
            {
                if (!compare(s1, s2, 0, i)) goto err;
                if (!compare(s1, s2, i, i)) goto err;
                if (!compare(s1, s2, i, 2*i)) goto err;
                if (!compare(s1, s2, 10 * i, 10 * i + 10)) goto err;
                if (!compare(s1, s2, i, 10000)) goto err;
            }
        }
    }

    if (bench)
    {
        // BENCH
        cout << "===BENCH===" << endl;
        {
            size_t num_queries = 0;
            size_t total_time_our = 0;
            size_t total_time_regex = 0;
            size_t total_time_regex_boost = 0;

            auto exec_bench = [&](
                string version, 
                std::function<size_t(string& t, decltype(index)& i, string s1, string s2, size_t min_gap, size_t max_gap, std::function<void(size_type a, size_type b)> callback)> match,
                size_t min_gap, 
                size_t max_gap, 
                size_t& total_time)
            {
                ctor_total = 0;
                size_t found = 0;
                auto start = timer::now();
                for (unsigned int tci = 0; tci < tcs.size(); tci++)
                    found += match(text, index, tcs[tci].first, tcs[tci].second, min_gap, max_gap, callback_nop);
                auto stop = timer::now();
                cout << min_gap << ".." << max_gap << " ; \x1b[1m" << version << "\x1b[0m: \x1b[1;33m" << found << "\x1b[0m occurrences found in " <<  duration_cast<milliseconds>(stop-start).count() - ctor_total << "ms ( + construct: " << ctor_total << "ms)" << endl;
                total_time += duration_cast<milliseconds>(stop-start).count();
                return found;
            };
            
            auto exec_bench_comparison = [&](size_t min_gap, size_t max_gap)
            {
                auto f1 = exec_bench("OUR", match_dfs<decltype(index)>, 
                    min_gap, max_gap, total_time_our);
                auto f2 = exec_bench("RXB", match_ref_boost<decltype(index)>, 
                    min_gap, max_gap, total_time_regex_boost);
                auto f3 = exec_bench("RX ", match_ref<decltype(index)>, 
                    min_gap, max_gap, total_time_regex);

                if (f1 != f2 
                ||  f2 != f3
                ) cout << " \x1b[31m PANIC: #occ mismatch! \x1b[1m" << endl;
            };

            for (int gap_size = 0; gap_size <= 10 * 1000; gap_size += 1000)
            {
                num_queries += tcs.size();

                exec_bench_comparison(0, gap_size);
                exec_bench_comparison(gap_size, gap_size + 10);
            }


            cout << "#queries = " << num_queries << endl;
            cout << "total ms OUR   = " << total_time_our << endl;
            cout << "total ms REGEX = " << total_time_regex << endl;
        }
    }

    // PROMPT
    { 
        size_t min1_gap;
        size_t max1_gap;
        //size_t min2_gap;
        //size_t max2_gap;
        string s1, s2, s3;
        cout << "PROMPT: " << endl;
        while (cin >> s1 >> s2 >> min1_gap >> max1_gap && max1_gap >= min1_gap)
            compare(s1, s2, min1_gap, max1_gap);
    }
    err:;
}

