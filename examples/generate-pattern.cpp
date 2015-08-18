#include <sdsl/suffix_trees.hpp>
#include <iostream>

using namespace sdsl;
using namespace std;

int main(int argc, char** argv)
{
    if (argc <= 3) 
    {
        cout << "Usage: " << argv[0] << " file count length [charset]" << endl;
        cout << " - file:       The data to extract the patterns from" << endl;
        cout << " - count:      Number of patterns to extract" << endl;
        cout << " - length:     Minimum length of patterns to extract" << endl;
        cout << " - charset:" << endl;
        cout << "   - 0 = alpha-numerical only" << endl;
        cout << "   - 1 = anything" << endl;
        cout << "   - default = length >= 10" << endl;
        cout << endl;
        cout << " Uses a CST to extract most common patterns/phrases matching provided criteria." << endl;
        cout << endl;
        return 1;
    }

    // parse args
    string file(argv[1]);
    size_t count = atoi(argv[2]);
    size_t length = atoi(argv[3]);
    int charset = argc <= 4 ? (length >= 10 ? 1 : 0) : atoi(argv[4]);
    
    // CHARSET
    auto is_special_char = [](char c) { return isspace(c) || ispunct(c); };

    std::function<bool(char c)> charset_predicates[] = {
        [](char c) { return isalnum(c) || c == '_'; },
        [](char c) { return c != '\n' && c != '\r' && c != '\t' && c != '\v'; }
    };
    auto charset_predicate = charset_predicates[charset];
    auto charset_string_predicate = [&](string s) { return all_of(s.begin(), s.end(), charset_predicate); };
    
    // load (cached) CST
    string cst_file = file + ".sdsl";
    cst_sct3<> cst;
    if (!load_from_file(cst, cst_file))
    {
        construct(cst, file, 1);
        store_to_file(cst, cst_file);
    }

    // traverse
    size_t size = cst.size();
    stack<pair<size_t, string>> found;
    auto min_occ = [&]() { return found.empty() ? 0 : found.top().first; };
    for (auto it=cst.begin(); it!=cst.end(); ++it) 
    {
        if (it.visit() == 1) 
        {
            auto v = *it;
            auto depth = cst.depth(v);
            if (depth == 0)
                continue;
            
            if (cst.size(v) <= min_occ())
            {
                it.skip_subtree();
                continue;
            }
            
            string phrase = extract(cst, v).substr(0, length);
                        
            // check charset
            if (!charset_string_predicate(phrase))
            {
                it.skip_subtree();
                continue;
            }
            
            // check length
            if (depth >= length)
            {
                // check for word boundary (approximation)
                if (!is_special_char(cst.csa.L[cst.rb(v)]) 
                 && !is_special_char(cst.csa.L[cst.lb(v)]))
                {
                    it.skip_subtree();
                    continue;
                }

                size_t occ = cst.size(v);
                stack<pair<size_t, string>> helper;
                while (!found.empty() && found.top().first < occ)
                {
                    helper.push(found.top());
                    found.pop();
                }
                found.emplace(occ, phrase);
                while (found.size() < count && !helper.empty())
                {
                    found.push(helper.top());
                    helper.pop();
                }
                it.skip_subtree();
            }
        }
    }
    
    // fill
    while (found.size() < count)
        found.emplace(0, string(length, 'x'));    

    // output
    while (!found.empty())
    {
        cout << found.top().second << endl;
        found.pop();
    }
}
