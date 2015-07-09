#include <sdsl/suffix_trees.hpp>
#include <iostream>

using namespace sdsl;
using namespace std;

int main(int argc, char** argv)
{
    if (argc <= 3) 
    {
        cout << "Usage: " << argv[0] << " file count length [mode [charset]]" << endl;
        cout << " - file:       The data to extract the patterns from" << endl;
        cout << " - count:      Number of patterns to extract" << endl;
        cout << " - length:     Minimum length of patterns to extract" << endl;
        cout << " - mode:" << endl;
        cout << "   - 0 = common (default)" << endl;
        cout << "   - 1 = rare" << endl;
        cout << "   - 2 = lexicographically" << endl;
        cout << " - charset:" << endl;
        cout << "   - 0 = alpha-numerical only (default)" << endl;
        cout << "   - 1 = anything" << endl;
        cout << endl;
        cout << " Uses a CST to extract patterns/phrases matching provided criteria." << endl;
        cout << endl;
        return 1;
    }

    // parse args
    string file(argv[1]);
    size_t count = atoi(argv[2]);
    size_t length = atoi(argv[3]);
    int mode = argc <= 4 ? 0 : atoi(argv[4]);
    int charset = argc <= 5 ? 0 : atoi(argv[5]);
    
    // CHARSET
    auto is_special_char = [](char c) { return isspace(c) || ispunct(c); };

    std::function<bool(char c)> charset_predicates[] = {
        [](char c) { return isalnum(c) || c == '_'; },
        [](char c) { (void)c; return true; }
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
    vector<pair<size_t, string>> found;
    for (auto it=cst.begin(); it!=cst.end(); ++it) 
    {
        if (it.visit() == 1) 
        {
            auto v = *it;
            if (cst.depth(v) == 0)
                continue;
            
            string phrase = extract(cst, v).substr(0, length);
            
            // check charset
            if (!charset_string_predicate(phrase))
            {
                it.skip_subtree();
                continue;
            }
            
            // check length
            if (phrase.size() == length)
            {
                // check for word boundary (approximation)
                if (!is_special_char(cst.csa.L[cst.rb(v)]) 
                 && !is_special_char(cst.csa.L[cst.lb(v)]))
                {
                    it.skip_subtree();
                    continue;
                }

                size_t occ = 1 + cst.rb(v) - cst.lb(v);
                found.emplace_back(occ, phrase);
                it.skip_subtree();
            }
        }
    }
    
    // sort
    switch (mode)
    {
        case 0: // common
            sort(found.rbegin(), found.rend());
            break;
        case 1: // rare
            sort(found.begin(), found.end());
            break;
        case 2: // lexicographically
            break;
    }
    
    // head
    if (found.size() > count)
        found.resize(count);
    
    // output
    for (auto result : found)
        cout << result.second << endl;
}
