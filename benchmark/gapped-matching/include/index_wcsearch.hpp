#pragma once

#include "utils.hpp"
#include "collection.hpp"
#include "sdsl/suffix_arrays.hpp"
#include "sdsl/matching.hpp"

// THIS IS JUST A PROXY until a consistent API for indices is figured out
// (or is the "collection"-approach final?)
class index_wcsearch
{
    private:
        sdsl::matching_index<> index;

    public:
        typedef sdsl::int_vector<0>::size_type size_type;
        std::string name() const
        {
            std::string index_name = IDXNAME;
            return "WCSEARCH-"+index_name;
        }

    public:
        index_wcsearch() { }
        index_wcsearch(collection& col)
        {
            sdsl::cache_config cc(false,".","WCSEARCH_TMP");
            sdsl::construct(index, col.file_map[consts::KEY_TEXT], cc, 0);
        }

        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v=NULL, std::string name="")const
        {
            return index.serialize(out, v, name);
        }

        void load(std::istream& in)
        {
            index.load(in);
        }

        void swap(index_wcsearch& ir)
        {
            if (this != &ir) {
                index.swap(ir.index);
            }
        }

        gapped_search_result
        search(const gapped_pattern& pat) const
        {
            gapped_search_result res;
            std::string s1;
            std::string s2;
            size_type min_gap = 0;
            size_type max_gap = 0;

            if (pat.subpatterns.size() == 1) {
                s1 = pat.subpatterns[0];
                s2 = pat.subpatterns[0];
            } else {
                s1 = pat.subpatterns[0];
                s2 = pat.subpatterns[1];
                min_gap = s1.size() + pat.gaps[0].first;
                max_gap = s1.size() + pat.gaps[0].second;
                // add "s1.size()" because "match2" currently requires word-beginning-relative gaps (which implicitly allows single-term matching by setting min/max_gap=0)
            }

            for (auto hit : index.match2(s1, sdsl::incremental_wildcard_pattern(s2, min_gap, max_gap))) {
                res.positions.push_back(hit.first);
            }
            return res;
        }
};
