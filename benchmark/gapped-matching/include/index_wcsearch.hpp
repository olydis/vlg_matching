#pragma once

#include "utils.hpp"
#include "collection.hpp"
#include "sdsl/suffix_arrays.hpp"
#include "sdsl/matching.hpp"

//size_t sdsl::PERFCTR_NUM_PROCESSED_WT_NODES;

// THIS IS JUST A PROXY until a consistent API for indices is figured out
// (or is the "collection"-approach final?)
class index_wcsearch
{
    private:
        typedef sdsl::matching_index<sdsl::csa_wt<sdsl::wt_int<>>, sdsl::wt_int<>, sdsl::rrr_vector<>> index_type;
        index_type index;

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

        std::string info(const gapped_pattern& pat) const 
        {
            // output SA-ranges (gives a good estimation about potential matches)
            index_type::size_type total_range = 0, sp = 0, ep = 0;

            std::string s1;
            std::string s2;

            for (size_t i = 0; i < pat.subpatterns.size(); ++i)
                total_range += backward_search(index.csa, 0, index.csa.size()-1, pat.subpatterns[i].begin(), pat.subpatterns[i].end(), sp, ep);

            return std::to_string(total_range);
        }
        void prepare(const gapped_pattern& pat) { (void)pat; }

        gapped_search_result
        search(const gapped_pattern& pat) const
        {
            gapped_search_result res;
            string_type s1;
            string_type s2;
            size_type min_gap;
            size_type max_gap;

            if (pat.subpatterns.size() == 1) {
                s1 = pat.subpatterns[0];
                s2 = pat.subpatterns[0];
                min_gap = 0;
                max_gap = 0;
            } else {
                s1 = pat.subpatterns[0];
                s2 = pat.subpatterns[1];
                min_gap = s1.size() + pat.gaps[0].first;
                max_gap = s1.size() + pat.gaps[0].second;
                // add "s1.size()" because "match2" currently requires word-beginning-relative gaps
                // (this is an important concept, as it allows single-term matching by setting min/max_gap=0)
            }

            sdsl::PERFCTR_NUM_PROCESSED_WT_NODES = 1;
            for (auto hit : index.match(s1, s2, min_gap, max_gap)) {
                res.positions.push_back(hit.first);
            }
            cout << "Processed nodes: " << sdsl::PERFCTR_NUM_PROCESSED_WT_NODES << endl;
            return res;
        }
};
