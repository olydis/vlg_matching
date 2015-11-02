#pragma once

#include "utils.hpp"
#include "collection.hpp"
#include "sdsl/suffix_arrays.hpp"
#include "sdsl/matching.hpp"

// implements double binary search
class index_sasearch
{
    private:
        typedef sdsl::matching_index<sdsl::wt_int<>, sdsl::rrr_vector<>> index_type;
        index_type index;
        typedef std::pair<index_type::size_type, index_type::size_type> range_type;
        typedef std::pair<range_type, range_type> double_range_type;

    public:
        typedef sdsl::int_vector<0>::size_type size_type;
        std::string name() const
        {
            std::string index_name = IDXNAME;
            return "SASEARCH-"+index_name;
        }

    public:
        index_sasearch() { }
        index_sasearch(collection& col)
        {
            sdsl::cache_config cc(false,".","SASEARCH_TMP");
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

        void swap(index_sasearch& ir)
        {
            if (this != &ir) {
                index.swap(ir.index);
            }
        }

        std::string info(const gapped_pattern& pat) const { (void)pat; return ""; }
        void prepare(const gapped_pattern& pat) { (void)pat; }

        gapped_search_result
        search(const gapped_pattern& pat) const
        {
            gapped_search_result res;
            string_type s1;
            string_type s2;
            size_type min_gap = 0;
            size_type max_gap = 0;

            if (pat.subpatterns.size() == 1) {
                s1 = pat.subpatterns[0];
                s2 = pat.subpatterns[0];
            } else {
                s1 = pat.subpatterns[0];
                s2 = pat.subpatterns[1];
                min_gap = pat.gaps[0].first + s1.size();
                max_gap = pat.gaps[0].second + s1.size();
            }

            // get ranges
            size_type sp1, ep1;
            size_type sp2, ep2;
            vector<size_type> range_a(forward_search(index.text.begin(), index.text.end(), index.wt, 0, index.wt.size()-1, s1.begin(), s1.end(), sp1, ep1));
            vector<size_type> range_b(forward_search(index.text.begin(), index.text.end(), index.wt, 0, index.wt.size()-1, s2.begin(), s2.end(), sp2, ep2));
            std::copy(index.wt.begin() + sp1, index.wt.begin() + ep1 + 1, range_a.begin());
            std::copy(index.wt.begin() + sp2, index.wt.begin() + ep2 + 1, range_b.begin());
            std::sort(range_a.begin(), range_a.end());
            std::sort(range_b.begin(), range_b.end());

            // linear search
            auto a_it = range_a.begin();
            auto b_it = range_b.begin();

            while (a_it != range_a.end()) {
                auto a_pos = *a_it;

                // enforcing min_gap
                bool b_valid;
                while ((b_valid = (b_it != range_b.end())) && a_pos + min_gap > *b_it)
                    ++b_it;
                if (!b_valid)
                    break;

                // check whether within max_gap
                auto b_pos = *b_it;
                if (a_pos + max_gap < b_pos) {
                    ++a_it;
                    continue;
                }

                // push greedy beyond max_gap
                ++b_it;
                while (b_it != range_b.end()) {
                    auto b_pos2 = *b_it;
                    if (a_pos + max_gap >= b_pos2)
                        b_pos = b_pos2;
                    else
                        break;
                    ++b_it;
                }

                res.positions.push_back(a_pos);

                // pull a beyond previous b (non-overlapping)
                b_pos += s2.size();
                while (a_it != range_a.end() && *a_it < b_pos)
                    ++a_it;
            }

            return res;
        }
};
