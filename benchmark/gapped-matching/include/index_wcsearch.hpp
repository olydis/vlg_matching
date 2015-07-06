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
            std::string s1;
            std::string s2;
            int min_gap;
            int max_gap;

            /* (1) parse pat (no error checking! parses: "<text>", "<text>.{<num>,<num>}<text>") */
            std::string regexp = pat.raw_regexp;
            int gap_decl_position = regexp.find(".{");
            if (gap_decl_position == -1)
            {
                s2 = s1 = regexp;
                min_gap = max_gap = 0;
            }
            else
            {
                s1 = regexp.substr(0, gap_decl_position);
                regexp = regexp.substr(gap_decl_position + 2);

                int gap_sepa_position = regexp.find(",");
                min_gap = std::stoi(regexp.substr(0, gap_sepa_position)) + s1.size();
                regexp = regexp.substr(gap_sepa_position + 1);

                int gap_end_position = regexp.find("}");
                max_gap = std::stoi(regexp.substr(0, gap_end_position)) + s1.size();
                s2 = regexp.substr(gap_end_position + 1);
            }

            /* (2) enum & copy the output */
            gapped_search_result res;
            for (auto hit : index.match2(s1, sdsl::incremental_wildcard_pattern(s2, min_gap, max_gap))) 
            {
                res.positions.push_back(hit.first - 1);
            }
            return res;
        }
};
