#pragma once

#include "collection.hpp"
#include "utils.hpp"
#include "sdsl/matching.hpp"
#include "sdsl/suffix_arrays.hpp"

#include <boost/regex.hpp>

template<class type_index>
class wild_card_match_iterator3 : public std::iterator<std::forward_iterator_tag, array<typename type_index::size_type, 3>>
{
    private:
        typedef typename type_index::node_type node_type;
        typedef typename type_index::size_type size_type;
        typedef typename type_index::wt_type   wt_type;
        typedef array<typename type_index::size_type, 3>     result_type;

        // (lex_range, node)
        vector<wavelet_tree_range_walker<type_index>> lex_ranges;

        size_t a = 0;
        size_t b = 0;
        size_t c = 0;

        size_t min_gap;
        size_t max_gap;

        result_type current;

        bool next()
        {
            // find next independent match
            while (valid()) {
                const auto& top0 = lex_ranges[0].current_node();
                const auto& top1 = lex_ranges[1].current_node();
                const auto& top2 = lex_ranges[2].current_node();
                if (top1->range_end_doc < top2->range_begin_doc)
                    lex_ranges[1].next();
                else if (top0->range_end_doc < top1->range_begin_doc)
                    lex_ranges[0].next();

                else if (top1->range_end + p2.max_gap < top2->range_begin)
                    lex_ranges[1].next();
                else if (top1->range_begin + p2.min_gap > top2->range_end)
                    lex_ranges[2].next();

                else if (top0->range_end + p1.max_gap < top1->range_begin)
                    lex_ranges[0].next();
                else if (top0->range_begin + p1.min_gap > top1->range_end)
                    lex_ranges[1].next();

                else if (top0->is_leaf && top1->is_leaf && top2->is_leaf) {
                    a = top0->range_begin;
                    auto doc = top0->range_end_doc;

                    if (c != 0 && a <= c) {
                        lex_ranges[0].next();
                        continue;
                    }

                    b = top1->range_begin;
                    c = top2->range_begin;

                    lex_ranges[1].next();
                    lex_ranges[2].next();

                    // push c forward
                    while (lex_ranges[2].has_more()
                           && b + p2.max_gap >= lex_ranges[2].current_node()->range_begin
                           && doc == lex_ranges[2].current_node()->range_begin_doc) {
                        auto leaf = lex_ranges[2].retrieve_leaf_and_traverse();
                        if (leaf != nullptr)
                            c = leaf->range_begin;
                    }

                    auto state1 = lex_ranges[1].save_state();
                    auto state2 = lex_ranges[2].save_state();

                    // push b forward
                    while (lex_ranges[1].has_more()
                           && a + p1.max_gap >= lex_ranges[1].current_node()->range_begin
                           && doc == lex_ranges[1].current_node()->range_begin_doc) {
                        auto leaf = lex_ranges[1].retrieve_leaf_and_traverse();
                        if (leaf != nullptr) {
                            auto b_temp = leaf->range_begin;

                            if (b_temp + p2.min_gap <= c)
                                b = b_temp;

                            // push c forward
                            while (lex_ranges[2].has_more()
                                   && b_temp + p2.max_gap >= lex_ranges[2].current_node()->range_begin
                                   && doc == lex_ranges[2].current_node()->range_begin_doc) {
                                auto leaf = lex_ranges[2].retrieve_leaf_and_traverse();
                                if (leaf != nullptr) {
                                    b = b_temp;
                                    c = leaf->range_begin;
                                }
                            }
                        }
                    }

                    lex_ranges[1].restore_state(state1);
                    lex_ranges[2].restore_state(state2);

                    // pull a forward
                    while (valid() &&
                           lex_ranges[0].current_node()->range_end <= c)
                        lex_ranges[0].next();


                    current = { a, b, c };
                    return true;
                } else
                    lex_ranges[top1->range_size() >= top0->range_size() ? (top2->range_size() >= top1->range_size() ? 2 : 1) : (top2->range_size() >= top0->range_size() ? 2 : 0)].split();
            }

            current = { 0, 0, 0 };
            return false;
        }

    public:
        wild_card_match_iterator3()
        {
            current = { 0, 0, 0 };
        }
        wild_card_match_iterator3(const type_index& index,
                                  string s1,
                                  string s2,
                                  string s3,
                                  size_t min_gap,
                                  size_t max_gap)
            : min_gap(min_gap), max_gap(max_gap)
        {
            auto root_node = make_shared<node_cache<type_index>>(index.wt.root(), index, nullptr, nullptr);
            size_type sp = 1, ep = 0;

            backward_search(index.csa, 0, index.csa.size()-1, s1.begin(), s1.end(), sp, ep);
            lex_ranges.emplace_back(index, range_type(sp, ep),root_node);

            backward_search(index.csa, 0, index.csa.size()-1, s2.begin(), s2.end(), sp, ep);
            lex_ranges.emplace_back(index, range_type(sp, ep),root_node);

            backward_search(index.csa, 0, index.csa.size()-1, s3.begin(), s3.end(), sp, ep);
            lex_ranges.emplace_back(index, range_type(sp, ep),root_node);

            next();
        }

        bool valid() const
        {
            return lex_ranges[0].has_more()
                   && lex_ranges[1].has_more()
                   && lex_ranges[2].has_more();
        }

        result_type operator*() const
        {
            return current;
        }
        result_type* operator->()
        {
            return &current;
        }

        wild_card_match_iterator3& operator++()
        {
            next();
            return *this;
        }

        friend bool operator==(
            const wild_card_match_iterator& a,
            const wild_card_match_iterator& b)
        {
            return a.current == b.current;
        }

        friend bool operator!=(
            const wild_card_match_iterator& a,
            const wild_card_match_iterator& b)
        {
            return !(a == b);
        }
};


class index_wcsearch3
{
    private:
        typedef sdsl::matching_index<sdsl::wt_int<>, sdsl::rrr_vector<>> index_type;
        index_type index;
        string text = "";

    public:
        typedef sdsl::int_vector<0>::size_type size_type;
        std::string name() const
        {
            std::string index_name = IDXNAME;
            return "WCSEARCH-"+index_name;
        }

    public:
        index_wcsearch3() { }
        index_wcsearch3(collection& col)
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

        void swap(index_wcsearch3& ir)
        {
            if (this != &ir) {
                index.swap(ir.index);
            }
        }

        std::string info(const gapped_pattern& pat) const
        {
            // output SA-ranges (gives a good estimation about potential matches)
            index_type::size_type total_range = 0, sp = 0, ep = 0;

            for (size_t i = 0; i < pat.subpatterns.size(); ++i)
                total_range += forward_search(index.text.begin(), index.text.end(), index.wt, 0, index.wt.size()-1, pat.subpatterns[i].begin(), pat.subpatterns[i].end(), sp, ep);

            return std::to_string(total_range);
        }
        void prepare(const gapped_pattern& pat)
        {
            (void)pat;
            if (index.text.width() <= 8)
                text = string(index.text.begin(), index.text.end());
        }

        gapped_search_result
        search(const gapped_pattern& pat) const
        {
            gapped_search_result res;
            string_type s1;
            string_type s2;
            string_type s3;
            size_type min_gap;
            size_type max_gap;

            if (pat.subpatterns.size() < 3) {
                s1 = pat.subpatterns[0];
                s2 = pat.subpatterns[0];
                s3 = pat.subpatterns[0];
                min_gap = 0;
                max_gap = 0;
            } else {
                s1 = pat.subpatterns[0];
                s2 = pat.subpatterns[1];
                s3 = pat.subpatterns[2];
                min_gap = s1.size() + pat.gaps[0].first;
                max_gap = s1.size() + pat.gaps[0].second;
                // add "s1.size()" because "match2" currently requires word-beginning-relative gaps
                // (this is an important concept, as it allows single-term matching by setting min/max_gap=0)
            }

            // smart scan
            auto container = matching_container<iterator>(
                                 wild_card_match_iterator3<index_type>(*this, s1, s2, s3, min_gap, max_gap),
                                 wild_card_match_iterator3<index_type>());
            for (auto hit : container) {
                res.positions.push_back(hit.first);
            }
            return res;
        }
};
