#pragma once

#include "collection.hpp"
#include "utils.hpp"
#include "sdsl/matching.hpp"
#include "sdsl/suffix_arrays.hpp"

#include <boost/regex.hpp>


template<class type_index>
class wavelet_tree_range_walker
{
    private:
        typedef shared_ptr<sdsl::node_cache<type_index>> node_type;
        const type_index& index;
        stack<pair<sdsl::range_type,node_type>> dfs_stack;

    public:
        typedef decltype(dfs_stack) state_type;
        wavelet_tree_range_walker(const type_index& index, sdsl::range_type initial_range, node_type root_node)
            : index(index)
        {
            dfs_stack.emplace(initial_range, root_node);
        }

        void restore_state(const state_type& state)
        {
            dfs_stack = state;
        }
        state_type save_state()
        {
            return dfs_stack;
        }

        bool has_more() const
        {
            return !dfs_stack.empty();
        }

        node_type current_node() const
        {
            return dfs_stack.top().second;
        }

        void skip_subtree()
        {
            dfs_stack.pop();
        }

        void expand()
        {
            auto top = dfs_stack.top(); dfs_stack.pop();
            auto& node = top.second;
            node->ensure_children();
            auto exp_range = index.wt.expand(node->node, top.first);
            if (!sdsl::empty(exp_range[1]))
                dfs_stack.emplace(exp_range[1], node->children.second);
            if (!sdsl::empty(exp_range[0]))
                dfs_stack.emplace(exp_range[0], node->children.first);
        }

        // performs one DFS step, retrieving a leaf or nullptr
        node_type retrieve_leaf_and_traverse()
        {
            node_type node = current_node();
            if (node->is_leaf) {
                skip_subtree();
                return node;
            }
            expand();
            return nullptr;
        }

        node_type next_leaf()
        {
            if (has_more() && current_node()->is_leaf)
                skip_subtree();
            while (has_more() && !current_node()->is_leaf)
                expand();
            return has_more() ? current_node() : nullptr;
        }
};

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

        size_t size3;

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
                if (top1->range_end + max_gap < top2->range_begin) {
                    lex_ranges[1].skip_subtree();
                } else if (top1->range_begin + min_gap > top2->range_end) {
                    lex_ranges[2].skip_subtree();
                }

                else if (top0->range_end + max_gap < top1->range_begin) {
                    lex_ranges[0].skip_subtree();
                } else if (top0->range_begin + min_gap > top1->range_end) {
                    lex_ranges[1].skip_subtree();
                }

                else if (top0->is_leaf && top1->is_leaf && top2->is_leaf) {
                    a = top0->range_begin;
                    b = top1->range_begin;
                    c = top2->range_begin;



                    // push c forward
                    while (lex_ranges[2].next_leaf() != nullptr
                           && b + max_gap >= lex_ranges[2].current_node()->range_begin) {
                        c = lex_ranges[2].current_node()->range_begin;
                    }

                    //----------
                    while (lex_ranges[1].next_leaf() != nullptr) {
                        auto b_pos2 = lex_ranges[1].current_node();
                        if (a + max_gap >= b_pos2->range_begin) {
                            b = b_pos2->range_begin;
                            // situation: found greedier ab ==> check c

                            // enforcing min_gap bc
                            bool c_valid;
                            if (b + min_gap > lex_ranges[2].current_node()->range_begin)
                                while ((c_valid = lex_ranges[2].next_leaf() != nullptr) && b + min_gap > lex_ranges[2].current_node()->range_begin)
                                    ;
                            if (!c_valid)
                                break;

                            // check whether within max_gap bc
                            if (b + max_gap < lex_ranges[2].current_node()->range_begin) {
                                continue;
                            }
                            c = lex_ranges[2].current_node()->range_begin;

                            // situation: bc still valid

                            // push c greedy beyond max_gap
                            lex_ranges[2].skip_subtree();
                            while (lex_ranges[2].next_leaf() != nullptr) {
                                auto c_pos2 = lex_ranges[2].current_node();
                                if (b + max_gap >= c_pos2->range_begin)
                                    c = c_pos2->range_begin;
                                else
                                    break;
                                lex_ranges[2].skip_subtree();
                            }
                        } else
                            break;
                        lex_ranges[1].skip_subtree();
                    }
                    //--------

                    current = { a, b, c };

                    // pull a forward
                    while (valid() &&
                           lex_ranges[0].current_node()->range_end <= c)
                        lex_ranges[0].skip_subtree();
                    while (lex_ranges[0].next_leaf() != nullptr && lex_ranges[0].current_node()->range_begin <= c + size3)
                        ;

                    return true;
                } else
                    lex_ranges[top1->range_size() >= top0->range_size() ? (top2->range_size() >= top1->range_size() ? 2 : 1) : (top2->range_size() >= top0->range_size() ? 2 : 0)].expand();
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
                                  string_type s1,
                                  string_type s2,
                                  string_type s3,
                                  size_t min_gap,
                                  size_t max_gap)
            : min_gap(min_gap), max_gap(max_gap), size3(s3.size())
        {
            auto root_node = make_shared<sdsl::node_cache<type_index>>(index.wt.root(), index, nullptr, nullptr);
            size_type sp = 1, ep = 0;

            forward_search(index.text.begin(), index.text.end(), index.wt, 0, index.wt.size()-1, s1.begin(), s1.end(), sp, ep);
            lex_ranges.emplace_back(index, sdsl::range_type(sp, ep),root_node);
            //std::cerr << "1 ::: " << sp << "-" << ep << std::endl;
            //std::cerr << "REGEX ::: " << (char)index.text[index.wt[sp]] << (char)index.text[index.wt[sp] + 1] << std::endl;
            //std::cerr << "REGEX ::: " << (char)index.text[index.wt[ep]] << (char)index.text[index.wt[ep] + 1] << std::endl;

            forward_search(index.text.begin(), index.text.end(), index.wt, 0, index.wt.size()-1, s2.begin(), s2.end(), sp, ep);
            lex_ranges.emplace_back(index, sdsl::range_type(sp, ep),root_node);
            //std::cerr << "2 ::: " << sp << "-" << ep << std::endl;
            //std::cerr << "REGEX ::: " << (char)index.text[index.wt[sp]] << (char)index.text[index.wt[sp] + 1] << std::endl;
            //std::cerr << "REGEX ::: " << (char)index.text[index.wt[ep]] << (char)index.text[index.wt[ep] + 1] << std::endl;

            forward_search(index.text.begin(), index.text.end(), index.wt, 0, index.wt.size()-1, s3.begin(), s3.end(), sp, ep);
            lex_ranges.emplace_back(index, sdsl::range_type(sp, ep),root_node);
            //std::cerr << "3 ::: " << sp << "-" << ep << std::endl;
            //std::cerr << "REGEX ::: " << (char)index.text[index.wt[sp]] << (char)index.text[index.wt[sp] + 1] << std::endl;
            //std::cerr << "REGEX ::: " << (char)index.text[index.wt[ep]] << (char)index.text[index.wt[ep] + 1] << std::endl;

            //std::cerr << "GAP ::: " << min_gap << "-" << max_gap << std::endl;

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
            const wild_card_match_iterator3& a,
            const wild_card_match_iterator3& b)
        {
            return a.current == b.current;
        }

        friend bool operator!=(
            const wild_card_match_iterator3& a,
            const wild_card_match_iterator3& b)
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

            std::cerr << "REGEX ::: " << pat.raw_regexp << std::endl;

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

            // linear scan?
            if (text.size() > 0) {
                index_type::size_type total_range = 0, sp = 0, ep = 0;
                for (size_t i = 0; i < pat.subpatterns.size(); ++i)
                    total_range += forward_search(index.text.begin(), index.text.end(), index.wt, 0, index.wt.size()-1, pat.subpatterns[i].begin(), pat.subpatterns[i].end(), sp, ep);

                // check for: |range| * log n > n
                //std::cout << "range = " << total_range << "; |text| = " << text.size() << std::endl;
                total_range *= sdsl::bits::hi(text.size());
                total_range *= CONST_LINEAR_THRESH;
                //std::cout << total_range << " > " << text.size() << " ==> " << (total_range > text.size()) << std::endl;

                if (total_range > text.size()) {
                    // linear scan
                    auto rx = boost::regex(pat.raw_regexp.begin(),pat.raw_regexp.end(),std::regex::ECMAScript);
                    auto matches_begin = boost::sregex_iterator(
                                             text.begin(),
                                             text.end(),
                                             rx,
                                             boost::regex_constants::match_flag_type::match_not_dot_newline);
                    auto matches_end = boost::sregex_iterator();

                    for (boost::sregex_iterator it = matches_begin; it != matches_end; ++it) {
                        res.positions.push_back(it->position());
                    }
                    return res;
                }
            }

            // smart scan
            auto container = sdsl::matching_container<wild_card_match_iterator3<index_type>>(
                                 wild_card_match_iterator3<index_type>(index, s1, s2, s3, min_gap, max_gap),
                                 wild_card_match_iterator3<index_type>());
            for (auto hit : container) {
                res.positions.push_back(std::get<0>(hit));
            }
            return res;
        }
};
