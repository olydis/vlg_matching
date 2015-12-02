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

        inline node_type current_node() const
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
class wild_card_match_iterator3 : public std::iterator<std::forward_iterator_tag, typename type_index::size_type>
{
    private:
        typedef typename type_index::node_type node_type;
        typedef typename type_index::size_type size_type;
        typedef typename type_index::wt_type   wt_type;
        typedef typename type_index::size_type result_type;

        // (lex_range, node)
        vector<wavelet_tree_range_walker<type_index>> lex_ranges;

        size_t min_gap;
        size_t max_gap;

        size_t size3;

        result_type current;

        bool next()
        {
            // find next independent match
            while (valid()) {
                bool skip = false;
                for (size_t i = 1; i < lex_ranges.size(); ++i) {
                    if (lex_ranges[i - 1].current_node()->range_end + max_gap < lex_ranges[i].current_node()->range_begin) {
                        lex_ranges[i - 1].skip_subtree();
                        skip = true;
                        if (!lex_ranges[i - 1].has_more())
                            break;
                    }
                    if (lex_ranges[i - 1].current_node()->range_begin + min_gap > lex_ranges[i].current_node()->range_end) {
                        lex_ranges[i].skip_subtree();
                        skip = true;
                        if (!lex_ranges[i].has_more())
                            break;
                    }
                }  

                if (skip)
                    continue;

                size_t r = 1;
                size_t j;
                skip = false;
                for (size_t i = 0; i < lex_ranges.size(); ++i) {
                    auto lr = lex_ranges[i].current_node()->range_size();
                    if (lr > r) {
                        r = lr;
                        j = i;
                        skip = true;
                    }
                }

                if (skip)
                    lex_ranges[j].expand();
                else {
                    current = lex_ranges[0].current_node()->range_begin;

                    auto x = lex_ranges[lex_ranges.size() - 1].current_node()->range_begin;

                    // pull a forward
                    while (valid() && lex_ranges[0].current_node()->range_end <= x) lex_ranges[0].skip_subtree();
                    while (lex_ranges[0].next_leaf() != nullptr && lex_ranges[0].current_node()->range_begin <= x + size3) ;

                    return true;
                }
            }

            current = (result_type)-1;
            return false;
        }

    public:
        wild_card_match_iterator3()
        {
            current = (result_type)-1;
        }
        wild_card_match_iterator3(const type_index& index,
                                  vector<string_type>& s,
                                  size_t min_gap,
                                  size_t max_gap)
            : min_gap(min_gap), max_gap(max_gap), size3(s[s.size() - 1].size())
        {
            auto root_node = make_shared<sdsl::node_cache<type_index>>(index.wt.root(), index, nullptr, nullptr);
            size_type sp = 1, ep = 0;

            for (auto sx : s) {
                forward_search(index.text.begin(), index.text.end(), index.wt, 0, index.wt.size()-1, sx.begin(), sx.end(), sp, ep);
                lex_ranges.emplace_back(index, sdsl::range_type(sp, ep), root_node);
std::cerr << std::string(sx.begin(), sx.end()) << ": " << sp << " " << ep << std::endl;
            }

            next();
        }

        bool valid() const
        {
            for (auto lr : lex_ranges)
                if (!lr.has_more())
                    return false;
            return true;
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
            vector<string_type> s;
            size_type min_gap;
            size_type max_gap;

            std::cerr << "REGEX ::: " << pat.raw_regexp << std::endl;

            s.push_back(pat.subpatterns[0]);
            s.push_back(pat.subpatterns[1]);
            for (size_t i = 2; i < NUM_PATTERNS; ++i)
                s.push_back(pat.subpatterns[1]);

            min_gap = s[0].size() + pat.gaps[0].first;
            max_gap = s[0].size() + pat.gaps[0].second;

            // smart scan
            auto container = sdsl::matching_container<wild_card_match_iterator3<index_type>>(
                                 wild_card_match_iterator3<index_type>(index, s, min_gap, max_gap),
                                 wild_card_match_iterator3<index_type>());
            for (auto hit : container) {
                res.positions.push_back(hit);
            }
            return res;
        }
};
