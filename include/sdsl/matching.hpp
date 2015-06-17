/* sdsl - succinct data structures library
    Copyright (C) 2011-2013 Simon Gog

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*! \file matching.hpp
    \brief matching.hpp contains methods for advanced matching into 
      indexed data.
    \author Johannes Bader
*/
#ifndef INCLUDED_SDSL_MATCHING
#define INCLUDED_SDSL_MATCHING

#include "sdsl_concepts.hpp"
#include "int_vector.hpp"
#include "sd_vector.hpp"// for standard initialisation of template parameters
#include "util.hpp"
#include "wt_huff.hpp"
#include <algorithm> // for std::swap
#include <vector>
#include <queue>
#include <iostream>

using namespace std;

//! Namespace for the succinct data structure library.
namespace sdsl
{

template<class t_csa, class t_wt, class t_bv>
class matching_index;

struct incremental_wildcard_pattern
{
    const string pattern;
    const size_t min_gap;
    const size_t max_gap;
    const string except;

    incremental_wildcard_pattern() : 
            pattern(""), 
            min_gap(0),
            max_gap(0),
            except("")
    { }

    incremental_wildcard_pattern(
        const string pattern,
        const size_t min_gap,
        const size_t max_gap,
        const string except = "") : 
            pattern(pattern),
            min_gap(min_gap),
            max_gap(max_gap),
            except(except)
    {
    }
};

template<class t_iter>
class matching_result
{
private:
    t_iter m_begin;
    t_iter m_end;

public:
    matching_result(t_iter begin, t_iter end) : m_begin(begin), m_end(end)
    {
    }

    t_iter begin()
    {
        return m_begin;
    } 

    t_iter end()
    {
        return m_end;
    }
};

template<class type_index>
struct node_cache
{
    typedef typename type_index::node_type node_type;
    typedef typename type_index::size_type size_type;
    
    const type_index& index;
    node_type node;
    pair<shared_ptr<node_cache>, shared_ptr<node_cache>> children;
    size_type range_begin;
    size_type range_end;
    size_type range_begin_doc;
    size_type range_end_doc;
    bool is_leaf;
    
    size_type range_size()
    {
        return range_end - range_begin;
    }
    
    node_cache(
        node_type node, 
        const type_index& index, 
        size_type* range_begin_doc, 
        size_type* range_end_doc)
        : index(index)
    {
        this->node = node;
        this->children = make_pair(nullptr, nullptr);
        auto range = index.wt.value_range(node);
        this->range_begin = get<0>(range);
        this->range_end = get<1>(range);
        this->range_begin_doc = 
            range_begin_doc == nullptr ? index.get_document_index(range_begin) : *range_begin_doc;
        this->range_end_doc = 
            range_end_doc == nullptr ? index.get_document_index(range_end) : *range_end_doc;
        this->is_leaf = index.wt.is_leaf(node);
    }
    
    void ensure_children()
    {
        if (children.first == nullptr)
        {
            size_type* range_begin_doc = &this->range_begin_doc;
            size_type* range_center_doc = this->range_begin_doc == this->range_end_doc ? range_begin_doc : nullptr;
            size_type* range_end_doc = &this->range_end_doc;
            
            auto children = index.wt.expand(node);
            this->children = make_pair(
                make_shared<node_cache>(children[0], index, range_begin_doc, range_center_doc), 
                make_shared<node_cache>(children[1], index, range_center_doc, range_end_doc));
        }
    }
};

template<class type_index>
class wavelet_tree_range_walker
{
private:
    typedef shared_ptr<node_cache<type_index>> node_type;
    const type_index& index;
    stack<pair<range_type,node_type>> dfs_stack;
        
public:
    wavelet_tree_range_walker(const type_index& index, range_type initial_range, node_type root_node)
        : index(index)
    {
        dfs_stack.emplace(initial_range, root_node);
    }
    
    bool has_more()
    {
        return !dfs_stack.empty();
    }
    
    node_type current_node()
    {
        return dfs_stack.top().second;
    }
    
    // moves to the next (adjacent) node in the iterator
    void next()
    {
        dfs_stack.pop();
    }

    // splits the current node into two adjacent nodes with the same value range
    void split()
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
        if (node->is_leaf)
        {
            next();
            return node;
        }
        split();
        return nullptr;
    }
};

template<class type_index>
class wild_card_match_iterator : public std::iterator<std::forward_iterator_tag, pair<typename type_index::size_type, typename type_index::size_type>>
{
private:
    typedef typename type_index::node_type node_type;
    typedef typename type_index::size_type size_type;
    typedef typename type_index::wt_type   wt_type;
    typedef pair<size_type, size_type>     result_type;

    // (lex_range, node)
    vector<wavelet_tree_range_walker<type_index>> lex_ranges;

    wt_type wts;
    size_t a;
    size_t a_doc;
    size_t b_idx = -1;
    deque<size_t> b_values;

    incremental_wildcard_pattern p1;

    result_type current;
    
    void adjust_b_range()
    {
        // shrink
        while (!b_values.empty() && a + p1.min_gap > b_values.front())
            b_values.pop_front();
        // expand
        while (lex_ranges[1].has_more() 
            && a + p1.max_gap >= lex_ranges[1].current_node()->range_begin
            && a_doc == lex_ranges[1].current_node()->range_begin_doc)
        {
            auto leaf = lex_ranges[1].retrieve_leaf_and_traverse();
            if (leaf != nullptr)
                b_values.push_back(leaf->range_begin);
        }
    }

    bool next_batch()
    {
        b_idx = 0;
        
        // find next connected match
        while (lex_ranges[0].has_more() 
            && !b_values.empty()
            && lex_ranges[0].current_node()->range_begin + p1.min_gap <= b_values.back())
        {
            auto leaf = lex_ranges[0].retrieve_leaf_and_traverse();
            if (leaf != nullptr)
            {
                a = leaf->range_begin;
                a_doc = leaf->range_begin_doc;
                adjust_b_range();
                if (!b_values.empty())
                    return true;
            }
        }
            
        // find next independent match
        while (lex_ranges[0].has_more() && lex_ranges[1].has_more())
        {
            const auto& top0 = lex_ranges[0].current_node();
            const auto& top1 = lex_ranges[1].current_node();
            if (top0->range_end_doc < top1->range_begin_doc)
                lex_ranges[0].next();
            else if (top0->range_end + p1.max_gap < top1->range_begin)
                lex_ranges[0].next();
            else if (top0->range_begin + p1.min_gap > top1->range_end)
                lex_ranges[1].next();
            else if (top0->is_leaf && top1->is_leaf)
            {
                a = top0->range_begin;
                a_doc = top0->range_begin_doc;
                b_values.push_back(top1->range_begin);
                lex_ranges[0].next();
                lex_ranges[1].next();
                adjust_b_range();
                return true;
            }
            else
                lex_ranges[top1->range_size() > top0->range_size() ? 1 : 0].split();
        }

        b_values.clear();
        return false;
    }

    void next()
    {
        ++b_idx;
        if (!valid() && !next_batch())
            return; // end

        current = make_pair(a, b_values[b_idx]);
    }

public:
    wild_card_match_iterator()
    {
    }
    wild_card_match_iterator(const type_index& index, 
        string s, 
        incremental_wildcard_pattern p1) 
        : wts(index.wt), p1(p1)
    {
        auto root_node = make_shared<node_cache<type_index>>(this->wts.root(), index, nullptr, nullptr);
        size_type sp = 1, ep = 0;
        
        backward_search(index.csa, 0, index.csa.size()-1, s.begin(), s.end(), sp, ep);
        lex_ranges.emplace_back(index, range_type(sp, ep),root_node);
        
        backward_search(index.csa, 0, index.csa.size()-1, p1.pattern.begin(), p1.pattern.end(), sp, ep);
        lex_ranges.emplace_back(index, range_type(sp, ep),root_node);
        
        if (p1.except != "")
        {
            backward_search(index.csa, 0, index.csa.size()-1, p1.except.begin(), p1.except.end(), sp, ep);
            lex_ranges.emplace_back(index, range_type(sp, ep),root_node);
        }

        next();
    }

    bool valid() const
    {
        return b_idx < b_values.size();
    }

    result_type operator*() const
    {
        return current;
    }
    result_type* operator->()
    {
        return &current;
    }

    wild_card_match_iterator& operator++()
    {
        next();
        return *this;
    }
    
    friend bool operator==(
        const wild_card_match_iterator& a,
        const wild_card_match_iterator& b)
    {
        if (!a.valid() && !b.valid())
            return true;
        if (!a.valid() || !b.valid())
            return false;
        if (a.a != b.a)
            return false;
        if (a.b_idx != b.b_idx)
            return false;
        return a.current == b.current;
    }
    
    friend bool operator!=(
        const wild_card_match_iterator& a,
        const wild_card_match_iterator& b)
    {
        return !(a == b);
    }
};

template<class t_csa=csa_wt<wt_huff<rrr_vector<63>>>, 
         class t_wt=wt_int<bit_vector, rank_support_v5<>, select_support_scan<1>, select_support_scan<0>>,
         class t_bv=rrr_vector<>>
class matching_index
{
    static_assert(std::is_same<typename index_tag<t_csa>::type, csa_tag>::value,
        "First template argument has to be a suffix array.");
    static_assert(std::is_same<typename index_tag<t_wt>::type, wt_tag>::value,
        "Second template argument has to be a wavelet tree.");
    static_assert(std::is_same<typename index_tag<t_bv>::type, bv_tag>::value,
        "Third template argument has to be a bitvector.");

private:
    typedef matching_index<t_csa, t_wt, t_bv> index_type;
public:
    typedef t_csa                         csa_type;
    typedef t_wt                          wt_type;
    typedef t_bv                          bv_type;
    typedef typename bv_type::rank_1_type rank_type;
    typedef typename wt_type::node_type   node_type;
    typedef typename csa_type::size_type  size_type;
    
    typedef wild_card_match_iterator<index_type> iterator;


private:
    const csa_type  m_csa;
    const wt_type   m_wt;
    const bv_type   m_dbs; // 1 marks the END of a document
    const rank_type m_dbs_rank;

public:
    const csa_type& csa = m_csa;
    const wt_type&  wt  = m_wt;
    
    matching_index(const csa_type csa, const wt_type wt, const bv_type dbs)
        : m_csa(csa), m_wt(wt), m_dbs(dbs), m_dbs_rank(&m_dbs)
    {
    }
    
    size_type get_document_index(size_type symbol_index) const
    {
        symbol_index = std::min(symbol_index, m_dbs.size());
        return m_dbs_rank.rank(symbol_index);
    }
    
    matching_result<iterator> match2(
        const string s,
        const incremental_wildcard_pattern p1
        ) const
    {
        return matching_result<iterator>(
            wild_card_match_iterator<index_type>(*this, s, p1),
            wild_card_match_iterator<index_type>());
    }
};

}// end namespace sdsl
#endif
