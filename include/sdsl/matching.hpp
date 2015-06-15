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

template<class t_csa, class t_wt>
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

template<class type_csa, class type_wt>
struct node_cache
{
    typedef typename type_wt::node_type  node_type;
    typedef typename type_csa::size_type size_type;
    
    type_wt wts;
    node_type node;
    pair<shared_ptr<node_cache>, shared_ptr<node_cache>> children;
    size_type range_begin;
    size_type range_end;
    bool is_leaf;
    
    size_type range_size()
    {
        return range_end - range_begin;
    }
    
    node_cache(node_type node, type_wt wts)
    {
        this->wts = wts;
        this->node = node;
        this->children = make_pair(nullptr, nullptr);
        auto range = wts.value_range(node);
        this->range_begin = get<0>(range);
        this->range_end = get<1>(range);
        this->is_leaf = wts.is_leaf(node);
    }
    
    void ensure_children()
    {
        if (children.first == nullptr)
        {
            auto children = wts.expand(node);
            this->children = make_pair(make_shared<node_cache>(children[0], wts), make_shared<node_cache>(children[1], wts));
        }
    }
};

template<class type_csa, class type_wt>
class wild_card_match_iterator : public std::iterator<std::forward_iterator_tag, pair<typename type_csa::size_type, typename type_csa::size_type>>
{
private:
    typedef typename type_wt::node_type  node_type;
    typedef typename type_csa::size_type size_type;
    typedef pair<size_type, size_type>   result_type;

    // (lex_range, node)
    array<stack<pair<range_type,shared_ptr<node_cache<type_csa, type_wt>>> >, 2> lex_ranges;

    type_wt wts;
    size_t a;
    size_t b_idx = -1;
    deque<size_t> b_values;

    incremental_wildcard_pattern p1;

    result_type current;

    void split_node(int t)
    {
        auto top = lex_ranges[t].top(); lex_ranges[t].pop();
        auto& node = top.second;
        node->ensure_children();
        auto exp_range = wts.expand(node->node, top.first);
        if (!empty(exp_range[1]))
            lex_ranges[t].emplace(exp_range[1], node->children.second);
        if (!empty(exp_range[0]))
            lex_ranges[t].emplace(exp_range[0], node->children.first);
    };
    
    void adjust_b_range()
    {
        // shrink
        while (!b_values.empty() && a + p1.min_gap > b_values.front())
            b_values.pop_front();
        // expand
        while (!lex_ranges[1].empty() && a + p1.max_gap >= lex_ranges[1].top().second->range_begin)
        {
            if (lex_ranges[1].top().second->is_leaf)
            {                    
                b_values.push_back(lex_ranges[1].top().second->range_begin);
                lex_ranges[1].pop();
            }
            else
                split_node(1);
        }
    }

    bool next_batch()
    {
        b_idx = 0;
        
        // find next connected match
        while (!lex_ranges[0].empty() 
            && !b_values.empty()
            && lex_ranges[0].top().second->range_begin + p1.min_gap <= b_values.back())
        {
            if (lex_ranges[0].top().second->is_leaf)
            {
                a = lex_ranges[0].top().second->range_begin;
                lex_ranges[0].pop();
                adjust_b_range();
                if (!b_values.empty())
                    return true;
            }
            else
                split_node(0);
        }
            
        // find next independent match
        while (!lex_ranges[0].empty() && !lex_ranges[1].empty())
        {
            const auto& top0 = lex_ranges[0].top().second;
            const auto& top1 = lex_ranges[1].top().second;
            if (top0->range_end + p1.max_gap < top1->range_begin)
                lex_ranges[0].pop();
            else if (top0->range_begin + p1.min_gap > top1->range_end)
                lex_ranges[1].pop();
            else if (top0->is_leaf && top1->is_leaf)
            {
                a = top0->range_begin;
                b_values.push_back(top1->range_begin);
                lex_ranges[0].pop();
                lex_ranges[1].pop();
                adjust_b_range();
                return true;
            }
            else
                split_node(top1->range_size() > top0->range_size() ? 1 : 0);
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
    wild_card_match_iterator(const matching_index<type_csa, type_wt>& index, 
        string s, 
        incremental_wildcard_pattern p1) 
        : wts(index.wt), p1(p1)
    {
        auto root_node = make_shared<node_cache<type_csa, type_wt>>(this->wts.root(), this->wts);
        size_type sp = 1, ep = 0;
        if (0 != backward_search(index.csa, 0, index.csa.size()-1, s.begin(), s.end(), sp, ep))
            lex_ranges[0].emplace(range_type(sp, ep),root_node);
        if (0 != backward_search(index.csa, 0, index.csa.size()-1, p1.pattern.begin(), p1.pattern.end(), sp, ep))
            lex_ranges[1].emplace(range_type(sp, ep),root_node);
        if (p1.except != "")
            if (0 != backward_search(index.csa, 0, index.csa.size()-1, p1.except.begin(), p1.except.end(), sp, ep))
                lex_ranges[2].emplace(range_type(sp, ep),root_node);

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
         class t_wt=wt_int<bit_vector, rank_support_v5<>, select_support_scan<1>, select_support_scan<0>>>
class matching_index
{
    static_assert(std::is_same<typename index_tag<t_csa>::type, csa_tag>::value,
        "First template argument has to be a suffix array.");
    static_assert(std::is_same<typename index_tag<t_wt>::type, wt_tag>::value,
        "Second template argument has to be a wavelet tree.");

private:
    typedef t_csa                        csa_type;
    typedef t_wt                         wt_type;
    typedef typename wt_type::node_type  node_type;
public:
    typedef typename csa_type::size_type size_type;
public:
    typedef wild_card_match_iterator<csa_type, wt_type> iterator;


private:
    const csa_type m_csa;
    const wt_type  m_wt;

public:
    const csa_type& csa = m_csa;
    const wt_type&  wt  = m_wt;
    
    matching_index(const csa_type csa, const wt_type wt)
        : m_csa(csa), m_wt(wt)
    {
    }   
    
    matching_result<iterator> match2(
        const string s,
        const incremental_wildcard_pattern p1
        ) const
    {
        return matching_result<iterator>(
            wild_card_match_iterator<csa_type, wt_type>(*this, s, p1),
            wild_card_match_iterator<csa_type, wt_type>());
    }
};

}// end namespace sdsl
#endif
