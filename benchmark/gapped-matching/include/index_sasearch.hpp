#pragma once

#include "utils.hpp"
#include "collection.hpp"
#include "sdsl/suffix_arrays.hpp"
#include "sdsl/matching.hpp"

class index_sasearch
{
    private:
        sdsl::int_vector<> m_text;
        sdsl::int_vector<> m_sa;

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
            sdsl::load_from_file(m_text, col.file_map[consts::KEY_TEXT]);

            sdsl::csa_wt<sdsl::wt_huff<>, 1> csa;
            sdsl::construct(csa, col.file_map[consts::KEY_TEXT], 0);
            m_sa = sdsl::int_vector<>(csa.size(), 0, sdsl::bits::hi(csa.size() - 1) + 1);
            std::copy(csa.begin(), csa.end(), m_sa.begin());
        }

        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v=NULL, std::string name="")const
        {
            auto size = m_sa.serialize(out, v, "sa");
            size += m_text.serialize(out, v, "text");
            return size;
        }

        void load(std::istream& in)
        {
            m_sa.load(in);
            m_text.load(in);
        }

        void swap(index_sasearch& ir)
        {
            if (this != &ir) {
                m_sa.swap(ir.m_sa);
                m_text.swap(ir.m_text);
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
            string_type s3;
            size_type min_gap = 0;
            size_type max_gap = 0;

            if (pat.subpatterns.size() <= 2) {
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
                vector<size_type> range_a(forward_search(m_text.begin(), m_text.end(), m_sa, 0, m_sa.size()-1, s1.begin(), s1.end(), sp1, ep1));
                vector<size_type> range_b(forward_search(m_text.begin(), m_text.end(), m_sa, 0, m_sa.size()-1, s2.begin(), s2.end(), sp2, ep2));
                std::copy(m_sa.begin() + sp1, m_sa.begin() + ep1 + 1, range_a.begin());
                std::copy(m_sa.begin() + sp2, m_sa.begin() + ep2 + 1, range_b.begin());
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
            } else {
                s1 = pat.subpatterns[0];
                s2 = pat.subpatterns[1];
                s3 = pat.subpatterns[2];
                min_gap = pat.gaps[0].first + s1.size();
                max_gap = pat.gaps[0].second + s1.size();

                // get ranges
                size_type sp1, ep1;
                size_type sp2, ep2;
                size_type sp3, ep3;
                vector<size_type> range_a(forward_search(m_text.begin(), m_text.end(), m_sa, 0, m_sa.size()-1, s1.begin(), s1.end(), sp1, ep1));
                vector<size_type> range_b(forward_search(m_text.begin(), m_text.end(), m_sa, 0, m_sa.size()-1, s2.begin(), s2.end(), sp2, ep2));
                vector<size_type> range_c(forward_search(m_text.begin(), m_text.end(), m_sa, 0, m_sa.size()-1, s3.begin(), s3.end(), sp3, ep3));
                std::copy(m_sa.begin() + sp1, m_sa.begin() + ep1 + 1, range_a.begin());
                std::copy(m_sa.begin() + sp2, m_sa.begin() + ep2 + 1, range_b.begin());
                std::copy(m_sa.begin() + sp3, m_sa.begin() + ep3 + 1, range_c.begin());
                std::sort(range_a.begin(), range_a.end());
                std::sort(range_b.begin(), range_b.end());
                std::sort(range_c.begin(), range_c.end());

                // linear search
                auto a_it = range_a.begin();
                auto b_it = range_b.begin();
                auto c_it = range_c.begin();

                while (a_it != range_a.end()) {
                    auto a_pos = *a_it;

                    // enforcing min_gap ab
                    bool b_valid;
                    while ((b_valid = (b_it != range_b.end())) && a_pos + min_gap > *b_it)
                        ++b_it;
                    if (!b_valid)
                        break;

                    // check whether within max_gap ab
                    auto b_pos = *b_it;
                    if (a_pos + max_gap < b_pos) {
                        ++a_it;
                        continue;
                    }

                    // enforcing min_gap bc
                    bool c_valid;
                    while ((c_valid = (c_it != range_c.end())) && b_pos + min_gap > *c_it)
                        ++c_it;
                    if (!c_valid)
                        break;

                    // check whether within max_gap bc
                    auto c_pos = *c_it;
                    if (b_pos + max_gap < c_pos) {
                        ++b_it;
                        continue;
                    }

                    // situation: VALID match, but LAZY

                    // push c greedy beyond max_gap
                    ++c_it;
                    while (c_it != range_c.end()) {
                        auto c_pos2 = *c_it;
                        if (b_pos + max_gap >= c_pos2)
                            c_pos = c_pos2;
                        else
                            break;
                        ++c_it;
                    }

                    // situation: bc greedy (c_it IN FRONT!)

                    // push b greedy beyond max_gap
                    ++b_it;
                    while (b_it != range_b.end()) {
                        auto b_pos2 = *b_it;
                        if (a_pos + max_gap >= b_pos2) {
                            b_pos = b_pos2;
                            // situation: found greedier ab ==> check c

                            // enforcing min_gap bc
                            bool c_valid;
                            while ((c_valid = (c_it != range_c.end())) && b_pos + min_gap > *c_it)
                                ++c_it;
                            if (!c_valid)
                                break;

                            // check whether within max_gap bc
                            if (b_pos + max_gap < *c_it) {
                                ++b_it;
                                continue;
                            }
                            c_pos = *c_it;

                            // situation: bc still valid

                            // push c greedy beyond max_gap
                            ++c_it;
                            while (c_it != range_c.end()) {
                                auto c_pos2 = *c_it;
                                if (b_pos + max_gap >= c_pos2)
                                    c_pos = c_pos2;
                                else
                                    break;
                                ++c_it;
                            }
                        } else
                            break;
                        ++b_it;
                    }

                    res.positions.push_back(a_pos);

                    // pull a beyond previous c (non-overlapping)
                    c_pos += s2.size();
                    while (a_it != range_a.end() && *a_it < c_pos)
                        ++a_it;
                }
            }
            return res;
        }
};
