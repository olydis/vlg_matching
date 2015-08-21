#pragma once

#include "utils.hpp"
#include "collection.hpp"
#include "sdsl/int_vector.hpp"

#include "bit_streams.hpp"
#include "eliasfano_skip_list.hpp"
#include "intersection.hpp"

#include <regex>

union qid_type {
    uint8_t u8id[8];
    uint64_t u64id;
};

template
<
    uint8_t t_q = 3,
    class t_list_type = eliasfano_skip_list<true,true,false>
    >
class index_qgram_regexp
{
        static_assert(t_q <= 8,"q-gram index only supported for q <= 8");

    public:
        enum { q = t_q };
        typedef sdsl::int_vector<0>::size_type size_type;
        typedef sdsl::int_vector<0>::value_type value_type;
        typedef std::string text_type;
        typedef t_list_type comp_list_type;
        std::string name() const
        {
            std::string index_name = IDXNAME;
            return "QGRAM-"+std::to_string(q)+"-"+index_name;
        }
    protected:
        std::regex rx;
        text_type m_text;
        std::unordered_map<uint64_t,uint64_t> m_qgram_lists;
        sdsl::bit_vector m_list_data;
        bit_istream m_list_strm;
    public:
        index_qgram_regexp() : m_list_strm(m_list_data) { }
        index_qgram_regexp(collection& col) : m_list_strm(m_list_data)
        {
            sdsl::int_vector_mapper<0> sdsl_text(col.file_map[consts::KEY_TEXT]);
            std::copy(sdsl_text.begin(),sdsl_text.end(),std::back_inserter(m_text));
            LOG(INFO) << "START QGRAM CONSTRUCTION!";
            // space inefficient construction for now!
            {
                std::unordered_map<uint64_t,std::vector<uint64_t>> tmp_qgram_lists;
                auto itr = sdsl_text.begin();
                auto end = sdsl_text.end() - (q-1); // TODO: problem when test.size() < q

                size_type cur_pos = 0;
                union qid_type qid;

                while (itr != end) {
                    qid.u64id = 0;
                    std::copy(itr,itr+q,std::begin(qid.u8id));
                    tmp_qgram_lists[qid.u64id].push_back(cur_pos);
                    cur_pos++; ++itr;
                }
                // compress the lists
                bit_ostream bvo(m_list_data);
                auto litr = tmp_qgram_lists.begin();
                while (litr != tmp_qgram_lists.end()) {
                    auto qid = litr->first;
                    auto& list = litr->second;
                    // LOG(INFO) << "BUILD LIST FOR QID=" << qid << " SIZE = " << list.size();
                    auto bv_offset = comp_list_type::create(bvo,list.begin(),list.end());
                    m_qgram_lists.emplace(qid,bv_offset);
                    litr = tmp_qgram_lists.erase(litr);
                }
            }
            m_list_strm.refresh(); // ugly but necessary for now
        }

        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v=NULL, std::string name="")const
        {
            sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            sdsl::int_vector<64> qids(m_qgram_lists.size()*2);
            size_t i = 0;
            for (const auto& pql : m_qgram_lists) {
                qids[i++] = pql.first;
                qids[i++] = pql.second;
            }
            size_type written_bytes = 0;
            written_bytes += qids.serialize(out,child,"qgram mapping");
            written_bytes += m_list_data.serialize(out,child,"qgram lists");
            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        void load(std::istream& in)
        {
            sdsl::int_vector<64> qids;
            qids.load(in);
            for (size_t i=0; i<qids.size(); i+=2) {
                m_qgram_lists.emplace(qids[i],qids[i+1]);
            }
            m_list_data.load(in);
            m_list_strm.refresh(); // ugly but necessary for now
        }

        void swap(index_qgram_regexp& ir)
        {
            if (this != &ir) {
                m_text.swap(ir.m_text);
                m_qgram_lists.swap(ir.m_qgram_lists);
                m_list_data.swap(ir.m_list_data);
            }
        }

        std::vector<uint64_t> str_to_qids(const std::string& str) const
        {
            std::vector<uint64_t> qids;
            union qid_type qid;
            auto itr = str.begin();
            auto end = str.end() - (q-1);
            while (itr != end) {
                qid.u64id = 0;
                std::copy(itr,itr+q,std::begin(qid.u8id));
                qids.push_back(qid.u64id);
                ++itr;
            }
            return qids;
        }

        void info(const gapped_pattern& pat) const { (void)pat; }

        void prepare(const gapped_pattern& pat) 
        { 
            /* (1) construct regexp */
            rx = std::regex(pat.raw_regexp.begin(),pat.raw_regexp.end(),REGEXP_TYPE);
        }

        //! Search for the k documents which contain the search term most frequent
        gapped_search_result
        search(const gapped_pattern& pat) const
        {
            std::cout << "search(" << pat.raw_regexp << ")" << std::endl;
            gapped_search_result res;

            if (pat.subpatterns.size() == 1) {
                // actually not a gapped pattern. SLOW FOR NOW!
                auto matches_begin = std::sregex_iterator(m_text.begin(),m_text.end(),rx);
                auto matches_end = std::sregex_iterator();
                for (std::sregex_iterator it = matches_begin; it != matches_end; ++it) {
                    res.positions.push_back(it->position());
                }
                return res;
            }

            /* (1) construct regexp */


            /* extract the different q-grams from the subpatterns */
            std::vector<uint64_t> potential_start_positions;
            size_t pat_start_offset = 0;
            for (size_t j=0; j<pat.subpatterns.size(); j++) {
                const auto& subp = pat.subpatterns[j];
                if (subp.size() < q) { // UNION over lists. TODO!
                    LOG(INFO) << "skip small q-gram for now: " << subp;
                } else { // Intersection over lists
                    auto qids = str_to_qids(subp);   // get vector of q-grams for subpattern subp
                    // LOG(INFO) << "subpattern = '" << subp << "' qgram-ids = " << qids;
                    /* for each qgram get the list */
                    std::vector<typename comp_list_type::list_type> plists;
                    std::vector<offset_proxy_list<typename comp_list_type::list_type>> lists;
                    size_type i = 0;
                    for (const auto& qid : qids) {
                        auto litr = m_qgram_lists.find(qid);
                        if (litr == m_qgram_lists.end()) {
                            // q-gram does not exist. no results possible -> return
                            return res;
                        } else {
                            auto list_offset = litr->second;
                            plists.emplace_back(comp_list_type::materialize(m_list_strm,list_offset));
                            lists.emplace_back(offset_proxy_list<typename comp_list_type::list_type>(plists.back(),i++));
                        }
                    }
                    if (lists.size() > 1) { // intersect lists of q-grams of subpatterns if there is more than on list (= qids.size() > 1)
                        auto ires = pos_intersect(lists);
                        if (potential_start_positions.empty() || ires.size() < potential_start_positions.size()) {
                            potential_start_positions.clear();
                            for (size_t l=0; l<ires.size(); l++) {
                                potential_start_positions.push_back(ires[l]-pat_start_offset);
                            }
                        }
                    } else { // no intersection required if there is only list (=qids.size()==1)

                        if (potential_start_positions.empty() || lists[0].size() < potential_start_positions.size()) {
                            potential_start_positions.clear();
                            auto itr = lists[0].begin();
                            auto end = lists[0].end();
                            while (itr != end) {
                                potential_start_positions.push_back(*itr - pat_start_offset);
                                ++itr;
                            }
                        }
                    }
                }
                pat_start_offset += subp.size();
                if (j != pat.gaps.size()) pat_start_offset += pat.gaps[j].second; // use min here?
                // LOG(INFO) << "pat_start_offset = " << pat_start_offset;
            }
            auto max_pattern_len = pat_start_offset;

            /* sort potential positions */
            std::sort(potential_start_positions.begin(),potential_start_positions.end());
            //LOG(INFO) << "potential_start_positions = " << potential_start_positions;

            /* (2) find all matching pos */
            if (potential_start_positions.empty()) { // case where we only have subpatterns smaller than q!
                auto matches_begin = std::sregex_iterator(m_text.begin(),m_text.end(),rx);
                auto matches_end = std::sregex_iterator();
                for (std::sregex_iterator it = matches_begin; it != matches_end; ++it) {
                    res.positions.push_back(it->position());
                }
            } else {
                for (const auto& start_pos : potential_start_positions) {
                    if (start_pos > m_text.size()) // bugfix by Johannes: start_pos seems to be a very very very large (= negative?) number sometimes => segfault in line below
                        continue;
                    auto matches_begin = std::sregex_iterator(m_text.begin()+start_pos,m_text.begin()+start_pos+max_pattern_len,rx);
                    auto matches_end = std::sregex_iterator();
                    for (std::sregex_iterator it = matches_begin; it != matches_end; ++it) {
                        res.positions.push_back(start_pos+it->position());
                    }
                }
            }

            return res;
        }
};
