#pragma once

#include "utils.hpp"
#include "collection.hpp"
#include "sdsl/int_vector.hpp"

#include "bit_streams.hpp"
#include "eliasfano_skip_list.hpp"

#include <regex>

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
        typedef t_list_type list_type;
        std::string name() const
        {
            std::string index_name = IDXNAME;
            return "QGRAM-"+std::to_string(q)+"-"+index_name;
        }
    protected:
        text_type m_text;
        std::unordered_map<uint64_t,uint64_t> m_qgram_lists;
        sdsl::bit_vector m_list_data;
        bit_istream m_list_strm;
    public:
        index_qgram_regexp() : m_list_strm(m_list_data) { }
        index_qgram_regexp(collection& col) : m_list_strm(m_list_data)
        {
            sdsl::int_vector_mapper<0> sdsl_text(col.file_map[consts::KEY_TEXT]);
            LOG(INFO) << "START QGRAM CONSTRUCTION!";
            // space inefficient construction for now!
            {
                std::unordered_map<uint64_t,std::vector<uint64_t>> tmp_qgram_lists;
                auto itr = sdsl_text.begin();
                auto end = sdsl_text.end() - (q-1);
                union qid_type {
                    uint8_t u8id[8];
                    uint64_t u64id;
                };
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
                    auto bv_offset = list_type::create(bvo,list.begin(),list.end());
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

        //! Search for the k documents which contain the search term most frequent
        gapped_search_result
        search(const gapped_pattern& pat) const
        {
            /* (1) construct regexp */
            std::regex rx(pat.raw_regexp.begin(),pat.raw_regexp.end(),REGEXP_TYPE);

            /* extract the different q-grams */

            /* (2) find all matching pos */
            auto matches_begin = std::sregex_iterator(m_text.begin(),m_text.end(),rx);
            auto matches_end = std::sregex_iterator();


            /* (3) copy the output */
            gapped_search_result res;
            for (std::sregex_iterator it = matches_begin; it != matches_end; ++it) {
                res.positions.push_back(it->position());
            }
            return res;
        }
};
