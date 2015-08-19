#pragma once

#include "utils.hpp"
#include "collection.hpp"
#include "sdsl/suffix_arrays.hpp"
#include "sdsl/matching.hpp"

class index_wcsearch_bfs
{
    private:
        sdsl::matching_index<> index;

    public:
        typedef sdsl::int_vector<0>::size_type size_type;
        std::string name() const
        {
            std::string index_name = IDXNAME;
            return "WCSEARCH-BFS-"+index_name;
        }

    public:
        index_wcsearch_bfs() { }
        index_wcsearch_bfs(collection& col)
        {
            sdsl::cache_config cc(false,".","WCSEARCH_BFS_TMP");
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

        void swap(index_wcsearch_bfs& ir)
        {
            if (this != &ir) {
                index.swap(ir.index);
            }
        }

        gapped_search_result
        search(const gapped_pattern& pat) const
        {
            auto csa = index.csa;
            auto wts = index.wt;
            using node_type = decltype(wts)::node_type;
            using size_type = decltype(csa)::size_type;
            using range_type = sdsl::range_type;

            gapped_search_result res;
            std::string s1;
            std::string s2;
            size_t min_gap = 0;
            size_t max_gap = 0;

            if (pat.subpatterns.size() == 1) {
                s1 = pat.subpatterns[0];
                s2 = pat.subpatterns[0];
            } else {
                s1 = pat.subpatterns[0];
                s2 = pat.subpatterns[1];
                min_gap = pat.gaps[0].first;
                max_gap = pat.gaps[0].second;
            }

            size_t cnt_nodes = 0;
            array<vector<pair<range_type,size_t>>, 2> lex_ranges;
            vector<node_type> nodes;
            nodes.push_back(wts.root());
            size_type sp = 1, ep = 0;
            if ( backward_search(csa, 0, csa.size()-1, s1.begin(), s1.end(), sp, ep) > 0 )
                lex_ranges[0].emplace_back(range_type(sp, ep),0);
            else
                nodes.clear();
            if ( backward_search(csa, 0, csa.size()-1, s2.begin(), s2.end(), sp, ep) > 0 )
                lex_ranges[1].emplace_back(range_type(sp,ep),0);
            else
                nodes.clear();

            cnt_nodes += nodes.size();

            auto _lb = [&lex_ranges, &nodes, &wts](size_type k, size_type i){
                return get<0>( wts.value_range(nodes[lex_ranges[k][i].second]) );
            };

            auto _rb = [&lex_ranges, &nodes, &wts](size_type k, size_type i){
                return get<1>( wts.value_range(nodes[lex_ranges[k][i].second]) );
            };
                
            while ( !nodes.empty() ) {
                   
                /*auto _output = [&nodes, &lex_ranges, &wts](size_t k){
                    cout<<">>>>>>>> k="<<k<<endl;
                    for (size_t i=0; i<lex_ranges[k].size(); ++i) {
                        auto pos_range = wts.value_range(nodes[lex_ranges[k][i].second]);
                        auto range = lex_ranges[k][i].first;
                        cout << "i = "<<i<<" range=["<<get<0>(range)<<","<<get<1>(range)<<"] "; 
                        cout <<  "pos_range=["<<get<0>(pos_range) << "," << get<1>(pos_range)<<"]"<<endl;
                    }
                };
                  
                cout<<"_nodes.size()="<<nodes.size()<<endl;
                cout<<"_lex_ranges[0].size()="<<lex_ranges[0].size()<<std::endl;
                _output(0);
                cout<<"_lex_ranges[1].size()="<<lex_ranges[1].size()<<std::endl;
                _output(1);*/



                if ( wts.is_leaf(nodes[0]) )
                    break;
                vector<node_type> _nodes;
                array<vector<pair<range_type,size_t>>, 2> _lex_ranges;
                std::array<decltype(lex_ranges[0].begin()),2> lex_range_it = {lex_ranges[0].begin(), lex_ranges[1].begin()};
                for (size_t i=0; i<nodes.size(); ++i){
                    auto exp_node = wts.expand(nodes[i]);
                    _nodes.push_back(exp_node[0]);
                    _nodes.push_back(exp_node[1]);
                    for (size_t j=0; j<2; ++j){
                        if (lex_range_it[j] != lex_ranges[j].end() and lex_range_it[j]->second == i){
                            auto exp_range = wts.expand(nodes[i], lex_range_it[j]->first);
                            ++lex_range_it[j];
                            for (size_t k=0; k<2; ++k){
                                if ( !sdsl::empty(exp_range[k]) ){
    //cout<<"k="<<k<<" exp_range["<<k<<"]=["<<get<0>(exp_range[k])<<","<<get<1>(exp_range[k])<<"]"<<endl;
                                    _lex_ranges[j].emplace_back(exp_range[k], _nodes.size()-2+k);    
                                }     
                            }
                        }
                    }
                }
                cnt_nodes += _nodes.size();
                nodes = std::move(_nodes);
                lex_ranges = std::move(_lex_ranges); 

                // filtering
                size_t n0 = 0, n1 = 0;
                for (size_t i=0,j=0; i<lex_ranges[0].size(); ++i){
                    while ( j < lex_ranges[1].size() and _rb(1,j) < _lb(0,i)+s1.size()+min_gap ) {
                        ++j;
                    }
                    if ( j < lex_ranges[1].size() and _lb(1,j) <= _rb(0,i)+s1.size()+max_gap ){
                        // lex_ranges[0][i] overlaps lex_ranges[1][j]
                        lex_ranges[0][n0++] = lex_ranges[0][i];
                    }
                }
                lex_ranges[0].resize(n0);
                for(size_t i=0,j=0; j<lex_ranges[1].size(); ++j){
                    while ( i < lex_ranges[0].size() and _lb(0,n0-i-1)+s1.size()+min_gap > _rb(1,lex_ranges[1].size()-j-1) ){
                        ++i;
                    }
                    if ( i < lex_ranges[0].size() and _lb(1,lex_ranges[1].size()-j-1) <= _rb(0,n0-i-1)+s1.size()+max_gap ) {
                        lex_ranges[1][lex_ranges[1].size()-n1-1] = lex_ranges[1][lex_ranges[1].size()-j-1];
                        ++n1;
                    }
                }
                for(size_t j=0; j<n1; ++j){
                    lex_ranges[1][j] = lex_ranges[1][lex_ranges[1].size()-n1+j];
                }
                lex_ranges[1].resize(n1);
                // remove unused nodes
                n0 = 0;
                for (size_t k=0, i=0, j=0; k<nodes.size(); ++k){
                    bool occurs = false;
                    if ( i < lex_ranges[0].size() and lex_ranges[0][i].second == k ){
                        occurs = true;
                        lex_ranges[0][i++].second = n0;
                    }
                    if ( j < lex_ranges[1].size() and lex_ranges[1][j].second == k ){
                        occurs = true;
                        lex_ranges[1][j++].second = n0;
                    }
                    if ( occurs ) {
                        nodes[n0++] = nodes[k];            
                    }
                }
                nodes.resize(n0);
            }

            for(size_t i = 0; i < lex_ranges[0].size(); ++i) {
                auto pos0 = _lb(0,i);
                for(size_t j = 0; j < lex_ranges[1].size(); ++j) {
                    auto pos1 = _lb(1,j);
                    if (pos0 + s1.size() + min_gap <= pos1
                     && pos0 + s1.size() + max_gap >= pos1)
                        res.positions.push_back(pos0);
                }
            }
            cout << "Processed nodes: " << cnt_nodes << " (" <<wts.size() << ")" << endl;

            return res;
        }
};
