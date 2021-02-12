#ifndef WRITER_HPP
#define WRITER_HPP

#include "data_format.hpp"

#include "TTree.h"

namespace sf_g{

template<class I, class O> struct writer {};

template<class I>
struct writer< I, TTree >{
    writer(data_stream<TTree>& tree_p, int channel_count_p  ) : 
        channel_count_m{ static_cast<size_t>( channel_count_p ) }, 
        data_mc{ channel_count_m },
        tree_mh{ tree_p.stream_h } {
            for( auto i{0}; i < channel_count_m ; ++i) {
                std::string name = "channel_" + std::to_string(i);
                tree_mh->Branch( name.c_str(), &data_mc[i] ); 
            } 
        }

    constexpr void operator()( std::vector<I>&& input_pc ) {
        for( auto i{0} ; i < channel_count_m ; ++i ) {
           data_mc[i] = std::move( input_pc[i] ) ;
//           std::cout << data_mc[i].fcr << '\n';
        }
        tree_mh->Fill();
    } 

    private:
    size_t const channel_count_m;
    std::vector<I> data_mc;
    TTree * tree_mh;
};


template<class ... Is>
struct writer< multi_input<Is...>, TTree > {
    writer( data_stream<TTree>& stream_p ) : 
        tree_mh{ stream_p.stream_h} { 
        std::size_t index{0};
        data_mc.apply_for_each(        
                [this, &index]( auto const& data_p ){ 
                    std::string name = "channel_" + std::to_string(index++);
                    tree_mh->Branch( name.c_str(), &data_p ); 
                    //some glue need to be fully compatible -> i.e. readable via root
                }
        ); 
   }

    constexpr void operator()( multi_output<Is...>&& input_pc ) {
        data_mc.apply_for_each(        
                [this, &input_pc]( auto& data_p ){ 
                    data_p = static_cast<decltype(data_p)>(input_pc);
                }
        ); 
        tree_mh->Fill();
    } 

private:
    TTree* tree_mh;
    multi_input<Is...> data_mc; 
};

} //namespace sf_g

#endif
