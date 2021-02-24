#ifndef WRITER_HPP
#define WRITER_HPP

#include "data_format.hpp"

#include "TTree.h"

namespace sf_g{

template<class I, class O> struct writer {};

template<class I>
struct writer< I, TTree >{
    writer(data_output<TTree>& tree_p, int channel_count_p  ) : 
        channel_count_m{ static_cast<size_t>( channel_count_p ) }, 
        data_mc{ channel_count_m },
        tree_m{ tree_p.output } {
            for( auto i{0}; i < channel_count_m ; ++i) {
                std::string name = "channel_" + std::to_string(i);
                tree_m.Branch( name.c_str(), &data_mc[i] ); 
            } 
        }

    constexpr void operator()( std::vector<I>&& input_pc ) {
        for( auto i{0} ; i < channel_count_m ; ++i ) {
           data_mc[i] = std::move( input_pc[i] ) ;
//           std::cout << data_mc[i].fcr << '\n';
        }
        tree_m.Fill();
    } 

    private:
    size_t const channel_count_m;
    std::vector<I> data_mc;
    TTree&  tree_m;
};


template<class ... Is>
struct writer< multi_input<Is...>, TTree > {
    writer( data_output<TTree>& output_p ) : 
        tree_m{ output_p.output } 
    { 
        std::size_t index{0};
        data_mc.apply_for_each(        
                [this, &index]( auto& data_p ){ 
                    std::string name = "channel_" + std::to_string(index++);
                    tree_m.Branch( name.c_str(), &data_p ); 
                }
        ); 
   }

    constexpr void operator()( multi_output<Is...>&& input_pc ) {
        data_mc.apply_for_each(        
                [this, &input_pc]( auto& data_p ){ 
                    data_p = static_cast<decltype(data_p)>(input_pc);
                }
        ); 
        tree_m.Fill();
    } 

private:
    TTree& tree_m;
    multi_input<Is...> data_mc; 
};

template<class I>
struct writer< multi_input<I>, TTree > {
    writer( data_output<TTree>& output_p ) : 
        tree_m{ output_p.output }
    { 
        tree_m.Branch( "channel_0", &data_m ); 
    }

    constexpr void operator()( multi_output<I>&& input_pc ) {
        data_m = static_cast<I>(input_pc);
//        data_m.value();
        tree_m.Fill();
    } 

private:
  int required_for_tree_filling{}; //apparently needed for root to fill the tree in the case of an lonely modifier
    TTree& tree_m;
    I data_m; 
};



template<class T> struct multi_input_deducer{};
template< template< class, class, class...> class M, 
          class I, 
          class ... Os, 
          class ... Ms > 
struct multi_input_deducer< M< I, multi_output<Os...>, Ms... >> {
    using type = multi_input<Os...> ;
};

template< class O, class M, class I = typename multi_input_deducer<M>::type>
writer<I, O> make_multi_writer( M const&, data_output<TTree>& output_p ) {return {output_p};}

} //namespace sf_g

#endif
