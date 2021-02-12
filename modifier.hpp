#ifndef MODIFIER_HPP
#define MODIFIER_HPP


#include "data_format.hpp"

#include <tuple>

namespace sf_g {
template< class I, class O, class ... Ms > struct modifier {};

template< class I > 
struct modifier< I, waveform > {
    modifier( metadata m_p) : metadata_m{ m_p } {}
    
    std::vector<waveform> operator()(std::vector<I>&& data_pc) const {
        std::vector<waveform> result_c{ static_cast<std::size_t>(metadata_m.channel_count) };
        for( auto& result : result_c) { 
            result.data = TH1D{ "", ";time;ADC", 1024, -metadata_m.sampling_period/2, (1024-1)*metadata_m.sampling_period+ metadata_m.sampling_period/2 };
        }
        for( auto k{0}; k < metadata_m.channel_count ; ++k ){ 
            for( auto i{0}; i < 1024 ; ++i ){
                result_c[k].data.SetBinContent( i+1, data_pc[k].sample_c[i] );
            }
        }
        return result_c;
    }

private:
    metadata const metadata_m;
};

namespace details{
    template< class ... Ts>
    struct pack{ constexpr static std::size_t size = sizeof...(Ts); };
} //namespace details    

template< class ... Modules > 
struct sub_modifier{
    using tuple_t = std::tuple< Modules... >;
    using output_t = composite< typename Modules::output_t... >;

    sub_modifier(Modules&&... ms_p) : module_mc{ std::make_tuple( std::move( ms_p... ) ) } {}
    output_t operator()( waveform&& input_p ) const {
        return call( std::move(input_p), std::make_index_sequence<sizeof...(Modules)>{} );
    }
    private:
    template<std::size_t ... Indices>
    output_t call( waveform&& input_p, std::index_sequence<Indices...> ) const {
        output_t output;
        int expander[] = { 0, ( static_cast< typename std::tuple_element_t<Indices, tuple_t>::output_t&>( output ) = std::get<Indices>(module_mc)( input_p ), void(), 0 ) ... }; 
        return output;
    }  

private:
    tuple_t const module_mc;
};

template< class ... Modules>
sub_modifier<Modules...> make_sub_modifier( Modules&& ... ms_p ) {return {std::move(ms_p)...};}



template< class ... Os, class ... Ms > 
struct modifier< waveform, multi_output<Os...>, Ms... > {
    using tuple_t = std::tuple<Ms...>;

    template< class T = details::pack<Os...>, 
              class Enabler = std::enable_if_t< T::size == sizeof...(Ms)> > 
    constexpr modifier(Ms&& ... sub_modifiers_p) : sub_modifier_mc{ std::make_tuple( sub_modifiers_p...) } {} 
    multi_output<Os...> operator()( std::vector<waveform>&& input_pc ) const {
        return call_impl( std::move(input_pc), std::make_index_sequence< sizeof...(Os)>{} ); 
    }

    private:
    template<std::size_t ... Indices>
    multi_output<Os...> call_impl( std::vector<waveform>&& input_pc, std::index_sequence<Indices...> ) const {
        multi_output<Os...> output;
        int expander[] = { 0, ( static_cast<Os&>( output ) = std::get<Indices>(sub_modifier_mc)( std::move(input_pc[Indices]) ), void(), 0 ) ... }; 
        return output;
    }
private:
    tuple_t const sub_modifier_mc;
};

template< class I, class O, class ... Ms>
modifier<I, O, Ms...> make_modifier( Ms&& ... sub_modifiers_p ) { return {sub_modifiers_p...};}

} //namespace sf_g
#endif
