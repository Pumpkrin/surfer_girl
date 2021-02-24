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

    sub_modifier(Modules&&... ms_p) : module_mc{ std::make_tuple( std::move( ms_p )...  ) } {}
    output_t operator()( waveform&& input_p ) const {
        return call( std::move(input_p), std::make_index_sequence<sizeof...(Modules)>{} );
    }
    private:
    template<std::size_t ... Indices>
    output_t call( waveform&& input_p, std::index_sequence<Indices...> ) const {
        output_t output;
        int expander[] = { 0, ( static_cast< typename std::tuple_element_t<Indices, tuple_t>::output_t&>( output ) = std::get<Indices>(module_mc)( std::move(input_p) ), void(), 0 ) ... }; 
        return output;
    }  

private:
    tuple_t const module_mc;
};

template< class M > 
struct sub_modifier<M>{
    using output_t = typename M::output_t;

    sub_modifier(M&& m_p) : module_m{ std::move(m_p) } {}
    output_t operator()( waveform&& input_p ) const {
        return module_m( std::move(input_p) );
    }

private:
    M const module_m;
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

template<class ... Ts>
struct multi_output_deducer{
    using type = multi_output< typename Ts::output_t...>; 
};


template< class I, class ... Ms, class O = typename multi_output_deducer<Ms...>::type>
modifier<I, O, Ms...> make_multi_modifier( Ms&& ... sub_modifiers_p ) { return {std::move(sub_modifiers_p)...};}


struct amplitude_finder {
    using output_t = amplitude;
    output_t operator()( waveform&& input_p ) const {
        double baseline{0}; 
        for( auto i{0}; i < 16 ; ++i ) { baseline += input_p.data.GetBinContent( i +1 ); }
        baseline /= 16;
        return {baseline - input_p.data.GetMinimum()}; 
    };
}; 

struct baseline_finder {
    using output_t = baseline;
    output_t operator()( waveform&& input_p ) const {
        double baseline{0};
        for( auto i{0}; i < 16 ; ++i ) { baseline += input_p.data.GetBinContent( i +1 ); }
        baseline /= 16;
        return {baseline};
    }
};


struct cfd_calculator {
    using output_t = cfd_time;
    output_t operator()(waveform&& input_p) const {
        double fraction = 0.4;
        std::size_t delay = 15;

        double baseline{0};
        for( auto i{0}; i < 16 ; ++i ) { baseline += input_p.data.GetBinContent( i +1 ); }
        baseline /= 16;

        double current{0}, last{0};
//        std::cout << "cfd:\n";
        for(auto i{1}; i < input_p.data.GetNbinsX()+1 - delay; ++i) {
            double value = input_p.data.GetBinContent(i) - baseline; 
            double df_value = fraction * ( input_p.data.GetBinContent(i + delay) - baseline );
//            std::cout << "value: " << value << " - df_value: " << df_value << '\n';
            
            current = value - df_value;  
            if( current * last < 0 ){ 
                double slope = (current - last)/( input_p.data.GetBinCenter(i) - input_p.data.GetBinCenter(i-1));
                double offset = current - input_p.data.GetBinCenter(i)*slope;
//                std::cout << "slope: " << slope << " - offset: " << offset << '\n';
//                std::cout << "time: " << -offset/slope << '\n';
                return {-offset/slope}; 
            }     
            last = current; 
        }
        return {-1};
    }
};

struct charge_integrator {
    using output_t = charge;
    output_t operator()(waveform&& input_p) const {

        double baseline{0};
        for( auto i{0}; i < 16 ; ++i ) { baseline += input_p.data.GetBinContent( i +1 ); }
        baseline /= 16;

        double waveform_length = input_p.data.GetBinCenter(input_p.data.GetXaxis()->GetLast());
//        std::cout << "rectangle_integral: " << baseline*waveform_length << '\n';
//        std::cout << "signal_integral: " << input_p.data.Integral("width") << '\n';
//        std::cout << "waveform_length: " << waveform_length << '\n';
        charge result{baseline*waveform_length - input_p.data.Integral("width")};  
        return result.charge > 0 ? result : charge{0};
    }
};

} //namespace sf_g
#endif
