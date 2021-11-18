#ifndef MODIFIER_HPP
#define MODIFIER_HPP


#include "data_structure.hpp"

#include <tuple>
#include <vector>
#include <string>
#include <regex>

namespace sf_g {
struct raw_modifier {
    raw_modifier( metadata m_p) : metadata_m{ m_p } {}
    
    std::vector<linked_waveform> operator()(std::vector<raw_waveform>&& data_pc) const {
        std::vector<linked_waveform> result_c{ static_cast<std::size_t>(metadata_m.channel_count) };
        for( auto k{0}; k < metadata_m.channel_count ; ++k ){ 
            result_c[k].data = TH1D{ "", ";time;ADC", 1024, -metadata_m.sampling_period/2, (1024-1)*metadata_m.sampling_period+ metadata_m.sampling_period/2 };
            result_c[k].channel_number = data_pc[k].channel_id;
            for( auto i{0}; i < 1024 ; ++i ){
                result_c[k].data.SetBinContent( i+1, data_pc[k].sample_c[i] );
            }
//            std::cout << "channel" << data_pc[k].channel_id << ": " << result_c[k].data.GetMinimum() << '\n';
        }
        return result_c;
    }

private:
    metadata const metadata_m;
};

struct cut_modifier{
    explicit cut_modifier( std::string const& input_file_p ) :
        channel_number_m{ retrieve_channel_number( input_file_p ) },
        event_mc{ retrieve_events( input_file_p ) }
    {}

    std::vector<linked_waveform> operator()( std::vector< linked_waveform >&& input_pc ) {
        if( counter_m++ != event_mc.front() ){ return {}; }  
        event_mc.erase( event_mc.begin() );
        return std::move( input_pc );            
    }

private:
    std::vector<std::size_t> retrieve_events( std::string const& input_file_p ) {
        std::ifstream input{ input_file_p, std::ios_base::binary };
        if( !input.is_open() ){ throw std::invalid_argument{ "unable to read the given file" }; }
        auto compute_length_l = [&input](){ 
            input.seekg(0, std::ios_base::end );
            auto result = input.tellg();
            input.seekg(0, std::ios_base::beg);
            return result;
                                          };
        auto const sentinel = compute_length_l();
        auto end_is_reached_l = [&input, &sentinel](){
            return input.eof() || ((sentinel - input.tellg()) == 0);
                                                     }; 
        
        std::vector< std::size_t > result_c;
        result_c.reserve( sentinel / sizeof(int) );
        while( !end_is_reached_l() ){
            int event;
            input.read(reinterpret_cast<char*>(&event), sizeof event );
            result_c.push_back( event);  
        } 
        return result_c;
    }

    std::size_t retrieve_channel_number( std::string const& input_file_p) const {
        std::smatch result;
        std::regex_search( input_file_p, result, std::regex{"[0-9](?=_!?[abcrtms])"});
        return static_cast<std::size_t>( std::stoi(result[0].str()) );
    }

private:
    std::size_t const channel_number_m;
    std::vector<std::size_t> event_mc;
    std::size_t counter_m{0};
};

namespace details{
    template< class ... Ts>
    struct pack{ constexpr static std::size_t size = sizeof...(Ts); };
} //namespace details    

template<class T> struct modifier_impl{};
template< class ... Modules > 
struct modifier_impl< details::pack< Modules...> >{
    using tuple_t = std::tuple< Modules... >;
    using output_t = composite< typename Modules::output_t... >;

    modifier_impl() = default;
    output_t operator()( linked_waveform const& input_p ) const {
        return modify_impl( input_p, std::make_index_sequence<sizeof...(Modules)>{} );
    }
    private:
    template<std::size_t ... Indices>
    output_t modify_impl( linked_waveform const& input_p, std::index_sequence<Indices...> ) const {
        output_t output;
        int expander[] = { 0, ( static_cast< typename std::tuple_element_t<Indices, tuple_t>::output_t&>( output ) = std::get<Indices>(module_mc)( input_p ), void(), 0 ) ... }; 
        return output;
    }  

private:
    tuple_t const module_mc;
};

template< class M > 
struct modifier_impl< details::pack<M> >{
    using output_t = typename M::output_t;

    output_t operator()( linked_waveform const& input_p ) const {
        return module_m( input_p );
    }

private:
    M const module_m{};
};

template< class Specifier >
using modifier = modifier_impl< typename Specifier::pack >;
 

//form specifier to list of Modules ?

struct amplitude_finder {
    using output_t = amplitude;
    output_t operator()( linked_waveform const& input_p ) const {
        double baseline{0}; 
        for( auto i{0}; i < 64 ; ++i ) { baseline += input_p.data.GetBinContent( i +1 ); }
        baseline /= 64;
//        std::cout << "modifer::amplitude: " << baseline - input_p.data.GetMinimum() << '\n';
        return {baseline - input_p.data.GetMinimum()}; 
    };
}; 

struct baseline_finder {
    using output_t = baseline;
    output_t operator()( linked_waveform const& input_p ) const {
        double baseline{0};
        for( auto i{0}; i < 64 ; ++i ) { baseline += input_p.data.GetBinContent( i +1 ); }
        baseline /= 64;
        return {baseline};
    }
};


struct cfd_calculator {
    using output_t = cfd_time;
    output_t operator()(linked_waveform const& input_p) const {
        double fraction = 0.4;
        std::size_t delay = 15;

        double baseline{0};
        for( auto i{0}; i < 64 ; ++i ) { baseline += input_p.data.GetBinContent( i +1 ); }
        baseline /= 64;
//        std::cout << "baseline: " << baseline << '\n'; 

        double current{0}, last{0};
//        std::cout << "cfd:\n";
        auto limit = input_p.data.GetMinimumBin();
        for(auto i{ limit }; i > 0 ; --i) {
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
    output_t operator()(linked_waveform const& input_p) const {

        double baseline{0};
        for( auto i{0}; i < 64 ; ++i ) { baseline += input_p.data.GetBinContent( i +1 ); }
        baseline /= 64;

        double waveform_length = input_p.data.GetBinCenter(input_p.data.GetXaxis()->GetLast());
//        std::cout << "rectangle_integral: " << baseline*waveform_length << '\n';
//        std::cout << "signal_integral: " << input_p.data.Integral("width") << '\n';
//        std::cout << "waveform_length: " << waveform_length << '\n';
        charge result{baseline*waveform_length - input_p.data.Integral("width")};  
//        return result.charge > 0 ? result : charge{0};
        return result;
    }
};

struct rise_time_calculator {
    using output_t = rise_time;
    output_t operator()( linked_waveform const& input_p ) const {
        double baseline{0}; 
        for( auto i{0}; i < 64 ; ++i ) { baseline += input_p.data.GetBinContent( i +1 ); }
        baseline /= 64;

        double amplitude = input_p.data.GetMinimum() - baseline;
        double low_cutoff = baseline + amplitude * 0.1 ;
        double high_cutoff = baseline + amplitude * 0.9 ; 

        double low_time{0}, high_time{0};
        for( auto i{ input_p.data.GetMinimumBin() - 1 } ; i > 0 ; --i ){
            if( high_time == 0 && input_p.data.GetBinContent(i) > high_cutoff ){ high_time = input_p.data.GetBinCenter(i); }                
            if( low_time == 0 && input_p.data.GetBinContent(i) > low_cutoff ){ low_time = input_p.data.GetBinCenter(i); break; }                
        }

        return { high_time - low_time }; 
    };
}; 

struct mean_finder {
    using output_t = mean;
    output_t operator()( linked_waveform const& input_p ) const {
        double result{0};
        for( auto i{1}; i < input_p.data.GetNbinsX()+1; ++i ){ result += input_p.data.GetBinContent( i ); }
        result /= input_p.data.GetNbinsX();
        return {result};
    }
};

struct sigma_finder {
    using output_t = sigma;
    output_t operator()( linked_waveform const& input_p ) const {
        double mean{0};
        double result{0};
        for( auto i{1}; i < input_p.data.GetNbinsX()+1; ++i ){ mean += input_p.data.GetBinContent( i ); }
        mean /= input_p.data.GetNbinsX();
        for( auto i{1}; i < input_p.data.GetNbinsX()+1; ++i ){ result += pow( input_p.data.GetBinContent( i ) - mean, 2); }
        result = sqrt( result/input_p.data.GetNbinsX() );
        return {result};
    }
};


} //namespace sf_g
#endif
