#ifndef READER_HPP
#define READER_HPP

#include "data_structure.hpp"
#include "data_stream.hpp"

//#include <fstream>
#include <regex>
#include <string>
#include <vector>


namespace sf_g{ 

template<class T>
struct mapped_element {
    T& value;
    constexpr void read(std::ifstream& stream_p){ stream_p.read( reinterpret_cast<char*>(&value), sizeof(T)); }
};

template<class T, std::size_t I>
struct mapped_element< std::array<T, I> > {
    std::array<T,I>& value;
    constexpr void read(std::ifstream& stream_p){ stream_p.read( reinterpret_cast<char*>(value.data()), I*sizeof(T)); }
};

template<class ... Ts>
struct mapper {
    using tuple_t = std::tuple<mapped_element<Ts>...>;

    mapper(mapped_element<Ts>&&... ts_p) : tuple_m{std::move(ts_p)...} {}
    constexpr void read(std::ifstream& stream_p){ read_impl<0>( stream_p ); }


private:
    template< std::size_t Index,
              typename std::enable_if_t< (Index < std::tuple_size<tuple_t>::value), std::nullptr_t > = nullptr >
    constexpr void read_impl(std::ifstream& stream_p)
    {
        std::get<Index>(tuple_m).read( stream_p );
        read_impl<Index+1>( stream_p );
    }
    
    template< std::size_t Index,
              typename std::enable_if_t< (Index >= std::tuple_size<tuple_t>::value), std::nullptr_t > = nullptr >
    constexpr void read_impl(std::ifstream& /*stream_p*/) {}

   
    private:
    tuple_t tuple_m;
};

template<class ... Ts>
auto make_mapper(Ts&... ts_p){
    return mapper<Ts...>{ mapped_element<Ts>{ts_p}... };
}

template<class ... Ts>
constexpr void read_into(std::ifstream& stream_p, mapper<Ts...>&& mapper_p){
   mapper_p.read( stream_p ); 
}


template<class I, class O> struct reader{};

template<>
struct reader< std::ifstream, raw_waveform > {
    using input_t = std::ifstream;
    using output_t = raw_waveform;
    
    reader( int size_p) : size_m{ static_cast<std::size_t>(size_p)} {}
    std::vector<output_t> operator()(data_input<input_t>& input_p) const {
        std::vector<output_t> result_c{ size_m };
        for( auto& result: result_c ){
            read_into( 
                input_p.input,
                make_mapper( 
                     result.channel_id,
                     result.event_id,
                     result.fcr,
                     result.baseline,
                     result.amplitude,
                     result.charge,
                     result.leading_edge,
                     result.trailing_edge,
                     result.rate_counter,
                     result.sample_c
                )
            );
        }
        return result_c;                   
    }

    private:
    std::size_t const size_m;
};

//you're not using a thing from data_source here -> why enforce it ?
//how will a root file be handled 
template<>
struct reader< std::ifstream, event_data> {
    using input_t = std::ifstream;
    using output_t = event_data;
    
    output_t operator()(data_input<input_t>& input_p) const {
        output_t result;
        read_into( 
            input_p.input,
            make_mapper( 
                result.event_id,
                result.epoch_time,
                result.date.year, 
                result.date.month,
                result.date.day,
                result.time.hour,
                result.time.minute,
                result.time.second,
                result.time.millisecond,
                result.tdc,
                result.corrected_tdc,
                result.channel_count
            )
        );
        return result;                   
    }
};

template<>
struct reader< std::ifstream, metadata> {
    using input_t = std::ifstream;
    using output_t = metadata;
    
    output_t operator()(data_input<input_t>& input_p) const {
        output_t result;
        std::string temp;
        std::getline(input_p.input, temp);
        std::getline(input_p.input, temp);
        std::getline(input_p.input, temp);
        std::getline(input_p.input, temp);
        std::size_t metadata_offset = temp.find( ':' );
        result.channel_count = std::stoi( temp.substr( metadata_offset +1 ) );
        metadata_offset = temp.find( ':', metadata_offset + 1);
        result.sampling_period = std::stod( temp.substr( metadata_offset +1) );
        return result;                   
    }
};

template<>
struct reader< TTree, linked_waveform> {
    using input_t = TTree;
    using output_t = linked_waveform;
    
    reader(data_input<input_t>& input_p) : 
        channel_count_m{ retrieve_channel_count(input_p) },
        channel_pairing_mc{ retrieve_channel_pairing( input_p ) },
        waveform_mc{channel_count_m},
        indirector_mc{channel_count_m} 
    {
        for( auto i{0} ; i < channel_count_m ; ++i ) { 
            waveform_mc[i] = std::make_unique<waveform>( waveform{} );
            indirector_mc[i] = waveform_mc[i].get();
            std::string name = "channel_" + std::to_string(channel_pairing_mc[i]) + ".";
            input_p.input_h->SetBranchAddress( name.c_str(), &indirector_mc[i] ); 
        }
    }

    std::vector< output_t > operator()( data_input<input_t>& input_p ) {
        std::vector< output_t > output_c{channel_count_m}; 
        input_p.load_entry();
//        std::cout << "reader::channel_number: " << waveform_mc.size() << " / " << channel_count_m << '\n';
        for( auto i{0} ; i < channel_count_m ; ++i ) { 
            output_c[i].data = waveform_mc[i]->data; 
//            std::cout << "reader::amplitude: " << waveform_mc[i]->data.GetMinimum() << " / " <<  output_c[i].data.GetMinimum() << '\n';
            output_c[i].channel_number =  channel_pairing_mc[i];
        }
        return output_c;
    }

    constexpr std::size_t channel_count() const { return channel_count_m; }
private:
    std::size_t retrieve_channel_count(data_input<input_t>& input_p) const {
        std::cout << "branches: " << input_p.input_h->GetListOfBranches()->GetEntries() << '\n';
        return input_p.input_h->GetListOfBranches()->GetEntries();
    }
    std::vector<std::size_t> retrieve_channel_pairing(data_input<input_t>& input_p) const {
        std::vector<std::size_t> result_c;
        auto * l_h = input_p.input_h->GetListOfBranches();
        result_c.reserve( channel_count_m );
        for( auto * b_h : *l_h){ 
            std::string name{ b_h->GetName() };
            std::smatch result;
            std::regex_search( name, result, std::regex{"[0-9]"});
            result_c.push_back( static_cast<std::size_t>( std::stoi(result[0].str() ) ) );
        }
        return result_c;
    }

private:
    std::size_t const channel_count_m;
    std::vector< std::size_t > channel_pairing_mc;
    std::vector< std::unique_ptr<waveform> > waveform_mc;
    std::vector< waveform* > indirector_mc;
};

} // namespace sf_g

#endif
