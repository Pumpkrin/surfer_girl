#include "reader.hpp"
#include "modifier.hpp"
#include "writer.hpp"
#include "data_format.hpp"

#include <iostream>

int main( int argc, char* argv[] ) {

    std::string input_file, output_file{"output.root"};

    if(argc < 3){ std::cerr << "surfer_girl should be called the following way: ./surfer_girl -in input_file.bin [-out output_file.root] \n"; return 1; }
    else{     
        for(auto i{0}; i < argc; ++i) {
            if( std::string( argv[i] ) == "-in" ){ input_file = std::string{ argv[++i] }; }
            if( std::string( argv[i] ) == "-out") { output_file = std::string{ argv[++i]}; }
        }
    }

    std::cout << output_file ;

    auto source = sf_g::data_input< std::ifstream >{ input_file };           
    auto const metadata = sf_g::reader<std::ifstream, sf_g::metadata>{}( source );
//    std::cout << "metadata: " << metadata.channel_count << " - " << metadata.sampling_period << "\n";
    auto const event_reader = sf_g::reader< std::ifstream, sf_g::event_data>{};
    auto const waveform_reader = sf_g::reader< std::ifstream, sf_g::raw_waveform>{metadata.channel_count}; 

    auto const m = sf_g::modifier<sf_g::raw_waveform, sf_g::waveform>{metadata};

    sf_g::data_output< TTree > sink{ output_file };
    auto w = sf_g::writer<sf_g::waveform, TTree>{ sink, metadata.channel_count };

//    std::cout << std::boolalpha << "end_is_reached: " << source.end_is_reached()  << '\n';
    while( !source.end_is_reached() ){
//         for( auto j{0} ; j < 1; ++j) {
        event_reader( source );
        auto data = waveform_reader( source );
//               std::cout << data[0].channel_id << " -- " << data[0].event_id << " -- " << data[0].fcr << " -- " << data[0].baseline << " -- " << data[0].amplitude << " -- " << data[0].charge << " -- " << data[0].leading_edge << " -- " << data[0].trailing_edge  << " -- " << data[0].rate_counter <<  '\n';
//               for( auto i{0}; i < 1024; ++i) {
//                   std::cout << data[0].sample_c[i] << " ";
//               }
        auto modified_data = m( std::move(data) );
//               std::cout << "\n";
        w( std::move(modified_data) ); 
    }
}

