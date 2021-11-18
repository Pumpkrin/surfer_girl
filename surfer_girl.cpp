#include "utilities.hpp"

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

    std::cout << "processing: " << input_file << '\n';
    std::cout << "outputing_to: " << output_file << '\n';

    auto source = sf_g::data_input< std::ifstream >{ input_file };           
    auto const metadata = sf_g::reader<std::ifstream, sf_g::metadata>{}( source );
    std::cout << "metadata: " << metadata.channel_count << " - " << metadata.sampling_period << "\n";
    auto const event_reader = sf_g::reader< std::ifstream, sf_g::event_data>{};
    auto const waveform_reader = sf_g::reader< std::ifstream, sf_g::raw_waveform>{metadata.channel_count}; 
    //add event id somewhere ? -> cross check with the waveforms
    auto const m = sf_g::raw_modifier{metadata};

    sf_g::data_output< TTree > sink{ output_file };
    auto w = sf_g::raw_writer{ sink, metadata.channel_count };

    while( !source.end_is_reached() ){
//         for( auto j{event_reader( source );0} ; j < 1; ++j) {
        event_reader( source );
//        auto timing = event_reader( source );
//        std::cout << "event: ";
//        std::cout << timing.event_id << " -- ";  
//        std::cout << timing.epoch_time << " -- ";
//        std::cout << timing.time.second << " - " << timing.time.millisecond << " -- ";
//        std::cout << "tdc: " <<  timing.tdc << " - " << corrected_tdc << '\n';
//        auto data = waveform_reader( source );
//        std::cout << "channel: " << data[0].channel_id << " -- " << data[0].event_id << " -- " << data[0].fcr << " -- " << data[0].baseline << " -- " << data[0].amplitude << " -- " << data[0].charge << " -- " << data[0].leading_edge << " -- " << data[0].trailing_edge  << " -- " << data[0].rate_counter <<  '\n';
//        for( auto i{0}; i < 1024 ; ++i){
//            std::cout << data[0].sample_c[i] << " ";
//        }
//        std::cout << '\n';
//        m( std::move(data) ) | w;
        waveform_reader( source) | m | w ;
    }
}

