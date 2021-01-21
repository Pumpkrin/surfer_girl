#include "policy.hpp"
#include "reader.hpp"
#include "writer.hpp"

#include <iostream>

int main( int argc, char* argv[] ) {

    std::string input_file, output_file{"output.root"};
    char format;

    if(argc < 5){ std::cerr << "surfer_girl should be called the following way: ./surfer_girl -in input_file.bin -format m(easurement)/w(aveform) [-out output_file.root] \n"; return 1; }
    else{     
        for(auto i{0}; i < argc; ++i) {
            if( std::string( argv[i]) == "-in" ){ input_file = std::string{ argv[++i] }; }
            if( std::string( argv[i] ) == "-format" ) { format = argv[++i][0]; } //no check on given value but not really necessary: switch will do the trick later on
            if( std::string( argv[i] ) == "-out") { output_file = std::string{ argv[++i]}; }
        }
    }

    std::cout << input_file << "_" << format << " -> " << output_file << '\n';

    switch(format){
        case 'w': {
           auto r = reader<waveform>{ input_file }; 
           r.metadata.channel_count = 3;
           writer<waveform> w{ r.metadata.channel_count,  output_file };
           while( !r.end_is_reached() ){
               auto data = r.read_event();
               std::cout << data[0].channel_id << " -- " << data[0].event_id << " -- " << data[0].fcr << " -- " << data[0].baseline << " -- " << data[0].amplitude << " -- " << data[0].charge << " -- " << data[0].leading_edge << " -- " << data[0].trailing_edge  << " -- " << data[0].rate_counter <<  '\n';
               for( auto i{0}; i < 1024; ++i) {
                   std::cout << data[0].sample_c[i] << " ";
               }
               std::cout << "\n";
               std::cout << data[1].channel_id << " -- " << data[1].event_id << " -- " << data[1].fcr << " -- " << data[1].baseline << " -- " << data[1].amplitude << " -- " << data[1].charge << " -- " << data[1].leading_edge << " -- " << data[1].trailing_edge  << " -- " << data[1].rate_counter <<  '\n';
               for( auto i{0}; i < 1024; ++i) {
                   std::cout << data[1].sample_c[i] << " ";
               }
               std::cout << "\n";
               std::cout << data[2].channel_id << " -- " << data[2].event_id << " -- " << data[2].fcr << " -- " << data[2].baseline << " -- " << data[2].amplitude << " -- " << data[2].charge << " -- " << data[2].leading_edge << " -- " << data[2].trailing_edge  << " -- " << data[2].rate_counter <<  '\n';
               for( auto i{0}; i < 1024; ++i) {
                   std::cout << data[2].sample_c[i] << " ";
               }
               w.write_data( std::move(data) ); 
           }
           w.save_data();
           break;
                  }
        case 'm': {
           auto r = reader<measurement>{ input_file }; 
           r.metadata.channel_count = 1;
           while( !r.end_is_reached() ){
               auto data = r.read_event();
           }
           break;
                  }
        default: {
           std::cerr << "The format given is unknown\n";
           return 1;
                 }
    }

}

