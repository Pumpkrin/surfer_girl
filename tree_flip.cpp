#include "utilities.hpp"
#include "flag_set.hpp"

#include <regex>
#include <iostream>

template<class T> struct td;

std::vector< std::string > regex_split( std::string const & text_p, std::regex regex_p ) {
   std::vector<std::string> result_c;
   result_c.reserve(10);
        
   std::sregex_iterator match_i{ text_p.begin(), text_p.end(), regex_p };
   auto end_i = std::sregex_iterator{};
        
   for(auto iterator = match_i; iterator != end_i ; ++iterator){
       result_c.push_back( iterator->str() );
   }
        
   return result_c;
}

int main( int argc, char* argv[] ) {

    std::string output_file;
    std::vector<std::string> input_file_c;
    bool input_provided{false};
    std::vector<std::string> module_list_c;

    if(argc < 5){ 
        std::cerr << "tree_flip is used to produce measurements from waveform, and can handle several input files.\nIf those are not provided, it will try to get them from std::cin.\nIt should be called the following way: ./tree_flip [-in input_file.root:...] -mod {module:...}[...] -out output_file.root \n"; return 1; }
    else{     
        for(auto i{0}; i < argc; ++i) {
            if( std::string( argv[i] ) == "-in" ){ input_file_c = regex_split(std::string{ argv[++i] }, std::regex{"(\\w|\\.)+([^:]|$)"}); input_provided = true; }
            if( std::string( argv[i] ) == "-out") { output_file = std::string{ argv[++i]}; }
            if( std::string( argv[i] ) == "-mod") {  module_list_c = regex_split(std::string{ argv[++i] }, std::regex{"\\{[^\\{\\}]*\\}"}); }
        }
    }

    std::size_t const channel_count = module_list_c.size();
    if(!input_provided){
        std::string input;
        std::cin >> input;
        input_file_c =  regex_split(input, std::regex{"\\w+.root"});
    }

    std::vector<uint8_t> opcode_c(channel_count);
    for( auto i{0}; i < channel_count ; ++i){ 
          opcode_c[i] = sf_g::make_single_opcode( module_list_c[i] ); 
    } 
    auto modifier_opcode = sf_g::make_combined_opcode( opcode_c );

    sf_g::data_output< TTree > sink{ output_file };
    for( auto& input_file : input_file_c ){ 

        std::cout << "processing: " << input_file << '\n'; 
        sf_g::data_input<TTree> source{ input_file };           
        auto r = sf_g::reader<TTree, sf_g::waveform>{ channel_count }; 
        
        printf( "modifier_opcode: %x\n", modifier_opcode);
        switch( modifier_opcode ) {
            using namespace sf_g;
            case to_opcode(flag_set<amplitude_flag>{}, flag_set<cfd_flag>{}) : {
            auto const sm1 = make_sub_modifier( amplitude_finder{} );
            auto const sm2 = make_sub_modifier( cfd_calculator{}  );
            auto const m = make_multi_modifier<waveform>( std::move(sm1), std::move(sm2) ); 
            auto w = make_multi_writer< TTree >( m, sink ); 
            while( !source.end_is_reached() ){
                r(source) | m | w;
            }
            break;
                                                                          }
            case to_opcode(flag_set<cfd_flag>{}, flag_set<amplitude_flag>{}) : {
            auto const sm1 = make_sub_modifier( cfd_calculator{} );
            auto const sm2 = make_sub_modifier( amplitude_finder{} );
            auto const m = make_multi_modifier<waveform>( std::move(sm1), std::move(sm2) ); 
            auto w = make_multi_writer< TTree >( m, sink ); 
            while( !source.end_is_reached() ){
                r(source) | m | w;
            }
            break;
                                                                          }
            case to_opcode(flag_set<amplitude_flag, cfd_flag>{}) : {
            auto const sm = make_sub_modifier( amplitude_finder{}, cfd_calculator{} );
            auto const m = make_multi_modifier<waveform>( std::move(sm) ); 
            auto w = make_multi_writer< TTree >( m, sink ); 
            while( !source.end_is_reached() ){
                r(source) | m | w;
            }
            break;
                                                                          }
            case to_opcode(flag_set<amplitude_flag>{}) : {
            auto const sm = make_sub_modifier( amplitude_finder{} );
            auto const m = make_multi_modifier<waveform>( std::move(sm) ); 
            auto w = make_multi_writer< TTree >( m, sink ); 
            while( !source.end_is_reached() ){
                r(source) | m | w;
            }
            break;
                                                                          }
            case to_opcode(flag_set<cfd_flag>{}) : {
            auto const sm = make_sub_modifier( cfd_calculator{} );
            auto const m = make_multi_modifier<waveform>( std::move(sm) ); 
            auto w = make_multi_writer< TTree >( m, sink ); 
            while( !source.end_is_reached() ){
                r(source) | m | w;
            }
            break;
                                                                          }
            case to_opcode(flag_set<baseline_flag>{}) : {
            auto const sm = make_sub_modifier( baseline_finder{} );
            auto const m = make_multi_modifier<waveform>( std::move(sm) ); 
            auto w = make_multi_writer< TTree >( m, sink ); 
            while( !source.end_is_reached() ){
                r(source) | m | w;
            }
            break;
                                                                          }
            default:{
                std::cerr << "This modifier configuration has not been implemented yet\n";
                break; 
                    }
        } 
    }

}
