#include "utilities.hpp"
#include "specifier.hpp"

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
        std::cerr << "tree_flip is used to produce measurements from waveform, and can handle several input files.\nIf those are not provided, it will try to get them from std::cin.\nIt should be called the following way: ./tree_flip [-in input_file.root:...] -mod {channel_number:modules:...}[...] -out output_file.root \n"; return 1; }
    else{     
        for(auto i{0}; i < argc; ++i) {
            if( std::string( argv[i] ) == "-in" ){ input_file_c = regex_split(std::string{ argv[++i] }, std::regex{"[^:]+"}); input_provided = true; }
            if( std::string( argv[i] ) == "-out") { output_file = std::string{ argv[++i]}; }
            if( std::string( argv[i] ) == "-mod") {  module_list_c = regex_split(std::string{ argv[++i] }, std::regex{"\\{[^\\{\\}]*\\}"}); }
        }
    }

    if(!input_provided){
        std::string input;
        std::cin >> input;
        input_file_c =  regex_split(input, std::regex{"\\w+.root"});
    }

    std::size_t const channel_count = module_list_c.size();
    std::vector<sf_g::specifier_pairing> sp_c(channel_count);
    for( auto i{0}; i < channel_count ; ++i){ 
          sp_c[i] = sf_g::make_specifier_pairing( module_list_c[i] ); 
//          printf("opcode: %x\n", opcode_c[i]);
    } 

    sf_g::data_output< TTree > sink{ output_file };
    sf_g::tree_editor te;
    te.reserve(channel_count);
    for( auto& sp : sp_c ){
        switch( sp.opcode ) {
        using namespace sf_g;
        case flag_set<amplitude_flag>{} : {
            te.add( sf_g::branch_editor< specifier<flag_set<amplitude_flag>{}> >{ sp.channel_number, sink } ); 
            break;
                                         }
        case flag_set<charge_flag>{} : {
            te.add( sf_g::branch_editor< specifier<flag_set<charge_flag>{}> >{ sp.channel_number, sink } ); 
            break;
                                         }
        case flag_set<rise_flag>{} : {
            te.add( sf_g::branch_editor< specifier<flag_set<rise_flag>{}> >{ sp.channel_number, sink } ); 
            break;
                                         }
        case flag_set<cfd_flag>{} : {
            te.add( sf_g::branch_editor< specifier<flag_set<cfd_flag>{}> >{ sp.channel_number, sink } ); 
            break;
                                         }
        case flag_set<baseline_flag>{} : {
            te.add( sf_g::branch_editor< specifier<flag_set<baseline_flag>{}> >{ sp.channel_number, sink } ); 
            break;
                                         }
        case flag_set<amplitude_flag, baseline_flag>{} : {
            te.add( sf_g::branch_editor< specifier<flag_set<amplitude_flag, baseline_flag>{}> >{ sp.channel_number, sink } ); 
            break;
                                         }
        case flag_set<amplitude_flag, charge_flag>{} : {
            te.add( sf_g::branch_editor< specifier<flag_set<amplitude_flag, charge_flag>{}> >{ sp.channel_number, sink } ); 
            break;
                                         }
        case flag_set<amplitude_flag, charge_flag, baseline_flag>{} : {
            te.add( sf_g::branch_editor< specifier<flag_set<amplitude_flag, charge_flag, baseline_flag>{}> >{ sp.channel_number, sink } ); 
            break;
                                         }
        case flag_set< baseline_flag, charge_flag, rise_flag >{} : {
            te.add( sf_g::branch_editor< specifier<flag_set<baseline_flag, charge_flag, rise_flag >{}> >{ sp.channel_number, sink } ); 
            break;
                                         }
        default: {
            std::cerr << "This configuration as not been implemented yet\n";
            break;
                 }
        }
    }
    for( auto& input_file : input_file_c ){ 
        std::cout << "processing: " << input_file << '\n'; 
        sf_g::data_input<TTree> source{ input_file };           
        auto r = sf_g::reader<TTree, sf_g::waveform>{ source }; //channel count can be deduced from number of keys in tree
        //on reader, always retrieve everything
        
        while( !source.end_is_reached() ){
            r(source) | te;
        }
    }

}
