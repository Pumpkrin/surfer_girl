#include "utilities.hpp"

#include <string>
#include <regex> 

int main( int argc, char* argv[]) {
    std::string cut_file;

    if(argc < 3){ 
        std::cerr << "tree_grubbing is used to apply a cut to a file output from surfer_girl. It will only contain the selected waveform and can then be reprocessed by tree_flip.\nIt should be called the following way: ./tree_grubbing -cut file.cut \n"; return 1; }
    else{     
        for(auto i{0}; i < argc; ++i) {
            if( std::string( argv[i] ) == "-cut") { cut_file = std::string{ argv[++i]}; }
        }
    }

    auto get_part_l = []( std::string const& filename_p, std::regex regex_p ){
                            std::smatch result;
                            std::regex_search( filename_p, result, regex_p);
                            return result[0].str();
                                                                             };
    auto waveform_repertory = get_part_l( cut_file, std::regex{"[^/]+"} ) + "/waveform/";
    auto output_file = waveform_repertory + get_part_l( cut_file, std::regex{"[^/]+(?=\\.cut)"} ) + ".root";
    std::cout << output_file << '\n';

    auto waveform_input = waveform_repertory + get_part_l( cut_file, std::regex{"[^/]+(?=_c)"} ) + ".root";
    std::cout << waveform_input << '\n';

    sf_g::data_input<TTree> source{ waveform_input };           
    auto r = sf_g::reader<TTree, sf_g::linked_waveform>{ source }; //channel count can be deduced from number of keys in tree

    auto m = sf_g::cut_modifier{ cut_file }; 

    sf_g::data_output< TTree > sink{ output_file };
    auto w = sf_g::raw_writer{ sink, 1 };
    while( !source.end_is_reached() ){
        r(source) | m | w;        
    }
}

