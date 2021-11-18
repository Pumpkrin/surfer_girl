#include <fstream>
#include <iostream>
#include <regex>
#include <string> 
        
#include "fit_tools.hpp"

#include "TGraphErrors.h"
#include "TFile.h"
#include "TFitResult.h"

int main( int argc, char* argv[]) {
    std::vector<std::string> file_c;
    auto get_part_l = []( std::string const& input_p , std::regex regex_p ){
                            std::smatch result;
                            std::regex_search( input_p, result, regex_p);
                            return result[0].str();
                      };
    auto get_result_path_l = [&get_part_l]( std::string const& input_p ){return input_p.substr(0, input_p.find( get_part_l( input_p, std::regex{"[0-9]+HV"} ) ) );};

    bool is_from_sources = true;
    std::string result_path;
    std::string curve_name;
    std::string fit_file ;
    std::string formulae;
    TF1 fit_function;
    for( auto i{1}; i < argc ; ++i ) { 
        std::string argument = argv[i];
        if( argument == std::string{"-func"} ){ 
            formulae = deduce_formulae( argv[++i] );
            fit_function = TF1{"f", formulae.c_str(), 0, 24}; 
        }
        if( argument == std::string{"-beam"} ){ is_from_sources = false; } 
        if( argument == std::string{"-folder"} ){
            result_path = get_result_path_l( argv[++i] );
            curve_name = "curve_" + get_part_l( argv[i], std::regex{"[0-9]+HV/\\w+"});
            curve_name.replace( curve_name.find( "/" ), 1, "_");
//            curve_name += "_" + formulae;
            std::string file_path = argv[i];
            fit_file = file_path.substr( 0, file_path.find_last_of( "/" ) ) + "_fit_parameters.txt";
            std::cout << fit_file << '\n';
            for( auto j{i}; j < argc ; ++j) {
                file_path = argv[j];
                auto file = file_path.substr( file_path.find_last_of( "/" ) + 1 );  
                std::cout << "currently extracting: " << file << '\n'; 
                file_c.push_back( file_path ) ;
            }
        }
    }  

    auto * g_h = new TGraphErrors{static_cast<int>(file_c.size())};
    g_h->SetName(curve_name.c_str());
    std::cout << "point_count: " << g_h->GetN() << '\n';

    auto counter{0};
    for( auto const& file : file_c ){
        std::cout << file << std::endl;
        auto energy = is_from_sources ? 
                extract_energy<source_tag>( file ) : 
                retrieve_value( 
                        result_path + "energy.map", 
                        [beam_energy = extract_energy<beam_tag>(file)](std::string const& buffer_p){ return std::stod(buffer_p) == beam_energy; },
                        [](std::string const& buffer_p){ 
                             value_and_error result;
                             result.value = find_value(buffer_p, ':') ;
                             result.error = find_value(buffer_p, '-') ; 
                             return result; 
                        }
                );
        std::cout << "energy: " << energy.value << " +/- " << energy.error << std::endl;

        auto mean = retrieve_value(
                    file,
                    [](std::string const& buffer_p){ return line_found(buffer_p, "Mean"); },
                    [](std::string const& buffer_p){ 
                        value_and_error result;
                        result.value = find_value(buffer_p, '='); 
                        result.error = find_value(buffer_p, '-');
                        return result; 
                    }
        );
        std::cout << "fit_result: " << mean.value << " -- " << mean.error << '\n';
//        g_h->SetPoint( counter , energy.value, mean.value );
        g_h->SetPoint( counter , mean.value, energy.value );
//        g_h->SetPointError( counter++ , energy.error, result.sigma.value );
//        g_h->SetPointError( counter++ , energy.error, mean.error );
        g_h->SetPointError( counter++ , mean.error, energy.error );
    }
       
//    TF1 f{"f", "[0]+x*[1]", 0, energy_c.back() * 1.2};
    if( file_c.size() > 2 ) { 
        std::cout << "writing fit results in: " << fit_file << '\n';
        auto fit = g_h->Fit(&fit_function, "S");
        std::ofstream output{ fit_file };
        auto* old_cout_buffer = std::cout.rdbuf();
        std::cout.rdbuf( output.rdbuf() );
        fit->Print();
        std::cout.rdbuf( old_cout_buffer );
    }
    std::string output = result_path + "results.root";
    std::cout << "storing curve in: " << output << '\n';
    TFile file( output.c_str() , "UPDATE" );
    g_h->Write();
}
