#include <fstream>
#include <iostream>
#include <regex>
#include <string> 

#include "fit_tools.hpp"

#include "TFile.h"
#include "TH1.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TFitResult.h"


int main( int argc, char* argv[]) {
    std::vector<std::string> file_c;
    auto get_part_l = []( std::string const& input_p , std::regex regex_p ){
                            std::smatch result;
                            std::regex_search( input_p, result, regex_p);
                            return result[0].str();
                      };
    auto get_result_path_l = [&get_part_l]( std::string const& input_p ){return input_p.substr(0, input_p.find( get_part_l( input_p, std::regex{"[0-9]+HT"} ) ) );};

    bool is_from_sources = true;
    std::string result_path;
    std::string curve_name;
    std::string fit_path ;
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
            curve_name = "energy_resolution_" + get_part_l( argv[i], std::regex{"[0-9]+HV/\\w+"});
            curve_name.replace( curve_name.find( "/" ), 1, "_");
//            curve_name += "_" + formulae;
            std::string file_path = argv[i];
            fit_path = file_path.substr( 0, file_path.find_last_of( "/" ) ) ;
            std::cout << fit_path << '\n';
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
        auto initial_energy = extract_energy<beam_tag>(file);
        auto energy = is_from_sources ? 
                extract_energy<source_tag>( file ) : 
                retrieve_value( 
                        result_path + "energy.map", 
                        [&initial_energy](std::string const& buffer_p){ return std::stod(buffer_p) == initial_energy; },
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
        auto sigma = retrieve_value(
                    file,
                    [](std::string const& buffer_p){ return line_found(buffer_p, "Sigma"); },
                    [](std::string const& buffer_p){ 
                        value_and_error result;
                        result.value = find_value(buffer_p, '='); 
                        result.error = find_value(buffer_p, '-');
                        return result; 
                    }
        );
        std::cout << "fit_result: " << mean.value << " -- " << sigma.value << '\n';
        if( !is_from_sources ){
            auto straggling = retrieve_value(
                    result_path + "energy.map",
                    [&initial_energy](std::string const& buffer_p){   return stod(buffer_p) == initial_energy; },
                    [](std::string const& buffer_p){ 
                        value_and_error result;
                        result.value = find_value(buffer_p, '|') ;
                        result.error = 0.005; //to be replaced later, but requires to run simulation in order to find error on the fit + initial values from beam
                        return result; 
                    }
            );
            std::cout << "uncalibrated_straggling: " << straggling.value << '\n';

            std::cout << "upper, lower: " << energy.value + straggling.value/2  << " - " << energy.value - straggling.value/2 << " -> ";
            auto straggling_lower_limit = apply_polynomial_calibration( fit_path + "_fit_parameters.txt", energy.value - straggling.value/2 );
            auto straggling_upper_limit = apply_polynomial_calibration( fit_path + "_fit_parameters.txt", energy.value + straggling.value/2 );
            std::cout << straggling_upper_limit << " - " << straggling_lower_limit << " -> " << straggling_upper_limit - straggling_lower_limit << '\n';
            straggling.value = straggling_upper_limit - straggling_lower_limit;

            std::cout << "straggling: " << straggling.value << '\n';
            sigma.error = sigma.value * sigma.error/sqrt( pow(sigma.value, 2) - pow(straggling.value, 2) );
            sigma.value = sqrt( pow(sigma.value, 2) - pow(straggling.value, 2) );  
            std::cout << "sigma: " << sigma.value << " - " << sigma.error << '\n';
        }
        auto resolution = compute_energy_resolution( mean, sigma );
        std::cout << "resolution: " << resolution.value << " -- " << resolution.error << '\n';
        
        g_h->SetPoint( counter , energy.value, resolution.value );
        g_h->SetPointError( counter++ , energy.error, resolution.error );
    }
       
    if( file_c.size() > 2 ) { 
        auto fit = g_h->Fit(&fit_function, "RS");
        std::ofstream output{ fit_path + "_energy_resolution_parameters.txt" };
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
