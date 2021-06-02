
#include <fstream>
#include <iostream>
#include <regex>
#include <string> 

#include "TFile.h"
#include "TH1.h"
#include "TGraphErrors.h"
#include "TF1.h"

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

struct fit_result {
    double mean;
    double sigma;
};


fit_result extract_fit_result( std::ifstream& input_p ) {
    fit_result result;
    std::string temp;
    auto find_value = [](std::string const& source_p ){
        auto offset = source_p.find( '=' );
        if( offset != std::string::npos ){
            return std::stod( source_p.substr( offset + 1 )); //stod will stop at first white space
        }
        return 0.;  
    };
    auto line_found = [](std::string const& source_p, std::string const& match_p){
        return source_p.find(match_p) != std::string::npos;
    };
    bool mean_filled{false};
    while( std::getline( input_p, temp )  ){
        if( !mean_filled && line_found( temp, "Mean") ){  result.mean = find_value(temp); mean_filled = true; continue; } 
        if( line_found( temp, "Sigma") ){ result.sigma = find_value(temp); break; } 
    }
    return result;
}

struct energy {
    double value;
    double error;
};

energy extract_energy( std::string const& file_p ) {
    auto get_part_l = []( std::string const& input_p , std::regex regex_p ){
                            std::smatch result;
                            std::regex_search( input_p, result, regex_p);
                            return result[0].str();
                                                              };
    energy result;
//    std::cout << get_part_l( file_p, std::regex{"[^_]+(?=pm)"} ) << " - " << get_part_l( file_p, std::regex{"[^pm]+(?=\\.cal)"} ) << std::endl; 
    result.value = std::stod( get_part_l( file_p, std::regex{"[^_]+(?=pm)"} ) );
    result.error = std::stod( get_part_l( file_p, std::regex{"[^pm]+(?=\\.cal)"} ));
    return result;
}

int main( int argc, char* argv[]) {
    std::vector<std::string> beam_file_c, source_file_c;
    std::string calibration_file;
    for( auto i{1}; i < argc ; ++i ) { 
        std::string full_path = argv[i];
        calibration_file = full_path.substr( 0, full_path.find_last_of( "/" ) ) + "_calibration.root";
        auto file = full_path.substr( full_path.find_last_of( "/" ) + 1 );  
        file.front() == 's' ? source_file_c.push_back( full_path ) : beam_file_c.push_back( full_path ) ;
    }  

    auto * g_h = new TGraphErrors{static_cast<int>(beam_file_c.size() + source_file_c.size())};
    g_h->SetName("curve");
    std::cout << "point_count: " << g_h->GetN() << '\n';

    auto counter{0};
    for( auto const& file : source_file_c ){
        std::cout << file << std::endl;
        auto energy = extract_energy( file );
        std::cout << "energy: " << energy.value << " +/- " << energy.error << std::endl;

        std::ifstream input{ file };
        if( !input.is_open() ){ std::cerr << "impossible to open file: " << file << '\n'; continue; }
        auto result = extract_fit_result( input ) ;
        std::cout << "fit_result: " << result.mean << " -- " << result.sigma << '\n';
        g_h->SetPoint( counter , energy.value, result.mean );
        g_h->SetPointError( counter++ , energy.error, result.sigma );
    }
       
//    TF1 f{"f", "[0]+x*[1]", 0, energy_c.back() * 1.2};
    TF1 f{"f", "[0]+x*[1]", 0, 24 };
    g_h->Fit(&f, "R");

    for( auto const& file : beam_file_c ) {
        std::cout << file << '\n';
        auto energy = extract_energy( file );
        std::cout << "energy: " << energy.value << " +/- " << energy.error << std::endl;

        std::ifstream input{ file };
        if( !input.is_open() ){ std::cerr << "impossible to open file: " << file << '\n'; continue; }
        auto result = extract_fit_result( input ) ;
        std::cout << "fit_result: " << result.mean << " -- " << result.sigma << '\n';
        g_h->SetPoint( counter , energy.value, result.mean );
        g_h->SetPointError( counter++ , energy.error, result.sigma );
    } 

    TFile file( calibration_file.c_str(), "RECREATE" );
    g_h->Write();
}
