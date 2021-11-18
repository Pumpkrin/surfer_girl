#ifndef FIT_TOOLS_HPP
#define FIT_TOOLS_HPP
#include <fstream>
#include <iostream>
#include <regex>
#include <string> 

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

struct value_and_error {
    double value;
    double error;
};


std::string get_part( std::string const& input_p , std::regex regex_p ){
    std::smatch result;
    std::regex_search( input_p, result, regex_p);
    return result[0].str();
}

double find_value(std::string const& source_p, const char delimiter_p  ){
    auto offset = source_p.find( delimiter_p );
    if( offset != std::string::npos ){
        return std::stod( source_p.substr( offset + 1 )); //stod will stop at first white space
    }
    return 0.;  
};

auto line_found(std::string const& source_p, std::string const& match_p){
    return source_p.find(match_p) != std::string::npos;
};

template<class P, class E>
value_and_error retrieve_value(std::string const& file_input_p, P predicate_p, E extractor_p ) {
    value_and_error result;
    std::string temp;
    std::ifstream input( file_input_p);
    if( !input.is_open() ){  std::cerr << "impossible to open file: " << file_input_p << '\n'; return result; }
    while( std::getline( input, temp )  ){
//        std::cout << temp << '\n';
        if( predicate_p( temp ) ){
            result = extractor_p(temp);
            break;
        }
    }
    return result;
}
    

struct source_tag{};
struct beam_tag{};

template<class T, typename std::enable_if_t< std::is_same<T, beam_tag>::value, std::nullptr_t> = nullptr> 
double extract_energy( std::string const& file_p ) {
    double result{};
    result = std::stod(get_part( file_p, std::regex{"[0-9]+\\.[0-9]+"} ));
    return result;
}

template<class T, typename std::enable_if_t< std::is_same<T, source_tag>::value, std::nullptr_t> = nullptr> 
value_and_error extract_energy( std::string const& file_p ) {
    value_and_error result;
    result.value = std::stod( get_part( file_p, std::regex{"[^/]+(?=pm)"} ) );
    result.error = std::stod( get_part( file_p, std::regex{"[^pm]+(?=\\.cal)"} ));
    return result;
}


double apply_polynomial_calibration( std::string const& file_p, double const& input_p ){
    double result{};
    std::vector<double> parameter_c;
    std::ifstream stream{ file_p };
    std::string buffer;
    if( !stream.is_open() ){  std::cerr << "impossible to open file: " << file_p << '\n'; return result; }
    while( std::getline(stream, buffer) ){
        if( line_found(buffer, "p") ){ parameter_c.push_back( find_value( buffer, '=' ) ); }
    }  
    for( auto i{0}; i < parameter_c.size() ; ++i ){
//        std::cout << "\ncurrent: "<<  parameter_c[i] << " " << pow(input_p, i) << '\n';
        result += pow(input_p, i) * parameter_c[i];
    }  
//    std::cout << result << '\n';
    return result;
}

                

value_and_error compute_energy_resolution( value_and_error const& mean_p, value_and_error const& sigma_p){
    value_and_error result;
    result.value = sigma_p.value * 2.35/mean_p.value;
    result.error = sqrt( pow(2.35/mean_p.value * sigma_p.error, 2) + 
                              pow( 2.35*sigma_p.value/pow(mean_p.value, 2)*mean_p.error, 2) );
    return result;
} 

std::string deduce_formulae( std::string const& func_p ){
    std::string formulae;
    switch( func_p[0] ){
    case 'p': {
        std::size_t order = std::stoi( std::string{func_p[1]} );
        formulae = "pol" + std::to_string(order) + "(0)";
        break;
    }
    case 'e': {
        formulae = "[0] + exp(-[1]*x) ";
        break;
              }
    case 's': {
        formulae = "[0]+[1]/sqrt(x)";
        break;
              }
    case 'b': {
//        formulae = "[0]+[1]/sqrt(x)+[2]*x";
//        formulae = "[0]+[1]/([2]*sqrt(x)+x)";
//        formulae = "[0]+[1]/(sqrt(x)+[2]*x)";
        formulae = "[0]+[1]/sqrt(x)+[2]/x";
//        formulae = "[0]+[1]/pow(x, 1./4)";
//        formulae = "sqrt([0]/x+[1])";
        break;
              }
    default: {std::cerr << "Unknown function provided \n"; break;}
    }
    return formulae;
}
#endif
