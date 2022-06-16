#ifndef FORMULAE_HPP
#define FORMULAE_HPP
#include <utility>
#include <fstream>
#include <array>

#include "fit_tools.hpp"

namespace sf_g{
template<class T, class Enabler = void> struct is_member_of_formulae_impl : std::false_type{};
template<class T> struct is_member_of_formulae_impl< T, decltype( adl_is_member_of_formulae( std::declval<T>() ) )> : std::true_type{};
template<class T> struct is_member_of_formulae : is_member_of_formulae_impl< T, void>{};

namespace formulae{
template<class T>
constexpr void adl_is_member_of_formulae(T &&);

template<std::size_t Order>
struct polynomial{
    static constexpr std::size_t order = Order+1;
    explicit constexpr polynomial( std::array< double, order>&& parameter_pc ) : parameter_c{std::move(parameter_pc)} {}
    constexpr double operator()(double value_p) const {
        double result{0};
        for( auto i{0}; i < parameter_c.size() ; ++i ){result += pow(value_p, i) * parameter_c[i];}
        return result;
    }
    std::array<double, order> parameter_c;
};

struct gamma_resolution{
    static constexpr std::size_t order = 2;
    explicit constexpr gamma_resolution( std::array< double, order>&& parameter_pc ) : parameter_c{std::move(parameter_pc)} {}
    double operator()(double value_p) const {
        return parameter_c[0] + parameter_c[1]/sqrt(value_p); 
    }
    std::array<double, order> parameter_c;
};

template<std::size_t Order>
std::array<double, Order> retrieve_parameters( std::string const& file_p ){
    std::array<double, Order> result_c;
    std::ifstream stream{ file_p };
    std::string buffer;
    std::size_t counter{0};
    if( !stream.is_open() ){  std::cerr << "impossible to open file: " << file_p << '\n'; return result_c; }
    while( std::getline(stream, buffer) ){
        if( line_found(buffer, "p") ){result_c[counter++]=find_value( buffer, '=' ); }
    }
    return result_c;
}  
} //namespace formulae
} // namespace formulae
#endif
