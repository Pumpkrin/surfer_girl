#ifndef SPECTRA_HPP
#define SPECTRA_HPP

#include "formulae.hpp"
#include "data_structure.hpp"

#include <fstream>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom3.h"

template<class C, class T> TH1D* extract(std::string);
template<class C, class T> double get_value(C*);
template<class T> class TD;

namespace sf_g{
template<class T, class Enabler = void> struct is_member_of_spectra_tag : std::false_type{};
template<class T> struct is_member_of_spectra_tag< T, decltype( adl_is_member_of_spectra_tag( std::declval<T>() ) )> : std::true_type{};

namespace spectra_tag{    
template<class T>
constexpr void adl_is_member_of_spectra_tag(T &&);
template<class C, class T> struct experiment{};
struct initial_energy{ 
    double operator()(gamma_response const* input_ph) const { return input_ph->gamma_energy;}
    double operator()(response const* input_ph) const { return input_ph->gamma_energy;}
};
struct deposited_energy{
    double operator()(gamma_response const* input_ph) const { return input_ph->deposited_energy;}
    double operator()(response const* input_ph) const { return input_ph->deposited_energy;}
};
struct total_deposited_energy{
    double operator()(response const* input_ph) const { return input_ph->total_deposited_energy;}
};
template<class D, class T> struct simulation{ using data_t = D;};

struct maximum{
    double operator()( TH1D const& histogram_p ) const { return histogram_p.GetMaximum(); }
};
struct integral{
    double operator()( TH1D const& histogram_p ) const { return histogram_p.Integral(); }
};

struct gold{
    struct parameter_t{
        std::size_t iteration;
    };
    gold(std::string const& file_p, parameter_t&& parameter_p) : file_m{file_p}, parameter_m{std::move(parameter_p)} {} 
    std::vector<std::size_t> apply(TH1D&& spectrum_p, std::vector<std::size_t> const& initial_pc) const {
        TFile input{file_m.c_str() };
        auto * response_matrix_transpose_h = input.Get<TH2D>("response_matrix_transpose");
        auto * toeplitz_matrix_h = input.Get<TH2D>("toeplitz_matrix");
        
       //TODO: cross-check dimensionality between input and response matrix  
        std::size_t const size = spectrum_p.GetNbinsX();

        std::vector<std::size_t> result_c{initial_pc};

//        auto integral{ spectrum_p.Integral() };
        std::vector<std::size_t> measured_c(size);
        for( auto row_index{0}; row_index < size ; ++row_index){
            for( auto k{1}; k < size+1 ; ++k){
                measured_c[row_index] += response_matrix_transpose_h->GetBinContent(k,row_index+1)*spectrum_p.GetBinContent(k);
            }
        }

        for( auto iteration{0}; iteration < parameter_m.iteration ; ++iteration ){
            auto last_result_c{result_c};
//            auto deconvoluted_integral{0.};
            for( auto k{0}; k < size ; ++k){
                auto denominator{0.};
                for( auto l{0}; l < size ; ++l ){
                    denominator += toeplitz_matrix_h->GetBinContent(l+1, k+1) * last_result_c[l]; 
                }
                result_c[k] = (measured_c[k] != 0) ? measured_c[k]/denominator * last_result_c[k] : 0; 
//                deconvoluted_integral += result_c[k]; 
            }
//            for( auto k{0}; k < size ; ++k){ result_c[k] = result_c[k] * integral/deconvoluted_integral;}
        }
        return result_c;
    }

    private:
    std::string file_m;
    parameter_t parameter_m;
};

struct richardson_lucy{
    struct parameter_t{
        std::size_t iteration;
    };
    richardson_lucy(std::string const& file_p, parameter_t&& parameter_p) : file_m{file_p}, parameter_m{std::move(parameter_p)} {} 
    std::vector<std::size_t> apply(TH1D&& spectrum_p, std::vector< std::size_t> const& initial_pc) const {
        TFile input{file_m.c_str() };
        auto * response_matrix_h = input.Get<TH2D>("response_matrix");
        
       //TODO: cross-check dimensionality between input and response matrix  
        std::size_t const size = spectrum_p.GetNbinsX();

        std::vector<std::size_t> result_c{ initial_pc };

        std::vector<std::size_t> measured_c(size);
        for( auto row_index{0}; row_index < size ; ++row_index){
            measured_c[row_index] = spectrum_p.GetBinContent(row_index+1);
        }

        for( auto iteration{0}; iteration < parameter_m.iteration ; ++iteration ){
            auto last_result_c{result_c};
            for( auto k{0}; k < size ; ++k){
                auto factor{0.};
                for( auto l{0}; l < size ; ++l ){
                    auto denominator{0.};
                    for( auto m{0}; m < size ; ++m ){denominator += response_matrix_h->GetBinContent(m+1, l+1) * last_result_c[m];} 
                    factor += measured_c[l] != 0 ? response_matrix_h->GetBinContent(k+1, l+1) * measured_c[l] / denominator : 0;    
                }
                result_c[k] = last_result_c[k] != 0 ? last_result_c[k] * factor : 0; 
            }
//            for( auto k{0}; k < size ; ++k){ result_c[k] = result_c[k] * integral/deconvoluted_integral;}
        }
        return result_c;
    }

    private:
    std::string file_m;
    parameter_t parameter_m;
};

struct map{
    struct parameter_t{
        std::size_t iteration;
    };
    map(std::string const& file_p, parameter_t&& parameter_p) : file_m{file_p}, parameter_m{std::move(parameter_p)} {} 
    std::vector<std::size_t> apply(TH1D&& spectrum_p, std::vector< std::size_t> const& initial_pc) const {
        TFile input{file_m.c_str() };
        auto * response_matrix_h = input.Get<TH2D>("response_matrix");
        
       //TODO: cross-check dimensionality between input and response matrix  
        std::size_t const size = spectrum_p.GetNbinsX();

        std::vector<std::size_t> result_c{ initial_pc };

        std::vector<std::size_t> measured_c(size);
        for( auto row_index{0}; row_index < size ; ++row_index){
            measured_c[row_index] = spectrum_p.GetBinContent(row_index+1);
//            std::cout << measured_c[row_index] << " ";
        }
//        std::cout << '\n';

        for( auto iteration{0}; iteration < parameter_m.iteration ; ++iteration ){
            auto last_result_c{result_c};
            for( auto k{0}; k < size ; ++k){
                auto exponent{0.};
                for( auto l{0}; l < size ; ++l ){
                    auto denominator{0.};
                    for( auto m{0}; m < size ; ++m ){denominator += response_matrix_h->GetBinContent(m+1, l+1) * last_result_c[m];} 
//                    std::cout << denominator << " ";
//                    exponent += measured_c[l] != 0 ? response_matrix_h->GetBinContent(k+1, l+1) * ( measured_c[l] / denominator -1 ) : -1 *response_matrix_h->GetBinContent(k+1, l+1);    
                    exponent += response_matrix_h->GetBinContent(k+1, l+1) * ( measured_c[l] / denominator -1 ) ;    
                }
//                std::cout << "\n" << last_result_c[k] << ", " << exp(exponent) << " -> " <<  exponent << '\n';
                result_c[k] = last_result_c[k] * exp(exponent); 
            }
//            for( auto k{0}; k < size ; ++k){ result_c[k] = result_c[k] * integral/deconvoluted_integral;}
        }
        return result_c;
    }

    private:
    std::string file_m;
    parameter_t parameter_m;
};
template< class M >
struct boosted{
    struct parameter_t{
        std::size_t iteration;
        std::size_t repetition;
        double boost_coefficient;
    };
    using method_parameter_t = typename M::parameter_t;
    boosted( std::string const& file_p, parameter_t&& parameter_p ) : method_m{file_p, method_parameter_t{parameter_p.iteration}}, parameter_m{std::move(parameter_p)} {}
    std::vector<std::size_t> apply(TH1D&& spectrum_p, std::vector<std::size_t> initial_pc) const {
        std::size_t const size = spectrum_p.GetNbinsX();

        std::vector<std::size_t> result_c{initial_pc};
        
        auto integral{ spectrum_p.Integral() };
         
        for( auto repetition{0}; repetition < parameter_m.repetition ; ++repetition ){
            result_c = method_m.apply( std::move( spectrum_p ), result_c ); 
//            auto boosted_integral{0.};
             for( auto k{0}; k < size ; ++k){
                 result_c[k] = pow( result_c[k], parameter_m.boost_coefficient );
//                 boosted_integral += result_c[k];
             }
//             for( auto k{0}; k < size ; ++k){result_c[k] = result_c[k] * integral/boosted_integral;}
        }
        auto boosted_integral{0.};
        for( auto k{0}; k < size ; ++k){boosted_integral += result_c[k];}
        for( auto k{0}; k < size ; ++k){result_c[k] = result_c[k] * integral/boosted_integral;}
        return result_c;
    }

private:
    parameter_t const parameter_m;
    M method_m;
};
struct boosted_gold{
    struct parameter_t{
        std::size_t iteration;
        std::size_t repetition;
        double boost_coefficient;
    };

    boosted_gold( std::string const& file_p, parameter_t&& parameter_p ) : file_m{file_p}, parameter_m{std::move(parameter_p)} {}
    std::vector<std::size_t> apply(TH1D&& spectrum_p, std::vector< std::size_t > const& initial_pc) const {
        TFile input{file_m.c_str()};
        auto * response_matrix_transpose_h = input.Get<TH2D>("response_matrix_transpose");
        auto * toeplitz_matrix_h = input.Get<TH2D>("toeplitz_matrix");

        std::size_t const size = spectrum_p.GetNbinsX();

        std::vector<std::size_t> result_c{initial_pc};
        
        auto integral{ spectrum_p.Integral() };
        std::vector<std::size_t> measured_c(size);
        for( auto row_index{0}; row_index < size ; ++row_index){
            for( auto k{1}; k < size+1 ; ++k){
                measured_c[row_index] += response_matrix_transpose_h->GetBinContent(k,row_index+1)*spectrum_p.GetBinContent(k);
            }
        }
        
        for( auto repetition{0}; repetition < parameter_m.repetition ; ++repetition ){
            for( auto iteration{0}; iteration < parameter_m.iteration ; ++iteration ){
                auto last_result_c{result_c};
                for( auto k{0}; k < size ; ++k){
                    auto denominator{0.};
                    for( auto l{0}; l < size ; ++l ){
                         denominator += toeplitz_matrix_h->GetBinContent(l+1, k+1) * last_result_c[l]; 
                     }
                result_c[k] = (measured_c[k] != 0) ? measured_c[k]/denominator * last_result_c[k] : 0; 
                }
             }
            auto boosted_integral{0.};
             for( auto k{0}; k < size ; ++k){
                 result_c[k] = pow( result_c[k], parameter_m.boost_coefficient );
                 boosted_integral += result_c[k];
             }
             for( auto k{0}; k < size ; ++k){result_c[k] = result_c[k] * integral/boosted_integral;}
        }
        return result_c;
    }

private:
    std::string const file_m;
    parameter_t const parameter_m;
};

//struct modified_boosted_gold{
//    struct parameter_t{
//        std::size_t iteration;
//        std::size_t repetition;
//    };
//    modified_boosted_gold( std::string const& file_p, parameter_t&& parameter_p ) : file_m{file_p}, parameter_m{std::move(parameter_p)} {}
//    std::vector<std::size_t> apply(TH1D&& spectrum_p) const {
//        TFile input{file_m.c_str() };
//        auto * response_matrix_transpose_h = input.Get<TH2D>("response_matrix_transpose");
//        auto * toeplitz_matrix_h = input.Get<TH2D>("toeplitz_matrix");
//
//        std::size_t const size = spectrum_p.GetNbinsX();
//        std::vector<std::size_t> result_c(size, 1);
//        
//        auto integral{ spectrum_p.Integral() };
//        std::vector<std::size_t> measured_c(size);
//        for( auto row_index{0}; row_index < size ; ++row_index){
//            for( auto k{1}; k < size+1 ; ++k){
//                measured_c[row_index] += response_matrix_transpose_h->GetBinContent(k,row_index+1)*spectrum_p.GetBinContent(k);
//            }
//        }
//       
//        for( auto repetition{0}; repetition < parameter_m.repetition ; ++repetition ){
//            for( auto iteration{0}; iteration < parameter_m.iteration ; ++iteration ){
//                auto last_result_c{result_c};
//                for( auto k{0}; k < size ; ++k){
//                    auto denominator{0.};
//                    for( auto l{0}; l < size ; ++l ){
//                         denominator += toeplitz_matrix_h->GetBinContent(l+1, k+1) * last_result_c[l]; 
//                     }
//                result_c[k] = (measured_c[k] != 0) ? measured_c[k]/denominator * last_result_c[k] : 0; 
//                }
//             }
//            auto boosted_integral{0.};
//             for( auto k{0}; k < size ; ++k){
//                 result_c[k] = pow( result_c[k], 1.05 + k*.9/size);
//                 boosted_integral += result_c[k];
//             }
//             for( auto k{0}; k < size ; ++k){result_c[k] = result_c[k] * integral/boosted_integral;}
//        }
//        return result_c;
//    }
//private:
//    std::string const file_m;
//    parameter_t const parameter_m;
//};
} //namespace spectra_tag
namespace details{
template< class T > struct is_value_processor : std::false_type{};
};//namespace details

template<class T> struct scaler{
    template<class Format>
    TH1D apply( TH1D&& histogram_p ) const {
        histogram_p.Scale( 1./T{}( histogram_p ) ); 
        return std::move( histogram_p );
    }
};

struct basic_optional{
    double value;
    bool has_value;
};
template<class T> struct extractor{};
template<class C, class T> 
struct extractor< spectra_tag::experiment<C, T> >{
    using composite_t = C;
    TH1D apply( std::string const& file_p ) const { return *extract<C, T>( file_p ); }
    basic_optional apply( C* composite_ph ) const { return {get_value<C,T>( composite_ph ), true}; } 
};

template<class D, class T>
struct extractor< spectra_tag::simulation<D,T> >{
    static_assert( is_member_of_spectra_tag< T >::value, "this function should be used with a spectra_tag class : i.e. spectra_tag::*****" );
    template<class Format>
    TH1D apply( std::string const& file_p ) const { 
        TFile file{ file_p.c_str() };
        auto* tree_h = file.Get<TTree>( "data" );
        std::unique_ptr<D> input_h{ new D{} }; 
        auto * result_h = input_h.get();
        tree_h->SetBranchAddress( "gamma_response", &result_h);

        auto histogram = Format{}(); 
        histogram.SetDirectory( gROOT );

        for(auto i{0}; i < tree_h->GetEntries(); ++i){
            tree_h->GetEntry(i);
            histogram.Fill( T{}(result_h) );
        }
        return histogram;
    }
    basic_optional apply(D const* value_p) const { return {T{}(value_p), true}; } 
};


template<class Formulae> 
struct calibrator{
    static_assert( is_member_of_formulae< Formulae >::value, "this function should be used with a formulae class : i.e. formulae::*****" );
    calibrator( std::string const& file_p ) : formulae_m{ formulae::retrieve_parameters<Formulae::order>(file_p) } {}
    basic_optional apply(basic_optional const& value_po) const {return {formulae_m(value_po.value), true};} 
        
private: 
    Formulae formulae_m;
};
template<class F> struct details::is_value_processor< calibrator<F>>: std::true_type{};

struct threshold{
    threshold(double value_p) : threshold_m{value_p} {}
    basic_optional apply(basic_optional const& value_po) const {return {value_po.value, value_po.value > threshold_m };} 
private:
    double threshold_m;
};
template<> struct details::is_value_processor< threshold >: std::true_type{};

struct gate{
    gate(double lower_limit_p, double upper_limit_p) : lower_limit_m{lower_limit_p}, upper_limit_m{upper_limit_p} {}
    template<class Format>
    TH1D apply( TH1D&& spectrum_p ) const {
        for( auto i{1}; i < spectrum_p.GetNbinsX() + 1; ++i){ 
            double center = spectrum_p.GetBinCenter(i);
            if( center < lower_limit_m || upper_limit_m < center){ spectrum_p.SetBinContent(i, 0); }
        }
        return std::move(spectrum_p);
    }

private:
    double lower_limit_m;
    double upper_limit_m;
};

template<class Formulae>
struct smearer{
    static_assert( is_member_of_formulae< Formulae >::value, "this function should be used with a formulae class : i.e. formulae::*****" );
    smearer( std::string const& resolution_input_p ) : dice_h{new TRandom3{}}, formulae_m{ formulae::retrieve_parameters<Formulae::order>( resolution_input_p ) }  {}
    smearer( smearer&& ) = default;
    basic_optional apply(basic_optional const& value_po) const {return { (*this)(value_po.value), true };} 
    double operator()( double value_p ) const {
        auto resolution = formulae_m( value_p );
        return dice_h->Gaus( value_p, resolution/2.35 * value_p );        
    }    
private:
    std::unique_ptr<TRandom3> dice_h;
    Formulae formulae_m;
};
template<class Formulae> struct details::is_value_processor< smearer<Formulae> >: std::true_type{};

namespace details{
template< class ... Ts> struct any_of : std::false_type {};
template< class T> struct any_of<T> : T {};
template< class T, class ... Ts>
struct any_of<T, Ts...> : std::conditional< bool(T::value), T, any_of<Ts...>>::type {};

template<class T, class U> struct is_same_meta : std::false_type {};
template< class T, class U, template<class> class TT>
struct is_same_meta< TT<T>, TT<U> > : std::true_type {};
template< class CT, class T, class CU, class U, template<class, class> class TT>
struct is_same_meta< TT<CT, T>, TT<CU, U> > : std::true_type {};
template< template<std::size_t> class TT, std::size_t T, std::size_t U>
struct is_same_meta< TT<T>, TT<U> > : std::true_type {};
template< template<std::size_t, std::size_t> class TT, std::size_t IU, std::size_t IT, std::size_t RU, std::size_t RT>
struct is_same_meta< TT<IT,RT>, TT<IU,RU> > : std::true_type {};
} //namespace details

template< class M >
struct unfolder{
    using parameter_t = typename M::parameter_t;
    static_assert( is_member_of_spectra_tag< M >::value, "this function should be used with a spectra_tag class : i.e. spectra_tag::*****" );
    unfolder( std::string const& response_matrix_file_p, parameter_t&& parameters_p ) : method_m{ response_matrix_file_p, std::move(parameters_p) } {}
    template<class Format>
    TH1D apply( TH1D&& spectrum_p ) const {
        gROOT->cd();
        TH1D result{Format{}()}; 
        std::size_t size = result.GetNbinsX();
        std::vector<std::size_t> initial_c(size, 1);
        auto unfolded_c = method_m.apply(std::move(spectrum_p), initial_c );
        for( std::size_t i{1}; i < result.GetNbinsX()+1; ++i ){ result.SetBinContent( i, unfolded_c[i-1]) ; } 
        return result;
    }
    private:
    M method_m;
};
//template<std::size_t Iteration, std::size_t BinCount>
//struct unfolder{
//    static constexpr std::size_t dimension = 512;
//    static constexpr std::size_t dimension = BinCount;
//    unfolder( std::string const& response_matrix_file_p ) { response_matrix_mc = retrieve_response_matrix( response_matrix_file_p ) ;}
//    template<class Format>
//    TH1D apply( TH1D&& spectrum_p ) const {
//        std::array<double, dimension> spectrum_c;
//        for( std::size_t i{1}; i < spectrum_p.GetNbinsX()+1; ++i ){ spectrum_c[i-1] = spectrum_p.GetBinContent( i ); } 
//        std::array<double, dimension> unfolded_c{ spectrum_c };
//        for( std::size_t i{0}; i < Iteration ; ++i ){
//            std::array<double, dimension> result_c;
//            for( std::size_t j{0} ; j < dimension ; ++j){ 
//                result_c[j] = inner_product( 
//                        unfolded_c.begin(), unfolded_c.end(),
//                        response_matrix_mc.begin() + j * dimension , 0.
//                );
//            }
//            for( std::size_t j{0} ; j < dimension ; ++j){unfolded_c[j] += spectrum_c[j] - result_c[j];} 
//        }
//        gROOT->cd();
//        TH1D result{Format{}()}; 
//        for( std::size_t i{1}; i < dimension+1; ++i ){ result.SetBinContent( i, unfolded_c[i-1]) ; } 
//        return result;
//    }
//    private:
//    std::array<double, dimension * dimension> retrieve_response_matrix( std::string const& response_matrix_file_p ) const {
//        TFile input{ response_matrix_file_p.c_str() };
//        auto * response_matrix_h = input.Get<TH2D>("response_matrix");
//        std::array<double, dimension * dimension> result_c;
//        for( std::size_t i{0} ; i < dimension*dimension ; ++i){
//            std::size_t row_index = i/dimension +1;
//            std::size_t column_index = i%dimension +1;
//            result_c[i] = response_matrix_h->GetBinContent( row_index, column_index );
//        }
//        return result_c;
//    }
//
//    private:
//    std::array<double, dimension*dimension> response_matrix_mc;
//};

template<std::size_t Iteration, std::size_t BinCount>
struct folder{
    static constexpr std::size_t dimension = BinCount;
    folder( std::string const& response_matrix_file_p ) { response_matrix_mc = retrieve_response_matrix( response_matrix_file_p ) ;}
    template<class Format>
    TH1D apply( TH1D&& spectrum_p ) const {
        std::array<double, dimension> spectrum_c;
        for( std::size_t i{1}; i < spectrum_p.GetNbinsX()+1; ++i ){ spectrum_c[i-1] = spectrum_p.GetBinContent( i ); } 
        std::array<double, dimension> folded_c{ 0 };
        std::array<double, dimension> unfolded_c{ spectrum_c };
        for( std::size_t i{0}; i < Iteration ; ++i ){
            std::cout << "unfolded: ";
            for( auto const& v : unfolded_c ){ std::cout << v << " "; } 
            std::cout << "\nmatrix:\n";
            for( std::size_t j{0} ; j < dimension ; ++j){ 
//                for( std::size_t i{0} ; i < dimension ; ++i ){
//                    std::cout << *(response_matrix_mc.begin() + j * dimension + i) << " ";
//                }
                auto result =  inner_product( 
                        unfolded_c.begin(), unfolded_c.end(),
                        response_matrix_mc.begin() + j * dimension , 0.
                );
                folded_c[j] = result;
//                std::cout << result << " ";
//                std::cout << '\n';
            }
            std::cout << '\n';
            for( std::size_t j{0} ; j < dimension ; ++j){unfolded_c[j] += spectrum_c[j] - folded_c[j];} 
        }
        gROOT->cd();
        TH1D result{Format{}()}; 
        for( std::size_t i{1}; i < dimension+1; ++i ){ result.SetBinContent( i, folded_c[i-1]) ; } 
        return result;
    }
    private:
    std::array<double, dimension * dimension> retrieve_response_matrix( std::string const& response_matrix_file_p ) const {
        TFile input{ response_matrix_file_p.c_str() };
        auto * response_matrix_h = input.Get<TH2D>("response_matrix");
        std::array<double, dimension * dimension> result_c;
        for( std::size_t i{0} ; i < dimension*dimension ; ++i){
            std::size_t row_index = i/dimension +1;
            std::size_t column_index = i%dimension +1;
            result_c[i] = response_matrix_h->GetBinContent( column_index, row_index );
            if( result_c[i] > 0){std::cout << "row_found: " << row_index-1 << '\n';}
        }
        return result_c;
    }

    private:
    std::array<double, dimension*dimension> response_matrix_mc;
};



template<class ... Ts>
struct processor_holder{
    using tuple_t = std::tuple<Ts...>;

    template< class T, class ... Us >
    constexpr processor_holder( processor_holder<Us...>&& ph_p, T&& t_p ) : processor_c{ construct_processors( std::move(ph_p), std::move(t_p), std::make_index_sequence<sizeof...(Us)>{} ) } {}   
private:
    template< class PH, class T, std::size_t ...Indices> 
    constexpr tuple_t construct_processors(PH&& ph_p, T&& t_p, std::index_sequence<Indices...>) const {
        return {std::move( std::get<Indices>(ph_p.processor_c) )..., std::move(t_p)}; 
    } 

public:
    tuple_t processor_c;
};
template<> struct processor_holder<>{};

template< class ... Ts > struct reconstruction_chain{};
template<> 
struct reconstruction_chain<>{
    template< class T>
    constexpr reconstruction_chain< extractor<T>, processor_holder<>, processor_holder<> > add_extractor() && { return { processor_holder<>{}, processor_holder<>{}}; }
};

reconstruction_chain<> make_reconstruction_chain(){ return {}; }

template< class D, class ... VPs, class ... GPs>
struct reconstruction_chain< extractor<D>, processor_holder<VPs...>, processor_holder<GPs...>>{

    reconstruction_chain(processor_holder<VPs...>&& vp_pc, processor_holder<GPs...>&& gp_pc) : vp_mc{ std::move(vp_pc)}, gp_mc{ std::move(gp_pc)} {}

    template< class T,
              typename std::enable_if_t< details::is_value_processor<T>::value && !details::any_of< details::is_same_meta<T, VPs>... >::value, std::nullptr_t> = nullptr>
    constexpr auto add() && -> reconstruction_chain< extractor<D>, processor_holder<VPs..., T >, processor_holder< GPs... > >{ return {{std::move(vp_mc), T{}},std::move(gp_mc)}; }
    template< class T,
              typename std::enable_if_t< !details::is_value_processor<T>::value && !details::any_of< details::is_same_meta<T, VPs>... >::value, std::nullptr_t> = nullptr>
    constexpr auto add() && -> reconstruction_chain< extractor<D>, processor_holder<VPs...>, processor_holder< GPs..., T > >{ return {{std::move(vp_mc)},{std::move(gp_mc), T{}}}; }
    template< class T, class ... Us,
              typename std::enable_if_t< details::is_value_processor<T>::value && !details::any_of< details::is_same_meta<T, VPs>... >::value, std::nullptr_t> = nullptr>
    constexpr auto add(Us&& ... us_p) && -> reconstruction_chain< extractor<D>, processor_holder<VPs..., T >, processor_holder< GPs... > >{ return {{std::move(vp_mc), T{std::forward<Us>(us_p)...}},std::move(gp_mc)}; }
    template< class T, class ... Us,
              typename std::enable_if_t< !details::is_value_processor<T>::value && !details::any_of< details::is_same_meta<T, VPs>... >::value, std::nullptr_t> = nullptr>
    constexpr auto add(Us&& ... us_p) && -> reconstruction_chain< extractor<D>, processor_holder<VPs...>, processor_holder< GPs..., T > >{ return {{std::move(vp_mc)},{std::move(gp_mc), T{std::forward<Us>(us_p)...}}}; }

public:
    template<std::size_t Channel, class Format, class D_ = D, 
             typename std::enable_if_t< !details::is_same_meta< D_, spectra_tag::simulation<double, double> >::value , std::nullptr_t > = nullptr>
    TH1D* apply( std::string const& file_p ){
        std::string specifier = file_p + ":channel_" + std::to_string(Channel) + ".";
        auto * result_h = new TH1D{ apply_chain_experiment<Format>( specifier, std::make_index_sequence<sizeof...(VPs)>{}, std::make_index_sequence<sizeof...(GPs)>{} ) }; 
        result_h->SetDirectory( gROOT );
        return result_h ;
    };
    template<class Format, class D_ = D, 
             typename std::enable_if_t< details::is_same_meta< D_, spectra_tag::simulation<double, double> >::value , std::nullptr_t > = nullptr>
    TH1D* apply( std::string const& file_p ){
        auto * result_h = new TH1D{ apply_chain_simulation<Format>( file_p , std::make_index_sequence<sizeof...(VPs)>{}, std::make_index_sequence<sizeof...(GPs)>{} ) }; 
        result_h->SetDirectory( gROOT );
        return result_h ;
    };
private:
    template<class Format, std::size_t ... IndicesVP, std::size_t ... IndicesGP,
             typename std::enable_if_t< sizeof...(IndicesVP) == 0 , std::nullptr_t > = nullptr>
    TH1D apply_chain_experiment( std::string const& specifier_p, std::index_sequence<IndicesVP...>, std::index_sequence<IndicesGP...> ){
        TH1D result = extractor<D>{}.apply( specifier_p );
        int expander[] = { 0, (result = std::get<IndicesGP>(gp_mc.processor_c).apply(std::move(result)), 0)...};
        return result;
    }
    template<class Format, std::size_t ... IndicesVP, std::size_t ... IndicesGP,
             typename std::enable_if_t< sizeof...(IndicesVP) != 0 , std::nullptr_t > = nullptr>
    TH1D apply_chain_experiment( std::string const& specifier_p, std::index_sequence<IndicesVP...>, std::index_sequence<IndicesGP...> ){
        auto get_part_l = []( std::string const& input_p , std::regex regex_p ){
                            std::smatch result;
                            std::regex_search( input_p, result, regex_p);
                            return result[0].str();
                                                              };
        auto file_regex = std::regex{"[^:]+"}; 
        auto branch_regex = std::regex{"[^:]+\\.$"};

        TFile file( get_part_l(specifier_p, std::move(file_regex)).c_str()  );
        TTree* tree_h = static_cast<TTree*>( file.Get("data") );
        using composite_t = typename extractor<D>::composite_t;
        std::unique_ptr< composite_t > input_h{ new composite_t{} };
        auto * composite_h = input_h.get();
        tree_h->SetBranchAddress( get_part_l( specifier_p, std::move(branch_regex) ).c_str(), &composite_h);

        gROOT->cd();
        TH1D result{ Format{}() }; //need to be defined outside of ile
        for(auto i{0}; i < tree_h->GetEntries(); ++i){
           tree_h->GetEntry(i);
           auto data_o = extractor<D>{}.apply(composite_h);
           int expander[] = {0, (data_o = std::get<IndicesVP>(vp_mc.processor_c).apply(std::move(data_o)), 0)...}; 
           if( data_o.has_value ){ result.Fill( data_o.value ); } 
        }
        int expander[] = { 0, (result = std::get<IndicesGP>(gp_mc.processor_c).template apply<Format>(std::move(result)), 0)...};
        return result;
    }
    template<class Format, std::size_t ... IndicesVP, std::size_t ... IndicesGP,
             typename std::enable_if_t< sizeof...(IndicesVP) == 0 , std::nullptr_t > = nullptr>
    TH1D apply_chain_simulation( std::string const& specifier_p, std::index_sequence<IndicesVP...>, std::index_sequence<IndicesGP...> ){
        TH1D result = extractor<D>{}.template apply<Format>( specifier_p );
        int expander[] = { 0, (result = std::get<IndicesGP>(gp_mc.processor_c).template apply<Format>(std::move(result)), 0)...};
        return result;
    }
    template<class Format, std::size_t ... IndicesVP, std::size_t ... IndicesGP,
             typename std::enable_if_t< sizeof...(IndicesVP) != 0 , std::nullptr_t > = nullptr>
    TH1D apply_chain_simulation( std::string const& file_p, std::index_sequence<IndicesVP...>, std::index_sequence<IndicesGP...> ){
        TFile file{ file_p.c_str() };
        auto* tree_h = file.Get<TTree>( "data" );
        using data_t = typename D::data_t;
        std::unique_ptr<data_t> input_h{ new data_t{} }; 
        auto * workaround_input_h = input_h.get();
        tree_h->SetBranchAddress( "gamma_response", &workaround_input_h);
        gROOT->cd();
        TH1D result{ Format{}() }; //need to be defined outside of ile
        for(auto i{0}; i < tree_h->GetEntries(); ++i){
           tree_h->GetEntry(i);
           auto data_o = extractor<D>{}.apply(workaround_input_h);
           int expander[] = {0, (data_o = std::get<IndicesVP>(vp_mc.processor_c).apply(std::move(data_o)), 0)...}; 
           if( data_o.has_value ){ result.Fill( data_o.value ); } 
        }
        int expander[] = { 0, (result = std::get<IndicesGP>(gp_mc.processor_c).template apply<Format>(std::move(result)), 0)...};
        return result;
    }
private:
    processor_holder<VPs...> vp_mc;
    processor_holder<GPs...> gp_mc;
};

template< std::size_t Lower, std::size_t Upper> struct range{};
template< std::size_t BinCount, class T> struct formatter{};
template< std::size_t BinCount, std::size_t Lower, std::size_t Upper > struct formatter< BinCount, range< Lower, Upper > >{
    TH1D operator()() const {
        return TH1D{ "h", ";energy (MeV);count", BinCount, Lower, Upper};
    }
};


}//namespace sf_g

#endif

