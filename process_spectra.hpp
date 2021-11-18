#ifndef SPECTRA_HPP
#define SPECTRA_HPP

#include "formulae.hpp"

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
};
struct deposited_energy{
    double operator()(gamma_response const* input_ph) const { return input_ph->deposited_energy;}
};
template<class T> struct simulation{};

struct maximum{
    double operator()( TH1D const& histogram_p ) const { return histogram_p.GetMaximum(); }
};
struct integral{
    double operator()( TH1D const& histogram_p ) const { return histogram_p.Integral(); }
};
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

template<class T>
struct extractor< spectra_tag::simulation<T> >{
    static_assert( is_member_of_spectra_tag< T >::value, "this function should be used with a spectra_tag class : i.e. spectra_tag::*****" );
    template<class Format>
    TH1D apply( std::string const& file_p ) const { 
        TFile file{ file_p.c_str() };
        auto* tree_h = file.Get<TTree>( "data" );
        std::unique_ptr<gamma_response> input_h{ new gamma_response{} }; 
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
    basic_optional apply(gamma_response const* value_p) const { return {T{}(value_p), true}; } 
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


template<std::size_t Iteration, std::size_t BinCount>
struct unfolder{
//    static constexpr std::size_t dimension = 512;
    static constexpr std::size_t dimension = BinCount;
    unfolder( std::string const& response_matrix_file_p ) { response_matrix_mc = retrieve_response_matrix( response_matrix_file_p ) ;}
    template<class Format>
    TH1D apply( TH1D&& spectrum_p ) const {
        std::array<double, dimension> spectrum_c;
        for( std::size_t i{1}; i < spectrum_p.GetNbinsX()+1; ++i ){ spectrum_c[i-1] = spectrum_p.GetBinContent( i ); } 
        std::array<double, dimension> unfolded_c{ spectrum_c };
        for( std::size_t i{0}; i < Iteration ; ++i ){
            std::array<double, dimension> result_c;
            for( std::size_t j{0} ; j < dimension ; ++j){ 
                result_c[j] = inner_product( 
                        unfolded_c.begin(), unfolded_c.end(),
                        response_matrix_mc.begin() + j * dimension , 0.
                );
            }
            for( std::size_t j{0} ; j < dimension ; ++j){unfolded_c[j] += spectrum_c[j] - result_c[j];} 
        }
        gROOT->cd();
        TH1D result{Format{}()}; 
        for( std::size_t i{1}; i < dimension+1; ++i ){ result.SetBinContent( i, unfolded_c[i-1]) ; } 
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
            result_c[i] = response_matrix_h->GetBinContent( row_index, column_index );
        }
        return result_c;
    }

    private:
    std::array<double, dimension*dimension> response_matrix_mc;
};

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
} //namespace details


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
    template< class T, class U,
              typename std::enable_if_t< details::is_value_processor<T>::value && !details::any_of< details::is_same_meta<T, VPs>... >::value, std::nullptr_t> = nullptr>
    constexpr auto add(U const& u_p) && -> reconstruction_chain< extractor<D>, processor_holder<VPs..., T >, processor_holder< GPs... > >{ return {{std::move(vp_mc), T{u_p}},std::move(gp_mc)}; }
    template< class T, class U,
              typename std::enable_if_t< !details::is_value_processor<T>::value && !details::any_of< details::is_same_meta<T, VPs>... >::value, std::nullptr_t> = nullptr>
    constexpr auto add(U const& u_p) && -> reconstruction_chain< extractor<D>, processor_holder<VPs...>, processor_holder< GPs..., T > >{ return {{std::move(vp_mc)},{std::move(gp_mc), T{u_p}}}; }

public:
    template<std::size_t Channel, class Format, class D_ = D, 
             typename std::enable_if_t< details::is_same_meta< D_, spectra_tag::simulation<double> >::value , std::nullptr_t > = nullptr>
    TH1D* apply( std::string const& file_p ){
        std::string specifier = file_p + ":channel_" + std::to_string(Channel) + ".";
        auto * result_h = new TH1D{ apply_chain_experiment<Format>( specifier, std::make_index_sequence<sizeof...(VPs)>{}, std::make_index_sequence<sizeof...(GPs)>{} ) }; 
        result_h->SetDirectory( gROOT );
        return result_h ;
    };
    template<class Format, class D_ = D, 
             typename std::enable_if_t< details::is_same_meta< D_, spectra_tag::simulation<double> >::value , std::nullptr_t > = nullptr>
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
        int expander[] = { 0, (result = std::get<IndicesGP>(gp_mc.processor_c).apply(std::move(result)), 0)...};
        return result;
    }
    template<class Format, std::size_t ... IndicesVP, std::size_t ... IndicesGP,
             typename std::enable_if_t< sizeof...(IndicesVP) != 0 , std::nullptr_t > = nullptr>
    TH1D apply_chain_simulation( std::string const& file_p, std::index_sequence<IndicesVP...>, std::index_sequence<IndicesGP...> ){
        TFile file{ file_p.c_str() };
        auto* tree_h = file.Get<TTree>( "data" );
        std::unique_ptr<gamma_response> input_h{ new gamma_response{} }; 
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

