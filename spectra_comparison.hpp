
#ifndef SPECTRA_COMPARISON_HPP
#define SPECTRA_COMPARISON_HPP

#include "fit_tools.hpp"
#include <fstream>

#include "TFile.h"
#include "TH1.h"

template<class C, class T> TH1D* extract(std::string);
template<class T> class TD;

namespace sf_g{
struct find_maximum{
    double operator()( TH1D const& histogram_p ) const { return histogram_p.GetMaximum(); }
};
struct find_integral{
    double operator()( TH1D const& histogram_p ) const { return histogram_p.Integral(); }
};
struct raw_output{
    void operator()( TH1D&& h1_p, TH1D&& h2_p) const {
        auto * c = new TCanvas{};
        auto * h1_h = new TH1D{ h1_p };
        h1_h->SetDirectory( gROOT );
        h1_h->SetLineColor(kBlue-3);
        h1_h->SetLineWidth(2);
        h1_h->Draw();
        auto * h2_h = new TH1D{ h2_p };
        h2_h->SetDirectory( gROOT );
        h2_h->SetLineColor(kRed+1);
        h2_h->SetLineWidth(2);
        h1_h->Draw("hist");
        h2_h->Draw("same hist");
    }
};
struct difference_output{
    void operator()(TH1D&& h1_p, TH1D&& h2_p) const {
        auto * c = new TCanvas{};
        auto * h_h = new TH1D{ h1_p }; 
        h_h->SetDirectory( gROOT );
        for( auto i{1}; i < h1_p.GetNbinsX(); ++i ){
            h_h->SetBinContent(i, h1_p.GetBinContent(i) - h2_p.GetBinContent(i) );
        }
        h_h->Draw("hist");
    }
};
struct pulled_difference_output{
    void operator()(TH1D&& h1_p, TH1D&& h2_p) const{};
};
namespace tag{
struct both{};
struct first{};
struct second{};
struct maximum{ using operation = find_maximum; };
struct integral{ using operation = find_integral; };
struct raw{ using operation = raw_output;};
struct difference{ using operation = difference_output;};
template<class T> struct pulled_difference{ using operation = pulled_difference_output;};
struct finished{};
} //namespace tag
template<class T> struct scaler{
    template<class U>
    TH1D apply( TH1D&& histogram_p ) const {
        using op = typename T::operation;
        histogram_p.Scale( 1./op{}( histogram_p ) ); 
        return std::move( histogram_p );
    }
};
template<class C, class T> struct extractor{
    template< class U >
    TH1D apply( std::string const& file_p ) const { return *extract<C, T>( file_p ); }
};

template<class T, class Enabler = void> struct is_member_of_formulae : std::false_type{};
template<class T> struct is_member_of_formulae< T, decltype( adl_is_member_of_formulae( std::declval<T>() ) )> : std::true_type{};
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

}; //namespace formuale
template<class Focus, class Formulae> 
struct calibrator{
    static_assert( is_member_of_formulae< Formulae >::value, "this function should be used with a formulae class : i.e. formulae::*****" );
    calibrator( std::string const& file_p ) : formulae_m{ retrieve_parameters(file_p) } {}
    template<class T, typename std::enable_if_t< std::is_same<Focus, tag::both>::value || std::is_same<T, Focus>::value , std::nullptr_t> = nullptr >
    TH1D apply(TH1D&& h_p) const {
        TH1D result{ "h", ";energy (MeV);count", h_p.GetNbinsX(), formulae_m( h_p.GetXaxis()->GetXmin() ), formulae_m( h_p.GetXaxis()->GetXmax() )}; 
        for(auto i{0} ; i < h_p.GetNbinsX(); ++i){ 
            result.SetBinContent( i+1, h_p.GetBinContent(i) );
            result.SetBinError( i+1, h_p.GetBinError(i) );
        }   
        return result;    
    }
    template<class T, typename std::enable_if_t< !std::is_same<Focus, tag::both>::value && !std::is_same<T, Focus>::value , std::nullptr_t> = nullptr >
    TH1D apply(TH1D&& h_p) const{ return std::move(h_p); }
private:
    std::array<double, Formulae::order> retrieve_parameters( std::string const& file_p ){
        std::array<double, Formulae::order> result_c;
        std::ifstream stream{ file_p };
        std::string buffer;
        std::size_t counter{0};
        if( !stream.is_open() ){  std::cerr << "impossible to open file: " << file_p << '\n'; return result_c; }
        while( std::getline(stream, buffer) ){
            if( line_found(buffer, "p") ){result_c[counter++]=find_value( buffer, '=' ); }
        }
        return result_c;
    }  
private: 
    Formulae formulae_m;
};

struct reconstruction_chain{
    struct eraser{
        virtual bool operator()(std::vector<I*> const&) const = 0;
        virtual ~eraser()=default;
    };

    template<class T>
    struct holder : eraser{
        constexpr holder() = default;
        constexpr holder(T&& t_p) : t_m{std::move(t_p)}{}
        bool operator()(std::vector<I*> const& input_pch) const override { return t_m(input_pch); }
        T t_m;
    };

    constexpr cut_vector() { erased_mch.reserve(5); input_mch.reserve(8); redirected_input_mch.reserve(8); }

    private:
    std::vector< std::unique_ptr< eraser> > erased_mch;

};

template<class F>
struct formatter{
    void apply(TH1D&& h1_p, TH1D&& h2_p) const {
        using op = typename F::operation;
//        F::operation{}(std::move(h1_p), std::move(h2_p));
        op{}(std::move(h1_p), std::move(h2_p));
    }
};

namespace details{
template< class ... Ts> struct any_of : std::false_type {};
template< class T> struct any_of<T> : T {};
template< class T, class ... Ts>
struct any_of<T, Ts...> : std::conditional< bool(T::value), T, any_of<Ts...>>::type {};

template<class T, class U> struct is_same_meta : std::false_type {};
template< class T, class U, template<class> class TT>
struct is_same_meta< TT<T>, TT<U> > : std::true_type {};
template< class T, template<class> class TT>
struct is_same_meta< TT<T>, TT<T> > : std::true_type {};
template< class CT, class T, class CU, class U, template<class, class> class TT>
struct is_same_meta< TT<CT, T>, TT<CU, U> > : std::true_type {};
template< class C, class T, class U, template<class, class> class TT>
struct is_same_meta< TT<C, T>, TT<C, U> > : std::true_type {};
template< class C, class T, template<class, class> class TT>
struct is_same_meta< TT<C, T>, TT<C, T> > : std::true_type {};

template<class T, class U> struct is_already_used : std::false_type {};
template<class FT, class FU> struct is_already_used< calibrator< tag::first, FT >, calibrator< tag::first, FU> > : std::true_type{};
template<class FT, class FU> struct is_already_used< calibrator< tag::first, FT >, calibrator< tag::both, FU> > : std::true_type{};
template<class FT, class FU> struct is_already_used< calibrator< tag::second, FT >, calibrator< tag::second, FU> > : std::true_type{};
template<class FT, class FU> struct is_already_used< calibrator< tag::second, FT >, calibrator< tag::both, FU> > : std::true_type{};

} //namespace details

template< class ... Ms > struct comparator{};
template<> 
struct comparator<>{
    template<class C, class T>
    constexpr comparator< extractor<  C, T> > add_extractor() && { return {extractor<C, T>{}}; }
};

comparator<> make_comparator(){ return {}; }

template< class C, class T, class ... Ms>
struct comparator< extractor< C, T>, Ms ...>{
    using tuple_t = std::tuple< extractor< C, T>, Ms ... >;
    static constexpr std::size_t size = std::tuple_size< tuple_t >::value;

    template<class ... Us> friend struct comparator;

    template<class ... Ts>
    constexpr comparator(Ts&& ... ts_p) : module_mc{std::make_tuple(std::move(ts_p)...)} {}

    template< class U,
              typename std::enable_if_t< !details::any_of< details::is_same_meta<scaler<U>, Ms>... >::value, std::nullptr_t> = nullptr>
    constexpr auto add_scaler() && -> comparator< extractor<C, T> , Ms ..., scaler<U>> { 
        return add_impl(scaler<U>{}, std::make_index_sequence<size>{});
    }
    template< class U >
    constexpr auto add_formatter() && -> comparator< tag::finished, extractor<C, T>, Ms..., formatter<U>>{ 
        auto final_comparator = add_impl(formatter<U>{}, std::make_index_sequence<size>{});
        return final_comparator.finish(std::make_index_sequence<size+1>{});
    }
    template< class U, class FU,
              typename std::enable_if_t< !details::any_of< details::is_already_used< calibrator<U, FU>, Ms> ...>::value, std::nullptr_t> = nullptr>
    constexpr auto add_calibrator(std::string const& file_p) && -> comparator< extractor<C, T>, Ms..., calibrator<U, FU>>{
        return add_impl(calibrator<U, FU>{file_p}, std::make_index_sequence<size>{});
    }
private:
    template<class M, std::size_t ... Indices>
    constexpr comparator< extractor<C, T>, Ms ..., M > add_impl( M&& module_p, std::index_sequence<Indices...> ) { return { std::move(std::get< Indices >( module_mc ) )..., std::move( module_p )}; }       
    template<std::size_t ... Indices>
    constexpr comparator< tag::finished, extractor<C, T>, Ms ...> finish( std::index_sequence<Indices...>) { return { std::move( std::get<Indices>(module_mc) ) ... }; }
private:
     tuple_t module_mc;
};

template<class ... Ms>
struct comparator< tag::finished,  Ms ... >{
    constexpr comparator(Ms&& ... ms_p) : module_mc{std::make_tuple(std::move(ms_p)...)} {}
    constexpr void operator()( std::string const& specifier1_p, std::string const& specifier2_p){
        auto histogram1 = apply_chain<tag::first>( specifier1_p, std::make_index_sequence< sizeof...(Ms)-2> {});
        auto histogram2 = apply_chain<tag::second>( specifier2_p, std::make_index_sequence< sizeof...(Ms)-2> {});
        std::get< sizeof...(Ms)-1 >( module_mc).apply( std::move(histogram1), std::move(histogram2) ); 
    };
private:
    template<class T, std::size_t ... Indices>
    TH1D apply_chain( std::string const& specifier_p, std::index_sequence<Indices...> ){
        TH1D result = std::get<0>(module_mc).template apply<T>( specifier_p );
        int expander[] = { 0, (result = std::get<Indices+1>(module_mc).template apply<T>(std::move(result)), 0)...};
        return result;
    }

private:
    std::tuple< Ms ... > module_mc;
};
}//namespace sf_g

#endif

