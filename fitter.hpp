#ifndef FITTER_HPP
#define FITTER_HPP

#include <regex>
#include <string> 

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TF1.h"

namespace sf_g {

template<class F, class T> struct get_value_impl{};
template<class F>
struct get_value_impl<F, charge > {
    constexpr double operator()( F* composite_ph ) const {return composite_ph->charge;}
};
template<class F>
struct get_value_impl<F, amplitude > {
    constexpr double operator()( F* composite_ph ) const {return composite_ph->amplitude;}
};


template<class F, class T>
constexpr double get_value( F* composite_ph){
    return get_value_impl<F, T>{}( composite_ph );
};

template<class T> struct new_histogram_impl{};
template<>
struct new_histogram_impl<charge>{
    TH1D * operator()() const { return new TH1D{"charge", ";charge (a.u.);count", 1000, 0, 200e6};} 
};        
template<>
struct new_histogram_impl<amplitude>{
    TH1D * operator()() const { return new TH1D{"amplitude", ";amplitude (a.u.);count", 1000, 0, 10000};} 
};        
template<class T>
TH1D* new_histogram(){
    return new_histogram_impl<T>{}();
}

template<class F, class T>    
TH1D* extract( std::string filename_p ) {
    auto get_part_l = []( std::string const& input_p , std::regex regex_p ){
                            std::smatch result;
                            std::regex_search( input_p, result, regex_p);
                            return result[0].str();
                                                              };
    auto file_regex = std::regex{"[^:]+"}; 
    auto branch_regex = std::regex{"[^:]+\\.$"};

    TFile file( get_part_l(filename_p, std::move(file_regex)).c_str()  );
    TTree* tree_h = static_cast<TTree*>( file.Get("data") );
    auto *composite_h = new F{};
    tree_h->SetBranchAddress( get_part_l( filename_p, std::move(branch_regex) ).c_str(), &composite_h);

    auto * h_h = new_histogram<T>();
    h_h->SetDirectory( gROOT );

    for(auto i{0}; i < tree_h->GetEntries(); ++i){
        tree_h->GetEntry(i);
        h_h->Fill( get_value<F, T>(composite_h) );
    }

    delete composite_h;
    return h_h;
}

}//namespace sf_g

template<class F, class T>
void fit(std::string filename_p) {
    TCanvas * c_h = new TCanvas{};
    auto * h_h = sf_g::extract<F, T>(filename_p) ;
    h_h->Draw("same");

    double peak = h_h->GetBinCenter( h_h->GetMaximumBin() ); 
    TF1 f{ "f", "gaus", 0.9 * peak, 1.1 * peak};
    h_h->Fit(&f, "R");
}



#endif

