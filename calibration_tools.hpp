#ifndef CALIBRATION_TOOLS_HPP
#define CALIBRATION_TOOLS_HPP

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
    TH1D * operator()( double maximum_p ) const { return new TH1D{"histogram", ";charge (a.u.);count", static_cast<int>(maximum_p * 1.1e-5/5), 0, maximum_p * 1.1};} 
};        
template<>
struct new_histogram_impl<amplitude>{
    TH1D * operator()( double maximum_p ) const { return new TH1D{"histogram", ";amplitude (a.u.);count", static_cast<int>(maximum_p * 1.1e-5/5) , 0, maximum_p * 1.1};} 
};        
template<class T>
TH1D* new_histogram( double maximum_p ){
    return new_histogram_impl<T>{}( maximum_p );
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
    double maximum{0};
    for(auto i{0}; i < tree_h->GetEntries(); ++i){
        tree_h->GetEntry(i);
        auto current_value = get_value<F, T>(composite_h);
        if( current_value > maximum ){ maximum = current_value; } 
    }

    auto * h_h = new_histogram<T>( maximum ); 
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
    auto fit = h_h->Fit(&f, "RS");
}

void store_fit( std::string output_p ) {
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,24,0)
    auto * current_directory_h = TDirectory::CurrentDirectory().load();
#else
    auto * current_directory_h = TDirectory::CurrentDirectory();
#endif
    auto& object_c = *current_directory_h->GetList();
    TH1D* h_h;
    for( auto * object_h : object_c ){
         if(std::string{ object_h->GetName() }== "histogram" ) { h_h = dynamic_cast<TH1D*>( object_h ); }
    }
    auto * f_h = h_h->GetFunction("f");
    auto fit = h_h->Fit(f_h, "RS");

    std::ofstream output{ output_p };
    auto* old_cout_buffer = std::cout.rdbuf();
    std::cout.rdbuf( output.rdbuf() );
    fit->Print();
    std::cout.rdbuf( old_cout_buffer );
}


template<class F, class T>
void extract_histogram(std::string filename_p) {
    TCanvas * c_h = new TCanvas{};
    auto * h_h = sf_g::extract<F, T>(filename_p) ;
    h_h->Draw("same");
}

void draw_calibration_curve( std::string filename_p) {
    TFile file(filename_p.c_str());
    auto * g_h = static_cast< TGraphErrors*>( file.Get("curve") );
    file.Close();

    auto * c_h = new TCanvas{};
    gStyle->SetOptStat(0);
    auto * energy_array_h = g_h->GetX();
    std::vector<double> energy_c;
    energy_c.assign( energy_array_h, energy_array_h + g_h->GetN());
    auto lower_limit = *std::min_element( energy_c.begin(), energy_c.end() );
    auto upper_limit = *std::max_element( energy_c.begin(), energy_c.end() );
    auto * frame_h = new TH1D{"frame", ";Energy (MeV); charge(a. u.)", 1, 0.8*lower_limit, 1.2*upper_limit};
    auto * mean_array_h = g_h->GetY();
    std::vector<double> mean_c;
    mean_c.assign( mean_array_h, mean_array_h + g_h->GetN());
    auto upper_limit_y = *std::max_element( mean_c.begin(), mean_c.end() );
    frame_h->GetYaxis()->SetRangeUser(0, 1.2* upper_limit_y);
    frame_h->Draw();
    g_h->SetMarkerStyle(21);
    g_h->Draw("ep same");
}

void look_at_waveforms( std::string filename_p ){
    auto get_part_l = []( std::string const& input_p , std::regex regex_p ){
                            std::smatch result;
                            std::regex_search( input_p, result, regex_p);
                            return result[0].str();
                                                              };
    auto file_regex = std::regex{"[^:]+"}; 
    auto branch_regex = std::regex{"[^:]+\\.$"};

    auto* file_h = new TFile( get_part_l(filename_p, std::move(file_regex)).c_str()  );
    TTree* tree_h = static_cast<TTree*>( file_h->Get("data") );
    auto *w_h = new sf_g::waveform{};
    auto branch_name = get_part_l( filename_p, std::move(branch_regex) );
    tree_h->SetBranchAddress( branch_name.c_str(), &w_h);
    tree_h->GetEntry(0);
    gROOT->cd();
    auto * h_h = new TH1D{ w_h->data };
    h_h->SetName( "h" ) ;
    h_h->SetTitle( std::string{ branch_name + "_0" }.c_str() ); 

    auto * c_h = new TCanvas{};
    c_h->Divide( 1, 2);
    auto * hist_pad_h = c_h->cd(1);
    hist_pad_h->SetPad( 0, 0.21, 1, 1 );
    h_h->Draw();
    delete w_h;    
    auto * slider_pad_h = c_h->cd(2);
    slider_pad_h->SetPad(0, 0, 1, 0.2);
    auto * slider_h = new TSlider("slider", "", 0.1, 0.3, 0.9, 0.7); 
    slider_h->SetRange(0., 0.1);
    slider_h->SetMethod( "output_waveform( )");
}

void output_waveform() {
    auto & canvas_c = *gROOT->GetListOfCanvases();
    auto * c_h = static_cast<TCanvas*>( canvas_c.At( canvas_c.GetLast() ) );
    auto * hist_pad_h = c_h->cd(1);
    auto * hist_h = static_cast<TH1D*>(hist_pad_h->FindObject( "h" ));

    auto& file_c = *gROOT->GetListOfFiles();
    auto * file_h = static_cast<TFile*>( file_c.At( file_c.GetLast() ) );
    TTree* tree_h = static_cast<TTree*>( file_h->Get("data") );
    auto *w_h = new sf_g::waveform{};
    auto branch_name = std::string{ hist_h->GetTitle() };
    branch_name.erase( branch_name.find_last_of( "_" ) ); 
    tree_h->SetBranchAddress( branch_name.c_str() , &w_h);

    auto * slider_h = static_cast<TSlider*>(c_h->FindObject("slider")); 
    auto size = slider_h->GetMaximum() - slider_h->GetMinimum();  
    auto entry = static_cast<int>( slider_h->GetMinimum() * (tree_h->GetEntries()-1)/(1-size) ); 
    tree_h->GetEntry( entry );
     
    gROOT->cd();
    hist_h = new TH1D{ w_h->data };
    hist_h->SetName("h");
    hist_h->SetTitle( std::string{ branch_name + "_" + entry }.c_str() );
    hist_pad_h->cd();
    hist_h->Draw();
    delete w_h;
}


#endif

