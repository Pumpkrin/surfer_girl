#ifndef TOOLS_HPP
#define TOOLS_HPP

#include "formulae.hpp"

#include <regex>
#include <string> 

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TSlider.h"


namespace sf_g {

template<class Composite, class T> struct get_value_impl{};
template<class Composite>
struct get_value_impl<Composite, charge > {
    constexpr double operator()( Composite* composite_ph ) const {return composite_ph->charge;}
};
template<class Composite>
struct get_value_impl<Composite, amplitude > {
    constexpr double operator()( Composite* composite_ph ) const {return composite_ph->amplitude;}
};
template<class Composite>
struct get_value_impl<Composite, cfd_time > {
    constexpr double operator()( Composite* composite_ph ) const {return composite_ph->time;}
};
template<class Composite>
struct get_value_impl<Composite, pile_up > {
    constexpr double operator()( Composite* composite_ph ) const {return composite_ph->pile_up;}
};


template<class Composite, class T>
constexpr double get_value( Composite* composite_ph){
    return get_value_impl<Composite, T>{}( composite_ph );
};

template<class T> struct new_histogram_impl{};
template<>
struct new_histogram_impl<charge>{
    TH1D * operator()( double maximum_p ) const { return new TH1D{"histogram", ";charge (a.u.);count", static_cast<int>(maximum_p * 1.1e-5/5), 0, maximum_p * 1.1};} 
};        
template<>
struct new_histogram_impl<amplitude>{
//    TH1D * operator()( double maximum_p ) const { return new TH1D{"histogram", ";amplitude (a.u.);count", static_cast<int>(maximum_p * 1.1/75) , 0, maximum_p * 1.1};} 
//    TH1D * operator()( double maximum_p ) const { return new TH1D{"histogram", ";amplitude (a.u.);count", maximum_p * 1.1/30 > 100 ? static_cast<int>(maximum_p * 1.1/30) : 100  , 0, maximum_p * 1.1};} 
    TH1D * operator()( double maximum_p ) const { return new TH1D{"histogram", ";amplitude (a.u.);count", 1000 ,100 , 7000};} 
};        
template<>
struct new_histogram_impl<pile_up>{
    TH1D * operator()( double maximum_p ) const { return new TH1D{"histogram", ";pile_up (a.u.);count", static_cast<int>(maximum_p * 1.1/30) , 0, maximum_p * 1.1};} 
};        
template<class T>
TH1D* new_histogram( double maximum_p ){
    return new_histogram_impl<T>{}( maximum_p );
}

template<class Composite, class T>    
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
    auto *composite_h = new Composite{};
    tree_h->SetBranchAddress( get_part_l( filename_p, std::move(branch_regex) ).c_str(), &composite_h);
    double maximum{0};
    for(auto i{0}; i < tree_h->GetEntries(); ++i){
        tree_h->GetEntry(i);
        auto current_value = get_value<Composite, T>(composite_h);
        if( current_value > maximum ){ maximum = current_value; } 
    }

    auto * h_h = new_histogram<T>( maximum ); 
//    auto * h_h = new_histogram<T>( 50e6 ); 
    h_h->SetDirectory( gROOT );

    for(auto i{0}; i < tree_h->GetEntries(); ++i){
        tree_h->GetEntry(i);
        h_h->Fill( get_value<Composite, T>(composite_h) );
    }

    delete composite_h;
    return h_h;
}

}//namespace sf_g

template<class Composite, class T>
void fit(std::string filename_p) {
    TCanvas * c_h = new TCanvas{};
    auto * h_h = sf_g::extract<Composite, T>(filename_p) ;
    h_h->Draw("same");

    double peak = h_h->GetBinCenter( h_h->GetMaximumBin() ); 
    TF1 f{ "f", "gaus", 0.9 * peak, 1.1 * peak};
    auto fit = h_h->Fit(&f, "RS");
}

void compute_integral_from_last_fit(std::string name_p) {
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,24,0)
    auto * current_directory_h = TDirectory::CurrentDirectory().load();
#else
    auto * current_directory_h = TDirectory::CurrentDirectory();
#endif
    auto& object_c = *current_directory_h->GetList();
    TH1D* h_h;
    for( auto * object_h : object_c ){
         if(std::string{ object_h->GetName() }==name_p  ) { h_h = dynamic_cast<TH1D*>( object_h ); }
    }
    if(!h_h){ std::cerr << "No histogram found\n"; return;}

    auto *list_h = h_h->GetListOfFunctions();
    auto * f_h = static_cast<TF1*>(list_h->Last());
//    auto * f_h = h_h->GetFunction("PrevFitTMP");
    if(!f_h){ std::cerr << "no function found\n";}
    auto fit_h = h_h->Fit(f_h, "S");
    auto lower_limit = fit_h->Parameter(1) - 3*fit_h->Parameter(2);
    auto upper_limit = fit_h->Parameter(1) + 3*fit_h->Parameter(2);
    auto width = h_h->GetBinWidth( h_h->FindBin( lower_limit ) );
    std::cout << "Mean: " << fit_h->Parameter(1) << " +/- " << fit_h->Error(1) << '\n';
    std::cout << "Integral: " << f_h->Integral(lower_limit, upper_limit)/width << " +/- " << f_h->IntegralError(lower_limit, upper_limit, fit_h->GetParams(), fit_h->GetCovarianceMatrix().GetMatrixArray())/width << "\n";
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


template<class Composite, class T>
void extract_histogram(std::string filename_p) {
    TCanvas * c_h = new TCanvas{};
    auto * h_h = sf_g::extract<Composite, T>(filename_p) ;
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


template< class Composite >
void generate_time_distribution( std::string filename ) {
    auto regex_split_l = []( std::string const & text_p, std::regex regex_p ) {
        std::vector<std::string> result_c;
        result_c.reserve(10);
        
        std::sregex_iterator match_i{ text_p.begin(), text_p.end(), regex_p };
        auto end_i = std::sregex_iterator{};
        for(auto iterator = match_i; iterator != end_i ; ++iterator){result_c.push_back( iterator->str() );}
        return result_c;
    };
    auto name_c = regex_split_l( filename, std::regex{"[^:]+"}); 

    auto* file_h = new TFile( name_c.front().c_str()  );
    TTree* tree_h = static_cast<TTree*>( file_h->Get("data") );
    auto *c1_h = new Composite{};
    auto *c2_h = new Composite{};
    tree_h->SetBranchAddress( name_c[1].c_str(), &c1_h);
    tree_h->SetBranchAddress( name_c[2].c_str(), &c2_h);
    
    gROOT->cd();
    auto * h_h =  new TH1D{"histogram", ";time difference (ps);count", 2000 , -2.5e6, 2.5e6};

    for(auto i{0}; i < tree_h->GetEntries(); ++i){
        tree_h->GetEntry(i);
        h_h->Fill( c2_h->time - c1_h->time);
    }

    h_h->Draw();
}

void draw_qqplot( TH1D* distribution1_p, TH1D* distribution2_p, double step_p ) {
    distribution1_p->Scale(1./distribution1_p->Integral());
    distribution2_p->Scale(1./distribution2_p->Integral());
    
    std::vector<double> abscissa_c;
    std::vector<double> ordinate_c;

    struct qq_data{
        double abscissa;
        double ordinate;
        double quantile;
        double quantile_abscissa;
        double quantile_ordinate;
    };

    std::vector<qq_data> qq_data_c;
    qq_data_c.reserve( 1./step_p);
    for(double quantile{step_p}; quantile < 1 ; quantile += step_p ){
        std::cout << "aiming for: " << quantile << '\n';

        int index1{1}, index2{1};
        if(qq_data_c.empty()){qq_data_c.push_back( qq_data{} );}
        else{ 
            qq_data_c.push_back( qq_data_c.back() ) ; 
            qq_data_c.back().quantile = 0;
            index1 = distribution1_p->FindBin( qq_data_c.back().abscissa) + 1; 
            index2 = distribution2_p->FindBin( qq_data_c.back().ordinate) + 1;
        }
        auto & qq_data = qq_data_c.back();

        while( qq_data.quantile_abscissa < quantile && index1 < distribution1_p->GetNbinsX() +1 ){qq_data.quantile_abscissa += distribution1_p->GetBinContent(index1++);} 
        while( qq_data.quantile_ordinate < quantile && index2 < distribution2_p->GetNbinsX() +1 ){qq_data.quantile_ordinate += distribution2_p->GetBinContent(index2++);} 
        std::cout << qq_data.quantile_abscissa << " -> " << qq_data.quantile_ordinate << '\n';

        if( qq_data.quantile_abscissa > quantile + step_p || qq_data.quantile_ordinate > quantile + step_p ){ std::cout << qq_data_c.size() << "erasing\n"; qq_data_c.erase( qq_data_c.end() -1) ;  continue; }

        qq_data.abscissa = distribution1_p->GetBinCenter(index1);
        qq_data.ordinate = distribution2_p->GetBinCenter(index2);
        qq_data.quantile = quantile;
        std::cout << quantile << ": " << qq_data.abscissa << " -> " << qq_data.ordinate <<  '\n';
    }

    qq_data_c.erase( std::remove_if( qq_data_c.begin(), qq_data_c.end(), [](auto const& qq_data_p){ return qq_data_p.quantile == 0; } ),
                     qq_data_c.end() );  
    std::for_each(qq_data_c.begin(), qq_data_c.end(), [&abscissa_c, &ordinate_c](auto const& v_p){ abscissa_c.push_back( v_p.abscissa ); ordinate_c.push_back( v_p.ordinate );  });

    for(auto i{0}; i < abscissa_c.size() ; ++i){
        std::cout << abscissa_c[i] << " -> " << ordinate_c[i] << '\n';
    }
    TGraph * g_h = new TGraph{ static_cast<int>(abscissa_c.size()), abscissa_c.data(), ordinate_c.data()};
    TCanvas * c_h = new TCanvas{};
    g_h->Draw("a*");
    TF1 * f_h = new TF1{ "f", "x",distribution1_p->GetBinCenter(1), distribution1_p->GetBinCenter( distribution1_p->GetNbinsX() ) };
    f_h->Draw("same");


}
void draw_reverse_qqplot( TH1D* distribution1_p, TH1D* distribution2_p, double step_p ) {
    distribution1_p->Scale(1./distribution1_p->Integral());
    distribution2_p->Scale(1./distribution2_p->Integral());
    
    std::vector<double> abscissa_c;
    std::vector<double> ordinate_c;

    struct qq_data{
        double abscissa;
        double ordinate;
        double quantile;
    };

//    std::vector<double> quantile_c{0.91, 1.22, 1.5, 1.6, 1.84, 2, 2.33, 2.6, 2.8, 3.1, 4.25, 4.65, 5.35, 5.9, 6.3, 6.8, 7.3};
    std::vector<double> quantile_c;
    int const size = distribution1_p->GetNbinsX();
    double range = distribution1_p->GetBinLowEdge( size ) + distribution1_p->GetBinWidth(size) - distribution1_p->GetBinLowEdge(1); 
    for(auto i{1}; i< range/step_p; ++i){ quantile_c.push_back( i * step_p); } 
    
    std::vector<qq_data> qq_data_c;
    qq_data_c.reserve(quantile_c.size());
    for(auto const& quantile : quantile_c){
        qq_data_c.push_back( qq_data{} );
        auto & qq_data = qq_data_c.back();

        auto index{0};
        while( distribution1_p->GetBinCenter(index) < quantile && index < distribution1_p->GetNbinsX() +1 ){
            qq_data.abscissa += distribution1_p->GetBinContent(index++);
        } 
        index = 0;
        while( distribution2_p->GetBinCenter(index) < quantile && index < distribution2_p->GetNbinsX() +1 ){
            qq_data.ordinate += distribution2_p->GetBinContent(index++);
        } 
        qq_data.quantile = quantile;
        std::cout << quantile << ": " << qq_data.abscissa << " -> " << qq_data.ordinate <<  '\n';
    }

    std::for_each(qq_data_c.begin(), qq_data_c.end(), [&abscissa_c, &ordinate_c](auto const& v_p){ abscissa_c.push_back( v_p.abscissa ); ordinate_c.push_back( v_p.ordinate );  });

    for(auto i{0}; i < abscissa_c.size() ; ++i){
        std::cout << abscissa_c[i] << " -> " << ordinate_c[i] << '\n';
    }
    TGraph * g_h = new TGraph{ static_cast<int>(abscissa_c.size()), abscissa_c.data(), ordinate_c.data()};
    TCanvas * c_h = new TCanvas{};
    g_h->Draw("a*");
    TF1 * f_h = new TF1{ "f", "x",distribution1_p->GetBinCenter(1), distribution1_p->GetBinCenter( distribution1_p->GetNbinsX() ) };
    f_h->Draw("same");

    TH1D * h_h = new TH1D{"h", "E (MeV)", static_cast<int>(abscissa_c.size()), 0, 10};
    for(auto i{0}; i < abscissa_c.size(); ++ i){
        h_h->Fill( quantile_c[i], sqrt(2*( pow(abscissa_c[i],2) + pow(ordinate_c[i],2))/2));
    }
    TCanvas * c2_h= new TCanvas{};
    h_h->Draw();
}

TH1D* compute_cumulative_distribution(TH1D * distribution_ph){
    int const size = distribution_ph->GetNbinsX();
    TH1D * value_h = new TH1D{ "h", "; E(MeV);", size, 0, 10};
    auto value{0.};
    for(auto i{1}; i < size +1 ; ++i){
        value += distribution_ph->GetBinContent(i);
        value_h->SetBinContent(i, value);
    }
    return value_h;
}

void draw_difference(TH1D * distribution1_ph,TH1D * distribution2_ph){
    int const size = distribution1_ph->GetNbinsX();
    TH1D * value_h = new TH1D{ "h", "; E(MeV);", size, 0, 10};
    for(auto i{1}; i < size +1 ; ++i){
        value_h->SetBinContent(i, distribution1_ph->GetBinContent(i) - distribution2_ph->GetBinContent(i));    
    }
    TCanvas * c_h = new TCanvas{};
    value_h->Draw();

} 

#endif

