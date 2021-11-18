#include "data_structure.hpp"
#include "process_spectra.hpp"
#include "formulae.hpp"

#include <string> 

#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TH2.h"
#include "TROOT.h"


int main(int argc, char* argv[]){
    std::string input_file;
    std::string output_file;
    std::size_t bin_count = 512;
    double lower_range = 0;
    double upper_range = 10;
    std::string resolution_file;
    for( auto i{1}; i < argc ; ++i ) { 
        std::string argument = argv[i];
        if( argument == std::string{"-in"} ){ input_file = argv[++i]; }
        if( argument == std::string{"-out"} ){ output_file = argv[++i]; }
        if( argument == std::string{"-resolution"} ){ resolution_file = argv[++i]; }
        if( argument == std::string{"-bin"} ){ bin_count = std::stoi(argv[++i]); }
        if( argument == std::string{"-range"} ){ lower_range = std::stod(argv[++i]); upper_range = std::stod(argv[++i]); }
    }  
    sf_g::smearer<sf_g::formulae::gamma_resolution> s{resolution_file};

    TFile input{ input_file.c_str() };
    auto* tree_h = input.Get<TTree>("data");
    std::unique_ptr<sf_g::gamma_response> holder_h{new sf_g::gamma_response{}};
    auto * response_h = holder_h.get();
    tree_h->SetBranchAddress("gamma_response", &response_h); 

    gROOT->cd();
    TH2D response_matrix{"response_matrix",";deposited_energy (MeV);gamma_energy (MeV)", static_cast<int>(bin_count), lower_range, upper_range, static_cast<int>(bin_count), lower_range, upper_range };
    TH1D spectrum{"spectrum",";deposited_energy (MeV);count", static_cast<int>(bin_count), lower_range, upper_range };

    for( std::size_t i{0}; i < tree_h->GetEntries() ; ++i ){
        tree_h->GetEntry(i);
//        if(response.deposited_energy > 0.7 && response.deposited_energy< 15){
            response_h->deposited_energy = s( response_h->deposited_energy );
            response_matrix.Fill( response_h->deposited_energy, response_h->gamma_energy );
            spectrum.Fill(response_h->deposited_energy);
            //fill result
//        }
    }
    input.Close();

    //rescale so to integral for 1 gamma energy
    std::vector<std::size_t> content_c( bin_count );
    for( std::size_t i{0}; i < bin_count * bin_count ; ++i ){
        content_c[i%bin_count] += response_matrix.GetBinContent( i%bin_count + 1, i/bin_count +1 );      
    }
    for( std::size_t i{0}; i < bin_count * bin_count ; ++i ){
        if(content_c[i%bin_count] > 0){ response_matrix.SetBinContent( i%bin_count+1, i/bin_count+1, response_matrix.GetBinContent( i%bin_count + 1, i/bin_count +1 )/content_c[i%bin_count] );      }
    }
    

    //fill output
    TFile output{ output_file.c_str(), "UPDATE" };
    spectrum.Write();
    response_matrix.Write();
    output.Close();
}
