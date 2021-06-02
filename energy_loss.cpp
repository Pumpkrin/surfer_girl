#include <fstream>
#include <string>
#include <iostream>

#include "TGraph.h"

struct value_and_error{
    double value;
    double error;
};
    
struct stopping_power{
    stopping_power( std::string const& input_p ) : data_m{ generate_bethe_bloch(input_p) } {}
    value_and_error evaluate( double energy_p ){
        value_and_error result;
        result.value = data_m.Eval(energy_p);
        result.error = compute_error( energy_p );
        return result;
    }

private:
    double compute_error( double energy_p ) const;
    TGraph generate_bethe_bloch( std::string const& input_p ) const;
private:
    TGraph data_m;
};

TGraph stopping_power::generate_bethe_bloch( std::string const& input_p ) const {
    std::ifstream input_flux( input_p );
    std::string buffer;
    
    Double_t kinetic_energy[132]={0};
    Double_t stopping_power[132]={0};
    
    if(input_flux){
        for(Int_t i=0;i<140;i++){
            if(i<8){std::getline(input_flux, buffer);}
            else{input_flux>>kinetic_energy[i-8]>>stopping_power[i-8];}
        }
    }
    input_flux.close();
    
    return TGraph(132,kinetic_energy,stopping_power);
}

double stopping_power::compute_error( double energy_p ) const{
    if(energy_p >= 10){ return 0.02 * data_m.Eval(energy_p);}
    else if(energy_p<10 && energy_p>=1){return 0.05*data_m.Eval(energy_p);}
    else if(energy_p<1 && energy_p>=0.1){return 0.1*data_m.Eval(energy_p);}
    else if(energy_p<0.1 && energy_p>=0.01){return 0.15*data_m.Eval(energy_p);}
    else{return 0.3*data_m.Eval(energy_p);}
}


int main( int argc, char* argv[] ) {
    
    double step = 10, position = 0, thickness = 0, density = 2.7 ;
    value_and_error energy;
    std::string input_file{ "data/p_aluminum.txt"};
    if(argc < 7){ std::cerr << "energy_loss should be called the following way: ./energy_loss -in input_file.txt -thick value -ene value \n"; return 1; }
    else{     
        for(auto i{0}; i < argc; ++i) {
            if( std::string( argv[i] ) == "-in" ){ input_file = std::string{ argv[++i] }; }
            if( std::string( argv[i] ) == "-thick") { thickness = std::stod( argv[++i]); }
            if( std::string( argv[i] ) == "-ene") { energy.value = std::stod( argv[++i]); }
            if( std::string( argv[i] ) == "-density") { density = std::stod( argv[++i]); }
        }
    }

    stopping_power sp{ input_file }; 
    
    while( position + step < thickness && energy.value > 0){
        auto current_sp = sp.evaluate( energy.value);
        energy.value -= current_sp.value * density * step * 1e-4;        
        energy.error += pow( current_sp.error * density * step * 1e-4, 2 ); 
        position+=step;
    }

    if(energy.value>0){
        auto current_sp = sp.evaluate( energy.value);
        energy.value -= current_sp.value * density * (thickness - position)  * 1e-4;        
        energy.error += pow( current_sp.error * density * (thickness - position) * 1e-4, 2 ); 
        energy.error = sqrt( energy.error);
        std::cout << "energy_remaining: " << energy.value << " +/- " << energy.error << " MeV\n";
    } else { std::cout << "no remaining energy\n"; }
        
}
