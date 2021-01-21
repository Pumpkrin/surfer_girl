#ifndef WRITER_HPP
#define WRITER_HPP

#include "policy.hpp"

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

template< class Policy >
struct root_converter{};

template<>
struct root_converter<waveform> {
    using input_data_t = typename waveform::data;
    struct data {
        int fcr;
        TH1D waveform;
    };
    data operator()(input_data_t& data_p) const {
        data result{};
        result.fcr = data_p.fcr;
        result.waveform = TH1D{ "", ";time;ADC", 1024, 0, 1 };
        for( auto i{0}; i < 1024 ; ++i ){
            result.waveform.SetBinContent( i+1, data_p.sample_c[i] );
        }
        return result;
    }
};

template<>
struct root_converter<measurement> {
    struct data {

    };

};


template<class Policy>
struct writer {
    using data_t = typename Policy::data;
    using root_data_t = typename root_converter<Policy>::data;

    writer( int channel_count_p, std::string output_file_p ) : 
        channel_count_m{ static_cast<size_t>( channel_count_p ) }, 
        ouput_m{ output_file_p.c_str(), "RECREATE" },
        data_mc{ std::vector<root_data_t>{ channel_count_m } },
        tree_m{ "data", "converted data from wavecatcher system" } {
            for( auto i{0}; i < channel_count_m ; ++i) {
                std::string name = "channel_" + std::to_string(i);
                tree_m.Branch( name.c_str(), &data_mc[i] ); 
            } 
        }

    constexpr void write_data( std::vector<data_t> data_pc ) {
        for( auto i{0} ; i < channel_count_m ; ++i ) {
           data_mc[i] = root_converter<Policy>{}( data_pc[i] ) ;
//           std::cout << data_mc[i].fcr << '\n';
        }
        tree_m.Fill();
    } 

    void save_data() {tree_m.Write();}

    private:
    size_t channel_count_m;
    TFile ouput_m;
    std::vector<root_data_t> data_mc;
    TTree tree_m;
};

#endif
