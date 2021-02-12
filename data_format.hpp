#ifndef DATA_FORMAT_HPP
#define DATA_FORMAT_HPP

#include <string>
#include <fstream>
#include <array>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

namespace sf_g{
    
//------------------------------------------data_stream----------------------------------------
template<class T>
struct data_stream{};

template<>
struct data_stream<std::ifstream> {
    data_stream(std::string input_file_p):
        stream{ check_validity( std::move(input_file_p) )},
        sentinel_m{ compute_length() } {}
    bool end_is_reached() { return stream.eof() || ((sentinel_m - stream.tellg()) == 0); }

private:
    std::ifstream check_validity( std::string input_file_p ) {
        std::ifstream stream{ input_file_p, std::ios_base::binary };
        if( stream.is_open() ){ return stream; }
        throw std::invalid_argument{ "unable to read the given file" };
    }
    std::ifstream::pos_type compute_length() {
        stream.seekg(0, std::ios_base::end );
        auto result = stream.tellg();
        stream.seekg(0, std::ios_base::beg);
        return result;
    }

public:
    std::ifstream stream;
private:
    std::ifstream::pos_type const sentinel_m;
};


template<>
struct data_stream<TTree> {
    data_stream(std::string filename_p):
        file_m{ filename_p.c_str(), "RECREATE" },
        stream_h{ retrieve_tree( )},
        sentinel_m{ static_cast<std::size_t>(stream_h->GetEntries()) } {}
    bool end_is_reached() { return current_entry_m != sentinel_m ; }
    void load_entry() {stream_h->GetEntry(current_entry_m++);}

    ~data_stream() { stream_h->Write(); }

private:
    TTree* retrieve_tree() {
        auto * current_directory_h = TDirectory::CurrentDirectory();
        auto& object_c = *current_directory_h->GetList();
        for( auto * object_h : object_c ){
            auto * tree_h = dynamic_cast<TTree*>( object_h );
            if(tree_h){ return tree_h; }
        }
        return new TTree{ "data", "data from wavecatcher"}; 
    }

private:
    TFile file_m;
    std::size_t current_entry_m{0};
public:
    TTree* stream_h;
private:    
    std::size_t const sentinel_m;
};


//------------------------------multi--------------------------------
template<class ... Ts>
struct multi_input {
    using tuple_t = std::tuple< Ts... >;
    constexpr multi_input(Ts&& ... t_p) : data_mc{std::make_tuple(std::move(t_p)...)} { }
    
private:
    template< std::size_t Index, class F,
              typename std::enable_if_t< (Index < std::tuple_size<tuple_t>::value), std::nullptr_t > = nullptr >
    constexpr void apply_for_each_impl(F&& f_p) 
    {
        f_p( std::get<Index>(data_mc) );
        apply_for_each_impl<Index+1>(std::forward<F>(f_p));
    }
    
    template< std::size_t Index, class F,
              typename std::enable_if_t< (Index >= std::tuple_size<tuple_t>::value), std::nullptr_t > = nullptr >
    constexpr void apply_for_each_impl(F&& /*f_p*/) {}

public: 
    template<class F>
    constexpr void apply_for_each(F&& f_p) 
    {
        apply_for_each_impl<0>( std::forward<F>(f_p) );
    }

private:
    std::tuple< Ts...> data_mc;
};

template<class ... Ts> struct multi_output : Ts... { };

template<class ... Ts> struct composite : Ts... {};

// ------------------------------raw----------------------------------
struct raw_waveform {
    int channel_id;
    int event_id;
    int fcr;
    float baseline;
    float amplitude;
    float charge;
    float leading_edge;
    float trailing_edge;
    float rate_counter;
    std::array<short, 1024> sample_c;
};

struct event_data {
    int event_id;
    double epoch_time;
    struct date_t{ 
        int year;
        int month;
        int day;
    } date;  
    struct time_t{
        int hour;
        int minute;
        int second;
        int millisecond;
    } time;
    int tdc;
    int corrected_tdc;
    int channel_count;
};

struct metadata {
    int channel_count;
    double sampling_period;
};

//-------------------------------------------transformed-------------------------------------------
struct waveform {
    TH1D data; 
};



}//namespace sf_g

#endif
