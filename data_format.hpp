#ifndef DATA_FORMAT_HPP
#define DATA_FORMAT_HPP

#include <string>
#include <fstream>
#include <array>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

#include <iostream>
namespace sf_g{
    
//------------------------------------------data_stream----------------------------------------
template<class T>
struct data_input{};

template<>
struct data_input<std::ifstream> {
    data_input(std::string input_file_p):
        input{ check_validity( std::move(input_file_p) )},
        sentinel_m{ compute_length() } {}
    bool end_is_reached() { return input.eof() || ((sentinel_m - input.tellg()) == 0); }

private:
    std::ifstream check_validity( std::string input_file_p ) {
        std::ifstream input{ input_file_p, std::ios_base::binary };
        if( input.is_open() ){ return input; }
        throw std::invalid_argument{ "unable to read the given file" };
    }
    std::ifstream::pos_type compute_length() {
        input.seekg(0, std::ios_base::end );
        auto result = input.tellg();
        input.seekg(0, std::ios_base::beg);
        return result;
    }

public:
    std::ifstream input;
private:
    std::ifstream::pos_type const sentinel_m;
};

template<>
struct data_input<TTree> {
    data_input(std::string filename_p):
        file_m{ filename_p.c_str(), "READ" },
        input_h{ retrieve_tree( )},
        sentinel_m{ static_cast<std::size_t>(input_h->GetEntries()) } {}
    bool end_is_reached() { return current_entry_m == sentinel_m ; }
    void load_entry() {input_h->GetEntry(current_entry_m++);}

private:
    TTree* retrieve_tree() {
        auto * current_directory_h = TDirectory::CurrentDirectory();
        auto& key_c = *current_directory_h->GetListOfKeys();
        for( auto * key_h : key_c ){
            auto * tree_h = dynamic_cast<TTree*>( file_m.Get( key_h->GetName() ));
            if(tree_h){ return tree_h; }
        }
        throw std::runtime_error{"no viable tree found in input file"};
    }

private:
    TFile file_m;
    std::size_t current_entry_m{0};
public:
    TTree* input_h;
private:    
    std::size_t const sentinel_m;
};

template<class T> struct data_output{};

template<>
struct data_output<TTree> {
    data_output(std::string filename_p):
        file_m{ filename_p.c_str(), "RECREATE" },
        output{ "data", "data from wavecatcher"} {}
    ~data_output() { file_m.cd(); output.Write(); }

private:
    TFile file_m;
public:
    TTree output;
};


//------------------------------multi--------------------------------
template<class ... Ts>
struct multi_input {
    using tuple_t = std::tuple< Ts... >;
    
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

template<class ... Ts> struct composite : Ts... {
    void value() { 
        int expander[] = {0, ( static_cast<Ts&>(*this).value(), void(), 0)...};
    }
};

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

struct amplitude { double amplitude; void value() const {  std::cout << amplitude << '\n';} };
struct baseline { double baseline; void value() const {  std::cout << baseline << '\n';} };
struct cfd_time { double time; void value() const { std::cout << time << '\n';} };
struct charge { double charge; void value() const { std::cout << charge << '\n';} };    

}//namespace sf_g

#endif
