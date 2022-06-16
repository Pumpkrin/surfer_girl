#ifndef DATA_STRUCTURE_HPP
#define DATA_STRUCTURE_HPP

#include <array>
#include <iostream>
#include "TH1.h"

namespace sf_g {

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

struct linked_waveform {
    TH1D data;
    std::size_t channel_number;
};

struct amplitude { double amplitude; void value() const {  std::cout << amplitude << '\n';} };
struct baseline { double baseline; void value() const {  std::cout << baseline << '\n';} };
struct cfd_time { double time; void value() const { std::cout << time << '\n';} };
struct charge { double charge; void value() const { std::cout << charge << '\n';} };    
struct rise_time { double rise_time; void value() const { std::cout << rise_time << '\n';} };    
struct fall_time { double fall_time; void value() const { std::cout << fall_time << '\n';} };    
struct mean { double mean; void value() const { std::cout << mean << '\n';} };    
struct sigma { double sigma; void value() const { std::cout << sigma << '\n';} };    
struct pile_up { double pile_up; void value() const { std::cout << pile_up << '\n';} };    


struct gamma_response{
    double gamma_energy;
    double deposited_energy;
};
struct response{
    double gamma_energy;
    double deposited_energy;
    double total_deposited_energy;
};

}
#endif
