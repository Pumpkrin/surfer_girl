#ifndef POLICY_HPP
#define POLICY_HPP

#include <tuple>
#include <array>
#include <fstream>
#include <iostream>

template<class T>
struct mapped_element {
    T& value;
    constexpr void read(std::ifstream& stream_p){ stream_p.read( reinterpret_cast<char*>(&value), sizeof(T)); }
};

template<class T, std::size_t I>
struct mapped_element< std::array<T, I> > {
    std::array<T,I>& value;
    constexpr void read(std::ifstream& stream_p){ stream_p.read( reinterpret_cast<char*>(value.data()), I*sizeof(T)); }
};

template<class ... Ts>
struct mapper {
    using tuple_type = std::tuple<mapped_element<Ts>...>;

    mapper(mapped_element<Ts>... ts_p) : tuple_m{std::move(ts_p)...} {}
    constexpr void read(std::ifstream& stream_p){ read_impl<0>( stream_p ); }


private:
    template< std::size_t Index,
              typename std::enable_if_t< (Index < std::tuple_size<tuple_type>::value), std::nullptr_t > = nullptr >
    constexpr void read_impl(std::ifstream& stream_p)
    {
        std::get<Index>(tuple_m).read( stream_p );
        read_impl<Index+1>( stream_p );
    }
    
    template< std::size_t Index,
              typename std::enable_if_t< (Index >= std::tuple_size<tuple_type>::value), std::nullptr_t > = nullptr >
    constexpr void read_impl(std::ifstream& /*stream_p*/) {}

   
    private:
    tuple_type tuple_m;
};

template<class ... Ts>
auto make_mapper(Ts&... ts_p){
    return mapper<Ts...>{ mapped_element<Ts>{ts_p}... };
}

template<class ... Ts>
constexpr void read_into(std::ifstream& stream_p, mapper<Ts...>&& mapper_p){
   mapper_p.read( stream_p ); 
}


struct measurement {
    struct metadata {
        int channel_count;
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
        int channel_count;
    };
    struct data {
        int channel_id;
        int event_id;
        float baseline;
        float amplitude;
        float charge;
        float leading_edge;
        float trailing_edge;
        float rate_counter;
    };

    metadata read_metadata(std::ifstream& stream_p){
        metadata result{};
        std::string temp;
        std::getline(stream_p, temp);
        std::getline(stream_p, temp);
        std::getline(stream_p, temp);
        std::getline(stream_p, temp);
        std::getline(stream_p, temp);
        return result;
    }
    void read_event_header(std::ifstream& stream_p) {
        event_data e{};
        read_into( stream_p, make_mapper( e.event_id, e.epoch_time, e.date.year, e.date.month, e.date.day, e.time.hour, e.time.minute, e.time.second, e.time.millisecond, e.tdc, e.channel_count ) );
        std::cout << e.event_id << " -- " << e.epoch_time << " -- ";
        std::cout << e.date.year << "-" << e.date.month << "-" << e.date.day << " -- ";
        std::cout << e.time.hour << "-" << e.time.minute << "-" << e.time.second << "-" << e.time.millisecond << " -- ";
        std::cout << e.tdc << " -- " << e.channel_count << '\n'; 
    }
    data read_data(std::ifstream& stream_p) {
        data result;
        read_into( stream_p, 
                   make_mapper( result.channel_id, 
                                result.event_id, 
                                result.baseline, 
                                result.amplitude, 
                                result.charge,
                                result.leading_edge,
                                result.trailing_edge,
                                result.rate_counter) );
        return result;
    }
};


struct waveform {
    struct metadata {
        int channel_count;
    };
    struct data {
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
        //int end_of_waveform;
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

    metadata read_metadata(std::ifstream& stream_p){
        metadata result{};
        std::string temp;
        std::getline(stream_p, temp);
        std::getline(stream_p, temp);
        std::getline(stream_p, temp);
        std::getline(stream_p, temp);
       //need to retrieve sampling_period and channel_count to fill in metadata 
        return result;
    }
    void read_event_header(std::ifstream& stream_p) {
        event_data e{};
        read_into( stream_p, 
                   make_mapper( e.event_id, 
                                e.epoch_time,
                                e.date.year, 
                                e.date.month,
                                e.date.day, 
                                e.time.hour,
                                e.time.minute, 
                                e.time.second, 
                                e.time.millisecond,
                                e.tdc, 
                                e.corrected_tdc,
                                e.channel_count ) );
        std::cout << e.event_id << " -- " << e.epoch_time << " -- ";
        std::cout << e.date.year << "-" << e.date.month << "-" << e.date.day << " -- ";
        std::cout << e.time.hour << "-" << e.time.minute << "-" << e.time.second << "-" << e.time.millisecond << " -- ";
        std::cout << e.tdc << " -- " << e.corrected_tdc << " -- " << e.channel_count << '\n';
    }
    data read_data(std::ifstream& stream_p) {
        data result;
        read_into( stream_p, 
                   make_mapper( result.channel_id, 
                                result.event_id, 
                                result.fcr,
                                result.baseline, 
                                result.amplitude, 
                                result.charge,
                                result.leading_edge,
                                result.trailing_edge,
                                result.rate_counter,
                                result.sample_c) );
        return result;
    }
    data read_data_with_offset(std::ifstream& stream_p) {
        data result;
        int i;
        read_into( stream_p, 
                   make_mapper( result.channel_id, 
                                result.event_id, 
                                result.fcr,
                                result.baseline, 
                                result.amplitude, 
                                result.charge,
                                result.leading_edge,
                                result.trailing_edge,
                                result.rate_counter,
                                result.sample_c,
                                i) );
//        std::cout << i << '\n';
        return result;
    }

};

#endif
