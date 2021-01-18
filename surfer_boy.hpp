#include <tuple>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

template<class T>
struct mapped_element {
    using type = T;
    T& value;
};

template<class ... Ts>
struct mapper {
    using tuple_type = std::tuple<mapped_element<Ts>...>;

    mapper(mapped_element<Ts>... ts_p) : tuple_m{std::move(ts_p)...} {}

    template< std::size_t Index, class F,
              typename std::enable_if_t< (Index < std::tuple_size<tuple_type>::value), std::nullptr_t > = nullptr >
    constexpr void apply_for_each_impl(F&& f_p)
    {
        f_p( std::get<Index>(tuple_m) );
        apply_for_each_impl<Index+1>(std::forward<F>(f_p));
    }
    
    template< std::size_t Index, class F,
              typename std::enable_if_t< (Index >= std::tuple_size<tuple_type>::value), std::nullptr_t > = nullptr >
    constexpr void apply_for_each_impl(F&& /*f_p*/) {}

    template<class F>
    constexpr void apply_for_each(F&& f_p)
    {
        apply_for_each_impl<0>( std::forward<F>(f_p) );
    }

    private:
    tuple_type tuple_m;
};

template<class ... Ts>
auto make_mapper(Ts&... ts_p){
    return mapper<Ts...>{ mapped_element<Ts>{ts_p}... };
}


namespace policy {
struct deserializer {
   template<class ... Ts>
   constexpr void read_into(std::ifstream& stream_p, mapper<Ts...>&& mapper_p){
       mapper_p.apply_for_each( 
           [&stream_p](auto& element_p){
              using element = typename std::decay<decltype(element_p)>::type;
              stream_p.read( reinterpret_cast<char*>(&element_p.value), sizeof( typename element::type ) );
           } 
       );
   }
};

} //namespace policy

struct measurement : policy::deserializer{
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


struct waveform : policy::deserializer {
    struct metadata {
       std::string value;
        double baseline;
        int channel_count;
    };
    struct data {
        int a;
        float b;
    };

    metadata read_metadata(std::ifstream& stream_p){
        metadata result{};
        std::getline(stream_p, result.value);
        read_into( stream_p, make_mapper(result.baseline, result.channel_count) );
        return result;
    }
    void read_event_header(std::ifstream& stream_p) {
    }
    data read_data(std::ifstream& stream_p) {
        data result;
        read_into( stream_p, make_mapper(result.a, result.b) );
        return result;
    }

};





template<class Policy>
struct reader {
    using data_t = typename Policy::data;
    using metadata_t = typename Policy::metadata;
    
    reader(std::ifstream&& stream_p) : 
        stream_m{ std::move(stream_p) }, 
        metadata{ Policy{}.read_metadata( stream_m  ) } {}

    std::vector<data_t> read_event() {
        std::vector<data_t> data_c;
        data_c.reserve(metadata.channel_count);
        Policy{}.read_event_header( stream_m );
        for(auto i{0}; i < metadata.channel_count ; ++i){
            data_c.push_back( Policy{}.read_data( stream_m ) );
        }
        return data_c;
    }

    private:
    std::ifstream stream_m;
    public:
    metadata_t metadata;
};

template<class Policy>
auto make_reader(std::ifstream stream_p) {
    return reader<Policy>{std::move(stream_p)};
}
