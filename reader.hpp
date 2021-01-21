#ifndef READER_HPP
#define READER_HPP

#include <fstream>
#include <string>
#include <vector>


template<class Policy>
struct reader {
    using data_t = typename Policy::data;
    using metadata_t = typename Policy::metadata;
    
    reader(std::string input_file_p) : 
        stream_m{ check_validity(std::move(input_file_p)) }, 
        end_m{ compute_length() }, 
        metadata{ Policy{}.read_metadata( stream_m  ) } {}

    std::vector<data_t> read_event() {
        std::vector<data_t> data_c;
        data_c.reserve(metadata.channel_count);
        Policy{}.read_event_header( stream_m );
        for(auto i{0}; i < metadata.channel_count ; ++i){
            switch( metada.channel_count) {
                case 1: {
                data_c.push_back( Policy{}.read_data_with_offset( stream_m ) );
                break;
                        }
                default:
                data_c.push_back( Policy{}.read_data( stream_m ) );
                break;
            }
        }
        return data_c;
    } 
    
    bool end_is_reached() { 
        auto current_position = stream_m.tellg();
        auto distance = end_m - current_position;
        return distance == 0 ; 
    }

    private:
    std::ifstream::pos_type compute_length() {
        stream_m.seekg(0, std::ios_base::end);
        auto result = stream_m.tellg();
        stream_m.seekg(0, std::ios_base::beg);
        return result;
    }
    std::ifstream check_validity( std::string input_file_p ) {
        std::ifstream stream{ input_file_p, std::ios_base::binary };
        if( stream.is_open() ){ return stream; }
        throw std::invalid_argument{ "unable to read the given input file"};
    }

    private:
    std::ifstream stream_m;
    const std::ifstream::pos_type end_m;
    public:
    metadata_t metadata;
};

template<class Policy>
auto make_reader(std::ifstream stream_p) {
    return reader<Policy>{std::move(stream_p)};
}

#endif
