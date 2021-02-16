//
//File      : flag_set.hpp
//Author    : Alexandre Sécher (alexandre.secher@iphc.cnrs.fr)
//Date      : 02/10/2020
//Framework : PhD thesis, CNRS/IPHC/DRHIM/Hadrontherapy, Strasbourg, France
//

#ifndef flag_set_h
#define flag_set_h

namespace sf_g {
namespace details{
    template<class ... Ts>
    struct any_of : std::false_type {};
    
    template<class T>
    struct any_of<T> : T  {};
    
    template< class T, class ... Ts >
    struct any_of<T, Ts ...> : std::conditional< bool(T::value), T, any_of<Ts...> >::type {};
} // namespace details


struct charge_flag{
    static constexpr uint8_t shift = 4;
};

struct fall_flag{
    static constexpr uint8_t shift = 3;
};

struct rise_flag{
    static constexpr uint8_t shift = 2;
};

struct cfd_flag{
    static constexpr uint8_t shift = 1;
};

struct amplitude_flag{
    static constexpr uint8_t shift = 0;
};

template< class ... Ts >
struct flag_set{
    
public:
    constexpr operator uint8_t() const { return compute_value(); }
    
private:
    flag_set() = default;
    
    constexpr uint8_t compute_value() const {
        uint8_t result = 0;
        int expander[] = { 0, ( result |= 1UL << Ts::shift , void(), 0) ... };
        return result;
        // return 1;
    }
};

template<>
struct flag_set<> {
    constexpr operator uint8_t() const {return 0;}
};

struct extended_flag_set{
    template<class ... Ts, class ... Us>
    friend constexpr extended_flag_set operator+(flag_set<Ts...> t_p, flag_set<Us...> u_p) {
        extended_flag_set result{};
        result.add( t_p );
        result.add( u_p );
        return result; 
    }
    template<class ... Ts>
    friend constexpr extended_flag_set to_final(flag_set<Ts...> t_p); 

public:
    constexpr operator uint32_t() const { return value; }

    template<class ... Ts>
    extended_flag_set& operator+( flag_set<Ts...> t_p ) {
        this->add(t_p);
        return *this;
    }

private:
    extended_flag_set() =default ;
    constexpr void add( uint8_t value_p ) {
        value += value_p << shift;
        shift -= 8;
    }
    uint32_t value{0};
    unsigned int shift{24};
};

uint8_t make_opcode( std::string& module_list_p ) {
    uint8_t opcode{0};
    for( auto& token : module_list_p ){ 
        switch(token){
        case 'a': { opcode |= sf_g::flag_set< sf_g::amplitude_flag >{} ; break ; }
        case ':': { break; }           
        case '{': { break; }           
        case '}': { break; }           
        default : { std::cerr << "modifier " << token << " not implemented yet\n"; break ; }
        }
    }  
    return opcode;
}

uint32_t make_opcode( std::vector<uint8_t> opcode_pc ) {
    uint32_t result{0};
    if( opcode_pc.size() > 4 ){ throw std::runtime_error{"this amount of channel is not supported yet"}; } 
    unsigned int shift = 24;
    for( auto opcode : opcode_pc ){ 
        result += opcode << shift;
        shift -= 8;   
    }
    return result;
}  

template<class ... Ts>
constexpr extended_flag_set to_final(flag_set<Ts...> t_p) {
   extended_flag_set result{};
   result.add( t_p );
   return result; 
}

} //namespace sf_g

#endif /* flag_set_h */
