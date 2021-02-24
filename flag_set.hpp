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

struct baseline_flag{
    static constexpr uint8_t shift = 5;
};

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

struct combined_flag_set{
    template<class ... Fs>
    friend constexpr combined_flag_set to_opcode(flag_set<Fs...> f_p); 
    template<class ... Ts>
    friend constexpr combined_flag_set to_opcode(Ts&&... ts_p); 

public:
    constexpr operator uint32_t() const { return value; }
    constexpr void add( uint8_t value_p ) {
        value += value_p << shift;
        shift -= 8;
    }

private:
    combined_flag_set() = default ;
    uint32_t value{0};
    unsigned int shift{24};
};

uint8_t make_single_opcode( std::string& module_list_p ) {
    uint8_t opcode{0};
    for( auto& token : module_list_p ){ 
        switch(token){
        case 'a': { opcode |= sf_g::flag_set< sf_g::amplitude_flag >{} ; break ; }
        case 't': { opcode |= sf_g::flag_set< sf_g::cfd_flag >{} ; break ; }
        case 'b': { opcode |= sf_g::flag_set< sf_g::baseline_flag >{} ; break ; }
        case 'c': { opcode |= sf_g::flag_set< sf_g::charge_flag >{} ; break ; }
        case ':': { break; }           
        case '{': { break; }           
        case '}': { break; }           
        default : { std::cerr << "modifier " << token << " not implemented yet\n"; break ; }
        }
    }  
    return opcode;
}

uint32_t make_combined_opcode( std::vector<uint8_t> opcode_pc ) {
    uint32_t result{0};
    if( opcode_pc.size() > 4 ){ throw std::runtime_error{"this amount of channel is not supported yet"}; } 
    unsigned int shift = 24;
    for( auto opcode : opcode_pc ){ 
        result += opcode << shift;
        shift -= 8;   
    }
    return result;
}  

template<class ... Fs>
constexpr combined_flag_set to_opcode(flag_set<Fs...> f_p) {
   combined_flag_set result{};
   result.add( f_p );
   return result; 
}

namespace details{
template<class ... Ts> struct to_opcode_impl;
template<class ... Fs, class ... Ts> 
struct to_opcode_impl< flag_set<Fs...>, Ts...>{
    constexpr void operator()(combined_flag_set& set_p, flag_set<Fs...> flag_p, Ts&&... ts_p) const {
        set_p.add( flag_p );
        to_opcode_impl<Ts...>{}( set_p, std::move(ts_p)... );
    }
};
template<class ...Fs>
struct to_opcode_impl< flag_set<Fs...> >{
    constexpr void operator()(combined_flag_set& set_p, flag_set<Fs...> flag_p) const {
        set_p.add(flag_p);
    }
};

}//namespace details
template<class ... Ts>
constexpr combined_flag_set to_opcode(Ts&&... ts_p){
    combined_flag_set result{};
    details::to_opcode_impl<Ts...>{}( result, std::move(ts_p)... );
    return result;
}

} //namespace sf_g

#endif /* flag_set_h */
