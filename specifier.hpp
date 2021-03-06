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
    static constexpr uint8_t shift = 4;
};

struct charge_flag{
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

struct mean_flag{
    static constexpr uint8_t shift = 5;
};

struct sigma_flag{
    static constexpr uint8_t shift = 6;
};

struct pile_up_flag{
    static constexpr uint8_t shift = 7;
};

template< class ... Ts >
struct flag_set{
    
public:
    constexpr operator uint8_t() const { return compute_value(); }
    
private:
    flag_set() = default;
    
    constexpr uint8_t compute_value() const {
        uint8_t result = 0;
        int expander[] = { 0, ( result |= 1UL << Ts::shift , 0) ... };
        return result;
        // return 1;
    }
};

template<>
struct flag_set<> {
    constexpr operator uint8_t() const {return 0;}
};


template< uint8_t Opcode > struct specifier{};
template<> struct specifier<0x01>{ using pack = details::pack< amplitude_finder>; using data_t = amplitude ; };
template<> struct specifier<0x02>{ using pack = details::pack< cfd_calculator> ; using data_t = cfd_time;  };
template<> struct specifier<0x04>{ using pack = details::pack< rise_time_calculator>; using data_t = rise_time; };
template<> struct specifier<0x08>{ using pack = details::pack< charge_integrator >; using data_t = charge; };
template<> struct specifier<0x10>{ using pack = details::pack< baseline_finder >; using data_t = baseline; };
template<> struct specifier<0x20>{ using pack = details::pack< mean_finder >; using data_t = mean; };
template<> struct specifier<0x80>{ using pack = details::pack< pile_up_finder >; using data_t = mean; };

template<> struct specifier<0x03>{ using pack = details::pack< amplitude_finder, cfd_calculator>; using data_t = composite<amplitude, cfd_time> ; };
template<> struct specifier<0x11>{ using pack = details::pack< amplitude_finder, baseline_finder >; using data_t = composite<amplitude, baseline>; };
template<> struct specifier<0x09>{ using pack = details::pack< amplitude_finder, charge_integrator >; using data_t = composite<amplitude, charge>; };
template<> struct specifier<0x18>{ using pack = details::pack< charge_integrator, baseline_finder >; using data_t = composite<charge, baseline>; };
template<> struct specifier<0x12>{ using pack = details::pack< cfd_calculator, baseline_finder >; using data_t = composite< cfd_time, baseline>; };

template<> struct specifier<0x19>{ using pack = details::pack< amplitude_finder, charge_integrator, baseline_finder >; using data_t = composite<amplitude, charge, baseline>; };
template<> struct specifier<0x1C>{ using pack = details::pack< baseline_finder, charge_integrator, rise_time_calculator >; using data_t = composite< baseline, charge, rise_time >; };
template<> struct specifier<0x1A>{ using pack = details::pack< charge_integrator, cfd_calculator, baseline_finder >; using data_t = composite< charge, cfd_time, baseline>; };
template<> struct specifier<0x13>{ using pack = details::pack< amplitude_finder, cfd_calculator, baseline_finder >; using data_t = composite< amplitude, cfd_time, baseline>; };
template<> struct specifier<0x91>{ using pack = details::pack< amplitude_finder, pile_up_finder, baseline_finder >; using data_t = composite<amplitude, pile_up, baseline>; };

template<> struct specifier<0x1B>{ using pack = details::pack< amplitude_finder, charge_integrator, cfd_calculator, baseline_finder >; using data_t = composite< amplitude, charge, cfd_time, baseline>; };
template<> struct specifier<0x2B>{ using pack = details::pack< amplitude_finder, charge_integrator, cfd_calculator, mean_finder >; using data_t = composite< amplitude, charge, cfd_time, mean>; };

template<> struct specifier<0x3B>{ using pack = details::pack< amplitude_finder, charge_integrator, cfd_calculator, mean_finder, baseline_finder >; using data_t = composite< amplitude, charge, cfd_time, mean, baseline>; };

template<> struct specifier<0x7B>{ using pack = details::pack< amplitude_finder, charge_integrator, cfd_calculator, mean_finder, sigma_finder, baseline_finder >; using data_t = composite< amplitude, charge, cfd_time, mean, sigma, baseline>; };

struct waveform_specifier{ using data_t = waveform; };

struct specifier_pairing{
    std::size_t channel_number;
    uint8_t opcode;
};    

specifier_pairing make_specifier_pairing( std::string& module_list_p ) {
    specifier_pairing sp{};
    for( auto& token : module_list_p ){ 
        switch(token){
        case 'a': { sp.opcode |= sf_g::flag_set< sf_g::amplitude_flag >{} ; break ; }
        case 't': { sp.opcode |= sf_g::flag_set< sf_g::cfd_flag >{} ; break ; }
        case 'b': { sp.opcode |= sf_g::flag_set< sf_g::baseline_flag >{} ; break ; }
        case 'c': { sp.opcode |= sf_g::flag_set< sf_g::charge_flag >{} ; break ; }
        case 'r': { sp.opcode |= sf_g::flag_set< sf_g::rise_flag >{} ; break ; }
        case 'm': { sp.opcode |= sf_g::flag_set< sf_g::mean_flag >{} ; break ; }
        case 's': { sp.opcode |= sf_g::flag_set< sf_g::sigma_flag >{} ; break ; }
        case 'p': { sp.opcode |= sf_g::flag_set< sf_g::pile_up_flag >{} ; break ; }
        case ':': { break; }           
        case '{': { break; }           
        case '}': { break; }           
        case '0': { sp.channel_number = 0 ; break; }
        case '1': { sp.channel_number = 1 ; break; }
        case '2': { sp.channel_number = 2 ; break; }
        case '3': { sp.channel_number = 3 ; break; }
        case '4': { sp.channel_number = 4 ; break; }
        case '5': { sp.channel_number = 5 ; break; }
        case '6': { sp.channel_number = 6 ; break; }
        case '7': { sp.channel_number = 7 ; break; }
        default : { std::cerr << "modifier " << token << " not implemented yet\n"; break ; }
        }
    }  
    return  sp;
}


} //namespace sf_g

#endif /* flag_set_h */
