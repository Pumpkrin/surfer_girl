#ifndef REPROCESS_WAVEFORM_HPP
#define REPROCESS_WAVEFORM_HPP

#include <fstream>

#include "TFile.h"
#include "TH1.h"


namespace sf_g{
struct tag{};    
struct amplitude_tag : tag { constexpr static const char name[] = "a"; };
constexpr const char amplitude_tag::name[];
struct baseline_tag : tag { constexpr static const char name[] = "b"; };
constexpr const char baseline_tag::name[];
struct cfd_time_tag : tag { constexpr static const char name[] = "t";};
constexpr const char cfd_time_tag::name[];
struct charge_tag : tag { constexpr static const char name[] = "c";};    
constexpr const char charge_tag::name[];
struct rise_time_tag : tag {constexpr static const char name[] = "r";};    
constexpr const char rise_time_tag::name[];
struct fall_time_tag : tag {constexpr static const char name[] = "f";};    
constexpr const char fall_time_tag::name[];

template<class T, class I> struct retrieve_value_impl{};
template<class I> 
struct retrieve_value_impl<baseline_tag, I>{
    constexpr double operator()(I const& input_p){
        return input_p.baseline;
    }
};
template<class I> 
struct retrieve_value_impl<charge_tag, I>{
    constexpr double operator()(I const& input_p){
        return input_p.charge;
    }
};

template<class T, class I> double retrieve_value( I const& input_p ) {
    static_assert( std::is_base_of< tag, T >::value, "this function should be used with a tag class : i.e. *****_tag" );
    return retrieve_value_impl<T, I>{}( input_p );
}

struct exclude_outer{
    constexpr bool operator()( 
        double const& lower_limit_p,
        double const& value_p,
        double const& upper_limit_p
                             ) {
        return lower_limit_p < value_p && value_p < upper_limit_p;
    }

    constexpr static const char name[] = "";
};
constexpr const char exclude_outer::name[];

struct exclude_inner{
    constexpr bool operator()( 
        double const& lower_limit_p,
        double const& value_p,
        double const& upper_limit_p
                             ) {
        return value_p < lower_limit_p || upper_limit_p < value_p;
    }
    constexpr static const char name[] = "!";
};

template<class T, class I, class S = exclude_outer>
struct cut{
    constexpr bool operator()(I const& input_p) const{
        double value = retrieve_value<T>(input_p);
        return value > lower_limit_m && value < upper_limit_m;
    }

    constexpr cut( double lower_limit_p, double upper_limit_p ) : lower_limit_m{lower_limit_p}, upper_limit_m{ upper_limit_p }{} 
    std::string name() const{ 
        std::string result = std::string{"_"} + S::name + T::name + "]";  
        auto const size = static_cast<std::size_t>( std::snprintf( nullptr, 0, "%.2e-%.2e[", lower_limit_m, upper_limit_m ) + 1);
        char buffer[ size ];
        std::snprintf( buffer, size, "%.2e-%.2e[", lower_limit_m, upper_limit_m );
        result += std::string{buffer, size -1};
        return result; 
    }

private:
    double const lower_limit_m;
    double const upper_limit_m;
};

template<class I>
struct cut_vector{
    struct eraser{
        virtual bool operator()(I const& input_p) const = 0;
        virtual std::string name() const = 0;
        virtual ~eraser()=default;
    };

    template<class T>
    struct holder : eraser{
        constexpr holder() = default;
        constexpr holder(T&& t_p) : t_m{std::move(t_p)}{}
        bool operator()(I const& input_p) const override { return t_m(input_p); }
        std::string name() const override { return t_m.name(); };
        T t_m;
    };

    constexpr cut_vector() = default;
    
    template<class T, class S = exclude_outer>
    constexpr void add_cut( double lower_limit_p, double upper_limit_p ){
        using cut_t = cut<T, I, S>;
        erased_mch.emplace_back( new holder<cut_t>{ cut_t{ lower_limit_p, upper_limit_p } }  );
    }

    constexpr bool operator()(I const& input_p) const {
        for( auto const& cut_h : erased_mch ){ if( (*cut_h)(input_p) ){return true;} }
        return false;
    }

    std::string name() const {
        std::string result;
        for( auto const& cut_h : erased_mch){ result +=  cut_h->name();  }
        return result;
    }

    private:
    std::vector< std::unique_ptr<eraser> > erased_mch;
};


}//namespace sf_g

//specification format is : path/to/file.root:branch_name.
template< class I >
void reprocess_waveform( std::string const& specification_p, sf_g::cut_vector<I> const& cv_p ) {
    auto get_part_l = []( std::string const& input_p , std::regex regex_p ){
                            std::smatch result;
                            std::regex_search( input_p, result, regex_p);
                            return result[0].str();
                                                              };
    auto file_regex = std::regex{"[^:]+"}; 
    auto branch_regex = std::regex{"[^:]+\\.$"};

    TFile file( get_part_l( specification_p, std::move(file_regex) ).c_str() );
    TTree* tree_h = static_cast<TTree*>( file.Get("data") );
    auto* composite_h = new I{};
    auto branch_name = get_part_l( specification_p, std::move(branch_regex));
    tree_h->SetBranchAddress( branch_name.c_str() , &composite_h);
    
    auto output_repertory = get_part_l( specification_p, std::regex{"\\w+(?=/)"} ) + "/cut/";
    auto channel_number = get_part_l(branch_name, std::regex{"[0-9](?=\\.$)"});
    auto output = get_part_l(specification_p, std::regex{"[^/]+(?=\\.root)"} ) + "_ch" + channel_number + cv_p.name() + ".cut" ;
    std::ofstream sink{ output_repertory + output, std::ios::binary }; 
    for(auto i{0}; i < tree_h->GetEntries(); ++i){
        tree_h->GetEntry(i);
        if( cv_p( *composite_h ) ){ sink.write( reinterpret_cast<char*>(&i), sizeof(i)); } 
    }

    delete composite_h;
} 

#endif

