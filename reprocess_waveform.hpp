#ifndef REPROCESS_WAVEFORM_HPP
#define REPROCESS_WAVEFORM_HPP

#include <fstream>

#include "TFile.h"
#include "TH1.h"


namespace sf_g{
template<class T, class Enabler = void> struct is_member_of_cut_tag : std::false_type{};
template<class T> struct is_member_of_cut_tag< T, decltype( adl_is_member_of_cut_tag( std::declval<T>() ) )> : std::true_type{};

namespace cut_tag{    
template<class T>
constexpr void adl_is_member_of_cut_tag(T &&);

template<std::size_t Index>
struct amplitude { 
    constexpr static const char name[] = "a"; 
    constexpr static std::size_t index = Index;
    template<class Input>
    constexpr double retrieve_value(std::vector<Input*> const& input_pch){ 
        return input_pch[Index]->amplitude;   
    }
};
template<std::size_t I> constexpr const char amplitude<I>::name[];

template<std::size_t Index>
struct baseline { 
    constexpr static const char name[] = "b"; 
    constexpr static std::size_t index = Index;
    template<class Input>
    constexpr double retrieve_value(std::vector<Input*> const& input_pch){ 
        return input_pch[Index]->baseline;   
    }
};
template<std::size_t I>
constexpr const char baseline<I>::name[];

template<std::size_t I1, std::size_t I2>
struct cfd_time { 
    constexpr static const char name[] = "t";
    constexpr static std::size_t index = I1 *10 + I2; //not so working hack -> 0 dissapears (of course)
    template<class Input>
    constexpr double retrieve_value(std::vector<Input*> const& input_pch){ 
        return input_pch[I2]->time - input_pch[I1]->time;   
    }
};
template<std::size_t I1, std::size_t I2>
constexpr const char cfd_time<I1,I2>::name[];

template<std::size_t I>
struct charge { constexpr static const char name[] = "c";};    
template<std::size_t I>
constexpr const char charge<I>::name[];
template<std::size_t I>
struct rise_time {constexpr static const char name[] = "r";};    
template<std::size_t I>
constexpr const char rise_time<I>::name[];
template<std::size_t I>
struct fall_time {constexpr static const char name[] = "f";};    
template<std::size_t I>
constexpr const char fall_time<I>::name[];
template<std::size_t Index>
struct mean { 
    constexpr static const char name[] = "m"; 
    constexpr static std::size_t index = Index;
    template<class Input>
    constexpr double retrieve_value(std::vector<Input*> const& input_pch){ 
        return input_pch[Index]->mean;   
    }
};
template<std::size_t I> constexpr const char mean<I>::name[];
template<std::size_t Index>
struct sigma { 
    constexpr static const char name[] = "s"; 
    constexpr static std::size_t index = Index;
    template<class Input>
    constexpr double retrieve_value(std::vector<Input*> const& input_pch){ 
        return input_pch[Index]->sigma;
    }
};
template<std::size_t I> constexpr const char sigma<I>::name[];
} //namespace tag

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
constexpr const char exclude_inner::name[];

template<class Type, class Input, class Selection = exclude_outer>
struct cut{
     constexpr bool operator()(std::vector< Input*> const& input_pch) const{
//        puts( __PRETTY_FUNCTION__ );
        double value = Type{}.template retrieve_value(input_pch);
        return Selection{}(lower_limit_m, value, upper_limit_m);
    }

    constexpr cut( double lower_limit_p, double upper_limit_p ) : lower_limit_m{lower_limit_p}, upper_limit_m{ upper_limit_p }{} 
    std::string name() const{ 
        std::string result = std::string{"_ch"} + std::to_string(Type::index) + "_" + Selection::name + Type::name + "]";  
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
        virtual bool operator()(std::vector<I*> const&) const = 0;
        virtual std::string name() const = 0;
        virtual ~eraser()=default;
    };

    template<class T>
    struct holder : eraser{
        constexpr holder() = default;
        constexpr holder(T&& t_p) : t_m{std::move(t_p)}{}
        bool operator()(std::vector<I*> const& input_pch) const override { return t_m(input_pch); }
        std::string name() const override { return t_m.name(); };
        T t_m;
    };

    constexpr cut_vector() { erased_mch.reserve(5); input_mch.reserve(8); redirected_input_mch.reserve(8); }
    
    //T becomes operator -> predicate ?
    template<class T, class S = exclude_outer>
    constexpr void add_cut( double lower_limit_p, double upper_limit_p ){
        static_assert( sf_g::is_member_of_cut_tag< T >::value, "this function should be used with a cut class : i.e. cut_tag::*****" );
        using cut_t = cut<T, I, S>;
        erased_mch.emplace_back( new holder<cut_t>{ cut_t{ lower_limit_p, upper_limit_p } }  );
    }

    void setup_tree(TTree * tree_ph) {
        for( auto i{0}; i < 8 ; ++i ){
            std::string branch_name = "channel_" + std::to_string(i) + ".";
            input_mch.emplace_back( new I{} );
            redirected_input_mch.emplace_back( input_mch.back().get() );
            if( tree_ph->GetBranch(branch_name.c_str()) ){ 
                std::cout << "setting_branch: " << branch_name << '\n'; 
                tree_ph->SetBranchAddress(branch_name.c_str() , &redirected_input_mch[i]);
            }
        }
    }
    //need to pass in tree pass to cut::whatever< S >, cut::time_distribution< I1, I2 >
    // should not go out after first cut
    constexpr bool operator()() const {
//        puts( __PRETTY_FUNCTION__ );
//        for( auto const& cut_h : erased_mch ){ auto is_in_cut = !(*cut_h)(redirected_input_mch); std::cout << "bool: " << std::boolalpha << is_in_cut << '\n'; if( is_in_cut ){ return false;} }
        for( auto const& cut_h : erased_mch ){ if( !(*cut_h)(redirected_input_mch) ){ return false;} }
        return true;
    }

    std::string name() const {
        std::string result;
        for( auto const& cut_h : erased_mch){ result +=  cut_h->name();  }
        return result;
    }

    void clear() { erased_mch.clear(); }
    private:
    std::vector< std::unique_ptr<eraser> > erased_mch;
    std::vector< std::unique_ptr<I> > input_mch;
    std::vector< I* > redirected_input_mch;
};


}//namespace sf_g

//specification format is : path/to/file.root:branch_name.
//add channel considered to cut vector, that way no modification is needed or even add it to cut ?
template< class I >
void reprocess_waveform( std::string const& file_p, sf_g::cut_vector<I>& cv_p ) {
    auto get_part_l = []( std::string const& input_p , std::regex regex_p ){
                            std::smatch result;
                            std::regex_search( input_p, result, regex_p);
                            return result[0].str();
                                                              };

    TFile file( file_p.c_str() );
    TTree* tree_h = static_cast<TTree*>( file.Get("data") );
    cv_p.setup_tree( tree_h );
    
    auto output_repertory = get_part_l( file_p, std::regex{"\\w+(?=/)"} ) + "/cut/";
    auto output = get_part_l(file_p, std::regex{"[^/]+(?=\\.root)"} ) + cv_p.name() + ".cut" ;

    std::ofstream sink{ output_repertory + output, std::ios::binary }; 
    std::cout << "before_cut: " << tree_h->GetEntries() << '\n';
    std::size_t counter{0};
    for(auto i{0}; i < tree_h->GetEntries(); ++i){
//        std::cout << "applying cut\n";
        tree_h->GetEntry(i);
        if( cv_p() ){ ++counter; sink.write( reinterpret_cast<char*>(&i), sizeof(i)); } 
    }
    std::cout<< "after_cut: " << counter << '\n';
} 

#endif

