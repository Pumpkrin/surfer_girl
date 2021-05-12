#ifndef WRITER_HPP
#define WRITER_HPP

#include "data_structure.hpp"
#include "data_stream.hpp"
#include "modifier.hpp"
#include "TTree.h"

namespace sf_g{

    //waveform_writer rather
struct raw_writer{
private:
    struct linker{
        linker( int channel_count_p ) : data_mc( channel_count_p ) { channel_number_mc.reserve( channel_count_p ); }
        waveform* add( std::size_t channel_number_p ){ 
            channel_number_mc.push_back( channel_number_p ); 
            return (*this)( channel_number_p ); 
        }
        waveform* operator()( std::size_t channel_number_p ){ 
            auto index = std::distance( 
                        channel_number_mc.begin(),
                        std::find( channel_number_mc.begin(), channel_number_mc.end(), channel_number_p )
                                       );
            return &data_mc[index];
        }
    private:
        std::vector<std::size_t> channel_number_mc;
        std::vector< waveform > data_mc;
    };

    void link_branches( std::vector<linked_waveform> && input_pc ){
//        for( auto && input : input_pc) {
//            auto * data_h = linker_m.add( input.channel_number );
//            std::string name = "channel_" + std::to_string(input.channel_number) + ".";
//            std::cout << name << '\n';
//            tree_m.branch( name.c_str(), data_h ); 
//            linker_m( input.channel_number )->data = std::move(input.data);
//        } 
        for( auto i{0}; i < input_pc.size() ; ++i) {
            auto * data_h = linker_m.add( i );
            std::string name = "channel_" + std::to_string(i) + ".";
            std::cout << name << '\n';
            tree_m.branch( name.c_str(), data_h ); 
            linker_m( i )->data = std::move(input_pc[i].data);
        } 
        current_function_h = &raw_writer::fill;
    };
    void fill( std::vector<linked_waveform> && input_pc) {
//        for( auto && input : input_pc) {
//           linker_m( input.channel_number )->data = std::move( input.data ) ;
//        }
        for( auto i{0}; i < input_pc.size() ; ++i ) {
           linker_m( i )->data = std::move( input_pc[i].data ) ;
        }
    };

public:
    raw_writer(data_output<TTree>& tree_p, int channel_count_p  ) : 
        current_function_h{ &raw_writer::link_branches },
        channel_count_m{ static_cast<size_t>( channel_count_p ) }, 
        linker_m{ channel_count_p },
        tree_m{ tree_p } {
            tree_m.register_writer();
        }

    constexpr void operator()( std::vector<linked_waveform>&& input_pc ) {
        if( !input_pc.empty() ){ 
            (this->*current_function_h)( std::move(input_pc) );
            tree_m.fill();
        }
    } 

    private:
    void (raw_writer::* current_function_h)( std::vector<linked_waveform> && input_p );
    size_t const channel_count_m;
    linker linker_m;
    data_output<TTree>&  tree_m;
};


template<class Specifier>
struct writer {
    using data_t = typename Specifier::data_t;

    writer( std::size_t associated_channel_p, data_output<TTree>& output_p ) : 
        tree_m{ output_p },
        data_mh{ new data_t{} }
    { 
        std::string name = "channel_" + std::to_string( associated_channel_p) +  ".";
        tree_m.register_branch_writer( name, data_mh.get() );
//        tree_m.branch( "test_" + std::to_string( associated_channel_p) ,&counter );
    }

public:
    constexpr void operator()( data_t&& input_p ) {
        *data_mh = input_p;
//        data_mh->value();
        tree_m.fill();
    } 

private:
    data_output<TTree>& tree_m;
    std::unique_ptr<data_t> data_mh; 
};


// ------------------------------------branch_editor-------------------------------------------

template<class Specifier>
struct branch_editor {
    branch_editor( std::size_t associated_channel_p, data_output<TTree>& sink_p) :
        associated_channel_m{associated_channel_p},
        w_m{ associated_channel_m, sink_p } {} 
    void operator()( std::vector<linked_waveform> const& input_pc ) {
        m_m( input_pc[associated_channel_m] ) | w_m;
    }
    std::size_t associated_channel() const { return associated_channel_m; }

private:
    std::size_t const associated_channel_m;
    
    modifier<Specifier> const m_m;
    writer<Specifier> w_m;
};

struct editor {
    struct eraser {
        virtual ~eraser() = default;
        virtual void apply_yourself( std::vector< linked_waveform > const& data_pc ) = 0; 
        virtual std::size_t associated_channel() const = 0 ;
    };
    template<class T>
    struct holder : eraser {
        constexpr holder() = default;
        constexpr holder( T t_p ) : t_m{ std::move(t_p)} {}
        void apply_yourself( std::vector< linked_waveform > const& data_pc ) override { t_m( data_pc ); }
        std::size_t associated_channel() const override { return t_m.associated_channel(); }
        T t_m;
    };

    constexpr editor() = default;
    editor( editor const& ) = delete;
    editor( editor&& ) = default;

    template< class T>
    constexpr editor(T t_p) : erased_mh{ new holder<T>{ std::move( t_p ) } } {}

    void operator()(std::vector<linked_waveform> const& data_pc){ erased_mh->apply_yourself( data_pc ); }
    std::size_t associated_channel() const { return erased_mh->associated_channel(); }

    private: 
    std::unique_ptr<eraser> erased_mh;
}; 

//--------------------------------tree_writer------------------------
struct tree_editor{
    void operator()(std::vector<linked_waveform> const&  input_pc){ 
        for( auto& e : e_mc){ e(input_pc); }
    }
    void reserve( std::size_t size_p ) { e_mc.reserve( size_p ); }
    editor& add( editor&& be_p ){
        if( unpaired_branch( be_p ) ){ e_mc.push_back( std::move(be_p) ); return e_mc.back(); }
        std::cerr << "This branch has been paired up with an other writer already\n";
        return *std::find_if( e_mc.begin(), e_mc.end(),
                             [&be_p](editor const& saved_editor_p){ return be_p.associated_channel() == saved_editor_p.associated_channel(); } );
    }
    private:
    bool unpaired_branch( editor const& be_p ) {
        return std::find_if( e_mc.begin(), e_mc.end(),
                             [&be_p](editor const& saved_editor_p){ return be_p.associated_channel() == saved_editor_p.associated_channel(); } )
                         == e_mc.end();
    }
private:
    std::vector<editor> e_mc;
};

} //namespace sf_g




#endif
