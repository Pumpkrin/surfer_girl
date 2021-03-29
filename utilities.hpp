#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include "data_format.hpp"
#include "reader.hpp"
#include "modifier.hpp"
#include "writer.hpp"

namespace sf_g{

template< class RO, 
          class S>
constexpr auto operator|( RO&& reader_output_p, modifier<S>& m_p) {
   return m_p( std::move(reader_output_p) );
} 
template< class RO, 
          class S>
constexpr auto operator|( RO&& reader_output_p, modifier<S> const& m_p) {
   return m_p( std::move(reader_output_p) );
} 

template< class RO > 
constexpr auto operator|( RO&& reader_output_p, raw_modifier& m_p) {
   return m_p( std::move(reader_output_p) );
} 
template< class RO >
constexpr auto operator|( RO&& reader_output_p, raw_modifier const& m_p) {
   return m_p( std::move(reader_output_p) );
} 

template< class MO,
          class S>
constexpr void operator|(MO&& modifier_output_p, writer<S>& w_p) {
    return w_p( std::move(modifier_output_p) );
}
template< class MO,
          class S>
constexpr void operator|(MO&& modifier_output_p, writer<S> const& w_p) {
    return w_p( std::move(modifier_output_p) );
}

template< class RO >
constexpr void operator|(RO const& reader_output_p, tree_editor const& te_p) {
    return te_p( reader_output_p );
}
template< class RO >
constexpr void operator|(RO const& reader_output_p, tree_editor& te_p) {
    return te_p( reader_output_p );
}

template< class RO >
constexpr void operator|(RO&& reader_output_p, raw_writer const& w_p) {
    return w_p( std::move(reader_output_p) );
}
template< class RO >
constexpr void operator|(RO&& reader_output_p, raw_writer& w_p) {
    return w_p( std::move(reader_output_p) );
}
} //namespace sf_g

#endif
