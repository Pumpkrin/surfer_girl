#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include "data_format.hpp"
#include "reader.hpp"
#include "modifier.hpp"
#include "writer.hpp"

namespace sf_g{
template< class RO, 
          class MI, class MO, class ... Ms>
constexpr auto operator|( RO&& reader_output_p, modifier<MI, MO, Ms...>& m_p) {
   return m_p( std::move(reader_output_p) );
} 
template< class RO, 
          class MI, class MO, class ... Ms>
constexpr auto operator|( RO&& reader_output_p, modifier<MI, MO, Ms...> const& m_p) {
   return m_p( std::move(reader_output_p) );
} 


template< class MO,
          class WI, class WO>
constexpr void operator|(MO&& modifier_output_p, writer<WI, WO>& w_p) {
    return w_p( std::move(modifier_output_p) );
}
template< class MO,
          class WI, class WO>
constexpr void operator|(MO&& modifier_output_p, writer<WI, WO> const& w_p) {
    return w_p( std::move(modifier_output_p) );
}
} //namespace sf_g

#endif
