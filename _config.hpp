// -----------------------------------------------------------------------
// This file is distributed under a 2-clause BSD like license
// Alternatively this file can be used under GPLv2 license
// See LICENSE.TXT for details.

#ifndef RS11_CONFIG_HPP
#define RS11_CONFIG_HPP

#if __cplusplus >= 201403L ||  _MSC_VER >= 1900
    #include  <gsl/span>
    #define RS11_HAVE_GSL
#endif

#ifndef __cpp_exceptions
    #define RS11_ERR_NOEXCEPT
#else
    #define RS11_ERR_NOEXCEPT   noexcept
#endif

#endif // RS11_CONFIG_HPP
