// -----------------------------------------------------------------------
// This file is distributed under a 2-clause BSD like license
// Alternatively this file can be used under GPLv2 license
// See LICENSE.TXT for details.

#ifndef RS11_CONFIG_HPP
#define RS11_CONFIG_HPP

#include <rs11_user_config.hpp>

// Check the version of the user configuration file.
#if RS11_USER_CONFIG_VERSION != 1
    #error "Version 1 of the RS11 user configuration is required."
#endif // RS11_USER_CONFIG_VERSION

#endif // RS11_CONFIG_HPP
