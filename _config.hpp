// For private use only.
// Copyright (c) Manuel Freiberger, 2016.

#ifndef RS11_CONFIG_HPP
#define RS11_CONFIG_HPP

#include <rs11_user_config.hpp>

// Check the version of the user configuration file.
#if RS11_USER_CONFIG_VERSION != 1
    #error "Version 1 of the RS11 user configuration is required."
#endif // RS11_USER_CONFIG_VERSION

#endif // RS11_CONFIG_HPP
