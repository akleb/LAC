/**
 * @file lac_Version.hpp
 * @author Alex Kleb (akleb@umich.edu)
 * @brief defines the function that holds the current compiled version
 * @version 0.1
 * @date 2023-06-09
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef _LAC_VERSION_
#define _LAC_VERSION_

#include <string>

#ifndef LAC_GIT_VERSION
#define LAC_GIT_VERSION "DID NOT FIND VERSION FROM GIT"
#endif

std::string lac_Version();

#endif
