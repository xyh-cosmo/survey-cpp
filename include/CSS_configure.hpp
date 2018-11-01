//	This class holds all the relevant configuration parameters that are needed
//	to perform the survey simulation. This class will be called by ALL other
//  modules to get parameters to be used!
//
//  Created @ 2018-11-1
//  Author: YOUHUA XU

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "dictionary.h"
#include "iniparser.h"

#ifndef _CSS_CONFIGURE_HPP_
#define _CSS_CONFIGURE_HPP_

class CSS_Configure{
    dictionary *ini_dict;
public:

//  configuration parameters:
    double survey_start_time;
    double survey_end_time;

    CSS_Configure();
    ~CSS_Configure();

    void read_configure( std::string& config_file_name );
    void save_configure( std::string& config_file_name );

    int get_int( std::string& key );
    double get_double( std::string& key );
    std::string get_string( std::string& key );
};


#endif //_CSS_CONFIGURE_HPP_