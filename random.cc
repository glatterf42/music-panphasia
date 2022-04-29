/*

 random.cc - This file is part of MUSIC -
 a code to generate multi-scale initial conditions
 for cosmological simulations

 Copyright (C) 2010  Oliver Hahn

 */

#include "random.hh"

// TODO: move all this into a plugin!!!

std::map< std::string, RNG_plugin_creator *>&
get_RNG_plugin_map()
{
  static std::map< std::string, RNG_plugin_creator* > RNG_plugin_map;
  return RNG_plugin_map;
}

void print_RNG_plugins()
{
  std::map< std::string, RNG_plugin_creator *>& m = get_RNG_plugin_map();
  std::map< std::string, RNG_plugin_creator *>::iterator it;
  it = m.begin();
  std::cout << " - Available random number generator plug-ins:\n";
  while( it!=m.end() )
    {
      if( (*it).second )
	std::cout << "\t\'" << (*it).first << "\'\n";
      ++it;
    }
}

RNG_plugin *select_RNG_plugin( config_file& cf )//, const refinement_hierarchy& refh )
{
	std::string rngname = cf.getValueSafe<std::string>( "random", "generator", "MUSIC" );

	RNG_plugin_creator *the_RNG_plugin_creator
	= get_RNG_plugin_map()[ rngname ];

	if( !the_RNG_plugin_creator )
	{
		std::cerr << " - Error: random number generator plug-in \'" << rngname << "\' not found." << std::endl;
		LOGERR("Invalid/Unregistered random number generator plug-in encountered : %s",rngname.c_str() );
		print_RNG_plugins();
		throw std::runtime_error("Unknown random number generator plug-in");

	}else
	{
		std::cout << " - Selecting random number generator plug-in \'" << rngname << "\'..." << std::endl;
		LOGUSER("Selecting random number generator plug-in  : %s",rngname.c_str() );
	}

	RNG_plugin *the_RNG_plugin
	= the_RNG_plugin_creator->create( cf );//, refh );

	return the_RNG_plugin;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
#pragma mark -
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////


//template class random_numbers<real_t>;
//template class random_number_generator< random_numbers<real_t>, real_t >;

/*template class random_numbers<float>;
template class random_numbers<double>;
template class random_number_generator< random_numbers<float>, float >;
template class random_number_generator< random_numbers<double>, double >;*/
