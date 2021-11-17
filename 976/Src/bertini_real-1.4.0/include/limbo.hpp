#pragma once

#include <map>

/**
 a templated function for looking up a value in a map by key, without accidentally creating a key-value pair when it didn't previously exist, and when it doesn't exist, gives back a default value.
 
 \return the value of the key, if it exists.  if nexists, returns default_value (which was input)
 \param mc_mapperson the map to look into.
 \param lookup_key the key to search for in the map.
 \param default_value If the map doesn't hold an entry for the lookup_key, returns this value.
 */
template <typename key_type, typename value_type>
value_type map_lookup_with_default(const  std::map <key_type,value_type> & mc_mapperson, const key_type & lookup_key, const value_type& default_value )
{
	typename std::map<key_type,value_type>::const_iterator it = mc_mapperson.find( lookup_key );
	if ( it == mc_mapperson.end() ) {
		return default_value;
	}
	else {
		return it->second;
	}
}


