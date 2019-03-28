#include "StringTools.h"

// Template specialization for std::strings - no conversion necessary
template<>
void StringTools::stringToValue<std::string>(const std::string& str, std::string& value)
{
	value = str;
}
