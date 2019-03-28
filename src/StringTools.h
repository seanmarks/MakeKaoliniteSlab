// Wrapper class with static methods for handling strings

#ifndef STRING_TOOLS_H
#define STRING_TOOLS_H

#include <algorithm>
#include <array>
#include <exception>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

class StringTools 
{
 public:
	// Prints a fixed-width string to the indicated stream
	// - Assumes any desired formatting flags have already been set (e.g. by std::fixed)
	// - The string is truncated as necessary according to the stream's 'adjustfield' flag
	//   - std::left     - trim from the right
	//   - std::right    - trim from the left
	//   - std::internal - unsupported (throws an exception)
	//   - (other)       - same as std::right
	// - If 'width' is non-negative, calls os.setw(width)
	//   - Else uses whatever value (if any) is already present in the stream
	template<typename T, typename C = char>
	static void printFixedWidthString(
		std::basic_ostream<C>& os,  // use for determining flags
		const T&   val,
		const int  width = -1   // default: don't limit width
	)
	{
		// Create a stringstream with the same formatting flags
		std::ios_base::fmtflags flags = os.flags();
		std::stringstream ss;
		ss.flags(flags);

		// Convert val to string
		ss << val;
		std::string raw_str( ss.str() );

		if ( width >= 0 ) {
			os << std::setw(width);
		}

		// Keep only part of string
		int len = raw_str.size();
		if ( width >= 0 and len > width ) {
			// Get 'adjustfield' flag, which determines text alignment
			std::ios_base::fmtflags adjust_field_flag = (flags & std::ios_base::adjustfield);

			// First character to keep
			int first = 0;
			if ( adjust_field_flag == std::ios_base::left ) {  // left-aligned
				first = 0;
			}
			else if ( adjust_field_flag == std::ios_base::internal ) {  // "internal" justification
				throw std::runtime_error("printFixedWidthString() doesn't support adjustfield 'internal'");
			}
			else {  // right-aligned (default)
				first = len - width;
			}

			os << raw_str.substr(first, width);
		}
		else {
			os << raw_str;
		}
	};

	// Convert a std::string to a value of the given type
	template<typename T> 
	static void stringToValue(const std::string& str, T& value) {
		std::stringstream ss( str );
		ss >> value;

		if ( ss.fail() ) {
			std::stringstream err_ss;
			err_ss << "unable to convert " << str << " to a number\n";
			throw std::runtime_error( err_ss.str() );
		}
	}

	template<typename T, std::size_t dim>
	static void stringsToArray( const std::vector<std::string>& strings, std::array<T,dim>& arr ) {
		if ( strings.size() != dim ) {
			throw std::runtime_error("size of vector of strings does not match output array size");
		}

		for ( unsigned i=0; i<dim; ++i ) {
			stringToValue(strings[i], arr[i]);
		}
	}

	template<typename T>
	static void stringsToVector( const std::vector<std::string>& strings, std::vector<T>& vec ) {
		// Ensure the output vector has the correct size
		int len = strings.size();
		vec.resize(len);

		for ( int i=0; i<len; ++i ) {
			stringToValue(strings[i], vec[i]);
		}
	}

	// Converts a string containing "yes/no", "true/false", etc. to the appropriate bool
	static bool stringToBool(const std::string& str)
	{
		std::string lowercase_str = StringTools::toLowercase(str);

		if ( lowercase_str == "yes" || lowercase_str == "true" || 
				 lowercase_str == "1"   || lowercase_str[0] == 'y' ) { 
			return true; 
		}
		else if ( lowercase_str == "no" || lowercase_str == "false" ||
							lowercase_str == "0"  || lowercase_str[0] == 'n' ) { 
			return false;
		}
		else {
			std::string errorMsg = "Error in StringTools::stringToBool - Could not interpret = \"" 
														 + str + "\" as true/yes/1 or false/no/0. Check your input.\n";
			throw std::runtime_error(errorMsg);
		}
	}

	// Convert a string to all lowercase
	static std::string toLowercase(const std::string& str) {
		std::string lowercase_str = str;
		std::transform(lowercase_str.begin(), lowercase_str.end(), lowercase_str.begin(), ::tolower);
		return lowercase_str;
	}
};

#endif /* STRING_TOOLS_H */
