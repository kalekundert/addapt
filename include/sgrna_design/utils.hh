#pragma once

#include <string>
#include <memory>
#include <utility>

#include <boost/format.hpp>

namespace sgrna_design {

using f = boost::format;
using std::string;
using std::vector;
using std::pair;

enum class ColorEnum {
	NORMAL = 0,
	BLACK = 30,
	RED = 31,
	GREEN = 32,
	YELLOW = 33,
	BLUE = 34,
	MAGENTA = 35,
	CYAN = 36,
	WHITE = 37
};

enum class StyleEnum {
	NORMAL = 0,
	BOLD = 1,
	REVERSE = 2,
};


string color(string, ColorEnum, StyleEnum);

// weighted_choice?

}
