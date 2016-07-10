#pragma once

#include <iostream>
#include <iterator>
#include <list>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <boost/format.hpp>

namespace sgrna_design {

using f = boost::format;
using std::make_shared;
using std::string;
using std::vector;
using std::list;
using std::pair;

/// @brief Do the indices refer to the items in the collection themselves, or 
/// to the spaces between the items?
enum class IndexEnum {
	ITEM,
	BETWEEN,
};

int
normalize_index(string, int, IndexEnum);

pair<int,int>
normalize_range(string, int, int, IndexEnum);

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

string
color(string, ColorEnum, StyleEnum);


}

namespace std {

template <typename T>
ostream& operator<< (ostream& out, const vector<T>& vec) {
  if ( !vec.empty() ) {
    out << '[';
    copy(vec.begin(), vec.end(), ostream_iterator<T>(out, ", "));
    out << "\b\b]";
  }
  return out;
}


}

