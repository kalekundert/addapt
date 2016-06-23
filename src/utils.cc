#include "utils.hh"

namespace sgrna {

string
color(
		string str,
		ColorEnum color=ColorEnum::NORMAL,
		StyleEnum style=StyleEnum::NORMAL) {

	auto fmt = boost::format("\033[%d;%dm%s\033[0;0m")
		% static_cast<int>(style)
		% static_cast<int>(color)
		% str;

	return fmt.str();
}

}
