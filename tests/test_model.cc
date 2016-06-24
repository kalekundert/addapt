#include <catch/catch.hpp>
#include "model.hh"

using namespace sgrna;

TEST_CASE("Test model constructor", "[model]") {
	Domain domain(
			"aptamer",
			"AUACCAGCCGAAAGGCCCUUGGCAG",
			ColorEnum::YELLOW,
			StyleEnum::BOLD);

	REQUIRE(domain.name() == "aptamer");
	REQUIRE(domain.seq() == "AUACCAGCCGAAAGGCCCUUGGCAG");
	REQUIRE(domain.color() == ColorEnum::YELLOW);
	REQUIRE(domain.style() == StyleEnum::BOLD);
}
