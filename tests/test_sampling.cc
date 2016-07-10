#include <catch/catch.hpp>
#include "model.hh"
#include "sampling.hh"

using namespace std;
using namespace sgrna_design;

TEST_CASE("Test the AutoPointMutation move", "[sampling]") {
	AutoPointMutation move;

	ConstructPtr construct = make_shared<Construct>();
	*construct += make_shared<Domain>("a", "NNn", "(..");
	*construct += make_shared<Domain>("b", "nNN", "..)");

	std::mt19937 rng(1);

	move.apply(construct, rng); CHECK(construct->seq() == "NUnnNN");
	move.apply(construct, rng); CHECK(construct->seq() == "NUnnUN");
	move.apply(construct, rng); CHECK(construct->seq() == "AUnnUU");
	move.apply(construct, rng); CHECK(construct->seq() == "UUnnUA");
	move.apply(construct, rng); CHECK(construct->seq() == "AUnnUU");
	move.apply(construct, rng); CHECK(construct->seq() == "CUnnUG");
	move.apply(construct, rng); CHECK(construct->seq() == "CUnnUG");
	move.apply(construct, rng); CHECK(construct->seq() == "CGnnUG");
	move.apply(construct, rng); CHECK(construct->seq() == "CUnnUG");
}

