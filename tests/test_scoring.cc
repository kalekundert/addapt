#include <set>
#include <catch/catch.hpp>
#include <sgrna_design/model.hh>
#include <sgrna_design/scoring.hh>

using namespace std;
using namespace sgrna_design;
using bp = pair<int,int>;

TEST_CASE("Test folding without an aptamer", "[scoring]") {
	// Make a construct that should fold into a hairpin.
	ConstructPtr hairpin = make_shared<Construct>();
	*hairpin += make_shared<Domain>("stem/a", "ACGU");
	*hairpin += make_shared<Domain>("loop", "GAAA");
	*hairpin += make_shared<Domain>("stem/b", "ACGU");

	// Indicate which base pairs I expect to find.
	set<bp> expected_base_pairs = {{0,11}, {1,10}, {2,9}, {3,8}};
	map<bp, double> base_pair_probs = {
		{{0,11}, 0.70},
		{{1,10}, 0.95},
		{{2, 9}, 0.95},
		{{3, 8}, 0.95}
	};

	RnaFold fold(hairpin);

	SECTION("the expected base pairs form") {
		for(int i = 0; i < hairpin->len(); i++) {
			for(int j = i; j < hairpin->len(); j++) {
				double p = fold.base_pair_prob(i, j);

				if(expected_base_pairs.count({i,j})) {
					double threshold = base_pair_probs[{i,j}];
					CHECK(p > threshold);
				}
				else {
					REQUIRE(p < 0.1);
				}
			}
		}
	}

	SECTION("the order of the indices doesn't matter") {
		for(int i = 0; i < hairpin->len(); i++) {
			for(int j = i; j < hairpin->len(); j++) {
				double ij = fold.base_pair_prob(i, j);
				double ji = fold.base_pair_prob(j, i);
				REQUIRE(ij == ji);
			}
		}
	}

	SECTION("out-of-bounds indices throw exceptions") {
		CHECK_THROWS(fold.base_pair_prob(0, 12));
		CHECK_THROWS(fold.base_pair_prob(0, -13));
	}
}

TEST_CASE("Test folding with an aptamer", "[scoring]") {
	// Make a construct consisting entirely of an aptamer.
	ConstructPtr aptamer = make_shared<Construct>();
	*aptamer += make_shared<Domain>("aptamer", "GAUACCAGCCGAAAGGCCCUUGGCAGC");

	// Indicate which base pairs I expect to form.
	set<bp> no_lig_base_pairs = {
		{4,22}, {5,21}, {6,20}, {7,16}, {8,15}, {9,14}
	};
	set<bp> lig_base_pairs = {
		{0,26}, {4,22}, {5,21}, {7,16}, {8,15}, {9,14}
	};

	RnaFold no_lig_fold(aptamer, LigandEnum::NONE);
	RnaFold lig_fold(aptamer, LigandEnum::THEO);

	SECTION("the expected base pairs form") {
		for(int i = 0; i < aptamer->len(); i++) {
			for(int j = i; j < aptamer->len(); j++) {
				double no_lig_bp = no_lig_fold.base_pair_prob(i,j);
				double lig_bp = lig_fold.base_pair_prob(i,j);

				if(no_lig_base_pairs.count({i,j})) {
					CHECK(no_lig_bp > 0.7);
				}
				else {
					REQUIRE(no_lig_bp < 0.3);
				}
				if(lig_base_pairs.count({i,j})) {
					CHECK(lig_bp > 0.7);
				}
				else {
					REQUIRE(lig_bp < 0.3);
				}
			}
		}
	}
}

