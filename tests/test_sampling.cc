#include <catch/catch.hpp>
#include "model.hh"
#include "sampling.hh"

using namespace std;
using namespace addapt;

TEST_CASE("Correctly identify mutable positions", "[sampling]") {

	SECTION("one macrostate") {
		ConstructPtr construct = make_shared<Construct>("UUUuuu");
		construct->add_macrostate("a", "(.)(.)");

		CHECK(can_be_mutated(construct, 0) == true);
		CHECK(can_be_mutated(construct, 1) == true);
		CHECK(can_be_mutated(construct, 2) == true);
		CHECK(can_be_mutated(construct, 3) == false);
		CHECK(can_be_mutated(construct, 4) == false);
		CHECK(can_be_mutated(construct, 5) == false);

		CHECK(can_be_freely_mutated(construct, 0) == true);
		CHECK(can_be_freely_mutated(construct, 1) == true);
		CHECK(can_be_freely_mutated(construct, 2) == false);
		CHECK(can_be_freely_mutated(construct, 3) == false);
		CHECK(can_be_freely_mutated(construct, 4) == false);
		CHECK(can_be_freely_mutated(construct, 5) == false);
	}

	SECTION("multiple macrostates") {
		ConstructPtr construct = make_shared<Construct>("UUUU");
		construct->add_macrostate("a", "().)");
		construct->add_macrostate("b", "(.))");

		CHECK(can_be_mutated(construct, 0) == true);
		CHECK(can_be_mutated(construct, 1) == true);
		CHECK(can_be_mutated(construct, 2) == true);
		CHECK(can_be_mutated(construct, 3) == true);

		CHECK(can_be_freely_mutated(construct, 0) == true);
		CHECK(can_be_freely_mutated(construct, 1) == false);
		CHECK(can_be_freely_mutated(construct, 2) == false);
		CHECK(can_be_freely_mutated(construct, 3) == false);
	}

}

TEST_CASE("Correctly mutate base-paired positions", "[sampling]") {
	ConstructPtr construct;

	struct Test {
		string sequence;
		vector<string> macrostates;
		string mutations;
		vector<string> expected_sequences;
	};

	// For each test, a construct will be created with the sequence from the 
	// first field and the macrostates from the second field.  The third field 
	// specifies which mutations will be made to that construct.  For example, 
	// "AG" specifies two tests.  In the first, the first position will be 
	// mutated to "A".  In the second, the second position will be mutated to 
	// "G".  The fourth field contains the expected sequence for each test.

	vector<Test> tests = {
		{"N",    {"."},    "A",    {"A"}},
		{"NN",   {"()"},   "AG",   {"AU", "CG"}},
		{"NNN",  {"(.)"},  "AGU",  {"ANU", "NGN", "ANU"}},

		{"NN",   {"()",
		          "()"},   "AG",   {"AU", "CG"}},

		{"NNN",  {"().",
		          ".()"},  "AGU",  {"AUA", "CGC", "UAU"}},

		{"NNN",  {"().",
		          "(.)"},  "AGU",  {"AUU", "CGG", "AUU"}},

		{"NNNN", {"()..",
		          "(())"}, "AGUC", {"AUAU", "CGCG", "UAUA", "GCGC"}},

		{"NNNN", {"(.).",
		          "(())"}, "AGUC", {"AAUU", "GGCC", "AAUU", "GGCC"}},
	};

	for(Test test: tests) {
		for(int i = 0; i < test.mutations.length(); i++) {
			construct = make_shared<Construct>(test.sequence);
			CAPTURE(test.sequence)

			for(int x = 0; x < test.macrostates.size(); x++) {
				CAPTURE(test.macrostates[x]);
				construct->add_macrostate(to_string(x), test.macrostates[x]);
			}

			CAPTURE(test.mutations);
			mutate_recursively(construct, i, test.mutations[i]);
			CHECK(construct->seq() == test.expected_sequences[i]);
		}
	}

	// Check some error conditions.
	construct = make_shared<Construct>("Nn");
	construct->add_macrostate("not mutable", "()");
	CHECK_THROWS(mutate_recursively(construct, 0, 'G'));

	construct = make_shared<Construct>("N");
	construct->add_macrostate("extra open", "(");
	CHECK_THROWS(mutate_recursively(construct, 0, 'G'));

	construct = make_shared<Construct>("N");
	construct->add_macrostate("extra close", ")");
	CHECK_THROWS(mutate_recursively(construct, 0, 'G'));


}
