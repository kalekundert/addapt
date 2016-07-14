#include <catch/catch.hpp>
#include "model.hh"
#include "sampling.hh"

using namespace std;
using namespace addapt;

TEST_CASE("Correctly identify mutable positions", "[sampling]") {

	SECTION("one macrostate") {
		DevicePtr device = make_shared<Device>("UUUuuu");
		device->add_macrostate("a", "(.)(.)");

		CHECK(can_be_mutated(device, 0) == true);
		CHECK(can_be_mutated(device, 1) == true);
		CHECK(can_be_mutated(device, 2) == true);
		CHECK(can_be_mutated(device, 3) == false);
		CHECK(can_be_mutated(device, 4) == false);
		CHECK(can_be_mutated(device, 5) == false);

		CHECK(can_be_freely_mutated(device, 0) == true);
		CHECK(can_be_freely_mutated(device, 1) == true);
		CHECK(can_be_freely_mutated(device, 2) == false);
		CHECK(can_be_freely_mutated(device, 3) == false);
		CHECK(can_be_freely_mutated(device, 4) == false);
		CHECK(can_be_freely_mutated(device, 5) == false);
	}

	SECTION("multiple macrostates") {
		DevicePtr device = make_shared<Device>("UUUU");
		device->add_macrostate("a", "().)");
		device->add_macrostate("b", "(.))");

		CHECK(can_be_mutated(device, 0) == true);
		CHECK(can_be_mutated(device, 1) == true);
		CHECK(can_be_mutated(device, 2) == true);
		CHECK(can_be_mutated(device, 3) == true);

		CHECK(can_be_freely_mutated(device, 0) == true);
		CHECK(can_be_freely_mutated(device, 1) == false);
		CHECK(can_be_freely_mutated(device, 2) == false);
		CHECK(can_be_freely_mutated(device, 3) == false);
	}

}

TEST_CASE("Correctly mutate base-paired positions", "[sampling]") {
	DevicePtr device;

	struct Test {
		string sequence;
		vector<string> macrostates;
		string mutations;
		vector<string> expected_sequences;
	};

	// For each test, a device will be created with the sequence from the 
	// first field and the macrostates from the second field.  The third field 
	// specifies which mutations will be made to that device.  For example, 
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
			device = make_shared<Device>(test.sequence);
			CAPTURE(test.sequence)

			for(int x = 0; x < test.macrostates.size(); x++) {
				CAPTURE(test.macrostates[x]);
				device->add_macrostate(to_string(x), test.macrostates[x]);
			}

			CAPTURE(test.mutations);
			mutate_recursively(device, i, test.mutations[i]);
			CHECK(device->seq() == test.expected_sequences[i]);
		}
	}

	// Check some error conditions.
	device = make_shared<Device>("Nn");
	device->add_macrostate("not mutable", "()");
	CHECK_THROWS(mutate_recursively(device, 0, 'G'));

	device = make_shared<Device>("N");
	device->add_macrostate("extra open", "(");
	CHECK_THROWS(mutate_recursively(device, 0, 'G'));

	device = make_shared<Device>("N");
	device->add_macrostate("extra close", ")");
	CHECK_THROWS(mutate_recursively(device, 0, 'G'));


}
