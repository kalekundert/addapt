#include <memory>
#include <catch/catch.hpp>
#include "model.hh"

using namespace std;
using namespace addapt;

TEST_CASE("Test the Device class", "[model]") {
	Device dummy("ACGU");
	CHECK(dummy.seq() == "ACGU");
	CHECK(dummy.len() == 4);

	CHECK(dummy.seq(0) == 'A');
	CHECK(dummy.seq(1) == 'C');
	CHECK(dummy.seq(2) == 'G');
	CHECK(dummy.seq(3) == 'U');
	CHECK_THROWS(dummy.seq(4));

	CHECK(dummy.seq(-1) == 'U');
	CHECK(dummy.seq(-2) == 'G');
	CHECK(dummy.seq(-3) == 'C');
	CHECK(dummy.seq(-4) == 'A');
	CHECK_THROWS(dummy.seq(-5));

	dummy.add_macrostate("a", "....");
	dummy.add_macrostate("b", "(())");
	CHECK(dummy.macrostate("a") == "....");
	CHECK(dummy.macrostate("b") == "(())");

	DevicePtr dummy_2 = dummy.copy();
	CHECK(dummy_2->seq() == "ACGU");
	CHECK(dummy_2->macrostate("a") == "....");
	CHECK(dummy_2->macrostate("b") == "(())");

	Device dummy_3("nnnn");
	CHECK(dummy_3.seq() == "nnnn");
	dummy_3.assign(dummy_2);
	CHECK(dummy_3.seq() == "ACGU");
	CHECK(dummy_3.macrostate("a") == "....");
	CHECK(dummy_3.macrostate("b") == "(())");
}

TEST_CASE("Test the Device class with contexts", "[model]") {
	Device dummy("C");
	dummy.add_macrostate("bp", "x");

	CHECK(dummy.context()->before() == "");
	CHECK(dummy.context()->after() == "");

	// Add a context and make sure it affects all the expected members.
	dummy.context(make_shared<Context>("A", "GU"));

	CHECK(dummy.len() == 4);
	CHECK(dummy.seq() == "ACGU");
	CHECK(dummy.seq(0) == 'A');
	CHECK(dummy.seq(1) == 'C');
	CHECK(dummy.seq(2) == 'G');
	CHECK(dummy.seq(3) == 'U');
	CHECK(dummy.raw_len() == 1);
	CHECK(dummy.raw_seq() == "C");
	CHECK(dummy.raw_seq(0) == 'C');
	CHECK(dummy.macrostate("bp") == ".x..");
	CHECK(dummy.context()->before() == "A");
	CHECK(dummy.context()->after() == "GU");

	for(auto pair: dummy.macrostates()) {
		CHECK(pair.first == "bp");
		CHECK(pair.second == ".x..");
	}

	// Remove the context and make sure all the expected members revert to their 
	// previous values.
	dummy.remove_context();

	CHECK(dummy.len() == 1);
	CHECK(dummy.seq() == "C");
	CHECK(dummy.seq(0) == 'C');
	CHECK(dummy.raw_len() == 1);
	CHECK(dummy.raw_seq() == "C");
	CHECK(dummy.raw_seq(0) == 'C');
	CHECK(dummy.macrostate("bp") == "x");
	CHECK(dummy.context()->before() == "");
	CHECK(dummy.context()->after() == "");

	for(auto pair: dummy.macrostates()) {
		CHECK(pair.first == "bp");
		CHECK(pair.second == "x");
	}
}

TEST_CASE("Test the Device::mutate method", "[model]") {
	Device dummy("AAAA");

	SECTION("positive indices count from the front") {
		dummy.mutate(0, 'U'); CHECK(dummy.seq() == "UAAA");
		dummy.mutate(1, 'U'); CHECK(dummy.seq() == "UUAA");
		dummy.mutate(2, 'U'); CHECK(dummy.seq() == "UUUA");
		dummy.mutate(3, 'U'); CHECK(dummy.seq() == "UUUU");
	}

	SECTION("negative indices count from the back") {
		dummy.mutate(-1, 'U'); CHECK(dummy.seq() == "AAAU");
		dummy.mutate(-2, 'U'); CHECK(dummy.seq() == "AAUU");
		dummy.mutate(-3, 'U'); CHECK(dummy.seq() == "AUUU");
		dummy.mutate(-4, 'U'); CHECK(dummy.seq() == "UUUU");
	}

	SECTION("out-of-bounds indices throw exceptions") {
		CHECK_THROWS(dummy.mutate(4, 'U'));
		CHECK_THROWS(dummy.mutate(-5, 'U'));
	}
}

TEST_CASE("Test the Aptamer class", "[model]") {
	Aptamer theo(
			"GAUACCAGCCGAAAGGCCCUUGGCAGC",
			"(...((.(((....)))....))...)",
			0.320);

	CHECK(theo.seq() == "GAUACCAGCCGAAAGGCCCUUGGCAGC");
	CHECK(theo.fold() == "(...((.(((....)))....))...)");
	CHECK(theo.affinity() == 0.320);
}


