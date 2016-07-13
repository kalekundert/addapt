#include <memory>
#include <catch/catch.hpp>
#include "model.hh"

using namespace std;
using namespace addapt;

TEST_CASE("Test the Construct class", "[model]") {
	Construct dummy("ACGU");
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

	ConstructPtr dummy_2 = dummy.copy();
	CHECK(dummy_2->seq() == "ACGU");
	CHECK(dummy_2->macrostate("a") == "....");
	CHECK(dummy_2->macrostate("b") == "(())");

	Construct dummy_3("nnnn");
	CHECK(dummy_3.seq() == "nnnn");
	dummy_3.assign(dummy_2);
	CHECK(dummy_3.seq() == "ACGU");
	CHECK(dummy_3.macrostate("a") == "....");
	CHECK(dummy_3.macrostate("b") == "(())");
}

TEST_CASE("Test the Construct::mutate method", "[model]") {
	Construct dummy("AAAA");

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


