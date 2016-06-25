#include <memory>
#include <catch/catch.hpp>
#include <sgrna_design/model.hh>

using namespace std;
using namespace sgrna_design;

TEST_CASE("Test domain constructor", "[model]") {
	Domain theo(
			"aptamer",
			"GAUACCAGCCGAAAGGCCCUUGGCAGC",
			ColorEnum::YELLOW,
			StyleEnum::BOLD);

	CHECK(theo.name() == "aptamer");
	CHECK(theo.seq() == "GAUACCAGCCGAAAGGCCCUUGGCAGC");
	CHECK(theo.color() == ColorEnum::YELLOW);
	CHECK(theo.style() == StyleEnum::BOLD);
}

TEST_CASE("Test Domain::mutate") {
	Domain dummy("dummy", "AAAA");

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

TEST_CASE("Test Domain::insert") {
	Domain dummy("dummy", "AAAA");

	SECTION("indices refer to positions between nucleotides") {
		dummy.seq("AAAA"); dummy.insert(0, "U"); CHECK(dummy.seq() == "UAAAA");
		dummy.seq("AAAA"); dummy.insert(1, "U"); CHECK(dummy.seq() == "AUAAA");
		dummy.seq("AAAA"); dummy.insert(2, "U"); CHECK(dummy.seq() == "AAUAA");
		dummy.seq("AAAA"); dummy.insert(3, "U"); CHECK(dummy.seq() == "AAAUA");
		dummy.seq("AAAA"); dummy.insert(4, "U"); CHECK(dummy.seq() == "AAAAU");
	}

	SECTION("negative indices count from the back") {
		dummy.seq("AAAA"); dummy.insert(-1, "U"); CHECK(dummy.seq() == "AAAAU");
		dummy.seq("AAAA"); dummy.insert(-2, "U"); CHECK(dummy.seq() == "AAAUA");
		dummy.seq("AAAA"); dummy.insert(-3, "U"); CHECK(dummy.seq() == "AAUAA");
		dummy.seq("AAAA"); dummy.insert(-4, "U"); CHECK(dummy.seq() == "AUAAA");
		dummy.seq("AAAA"); dummy.insert(-5, "U"); CHECK(dummy.seq() == "UAAAA");
	}

	SECTION("out-of-bounds indices throw exceptions") {
		CHECK_THROWS(dummy.insert(5, "U"));
		CHECK_THROWS(dummy.insert(-6, "U"));
	}
}

TEST_CASE("Test the construct class", "[model]") {
	Construct rna;

	rna += make_shared<Domain>("stem/5", "ACGU");
	rna += make_shared<Domain>("tetraloop", "GAAA");
	rna += make_shared<Domain>("stem/3", "ACGU");

	REQUIRE(rna.seq() == "ACGUGAAAACGU");

	REQUIRE(rna["stem/5"]->seq() == "ACGU");
	REQUIRE(rna["tetraloop"]->seq() == "GAAA");
	REQUIRE(rna["stem/3"]->seq() == "ACGU");

	REQUIRE(rna.index_5("stem/5") == 0);
	REQUIRE(rna.index_3("stem/5") == 3);
	REQUIRE(rna.index_5("tetraloop") == 4);
	REQUIRE(rna.index_3("tetraloop") == 7);
	REQUIRE(rna.index_5("stem/3") == 8);
	REQUIRE(rna.index_3("stem/3") == 11);

	ConstructPtr rna_copy = rna.copy();

	REQUIRE(rna_copy->seq() == "ACGUGAAAACGU");

	rna["stem/3"]->seq("UUUU");

	REQUIRE(rna.seq() == "ACGUGAAAUUUU");
	REQUIRE(rna_copy->seq() == "ACGUGAAAACGU");
}
