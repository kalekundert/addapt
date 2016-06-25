#include <memory>
#include <catch/catch.hpp>
#include <sgrna_design/model.hh>

using namespace std;
using namespace sgrna_design;

TEST_CASE("Test domain constructor", "[model]") {
	Domain domain(
			"aptamer",
			"GAUACCAGCCGAAAGGCCCUUGGCAGC",
			ColorEnum::YELLOW,
			StyleEnum::BOLD);

	CHECK(domain.name() == "aptamer");
	CHECK(domain.seq() == "GAUACCAGCCGAAAGGCCCUUGGCAGC");
	CHECK(domain.color() == ColorEnum::YELLOW);
	CHECK(domain.style() == StyleEnum::BOLD);
}

TEST_CASE("Test Domain::mutate") {
	Domain domain("dummy", "AAAA");

	SECTION("indices refer to the right positions") {
		domain.mutate(0, 'U');
		CHECK(domain.seq() == "UAAA");
		domain.mutate(1, 'U');
		CHECK(domain.seq() == "UUAA");
		domain.mutate(2, 'U');
		CHECK(domain.seq() == "UUUA");
		domain.mutate(3, 'U');
		CHECK(domain.seq() == "UUUU");
	}

	SECTION("negative indices count from the back") {
		domain.mutate(-1, 'U');
		CHECK(domain.seq() == "AAAU");
		domain.mutate(-2, 'U');
		CHECK(domain.seq() == "AAUU");
		domain.mutate(-3, 'U');
		CHECK(domain.seq() == "AUUU");
		domain.mutate(-4, 'U');
		CHECK(domain.seq() == "UUUU");
	}

	SECTION("out-of-bounds indices throw exceptions") {
		CHECK_THROWS(domain.mutate(4, 'U'));
		CHECK_THROWS(domain.mutate(-5, 'U'));
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
