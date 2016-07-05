#include <memory>
#include <catch/catch.hpp>
#include <sgrna_design/model.hh>

using namespace std;
using namespace sgrna_design;

TEST_CASE("Test the Domain constructor", "[model]") {
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

TEST_CASE("Test the Domain::mutate method", "[model]") {
	Domain dummy("dummy");

	SECTION("positive indices count from the front") {
		dummy.seq("AAAA"); dummy.mutate(0, 'U'); CHECK(dummy.seq() == "UAAA");
		dummy.seq("AAAA"); dummy.mutate(1, 'U'); CHECK(dummy.seq() == "AUAA");
		dummy.seq("AAAA"); dummy.mutate(2, 'U'); CHECK(dummy.seq() == "AAUA");
		dummy.seq("AAAA"); dummy.mutate(3, 'U'); CHECK(dummy.seq() == "AAAU");
	}

	SECTION("negative indices count from the back") {
		dummy.seq("AAAA"); dummy.mutate(-1, 'U'); CHECK(dummy.seq() == "AAAU");
		dummy.seq("AAAA"); dummy.mutate(-2, 'U'); CHECK(dummy.seq() == "AAUA");
		dummy.seq("AAAA"); dummy.mutate(-3, 'U'); CHECK(dummy.seq() == "AUAA");
		dummy.seq("AAAA"); dummy.mutate(-4, 'U'); CHECK(dummy.seq() == "UAAA");
	}

	SECTION("out-of-bounds indices throw exceptions") {
		CHECK_THROWS(dummy.mutate(4, 'U'));
		CHECK_THROWS(dummy.mutate(-5, 'U'));
	}
}

TEST_CASE("Test the Domain::insert method", "[model]") {
	Domain dummy("dummy");

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

	SECTION("multiple nucleotides can be inserted") {
		dummy.seq("AAAA"); dummy.insert(2, "");      CHECK(dummy.seq() == "AAAA");
		dummy.seq("AAAA"); dummy.insert(2, "U");     CHECK(dummy.seq() == "AAUAA");
		dummy.seq("AAAA"); dummy.insert(2, "UU");    CHECK(dummy.seq() == "AAUUAA");
		dummy.seq("AAAA"); dummy.insert(2, "UUU");   CHECK(dummy.seq() == "AAUUUAA");
		dummy.seq("AAAA"); dummy.insert(2, "UUUU");  CHECK(dummy.seq() == "AAUUUUAA");
		dummy.seq("AAAA"); dummy.insert(2, "UUUUU"); CHECK(dummy.seq() == "AAUUUUUAA");
	}

	SECTION("chars can also be inserted") {
		dummy.seq("AAAA"); dummy.insert(2, 'U'); CHECK(dummy.seq() == "AAUAA");
	}

	SECTION("out-of-bounds indices throw exceptions") {
		CHECK_THROWS(dummy.insert(5, "U"));
		CHECK_THROWS(dummy.insert(-6, "U"));
	}
}

TEST_CASE("Test the Domain::remove(int) method", "[model]") {
	Domain dummy("dummy");

	SECTION("the index refers to the nucleotides, counting from 0") {
		dummy.seq("ACGU"); dummy.remove(0); CHECK(dummy.seq() == "CGU");
		dummy.seq("ACGU"); dummy.remove(1); CHECK(dummy.seq() == "AGU");
		dummy.seq("ACGU"); dummy.remove(2); CHECK(dummy.seq() == "ACU");
		dummy.seq("ACGU"); dummy.remove(3); CHECK(dummy.seq() == "ACG");
	}

	SECTION("negative indices count from the back") {
		dummy.seq("ACGU"); dummy.remove(-1); CHECK(dummy.seq() == "ACG");
		dummy.seq("ACGU"); dummy.remove(-2); CHECK(dummy.seq() == "ACU");
		dummy.seq("ACGU"); dummy.remove(-3); CHECK(dummy.seq() == "AGU");
		dummy.seq("ACGU"); dummy.remove(-4); CHECK(dummy.seq() == "CGU");
	}

	SECTION("out-of-bounds indices throw exceptions") {
		CHECK_THROWS(dummy.remove(3));
		CHECK_THROWS(dummy.remove(-4));
	}
}

TEST_CASE("Test the Domain::remove(int, int) method", "[model]") {
	Domain dummy("dummy");

	SECTION("two indices: refer to positions between nucleotides") {
		dummy.seq("ACGU"); dummy.remove(0, 0); CHECK(dummy.seq() == "ACGU");
		dummy.seq("ACGU"); dummy.remove(0, 1); CHECK(dummy.seq() == "CGU");
		dummy.seq("ACGU"); dummy.remove(0, 2); CHECK(dummy.seq() == "GU");
		dummy.seq("ACGU"); dummy.remove(0, 3); CHECK(dummy.seq() == "U");
		dummy.seq("ACGU"); dummy.remove(0, 4); CHECK(dummy.seq() == "");
	}

	SECTION("negative indices count from the back") {
		dummy.seq("ACGU"); dummy.remove(-1, -1); CHECK(dummy.seq() == "ACGU");
		dummy.seq("ACGU"); dummy.remove(-1, -2); CHECK(dummy.seq() == "ACG");
		dummy.seq("ACGU"); dummy.remove(-1, -3); CHECK(dummy.seq() == "AC");
		dummy.seq("ACGU"); dummy.remove(-1, -4); CHECK(dummy.seq() == "A");
		dummy.seq("ACGU"); dummy.remove(-1, -5); CHECK(dummy.seq() == "");
	}

	SECTION("positive and negative indices can be mixed") {
		dummy.seq("ACGU"); dummy.remove(0, -1); CHECK(dummy.seq() == "");
		dummy.seq("ACGU"); dummy.remove(0, -2); CHECK(dummy.seq() == "U");
		dummy.seq("ACGU"); dummy.remove(0, -3); CHECK(dummy.seq() == "GU");
		dummy.seq("ACGU"); dummy.remove(0, -4); CHECK(dummy.seq() == "CGU");
		dummy.seq("ACGU"); dummy.remove(0, -5); CHECK(dummy.seq() == "ACGU");

		dummy.seq("ACGU"); dummy.remove(-1, 0); CHECK(dummy.seq() == "");
		dummy.seq("ACGU"); dummy.remove(-1, 1); CHECK(dummy.seq() == "A");
		dummy.seq("ACGU"); dummy.remove(-1, 2); CHECK(dummy.seq() == "AC");
		dummy.seq("ACGU"); dummy.remove(-1, 3); CHECK(dummy.seq() == "ACG");
		dummy.seq("ACGU"); dummy.remove(-1, 4); CHECK(dummy.seq() == "ACGU");
	}

	SECTION("the order of the indices doesn't matter") {
		dummy.seq("ACGU"); dummy.remove( 1,  3); CHECK(dummy.seq() == "AU");
		dummy.seq("ACGU"); dummy.remove( 3,  1); CHECK(dummy.seq() == "AU");
		dummy.seq("ACGU"); dummy.remove(-2, -4); CHECK(dummy.seq() == "AU");
		dummy.seq("ACGU"); dummy.remove(-4, -2); CHECK(dummy.seq() == "AU");
		dummy.seq("ACGU"); dummy.remove( 1, -2); CHECK(dummy.seq() == "AU");
		dummy.seq("ACGU"); dummy.remove(-2,  1); CHECK(dummy.seq() == "AU");
		dummy.seq("ACGU"); dummy.remove( 3, -4); CHECK(dummy.seq() == "AU");
		dummy.seq("ACGU"); dummy.remove(-4,  3); CHECK(dummy.seq() == "AU");
	}

	SECTION("out-of-bounds indices throw exceptions") {
		CHECK_THROWS(dummy.remove(0, 5));
		CHECK_THROWS(dummy.remove(0, -6));
	}
}

TEST_CASE("Test the Domain::replace method", "[model]") {
	Domain dummy("dummy");

	SECTION("indices refer to positions between nucleotides") {
		dummy.seq("AAAA"); dummy.replace(0, 0, "U"); CHECK(dummy.seq() == "UAAAA");
		dummy.seq("AAAA"); dummy.replace(0, 1, "U"); CHECK(dummy.seq() == "UAAA");
		dummy.seq("AAAA"); dummy.replace(0, 2, "U"); CHECK(dummy.seq() == "UAA");
		dummy.seq("AAAA"); dummy.replace(0, 3, "U"); CHECK(dummy.seq() == "UA");
		dummy.seq("AAAA"); dummy.replace(0, 4, "U"); CHECK(dummy.seq() == "U");
	}

	SECTION("negative indices count from the back") {
		dummy.seq("AAAA"); dummy.replace(-1, -1, "U"); CHECK(dummy.seq() == "AAAAU");
		dummy.seq("AAAA"); dummy.replace(-1, -2, "U"); CHECK(dummy.seq() == "AAAU");
		dummy.seq("AAAA"); dummy.replace(-1, -3, "U"); CHECK(dummy.seq() == "AAU");
		dummy.seq("AAAA"); dummy.replace(-1, -4, "U"); CHECK(dummy.seq() == "AU");
		dummy.seq("AAAA"); dummy.replace(-1, -5, "U"); CHECK(dummy.seq() == "U");
	}

	SECTION("positive and negative indices can be mixed") {
		dummy.seq("AAAA"); dummy.replace(0, -1, "U"); CHECK(dummy.seq() == "U");
		dummy.seq("AAAA"); dummy.replace(0, -2, "U"); CHECK(dummy.seq() == "UA");
		dummy.seq("AAAA"); dummy.replace(0, -3, "U"); CHECK(dummy.seq() == "UAA");
		dummy.seq("AAAA"); dummy.replace(0, -4, "U"); CHECK(dummy.seq() == "UAAA");
		dummy.seq("AAAA"); dummy.replace(0, -5, "U"); CHECK(dummy.seq() == "UAAAA");

		dummy.seq("AAAA"); dummy.replace(-1, 0, "U"); CHECK(dummy.seq() == "U");
		dummy.seq("AAAA"); dummy.replace(-1, 1, "U"); CHECK(dummy.seq() == "AU");
		dummy.seq("AAAA"); dummy.replace(-1, 2, "U"); CHECK(dummy.seq() == "AAU");
		dummy.seq("AAAA"); dummy.replace(-1, 3, "U"); CHECK(dummy.seq() == "AAAU");
		dummy.seq("AAAA"); dummy.replace(-1, 4, "U"); CHECK(dummy.seq() == "AAAAU");
	}

	SECTION("the order of the indices doesn't matter") {
		dummy.seq("AAAA"); dummy.replace( 1,  3, "U"); CHECK(dummy.seq() == "AUA");
		dummy.seq("AAAA"); dummy.replace( 3,  1, "U"); CHECK(dummy.seq() == "AUA");
		dummy.seq("AAAA"); dummy.replace(-2, -4, "U"); CHECK(dummy.seq() == "AUA");
		dummy.seq("AAAA"); dummy.replace(-4, -2, "U"); CHECK(dummy.seq() == "AUA");
		dummy.seq("AAAA"); dummy.replace( 1, -2, "U"); CHECK(dummy.seq() == "AUA");
		dummy.seq("AAAA"); dummy.replace(-2,  1, "U"); CHECK(dummy.seq() == "AUA");
		dummy.seq("AAAA"); dummy.replace( 3, -4, "U"); CHECK(dummy.seq() == "AUA");
		dummy.seq("AAAA"); dummy.replace(-4,  3, "U"); CHECK(dummy.seq() == "AUA");
	}

	SECTION("multiple nucleotides can be replaceed") {
		dummy.seq("AAAA"); dummy.replace(1, 3, "");      CHECK(dummy.seq() == "AA");
		dummy.seq("AAAA"); dummy.replace(1, 3, "U");     CHECK(dummy.seq() == "AUA");
		dummy.seq("AAAA"); dummy.replace(1, 3, "UU");    CHECK(dummy.seq() == "AUUA");
		dummy.seq("AAAA"); dummy.replace(1, 3, "UUU");   CHECK(dummy.seq() == "AUUUA");
		dummy.seq("AAAA"); dummy.replace(1, 3, "UUUU");  CHECK(dummy.seq() == "AUUUUA");
		dummy.seq("AAAA"); dummy.replace(1, 3, "UUUUU"); CHECK(dummy.seq() == "AUUUUUA");
	}

	SECTION("out-of-bounds indices throw exceptions") {
		CHECK_THROWS(dummy.replace(0, 5, "U"));
		CHECK_THROWS(dummy.replace(0, -6, "U"));
	}
}

TEST_CASE("Test the Construct class", "[model]") {
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
