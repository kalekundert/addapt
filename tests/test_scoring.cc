#include <set>
#include <vector>
#include <catch/catch.hpp>
#include <sgrna_design/model.hh>
#include <sgrna_design/scoring.hh>
#include <sgrna_design/utils.hh>

using namespace std;
using namespace sgrna_design;
using bp = pair<int,int>;

class DummyRnaFold : public RnaFold {

public:

	double &
	operator[](bp key) {
		key = {
			std::min(key.first, key.second),
			std::max(key.first, key.second)};
		return my_base_pair_probs[key];
	}

	double
	base_pair_prob(int a, int b) const {
		bp key = {std::min(a,b), std::max(a,b)};
		auto it = my_base_pair_probs.find(key);
		return (it != my_base_pair_probs.end())? it->second : 0.0;
	}

private:
	map<bp,double> my_base_pair_probs;

};

ConstructConstPtr
build_rhf_6_construct() {
	ConstructPtr rhf_6 = make_shared<Construct>();

	// 0....,....1....,....2....,....3....,....4....,....5....,....6....,....7....,....8....,....9....,....0.
	// GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCCUUUUCGCCGAUACCAGCCGAAAGGCCCUUGGCAGCGACGGCACCGAGUCGGUGCUUUUUU
	// (((((((.((((....))))...))))))).,,({{,{..|||{{(,((,{....,.||{}}}})),..,}))).,,||.(((((((...)))))))..... [ΔG=-29.58 kcal/mol, apo]
	// (((((((.((((....))))...))))))){(((......)))}..{{.((...((.(((....)))....))...)).}|((((((...)))))),..... [ΔG=-33.82 kcal/mol, holo]

	*rhf_6 += make_shared<Domain>("lower_stem/a", "guuuua");
	*rhf_6 += make_shared<Domain>("upper_stem", "gagcuagaaauagcaagu");
	*rhf_6 += make_shared<Domain>("lower_stem/b", "uaaaau");
	*rhf_6 += make_shared<Domain>("nexus/x", "aagg");
	*rhf_6 += make_shared<Domain>("nexus/y", "cuagu");
	*rhf_6 += make_shared<Domain>("nexus/z", "ccCuU");
	*rhf_6 += make_shared<Domain>("ruler", "UUC");
	*rhf_6 += make_shared<Domain>("hairpin/u", "GCC");
	*rhf_6 += make_shared<Domain>("aptamer", "gauaccagccgaaaggcccuuggcagc");
	*rhf_6 += make_shared<Domain>("hairpin/v", "GAC");
	*rhf_6 += make_shared<Domain>("tail", "ggcaccgagucggugcuuuuuu");

	return rhf_6;
}


TEST_CASE("Test the DummyRnaFold helper class") {
	DummyRnaFold fold;

	// Base pairs have 0 probability unless otherwise specified.
	CHECK(fold.base_pair_prob(1,2) == 0.0);

	// Base pair probabilities can be set.  The order of the base pair indices 
	// doesn't matter.
	fold[{1,2}] = 0.75;
	CHECK(fold.base_pair_prob(1,2) == 0.75);
	CHECK(fold.base_pair_prob(2,1) == 0.75);

	fold[{2,1}] = 0.50;
	CHECK(fold.base_pair_prob(1,2) == 0.50);
	CHECK(fold.base_pair_prob(2,1) == 0.50);
}

TEST_CASE("Test folding a hairpin without an aptamer", "[scoring]") {
	// Make a construct that should fold into a hairpin:
	
	// ACGUGAAAACGU
	// ((((....)))) [ΔG=-2.20 kcal/mol]

	ConstructPtr hairpin = make_shared<Construct>();
	*hairpin += make_shared<Domain>("stem/a", "ACGU");
	*hairpin += make_shared<Domain>("loop", "GAAA");
	*hairpin += make_shared<Domain>("stem/b", "ACGU");

	// Indicate which base pairs I expect to find.
	set<bp> expected_base_pairs = {{0,11}, {1,10}, {2,9}, {3,8}};
	map<bp,double> base_pair_probs = {
		{{0,11}, 0.70},
		{{1,10}, 0.95},
		{{2, 9}, 0.95},
		{{3, 8}, 0.95}
	};

	ViennaRnaFold fold(hairpin);

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

TEST_CASE("Test folding a hairpin with an aptamer", "[scoring]") {
	// Make a construct consisting entirely of an aptamer:

	// GAUACCAGCCGAAAGGCCCUUGGCAGC
	// ....((((((....)))...))).... [ΔG=-6.20 kcal/mol, apo]
	// (...((.(((....)))....))...) [ΔG=-9.22 kcal/mol, holo]

	ConstructPtr hairpin = make_shared<Construct>();
	*hairpin += make_shared<Domain>("aptamer", "GAUACCAGCCGAAAGGCCCUUGGCAGC");

	// Indicate which base pairs I expect to form.
	set<bp> apo_base_pairs = {
		{4,22}, {5,21}, {6,20}, {7,16}, {8,15}, {9,14}
	};
	set<bp> holo_base_pairs = {
		{0,26}, {4,22}, {5,21}, {7,16}, {8,15}, {9,14}
	};

	// Make sure the expected base pairs form.
	ViennaRnaFold apo_fold(hairpin, LigandEnum::NONE);
	ViennaRnaFold holo_fold(hairpin, LigandEnum::THEO);

	for(int i = 0; i < hairpin->len(); i++) {
		for(int j = i; j < hairpin->len(); j++) {
			double apo_bp = apo_fold.base_pair_prob(i,j);
			double holo_bp = holo_fold.base_pair_prob(i,j);

			if(apo_base_pairs.count({i,j})) {
				CHECK(apo_bp > 0.7);
			}
			else {
				REQUIRE(apo_bp < 0.3);
			}
			if(holo_base_pairs.count({i,j})) {
				CHECK(holo_bp > 0.7);
			}
			else {
				REQUIRE(holo_bp < 0.3);
			}
		}
	}
}

TEST_CASE("Test folding rhf(6)") {
	ConstructConstPtr rhf_6 = build_rhf_6_construct();

	// Indicate which base pairs I expect to form.
	map<bp,double> constitutive_base_pairs = {
		{{ 0,29}, 0.75}, {{ 1,28}, 0.80}, {{ 2,27}, 0.85}, {{ 3,26}, 0.85},
		{{ 4,25}, 0.85}, {{ 5,24}, 0.85}, {{ 6,23}, 0.80}, {{ 8,19}, 0.95},
		{{ 9,18}, 0.95}, {{10,17}, 0.95}, {{11,16}, 0.95}, {{81,95}, 0.90},
		{{82,94}, 0.90}, {{83,93}, 0.90}, {{84,92}, 0.90}, {{85,91}, 0.90},
		{{86,90}, 0.75}
	};
	map<bp,double> apo_base_pairs = {
		{{33,73}, 0.30}, {{34,72}, 0.30}, {{35,71}, 0.30}, {{36,70}, 0.25},
		{{37,69}, 0.10}, {{40,65}, 0.35}, {{41,64}, 0.35}, {{42,63}, 0.35},
		{{43,62}, 0.35}, {{44,61}, 0.30}, {{45,60}, 0.25}, {{47,58}, 0.35},
		{{48,57}, 0.40}
	};
	map<bp,pair<double,double> > holo_base_pairs = {
		{{30,43}, {0.55,0.25}}, {{31,42}, {0.75,0.35}}, {{32,41}, {0.75,0.35}},
		{{33,40}, {0.75,0.35}}, {{46,80}, {0.50,0.05}}, {{47,79}, {0.60,0.05}},
		{{49,77}, {0.85,0.05}}, {{50,76}, {0.95,0.05}}, {{54,72}, {0.95,0.20}},
		{{55,71}, {0.95,0.20}}, {{57,66}, {0.95,0.50}}, {{58,65}, {0.95,0.50}},
		{{59,64}, {0.95,0.50}}
	};

	// Make sure the expected base pairs form in the expected conditions.
	ViennaRnaFold apo_fold(rhf_6, LigandEnum::NONE);
	ViennaRnaFold holo_fold(rhf_6, LigandEnum::THEO);

	for(auto base_pair: constitutive_base_pairs) {
		int i = base_pair.first.first;
		int j = base_pair.first.second;
		double threshold = base_pair.second;

		REQUIRE(apo_fold.base_pair_prob(i,j) > threshold);
		REQUIRE(holo_fold.base_pair_prob(i,j) > threshold);
	}

	for(auto base_pair: apo_base_pairs) {
		int i = base_pair.first.first;
		int j = base_pair.first.second;
		double threshold = base_pair.second;

		REQUIRE(apo_fold.base_pair_prob(i,j) > threshold);
		REQUIRE(holo_fold.base_pair_prob(i,j) < 1e-3);
	}

	for(auto base_pair: holo_base_pairs) {
		int i = base_pair.first.first;
		int j = base_pair.first.second;
		double apo_threshold = base_pair.second.second;
		double holo_threshold = base_pair.second.first;

		CHECK(apo_fold.base_pair_prob(i,j) < apo_threshold);
		CHECK(holo_fold.base_pair_prob(i,j) > holo_threshold);
	}
}

TEST_CASE("Test the 'ligand sensitivity' score term", "[scoring]") {
	// Make a construct with individually indexable nucleotides.
	ConstructPtr dummy_construct = make_shared<Construct>();
	*dummy_construct += make_shared<Domain>("a", "UU");
	*dummy_construct += make_shared<Domain>("b", "U");
	*dummy_construct += make_shared<Domain>("c", "U");

	// Define fake base-pairing probabilities for this construct.
	DummyRnaFold dummy_apo_fold;
	DummyRnaFold dummy_holo_fold;

	dummy_apo_fold[{0,3}] = 0.10;
	dummy_holo_fold[{0,3}] = 0.40;

	dummy_apo_fold[{1,2}] = 0.30;
	dummy_holo_fold[{1,2}] = 0.10;

	// Define expected scores for various combinations of nucleotides.
	struct Test {
		vector<string> selection;
		double expected_score;
	};

	vector<Test> tests = {
		// The score is 0 if there aren't any base pairs.
		{{"a"}, 0},
		{{"b"}, 0},
		{{"c"}, 0},

		// The score increases with the number of ligand-sensitive base pairs.
		{{"a","b"}, (0.3-0.1)*log(3)/3},
		{{"a","c"}, (0.4-0.1)*log(4)/3},
		{{"a","b","c"}, (0.3-0.1)*log(3)/4 + (0.4-0.1)*log(4)/4}
	};

	// Make sure each scenario is scored correctly.

	for(Test test: tests) {
		LigandSensitivityTerm term(test.selection);
		double score = term.evaluate(
				dummy_construct, dummy_apo_fold, dummy_holo_fold);

		CAPTURE(test.selection);
		CHECK(score == Approx(test.expected_score));
	}

	// You can't score empty selections.

	CHECK_THROWS(LigandSensitivityTerm term({}));

}

TEST_CASE("Test the 'specific ligand sensitivity' score term", "[scoring]") {
	// Make a construct with individually indexable nucleotides.
	ConstructPtr dummy_construct = make_shared<Construct>();
	*dummy_construct += make_shared<Domain>("a", "U");
	*dummy_construct += make_shared<Domain>("b", "U");
	*dummy_construct += make_shared<Domain>("B", "U");
	*dummy_construct += make_shared<Domain>("A", "U");

	// Define fake base-pairing probabilities for this construct.
	DummyRnaFold dummy_apo_fold;
	DummyRnaFold dummy_holo_fold;

	dummy_apo_fold[{0,3}] = 0.10;
	dummy_holo_fold[{0,3}] = 0.40;

	dummy_apo_fold[{1,2}] = 0.30;
	dummy_holo_fold[{1,2}] = 0.10;

	// Define expected scores for various combinations of nucleotides.
	struct Test {
		ConditionEnum condition;
		vector<string> selection;
		vector<string> targets;
		double expected_score;
	};

	vector<Test> tests = {
		// The score is 0 if there aren't any base pairs.
		{ConditionEnum::APO,  {"a"}, {}, 0},
		{ConditionEnum::HOLO, {"a"}, {}, 0},
		{ConditionEnum::APO,  {"a"}, {"a"}, 0},
		{ConditionEnum::HOLO, {"a"}, {"a"}, 0},
		{ConditionEnum::APO,  {"a"}, {"B"}, 0},
		{ConditionEnum::HOLO, {"a"}, {"B"}, 0},
		{ConditionEnum::APO,  {"b"}, {}, 0},
		{ConditionEnum::HOLO, {"b"}, {}, 0},
		{ConditionEnum::APO,  {"b"}, {"b"}, 0},
		{ConditionEnum::HOLO, {"b"}, {"b"}, 0},
		{ConditionEnum::APO,  {"b"}, {"A"}, 0},
		{ConditionEnum::HOLO, {"b"}, {"A"}, 0},

		// The score is proportional to the ligand sensitivity of the base pair.
		{ConditionEnum::APO,  {"a"}, {"A"}, -0.1*log(4)},
		{ConditionEnum::HOLO, {"a"}, {"A"},  0.4*log(4)},
		{ConditionEnum::APO,  {"b"}, {"B"},  0.3*log(3)},
		{ConditionEnum::HOLO, {"b"}, {"B"}, -0.1*log(3)},

		// The score is normalized by the number of nucleotides in the selection.
		{ConditionEnum::APO,  {"a"},     {"A","B"}, -0.1*log(4)/1},
		{ConditionEnum::HOLO, {"a"},     {"A","B"},  0.4*log(4)/1},
		{ConditionEnum::APO,  {"b"},     {"A","B"},  0.3*log(3)/1},
		{ConditionEnum::HOLO, {"b"},     {"A","B"}, -0.1*log(3)/1},
		{ConditionEnum::APO,  {"a","b"}, {"A"},     -0.1*log(4)/2},
		{ConditionEnum::HOLO, {"a","b"}, {"A"},      0.4*log(4)/2},
		{ConditionEnum::APO,  {"a","b"}, {"B"},      0.3*log(3)/2},
		{ConditionEnum::HOLO, {"a","b"}, {"B"},     -0.1*log(3)/2},
		{ConditionEnum::APO,  {"a","b"}, {"A","B"}, -0.1*log(4)/2 +0.3*log(3)/2},
		{ConditionEnum::HOLO, {"a","b"}, {"A","B"},  0.4*log(4)/2 -0.1*log(3)/2},
	};

	// Make sure each scenario is scored correctly.

	for(Test test: tests) {
		SpecificLigandSensitivityTerm term(
				test.condition, test.selection, test.targets);
		double score = term.evaluate(
				dummy_construct, dummy_apo_fold, dummy_holo_fold);

		CAPTURE(test.condition);
		CAPTURE(test.selection);
		CAPTURE(test.targets);

		CHECK(score == Approx(test.expected_score));
	}

	// You can't score empty selections.

	CHECK_THROWS(
			SpecificLigandSensitivityTerm term(ConditionEnum::APO, {}, {"a"}));

}


/*

TEST_CASE("Test the BasePairingTerm::evaluate method", "[scoring]") {
	// Make a construct that will fold very differently with and without ligand.  
	// This construct has an aptamer on the 5' end and a decoy sequence (that 
	// pairs internally with the aptamer) on the 3' end.  The decoy stem should 
	// form in the absence of ligand, while the ends of the aptamer should base 
	// pair with ligand:

	// G AUACCA GCC GAAA GGC CCUUGGCA GCC UUUC
	// . ...... ... (((( ((( ........ ))) ))))  [ΔG=-8.80 kcal/mol, apo]
	// ( ...((. ((( .... ))) ....)).. .). ....  [ΔG=-11.81 kcal/mol, holo]

	ConstructPtr sensor = make_shared<Construct>();
	*sensor += make_shared<Domain>("target", "G");
	*sensor += make_shared<Domain>("aptamer/a", "AUACCA");
	*sensor += make_shared<Domain>("aptamer/b", "GCC");
	*sensor += make_shared<Domain>("aptamer/c", "GAAA");
	*sensor += make_shared<Domain>("aptamer/d", "GGC");
	*sensor += make_shared<Domain>("aptamer/e", "CCUUGGCA");
	*sensor += make_shared<Domain>("decoy/a", "GCC");
	*sensor += make_shared<Domain>("decoy/b", "UUUC");

	// Create the fold objects the score terms will need.

	RnaFold no_lig_fold(sensor, LigandEnum::NONE);
	RnaFold lig_fold(sensor, LigandEnum::THEO);

	// Evaluate the base-pairing of the ligand-unbound state.

	SECTION("evaluate the decoy stem") {
		BasePairingTerm c_term(
				ConditionEnum::APO,
				{"aptamer/c"},
				{"decoy/b"}
		);
		BasePairingTerm d_term(
				ConditionEnum::APO,
				{"aptamer/d"},
				{"decoy/a"}
		);
		BasePairingTerm cd_term(
				ConditionEnum::APO,
				{"aptamer/c", "aptamer/d"},
				{"decoy/a", "decoy/b"}
		);

		double c_eval = c_term.evaluate(sensor, no_lig_fold, lig_fold);
		double d_eval = d_term.evaluate(sensor, no_lig_fold, lig_fold);
		double cd_eval = cd_term.evaluate(sensor, no_lig_fold, lig_fold);

		CHECK(c_eval > 4 * 0.9);
		CHECK(d_eval > 3 * 0.9);
		CHECK(cd_eval == Approx(c_eval + d_eval));
	}

	// Evaluate the base-pairing of the ligand-bound state.

	SECTION("evaluate the target stem") {
		BasePairingTerm t_term(
				ConditionEnum::HOLO,
				{"target"},
				{"decoy/a"}
		);
		BasePairingTerm a_term(
				ConditionEnum::HOLO,
				{"aptamer/a"},
				{"aptamer/e"}
		);
		BasePairingTerm b_term(
				ConditionEnum::HOLO,
				{"aptamer/b"},
				{"aptamer/d"}
		);
		BasePairingTerm tab_term(
				ConditionEnum::HOLO,
				{"target", "aptamer/a", "aptamer/b"},
				{"aptamer/d", "aptamer/e", "decoy/a"}
		);

		double t_eval = t_term.evaluate(sensor, no_lig_fold, lig_fold);
		double a_eval = a_term.evaluate(sensor, no_lig_fold, lig_fold);
		double b_eval = b_term.evaluate(sensor, no_lig_fold, lig_fold);
		double tab_eval = tab_term.evaluate(sensor, no_lig_fold, lig_fold);

		CHECK(t_eval > 1 * 0.9);
		CHECK(a_eval > 2 * 0.9);
		CHECK(b_eval > 3 * 0.9);
		CHECK(tab_eval == Approx(t_eval + a_eval + b_eval).epsilon(1e-3));
	}

	// Evaluate the base pairing of a state that shouldn't exist.
	
	SECTION("evaluate non-existent base pairs") {
		BasePairingTerm cd_term(
				ConditionEnum::HOLO,	// switch condition
				{"aptamer/c", "aptamer/d"},
				{"decoy/a", "decoy/b"}
		);
		BasePairingTerm tab_term(
				ConditionEnum::APO,	// switch condition
				{"target", "aptamer/a", "aptamer/b"},
				{"aptamer/d", "aptamer/e", "decoy/b"}
		);
		BasePairingTerm c_term(
				ConditionEnum::APO,
				{"aptamer/c"},
				{"aptamer/d"}
		);

		double cd_eval = cd_term.evaluate(sensor, no_lig_fold, lig_fold);
		double tab_eval = tab_term.evaluate(sensor, no_lig_fold, lig_fold);
		double c_eval = c_term.evaluate(sensor, no_lig_fold, lig_fold);

		CHECK(cd_eval < 0.1);
		CHECK(tab_eval < 0.12);
		CHECK(c_eval < 0.1);
	}
}

TEST_CASE("Test the 'ligand sensitivity' term", "[scoring]") {
	// Make a term that will score the "eval/a" and "eval/b" domains.  All the 
	// tests in this case are constructed to have these domains.
	LigandSensitivityTerm term({"eval/a", "eval/b"});

	SECTION("sequences that are always base-paired score poorly") {
		ConstructPtr rna = make_shared<Construct>();
		*rna += make_shared<Domain>("eval/a", "ACGU");
		*rna += make_shared<Domain>("ignore", "GAAA");
		*rna += make_shared<Domain>("eval/b", "ACGU");

		RnaFold no_lig_fold(rna, LigandEnum::NONE);
		RnaFold lig_fold(rna, LigandEnum::THEO);

		CHECK(term.evaluate(rna, no_lig_fold, lig_fold) == Approx(0.0));
	}

	SECTION("sequences that are never base-paired score poorly") {
		ConstructPtr rna = make_shared<Construct>();
		*rna += make_shared<Domain>("eval/a", "UUUUUU");
		*rna += make_shared<Domain>("eval/b", "UUUUUU");

		RnaFold no_lig_fold(rna, LigandEnum::NONE);
		RnaFold lig_fold(rna, LigandEnum::THEO);

		CHECK(term.evaluate(rna, no_lig_fold, lig_fold) == Approx(0.0));
	}

	SECTION("sequences that are affected by the ligand score well") {
		ConstructPtr rna = make_shared<Construct>();
		*rna += make_shared<Domain>("eval/a", "G");
		*rna += make_shared<Domain>("aptamer", "AUACCAGCCGAAAGGCCCUUGGCAG");
		*rna += make_shared<Domain>("eval/b", "C");

		RnaFold no_lig_fold(rna, LigandEnum::NONE);
		RnaFold lig_fold(rna, LigandEnum::THEO);

		CHECK(term.evaluate(rna, no_lig_fold, lig_fold) > log(1 * 50));
	}

	SECTION("rhf(6) scores well") {
		ConstructPtr rna = make_shared<Construct>();
		*rna += make_shared<Domain>("lower_stem/a", "guuuua");
		*rna += make_shared<Domain>("upper_stem", "gagcuagaaauagcaagu");
		*rna += make_shared<Domain>("lower_stem/b", "uaaaau");
		*rna += make_shared<Domain>("nexus/x", "aagg");
		*rna += make_shared<Domain>("nexus/y", "cuagu");
		*rna += make_shared<Domain>("nexus/z", "ccCuU");
		*rna += make_shared<Domain>("ruler", "UUC");
		*rna += make_shared<Domain>("hairpin/u", "GCC");
		*rna += make_shared<Domain>("aptamer", "gauaccagccgaaaggcccuuggcagc");
		*rna += make_shared<Domain>("hairpin/v", "GAC");
		*rna += make_shared<Domain>("tail", "ggcaccgagucggugcuuuuuu");

		RnaFold no_lig_fold(rna, LigandEnum::NONE);
		RnaFold lig_fold(rna, LigandEnum::THEO);

		cout << *rna << endl;
		cout << no_lig_fold.base_pair_string() << endl;
		cout << lig_fold.base_pair_string() << endl;

		LigandSensitivityTerm xz_term({"nexus/x", "nexus/z"});
		CHECK(xz_term.evaluate(rna, no_lig_fold, lig_fold) > log(4 * 2.5));
	}
}


*/
