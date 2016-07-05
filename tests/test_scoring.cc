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

TEST_CASE("Test folding rhf(6)", "[scoring]") {
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

TEST_CASE("Test the score function class", "[scoring]") {
	ScoreFunction scorefxn;
	ConstructPtr dummy_construct = make_shared<Construct>();
	*dummy_construct += make_shared<Domain>("dummy", "UUUU");

	class DummyTerm : public ScoreTerm {

	public:

		DummyTerm(double score, double weight):
			ScoreTerm("dummy", weight), my_score(score) {}

		double
		evaluate(ConstructConstPtr, RnaFold const &, RnaFold const &) const {
			return my_score;
		}

	private:
		double my_score;

	};

	CHECK(scorefxn.evaluate(dummy_construct) == 0);

	SECTION("single score terms are summed correctly") {
		scorefxn += make_shared<DummyTerm>(10, 1);
		CHECK(scorefxn.evaluate(dummy_construct) == Approx(10));
	}

	SECTION("multiple score terms are summed correctly") {
		scorefxn += make_shared<DummyTerm>(10, 1);
		scorefxn += make_shared<DummyTerm>(5, 1);
		CHECK(scorefxn.evaluate(dummy_construct) == Approx(15));
	}

	SECTION("weights are applied correctly") {
		scorefxn += make_shared<DummyTerm>(10, 1);
		scorefxn += make_shared<DummyTerm>(10, 0.5);
		CHECK(scorefxn.evaluate(dummy_construct) == Approx(15));
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

	dummy_apo_fold[{0,3}] = 0.20;
	dummy_holo_fold[{0,3}] = 0.40;

	dummy_apo_fold[{1,2}] = 0.30;
	dummy_holo_fold[{1,2}] = 0.10;

	// Define expected scores for various combinations of nucleotides.
	struct Test {
		vector<string> selection;
		double expected_score;
	};

	vector<Test> tests = {
		// The score increases with the number of ligand-sensitive base pairs.
		{{"a","b"}, (0.3-0.1)*log(3)/3},
		{{"a","c"}, (0.4-0.2)*log(2)/3},
		{{"a","b","c"}, (0.3-0.1)*log(3)/4 + (0.4-0.2)*log(2)/4},

		// The score is 0 if there aren't any base pairs.
		{{},    0},
		{{"a"}, 0},
		{{"b"}, 0},
		{{"c"}, 0}
	};

	// Each scenario is scored correctly.
	for(Test test: tests) {
		LigandSensitivityTerm term("dummy", test.selection);
		double score = term.evaluate(
				dummy_construct, dummy_apo_fold, dummy_holo_fold);

		CAPTURE(test.selection);
		CHECK(score == Approx(test.expected_score));
	}
}

TEST_CASE("Test the 'conditionally paired' score term", "[scoring]") {
	// Make a construct with individually indexable nucleotides.
	ConstructPtr dummy_construct = make_shared<Construct>();
	*dummy_construct += make_shared<Domain>("a", "U");
	*dummy_construct += make_shared<Domain>("A", "U");
	*dummy_construct += make_shared<Domain>("b", "U");
	*dummy_construct += make_shared<Domain>("B", "U");

	// Define fake base-pairing probabilities for this construct.
	DummyRnaFold dummy_apo_fold;
	DummyRnaFold dummy_holo_fold;

	dummy_apo_fold[{0,1}] = 0.20;
	dummy_holo_fold[{0,1}] = 0.40;

	dummy_apo_fold[{2,3}] = 0.30;
	dummy_holo_fold[{2,3}] = 0.10;

	// Define expected scores for various combinations of nucleotides.
	struct Test {
		vector<string> selection;
		vector<string> targets;
		double expected_score;
	};

	vector<Test> tests = {
		// The score is proportional to the ligand sensitivity of the base pair.
		{{"a"}, {"A"},  0.6*log(2)},
		{{"b"}, {"B"}, -0.4*log(3)},

		// The score is normalized by the number of nucleotides in the selection.
		{{"a"},     {"A","B"},  0.6*log(2)/1},
		{{"b"},     {"A","B"}, -0.4*log(3)/1},
		{{"a","b"}, {"A"},      0.6*log(2)/2},
		{{"a","b"}, {"B"},     -0.4*log(3)/2},
		{{"a","b"}, {"A","B"},  0.6*log(2)/2 -0.4*log(3)/2},

		// The score is 0 if there aren't any base pairs.
		{{},    {},    0},
		{{"a"}, {},    0},
		{{"a"}, {"a"}, 0},
		{{"a"}, {"B"}, 0},
		{{"b"}, {},    0},
		{{"b"}, {"b"}, 0},
		{{"b"}, {"A"}, 0},
	};

	// Each scenario is scored correctly in both the apo and holo conditions.
	for(Test test: tests) {
		ConditionallyPairedTerm apo_term(
				"dummy", ConditionEnum::APO, test.selection, test.targets);
		ConditionallyPairedTerm holo_term(
				"dummy", ConditionEnum::HOLO, test.selection, test.targets);

		double apo_score = apo_term.evaluate(
				dummy_construct, dummy_apo_fold, dummy_holo_fold);
		double holo_score = holo_term.evaluate(
				dummy_construct, dummy_apo_fold, dummy_holo_fold);

		CAPTURE(test.selection);
		CAPTURE(test.targets);

		CHECK(apo_score == Approx(-test.expected_score));
		CHECK(holo_score == Approx(test.expected_score));
	}
}

TEST_CASE("Test the 'conditionally unpaired' score term", "[scoring]") {
	// Make a construct with individually indexable nucleotides.
	ConstructPtr dummy_construct = make_shared<Construct>();
	*dummy_construct += make_shared<Domain>("a", "U");
	*dummy_construct += make_shared<Domain>("A", "U");
	*dummy_construct += make_shared<Domain>("b", "U");
	*dummy_construct += make_shared<Domain>("B", "U");

	// Define fake base-pairing probabilities for this construct.
	DummyRnaFold dummy_apo_fold;
	DummyRnaFold dummy_holo_fold;

	dummy_apo_fold[{0,1}]  = 0.8; // p_unpaired(0, apo)  = 0.2
	dummy_holo_fold[{0,1}] = 0.6; // p_unpaired(0, holo) = 0.4

	dummy_apo_fold[{2,3}]  = 0.6; // p_unpaired(2, apo)  = 0.3
	dummy_apo_fold[{2,1}]  = 0.1;
	dummy_holo_fold[{2,3}] = 0.6; // p_unpaired(2, holo) = 0.1
	dummy_holo_fold[{2,1}] = 0.3;

	// Define expected scores for various combinations of nucleotides.
	struct Test {
		vector<string> selection;
		double expected_score;
	};

	vector<Test> tests = {
		{{},         0},
		{{"a"},      0.6*log(2)},
		{{"b"},     -0.4*log(3)},
		{{"a","b"},  0.6*log(2)/2 -0.4*log(3)/2}
	};

	// Each scenario is scored correctly in both the apo and holo conditions.
	for(Test test: tests) {
		ConditionallyUnpairedTerm apo_term(
				"dummy", ConditionEnum::APO, test.selection);
		ConditionallyUnpairedTerm holo_term(
				"dummy", ConditionEnum::HOLO, test.selection);

		double apo_score = apo_term.evaluate(
				dummy_construct, dummy_apo_fold, dummy_holo_fold);
		double holo_score = holo_term.evaluate(
				dummy_construct, dummy_apo_fold, dummy_holo_fold);

		CAPTURE(test.selection);

		CHECK(apo_score == Approx(-test.expected_score));
		CHECK(holo_score == Approx(test.expected_score));
	}
}

TEST_CASE("Test the 'always paired' score term", "[scoring]") {
	// Make a construct with individually indexable nucleotides.
	ConstructPtr dummy_construct = make_shared<Construct>();
	*dummy_construct += make_shared<Domain>("a", "U");
	*dummy_construct += make_shared<Domain>("A", "U");
	*dummy_construct += make_shared<Domain>("b", "U");
	*dummy_construct += make_shared<Domain>("B", "U");
	*dummy_construct += make_shared<Domain>("c", "U");
	*dummy_construct += make_shared<Domain>("C", "U");

	// Define fake base-pairing probabilities for this construct.
	DummyRnaFold dummy_apo_fold;
	DummyRnaFold dummy_holo_fold;

	dummy_apo_fold[{0,1}] = 0.50;
	dummy_holo_fold[{0,1}] = 0.25;

	dummy_apo_fold[{2,3}] = 0.80;
	dummy_holo_fold[{2,3}] = 0.90;

	dummy_apo_fold[{3,4}] = 1.00;
	dummy_holo_fold[{3,4}] = 1.00;

	// Define expected scores for various combinations of nucleotides.
	struct Test {
		vector<string> selection;
		vector<string> targets;
		double expected_score;
	};

	vector<Test> tests = {
		// The score indicates how well the base pair forms in both conditions.
		{{"a"}, {"A"}, -log(1*3)/2},
		{{"b"}, {"B"}, +log(4*9)/2},

		// The score is normalized by the number of nucleotides in the selection.
		{{"a"},     {"A","B"}, -log(1*3)/2},
		{{"b"},     {"A","B"}, +log(4*9)/2},
		{{"a","b"}, {"A"},     -log(1*3)/4},
		{{"a","b"}, {"B"},     +log(4*9)/4},
		{{"a","b"}, {"A","B"}, -log(1*3)/4 +log(4*9)/4},

		// The score is 0 if there aren't any base pairs.
		{{},    {},    0},
		{{"a"}, {},    0},
		{{"a"}, {"a"}, 0},
		{{"a"}, {"B"}, 0},
		{{"b"}, {},    0},
		{{"b"}, {"b"}, 0},
		{{"b"}, {"A"}, 0},

		// Counter-intuitively, the score is also zero if no off-targets are 
		// possible, because otherwise it would be infinite in this case.
		{{"c"}, {"C"}, 0}
	};

	// Each scenario is scored correctly in both the apo and holo conditions.
	for(Test test: tests) {
		AlwaysPairedTerm term(
				"dummy", test.selection, test.targets);

		double score = term.evaluate(
				dummy_construct, dummy_apo_fold, dummy_holo_fold);

		CAPTURE(test.selection);
		CAPTURE(test.targets);

		CHECK(score == Approx(test.expected_score));
	}
}

TEST_CASE("Test the 'always unpaired' score term", "[scoring]") {
	// Make a construct with individually indexable nucleotides.
	ConstructPtr dummy_construct = make_shared<Construct>();
	*dummy_construct += make_shared<Domain>("a", "U");
	*dummy_construct += make_shared<Domain>("A", "U");
	*dummy_construct += make_shared<Domain>("b", "U");
	*dummy_construct += make_shared<Domain>("B", "U");

	// Define fake base-pairing probabilities for this construct.
	DummyRnaFold dummy_apo_fold;
	DummyRnaFold dummy_holo_fold;

	dummy_apo_fold[{0,1}]  = 0.50; // p_unpaired(0, apo)  = 0.50
	dummy_holo_fold[{0,1}] = 0.75; // p_unpaired(0, holo) = 0.25

	dummy_apo_fold[{2,3}]  = 0.10; // p_unpaired(2, apo)  = 0.80
	dummy_apo_fold[{2,1}]  = 0.10;
	dummy_holo_fold[{2,3}] = 0.05; // p_unpaired(2, holo) = 0.90
	dummy_holo_fold[{2,1}] = 0.05;

	// Define expected scores for various combinations of nucleotides.
	struct Test {
		vector<string> selection;
		double expected_score;
	};

	vector<Test> tests = {
		{{},         0},
		{{"a"},     -log(1*3)/2},
		{{"b"},     +log(4*9)/2},
		{{"a","b"}, -log(1*3)/4 +log(4*9)/4}
	};

	// Each scenario is scored correctly in both the apo and holo conditions.
	for(Test test: tests) {
		AlwaysUnpairedTerm apo_term("dummy", test.selection);
		double score = apo_term.evaluate(
				dummy_construct, dummy_apo_fold, dummy_holo_fold);
		CAPTURE(test.selection);
		CHECK(score == Approx(test.expected_score));
	}
}


