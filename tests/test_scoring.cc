#include <set>
#include <vector>
#include <catch/catch.hpp>
#include "model.hh"
#include "scoring.hh"
#include "utils.hh"

using namespace std;
using namespace addapt;
using bp = pair<int,int>;

class DummyRnaFold : public RnaFold {

public:

	DummyRnaFold(double p=0): my_macrostate_prob(p) {}

	double &
	operator[](bp key) {
		key = {
			std::min(key.first, key.second),
			std::max(key.first, key.second)};
		return my_base_pair_probs[key];
	}

	double
	base_pair_prob(int a, int b) const {
		bp key = {std::min(a,b), std::max(a,b)}; // order doesn't matter
		auto it = my_base_pair_probs.find(key);
		return (it != my_base_pair_probs.end())? it->second : 0.0;
	}

	double
	macrostate_prob(string) const {
		return my_macrostate_prob;
	}

private:

	map<bp,double> my_base_pair_probs;
	double my_macrostate_prob;

};

DeviceConstPtr
build_rhf_6_device() {
	// 0....,....1....,....2....,....3....,....4....,....5....,....6....,....7....,....8....,....9....,....0.
	// GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCCUUUUCGCCGAUACCAGCCGAAAGGCCCUUGGCAGCGACGGCACCGAGUCGGUGCUUUUUU
	// (((((((.((((....))))...))))))).,,({{,{..|||{{(,((,{....,.||{}}}})),..,}))).,,||.(((((((...)))))))..... [ΔG=-29.58 kcal/mol, apo]
	// (((((((.((((....))))...))))))){(((......)))}..{{.((...((.(((....)))....))...)).}|((((((...)))))),..... [ΔG=-33.82 kcal/mol, holo]

	string seq = "guuuuagagcuagaaauagcaaguuaaaauaaggcuaguccCuUUUCGCCgauaccagccgaaaggcccuuggcagcGACggcaccgagucggugcuuuuuu";
	string cst = "(............................)xx..xxxxx..xxxxxx(...............................).(.............)......";

	DevicePtr rhf_6 = make_shared<Device>(seq);
	rhf_6->add_macrostate("active", cst);
	return rhf_6;
}

AptamerConstPtr THEO_APTAMER = make_shared<Aptamer>(
  "GAUACCAGCCGAAAGGCCCUUGGCAGC",
  "(...((.(((....)))....))...)",
  0.32);

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
	// Make a device that should fold into a hairpin:
	
	// ACGUGAAAACGU
	// ((((....)))) [ΔG=-2.20 kcal/mol]

	DevicePtr hairpin = make_shared<Device>("ACGUGAAAACGU");

	// Indicate which base pairs I expect to find.
	set<bp> expected_base_pairs = {{0,11}, {1,10}, {2,9}, {3,8}};
	map<bp,double> base_pair_probs = {
		{{0,11}, 0.70},
		{{1,10}, 0.95},
		{{2, 9}, 0.95},
		{{3, 8}, 0.95}
	};

	ViennaRnaFold fold(hairpin);

	SECTION("the expected base pairs probabilities are correct") {
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

	SECTION("the order of the base pairs doesn't matter") {
		for(int i = 0; i < hairpin->len(); i++) {
			for(int j = i; j < hairpin->len(); j++) {
				double ij = fold.base_pair_prob(i, j);
				double ji = fold.base_pair_prob(j, i);
				REQUIRE(ij == ji);
			}
		}
	}

	SECTION("out-of-bounds base pair indices throw exceptions") {
		CHECK_THROWS(fold.base_pair_prob(0, 12));
		CHECK_THROWS(fold.base_pair_prob(0, -13));
	}

	SECTION("the macrostate probabilities are correct") {
		CHECK(fold.macrostate_prob("...(....)...") > 0.95);
		CHECK(fold.macrostate_prob("..((....))..") > 0.85);
		CHECK(fold.macrostate_prob(".(((....))).") > 0.75);
		CHECK(fold.macrostate_prob("((((....))))") > 0.65);
		CHECK(fold.macrostate_prob("(((......)))") > 0.65);
		CHECK(fold.macrostate_prob("((........))") > 0.65);
		CHECK(fold.macrostate_prob("(..........)") > 0.65);

		CHECK(fold.macrostate_prob("xxxxxxxxxxxx") < 0.05);
		CHECK(fold.macrostate_prob("xxxx........") < 0.05);
		CHECK(fold.macrostate_prob("....xxxx....") > 0.95);
		CHECK(fold.macrostate_prob("........xxxx") < 0.05);
	}
}

TEST_CASE("Test folding a hairpin with an aptamer", "[scoring]") {
	// Make a device consisting entirely of an aptamer:

	// GAUACCAGCCGAAAGGCCCUUGGCAGC
	// ....((((((....)))...))).... [ΔG=-6.20 kcal/mol, apo]
	// (...((.(((....)))....))...) [ΔG=-9.22 kcal/mol, holo]

	DevicePtr hairpin = make_shared<Device>("GAUACCAGCCGAAAGGCCCUUGGCAGC");

	// Indicate which base pairs I expect to form.
	set<bp> apo_base_pairs = {
		{4,22}, {5,21}, {6,20}, {7,16}, {8,15}, {9,14}
	};
	set<bp> holo_base_pairs = {
		{0,26}, {4,22}, {5,21}, {7,16}, {8,15}, {9,14}
	};

	// Make sure the expected base pairs form.
	ViennaRnaFold apo_fold(hairpin);
	ViennaRnaFold holo_fold(hairpin, THEO_APTAMER);

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

	// Make sure the various macrostates have the expected populations.
	CHECK(apo_fold.macrostate_prob("....((((((....)))...)))....") > 0.65);
	CHECK(apo_fold.macrostate_prob("(.........................)") < 0.05);
	CHECK(holo_fold.macrostate_prob("....((((((....)))...)))....") < 0.01);
	CHECK(holo_fold.macrostate_prob("(.........................)") > 0.95);
}

TEST_CASE("Test folding rhf(6)", "[scoring]") {
	DeviceConstPtr rhf_6 = build_rhf_6_device();

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
	ViennaRnaFold apo_fold(rhf_6);
	ViennaRnaFold holo_fold(rhf_6, THEO_APTAMER);

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
	
	// Make sure the various macrostates have the expected populations.
	CAPTURE(rhf_6->seq());
	CAPTURE(rhf_6->macrostate("active"));
	CHECK(apo_fold.macrostate_prob(rhf_6->macrostate("active")) < 7e-5);
	CHECK(holo_fold.macrostate_prob(rhf_6->macrostate("active")) > 4e-3);
}

TEST_CASE("Test the score function class", "[scoring]") {
	ScoreFunction scorefxn;
	DevicePtr dummy_device = make_shared<Device>("UUUU");

	class DummyTerm : public ScoreTerm {

	public:

		DummyTerm(double score, double weight):
			ScoreTerm("dummy", weight), my_score(score) {}

		double
		evaluate(DeviceConstPtr, RnaFold const &, RnaFold const &) const {
			return my_score;
		}

	private:
		double my_score;

	};

	CHECK(scorefxn.evaluate(dummy_device) == 0);

	SECTION("single score terms are summed correctly") {
		scorefxn += make_shared<DummyTerm>(10, 1);
		CHECK(scorefxn.evaluate(dummy_device) == Approx(10));
	}

	SECTION("multiple score terms are summed correctly") {
		scorefxn += make_shared<DummyTerm>(10, 1);
		scorefxn += make_shared<DummyTerm>(5, 1);
		CHECK(scorefxn.evaluate(dummy_device) == Approx(15));
	}

	SECTION("weights are applied correctly") {
		scorefxn += make_shared<DummyTerm>(10, 1);
		scorefxn += make_shared<DummyTerm>(10, 0.5);
		CHECK(scorefxn.evaluate(dummy_device) == Approx(15));
	}
}

TEST_CASE("Test the 'macrostate prob' score term", "[scoring]") {
	DevicePtr dummy_device = make_shared<Device>("");
	dummy_device->add_macrostate("dummy", "");

	struct Test {
		string name;
		ConditionEnum condition;
		FavorableEnum	favorable;
		double apo_prob;
		double holo_prob;
		double expected_score;
	};

	vector<Test> tests = {
		{"apo: not dummy",  ConditionEnum::APO,  FavorableEnum::NO,  0.2, 0.2, log(0.8)},
		{"apo: not dummy",  ConditionEnum::APO,  FavorableEnum::NO,  0.9, 0.2, log(0.1)},
		{"apo: not dummy",  ConditionEnum::APO,  FavorableEnum::NO,  0.2, 0.9, log(0.8)},
		{"apo: not dummy",  ConditionEnum::APO,  FavorableEnum::NO,  0.9, 0.9, log(0.1)},

		{"apo: dummy",      ConditionEnum::APO,  FavorableEnum::YES, 0.2, 0.2, log(0.2)},
		{"apo: dummy",      ConditionEnum::APO,  FavorableEnum::YES, 0.9, 0.2, log(0.9)},
		{"apo: dummy",      ConditionEnum::APO,  FavorableEnum::YES, 0.2, 0.9, log(0.2)},
		{"apo: dummy",      ConditionEnum::APO,  FavorableEnum::YES, 0.9, 0.9, log(0.9)},

		{"holo: not dummy", ConditionEnum::HOLO, FavorableEnum::NO,  0.2, 0.2, log(0.8)},
		{"holo: not dummy", ConditionEnum::HOLO, FavorableEnum::NO,  0.9, 0.2, log(0.8)},
		{"holo: not dummy", ConditionEnum::HOLO, FavorableEnum::NO,  0.2, 0.9, log(0.1)},
		{"holo: not dummy", ConditionEnum::HOLO, FavorableEnum::NO,  0.9, 0.9, log(0.1)},

		{"holo: dummy",     ConditionEnum::HOLO, FavorableEnum::YES, 0.2, 0.2, log(0.2)},
		{"holo: dummy",     ConditionEnum::HOLO, FavorableEnum::YES, 0.9, 0.2, log(0.2)},
		{"holo: dummy",     ConditionEnum::HOLO, FavorableEnum::YES, 0.2, 0.9, log(0.9)},
		{"holo: dummy",     ConditionEnum::HOLO, FavorableEnum::YES, 0.9, 0.9, log(0.9)},
	};

	for(Test test: tests) {
		DummyRnaFold dummy_apo_fold(test.apo_prob);
		DummyRnaFold dummy_holo_fold(test.holo_prob);

		CAPTURE(test.condition)
		CAPTURE(test.favorable)
		CAPTURE(test.apo_prob)
		CAPTURE(test.holo_prob)

		MacrostateProbTerm term(
				"dummy", test.condition, test.favorable);
		double score = term.evaluate(
				dummy_device, dummy_apo_fold, dummy_holo_fold);

		CHECK(term.name() == test.name);
		CHECK(score == Approx(test.expected_score));
	}
}

