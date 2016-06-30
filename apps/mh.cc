#include <cmath>
#include <iostream>
#include <string>

extern "C" {
#include <ViennaRNA/utils.h>
#include <ViennaRNA/fold.h>
}

#include <docopt/docopt.h>

#include <sgrna_design/model.hh>
#include <sgrna_design/sampling.hh>
#include <sgrna_design/scoring.hh>
#include <sgrna_design/utils.hh>

using namespace std;
using namespace sgrna_design;

static char const USAGE[] = R"""(
Insert the theophylline aptamer into the first hairpin of the sgRNA and run a 
Monte Carlo design simulation to optimize variable regions in the nexus, the 
ruler, and the rest of the first hairpin.  The design goal is to only form the 
wildtype nexus and hairpin base pairs when theophylline is bound.

Usage:
  mh [options]

Options:
  -n --num-moves <num>   [default=100]
    The number of moves to attempt in the Monte Carlo simulation.  I haven't 
    yet determined how many moves are required to reach convergence.

  -t --kt <beta>         [default=1.0]
    The likelihood of accepting a negative move.  If kT=0, only positive moves 
    will be accepted.  In the limit that kT=inf, every move will be accepted.  
    This parameter should be tuned to the magnitude of the score function.

  -v --version
    Display the version of ``mh`` being used.

  -h --help
    Display this usage information.
)""";

enum class ScorefxnEnum {
	GENERAL,
	SPECIFIC,
};


ConstructPtr
build_rhf_6_sgrna() {
	ColorEnum GREEN = ColorEnum::GREEN;
	ColorEnum RED = ColorEnum::RED;
	ColorEnum MAGENTA = ColorEnum::MAGENTA;
	ColorEnum BLUE = ColorEnum::BLUE;
	ColorEnum YELLOW = ColorEnum::YELLOW;
	StyleEnum BOLD = StyleEnum::BOLD;

	ConstructPtr sgrna = make_shared<Construct>();

	*sgrna += make_shared<Domain>("spacer", "");
	*sgrna += make_shared<Domain>("lower_stem/a", "guuuua", GREEN);
	*sgrna += make_shared<Domain>("upper_stem", "gagcuagaaauagcaagu", GREEN);
	*sgrna += make_shared<Domain>("lower_stem/b", "uaaaau", GREEN);
	*sgrna += make_shared<Domain>("nexus/a", "aagg", RED);
	*sgrna += make_shared<Domain>("nexus/b", "cuagu", RED, BOLD);
	*sgrna += make_shared<Domain>("nexus/c", "cc", RED);
	*sgrna += make_shared<Domain>("ruler", "CuUUUC", MAGENTA, BOLD);
	*sgrna += make_shared<Domain>("hairpin/a", "GCC", BLUE, BOLD);
	*sgrna += make_shared<Domain>("aptamer", "gauaccagccgaaaggcccuuggcagc", YELLOW);
	*sgrna += make_shared<Domain>("hairpin/b", "GAC", BLUE, BOLD);
	*sgrna += make_shared<Domain>("tail", "ggcaccgagucggugcuuuuuu", BLUE);

	return sgrna;
}

ConstructPtr
build_mh_sgrna() {
	ColorEnum GREEN = ColorEnum::GREEN;
	ColorEnum RED = ColorEnum::RED;
	ColorEnum MAGENTA = ColorEnum::MAGENTA;
	ColorEnum BLUE = ColorEnum::BLUE;
	ColorEnum YELLOW = ColorEnum::YELLOW;
	StyleEnum BOLD = StyleEnum::BOLD;

	ConstructPtr sgrna = make_shared<Construct>();

	*sgrna += make_shared<Domain>("spacer", "");
	*sgrna += make_shared<Domain>("lower_stem/a", "guuuua", GREEN);
	*sgrna += make_shared<Domain>("upper_stem", "gagcuagaaauagcaagu", GREEN);
	*sgrna += make_shared<Domain>("lower_stem/b", "uaaaau", GREEN);
	*sgrna += make_shared<Domain>("nexus/a", "aagg", RED);
	*sgrna += make_shared<Domain>("nexus/b", "CUAGU", RED, BOLD);
	*sgrna += make_shared<Domain>("nexus/c", "cc", RED);
	*sgrna += make_shared<Domain>("ruler", "GUUAUCA", MAGENTA, BOLD);
	*sgrna += make_shared<Domain>("hairpin/a", "ACU", BLUE, BOLD);
	*sgrna += make_shared<Domain>("aptamer", "gauaccagccgaaaggcccuuggcagc", YELLOW);
	*sgrna += make_shared<Domain>("hairpin/b", "AGU", BLUE, BOLD);
	*sgrna += make_shared<Domain>("tail", "ggcaccgagucggugcuuuuuu", BLUE);

	return sgrna;
}

ScoreFunctionPtr
build_mh_scorefxn(ScorefxnEnum style=ScorefxnEnum::SPECIFIC) {
	ScoreFunctionPtr scorefxn = make_shared<ScoreFunction>();

	switch(style) {

		case ScorefxnEnum::GENERAL:
			// Fold differently with the ligand than without it.

			*scorefxn += ScoreTermPtr(new LigandSensitivityTerm(
					{"nexus/a", "nexus/b", "nexus/c", "ruler", "hairpin/a", "hairpin/b"},
					1.0
					));

			break;

		case ScorefxnEnum::SPECIFIC:
			// Desired base pairs *without* ligand.

			*scorefxn += ScoreTermPtr(new SpecificLigandSensitivityTerm(
					ConditionEnum::APO,
					{"nexus/a", "nexus/b", "nexus/c"},
					{"ruler", "hairpin/a", "aptamer", "hairpin/b"}
			));

			*scorefxn += ScoreTermPtr(new SpecificLigandSensitivityTerm(
					ConditionEnum::APO,
					{"hairpin/a"},
					{"nexus/a", "nexus/b", "nexus/c", "ruler"}
			));

			*scorefxn += ScoreTermPtr(new SpecificLigandSensitivityTerm(
					ConditionEnum::APO,
					{"hairpin/b"},
					{"aptamer"}
			));

			// Desired base pairs *with* ligand.

			*scorefxn += ScoreTermPtr(new SpecificLigandSensitivityTerm(
					ConditionEnum::HOLO,
					{"hairpin/a"},
					{"hairpin/b"}
			));

			// Desired base pairs *with and without* ligand.
			//
			//*scorefxn += ScoreTermPtr(new RequiredBasePairing(
			//		{"lower_stem/a"},
			//		{"lower_stem/b"},
			//		1.0 / 6.0
			//));

			break;
	}

	//*scorefxn += StableFoldTerm();

	return scorefxn;
}

MonteCarloPtr
build_mh_sampler(ConstructPtr wt) {
	MonteCarloPtr sampler = make_shared<MonteCarlo>();
	vector<string> mutable_domains = {
		"nexus/b", "ruler", "hairpin/a", "hairpin/b"
	};
	
	*sampler += make_shared<MakePointMutation>(mutable_domains);
	*sampler += make_shared<ChangeDomainLength>("ruler", 4, 7);
	//*sampler += make_shared<MakeWtReversion>(mutable_domains, wt);

	return sampler;
}


int main(int argc, char **argv) {
	map<string, docopt::value> args = docopt::docopt(
			USAGE+1, {argv + 1, argv + argc}, true, "0.0");

	ConstructPtr mh = build_mh_sgrna();
	ConstructPtr rhf_6 = build_rhf_6_sgrna();
	ConstructPtr wt = mh->copy();
	ScoreFunctionPtr scorefxn = build_mh_scorefxn();
	MonteCarloPtr sampler = build_mh_sampler(wt);

	sampler->scorefxn(scorefxn);
	sampler->num_steps(stoi(args["--num-moves"].asString()));
	//sampler->beta(stod(args["--kt"].asString()));
	sampler->beta(10);
	mh = sampler->apply(mh);

	cout << f("wt:     %s (%.4f)") % *wt % scorefxn->evaluate(wt) << endl;
	cout << f("rhf(6): %s (%.4f)") % *rhf_6 % scorefxn->evaluate(rhf_6) << endl;
	cout << f("mh:     %s (%.4f)") % *mh % scorefxn->evaluate(mh) << endl;

	return 0;
}


