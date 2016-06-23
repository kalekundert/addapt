#include <cmath>
#include <iostream>

extern "C" {
#include <ViennaRNA/utils.h>
#include <ViennaRNA/fold.h>
}

#include "model.hh"
#include "sampling.hh"
#include "scoring.hh"
#include "utils.hh"

using namespace std;
using namespace sgrna;

ColorEnum GREEN = ColorEnum::GREEN;
ColorEnum RED = ColorEnum::RED;
ColorEnum MAGENTA = ColorEnum::MAGENTA;
ColorEnum BLUE = ColorEnum::BLUE;
ColorEnum YELLOW = ColorEnum::YELLOW;
StyleEnum BOLD = StyleEnum::BOLD;

ConstructPtr
build_mh_sgrna() {
	ConstructPtr sgrna = make_shared<Construct>();

	sgrna += make_shared<Domain>("spacer", "");
	sgrna += make_shared<Domain>("stem", "GUUUUAGAGCUAGAAAUAGCAAGUUAAAAU", GREEN);
	sgrna += make_shared<Domain>("nexus/a", "AAGG", RED);
	sgrna += make_shared<Domain>("nexus/b", "CUAGU", RED, BOLD);
	sgrna += make_shared<Domain>("nexus/c", "CC", RED);
	sgrna += make_shared<Domain>("ruler", "GUUAUCA", MAGENTA, BOLD);
	sgrna += make_shared<Domain>("hairpin/a", "ACUU", BLUE, BOLD);
	sgrna += make_shared<Domain>("aptamer", "AUACCAGCCGAAAGGCCCUUGGCAG", YELLOW);
	sgrna += make_shared<Domain>("hairpin/b", "AAGU", BLUE, BOLD);
	sgrna += make_shared<Domain>("tail", "GGCACCGAGUCGGUGCUUUUUU", BLUE);

	return sgrna;
}

ScoreFunctionPtr
build_mh_scorefxn() {
	ScoreFunctionPtr scorefxn = make_shared<ScoreFunction>();
	
	// Add score terms...

	return scorefxn;
}

MonteCarloPtr
build_mh_sampler(ConstructPtr wt) {
	MonteCarloPtr sampler = make_shared<MonteCarlo>();
	
	// How to decide which domains?
	sampler += MakePointMutation();
	// Add moves...

	return sampler;
}


int main(int argc, char **argv) {
	ConstructPtr mh = build_mh_sgrna();
	ConstructPtr wt = mh->copy();
	ScoreFunctionPtr scorefxn = build_mh_scorefxn();
	MonteCarloPtr sampler = build_mh_sampler(wt);

	sampler->scorefxn(scorefxn);
	sampler->apply(mh);

	cout << f("wt: %1%") % *wt << endl;
	cout << f("mh: %1%") % *mh << endl;

	return 0;
}


