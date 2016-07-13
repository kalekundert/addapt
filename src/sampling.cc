#include <algorithm>
#include <cassert>
#include <cctype>
#include <cmath>
#include <iostream>
#include <regex>

#include "sampling.hh"
#include "utils.hh"

using namespace std;

namespace addapt {

MonteCarlo::MonteCarlo(): 
	my_steps(0),
	my_thermostat(std::make_shared<FixedThermostat>(1)),
	my_scorefxn(std::make_shared<ScoreFunction>()),
	my_moves(),
	my_reporters() {}

ConstructPtr
MonteCarlo::apply(ConstructPtr construct, std::mt19937 &rng) const {
	if (my_moves.empty()) {
		return construct;
	}

	// Setup to data structure that will hold all the information about each 
	// step.  The purpose of this structure is to support external logging 
	// methods and thermostats.
	MonteCarloStep step;
	step.num_steps = my_steps; step.i = -1;
	step.current_construct = construct;

	// Create random number distributions for later use.
	auto random = std::bind(std::uniform_real_distribution<>(), rng);
	auto randmove = std::bind(std::uniform_int_distribution<>(0, my_moves.size()-1), rng);

	// Get an initial score.
	step.current_score = my_scorefxn->evaluate(step.current_construct, step.score_table);
	step.proposed_score = step.current_score;

	// Initialize the counters that will keep track of how often moves are 
	// accepted and rejected.
	step.outcome_counters[OutcomeEnum::REJECT] = 0;
	step.outcome_counters[OutcomeEnum::ACCEPT_WORSENED] = 0;
	step.outcome_counters[OutcomeEnum::ACCEPT_UNCHANGED] = 0;
	step.outcome_counters[OutcomeEnum::ACCEPT_IMPROVED] = 0;

	// Initialize the reporters
	for(auto reporter: my_reporters) {
		reporter->start(step);
	}

	for(step.i = 0; step.i < step.num_steps; step.i++) {
		// Get the temperature for the Metropolis criterion.  This has to be done 
		// every iteration, even if no accept/reject decision needs to be made.

		step.temperature = my_thermostat->adjust(step);

		// Copy the sgRNA so we can easily undo the move.
		step.proposed_construct = step.current_construct->copy();

		// Randomly pick a move to apply.
		step.move = my_moves[randmove()];
		step.move->apply(step.proposed_construct, rng);

		// Skip the score function evaluation if the sequence didn't change.
		if(step.current_construct->seq() == step.proposed_construct->seq()) {
			step.outcome = OutcomeEnum::ACCEPT_UNCHANGED;
		}

		// Score the proposed move, then either accept or reject it.
		else {
			step.proposed_score = my_scorefxn->evaluate(step.proposed_construct, step.score_table);
			step.score_diff = step.proposed_score - step.current_score;
			step.metropolis_criterion = std::exp(step.score_diff / step.temperature);
			step.random_threshold = random();

			if(step.metropolis_criterion < step.random_threshold) {
				step.outcome = OutcomeEnum::REJECT;
			}
			else{
				step.outcome = (step.score_diff > 0)?
					OutcomeEnum::ACCEPT_IMPROVED : OutcomeEnum::ACCEPT_WORSENED;

				step.current_construct = step.proposed_construct;
				step.current_score = step.proposed_score;
			}
		}

		// Update the accept/reject statistics.
		step.outcome_counters[step.outcome] += 1;

		// Give the reporters a chance to react to the move.
		for(auto reporter: my_reporters) {
			reporter->update(step);
		}
	}

	// Give the reporters one last chance to report things.
	for(auto reporter: my_reporters) {
		reporter->finish(step);
	}

	return step.current_construct;
}

int
MonteCarlo::num_steps() const {
	return my_steps;
}

void
MonteCarlo::num_steps(int num_steps) {
	my_steps = num_steps;
}

ThermostatPtr
MonteCarlo::thermostat() const {
	return my_thermostat;
}

void
MonteCarlo::thermostat(ThermostatPtr thermostat) {
	my_thermostat = thermostat;
}

ScoreFunctionPtr
MonteCarlo::scorefxn() const {
	return my_scorefxn;
}

void
MonteCarlo::scorefxn(ScoreFunctionPtr scorefxn) {
	my_scorefxn = scorefxn;
}

MoveList
MonteCarlo::moves() const {
	return my_moves;
}

void
MonteCarlo::add_move(MovePtr move) {
	my_moves.push_back(move);
}

void
MonteCarlo::operator+=(MovePtr move) {
	add_move(move);
}

ReporterList
MonteCarlo::reporters() const {
	return my_reporters;
}

void
MonteCarlo::add_reporter(ReporterPtr reporter) {
	my_reporters.push_back(reporter);
}

void
MonteCarlo::operator+=(ReporterPtr reporter) {
	add_reporter(reporter);
}


bool
can_be_mutated(ConstructConstPtr construct, int position) {
	// Only mutate positions that are upper case.  This is a simple way for the 
	// user to indicate which positions should be mutable.
	return isupper(construct->seq()[position]);
}

bool
can_be_freely_mutated(ConstructConstPtr construct, int position) {
	if(not can_be_mutated(construct, position)) {
			return false;
	}

	// Don't mutate positions that are constrained to be the 3' end of a base 
	// pair.  We don't want to over-sample these positions, and they'll be 
	// mutated as a unit with their partner.
	for(auto /*pair*/ macrostate: construct->macrostates()) {
		if(macrostate.second[position] == ')') {
			return false;
		}
	}

	return true;
}

void
mutate_recursively(
		ConstructPtr construct,
		int const position,
		char const mutation) {

	vector<bool> already_mutated(construct->len(), false);
	mutate_recursively(construct, position, mutation, already_mutated);
}

void
mutate_recursively(
		ConstructPtr construct,
		int const position,
		char const mutation,
		vector<bool> &already_mutated) {

	assert(already_mutated.size() == construct->len());
	assert(not already_mutated[position]);
	assert(can_be_mutated(construct, position));

	// Make the indicated mutation to the construct.
	construct->mutate(position, mutation);
	already_mutated[position] = true;

	// Recursively make complementary mutations in any position constrained to be 
	// base paired to this one.  GU base pairs are not considered because they 
	// would necessarily cause a bias in the final nucleotide distribution.
	for(auto item: construct->macrostates()) {
		string macrostate = item.second;
		char open_symbol, close_symbol;
		int step;

		// If this position is base-paired in this macrostate, setup the variables 
		// that will be used to search either forward or backward for its partner.
		if(macrostate[position] == '(') {
			open_symbol = '(';
			close_symbol = ')';
			step = 1;
		}

		else if(macrostate[position] == ')') {
			open_symbol = ')';
			close_symbol = '(';
			step = -1;
		}

		// Otherwise, if this position is not base-paired in this macrostate, move 
		// on to the next macrostate without doing anything.
		else {
			continue;
		}

		// Find the partner that this position is base-paired with in this 
		// macrostate.
		int level = 1;
		int partner = position;
		
		while(level != 0) {
			partner += step;
			if(partner < 0 or partner >= macrostate.length()) {
				throw (f("mismatched base-pair in '%s' macrostate: '%s'") % item.first % macrostate).str();
			}

			level += (macrostate[partner] == open_symbol);
			level -= (macrostate[partner] == close_symbol);
		}

		// Make sure the partner is mutable.
		if(not can_be_mutated(construct, partner)) {
			throw (f("position '%d' can be mutated, but it's base-paired to position '%d' which can't be.") % position % partner).str();
		}

		// Recursively mutate the partner such that every base-pairing constraint 
		// in every macrostate can be satisfied.
		char complementary_mutation = COMPLEMENTARY_NUCS.at(mutation);
		if(not already_mutated[partner]) {
			mutate_recursively(
					construct, partner, complementary_mutation, already_mutated);
		}
		else if(construct->seq()[partner] != complementary_mutation) {
			// I can't think of a way to trigger this condition (programming errors 
			// excluded), but I'm not yet convinced that there isn't a way.  If I 
			// were, I would've made this an assertion rather than an exception.
			throw "no way to satisfy all base pairing constraints.";
		}
	}
}


UnbiasedMutationMove::UnbiasedMutationMove() {}

void
UnbiasedMutationMove::apply(ConstructPtr construct, std::mt19937 &rng) const {
	// Make a list of the positions that can be mutated.
	vector<int> mutable_positions;
	for(int i = 0; i < construct->len(); i++) {
		if(can_be_freely_mutated(construct, i)) {
			mutable_positions.push_back(i);
		}
	}

	// Mutate a randomly chosen position to a randomly chosen base.
	int random_i = mutable_positions[
		std::uniform_int_distribution<>(0, mutable_positions.size()-1)(rng)];
	char random_acgu = "ACGU"[std::uniform_int_distribution<>(0, 3)(rng)];

	mutate_recursively(construct, random_i, random_acgu);
}

FixedThermostat::FixedThermostat(double temperature):
	my_temperature(temperature) {}

double
FixedThermostat::adjust(MonteCarloStep const & step) {
	return my_temperature;
}

double
FixedThermostat::temperature() const {
	return my_temperature;
}

void
FixedThermostat::temperature(double value) {
	my_temperature = value;
}


AnnealingThermostat::AnnealingThermostat(
		int cycle_len, double max_temperature, double min_temperature):

	my_cycle_len(cycle_len),
	my_max_temperature(max_temperature),
	my_min_temperature(min_temperature) {}

double
AnnealingThermostat::adjust(MonteCarloStep const &step) {
	int N = my_cycle_len;
	int i = step.i;
	double T_hi = my_max_temperature;
	double T_lo = my_min_temperature;
	return ((T_lo - T_hi) / N) * (i % N) + T_hi;
}

int
AnnealingThermostat::cycle_len() const {
	return my_cycle_len;
}

void
AnnealingThermostat::cycle_len(int value) {
	my_cycle_len = value;
}

double
AnnealingThermostat::min_temperature() const {
	return my_min_temperature;
}

void
AnnealingThermostat::min_temperature(double value) {
	my_min_temperature = value;
}

double
AnnealingThermostat::max_temperature() const {
	return my_max_temperature;
}

void
AnnealingThermostat::max_temperature(double value) {
	my_max_temperature = value;
}


AutoScalingThermostat::AutoScalingThermostat(
		double target_acceptance_rate,
		unsigned training_period,
		double initial_temperature):

	my_temperature(initial_temperature),
	my_target_acceptance_rate(target_acceptance_rate),
	my_training_period(training_period),
	my_training_set() {}

double
AutoScalingThermostat::adjust(MonteCarloStep const &step) {
	// Add the new score difference to the training set.
	my_training_set.push_back(step.score_diff);

	// Once we've acquired a certain amount of training data, calculate a new 
	// temperature by finding the median score difference and solving for the 
	// temperature that would give the target acceptance rate.
	if(my_training_set.size() >= my_training_period) {
		int n = my_training_set.size() / 2;
		nth_element(
				my_training_set.begin(),
				my_training_set.begin() + n,
				my_training_set.end());
		my_temperature = std::max(
				my_training_set[n] / log(my_target_acceptance_rate), 0.0);
		my_training_set.clear();
	}

	return my_temperature;
}

void
AutoScalingThermostat::initial_temperature(double value) {
	my_temperature = value;
}

void
AutoScalingThermostat::training_period(unsigned value) {
	my_training_period = value;
}


void
ProgressReporter::update(MonteCarloStep const &step) {
	// Print a progress bar if the program is running in a TTY.
	if(isatty(fileno(stdout))) {
		using namespace std;
		string const clear_line = "\033[2K\r";
		cout << clear_line << f("[%d/%d]") % (step.i + 1) % step.num_steps;
		if(step.i + 1 == step.num_steps) cout << endl;
		else cout << flush;
	}
}


TsvTrajectoryReporter::TsvTrajectoryReporter(string path, int interval):
	my_path(path), my_interval(interval) {}

void
TsvTrajectoryReporter::start(MonteCarloStep const & step) {
	my_tsv.open(my_path);
	if(not my_tsv.is_open()) {
		throw (f("couldn't open '%s' for writing") % my_path).str();
	}

	// Record parameters that apply to the whole trajectory.  All the lines in 
	// this section are prefixed with '#' so they won't be parsed by pandas as 
	// part of the trajectory.  This is a hack.  It'd be better to use HDF5.
	my_tsv << "#\t" << "initial_seq\t" << step.current_construct->seq() << "\n";

	// Write the column headers.
	my_tsv << "step\t";
	my_tsv << "num_steps\t";
	my_tsv << "current_score\t";
	my_tsv << "proposed_score\t";

	for(auto row: step.score_table) {
		my_tsv << f("term_weight[%s]") % row.name << "\t";
		my_tsv << f("term_value[%s]") % row.name << "\t";
	}

	my_tsv << "score_diff\t";
	my_tsv << "temperature\t";
	my_tsv << "metropolis_criterion\t";
	my_tsv << "random_threshold\t";
	my_tsv << "move\t";
	my_tsv << "outcome\t";
	my_tsv << "current_seq\t";
	my_tsv << "proposed_seq\t";
	
	my_tsv << std::endl;
}

void
TsvTrajectoryReporter::update(MonteCarloStep const &step) {
	if(step.i % my_interval == 0) {
		my_tsv << step.i << "\t";
		my_tsv << step.num_steps << "\t";
		my_tsv << step.current_score << "\t";
		my_tsv << step.proposed_score << "\t";

		for(auto row: step.score_table) {
			my_tsv << row.weight << "\t";
			my_tsv << row.term << "\t";
		}

		my_tsv << step.score_diff << "\t";
		my_tsv << step.temperature << "\t";
		my_tsv << step.metropolis_criterion << "\t";
		my_tsv << step.random_threshold << "\t";
		my_tsv << step.move->name() << "\t";
		my_tsv << step.outcome << "\t";
		my_tsv << step.current_construct->seq() << "\t";
		my_tsv << step.proposed_construct->seq() << "\t";

		my_tsv << std::endl;
	}
}

void
TsvTrajectoryReporter::finish(MonteCarloStep const &step) {
	my_tsv.close();
}


}

namespace std {

ostream &
operator<<(ostream& out, const addapt::OutcomeEnum& condition) {
	switch(condition) {
		case addapt::OutcomeEnum::REJECT: out << "REJECT"; break;
		case addapt::OutcomeEnum::ACCEPT_WORSENED: out << "ACCEPT_WORSENED"; break;
		case addapt::OutcomeEnum::ACCEPT_UNCHANGED: out << "ACCEPT_UNCHANGED"; break;
		case addapt::OutcomeEnum::ACCEPT_IMPROVED: out << "ACCEPT_IMPROVED"; break;
	}
	return out;
}

}

