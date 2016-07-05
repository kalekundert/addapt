#include <algorithm>
#include <cmath>
#include <iostream>
#include <regex>

#include <sgrna_design/sampling.hh>
#include <sgrna_design/utils.hh>

namespace sgrna_design {

MonteCarlo::MonteCarlo(): 
	my_steps(0),
	my_thermostat(std::make_shared<FixedThermostat>(1)),
	my_scorefxn(std::make_shared<ScoreFunction>()),
	my_moves(),
	my_reporters() {}

ConstructPtr
MonteCarlo::apply(ConstructPtr sgrna, std::mt19937 &rng) const {
	if (my_moves.empty()) {
		return sgrna;
	}

	// Setup to data structure that will hold all the information about each 
	// step.  The purpose of this structure is to support external logging 
	// methods and thermostats.
	MonteCarloStep step;
	step.num_steps = my_steps; step.i = -1;
	step.current_sgrna = sgrna;

	// Create random number distributions for later use.
	auto random = std::bind(std::uniform_real_distribution<>(), rng);
	auto randmove = std::bind(std::uniform_int_distribution<>(0, my_moves.size()-1), rng);

	// Get an initial score.
	step.current_score = my_scorefxn->evaluate(step.current_sgrna, step.score_table);
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
		step.proposed_sgrna = step.current_sgrna->copy();

		// Randomly pick a move to apply.
		step.move = my_moves[randmove()];
		step.move->apply(step.proposed_sgrna, rng);

		// Skip the score function evaluation if the sequence didn't change.
		if(step.current_sgrna->seq() == step.proposed_sgrna->seq()) {
			step.outcome = OutcomeEnum::ACCEPT_UNCHANGED;
		}

		// Score the proposed move, then either accept or reject it.
		else {
			step.proposed_score = my_scorefxn->evaluate(step.proposed_sgrna, step.score_table);
			step.score_diff = step.proposed_score - step.current_score;
			step.metropolis_criterion = std::exp(step.score_diff / step.temperature);
			step.random_threshold = random();

			if(step.metropolis_criterion < step.random_threshold) {
				step.outcome = OutcomeEnum::REJECT;
			}
			else{
				step.outcome = (step.score_diff > 0)?
					OutcomeEnum::ACCEPT_IMPROVED : OutcomeEnum::ACCEPT_WORSENED;

				step.current_sgrna = step.proposed_sgrna;
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

	return step.current_sgrna;
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


ThermostatPtr
build_thermostat(string spec) {
	std::regex fixed_pattern(
			"([0-9.e+-]+)" 	    // A floating-point number (the temperature).
	);
	std::regex annealing_pattern(
			"(?:"               // Optional argument.
			"([0-9]+)x"         // An integer followed by 'x' (the number of cycles).
			"\\s+"
			")?"
			"([0-9.e+-]+)"      // A floating point number (the high temperature).
			"\\s*"
			"=>"                // An arrow.
			"\\s*"
			"([0-9.e+-]+)"      // A floating point number (the low temperature).
	);
	std::regex auto_scaling_pattern(
			"auto"
			"(?:"						    // Optional argument.
			"\\s+"
			"([0-9.]+)%"        // A percentage (the target acceptance rate).
				"(?:"         
				"\\s+"
				"([0-9]+)"        // An integer (the training period).
					"(?:"
					"\\s+"
					"([0-9.e+-]+)"  // A floating point number (the initial temperature).
					")?"
				")?"
			")?"
	);

	std::smatch match;

	if(std::regex_match(spec, match, fixed_pattern)) {
		double temperature = stod(match[1]);
		return make_shared<FixedThermostat>(temperature);
	}

	if(std::regex_match(spec, match, annealing_pattern)) {
		int num_cycles = stoi(match[1].length()? match[1].str() : "1");
		double high_temp = stod(match[2]);
		double low_temp = stod(match[3]);
		return make_shared<AnnealingThermostat>(num_cycles, high_temp, low_temp);
	}

	if(std::regex_match(spec, match, auto_scaling_pattern)) {
		double accept_rate = stod(match[1].length()? match[1].str() : "50") / 100;
		int training_period = stod(match[2].length()? match[2].str() : "100");
		double initial_temp = stod(match[3].length()? match[3].str() : "1");
		return make_shared<AutoScalingThermostat>(
				accept_rate, training_period, initial_temp);
	}

	throw (f("can't make a thermostat from '%s'") % spec).str();
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
		int num_cycles, double max_temperature, double min_temperature):

	my_num_cycles(num_cycles),
	my_max_temperature(max_temperature),
	my_min_temperature(min_temperature) {}

double
AnnealingThermostat::adjust(MonteCarloStep const &step) {
	int N = step.num_steps / my_num_cycles;
	int i = step.i;
	double T_hi = my_max_temperature;
	double T_lo = my_min_temperature;
	return ((T_lo - T_hi) / N) * (i % N) + T_hi;
}

int
AnnealingThermostat::num_cycles() const {
	return my_num_cycles;
}

void
AnnealingThermostat::num_cycles(int value) {
	my_num_cycles = value;
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


TsvTrajectoryReporter::TsvTrajectoryReporter(string path):
	my_path(path) {}

void
TsvTrajectoryReporter::start(MonteCarloStep const & step) {
	my_tsv.open(my_path);

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

	for(auto domain: step.current_sgrna->domains()) {
		my_tsv << f("current_domain[%s]") % domain->name() << "\t";
	}
	for(auto domain: step.current_sgrna->domains()) {
		my_tsv << f("proposed_domain[%s]") % domain->name() << "\t";
	}
	
	my_tsv << std::endl;
}

void
TsvTrajectoryReporter::update(MonteCarloStep const &step) {
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

	for(auto domain: step.current_sgrna->domains()) {
		my_tsv << domain->seq() << "\t";
	}
	for(auto domain: step.proposed_sgrna->domains()) {
		my_tsv << domain->seq() << "\t";
	}

	my_tsv << std::endl;
}

void
TsvTrajectoryReporter::finish(MonteCarloStep const &step) {
	my_tsv.close();
}


DomainMove::DomainMove(std::vector<string> domain_names):
	my_domain_names(domain_names) {}

std::vector<string>
DomainMove::domain_names() const {
	return my_domain_names;
}

string
DomainMove::random_domain_name(std::mt19937 &rng) const {
	auto dist = std::uniform_int_distribution<>(0, my_domain_names.size()-1);
	return my_domain_names[dist(rng)];
}


MakePointMutation::MakePointMutation(std::vector<string> domain_names):
	DomainMove(domain_names) {}

void
MakePointMutation::apply(ConstructPtr sgrna, std::mt19937 &rng) const {
	DomainPtr random_domain = (*sgrna)[random_domain_name(rng)];
	int random_i = std::uniform_int_distribution<>(0, random_domain->len()-1)(rng);
	char random_acgu = "ACGU"[std::uniform_int_distribution<>(0, 3)(rng)];
	random_domain->mutate(random_i, random_acgu);
}


MakeWtReversion::MakeWtReversion(
		std::vector<string> domain_names,
		ConstructConstPtr wt):

	DomainMove(domain_names),
  my_wt(wt) {}

void
MakeWtReversion::apply(ConstructPtr sgrna, std::mt19937 &rng) const {
	DomainPtr random_domain = (*sgrna)[random_domain_name(rng)];
	int random_i = std::uniform_int_distribution<>(0, random_domain->len()-1)(rng);
	char wt_acgu = (*my_wt)[random_domain]->seq(random_i);
	random_domain->mutate(random_i, wt_acgu);
}


ChangeDomainLength::ChangeDomainLength(
		string domain_name, int min_length, int max_length):

	my_domain_name(domain_name),
	my_max_len(max_length),
	my_min_len(min_length) {}


void
ChangeDomainLength::apply(ConstructPtr sgrna, std::mt19937 &rng) const {
	DomainPtr domain = (*sgrna)[my_domain_name];

	// Either try to increment or decrement the length of the domain by one, 
	// being careful to keep the length of domain within the prescribed limits.  
	// Note that the rest of the function is written generally enough to 
	// accommodate any algorithm for picking random lengths.

	int random_length = std::min(std::max(
			domain->len() + (std::bernoulli_distribution()(rng) ? 1 : -1),
			my_min_len), my_max_len);

	// Insert or remove nucleotides as necessary to bring the domain to the 
	// randomly chosen length.

	while(domain->len() < random_length) {
		int random_i = std::uniform_int_distribution<>(0, domain->len()-1)(rng);
		char random_acgu = "ACGU"[std::uniform_int_distribution<>(0, 3)(rng)];
		domain->insert(random_i, random_acgu);
	}

	while(domain->len() > random_length) {
		int random_i = std::uniform_int_distribution<>(0, domain->len()-1)(rng);
		domain->remove(random_i);
	}
}


}

namespace std {

ostream &
operator<<(ostream& out, const sgrna_design::OutcomeEnum& condition) {
	switch(condition) {
		case sgrna_design::OutcomeEnum::REJECT: out << "REJECT"; break;
		case sgrna_design::OutcomeEnum::ACCEPT_WORSENED: out << "ACCEPT_WORSENED"; break;
		case sgrna_design::OutcomeEnum::ACCEPT_UNCHANGED: out << "ACCEPT_UNCHANGED"; break;
		case sgrna_design::OutcomeEnum::ACCEPT_IMPROVED: out << "ACCEPT_IMPROVED"; break;
	}
	return out;
}

}

