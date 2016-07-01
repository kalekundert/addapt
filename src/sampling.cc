#include <algorithm>
#include <cmath>
#include <iostream>

#include <sgrna_design/sampling.hh>

namespace sgrna_design {

MonteCarlo::MonteCarlo(double beta, int num_steps, ScoreFunctionPtr scorefxn):
	my_beta(beta),
	my_steps(num_steps),
	my_scorefxn(scorefxn),
	my_reporters({std::make_shared<ProgressReporter>()}) {}

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
	step.beta = my_beta;

	// Create random number distributions for later use.
	auto random = std::bind(std::uniform_real_distribution<>(), rng);
	auto randmove = std::bind(std::uniform_int_distribution<>(0, my_moves.size()-1), rng);

	// Get an initial score.
	step.current_score = my_scorefxn->evaluate(sgrna);
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
			step.proposed_score = my_scorefxn->evaluate(step.proposed_sgrna);
			step.score_diff = step.proposed_score - step.current_score;
			step.metropolis_criterion = std::exp(step.beta * step.score_diff);
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

double
MonteCarlo::beta() const {
	return my_beta;
}

void
MonteCarlo::beta(double beta) {
	my_beta = beta;
}

int
MonteCarlo::num_steps() const {
	return my_steps;
}

void
MonteCarlo::num_steps(int num_steps) {
	my_steps = num_steps;
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


CsvTrajectoryReporter::CsvTrajectoryReporter(string path):
	my_path(path) {}

void
CsvTrajectoryReporter::start(MonteCarloStep const & step) {
	my_csv.open(my_path);

	my_csv << "step\t";
	my_csv << "num_steps\t";
	my_csv << "current_score\t";
	my_csv << "proposed_score\t";
	my_csv << "score_diff\t";
	my_csv << "beta\t";
	my_csv << "metropolis_criterion\t";
	my_csv << "random_threshold\t";
	my_csv << "move\t";
	my_csv << "outcome\t";

	for(auto domain: step.current_sgrna->domains()) {
		my_csv << "current_" << domain->name() << "\t";
	}
	for(auto domain: step.current_sgrna->domains()) {
		my_csv << "proposed_" << domain->name() << "\t";
	}
	
	my_csv << std::endl;
}

void
CsvTrajectoryReporter::update(MonteCarloStep const &step) {
	my_csv << step.i << "\t";
	my_csv << step.num_steps << "\t";
	my_csv << step.current_score << "\t";
	my_csv << step.proposed_score << "\t";
	my_csv << step.score_diff << "\t";
	my_csv << step.beta << "\t";
	my_csv << step.metropolis_criterion << "\t";
	my_csv << step.random_threshold << "\t";
	my_csv << step.move->name() << "\t";
	my_csv << step.outcome << "\t";

	for(auto domain: step.current_sgrna->domains()) {
		my_csv << domain->seq() << "\t";
	}
	for(auto domain: step.proposed_sgrna->domains()) {
		my_csv << domain->seq() << "\t";
	}

	my_csv << std::endl;
}

void
CsvTrajectoryReporter::finish(MonteCarloStep const &step) {
	my_csv.close();
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

