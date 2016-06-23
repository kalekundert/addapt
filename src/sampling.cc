#include <cmath>

#include "sampling.hh"

namespace sgrna {

MonteCarlo::MonteCarlo(double beta, int num_steps, ScoreFunctionPtr scorefxn):
	my_beta(beta),
	my_steps(num_steps),
	my_scorefxn(scorefxn) {}

ConstructPtr
MonteCarlo::apply(ConstructPtr sgrna) const {
	if (my_moves.empty()) {
		return sgrna;
	}

	// Create random number generators that use the Mersenne twister algorithm 
	// with a state size of 19937 bits.
	std::mt19937 rng;
	auto random = std::bind(std::uniform_real_distribution<>(), rng);
	auto randmove = std::bind(std::uniform_int_distribution<>(0, my_moves.size()), rng);

	// Get an initial score.
	double score = my_scorefxn->evaluate(sgrna);

	for(int i = 0; i <= my_steps; i++) {
		// Copy the sgRNA so we can easily undo the move.
		ConstructPtr trial_sgrna = sgrna->copy();

		// Randomly pick a move to apply.
		MovePtr move = my_moves[randmove()];
		move->apply(trial_sgrna, rng);

		// Decide whether or not to accept the move.
		double trial_score = my_scorefxn->evaluate(trial_sgrna);
		double score_diff = score - trial_score;
		double metropolis_criterion = std::exp(-my_beta * score_diff);
		double random_threshold = random();

		if(metropolis_criterion > random_threshold) {
			sgrna = trial_sgrna;
			score = trial_score;
		}
	}

	return sgrna;
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
operator+=(MonteCarloPtr monte_carlo, MovePtr move) {
	monte_carlo->add_move(move);
}


/*
void
MakePointMutation::apply(ConstructPtr sgrna, std::mt19937 rng) {
	// Pick domain
	DomainPtr random_domain = my_domains[
			std::uniform_int_distribution<>(0, my_domains.size())(rng)];
	int random_i = std::uniform_int_distribution(0, random_domain.len())(rng);
	char random_acgu = "ACGU"[std::uniform_int_distribution(0, 4)(rng)];

	random_domain.mutate(random_i, random_acgu);
}

MakeWtReversion::MakeWtReversion(ConstructPtr wt) : my_wt(wt) {}

void
MakeWtReversion::apply(ConstructPtr agrna, std:mt19937 rng) {
}
*/










}
