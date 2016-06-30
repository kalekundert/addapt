#pragma once

#include <iostream>
#include <fstream>
#include <map>
#include <memory>
#include <vector>

#include <sgrna_design/model.hh>
#include <sgrna_design/scoring.hh>

namespace sgrna_design {

class MonteCarlo;
using MonteCarloPtr = std::shared_ptr<MonteCarlo>;

class Move;
using MovePtr = std::shared_ptr<Move>;
using MoveList = std::vector<MovePtr>;

class Reporter;
using ReporterPtr = std::shared_ptr<Reporter>;
using ReporterList = std::vector<ReporterPtr>;

class MonteCarlo {

public:

	/// @brief Default constructor with optional beta, num_steps, and score 
	/// function arguments.
	MonteCarlo(
			double beta=1,
			int num_steps=1000,
			ScoreFunctionPtr=std::make_shared<ScoreFunction>());

	/// @brief Perform the Monte Carlo design simulation.
	ConstructPtr apply(ConstructPtr, std::mt19937 &) const;

	/// @brief Return the factor that controls how often moves are accepted.
	double beta() const;

	/// @brief Set the factor that controls how often moves are accepted.
	void beta(double);

	/// @brief Return the number of moves that will be tried during the 
	/// simulation.
	int num_steps() const;

	/// @brief Set the number of moves that will be tried during the simulation.
	void num_steps(int);

	/// @brief Return the score function.
	ScoreFunctionPtr scorefxn() const;

	/// @brief Set the score function.
	void scorefxn(ScoreFunctionPtr);

	/// @brief Return the list of possible moves.
	MoveList moves() const;

	/// @brief Add a move.
	void add_move(MovePtr);

	/// @brief Add a move.
	void operator+=(MovePtr);

	/// @brief Return the reporters being used in this simulation.
	ReporterList reporters() const;

	/// @brief Add a reporter.
	void add_reporter(ReporterPtr);

	/// @brief Add a reporter.
	void operator+=(ReporterPtr);

	private:

		double my_beta;
		int my_steps;
		ScoreFunctionPtr my_scorefxn;
		MoveList my_moves;
		ReporterList my_reporters;

	};

enum class OutcomeEnum {
	REJECT,
	ACCEPT_WORSENED,
	ACCEPT_UNCHANGED,
	ACCEPT_IMPROVED,
};

struct MonteCarloStep {
	int i, num_steps;
	ConstructPtr current_sgrna, proposed_sgrna;
	MovePtr move;
	double current_score, proposed_score, score_diff;
	double beta, metropolis_criterion, random_threshold;
	OutcomeEnum outcome;
	std::map<OutcomeEnum,int> outcome_counters;
};


class Reporter {

public:

	virtual void start(MonteCarloStep const &) {};
	virtual void update(MonteCarloStep const &) {};
	virtual void finish(MonteCarloStep const &) {};

};

class ProgressReporter : public Reporter {

public:

	void update(MonteCarloStep const &);

};

/// @details Columns will become misaligned if domains are added or removed 
/// during the simulation.
class CsvTrajectoryReporter : public Reporter {

public:

	CsvTrajectoryReporter(string);

	void start(MonteCarloStep const &);
	void update(MonteCarloStep const &);
	void finish(MonteCarloStep const &);

private:

	string my_path;
	std::ofstream my_csv;
};


class Thermostat {

public:

	virtual double beta(MonteCarloStep const &) = 0;

};


class Move {

public:

	virtual string name() const = 0;

	virtual void apply(ConstructPtr, std::mt19937 &) const = 0;

};

class DomainMove : public Move {

public:

	DomainMove(std::vector<string>);

	std::vector<string> domain_names() const;

protected:

	string random_domain_name(std::mt19937 &) const;

private:

	std::vector<string> my_domain_names;

};

class MakePointMutation : public DomainMove {

public:

	MakePointMutation(std::vector<string>);

	string name() const { return "MakePointMutation"; }

	void apply(ConstructPtr, std::mt19937 &) const;

};

class MakeWtReversion : public DomainMove {

public:

	MakeWtReversion(std::vector<string>, ConstructConstPtr);

	string name() const { return "MakeWtReversion"; }

	void apply(ConstructPtr, std::mt19937 &) const;

private:

	ConstructConstPtr my_wt;
};

class ChangeDomainLength : public Move {

public:

	ChangeDomainLength(string, int, int);

	string name() const { return "ChangeDomainLength"; }

	void apply(ConstructPtr, std::mt19937 &) const;

private:

	string my_domain_name;
	int my_max_len, my_min_len;
};

}

namespace std {

ostream&
operator<<(ostream&, const sgrna_design::OutcomeEnum&);

}



