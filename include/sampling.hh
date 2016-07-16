#pragma once

#include <iostream>
#include <fstream>
#include <map>
#include <memory>
#include <vector>

#include "model.hh"
#include "scoring.hh"
#include "utils.hh"

namespace addapt {

class MonteCarlo;
using MonteCarloPtr = std::shared_ptr<MonteCarlo>;

class Move;
using MovePtr = std::shared_ptr<Move>;
using MoveList = std::vector<MovePtr>;

class Thermostat;
using ThermostatPtr = std::shared_ptr<Thermostat>;

class Reporter;
using ReporterPtr = std::shared_ptr<Reporter>;
using ReporterList = std::vector<ReporterPtr>;

class MonteCarlo {

public:

	/// @brief Default deviceor.
	MonteCarlo();

	/// @brief Perform the Monte Carlo design simulation.
	DevicePtr apply(DevicePtr, std::mt19937 &) const;

	/// @brief Return the number of moves that will be tried during the 
	/// simulation.
	int num_steps() const;

	/// @brief Set the number of moves that will be tried during the simulation.
	void num_steps(int);

	/// @brief Return the object responsible for setting the "temperature" of the 
	/// Metropolis criterion.
	ThermostatPtr thermostat() const;

	/// @brief Set the object responsible for setting the "temperature" of the 
	/// Metropolis criterion.
	void thermostat(ThermostatPtr);

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

		int my_steps;
		ThermostatPtr my_thermostat;
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
	DevicePtr current_device, proposed_device;
	MovePtr move;
	EvaluatedScoreFunction score_table;
	double current_score, proposed_score, score_diff;
	double temperature, metropolis_criterion, random_threshold;
	OutcomeEnum outcome;
	std::map<OutcomeEnum,int> outcome_counters;
};


map<char, char> const COMPLEMENTARY_NUCS = {
	{'A','U'},{'G','C'},{'C','G'},{'U','A'}};

bool
can_be_mutated(DeviceConstPtr, int);

bool
can_be_freely_mutated(DeviceConstPtr, int);

void
mutate_recursively(DevicePtr, int const, char const);

void
mutate_recursively(DevicePtr, int const, char const, vector<bool> &);


class Move {

public:

	virtual string name() const = 0;

	virtual void apply(DevicePtr, std::mt19937 &) const = 0;

};

class UnbiasedMutationMove : public Move {

public:

	UnbiasedMutationMove();

	string name() const { return "UnbiasedMutation"; }

	void apply(DevicePtr, std::mt19937 &) const;

};


class Thermostat {

public:

	virtual double adjust(MonteCarloStep const &) = 0;

};

class FixedThermostat : public Thermostat {

public:

	FixedThermostat(double);

	double adjust(MonteCarloStep const &);

	double temperature() const;

	void temperature(double);

private:

	double my_temperature;

};

class AnnealingThermostat : public Thermostat {

public:

	AnnealingThermostat(int, double, double);

	double adjust(MonteCarloStep const &);

	int cycle_len() const;

	void cycle_len(int);

	double max_temperature() const;

	void max_temperature(double);

	double min_temperature() const;

	void min_temperature(double);

private:

	int my_cycle_len;
	double my_max_temperature;
	double my_min_temperature;

};

class AutoScalingThermostat : public Thermostat {

public:

	AutoScalingThermostat(double=0.5, unsigned=100, double=1.0);

	double adjust(MonteCarloStep const &);

	void initial_temperature(double);

	void target_acceptance_rate(double);

	void training_period(unsigned);

private:

	double my_temperature;
	double my_target_acceptance_rate;
	unsigned my_training_period;
	vector<double> my_training_set;

};

class AutoAnnealingThermostat : public Thermostat {

public:

	AutoAnnealingThermostat(int, double, double, double=1.0, double=0.0);

	double adjust(MonteCarloStep const &);

private:

	int my_cycle_len;
	double my_high_acceptance_rate, my_low_acceptance_rate;
	double my_initial_high_temperature, my_initial_low_temperature;
	double my_sum_score_diffs;
	int my_num_score_diffs;
	double my_high_temperature, my_low_temperature;

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
class TsvTrajectoryReporter : public Reporter {

public:

	TsvTrajectoryReporter(string, int);

	void start(MonteCarloStep const &);
	void update(MonteCarloStep const &);
	void finish(MonteCarloStep const &);

private:

	string my_path;
	int my_interval;
	std::ofstream my_tsv;
};


}

namespace std {

ostream&
operator<<(ostream&, const addapt::OutcomeEnum&);

}



