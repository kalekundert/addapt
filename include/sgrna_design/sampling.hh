#pragma once

#include <iostream>
#include <fstream>
#include <map>
#include <memory>
#include <vector>

#include <sgrna_design/model.hh>
#include <sgrna_design/scoring.hh>
#include <sgrna_design/utils.hh>

namespace sgrna_design {

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

	/// @brief Default constructor.
	MonteCarlo();

	/// @brief Perform the Monte Carlo design simulation.
	ConstructPtr apply(ConstructPtr, std::mt19937 &) const;

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
	ConstructPtr current_sgrna, proposed_sgrna;
	MovePtr move;
	double current_score, proposed_score, score_diff;
	double temperature, metropolis_criterion, random_threshold;
	OutcomeEnum outcome;
	std::map<OutcomeEnum,int> outcome_counters;
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

	int num_cycles() const;

	void num_cycles(int);

	double max_temperature() const;

	void max_temperature(double);

	double min_temperature() const;

	void min_temperature(double);

private:

	int my_num_cycles;
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

ThermostatPtr
build_thermostat(string spec);


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

	TsvTrajectoryReporter(string);

	void start(MonteCarloStep const &);
	void update(MonteCarloStep const &);
	void finish(MonteCarloStep const &);

private:

	string my_path;
	std::ofstream my_tsv;
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



