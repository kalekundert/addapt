#pragma once

#include <memory>
#include <vector>

#include <sgrna_design/model.hh>
#include <sgrna_design/scoring.hh>

namespace sgrna_design {

class Move;
using MovePtr = std::shared_ptr<Move>;
using MoveList = std::vector<MovePtr>;

class MonteCarlo;
using MonteCarloPtr = std::shared_ptr<MonteCarlo>;

class MonteCarlo {

public:

	/// @brief Default constructor with optional beta, num_steps, and score 
	/// function arguments.
	MonteCarlo(
			double beta=1,
			int num_steps=1000,
			ScoreFunctionPtr=std::make_shared<ScoreFunction>());

	/// @brief Perform the Monte Carlo design simulation.
	ConstructPtr apply(ConstructPtr) const;

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

	private:

		double my_beta;
		int my_steps;
		ScoreFunctionPtr my_scorefxn;
		MoveList my_moves;

	};

class Move {

public:

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

	void apply(ConstructPtr, std::mt19937 &) const;

};

class MakeWtReversion : public DomainMove {

public:

	MakeWtReversion(std::vector<string>, ConstructConstPtr);

	void apply(ConstructPtr, std::mt19937 &) const;

private:

	ConstructConstPtr my_wt;
};

}




