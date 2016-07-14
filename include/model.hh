#pragma once

#include <iostream>
#include <iterator>
#include <memory>
#include <vector>

#include "utils.hh"

namespace addapt {

class Device;
using DevicePtr = std::shared_ptr<Device>;
using DeviceConstPtr = std::shared_ptr<Device const>;

class Aptamer;
using AptamerConstPtr = std::shared_ptr<Aptamer const>;

class Context;
using ContextConstPtr = std::shared_ptr<Context const>;

class Device {

public:

	/// @brief Default deviceor.
	Device(string);

	/// @brief Return the length of this device.
	int len() const;

	/// @brief Return the sequence of this device.
	string seq() const;

	/// @brief Return the nucleotide at the given position of this device.
	char seq(int) const;

	/// @brief Return the length of this device, neglecting the current context.
	int raw_len() const;

	/// @brief Return the sequence of this device, neglecting the current 
	/// context.
	string raw_seq() const;

	/// @brief Return the nucleotide at the given position of this device, 
	/// neglecting the current context.
	char raw_seq(int) const;

	/// @brief Return constraints that define the given macrostate.
	string macrostate(string) const;

	/// @brief Add constraints that define a particular macrostate.
	void add_macrostate(string, string);

	/// @brief Return the context that this device is operating in.
	ContextConstPtr context() const;

	/// @brief Set the context that this device is operating in.
	void context(ContextConstPtr);

	/// @brief Simulate this device by itself, without context.
	void remove_context();

	/// @brief Make a point mutation in this device.
	void mutate(int, char const);

	/// @brief Return a deep-copy of this device.
	DevicePtr copy() const;

	/// @brief Make this device equivalent to the given one.
	void assign(DevicePtr);

public:

	/// @brief
	class macrostate_iterator :
		public std::iterator<std::forward_iterator_tag, std::pair<string,string> > {

	using internal_iterator = unordered_map<string,string>::const_iterator;
	internal_iterator my_it;
	Device const &my_device;

	public:
		macrostate_iterator(Device const &dev, internal_iterator it): my_device(dev), my_it(it) {}
		macrostate_iterator(const macrostate_iterator &other): my_device(other.my_device), my_it(other.my_it) {}
		macrostate_iterator &operator++() { ++my_it; return *this; }
		macrostate_iterator operator++(int) { auto tmp(*this); operator++(); return tmp; }
		bool operator==(macrostate_iterator const &other) { return my_it == other.my_it; }
		bool operator!=(macrostate_iterator const &other) { return my_it != other.my_it; }
		value_type operator*() { return value_type(my_it->first, my_device.macrostate(my_it->first)); }
	};

	/// @brief 
	class macrostate_view {
	Device const &my_device;
	public:
		macrostate_view(Device const &dev): my_device(dev) {}
		macrostate_iterator begin() { return macrostate_iterator(my_device, my_device.my_macrostates.begin()); }
		macrostate_iterator end() { return macrostate_iterator(my_device, my_device.my_macrostates.end()); }
	};

	/// @brief Return all the macrostates associated with this device.
	macrostate_view macrostates() const;

public:

	string my_seq;
	unordered_map<string,string> my_macrostates;
	ContextConstPtr my_context;

};

class Aptamer {

public:

	/// @brief Specify the sequence, fold, and binding affinity (μM) of the 
	/// aptamer.
	Aptamer(string, string, double);

	/// @brief Return the sequence of this aptamer.
	string seq() const;

	/// @brief Return a pseudo-dot-bracket string specifying the base pairs 
	/// formed by this aptamer in the holo condition.
	string fold() const;

	/// @brief Return the affinity of this aptamer for its ligand (i.e. its  
	/// dissociation constant in units of μM).
	double affinity() const;

private:
	string my_seq;
	string my_fold;
	double my_affinity;

};

class Context {

public:

	/// @brief Specify sequences that will surround the device.
	Context(string="", string="");

	/// @brief Return the sequence that will be 5' of the device.
	string before() const;

	/// @brief Return the sequence that will be 3' of the device.
	string after() const;

private:
	string my_before;
	string my_after;

};


}
