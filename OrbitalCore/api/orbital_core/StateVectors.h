#ifndef STATEVECTORS_H
#define STATEVECTORS_H
#include "orbital_core/Vec3D.h"
#include <iostream>
#include <utility>

class StateVectors : public std::pair<Vec3D, Vec3D> {
   public:
	StateVectors() : std::pair<Vec3D, Vec3D>() {}

	StateVectors(const Vec3D& first, const Vec3D& second) : std::pair<Vec3D, Vec3D>(first, second) {}

	StateVectors(const StateVectors& other) : std::pair<Vec3D, Vec3D>(other) {}

	StateVectors(StateVectors&& other) : std::pair<Vec3D, Vec3D>(std::move(other)) {}

	StateVectors& operator=(const StateVectors& other)
	{
		first = other.first;
		second = other.second;
		return *this;
	}

	StateVectors& operator=(StateVectors&& other)
	{
		first = std::move(other.first);
		second = std::move(other.second);
		return *this;
	}

	~StateVectors() = default;

	friend std::ostream& operator<<(std::ostream& os, const StateVectors& sv)
	{
		os << "R: " << sv.first << std::endl << "V: " << sv.second;
		return os;
	}

	bool isValid() const { return first.magnitude() != 0 && second.magnitude() != 0; }

	void print(std::ostream& os = std::cout) const { os << *this; }
};

#endif  // STATEVECTORS_H
