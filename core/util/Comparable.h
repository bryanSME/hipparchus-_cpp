#pragma once


#include <complex>

// this is a simple virtual class to compare two classes
template<typename T>
class Comparable
{
public:
	virtual double compare(const std::complex<double>& o1, const std::complex<double>& o2) = 0;
	//virtual bool operator==(const Comparable& other) = 0;
};