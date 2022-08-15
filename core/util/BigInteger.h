#pragma once

#include <iostream>
#include <string>
#include <vector>


class BigInteger
{
    std::string digits;
public:

    //Constructors:
    BigInteger(unsigned long long n = 0);
    BigInteger(std::string&);
    BigInteger(const char*);
    BigInteger(BigInteger&);

    //Helper Functions:
    friend void divide_by_2(BigInteger& a);
    friend bool Null(const BigInteger&);
    friend int Length(const BigInteger&);
    int operator[](const int)const;

    /* * * * Operator Overloading * * * */

//Direct assignment
    BigInteger& operator=(const BigInteger&);

    //Post/Pre - Incrementation
    BigInteger& operator++();
    BigInteger operator++(int temp);
    BigInteger& operator--();
    BigInteger operator--(int temp);

    //Addition and Subtraction
    friend BigInteger& operator+=(BigInteger&, const BigInteger&);
    friend BigInteger operator+(const BigInteger&, const BigInteger&);
    friend BigInteger operator-(const BigInteger&, const BigInteger&);
    friend BigInteger& operator-=(BigInteger&, const BigInteger&);

    //Comparison operators
    friend bool operator==(const BigInteger&, const BigInteger&);
    friend bool operator!=(const BigInteger&, const BigInteger&);

    friend bool operator>(const BigInteger&, const BigInteger&);
    friend bool operator>=(const BigInteger&, const BigInteger&);
    friend bool operator<(const BigInteger&, const BigInteger&);
    friend bool operator<=(const BigInteger&, const BigInteger&);

    //Multiplication and Division
    friend BigInteger& operator*=(BigInteger&, const BigInteger&);
    friend BigInteger operator*(const BigInteger&, const BigInteger&);
    friend BigInteger& operator/=(BigInteger&, const BigInteger&);
    friend BigInteger operator/(const BigInteger&, const BigInteger&);

    //Modulo
    friend BigInteger operator%(const BigInteger&, const BigInteger&);
    friend BigInteger& operator%=(BigInteger&, const BigInteger&);

    //Power Function
    friend BigInteger& operator^=(BigInteger&, const BigInteger&);
    friend BigInteger operator^(BigInteger&, const BigInteger&);

    //Square Root Function
    friend BigInteger sqrt(BigInteger& a);

    //Read and Write
    friend std::ostream& operator<<(std::ostream&, const BigInteger&);
    friend std::istream& operator>>(std::istream&, BigInteger&);

    //Others
    friend BigInteger NthCatalan(int n);
    friend BigInteger NthFibonacci(int n);
    friend BigInteger Factorial(int n);
};