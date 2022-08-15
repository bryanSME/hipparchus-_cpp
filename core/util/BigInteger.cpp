#pragma once

#include "BigInteger.h"
#include <iostream>
#include <string>
#include <vector>

BigInteger::BigInteger(std::string& s)
{
    digits = "";
    int n = s.size();
    for (int i = n - 1; i >= 0; i--)
    {
        if (!isdigit(s[i]))
        {
            throw("ERROR");
        }
        digits.push_back(s[i] - '0');
    }
}

BigInteger::BigInteger(unsigned long long nr)
{
    do
    {
        digits.push_back(nr % 10);
        nr /= 10;
    } while (nr);
}

BigInteger::BigInteger(const char* s)
{
    digits = "";
    for (int i = strlen(s) - 1; i >= 0; i--)
    {
        if (!isdigit(s[i]))
        {
            throw("ERROR");
        }
        digits.push_back(s[i] - '0');
    }
}

BigInteger::BigInteger(BigInteger& a)
{
    digits = a.digits;
}

bool Null(const BigInteger& a)
{
    return (a.digits.size() == 1 && a.digits[0] == 0);
}

int Length(const BigInteger& a)
{
    return a.digits.size();
}

int BigInteger::operator[](const int index) const
{
    if (digits.size() <= index || index < 0)
    {
        throw("ERROR");
    }
    return digits[index];
}

bool operator==(const BigInteger& a, const BigInteger& b)
{
    return a.digits == b.digits;
}

bool operator!=(const BigInteger& a, const BigInteger& b) {
    return !(a == b);
}

bool operator<(const BigInteger& a, const BigInteger& b)
{
    int n = Length(a), m = Length(b);
    if (n != m)
    {
        return n < m;
    }
    while (n--)
    {
        if (a.digits[n] != b.digits[n])
        {
            return a.digits[n] < b.digits[n];
        }
    }
    return false;
}

bool operator>(const BigInteger& a, const BigInteger& b)
{
    return b < a;
}

bool operator>=(const BigInteger& a, const BigInteger& b)
{
    return !(a < b);
}

bool operator<=(const BigInteger& a, const BigInteger& b)
{
    return !(a > b);
}

BigInteger& BigInteger::operator=(const BigInteger& a)
{
    digits = a.digits;
    return *this;
}

BigInteger& BigInteger::operator++()
{
    int i;
    int n = digits.size();
    for (i = 0; i < n && digits[i] == 9; i++)
    {
        digits[i] = 0;
    }
    if (i == n)
    {
        digits.push_back(1);
    }
    else
    {
        digits[i]++;
    }
    return *this;
}

BigInteger BigInteger::operator++(int temp)
{
    BigInteger aux;
    aux = *this;
    ++(*this);
    return aux;
}

BigInteger& BigInteger::operator--()
{
    if (digits[0] == 0 && digits.size() == 1)
    {
        throw("UNDERFLOW");
    }
    int i;
    int n = digits.size();
    for (i = 0; digits[i] == 0 && i < n; i++)
        digits[i] = 9;
    digits[i]--;
    if (n > 1 && digits[n - 1] == 0)
    {
        digits.pop_back();
    }
    return *this;
}

BigInteger BigInteger::operator--(int temp)
{
    BigInteger aux;
    aux = *this;
    --(*this);
    return aux;
}

BigInteger& operator+=(BigInteger& a, const BigInteger& b)
{
    int t{};
    int s;
    int i;
    int n = Length(a);
    int m = Length(b);
    if (m > n)
    {
        a.digits.append(m - n, 0);
    }
    n = Length(a);
    for (i = 0; i < n; i++)
    {
        if (i < m)
        {
            s = (a.digits[i] + b.digits[i]) + t;
        }
        else
        {
            s = a.digits[i] + t;
        }
        t = s / 10;
        a.digits[i] = (s % 10);
    }
    if (t)
    {
        a.digits.push_back(t);
    }
    return a;
}

BigInteger operator+(const BigInteger& a, const BigInteger& b)
{
    BigInteger temp;
    temp = a;
    temp += b;
    return temp;
}

BigInteger& operator-=(BigInteger& a, const BigInteger& b)
{
    if (a < b)
    {
        throw("UNDERFLOW");
    }
    int n = Length(a);
    int m = Length(b);
    int i;
    int t{};
    int s;
    for (i = 0; i < n; i++)
    {
        if (i < m)
        {
            s = a.digits[i] - b.digits[i] + t;
        }
        else
        {
            s = a.digits[i] + t;
        }
        if (s < 0)
        {
            s += 10;
            t = -1;
        }
        else
        {
            t = 0;
        }
        a.digits[i] = s;
    }
    while (n > 1 && a.digits[n - 1] == 0)
    {
        a.digits.pop_back();
        n--;
    }
    return a;
}

BigInteger operator-(const BigInteger& a, const BigInteger& b)
{
    BigInteger temp;
    temp = a;
    temp -= b;
    return temp;
}

BigInteger& operator*=(BigInteger& a, const BigInteger& b)
{
    if (Null(a) || Null(b))
    {
        a = BigInteger();
        return a;
    }
    int n = a.digits.size();
    int m = b.digits.size();
    std::vector<int> v(n + m, 0);
    for (int i{}; i < n; i++)
    {
        for (int j{}; j < m; j++)
        {
            v[i + j] += (a.digits[i]) * (b.digits[j]);
        }
    }
    n += m;
    a.digits.resize(v.size());
    for (int s, i = 0, t = 0; i < n; i++)
    {
        s = t + v[i];
        v[i] = s % 10;
        t = s / 10;
        a.digits[i] = v[i];
    }
    for (int i = n - 1; i >= 1 && !v[i]; i--)
    {
        a.digits.pop_back();
    }
    return a;
}

BigInteger operator*(const BigInteger& a, const BigInteger& b)
{
    BigInteger temp;
    temp = a;
    temp *= b;
    return temp;
}

BigInteger& operator/=(BigInteger& a, const BigInteger& b)
{
    if (Null(b))
    {
        throw("Arithmetic Error: Division By 0");
    }
    if (a < b)
    {
        a = BigInteger();
        return a;
    }
    if (a == b)
    {
        a = BigInteger(1);
        return a;
    }
    int i;
    int lgcat{};
    int cc;
    int n = Length(a);
    int m = Length(b);
    std::vector<int> cat(n, 0);
    BigInteger t;
    for (i = n - 1; t * 10 + a.digits[i] < b; i--)
    {
        t *= 10;
        t += a.digits[i];
    }
    for (; i >= 0; i--)
    {
        t = t * 10 + a.digits[i];
        for (cc = 9; cc * b > t; cc--);
        t -= cc * b;
        cat[lgcat++] = cc;
    }
    a.digits.resize(cat.size());
    for (i = 0; i < lgcat; i++)
    {
        a.digits[i] = cat[lgcat - i - 1];
    }
    a.digits.resize(lgcat);
    return a;
}

BigInteger operator/(const BigInteger& a, const BigInteger& b)
{
    BigInteger temp;
    temp = a;
    temp /= b;
    return temp;
}

BigInteger& operator%=(BigInteger& a, const BigInteger& b)
{
    if (Null(b))
    {
        throw("Arithmetic Error: Division By 0");
    }
    if (a < b)
    {
        a = BigInteger();
        return a;
    }
    if (a == b)
    {
        a = BigInteger(1);
        return a;
    }
    int i, lgcat = 0;
    int cc;
    int n = Length(a);
    int m = Length(b);
    std::vector<int> cat(n, 0);
    BigInteger t;
    for (i = n - 1; t * 10 + a.digits[i] < b; i--)
    {
        t *= 10;
        t += a.digits[i];
    }
    for (; i >= 0; i--)
    {
        t = t * 10 + a.digits[i];
        for (cc = 9; cc * b > t; cc--);
        t -= cc * b;
        cat[lgcat++] = cc;
    }
    a = t;
    return a;
}

BigInteger operator%(const BigInteger& a, BigInteger& b)
{
    BigInteger temp;
    temp = a;
    temp %= b;
    return temp;
}

BigInteger& operator^=(BigInteger& a, const BigInteger& b)
{
    BigInteger Exponent, Base(a);
    Exponent = b;
    a = 1;
    while (!Null(Exponent))
    {
        if (Exponent[0] & 1)
        {
            a *= Base;
        }
        Base *= Base;
        divide_by_2(Exponent);
    }
    return a;
}

BigInteger operator^(BigInteger& a, BigInteger& b)
{
    BigInteger temp(a);
    temp ^= b;
    return temp;
}

void divide_by_2(BigInteger& a)
{
    int add{};
    for (int i = a.digits.size() - 1; i >= 0; i--)
    {
        int digit = (a.digits[i] >> 1) + add;
        add = ((a.digits[i] & 1) * 5);
        a.digits[i] = digit;
    }
    while (a.digits.size() > 1 && !a.digits.back())
    {
        a.digits.pop_back();
    }
}

BigInteger sqrt(BigInteger& a)
{
    BigInteger left(1);
    BigInteger right(a);
    BigInteger v(1);
    BigInteger mid;
    BigInteger prod;
    divide_by_2(right);
    while (left <= right)
    {
        mid += left;
        mid += right;
        divide_by_2(mid);
        prod = (mid * mid);
        if (prod <= a)
        {
            v = mid;
            ++mid;
            left = mid;
        }
        else
        {
            --mid;
            right = mid;
        }
        mid = BigInteger();
    }
    return v;
}

BigInteger NthCatalan(int n)
{
    BigInteger a(1), b;
    for (int i = 2; i <= n; i++)
    {
        a *= i;
    }
    b = a;
    for (int i = n + 1; i <= 2 * n; i++)
    {
        b *= i;
    }
    a *= a;
    a *= (n + 1);
    b /= a;
    return b;
}

BigInteger NthFibonacci(int n)
{
    BigInteger a(1), b(1), c;
    if (!n)
        return c;
    n--;
    while (n--)
    {
        c = a + b;
        b = a;
        a = c;
    }
    return b;
}

BigInteger Factorial(int n)
{
    BigInteger f(1);
    for (int i = 2; i <= n; i++)
    {
        f *= i;
    }
    return f;
}

std::istream& operator>>(std::istream& in, BigInteger& a)
{
    std::string s;
    in >> s;
    int n = s.size();
    for (int i = n - 1; i >= 0; i--)
    {
        if (!isdigit(s[i]))
        {
            throw("INVALID NUMBER");
        }
        a.digits[n - i - 1] = s[i];
    }
    return in;
}

std::ostream& operator<<(std::ostream& out, const BigInteger& a)
{
    for (int i = a.digits.size() - 1; i >= 0; i--)
    {
        std::cout << (short)a.digits[i];
    }
    return std::cout;
}

//Driver code with some examples
int main()
{
    BigInteger first("12345");
    std::cout
        << "The number of digits"
        << " in first big integer = "
        << Length(first) << '\n';
    BigInteger second(12345);
    if (first == second)
    {
        std::cout << "first and second are equal!\n";
    }
    else
    {
        std::cout << "Not equal!\n";
    }
    BigInteger third("10000");
    BigInteger fourth("100000");
    if (third < fourth)
    {
        std::cout << "third is smaller than fourth!\n";
    }
    BigInteger fifth("10000000");
    if (fifth > fourth)
    {
        std::cout << "fifth is larger than fourth!\n";
    }

    // Printing all the numbers
    std::cout << "first = " << first << '\n';
    std::cout << "second = " << second << '\n';
    std::cout << "third = " << third << '\n';
    std::cout << "fourth = " << fourth << '\n';
    std::cout << "fifth = " << fifth << '\n';

    // Incrementing the value of first
    first++;
    std::cout
        << "After incrementing the"
        << " value of first is : ";
    std::cout << first << '\n';
    BigInteger sum;
    sum = (fourth + fifth);
    std::cout
        << "Sum of fourth and fifth = "
        << sum << '\n';
    BigInteger product;
    product = second * third;
    std::cout
        << "Product of second and third = "
        << product << '\n';

    // Print the fibonaccii number from 1 to 100
    std::cout
        << "-------------------------Fibonacci"
        << "------------------------------\n";
    for (int i{}; i <= 100; i++)
    {
        BigInteger Fib;
        Fib = NthFibonacci(i);
        std::cout << "Fibonacci " << i << " = " << Fib << '\n';
    }
    std::cout
        << "-------------------------Catalan"
        << "------------------------------\n";
    for (int i{}; i <= 100; i++)
    {
        BigInteger Cat;
        Cat = NthCatalan(i);
        std::cout << "Catalan " << i << " = " << Cat << '\n';
    }

    // Calculating factorial of from 1 to 100
    std::cout
        << "-------------------------Factorial"
        << "------------------------------\n";
    for (int i = 0; i <= 100; i++)
    {
        BigInteger fact;
        fact = Factorial(i);
        std::cout
            << "Factorial of "
            << i << " = ";
        std::cout << fact << '\n';
    }
    // This code is contributed
    // by Gatea David
}