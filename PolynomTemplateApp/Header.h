#pragma once
#pragma warning(disable: 4996)
#include<iostream>
#include<sstream>
#include<string>
using namespace std;

template<typename T>
class Polinom {
private: int grad;
	   T* coef;
public:
	Polinom(const int = 0);
	Polinom(int, T[]);
	Polinom(int, T);
	Polinom(const Polinom<T>&);
	~Polinom();
	T daGrad() { return grad; }

	Polinom<T>& operator=(const Polinom<T>&);
	Polinom<T> operator+(const Polinom<T>&);
	Polinom<T> operator-(const Polinom<T>&);
	Polinom<T> operator*(const Polinom<T>&);
	Polinom<T> operator/(const Polinom<T>&);
	T operator==(const Polinom<T>&);
	T operator!=(const Polinom<T>&);
	Polinom<T> operator-();
	Polinom<T>& operator+=(const Polinom<T>&);
	Polinom<T>& operator-=(const Polinom<T>&);
	Polinom<T>& operator*=(const Polinom<T>&);
	Polinom<T>& operator/=(const Polinom<T>&);
	T& operator[](int);
	T operator()(const T&);
	template <typename T> friend ostream& operator<<(ostream&, const Polinom<T>&);
	template <typename T> friend istream& operator>>(istream&, Polinom<T>&);
};