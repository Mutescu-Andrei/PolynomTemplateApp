#include"Header.h"
template <typename T>

Polinom<T>::Polinom(const int i) {
	grad = 0;
	coef = new T[grad + 1];
	coef[0] = i;
}

template <typename T>
Polinom<T>::Polinom(int grad, T coef[])
{
	this->grad = grad;
	this->coef = new T[this->grad + 1];
	for (int i = 0; i <= this->grad; i++)
		this->coef[i] = coef[i];
}

template <typename T>
Polinom<T>::Polinom(int grad, T c) {
	this->grad = grad;
	coef = new T[this->grad + 1];
	for (int i = 0; i < this->grad; i++)
		this->coef[i] = (T)0;
	this->coef[this->grad] = c;
}

template <typename T>
Polinom<T>::Polinom(const Polinom<T>& P)
{
	this->grad = P.grad;
	this->coef = new T[this->grad + 1];
	for (int i = 0; i <= this->grad; i++)
		this->coef[i] = P.coef[i];
}

template <typename T>
Polinom<T>::~Polinom()
{
	if (coef)
		delete[] coef;
}

template <typename T>
Polinom<T>& Polinom<T>::operator=(const Polinom<T>& P)
{
	if (this != &P)
	{
		if (coef)
			delete[] coef;
		grad = P.grad;
		coef = new T[grad + 1];
		for (int i = 0; i <= P.grad; i++)
			coef[i] = P.coef[i];
	}
	return *this;
}

template <typename T>
Polinom<T> Polinom<T>::operator+(const Polinom<T>& P)
{
	int gr = grad > P.grad ? grad : P.grad;
	T* cf = new T[gr + 1];
	for (int i = 0; i <= gr; i++)
		cf[i] = (T)0;
	for (int i = 0; i <= grad; i++)
		cf[i] = cf[i] + coef[i];
	for (int i = 0; i <= P.grad; i++)
		cf[i] = cf[i] + P.coef[i];
	while (gr > 0 && cf[gr] == (T)0)
		gr--;
	Polinom p(gr, cf);
	delete[] cf;
	return p;
}

template <typename T>
Polinom<T> Polinom<T>::operator-(const Polinom<T>& P)
{
	int gr = grad > P.grad ? grad : P.grad;
	T* cf = new T[gr + 1];
	for (int i = 0; i <= gr; i++)
		cf[i] = (T)0;
	for (int i = 0; i <= grad; i++)
		cf[i] = cf[i] + coef[i];
	for (int i = 0; i <= P.grad; i++)
		cf[i] = cf[i] - P.coef[i];
	while (gr > 0 && cf[gr] == (T)0)gr--;
	Polinom p(gr, cf);
	delete[] cf;
	return p;
}

template <typename T>
Polinom<T> Polinom<T>::operator*(const Polinom<T>& P)
{
	int gr = grad + P.grad;
	T* cf = new T[gr + 1];
	for (int i = 0; i <= gr; i++)
		cf[i] = (T)0;
	for (int i = 0; i <= grad; i++)
		for (int j = 0; j <= P.grad; j++)
			cf[i + j] += coef[i] * P.coef[j];
	while (gr > 0 && cf[gr] == (T)0)gr--;
	Polinom p(gr, cf);
	delete[] cf;
	return p;
}

template <typename T>
Polinom<T> Polinom<T>::operator/(const Polinom<T>& P)
{
	if (grad - P.grad < 0)return Polinom();
	int gr = -1;
	T* cf = new T[grad - P.grad + 1];
	Polinom Q(*this);
	while (Q.grad >= P.grad) {
		cf[grad - P.grad - ++gr] = Q.coef[Q.grad] / P.coef[P.grad];
		Q -= Polinom(Q.grad - P.grad, cf[grad - P.grad - gr]) * P;
	}
	Polinom p(gr, cf);
	delete[] cf;
	return p;
}

template <typename T>
T Polinom<T>::operator==(const Polinom<T>& P) {
	if (grad != P.grad)return 0;
	for (int i = 0; i <= grad; i++) {
		if (coef[i] != P.coef[i])return 0;
	}
	return 1;
}

template <typename T>
T Polinom<T>::operator!=(const Polinom<T>& P) {
	return !(*this == P);
}

template <typename T>
Polinom<T> Polinom<T>::operator-() {
	Polinom P = *this;
	for (int i = 0; i <= P.grad; i++)
		P.coef[i] = -P.coef[i];
	return P;
}

template <typename T>
Polinom<T>& Polinom<T>::operator+=(const Polinom<T>& P) {
	*this = *this + P;
	return *this;
}

template <typename T>
Polinom<T>& Polinom<T>::operator-=(const Polinom<T>& P) {
	*this = *this - P;
	return *this;
}

template <typename T>
Polinom<T>& Polinom<T>::operator*=(const Polinom<T>& P) {
	*this = *this * P;
	return *this;
}

template <typename T>
Polinom<T>& Polinom<T>::operator/=(const Polinom<T>& P) {
	*this = *this / P;
	return *this;
}

template <typename T>
T& Polinom<T>::operator[](int i) {
	return coef[i];
}

template <typename T>
T Polinom<T>::operator()(const T& c) {
	T rez = (T)0;

	for (int i = grad; i >= 0; i--)
		rez = rez * c + coef[i];

	return rez;
}

template <typename T>
ostream& operator<<(ostream& out, const Polinom<T>& P) {
	if (P.grad == 0)
		out << P.coef[0];
	else if (P.grad == 1) {
		out << "(" << P.coef[1] << ")X";
		if (P.coef[0] != (T)0)
			out << " + " << P.coef[0];
	}
	else {
		out << "(" << P.coef[P.grad] << ")X^" << P.grad;

		for (int i = P.grad - 1; i >= 2; i--) {
			if (P.coef[i] != (T)0)
				out << " + (" << P.coef[i] << ")X^" << i;
		}
		if (P.coef[1] != (T)0)
			out << " + (" << P.coef[1] << ")X";
		if (P.coef[0] != (T)0)
			out << " + " << P.coef[0];
	}
	return out;
}

template <typename T>
istream& operator>>(istream& in, Polinom<T>& P) {

	char linie[1000];
	in.getline(linie, 1000);
	int gr = -1;
	T cf[100];
	char* pt;
	pt = strtok(linie, " ");
	while (pt) {//citesc totul ca un string si apoi il impart

		std::istringstream ss(pt);
		ss >> cf[++gr];


		pt = strtok(NULL, " ");

	}
	P = Polinom<T>(gr, cf);

	return in;
}

int main()
{
	Polinom<int> P, R;
	cout << "Dati parametrii pentru primul polinom ";
	cin >> P;
	cout << "P: " << P << endl;
	cout << "Dati parametrii pentru al doilea polinom ";
	cin >> R;
	cout << "R: " << R << endl;


	cout << "P*R: " << P * R << endl;
	cout << "P/R: " << P / R << endl;
	cout << "P+R: " << P + R << endl;
	cout << "P-R: " << P - R << endl;
	cout << endl;

	Polinom<int>T = P + R;
	cout << T << endl;
	system("pause");
	return 0;
}