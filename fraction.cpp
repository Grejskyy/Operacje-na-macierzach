#include <boost/multiprecision/gmp.hpp>

using namespace boost::multiprecision;
using namespace std;

mpz_int gcd(mpz_int a, mpz_int b) {
    return b == 0 ? a : gcd(b, a % b);
}

class Fraction {
    public:
        mpz_int numerator, denominator;

        Fraction() {
            numerator = 0;
            denominator = 1;
        }

        Fraction(mpz_int n, mpz_int d) {
            if (d==0) {
                cerr << "Denominator may not be 0." << endl;
                exit(0);
            } else if (n == 0) {
                numerator = 0;
                denominator = 1;
            } else {
                int sign = 1;
                if (n < 0) {
                    sign *= -1;
                    n *= -1;
                }
                if (d < 0) {
                    sign *= -1;
                    d *= -1;
                }

                mpz_int tmp = gcd(n, d);
                numerator = n/tmp*sign;
                denominator = d/tmp;
            }
        }

        operator int() {return numerator.convert_to<int>()/denominator.convert_to<int>();}
        operator float() {return numerator.convert_to<float>()/denominator.convert_to<float>();}
        operator double() {return numerator.convert_to<double>()/denominator.convert_to<double>();}
        operator long double() {return numerator.convert_to<long double>()/denominator.convert_to<long double>();}
};
template <typename T>
T giveFloat(vector <T> data, int a, int b){
    Fraction temp(a,b);
    return (T)temp;
}

Fraction gaussMath(Fraction &m, Fraction division, Fraction temp){
    Fraction tmp(m.numerator * division.denominator * temp.denominator - division.numerator * temp.numerator * m.denominator, m.denominator * (division.denominator * temp.denominator));
    return tmp;
}
template<typename T>
T gaussMath(T m, T division, T temp){
    return m-division*temp;
}

mpf_float CheckErrorValue(mpf_float eigen, Fraction fraction){
    mpf_float ourFunction = fraction.numerator.convert_to<mpf_float>()/fraction.denominator.convert_to<mpf_float>();
    mpf_float value = eigen.convert_to<mpf_float>() - ourFunction;
    return value < 0 ? value *= -1 : value;
}
template<typename T>
mpf_float CheckErrorValue(mpf_float eigen, T dValue){;
    mpf_float value = eigen - dValue;
    return value < 0 ? value *= -1 : value;
}

Fraction abs(const Fraction& obj){
    Fraction temp(abs(obj.numerator),abs(obj.denominator));
    return temp;
}

Fraction operator+(const Fraction& lhs, const Fraction& rhs) {
    Fraction tmp(lhs.numerator*rhs.denominator+rhs.numerator*lhs.denominator,
                lhs.denominator*rhs.denominator);
    return tmp;
}

Fraction operator+=(Fraction& lhs, const Fraction& rhs) {
    Fraction tmp(lhs.numerator*rhs.denominator
                +rhs.numerator*lhs.denominator,
                lhs.denominator*rhs.denominator);
    lhs = tmp;
    return lhs;
}

Fraction operator-(const Fraction& lhs, const Fraction& rhs) {
    Fraction tmp(lhs.numerator*rhs.denominator
                -rhs.numerator*lhs.denominator,
                lhs.denominator*rhs.denominator);
    return tmp;
}

Fraction operator-=(Fraction& lhs, const Fraction& rhs) {
    Fraction tmp(lhs.numerator*rhs.denominator
                -rhs.numerator*lhs.denominator,
                lhs.denominator*rhs.denominator);
    lhs = tmp;
    return lhs;
}

Fraction operator*(const Fraction& lhs, const Fraction& rhs) {
    Fraction tmp(lhs.numerator*rhs.numerator,
               lhs.denominator*rhs.denominator);
    return tmp;
}

Fraction operator*=(Fraction& lhs, const Fraction& rhs) {
    Fraction tmp(lhs.numerator*rhs.numerator,
               lhs.denominator*rhs.denominator);
    lhs = tmp;
    return lhs;
}

Fraction operator*(int lhs, const Fraction& rhs) {
    Fraction tmp(lhs*rhs.numerator,rhs.denominator);
    return tmp;
}

Fraction operator*(const Fraction& rhs, int lhs) {
    Fraction tmp(lhs*rhs.numerator,rhs.denominator);
    return tmp;
}

Fraction operator/(const Fraction& lhs, const Fraction& rhs) {
    Fraction tmp(lhs.numerator*rhs.denominator,
                 lhs.denominator*rhs.numerator);
    return tmp;
}

bool operator>(const Fraction& lhs, const Fraction& rhs) {
    Fraction temp1 = lhs;
    Fraction temp2 = rhs;
    return (long double)temp1 > (long double)temp2;
}

std::ostream& operator<<(std::ostream &strm, const Fraction &a) {
    if (a.denominator == 1) {
        strm << a.numerator;
    } else {
        strm << a.numerator << "/" << a.denominator;
    }
    return strm;
}