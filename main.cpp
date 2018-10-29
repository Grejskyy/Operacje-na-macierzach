// Trzeba ogarnąć fraction, próbowałem cpp.int ale: 
// -int128 - poprawnie liczy tylko przy GaussC
// -int256 - wolniejszy od mpz_int
// Rzutowanie chyba okej
// Eigen ale dopiero po ogarnięciu fraction, nwm może Kuszner coś kojarzy dx
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <ctime>
#include "fraction.cpp"

using namespace boost::multiprecision;
using namespace std;
 
template<typename T>
class Matrix {
	protected:
		
		
	public:
		unsigned rows, cols;
		vector <T> data; 
		Matrix(unsigned rows, unsigned cols)
			: rows(rows)
			, cols(cols)
			{
				random_device rd;
    			mt19937 gen(rd());
    			uniform_int_distribution<> dis(-65536, 65535);
				for (int i = 0; i < rows*cols; ++i)
				{
<<<<<<< HEAD
					data.push_back(giveFloat(data,dis(gen),65536));
=======
					data.push_back(giveFraction(dis(gen),65536));
					//cout << i << endl;
>>>>>>> 7e56cfc8aa800353fb4ee970e1712911bc69a959
				}
			}
		Matrix(const Matrix &obj)
		: rows(obj.rows)
		, cols(obj.cols)
		{
			data = obj.data;
		}
		Matrix(Matrix <Fraction> &obj)
		: rows(obj.rows)
		, cols(obj.cols)
		{
			for (int i = 0; i < rows*cols; ++i)
		 	{
			data.push_back((T)obj.data[i]);
		 	}
		}
		~Matrix()
		{
		}
		T& operator() (unsigned row, unsigned col){
			return data[cols*row + col];
		}
		T& operator=(const Matrix &obj)
		{
			rows = obj.rows;
		 	cols = obj.cols;
		 	for (int i = 0; i < rows*cols; ++i)
		 	{
			data[i] = (T)obj.data[i];
		 	}
		}
		unsigned getRows(){
			return rows;
		}
		unsigned getCols(){
			return cols;
		}
};
template<typename T>
void Gauss(Matrix <T> &m){
	Matrix <T> temp(m);
	T division;
  	for (int n = 0; n < m.getRows(); ++n)
  	{
		for (int i = n+1; i < m.getRows(); ++i)
		{	
<<<<<<< HEAD
			division = temp(i,n)/temp(n,n);
=======
>>>>>>> 7e56cfc8aa800353fb4ee970e1712911bc69a959
			for (int j = n; j < m.getCols(); ++j)
			{
				m(i,j) = gaussMath(m(i,j),division,temp(n,j));
			}
		}
		cout << "Test n = " << n << endl;
<<<<<<< HEAD
=======
		n++;
		temp = m;
	}
	printX(m);
}

template<typename T>
void GaussDouble(Matrix <T> &m){
	Matrix <T> temp(m);
  	int n = 0;
	double div = 0;
  	while(n < m.getRows()){
		for (int i = n+1; i < m.getRows(); ++i)
		{	
			div = temp(i,n)/temp(n,n);
			for (int j = n; j < m.getCols(); ++j)
			{
				m(i,j) -= div*temp(n,j);
			}
		}
		cout << "Test n = " << n << endl;
		n++;
		temp = m;
	}
	printX(m);
}

template<typename T>
void GaussFr(Matrix <T> &m){
	Matrix <T> temp(m);
  	int n = 0;
	Fraction division(1,1);
  	while(n < m.getRows()){
		for (int i = n+1; i < m.getRows(); ++i)
		{	
			division = temp(i,n)/temp(n,n);
			for (int j = n; j < m.getCols(); ++j)
			{
				m(i,j) -= division*temp(n,j);
			}
		}
		cout << "Test n = " << n << endl;
		n++;
		temp = m;
	}
	printX(m);
}

template<typename T>
void GaussFr2(Matrix <T> &m){
	Matrix <T> temp(m);
  	int n = 0;
	Fraction division(1,1);
  	while(n < m.getRows()){
		for (int i = n+1; i < m.getRows(); ++i)
		{	
			division = temp(i,n)/temp(n,n);
			for (int j = n; j < m.getCols(); ++j)
			{
				m(i,j) = gaussMath(m(i,j),division,temp(n,j));
			}
		}
		cout << "Test n = " << n << endl;
		n++;
>>>>>>> 7e56cfc8aa800353fb4ee970e1712911bc69a959
		temp = m;
	}
	printX(m);
}

template <typename T>
void Gaussb(Matrix <T> &m){
	Matrix <T> temp(m);
  	int n = 0;
  	int hlp=0;
  	T division;
  	while(n < m.getRows()){
  		hlp=0;
  		for (int i = n; i < m.getRows(); ++i)
  			if(m(i,n)>m(hlp,n))
				hlp=i;
  		if(hlp!=0){
  			for (int i = 0; i < m.getCols(); ++i)
  			{
  				m(n,i) = temp(hlp,i);
  				m(hlp,i) = temp(n,i);
				  		
  			}
  			temp = m;
  		}
		for (int i = n+1; i < m.getRows(); ++i)
		{	
			division = temp(i,n)/temp(n,n);
			for (int j = n; j < m.getCols(); ++j)
			{
				m(i,j) = gaussMath(m(i,j),division,temp(n,j));
			}
		}
		n++;
		temp = m;
	}
	printX(m);
}

template <typename T>
void Gaussc(Matrix <T> &m){
	Matrix <T> temp(m);
	
  	int n = 0;
  	int hlp=0;
  	int hlp2=0;
  	vector<unsigned> colPosition;
  	for (int i = 0; i <= m.getRows(); ++i)
  	{
  		colPosition.push_back(i);
  	}
  	T division;
  	while(n < m.getRows()){
  		hlp=0;
  		hlp2=0;
  		for (int i = n; i < m.getRows(); ++i)
  			for (int j = n; j < m.getCols()-1; ++j)
  				if(abs(m(i,j))>abs(m(hlp,hlp2))){
					hlp = i;
					hlp2 = j;
  				}
  		if(hlp2!=0){
  			int a;
  			a= colPosition[hlp2];
  			colPosition[hlp2] = colPosition[n];
  			colPosition[n] = a; 
  			for (int i = 0; i < m.getRows(); ++i)
  			{
  				m(i,n) = temp(i,hlp2);
  				m(i,hlp2) = temp(i,n);	
  				
  			}
  			temp = m;
  		}
  		if(hlp!=0){
  			for (int i = 0; i < m.getCols(); ++i)
  			{
  				m(n,i) = temp(hlp,i);
  				m(hlp,i) = temp(n,i);	
  			}
  			temp = m;
  		}
		for (int i = n+1; i < m.getRows(); ++i)
		{	
			division = temp(i,n)/temp(n,n);
			for (int j = n; j < m.getCols(); ++j)
			{
				m(i,j) = gaussMath(m(i,j),division,temp(n,j));
			}
		}
		n++;
		temp = m;
	}
  	cout << endl;
	printX(m, colPosition);
}

template <typename T>
void printX(Matrix <T> &m){
	vector <T> x;
	typename vector <T>::iterator it; 
	unsigned row = m.getRows();
	unsigned col = m.getCols();
	T value;
	for (int i = 1; i <=row; ++i)
	{
		value = m(row - i,col - 1);
		for (int j = 2; j <= i; ++j)
		{
			value -= m(row - i, col - j)*x[j-2];
		}
		value = value / m(row - i, col - (i+1));
		x.push_back(value);
	}
	reverse(x.begin(), x.end());
	unsigned iter = 1;
	for (it=x.begin(); it!=x.end(); ++it){
    	cout<< "x"  << iter << ": " << (long double)*it << endl;
    	iter++;
	}
}

template <typename T>
void printX(Matrix <T> &m, vector<unsigned> colPosition){
	vector <T> x;
	typename vector <T>::iterator it; 
	unsigned row = m.getRows();
	unsigned col = m.getCols();
	T value;
	for (int i = 1; i <=row; ++i)
	{
		value = m(row - i,col - 1);
		for (int j = 2; j <= i; ++j)
		{
			value -= m(row - i, col - j)*x[j-2];
		}
		value = value / m(row - i, col - (i+1));
		x.push_back(value);
	}
	reverse(x.begin(), x.end());
	for (int i = 0; i < x.size(); ++i)
	{
		for (int j = 0; j < x.size(); ++j)
		{
			if(colPosition[j] == i){
					cout<< "x"  << i+1 << ": " << (long double)x[j] << endl;
					break; 
			}
		}
	}
}

int main(int argc, char const *argv[])
{
<<<<<<< HEAD
	int size = 40;
	time_t start, end;
    Matrix <Fraction> U1(size,size+1);
    Matrix <Fraction> U2(U1);
    Matrix <Fraction> U3(U2);
    Matrix <double> D1(U1);
    Matrix <double> D2(D1);
    Matrix <double> D3(D2);
    Matrix <float> F1(U1);
    Matrix <float> F2(F1);
    Matrix <float> F3(F2);
	start = clock();
	Gauss(U1);
	end = clock();
	double U1Time = (static_cast <double>(end - start) / CLOCKS_PER_SEC);
	start = clock();
	Gaussb(U2);
	end = clock();
	double U2Time = (static_cast <double>(end - start) / CLOCKS_PER_SEC);
	start = clock();
	Gaussc(U3);
	end = clock();
	double U3Time = (static_cast <double>(end - start) / CLOCKS_PER_SEC);
	start = clock();
	Gauss(D1);
	end = clock();
	double D1Time = (static_cast <double>(end - start) / CLOCKS_PER_SEC);
	start = clock();
	Gaussb(D2);
	end = clock();
	double D2Time = (static_cast <double>(end - start) / CLOCKS_PER_SEC);
	start = clock();
	Gaussc(D3);
	end = clock();
	double D3Time = (static_cast <double>(end - start) / CLOCKS_PER_SEC);
	start = clock();
	Gauss(F1);
	end = clock();
	double F1Time = (static_cast <double>(end - start) / CLOCKS_PER_SEC);
	start = clock();
	Gaussb(F2);
	end = clock();
	double F2Time = (static_cast <double>(end - start) / CLOCKS_PER_SEC);
	start = clock();
	Gaussc(F3);
	end = clock();
	double F3Time = (static_cast <double>(end - start) / CLOCKS_PER_SEC);
	start = clock();
	cout << "Matrix size: " << size << "x" << size << endl;
	cout << "Fraction Gauss1: " << U1Time << endl;
	cout << "Fraction Gauss2: " << U2Time << endl;
	cout << "Fraction Gauss3: " << U3Time << endl;
	cout << "Double Gauss1: " << D1Time << endl;
	cout << "Double Gauss2: " << D2Time << endl;
	cout << "Double Gauss3: " << D3Time << endl;
	cout << "Float Gauss1: " << F1Time << endl;
	cout << "Float Gauss2: " << F2Time << endl;
	cout << "Float Gauss3: " << F3Time << endl;
=======
	time_t start, end;
    Matrix <Fraction> m(500,501);
					// m.data.push_back(giveFraction(-1,1));
					// m.data.push_back(giveFraction(2,1));
					// m.data.push_back(giveFraction(1,1));
					// m.data.push_back(giveFraction(-1,1));
					// m.data.push_back(giveFraction(1,1));
					// m.data.push_back(giveFraction(-3,1));
					// m.data.push_back(giveFraction(-2,1));
					// m.data.push_back(giveFraction(-1,1));
					// m.data.push_back(giveFraction(3,1));
					// m.data.push_back(giveFraction(-1,1));
					// m.data.push_back(giveFraction(-1,1));
					// m.data.push_back(giveFraction(4,1));
//CHWILOWE RZUTOWANIE NA ODPIERDOL					
    Matrix <double> n(0,0);
     n.rows = m.rows;
	n.cols = m.cols;
    for (int i = 0; i < m.rows*m.cols; ++i)
    {
    	n.data.push_back((double)m.data[i]);
    }
	cout << "Gauss:" << endl;
	/*start = clock();
	Gauss(n);
	end = clock();
	double db1time = (static_cast <double>(end - start) / CLOCKS_PER_SEC);*/

	start = clock();
	GaussDouble(n);
	end = clock();
	double db2time = (static_cast <double>(end - start) / CLOCKS_PER_SEC);

	/*start = clock();
	GaussFr(m);
	end = clock();
	double fr1time = (static_cast <double>(end - start) / CLOCKS_PER_SEC);*/

	start = clock();
	GaussFr2(m);
	end = clock();
	//cout << "Double db1: " << db1time << endl;
	cout << "Double db2: " << db2time << endl;
	//cout << "Fraction Fr: " << fr1time << endl;
	cout << "Fraction Fr2: " <<(static_cast <double>(end - start) / CLOCKS_PER_SEC) << endl;
>>>>>>> 7e56cfc8aa800353fb4ee970e1712911bc69a959
    return 0;
}