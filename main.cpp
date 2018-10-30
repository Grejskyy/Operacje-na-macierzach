#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <sstream>
#include <fstream>
#include <string>
#include <ctime>
#include <Eigen/Dense>
#include "fraction.cpp"

using namespace boost::multiprecision;
using namespace std;

mpf_float avError = 0;
 
template<typename T>
class MyMatrix {
	protected:
		
		
	public:
		unsigned rows, cols;
		vector <T> data; 
		MyMatrix(unsigned rows, unsigned cols)
			: rows(rows)
			, cols(cols)
			{
				random_device rd;
    			mt19937 gen(rd());
    			uniform_int_distribution<> dis(-65536, 65535);
				for (int i = 0; i < rows*cols; ++i)
				{
					data.push_back(giveFloat(data,dis(gen),65536));
				}
			}
		MyMatrix(const MyMatrix &obj)
		: rows(obj.rows)
		, cols(obj.cols)
		{
			data = obj.data;
		}
		MyMatrix(MyMatrix <Fraction> &obj)
		: rows(obj.rows)
		, cols(obj.cols)
		{
			for (int i = 0; i < rows*cols; ++i)
		 	{
			data.push_back((T)obj.data[i]);
		 	}
		}
		~MyMatrix()
		{
		}
		T& operator() (unsigned row, unsigned col){
			return data[cols*row + col];
		}
		T& operator=(const MyMatrix &obj)
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
void Gauss(MyMatrix <T> &m, Eigen::MatrixXd x){
	MyMatrix <T> temp(m);
	T division;
  	for (int n = 0; n < m.getRows(); ++n)
  	{
		for (int i = n+1; i < m.getRows(); ++i)
		{	
			division = temp(i,n)/temp(n,n);
			for (int j = n; j < m.getCols(); ++j)
			{	
				m(i,j) = gaussMath(m(i,j),division,temp(n,j));
			}
		}
		temp = m;
	}
	checkErrorX(m,x);
}

template <typename T>
void Gaussb(MyMatrix <T> &m, Eigen::MatrixXd x){
	MyMatrix <T> temp(m);
  	int n = 0;
  	int hlp=0;
  	T division;
  	while(n < m.getRows()){
  		hlp=0;
  		for (int i = n; i < m.getRows(); ++i)
  			if(abs(m(i,n))>abs(m(hlp,n)))
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
	checkErrorX(m,x);
}

template <typename T>
void Gaussc(MyMatrix <T> &m, Eigen::MatrixXd x){
	MyMatrix <T> temp(m);
	
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
	checkErrorX(m, colPosition, x);
}

template <typename T>
void checkErrorX(MyMatrix <T> &m, Eigen::MatrixXd y){
	vector <T> x;
	typename vector <T>::iterator it; 
	unsigned row = m.getRows();
	unsigned col = m.getCols();
	T value;
	avError = 0;
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
	unsigned i = 0;
	for (it=x.begin(); it!=x.end(); ++it){
		avError += CheckErrorValue(y(i),*it);
    	i++;
	}
}

template <typename T>
void checkErrorX(MyMatrix <T> &m, vector<unsigned> colPosition, Eigen::MatrixXd y){
	vector <T> x;
	typename vector <T>::iterator it; 
	unsigned row = m.getRows();
	unsigned col = m.getCols();
	T value;
	avError = 0;
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
				avError += CheckErrorValue(y(i),x[j]);
				break; 
			}
		}
	}
}

int main(int argc, char const *argv[])
{
	istringstream ss(argv[1]);
	int size;
	ss >> size;
	time_t start, end;
    MyMatrix <Fraction> U1(size,size+1);
    MyMatrix <Fraction> U2(U1);
    MyMatrix <Fraction> U3(U2);
    MyMatrix <double> D1(U1);
    MyMatrix <double> D2(D1);
    MyMatrix <double> D3(D2);
    MyMatrix <float> F1(U1);
    MyMatrix <float> F2(F1);
    MyMatrix <float> F3(F2);
    Eigen::MatrixXd eig(size, size);
    Eigen::VectorXd vec(size);

    for (int i = 0; i < size; ++i)
    {
    	for (int j = 0; j < size; ++j)
    	{
    		eig(i,j) = D1(i,j);
    	}
    	vec(i) = D1(i,size);
    }
	
    Eigen::MatrixXd x = eig.colPivHouseholderQr().solve(vec);

	start = clock();
	Gauss(F1,x);
	end = clock();
	mpf_float GaussFlAvError = avError/size;
	double F1Time = (static_cast <double>(end - start) / CLOCKS_PER_SEC);
	start = clock();
	Gaussb(F2,x);
	end = clock();
	mpf_float Gauss2FlAvError = avError/size;
	double F2Time = (static_cast <double>(end - start) / CLOCKS_PER_SEC);
	start = clock();
	Gaussc(F3,x);
	end = clock();
	mpf_float Gauss3FlAvError = avError/size;
	double F3Time = (static_cast <double>(end - start) / CLOCKS_PER_SEC);

	start = clock();
	Gauss(D1,x);
	end = clock();
	mpf_float GaussDbAvError = avError/size;
	double D1Time = (static_cast <double>(end - start) / CLOCKS_PER_SEC);
	start = clock();
	Gaussb(D2,x);
	end = clock();
	mpf_float Gauss2DbAvError = avError/size;
	double D2Time = (static_cast <double>(end - start) / CLOCKS_PER_SEC);
	start = clock();
	Gaussc(D3,x);
	end = clock();
	mpf_float Gauss3DbAvError = avError/size;
	double D3Time = (static_cast <double>(end - start) / CLOCKS_PER_SEC);

	start = clock();
	Gauss(U1,x);
	end = clock();
	mpf_float GaussFrAvError = avError/size;
	double U1Time = (static_cast <double>(end - start) / CLOCKS_PER_SEC);
	start = clock();
	Gaussb(U2,x);
	end = clock();
	mpf_float Gauss2FrAvError = avError/size;
	double U2Time = (static_cast <double>(end - start) / CLOCKS_PER_SEC);
	start = clock();
	Gaussc(U3,x);
	end = clock();
	mpf_float Gauss3FrAvError = avError/size;
	double U3Time = (static_cast <double>(end - start) / CLOCKS_PER_SEC);
	
	string fileName = "matrix" + to_string(size) + ".csv";
	ofstream csvFile;
    csvFile.open (fileName);
    csvFile << "Float Gauss1 Time;Float Gauss2 Time;Float Gauss3 Time;Double Gauss1 Time;Double Gauss2 Time;Double Gauss3 Time;Fraction Gauss1 Time;Fraction Gauss2 Time;Fraction Gauss3 Time;Size "<< size << "x" << size << ";\n";
    csvFile << F1Time << ";" << F2Time << ";"<< F3Time << ";" << D1Time << ";" << D2Time << ";"<< D3Time << ";" << U1Time << ";" << U2Time << ";"<< U3Time << ";" << "\n";
    csvFile << "\n";
    csvFile << "Float Gauss1 Average Error;Float Gauss2 Average Error;Float Gauss3 Average Error;Double Gauss1 Average Error;Double Gauss2 Average Error;Double Gauss3 Average Error;Fraction Gauss1 Average Error;Fraction Gauss2 Average Error;Fraction Gauss3 Average Error;\n";
    csvFile << GaussFlAvError << ";" << Gauss2FlAvError << ";"<< Gauss3FlAvError << ";" << GaussDbAvError << ";" << Gauss2DbAvError << ";"<< Gauss3DbAvError << ";" << GaussFrAvError << ";" << Gauss2FrAvError << ";"<< Gauss3FrAvError << ";" << "\n";
    csvFile.close();
    return 0;
}