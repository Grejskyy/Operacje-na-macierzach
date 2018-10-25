#include <iostream>
#include <vector>
#include <algorithm>
#include <random>

using namespace std;

template<typename T>
class Matrix {
	protected:
		unsigned rows, cols;
		vector <T> data; 
	public:
		Matrix(unsigned rows, unsigned cols)
			: rows(rows)
			, cols(cols)
			{
				random_device rd;  
    			mt19937 gen(rd()); 
    			uniform_real_distribution<> dis(-1, 1);
				for (int i = 0; i < rows*cols; ++i)
				{
					data.push_back(dis(gen));
				}
			}
		Matrix(const Matrix &obj)
		: rows(obj.rows)
		, cols(obj.cols)
		{
			data = obj.data;
		}
		~Matrix()
		{
		}
		T& operator() (unsigned row, unsigned col){
			return data[cols*row + col];
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
  	int n = 0;
  	while(n < m.getRows()){
		for (int i = n+1; i < m.getRows(); ++i)
		{
			for (int j = n; j < m.getCols(); ++j)
			{
				m(i,j) -= temp(i,n)/temp(n,n)*temp(n,j);
			}
		}
		n++;
		temp = m;
	}
	printX(m);
}
template <typename T>
void Gaussb(Matrix <T> &m){
	Matrix <T> temp(m);
	
  	int n = 0;
  	int hlp=0;
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
			for (int j = n; j < m.getCols(); ++j)
			{
				m(i,j) -= temp(i,n)/temp(n,n)*temp(n,j);
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
  	while(n < m.getRows()){
  		hlp=0;
  		hlp2=0;
  		for (int i = n; i < m.getRows(); ++i)
  			for (int j = n; j < m.getCols()-1; ++j)
  				if(m(i,j)>m(hlp,hlp2)){
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
			for (int j = n; j < m.getCols(); ++j)
			{
				m(i,j) -= temp(i,n)/temp(n,n)*temp(n,j);
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
		value /= m(row - i, col - (i+1));
		x.push_back(value);
	}
	reverse(x.begin(), x.end());
	unsigned iter = 1;
	for (it=x.begin(); it!=x.end(); ++it){
    	cout<< "x"  << iter << ": " << *it << endl;
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
		value /= m(row - i, col - (i+1));
		x.push_back(value);
	}
	reverse(x.begin(), x.end());
	for (int i = 0; i < x.size(); ++i)
	{
		for (int j = 0; j < x.size(); ++j)
		{
			if(colPosition[j] == i){
					cout<< "x"  << i+1 << ": " << x[j] << endl;
					break; 
			}
		}
	}
}
int main(int argc, char const *argv[])
{
	Matrix <float> m(1500,1501);
 	Matrix <float> n=m;
	Gauss(m);
	Gaussc(n);        
	return 0;
}
