#include<iostream>
#include<sstream>
#include<string>
#include<iomanip>
#include<cmath>
#define Pi 3.141592653589793
typedef long T;
using namespace std;

class BigFloat;
struct complex //定义复数结构体
{
	double re, im;
	complex(double r = 0.0, double i = 0.0) { re = r, im = i; }  //初始化
	 //定义三种运算 
	complex operator +(complex com)
	{
		return complex(re + com.re, im + com.im);
	}
	complex operator -(complex com)
	{
		return complex(re - com.re, im - com.im);
	}
	complex operator *(complex com)
	{
		return complex(re * com.re - im * com.im, re * com.im + im * com.re);
	}
	friend  void BRC(complex* y, int len);
	friend	void FFT(complex* y, int len, double on);
};
void BRC(complex* y, int len)//二进制反转倒置
{
	int i, j, k;
	for (i = 1, j = len / 2; i < len - 1; i++)
	{
		if (i < j)swap(y[i], y[j]);//i<j保证只交换一次
		k = len >> 1;
		while (j >= k)
		{
			j -= k; k = k >> 1;
		}
		if (j < k)j += k;
	}
}

void FFT(complex* y, int len, double on)//on=1表示顺，-1表示逆
{
	int i, j, k, h;
	complex u, t;
	BRC(y, len);
	for (h = 2; h <= len; h <<= 1)//控制层数
	{
		//初始化单位复根
		complex wn(cos(on * 2 * Pi / h), sin(on * 2 * Pi / h));
		for (j = 0; j < len; j += h)//控制起始下标
		{
			//初始化螺旋因子
			complex w(1, 0);
			for (k = j; k < j + h / 2; k++)
			{
				u = y[k];
				t = w * y[k + h / 2];
				y[k] = u + t;
				y[k + h / 2] = u - t;
				w = w * wn;//更新螺旋因子
			}
		}
	}
	if (on == -1)
		for (i = 0; i < len; i++) //逆FFT(IDFT)
			y[i].re /= len;

}
/****************************************
*BigInt是实现任意精度的整数的四则运算类 *
*copyright：南京工程学院                *
*author:网络工程181马新尧               *
****************************************/
class BigInt
{
public:
	enum myconst { SPACE_CAPACITY = 4, Int_Size = 4, num_Weight = 10000 };
	//SPACE_CAPACITY冗余容量；Int_Size位权位数；num_Weight位权   
	BigInt(string numstr); //构造函数：用字符串初始化
	BigInt(const BigInt& bi);    //复制构造函数：用BigInt对象初始化
	BigInt(T n = 0, T num_Capacity = SPACE_CAPACITY);//构造函数：容量为num_Capacity值为0
	~BigInt() { delete[] num; }
	static int BigIntcmp(BigInt& bi1, BigInt& bi2);//比较绝对值大小，如果bi1>bi2返回1；如果bi1=bi2返回0；如果bi1<bi2返回-1
	BigInt operator =(BigInt bi);
	BigInt operator +=(BigInt bi) { return *this = *this + bi; }
	BigInt operator -=(BigInt bi) { return *this = *this - bi; }
	BigInt operator *=(BigInt bi) { return *this = *this * bi; }
	BigInt operator /=(BigInt bi) { return *this = *this / bi; }
	BigInt operator -() { BigInt temp = *this; temp.num[0] = (temp.num[0] + 1) % 2; return temp; }
	BigInt& operator ++() { *this = *this + (BigInt)1; return *this; }
	BigInt& operator --() { *this = *this - (BigInt)1; return *this; }
	BigInt operator ++(int) { BigInt temp = *this; *this = *this + (BigInt)1; return temp; }
	BigInt operator --(int) { BigInt temp = *this; *this = *this - (BigInt)1; return temp; }
	friend BigInt operator >>(BigInt bi, T a);
	friend BigInt operator <<(BigInt& bi, T a);
	friend istream& operator >>(istream& input, BigInt& bi);
	friend ostream& operator <<(ostream& output, BigInt& bi);
	friend BigInt operator +(BigInt bi1, BigInt bi2);
	friend BigInt operator -(BigInt bi1, BigInt bi2);
	friend BigInt operator *(BigInt bi1, BigInt bi2);
	friend BigInt operator /(BigInt bi1, BigInt bi2);
	friend BigInt operator %(BigInt bi1, BigInt bi2);
	friend bool operator == (BigInt& bi1, BigInt bi2);
	friend bool operator != (BigInt& bi1, BigInt bi2) { return !operator == (bi1, bi2); }
	friend bool operator > (BigInt& bi1, BigInt bi2);
	friend bool operator < (BigInt& bi1, BigInt bi2);
	friend bool operator >= (BigInt& bi1, BigInt bi2) { return !operator < (bi1, bi2); }
	friend bool operator <= (BigInt& bi1, BigInt bi2) { return !operator > (bi1, bi2); }
	friend BigFloat;
	//void display();        //输出整数，初期函数，末期要重载<<运算符
	BigInt add(BigInt& bi1, BigInt& bi2);//加法函数：实现绝对值加法，末期要重载+运算符
	BigInt sub(BigInt& bi1, BigInt& bi2);//减法函数：实现绝对值大数减小数，末期要重载-运算符
	BigInt div(BigInt bi1, BigInt bi2, char flag);
private:
	void numstr_convert(string numstr);//将字符串转化为数组的函数，BigInt(string numstr)构造函数中调用
	void numstr_test(string numstr);//检测字符串是否是合法的整数
private:
	T* num;//指向整形数组的指针
	long num_Size;//整形数组中表示任意精度整数的元素个数
	long num_Capacity;//整形数组的大小
 /******************************************************
 字符串“-123456789”，在对象中
 表示为num[0]=1、num[1]=6789、num[2]=2345、num[3]=1，
 num[0]表示整数的符号，0表示正数、1表示负数，num_Size=3
 ******************************************************/
};
BigInt::BigInt(string numstr)
{
	long len = numstr.size();
	numstr_test(numstr);
	num_Capacity = len / Int_Size + SPACE_CAPACITY;
	num = new T[num_Capacity];
	if (numstr[0] == '-')
	{
		num[0] = 1;
		numstr_convert(numstr.substr(1));
	}
	else
	{
		num[0] = 0;
		numstr_convert(numstr);
	}
}

BigInt::BigInt(const BigInt& bi)
{
	num_Capacity = bi.num_Capacity;
	num_Size = bi.num_Size;
	num = new T[num_Capacity];
	for (int i = 0; i <= num_Size; i++)
		num[i] = bi.num[i];
}
BigInt::BigInt(T n, T num_Capacity)
{
	int i;
	this->num_Capacity = num_Capacity;
	num = new T[num_Capacity];
	n >= 0 ? num[0] = 0 : (num[0] = 1, n = -n);
	for (i = 1; i == 1 || n > 0; i++)
	{
		num[i] = n % num_Weight;
		n = n / num_Weight;
	}
	num_Size = i - 1;
}

BigInt BigInt::operator =(BigInt bi)
{
	num_Capacity = bi.num_Capacity;
	num_Size = bi.num_Size;
	delete[]num;//释放以前的内存，很重要，没有加时导致程序崩溃
	num = new T[num_Capacity];
	for (int i = 0; i <= num_Size; i++)
		num[i] = bi.num[i];
	return *this;
}
BigInt operator >>(BigInt bi, T a)
{
	if (a<0 || a>bi.num_Size)
	{
		cout << "超出移位范围！" << endl;
		return bi;
	}
	if (a == bi.num_Size)
	{
		bi.num[0] = 0;
		bi.num[1] = 0;
		bi.num_Size = 1;
		return bi;
	}
	bi.num_Size -= a;
	for (T i = 1; i <= bi.num_Size; i++)
		bi.num[i] = bi.num[i + a];

	return bi;
}
BigInt operator <<(BigInt& bi, T a)
{
	T i;
	if (a < 0)
	{
		cout << "超出移位范围！" << endl;
		return bi;
	}
	BigInt r(0, bi.num_Size + a + BigInt::SPACE_CAPACITY);
	r.num_Size = bi.num_Size + a;
	r.num[0] = bi.num[0];
	for (i = 1; i <= a; i++)
		r.num[i] = 0;
	for (; i <= r.num_Size; i++)
		r.num[i] = bi.num[i - a];

	return r;
}
istream& operator >>(istream& input, BigInt& bi)
{
	string st;
	input >> st;
	bi = BigInt(st);
	return input;
}
ostream& operator <<(ostream& output, BigInt& bi)
{
	if (bi.num[0] == 1) cout << '-';
	output << bi.num[bi.num_Size];
	for (int i = bi.num_Size - 1; i >= 1; i--)
		output << setw(bi.Int_Size) << setfill('0') << bi.num[i];
	return output;
}
BigInt operator +(BigInt bi1, BigInt bi2)
{
	int m = bi1.num_Size > bi2.num_Size ? bi1.num_Size : bi2.num_Size;
	BigInt bi(0, m + BigInt::SPACE_CAPACITY + 1);
	if (bi1.num[0] == bi2.num[0])
	{
		bi = bi.add(bi1, bi2);
		bi.num[0] = bi1.num[0];
	}
	else
	{
		if (BigInt::BigIntcmp(bi1, bi2) == 0)
		{
			bi.num[0] = 0;
			bi.num[1] = 0;
			bi.num_Size = 1;
		}
		if (BigInt::BigIntcmp(bi1, bi2) == 1)
		{
			bi = bi.sub(bi1, bi2);
			bi.num[0] = bi1.num[0];
		}
		if (BigInt::BigIntcmp(bi1, bi2) == -1)
		{
			bi = bi.sub(bi2, bi1);
			bi.num[0] = bi2.num[0];
		}
	}
	return bi;
}
BigInt operator -(BigInt bi1, BigInt bi2)
{
	bi2.num[0] = (bi2.num[0] + 1) % 2;

	return bi1 + bi2;
}
BigInt operator *(BigInt bi1, BigInt bi2)
{
	int len1 = bi1.num_Size, len2 = bi2.num_Size, len = 1;
	int i;
	if (bi1.num_Size == 1 && bi1.num[1] == 0) return bi1;//bi1==0 
	if (bi2.num_Size == 1 && bi2.num[1] == 0) return bi2;//bi2==0

	while (len < 2 * len1 || len < 2 * len2)len <<= 1;//右移相当于len=len*2
	complex* x1 = new complex[len + 1];
	complex* x2 = new complex[len + 1];
	T* sum = new T[len + 1];
	//倒置存储
	for (i = 1; i <= len1; i++)
	{
		x1[i].re = bi1.num[i]; x1[i].im = 0.0;
	}
	for (; i <= len; i++)  //多余次数界初始化为0
	{
		x1[i].re = x1[i].im = 0.0;
	}
	for (i = 1; i <= len2; i++)
	{
		x2[i].re = bi2.num[i];; x2[i].im = 0.0;
	}
	for (; i <= len; i++)  //多余次数界初始化为0
	{
		x2[i].re = x2[i].im = 0.0;
	}
	//FFT求值
	FFT(x1 + 1, len, 1);//FFT(a) 1表示顺 -1表示逆
	FFT(x2 + 1, len, 1);//FFT(b)
//点乘，结果存入x1
	for (i = 1; i <= len; i++)
		x1[i] = x1[i] * x2[i];
	//插值，逆FFT（IDTF）
	FFT(x1 + 1, len, -1);

	//细节处理
	for (i = 1; i <= len; i++)
		sum[i] = x1[i].re + 0.5;//四舍五入
	for (i = 1; i < len; i++)     //进位
	{
		sum[i + 1] += sum[i] / BigInt::num_Weight;
		sum[i] %= BigInt::num_Weight;
	}
	//输出

	BigInt result;
	result.num_Capacity = len + 1;
	len = len1 + len2;
	while (sum[len] == 0 && len > 0)len--;//检索最高位
	result.num_Size = len;

	result.num = sum;
	result.num[0] = (bi1.num[0] + bi2.num[0]) % 2;
	/*
		 vector<unsigned short> multResult;

		 for(i=0;i<=len;i++)
		  multResult.push_back(sum[i]);
		BigInt result(multResult);
		*/
	delete[]x1;
	delete[]x2;
	return result;
}
BigInt operator /(BigInt bi1, BigInt bi2)
{
	BigInt r;
	r = r.div(bi1, bi2, 'q');
	r.num[0] = (bi1.num[0] + bi2.num[0]) % 2;
	return r;
}
BigInt operator %(BigInt bi1, BigInt bi2)
{
	BigInt r;
	r = r.div(bi1, bi2, 'r');
	r.num[0] = (bi1.num[0] + bi2.num[0]) % 2;
	return r;
}
bool operator == (BigInt& bi1, BigInt bi2)
{
	if (bi1.num[0] != bi2.num[0]) return false;
	if (BigInt::BigIntcmp(bi1, bi2) == 0) return true;
	return false;
}
bool operator >(BigInt& bi1, BigInt bi2)
{
	bool r;
	if (bi1.num[0] > bi2.num[0]) return false;
	if (bi1.num[0] < bi2.num[0]) return true;
	int temp = BigInt::BigIntcmp(bi1, bi2);
	if (temp == 0) return false;
	if (temp == 1)
		r = true;
	else
		r = false;
	if (bi1.num[0] == 1)
		r = !r;
	return r;
}
bool operator <(BigInt& bi1, BigInt bi2)
{
	bool r;
	if (bi1.num[0] > bi2.num[0]) return true;
	if (bi1.num[0] < bi2.num[0]) return false;
	int temp = BigInt::BigIntcmp(bi1, bi2);
	if (temp == 0) return false;
	if (temp == 1)
		r = false;
	else
		r = true;
	if (bi1.num[0] == 1)
		r = !r;
	return r;
}
BigInt BigInt::add(BigInt& bi1, BigInt& bi2)
{
	int i;
	BigInt* mx, * mn;
	int carry = 0;
	if (bi1.num_Size >= bi2.num_Size)
		mx = &bi1, mn = &bi2;
	else
		mx = &bi2, mn = &bi1;
	BigInt bi(0, mx->num_Size + SPACE_CAPACITY);
	for (i = 1; i <= mn->num_Size; i++)
	{
		bi.num[i] = mx->num[i] + mn->num[i] + carry;
		if (bi.num[i] >= num_Weight)
		{
			bi.num[i] -= num_Weight;
			carry = 1;
		}
		else
		{
			carry = 0;
		}
	}
	for (; i <= mx->num_Size; i++)
	{
		bi.num[i] = mx->num[i] + carry;
		if (bi.num[i] >= num_Weight)
		{
			bi.num[i] -= num_Weight;
			carry = 1;
		}
		else
		{
			carry = 0;
		}
	}
	bi.num_Size = mx->num_Size;
	if (carry == 1)
		bi.num[++bi.num_Size] = 1;
	// bi.display();
	return bi;

}
BigInt BigInt::div(BigInt bi1, BigInt bi2, char flag)
{
	BigInt quo(0, bi1.num_Capacity), rem(0, bi2.num_Capacity);
	quo.num_Size = bi1.num_Size - bi2.num_Size + 1;
	quo.num_Capacity = bi1.num_Capacity;
	rem.num_Size = bi2.num_Size;
	rem.num_Capacity = bi2.num_Capacity;
	if (bi2.num_Size == 1)
	{
		if (bi2.num[1] == 0)
		{
			cout << "错误：除数为0！" << endl;
			exit(1);
		}
		T sor = bi2.num[1], i = bi1.num_Size, dend = bi1.num[bi1.num_Size];
		if (dend < sor)
		{
			quo.num_Size--;
			i--;
			dend = dend * num_Weight + bi1.num[i];
		}
		for (; i > 1; i--)
		{
			quo.num[i] = dend / sor;
			dend = dend % sor * num_Weight + bi1.num[i - 1];
		}
		quo.num[1] = dend / sor;
		dend = dend % sor;
		if (flag == 'q')
			return quo;
		else
		{
			rem.num[1] = dend; return rem;
		}


	}
	if (bi2.num_Size > 1)
	{
		int i;
		T d = num_Weight / (bi2.num[bi2.num_Size] + 1);
		BigInt U = bi1 * d, v = bi2 * d;
		U.num[0] = 0;
		v.num[0] = 0;
		BigInt u(0, v.num_Size + SPACE_CAPACITY + 1);
		u.num[0] = 0;
		u.num_Size = v.num_Size;
		for (i = u.num_Size; i >= 1; i--)
		{
			u.num[i] = U.num[U.num_Size - u.num_Size + i];
		}
		T q, r, n = v.num_Size, j;
		T v1 = v.num[v.num_Size], v2 = v.num[v.num_Size - 1], temp;
		//cout<<u<<endl;
		//cout<<v<<endl;
		for (j = U.num_Size - v.num_Size + 1; j >= 1; j--)
		{
			if (u.num_Size < n)
			{
				quo.num[j] = q;
				u = u * num_Weight + U.num[j - 1];
				continue;
			}
			if (u.num_Size == n)
			{
				temp = u.num[u.num_Size];
				q = temp / v1;
				r = temp % v1;
				while (q * v2 > num_Weight * r + u.num[u.num_Size - 1] && r < num_Weight)
				{
					q--;
					r += v1;
				}


			}
			else
			{
				temp = u.num[u.num_Size] * num_Weight + u.num[u.num_Size - 1];
				q = temp / v1;
				r = temp % v1;
				while (q * v2 > num_Weight * r + u.num[u.num_Size - 2] && r < num_Weight)
				{
					q--;
					r += v1;
				}
			}
			u = u - q * v;
			if (u.num[0] == 1)
			{
				u = u + v; q--;
			}
			quo.num[j] = q;
			u = u * num_Weight + U.num[j - 1];
		}


		if (flag == 'q')
		{
			for (i = quo.num_Size; i > 1; i--)
				if (quo.num[i] == 0)
					quo.num_Size--;
				else
					break;
			return quo;
		}
		else
		{
			for (i = u.num_Size; i > 1; i--)
				if (u.num[i] == 0)
					u.num_Size--;
				else
					break;
			u = div(u, num_Weight, 'q');
			u = div(u, d, 'q');
			return u;
		}

	}

}
void BigInt::numstr_convert(string numstr)
{
	long len = numstr.size(), lenMint_s;
	int i, j;
	lenMint_s = len % Int_Size;
	num_Size = len / Int_Size;
	j = num_Size;
	if (lenMint_s != 0)
	{
		num_Size++;
		j = num_Size;

		for (i = 0, num[j] = 0; i < lenMint_s; i++)
		{
			num[j] = num[j] * 10 + numstr[i] - 48;
		}
		j--;
	}
	else
		lenMint_s = Int_Size;
	for (; j >= 1; j--)
	{
		num[j] = 0;
		for (i = 0; i < Int_Size; i++)
		{
			num[j] = num[j] * 10 + numstr[lenMint_s + (num_Size - j - 1) * Int_Size + i] - 48;
		}
	}
}

void BigInt::numstr_test(string numstr)
{//检测字符串是否是合法的整数


}
int BigInt::BigIntcmp(BigInt& bi1, BigInt& bi2)
{//比较绝对值大小，如果bi1>bi2返回1；如果bi1=bi2返回0；如果bi1<bi2返回-1
	if (bi1.num_Size > bi2.num_Size) return 1;
	if (bi1.num_Size < bi2.num_Size) return -1;
	for (int i = bi1.num_Size; i >= 1; i--)
	{
		if (bi1.num[i] > bi2.num[i]) return 1;
		if (bi1.num[i] < bi2.num[i]) return -1;
	}
	return 0;
}
BigInt BigInt::sub(BigInt& bi1, BigInt& bi2)
{//减法函数：实现绝对值大数减小数，末期要重载-运算符
	int i;
	int carry = 0;
	BigInt bi(0, bi1.num_Size + SPACE_CAPACITY);
	for (i = 1; i <= bi2.num_Size; i++)
	{
		bi.num[i] = bi1.num[i] - bi2.num[i] + carry;
		if (bi.num[i] < 0)
		{
			bi.num[i] += num_Weight;
			carry = -1;
		}
		else
		{
			carry = 0;
		}
	}
	for (; i <= bi1.num_Size; i++)
	{
		bi.num[i] = bi1.num[i] + carry;
		if (bi.num[i] < 0)
		{
			bi.num[i] += num_Weight;
			carry = -1;
		}
		else
		{
			carry = 0;
		}
	}

	for (i = bi1.num_Size; bi.num[i] == 0 && i > 1; i--);
	bi.num_Size = i;
	return bi;
}

class BigFloat
{
public:
	BigFloat(string str);
	BigFloat(BigInt bi) { this->bi = bi; this->order = 0; }
	BigFloat(double d = 0) { string str; ostringstream oss; oss << d; str = oss.str(); *this = (BigFloat)str; }
	/*static setprecision(T p) { precision = p; }
	static showprecision() { cout << precision << endl; }*/
	void simple();
	BigFloat operator =(BigFloat bf) { this->bi = bf.bi; this->order = bf.order; return *this; }
	BigFloat operator -() { BigFloat temp = *this; temp.bi = -temp.bi; return temp; }
	BigFloat operator +=(BigFloat bf) { return *this = *this + bf; }
	BigFloat operator -=(BigFloat bf) { return *this = *this - bf; }
	BigFloat operator *=(BigFloat bf) { return *this = *this * bf; }
	BigFloat operator /=(BigFloat bf) { return *this = *this / bf; }
	friend istream& operator >>(istream& input, BigFloat& bf);
	friend ostream& operator <<(ostream& output, BigFloat& bf);
	friend BigFloat operator +(BigFloat bf1, BigFloat bf2);
	friend BigFloat operator -(BigFloat& bf1, BigFloat& bf2);
	friend BigFloat operator *(const BigFloat& bf1, const BigFloat& bf2);
	friend BigFloat operator /(BigFloat bf1, BigFloat bf2);
	friend bool operator ==(BigFloat bf1, BigFloat bf2);
	friend bool operator !=(BigFloat bf1, BigFloat bf2);
	friend bool operator >(BigFloat bf1, BigFloat bf2);
	friend bool operator >=(BigFloat bf1, BigFloat bf2);
	friend bool operator <(BigFloat bf1, BigFloat bf2);
	friend bool operator <=(BigFloat bf1, BigFloat bf2);
private:
	BigInt bi;
	T order;
	static T precision;


};
T BigFloat::precision = 100;
BigFloat::BigFloat(string strS)
{
	int ee = 0;
	string str;
	T loc = strS.find('e', 0);
	if (loc == string::npos)
		loc = strS.find('E', 0);
	if (loc == string::npos)
	{
		str = strS;
	}
	else
	{
		str = strS.substr(0, loc);
		stringstream ss;
		ss << strS.substr(loc + 1, strS.size());//strS.length
		ss >> ee;
	}
	T sign = 0;
	if (str[0] == '-')
	{
		sign = 1;
		str.erase(0, 1);
	}
	loc = str.find(".", 0);
	if (loc == string::npos)
	{
		order = 0;
	}
	else
	{
		order = loc - str.length() + 1;
		str.erase(loc, 1);
		while (str[0] == '0')
			str.erase(0, 1);
		str.append(order % BigInt::Int_Size + BigInt::Int_Size, '0');
		order -= order % BigInt::Int_Size + BigInt::Int_Size;

	}
	if (sign == 1)
		str = '-' + str;
	bi = BigInt(str);
	order += ee;
}
void BigFloat::simple()
{
	int i;
	if (-order > precision)
	{
		if (bi.num_Size <= (-order - precision) / 4)
		{
			bi.num[1] = 0;
			bi.num[0] = 0;
			order = 0;
			return;
		}
		bi = bi >> (-order - precision) / 4;
		order = -precision;
	}

	for (i = 1; i < bi.num_Size; i++)
		if (bi.num[i] != 0)
			break;
		else
			order += 4;
	bi = bi >> (--i);
}
//ostream& operator <<(ostream& output, BigFloat& bf)
//{
//	output = operator <<(output, bf.bi);
//	if (bf.order != 0)
//		output << "e" << bf.order;
//	return output;
//}
istream& operator >>(istream& input, BigFloat& bf)
{
	string st;
	input >> st;
	bf = BigFloat(st);
	return input;
}
BigFloat operator +(BigFloat bf1, BigFloat bf2)
{
	if (bf1.order > bf2.order)
	{
		bf1.bi = bf1.bi << (bf1.order - bf2.order) / BigInt::Int_Size;
		bf1.order = bf2.order;
	}
	if (bf1.order < bf2.order)
	{
		bf2.bi = bf2.bi << (bf2.order - bf1.order) / BigInt::Int_Size;
		bf2.order = bf1.order;
	}
	bf1.bi += bf2.bi;
	bf1.simple();
	return bf1;
}
BigFloat operator -(BigFloat& bf1, BigFloat& bf2)
{

	return bf1 + (-bf2);
}
BigFloat operator *(const BigFloat& bf1, const BigFloat& bf2)
{
	BigFloat r;
	r.bi = bf1.bi * bf2.bi;
	r.order = bf1.order + bf2.order;
	r.simple();
	return r;
}
BigFloat operator /(BigFloat bf1, BigFloat bf2)
{
	if (bf2.order - bf1.order < BigFloat::precision)
	{
		bf1.bi = bf1.bi << (BigFloat::precision + bf1.order - bf2.order) / BigInt::Int_Size;
		bf1.order -= (BigFloat::precision + bf1.order - bf2.order);
	}
	bf1.bi /= bf2.bi;
	bf1.order = bf1.order - bf2.order;
	bf1.simple();
	return bf1;
}
bool operator ==(BigFloat bf1, BigFloat bf2)
{
	if (bf1.order > bf2.order) bf1.bi << (bf1.order - bf2.order) / BigInt::Int_Size;
	if (bf1.order < bf2.order) bf2.bi << (bf2.order - bf2.order) / BigInt::Int_Size;
	return bf1.bi == bf2.bi;
}
bool operator !=(BigFloat bf1, BigFloat bf2)
{
	if (bf1.order > bf2.order) bf1.bi << (bf1.order - bf2.order) / BigInt::Int_Size;
	if (bf1.order < bf2.order) bf2.bi << (bf2.order - bf2.order) / BigInt::Int_Size;
	return bf1.bi != bf2.bi;
}
bool operator >(BigFloat bf1, BigFloat bf2)
{
	if (bf1.order > bf2.order) bf1.bi << (bf1.order - bf2.order) / BigInt::Int_Size;
	if (bf1.order < bf2.order) bf2.bi << (bf2.order - bf2.order) / BigInt::Int_Size;
	return bf1.bi > bf2.bi;
}
bool operator >=(BigFloat bf1, BigFloat bf2)
{
	if (bf1.order > bf2.order) bf1.bi << (bf1.order - bf2.order) / BigInt::Int_Size;
	if (bf1.order < bf2.order) bf2.bi << (bf2.order - bf2.order) / BigInt::Int_Size;
	return bf1.bi >= bf2.bi;
}
bool operator <(BigFloat bf1, BigFloat bf2)
{
	if (bf1.order > bf2.order) bf1.bi << (bf1.order - bf2.order) / BigInt::Int_Size;
	if (bf1.order < bf2.order) bf2.bi << (bf2.order - bf2.order) / BigInt::Int_Size;
	return bf1.bi < bf2.bi;
}
bool operator <=(BigFloat bf1, BigFloat bf2)
{
	if (bf1.order > bf2.order) bf1.bi << (bf1.order - bf2.order) / BigInt::Int_Size;
	if (bf1.order < bf2.order) bf2.bi << (bf2.order - bf2.order) / BigInt::Int_Size;
	return bf1.bi <= bf2.bi;
}
int main()
{
	/*BigFloat x0 = 2, x1, es = 1e-200;
	BigFloat::setprecision(600);
	for (int i = 0; i <= 300; i++)
	{
		x1 = 0.5 * (x0 + 2 / x0);
		if (x0 - x1 < es) break;
		x0 = x1;
	}

	cout << x1 << " " << i;*/

	//1.4142135623730950488016887242096980785696718
	//75376948073176679737990732478462107038850387534327641573
	BigInt a = BigInt("124123671517634124312");
	cout << a << endl;
	return 0;
}