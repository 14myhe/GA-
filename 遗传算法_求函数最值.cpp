#include<iostream>
#include<math.h>
#include<stdlib.h>
#include<time.h> 
#include<string.h> 
#include<cstring>
#include<algorithm>
#include<memory.h>
#include<stdlib.h> 
using namespace std;
//适应性函数

#define pc 0.7    //交叉概率
#define pm 0.001  //突变概率
//#define each_d  0.0000007152557373046875  
#define each_d 6.0/(pow(2,23)-1)
#define N 500 //初始化染色体种群数量
#define times 500//设置遗传迭代的次数
#define mutation_num  round(46*0.001*N)


/*
1.明确问题，有两个变量x1,x2(用一个46位的数组来表示,后续在分割为两个数组);所以一个46位的数组是一个自变量
根据精度0.000001，计算出至少需要23位二进制来表示）
2.有两个类一个是染色体类，一个是群体类
chromosomes（染色体类）成员属性有染色体01数组、适应值、适应比率等等
成员方法有
void transform(int* i);//将染色体分割为x1,x2
double cal_adapt_value(double x1,double x2); //计算适应值 
群体类主要是包括多个染色体对象 
成员方法有：
void init();//初始化 
void cal_factor(chromosomes* c);//计算染色体元素的各个属性值 
int  current_greatest();//取出适应性比率最高的染色体 
int select_two();//随机取出一条染色体 
int  return_index(double x);//返回该染色体对应在父种群/子种群的下标 
void recombine();//执行交叉
void mutation();//执行突变
void create_new();//形成新的子群 
void refresh();//在父、子群中选出最优秀的N个个体 
void rebuild();//重复多次直到得到最优解
3.先初始化一个种群，计算个体的属性值（包括适应值、适应性比率等等） 
4. 随机选取N*pm个个体来进行交叉，（达到适应性比率高的被抽中的概率大[主要是通过设置begin,end来表明每个个体可能被抽中的数值概率]） 
5.随机选取N*pc*46个基因位突变
6.选取父、子种群中最优的N个个体形成新的"父种群"；
7.一直迭代直到适应值符合要求

*/


//染色体类
class chromosomes{
      public:
        int value_str[46];//染色体串
        double x1=0.0;
        double x2=0.0;
        double adapt_value=0.0;//适应性
        double adapt_rate=0.0;//适应性比率（%）
        double begin,end;
        int no=0;
      public:
        void transform(int* i);//将染色体分割为x1,x2
        double cal_adapt_value(double x1,double x2);

};
//分割为两部分，计算出x1,x2的值（在区间[0.6]） 
void chromosomes::transform(int* binarynum)
{
   int x=0,y=0;
     for(int i=0;i<23;i++)
     {
        x+=binarynum[i]*pow(2,i);
     }
     for(int i=23;i<46;i++)
     {
        y+=binarynum[i]*pow(2,i-23);
     } 
     x1=x*each_d;
     x2=y*each_d;
}
//计算适应性函数
double chromosomes::cal_adapt_value(double x1,double x2)//x1,x2是弧度
{
    double res=0.0;
    double j_x1=sin(2*x1);
    double j_x2=sin(2*x2);
    res=j_x1*j_x1+j_x2*j_x2;
    return res;
}

//染色体组类
class c_group{
  public:
    chromosomes c[N];
    chromosomes newc[N];//繁殖池
    chromosomes muc[N];
  public:
    void init();//初始化 
    void cal_factor(chromosomes* c);//
    int  current_greatest();
    int select_two();//取出一条染色体 
    int  return_index(double x);
    void recombine();//交叉
    void mutation();//突变
    void create_new();
    void refresh();
	void rebuild();//重构
    void train();
    
}; 

void c_group::init()//建立染色体的初始化种群
{
     //重新启动 
  srand(time(0));
  for(int i=0;i<N;i++)
  {
    for(int j=0;j<46;j++)
    {
       c[i].value_str[j]=rand()%2;
    }
  }
  cal_factor(c);
  create_new();

cout<<"染色体串："<<endl;
   for(int  i=0;i<N;i++)
   {
     for(int j=0;j<46;j++)
      {
           cout<<c[i].value_str[j]<<",";
      } 
      cout<<endl;
   } 
}
void c_group::cal_factor(chromosomes* c)
{
  //计算x1，x2实际值，以及适应性
  double total=0.0;
  for(int i=0;i<N;i++)
  {
        c[i].transform(c[i].value_str);
        c[i].adapt_value=c[i].cal_adapt_value(c[i].x1,c[i].x2);//计算染色体适应值 
        total+=c[i].adapt_value;
  }
  //计算适应性比率
  for(int i=0;i<N;i++)
  {
      c[i].adapt_rate=c[i].adapt_value/total;//相对适应值
      c[i].no=round(c[i].adapt_rate*N);//计算个数 
  } 
  //构造适应性比率大被抽中的几率高
  double first=0.0;
  for(int i=0;i<N;i++)
  {
      c[i].begin=first;
      c[i].end=c[i].begin+c[i].adapt_rate;
      first=c[i].end;
  }
}
//子代 
void c_group::create_new()
{
	int index=0;
	for(int i=0;i<N&&index<200;i++)
	{
		for(int j=0;j<c[i].no&&index<200;j++)
		{
			for(int k=0;k<46;k++)
			{
				newc[index].value_str[k]=c[i].value_str[k];
			}
			index++;
		}
	}
	//赋值 
	for(int i=0;i<index-1;i++)
	{
		for(int j=0;j<46;j++)
		{
			muc[i].value_str[j]=newc[i].value_str[j];
		}
	}
	cal_factor(newc);//重新计算 
	cal_factor(muc); 

}
int c_group::select_two()
{
  //选出一条染色体
     double ran1;
     int rindex;
     ran1=(double)(rand()%10001*0.0001);
     rindex=return_index(ran1);
     return rindex;
}
int c_group::return_index(double x)
{

    for(int i=0;i<N;i++)
    {
        if((x>=c[i].begin)&&(x<c[i].end))
        {
            return i;
        }
    }
    return 0; 
}
//交叉
void c_group::recombine()
{
    //产生一个随机数来判断需不需要真正执行重组
  //选中两条染色体L1,L2
    int l1,l2;
    for(int i=0;i<(N*pc)/2;i++)
    {
		l1=select_two();
        l2=select_two();
     //交叉，实现交叉需要有两个随机数，一个是交换的位数length，一个是起始位置begin_index;
        int begin_index=rand()%46;
        int length=rand()%(46-begin_index);
        int temp;
        for(int i=0;i<=length;i++)
        {
            temp=newc[l1].value_str[begin_index+i];
		    newc[l1].value_str[begin_index+i]=newc[l2].value_str[begin_index+i];
            newc[l2].value_str[begin_index+i]=temp;
        } 
    }
    cal_factor(newc);//重新计算 
  
}

//突变
void c_group::mutation()
{
    
    for(int k=0;k<mutation_num;k++)
    {
        int i=select_two();
        int m_index1=rand()%46;
          
        if(muc[i].value_str[m_index1]==0)
        {
            muc[i].value_str[m_index1]=1;
        }
        else
        {
            muc[i].value_str[m_index1]=0;
        }
    }
     cal_factor(muc);//重新计算 
    
}

int cmp(chromosomes a,chromosomes b)
{
	return a.adapt_rate>b.adapt_rate;
}
void c_group::refresh()
{
	recombine();
	mutation();
//筛选出父代和所有子代中最优秀的N个 
     int sum=3*N;
     chromosomes total[sum];
     for(int i=0;i<N;i++)
     {	
	    total[i].adapt_rate=c[i].adapt_rate; 
     	for(int j=0;j<46;j++)
     	{
     		total[i].value_str[j]=c[i].value_str[j];
     	
     	}
     }
     for(int i=0;i<N;i++)
     {
     	total[i+N].adapt_rate=newc[i].adapt_rate; 
     	for(int j=0;j<46;j++)
     	{
     		total[i+N].value_str[j]=newc[i].value_str[j];
     	}
     }
     int a=2*N; 
     for(int i=0;i<N;i++)
     {
     	total[i+a].adapt_rate=muc[i].adapt_rate; 
     	for(int j=0;j<46;j++)
     	{
     		total[i+a].value_str[j]=muc[i].value_str[j];
     	}
     }
     sort(total,total+sum,cmp);
    //取前N个 
     for(int i=0;i<N;i++)
     {
     	for(int j=0;j<46;j++)
     	{
     		c[i].value_str[j]=total[i].value_str[j];
     	
     	}
     }
     cal_factor(c);//重新计算 

}
int c_group::current_greatest()
{
  //当前适应性比最高的染色体下标
   double r=c[0].adapt_rate;
   int greatest_index=0; 
   for(int i=1;i<N;i++)
   { 
       if(c[i].adapt_rate>r)
       {
           greatest_index=i;
       }
   }
   return greatest_index;
}

void c_group::rebuild()
{
	
    int count=0;
    int flag=0;
    while(count<times)
    {
        refresh();
        int i=current_greatest();
        if(abs(2-c[i].cal_adapt_value(c[i].x1,c[i].x2))<=0.000001)
        {
           flag=1;
           cout<<"符合要求此时x1="<<c[i].x1<<"  x2="<<c[i].x2<<" sin(2*x1)*sin(2*x1)="<<sin(2*c[i].x1)*sin(2*c[i].x1)<<"  sin(2*x2)*sin(2*x2)="<<sin(2*c[i].x2)*sin(2*c[i].x2)<<endl;
           cout<<"res="<<3-c[i].cal_adapt_value(c[i].x1,c[i].x2)<<endl;
           cout<<"迭代次数为："<<count<<endl; 
        }
    count++;
    }
   if(flag==0)
   {
   	   init();
       rebuild();
   } 

}

void c_group::train()
{
   init();
   rebuild();
     
}

 
int main()
{
     c_group cn;//初始化的总群 
     cn.init();//初始化
     cn.rebuild();//训练
     system("pause");
}
