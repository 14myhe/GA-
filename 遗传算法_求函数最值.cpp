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
//��Ӧ�Ժ���

#define pc 0.7    //�������
#define pm 0.001  //ͻ�����
//#define each_d  0.0000007152557373046875  
#define each_d 6.0/(pow(2,23)-1)
#define N 500 //��ʼ��Ⱦɫ����Ⱥ����
#define times 500//�����Ŵ������Ĵ���
#define mutation_num  round(46*0.001*N)


/*
1.��ȷ���⣬����������x1,x2(��һ��46λ����������ʾ,�����ڷָ�Ϊ��������);����һ��46λ��������һ���Ա���
���ݾ���0.000001�������������Ҫ23λ����������ʾ��
2.��������һ����Ⱦɫ���࣬һ����Ⱥ����
chromosomes��Ⱦɫ���ࣩ��Ա������Ⱦɫ��01���顢��Ӧֵ����Ӧ���ʵȵ�
��Ա������
void transform(int* i);//��Ⱦɫ��ָ�Ϊx1,x2
double cal_adapt_value(double x1,double x2); //������Ӧֵ 
Ⱥ������Ҫ�ǰ������Ⱦɫ����� 
��Ա�����У�
void init();//��ʼ�� 
void cal_factor(chromosomes* c);//����Ⱦɫ��Ԫ�صĸ�������ֵ 
int  current_greatest();//ȡ����Ӧ�Ա�����ߵ�Ⱦɫ�� 
int select_two();//���ȡ��һ��Ⱦɫ�� 
int  return_index(double x);//���ظ�Ⱦɫ���Ӧ�ڸ���Ⱥ/����Ⱥ���±� 
void recombine();//ִ�н���
void mutation();//ִ��ͻ��
void create_new();//�γ��µ���Ⱥ 
void refresh();//�ڸ�����Ⱥ��ѡ���������N������ 
void rebuild();//�ظ����ֱ���õ����Ž�
3.�ȳ�ʼ��һ����Ⱥ��������������ֵ��������Ӧֵ����Ӧ�Ա��ʵȵȣ� 
4. ���ѡȡN*pm�����������н��棬���ﵽ��Ӧ�Ա��ʸߵı����еĸ��ʴ�[��Ҫ��ͨ������begin,end������ÿ��������ܱ����е���ֵ����]�� 
5.���ѡȡN*pc*46������λͻ��
6.ѡȡ��������Ⱥ�����ŵ�N�������γ��µ�"����Ⱥ"��
7.һֱ����ֱ����Ӧֵ����Ҫ��

*/


//Ⱦɫ����
class chromosomes{
      public:
        int value_str[46];//Ⱦɫ�崮
        double x1=0.0;
        double x2=0.0;
        double adapt_value=0.0;//��Ӧ��
        double adapt_rate=0.0;//��Ӧ�Ա��ʣ�%��
        double begin,end;
        int no=0;
      public:
        void transform(int* i);//��Ⱦɫ��ָ�Ϊx1,x2
        double cal_adapt_value(double x1,double x2);

};
//�ָ�Ϊ�����֣������x1,x2��ֵ��������[0.6]�� 
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
//������Ӧ�Ժ���
double chromosomes::cal_adapt_value(double x1,double x2)//x1,x2�ǻ���
{
    double res=0.0;
    double j_x1=sin(2*x1);
    double j_x2=sin(2*x2);
    res=j_x1*j_x1+j_x2*j_x2;
    return res;
}

//Ⱦɫ������
class c_group{
  public:
    chromosomes c[N];
    chromosomes newc[N];//��ֳ��
    chromosomes muc[N];
  public:
    void init();//��ʼ�� 
    void cal_factor(chromosomes* c);//
    int  current_greatest();
    int select_two();//ȡ��һ��Ⱦɫ�� 
    int  return_index(double x);
    void recombine();//����
    void mutation();//ͻ��
    void create_new();
    void refresh();
	void rebuild();//�ع�
    void train();
    
}; 

void c_group::init()//����Ⱦɫ��ĳ�ʼ����Ⱥ
{
     //�������� 
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

cout<<"Ⱦɫ�崮��"<<endl;
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
  //����x1��x2ʵ��ֵ���Լ���Ӧ��
  double total=0.0;
  for(int i=0;i<N;i++)
  {
        c[i].transform(c[i].value_str);
        c[i].adapt_value=c[i].cal_adapt_value(c[i].x1,c[i].x2);//����Ⱦɫ����Ӧֵ 
        total+=c[i].adapt_value;
  }
  //������Ӧ�Ա���
  for(int i=0;i<N;i++)
  {
      c[i].adapt_rate=c[i].adapt_value/total;//�����Ӧֵ
      c[i].no=round(c[i].adapt_rate*N);//������� 
  } 
  //������Ӧ�Ա��ʴ󱻳��еļ��ʸ�
  double first=0.0;
  for(int i=0;i<N;i++)
  {
      c[i].begin=first;
      c[i].end=c[i].begin+c[i].adapt_rate;
      first=c[i].end;
  }
}
//�Ӵ� 
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
	//��ֵ 
	for(int i=0;i<index-1;i++)
	{
		for(int j=0;j<46;j++)
		{
			muc[i].value_str[j]=newc[i].value_str[j];
		}
	}
	cal_factor(newc);//���¼��� 
	cal_factor(muc); 

}
int c_group::select_two()
{
  //ѡ��һ��Ⱦɫ��
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
//����
void c_group::recombine()
{
    //����һ����������ж��費��Ҫ����ִ������
  //ѡ������Ⱦɫ��L1,L2
    int l1,l2;
    for(int i=0;i<(N*pc)/2;i++)
    {
		l1=select_two();
        l2=select_two();
     //���棬ʵ�ֽ�����Ҫ�������������һ���ǽ�����λ��length��һ������ʼλ��begin_index;
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
    cal_factor(newc);//���¼��� 
  
}

//ͻ��
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
     cal_factor(muc);//���¼��� 
    
}

int cmp(chromosomes a,chromosomes b)
{
	return a.adapt_rate>b.adapt_rate;
}
void c_group::refresh()
{
	recombine();
	mutation();
//ɸѡ�������������Ӵ����������N�� 
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
    //ȡǰN�� 
     for(int i=0;i<N;i++)
     {
     	for(int j=0;j<46;j++)
     	{
     		c[i].value_str[j]=total[i].value_str[j];
     	
     	}
     }
     cal_factor(c);//���¼��� 

}
int c_group::current_greatest()
{
  //��ǰ��Ӧ�Ա���ߵ�Ⱦɫ���±�
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
           cout<<"����Ҫ���ʱx1="<<c[i].x1<<"  x2="<<c[i].x2<<" sin(2*x1)*sin(2*x1)="<<sin(2*c[i].x1)*sin(2*c[i].x1)<<"  sin(2*x2)*sin(2*x2)="<<sin(2*c[i].x2)*sin(2*c[i].x2)<<endl;
           cout<<"res="<<3-c[i].cal_adapt_value(c[i].x1,c[i].x2)<<endl;
           cout<<"��������Ϊ��"<<count<<endl; 
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
     c_group cn;//��ʼ������Ⱥ 
     cn.init();//��ʼ��
     cn.rebuild();//ѵ��
     system("pause");
}
