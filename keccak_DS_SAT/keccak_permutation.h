#include<bits/stdc++.h>
#include "table.h"
using namespace std;

int mod(int a,int b)
{
    while(a<0)
    {
        a+=b;
    }
    return a%b;
}

void bin_to_intvector(int x,int length,int* result)
{
    for(int i=0;i<length;i++)
    {
        uint8_t b=x/int(pow(2,length-1));
        result[i]=b;
        x=(x<<1)%int(pow(2,length));
    }
}

void bin_to_vector(int x,int length,uint8_t* result)
{
    for(int i=0;i<length;i++)
    {
        uint8_t b=x/int(pow(2,length-1));
        result[i]=b;
        x=(x<<1)%int(pow(2,length));
    }
}

void bin64_to_vector(uint64_t x,uint8_t* result)
{
    for(int i=0;i<64;i++)
    {
        uint8_t b=x/uint64_t(pow(2,64-1));
        result[i]=b;
        x=(x<<1);
    }
}

int comp_vector(uint8_t* a,uint8_t* b,int len)
{
    for(int i=0;i<len;i++)
    {
        if(a[i]!=b[i])
        {
            return 0;
        }
    }
    return 1;
}

int get_locate(int x,int y,int z)
{
    return 64*(5*y+x)+z;
}

int pusai0(int i)
{
    return mod((i/64),5);
}

int pusai1(int i)
{
    return i/320;
}

int pusai2(int i)
{
    return mod(i,64);
}

int phi0(int i)
{
    return 64*pusai0(i)+pusai2(i);
}

int phi1(int i)
{
    return 64*(mod(pusai0(i)-1,5))+pusai2(i);
}

int phi2(int i)
{
    return 64*(mod(pusai0(i)+1,5))+mod(pusai2(i)-1,64);
}

void theta(uint8_t* a,uint8_t* b)
{
    for(int x=0;x<5;x++)
    {
        for(int y=0;y<5;y++)
        {
            for(int z=0;z<64;z++)
            {
                uint8_t _sum=0;
                for(int y_p=0;y_p<5;y_p++)
                {
                    _sum=_sum^a[get_locate(mod(x-1,5),y_p,z)]^a[get_locate(mod(x+1,5),y_p,mod(z-1,64))];
                }
                b[get_locate(x,y,z)]=a[get_locate(x,y,z)]^_sum;
            }
        }
    }
}

void rho(uint8_t* a,uint8_t* b)
{
    for(int x=0;x<5;x++)
    {
        for(int y=0;y<5;y++)
        {
            for(int z=0;z<64;z++)
            {
                b[get_locate(x,y,z)]=a[get_locate(x,y,mod(z-t_step[x][y],64))];
            }
        }
    }
}

void pi(uint8_t*a,uint8_t* b)
{
    for(int x=0;x<5;x++)
    {
        for(int y=0;y<5;y++)
        {
            for(int z=0;z<64;z++)
            {
                b[get_locate(x,y,z)]=a[get_locate(pi_step[x][y][0],pi_step[x][y][1],z)];
            }
        }
    }
}

void chi(uint8_t* a,uint8_t* b)
{
    for(int x=0;x<5;x++)
    {
        for(int y=0;y<5;y++)
        {
            for(int z=0;z<64;z++)
            {
                b[get_locate(x,y,z)]=a[get_locate(x,y,z)]^(a[get_locate(mod(x+1,5),y,z)]&a[get_locate(mod(x+2,5),y,z)])^a[get_locate(mod(x+2,5),y,z)];
            }
        }
    }
}

void iota(uint8_t* a,uint8_t* b,int ir)
{
    memcpy(b,a,1600);
    for(int z=0;z<64;z++)
    {
        b[get_locate(0,0,z)]=a[get_locate(0,0,z)]^RC[ir][z];
    }
}

int sigma(int i)
{
    int x=pusai0(i);
    int y=pusai1(i);
    int z=pusai2(i);
    int z_l=mod(z+t_step[x][y],64);
    int x_l=pi_inverse[x][y][0];
    int y_l=pi_inverse[x][y][1];
    return get_locate(x_l,y_l,z_l);
}

int sigma_inverse(int i)
{
    int x=pusai0(i);
    int y=pusai1(i);
    int z=pusai2(i);
    int x_l=pi_step[x][y][0];
    int y_l=pi_step[x][y][1];
    int z_l=mod(z-t_step[x_l][y_l],64);
    return get_locate(x_l,y_l,z_l);
}

void f_4(uint8_t* a,uint8_t* b)
{
    uint8_t tmp1[1600]={0};
    uint8_t tmp2[1600]={0};
    memcpy(tmp1,a,1600);
    for(int i=0;i<4;i++)
    {
        theta(tmp1,tmp2);
        rho(tmp2,tmp1);
        pi(tmp1,tmp2);
        chi(tmp2,tmp1);
        iota(tmp1,tmp2,i);
        memcpy(tmp1,tmp2,1600);
    }
    memcpy(b,tmp2,1600);
}

set<int> get_S0()
{
    set<int> S0;
    set<int> S1;
    for(int i=0;i<1600;i++)
    {
        int x=pusai0(i);
        int y=pusai1(i);
        int z=pusai2(i);
        int deta_out=0;
        for(int x_=0;x_<5;x_++)
        {
            deta_out+=((a1[get_locate(x_,y,z)])<<(4-x_));
        }
        if(deta_out==0)
        {
            S0.insert(i);
        }
    }
    return S0;
}

set<int> get_S1()
{
    set<int> S0;
    set<int> S1;
    for(int i=0;i<1600;i++)
    {
        int x=pusai0(i);
        int y=pusai1(i);
        int z=pusai2(i);
        int deta_out=0;
        for(int x_=0;x_<5;x_++)
        {
            deta_out+=((a1[get_locate(x_,y,z)])<<(4-x_));
        }
        if((deta_out==1 || deta_out==2 || deta_out==4 || deta_out==8 || deta_out==16) && a1[i]==1)
        {
            S1.insert(i);
        }
    }
    return S1;
}

uint8_t discard_by_more_condition(uint8_t* x)
{
    uint8_t condition1=(x[1150]^x[830])&(x[1563]^x[923]^1);
    uint8_t condition2=(x[1149]^x[829]^1)&(x[1491]^x[1171]^1);
    uint8_t condition3=(x[1151]^x[831])&(x[1564]^x[924]^1);
    uint8_t condition4=(x[1487]^x[847]^1)&(x[1248]^x[928]^1);
    uint8_t condition5=(x[1490]^x[850])&(x[1251]^x[931]^1)&(x[1305]^x[985]);
    uint8_t condition6=(x[1490]^x[850]^1)&(x[1251]^x[931]^1)&(x[1305]^x[985]^1);
    uint8_t condition7=(x[1253]^x[933]^1)&(x[1307]^x[987]^1);
    uint8_t condition8=(x[1493]^x[853])&(x[1254]^x[1574]^1)&(x[1308]^x[988]^1);
    uint8_t condition9=(x[1494]^x[854])&(x[1575]^x[1255]^1)&(x[1309]^x[989]^1);
    uint8_t condition10=(x[1491]^x[851])&(x[1306]^x[986]^1);
    return condition1|condition2|condition3|condition4|condition5|condition6|condition7|condition8|condition9|condition10;
}

void INITIAL_41bit(uint8_t* M1,uint8_t* M1_p,set<int> S0,set<int> S1,
uint8_t* A,uint8_t* A_p,uint8_t* As,uint8_t* As_p,uint8_t* B,uint8_t* B_p,uint8_t* Bs,
uint8_t* Bs_p,uint8_t* a0,uint8_t* b0,uint8_t* a0_s,uint8_t* b0_s,uint8_t* Sigma,uint8_t* Sigma_s,int* MarkedBit)
{
    f_4(M1,A);f_4(M1_p,A_p);
    //modify message
    for(int i=0;i<mat_rank;i++)
    {
        int index[4]={-1,-1,-1,-1};
        int index_cnt=0;
        for(int j=0;j<1600-828;j++)
        {
            if(Mat2[i][j]==1)
            {
                index[index_cnt]=828+j;
                index_cnt++;
            }
        }
        if(index_cnt==2)
        {
            A_p[index[0]]=A[index[0]]^A[index[1]]^A_p[index[1]]^Mat2[i][1600-828];
        }
        else
        {
            A_p[index[0]]=A[index[0]]^A[index[1]]^A_p[index[1]]^A[index[2]]^A_p[index[2]]^A[index[3]]^A_p[index[3]]^Mat2[i][1600-828];
        }
    }
    for(int i=0;i<1600;i++)
    {
        if(i>=828)
        {
            a0[i]=A[i]^A_p[i];
            a0_s[i]=1;
            As[i]=1;
            As_p[i]=1;
        }
        else
        {
            a0[i]=0;
            a0_s[i]=0;
            A[i]=0;
            A_p[i]=0;
            As[i]=0;
            As_p[i]=0;
        }
        auto it0=S0.find(i);
        auto it1=S1.find(i);
        if(it0!=S0.end())
        {
            b0[i]=0;
            b0_s[i]=1;
        }
        else if(it1!=S1.end())
        {
            b0[i]=1;
            b0_s[i]=1;
        }
        else
        {
            b0_s[i]=0;
        }
        B[i]=0;
        B_p[i]=0;
        Bs[i]=0;
        Bs_p[i]=0;
    }
    for(int i=0;i<320;i++)
    {
        Sigma[i]=0;
        Sigma_s[i]=0;
        MarkedBit[i]=-1;
    }
}

void INITIAL_39bit(uint8_t* M1,uint8_t* M1_p,set<int> S0,set<int> S1,
uint8_t* A,uint8_t* A_p,uint8_t* As,uint8_t* As_p,uint8_t* B,uint8_t* B_p,uint8_t* Bs,
uint8_t* Bs_p,uint8_t* a0,uint8_t* b0,uint8_t* a0_s,uint8_t* b0_s,uint8_t* Sigma,uint8_t* Sigma_s,int* MarkedBit)
{
    f_4(M1,A);f_4(M1_p,A_p);
    //modify message
    for(int i=0;i<mat3_rank;i++)
    {
        int index[4]={-1,-1,-1,-1};
        int index_cnt=0;
        for(int j=0;j<1600-828;j++)
        {
            if(Mat3[i][j]==1)
            {
                index[index_cnt]=828+j;
                index_cnt++;
            }
        }
        if(index_cnt==2)
        {
            A_p[index[0]]=A[index[0]]^A[index[1]]^A_p[index[1]]^Mat3[i][1600-828];
        }
        else
        {
            A_p[index[0]]=A[index[0]]^A[index[1]]^A_p[index[1]]^A[index[2]]^A_p[index[2]]^A[index[3]]^A_p[index[3]]^Mat3[i][1600-828];
        }
    }
    //modify done
    for(int i=0;i<1600;i++)
    {
        if(i>=828)
        {
            a0[i]=A[i]^A_p[i];
            a0_s[i]=1;
            As[i]=1;
            As_p[i]=1;
        }
        else
        {
            a0[i]=0;
            a0_s[i]=0;
            A[i]=0;
            A_p[i]=0;
            As[i]=0;
            As_p[i]=0;
        }
        auto it0=S0.find(i);
        auto it1=S1.find(i);
        if(it0!=S0.end())
        {
            b0[i]=0;
            b0_s[i]=1;
        }
        else if(it1!=S1.end())
        {
            b0[i]=1;
            b0_s[i]=1;
        }
        else
        {
            b0_s[i]=0;
        }
        B[i]=0;
        B_p[i]=0;
        Bs[i]=0;
        Bs_p[i]=0;
    }
    for(int i=0;i<320;i++)
    {
        Sigma[i]=0;
        Sigma_s[i]=0;
        MarkedBit[i]=-1;
    }
}

uint8_t INITIAL_discard(uint8_t* M1,uint8_t* M1_p,set<int> S0,set<int> S1,
uint8_t* A,uint8_t* A_p,uint8_t* As,uint8_t* As_p,uint8_t* B,uint8_t* B_p,uint8_t* Bs,
uint8_t* Bs_p,uint8_t* a0,uint8_t* b0,uint8_t* a0_s,uint8_t* b0_s,uint8_t* Sigma,uint8_t* Sigma_s,int* MarkedBit)
{
    f_4(M1,A);f_4(M1_p,A_p);
    //modify message
    for(int i=0;i<mat3_rank;i++)
    {
        int index[4]={-1,-1,-1,-1};
        int index_cnt=0;
        for(int j=0;j<1600-828;j++)
        {
            if(Mat3[i][j]==1)
            {
                index[index_cnt]=828+j;
                index_cnt++;
            }
        }
        if(index_cnt==2)
        {
            A_p[index[0]]=A[index[0]]^A[index[1]]^A_p[index[1]]^Mat3[i][1600-828];
        }
        else
        {
            A_p[index[0]]=A[index[0]]^A[index[1]]^A_p[index[1]]^A[index[2]]^A_p[index[2]]^A[index[3]]^A_p[index[3]]^Mat3[i][1600-828];
        }
    }
    //modify done
    A[829]^=1;A[830]^=1;A[831]^=1;A_p[829]^=1;A_p[830]^=1;A_p[831]^=1;
    for(int i=0;i<1600;i++)
    {
        if(i>=828)
        {
            a0[i]=A[i]^A_p[i];
            a0_s[i]=1;
            As[i]=1;
            As_p[i]=1;
        }
        else
        {
            a0[i]=0;
            a0_s[i]=0;
            A[i]=0;
            A_p[i]=0;
            As[i]=0;
            As_p[i]=0;
        }
    }
    uint8_t flag=discard_by_more_condition(a0);
    if(flag==1)
    {
        return flag;
    }
    for(int i=0;i<1600;i++)
    {
        auto it0=S0.find(i);
        auto it1=S1.find(i);
        if(it0!=S0.end())
        {
            b0[i]=0;
            b0_s[i]=1;
        }
        else if(it1!=S1.end())
        {
            b0[i]=1;
            b0_s[i]=1;
        }
        else
        {
            b0_s[i]=0;
        }
        B[i]=0;
        B_p[i]=0;
        Bs[i]=0;
        Bs_p[i]=0;
    }
    for(int i=0;i<320;i++)
    {
        Sigma[i]=0;
        Sigma_s[i]=0;
        MarkedBit[i]=-1;
    }
    return flag;
}

uint8_t INITIAL_inversestate(uint8_t* M1,uint8_t* M1_p,set<int> S0,set<int> S1,
uint8_t* A,uint8_t* A_p,uint8_t* As,uint8_t* As_p,uint8_t* B,uint8_t* B_p,uint8_t* Bs,
uint8_t* Bs_p,uint8_t* a0,uint8_t* b0,uint8_t* a0_s,uint8_t* b0_s,uint8_t* Sigma,uint8_t* Sigma_s,int* MarkedBit)
{
    memcpy(A,M1,1600);memcpy(A_p,M1_p,1600);
    for(int i=0;i<1600;i++)
    {
        if(i>=828)
        {
            a0[i]=A[i]^A_p[i];
            a0_s[i]=1;
            As[i]=1;
            As_p[i]=1;
        }
        else
        {
            a0[i]=0;
            a0_s[i]=0;
            A[i]=0;
            A_p[i]=0;
            As[i]=0;
            As_p[i]=1;
        }
    }
    uint8_t flag=discard_by_more_condition(a0);
    if(flag==1)
    {
        return flag;
    }
    for(int i=0;i<1600;i++)
    {
        auto it0=S0.find(i);
        auto it1=S1.find(i);
        if(it0!=S0.end())
        {
            b0[i]=0;
            b0_s[i]=1;
        }
        else if(it1!=S1.end())
        {
            b0[i]=1;
            b0_s[i]=1;
        }
        else
        {
            b0_s[i]=0;
        }
        B[i]=0;
        B_p[i]=0;
        Bs[i]=0;
        Bs_p[i]=0;
    }
    for(int i=0;i<320;i++)
    {
        Sigma[i]=0;
        Sigma_s[i]=0;
        MarkedBit[i]=-1;
    }
    return flag;
}

int CPKERNEL(uint8_t* a0,uint8_t* b0,uint8_t* a0_s,uint8_t* b0_s,int column_index)
{
    int i[5]={0};
    int x=column_index/64;
    int z=mod(column_index,64);
    i[0]=get_locate(x,0,z);i[1]=get_locate(x,1,z);i[2]=get_locate(x,2,z);i[3]=get_locate(x,3,z);i[4]=get_locate(x,4,z);
    uint8_t flag=0;
    uint8_t SUM=0;
    int bit_cnt=0;
    for(int j=0;j<5;j++)
    {
        if(a0_s[i[j]]==1 and b0_s[sigma(i[j])]==1)
        {
            flag=1;
            SUM=a0[i[j]]^b0[sigma(i[j])];
            break;
        }
    }
    if(flag==1)
    {
        for(int j=0;j<5;j++)
        {
            if(a0_s[i[j]]==1 and b0_s[sigma(i[j])]==0)
            {
                b0[sigma(i[j])]=SUM^a0[i[j]];
                b0_s[sigma(i[j])]=1;
                bit_cnt++;
            }
            else if(a0_s[i[j]]==0 and b0_s[sigma(i[j])]==1)
            {
                a0[i[j]]=SUM^b0[sigma(i[j])];
                a0_s[i[j]]=1;
                bit_cnt++;
            }
        }
    }
    return bit_cnt;
}

void judge_state(uint8_t* a0,uint8_t* b0,uint8_t* a0_s,uint8_t* b0_s,int column_index,uint8_t &column_sum,int &state,int &marked_bit)
{
    int i[5]={0};
    int j[5]={0};
    int x=column_index/64;
    int z=mod(column_index,64);
    i[0]=get_locate(x,0,z);i[1]=get_locate(x,1,z);i[2]=get_locate(x,2,z);i[3]=get_locate(x,3,z);i[4]=get_locate(x,4,z);
    j[0]=sigma(i[0]);j[1]=sigma(i[1]);j[2]=sigma(i[2]);j[3]=sigma(i[3]);j[4]=sigma(i[4]);
    if( a0_s[i[0]]==1 and a0_s[i[1]]==1 and a0_s[i[2]]==1 and a0_s[i[3]]==1 and a0_s[i[4]]==1)
        {column_sum= a0[i[0]]^a0[i[1]]^a0[i[2]]^a0[i[3]]^a0[i[4]];state=1;marked_bit=-1;return;}
    if( b0_s[j[0]]==1 and b0_s[j[1]]==1 and a0_s[i[2]]==1 and a0_s[i[3]]==1 and a0_s[i[4]]==1)
        {column_sum= b0[j[0]]^b0[j[1]]^a0[i[2]]^a0[i[3]]^a0[i[4]];state=1;marked_bit=-1;return;}
    if( b0_s[j[0]]==1 and a0_s[i[1]]==1 and b0_s[j[2]]==1 and a0_s[i[3]]==1 and a0_s[i[4]]==1)
        {column_sum= b0[j[0]]^a0[i[1]]^b0[j[2]]^a0[i[3]]^a0[i[4]];state=1;marked_bit=-1;return;}
    if( b0_s[j[0]]==1 and a0_s[i[1]]==1 and a0_s[i[2]]==1 and b0_s[j[3]]==1 and a0_s[i[4]]==1)
        {column_sum= b0[j[0]]^a0[i[1]]^a0[i[2]]^b0[j[3]]^a0[i[4]];state=1;marked_bit=-1;return;}
    if( b0_s[j[0]]==1 and a0_s[i[1]]==1 and a0_s[i[2]]==1 and a0_s[i[3]]==1 and b0_s[j[4]]==1)
        {column_sum= b0[j[0]]^a0[i[1]]^a0[i[2]]^a0[i[3]]^b0[j[4]];state=1;marked_bit=-1;return;}
    if( a0_s[i[0]]==1 and b0_s[j[1]]==1 and b0_s[j[2]]==1 and a0_s[i[3]]==1 and a0_s[i[4]]==1)
        {column_sum= a0[i[0]]^b0[j[1]]^b0[j[2]]^a0[i[3]]^a0[i[4]];state=1;marked_bit=-1;return;}
    if( a0_s[i[0]]==1 and b0_s[j[1]]==1 and a0_s[i[2]]==1 and b0_s[j[3]]==1 and a0_s[i[4]]==1)
        {column_sum= a0[i[0]]^b0[j[1]]^a0[i[2]]^b0[j[3]]^a0[i[4]];state=1;marked_bit=-1;return;}
    if( a0_s[i[0]]==1 and b0_s[j[1]]==1 and a0_s[i[2]]==1 and a0_s[i[3]]==1 and b0_s[j[4]]==1)
        {column_sum= a0[i[0]]^b0[j[1]]^a0[i[2]]^a0[i[3]]^b0[j[4]];state=1;marked_bit=-1;return;}
    if( a0_s[i[0]]==1 and a0_s[i[1]]==1 and b0_s[j[2]]==1 and b0_s[j[3]]==1 and a0_s[i[4]]==1)
        {column_sum= a0[i[0]]^a0[i[1]]^b0[j[2]]^b0[j[3]]^a0[i[4]];state=1;marked_bit=-1;return;}
    if( a0_s[i[0]]==1 and a0_s[i[1]]==1 and b0_s[j[2]]==1 and a0_s[i[3]]==1 and b0_s[j[4]]==1)
        {column_sum= a0[i[0]]^a0[i[1]]^b0[j[2]]^a0[i[3]]^b0[j[4]];state=1;marked_bit=-1;return;}
    if( a0_s[i[0]]==1 and a0_s[i[1]]==1 and a0_s[i[2]]==1 and b0_s[j[3]]==1 and b0_s[j[4]]==1)
        {column_sum= a0[i[0]]^a0[i[1]]^a0[i[2]]^b0[j[3]]^b0[j[4]];state=1;marked_bit=-1;return;}
    //state 2
    //左侧4个1
    if( a0_s[i[0]]==1 and a0_s[i[1]]==1 and a0_s[i[2]]==1 and a0_s[i[3]]==1)
        {column_sum= a0[i[0]]^a0[i[1]]^a0[i[2]]^a0[i[3]];state=2;marked_bit=i[4];return;}
    if( a0_s[i[0]]==1 and a0_s[i[1]]==1 and a0_s[i[2]]==1 and a0_s[i[4]]==1)
        {column_sum= a0[i[0]]^a0[i[1]]^a0[i[2]]^a0[i[4]];state=2;marked_bit=i[3];return;}
    if( a0_s[i[0]]==1 and a0_s[i[1]]==1 and a0_s[i[3]]==1 and a0_s[i[4]]==1)
        {column_sum= a0[i[0]]^a0[i[1]]^a0[i[3]]^a0[i[4]];state=2;marked_bit=i[2];return;}
    if( a0_s[i[0]]==1 and a0_s[i[2]]==1 and a0_s[i[3]]==1 and a0_s[i[4]]==1)
        {column_sum= a0[i[0]]^a0[i[2]]^a0[i[3]]^a0[i[4]];state=2;marked_bit=i[1];return;}
    if( a0_s[i[1]]==1 and a0_s[i[2]]==1 and a0_s[i[3]]==1 and a0_s[i[4]]==1)
        {column_sum= a0[i[1]]^a0[i[2]]^a0[i[3]]^a0[i[4]];state=2;marked_bit=i[0];return;}//这里做了修改
    //右侧4个1
    if( b0_s[j[0]]==1 and b0_s[j[1]]==1 and b0_s[j[2]]==1 and b0_s[j[3]]==1)
        {column_sum= b0[j[0]]^b0[j[1]]^b0[j[2]]^b0[j[3]];state=2;marked_bit=i[4];return;}
    if( b0_s[j[0]]==1 and b0_s[j[1]]==1 and b0_s[j[2]]==1 and b0_s[j[4]]==1)
        {column_sum= b0[j[0]]^b0[j[1]]^b0[j[2]]^b0[j[4]];state=2;marked_bit=i[3];return;}
    if( b0_s[j[0]]==1 and b0_s[j[1]]==1 and b0_s[j[4]]==1 and b0_s[j[3]]==1)
        {column_sum= b0[j[0]]^b0[j[1]]^b0[j[4]]^b0[j[3]];state=2;marked_bit=i[2];return;}
    if( b0_s[j[0]]==1 and b0_s[j[4]]==1 and b0_s[j[2]]==1 and b0_s[j[3]]==1)
        {column_sum= b0[j[0]]^b0[j[4]]^b0[j[2]]^b0[j[3]];state=2;marked_bit=i[1];return;}
    if( b0_s[j[4]]==1 and b0_s[j[1]]==1 and b0_s[j[2]]==1 and b0_s[j[3]]==1)
        {column_sum= b0[j[4]]^b0[j[1]]^b0[j[2]]^b0[j[3]];state=2;marked_bit=i[0];return;}
    //左2右2
    if( a0_s[i[0]]==1 and a0_s[i[1]]==1 and b0_s[j[2]]==1 and b0_s[j[3]]==1)
        {column_sum= a0[i[0]]^a0[i[1]]^b0[j[2]]^b0[j[3]];state=2;marked_bit=i[4];return;}
    if( a0_s[i[0]]==1 and a0_s[i[2]]==1 and b0_s[j[1]]==1 and b0_s[j[3]]==1)
        {column_sum= a0[i[0]]^a0[i[2]]^b0[j[1]]^b0[j[3]];state=2;marked_bit=i[4];return;}
    if( a0_s[i[0]]==1 and a0_s[i[3]]==1 and b0_s[j[1]]==1 and b0_s[j[2]]==1)
        {column_sum= a0[i[0]]^a0[i[3]]^b0[j[1]]^b0[j[2]];state=2;marked_bit=i[4];return;}
    if( a0_s[i[1]]==1 and a0_s[i[2]]==1 and b0_s[j[0]]==1 and b0_s[j[3]]==1)
        {column_sum= a0[i[1]]^a0[i[2]]^b0[j[0]]^b0[j[3]];state=2;marked_bit=i[4];return;}
    if( a0_s[i[1]]==1 and a0_s[i[3]]==1 and b0_s[j[0]]==1 and b0_s[j[2]]==1)
        {column_sum= a0[i[1]]^a0[i[3]]^b0[j[0]]^b0[j[2]];state=2;marked_bit=i[4];return;}
    if( a0_s[i[2]]==1 and a0_s[i[3]]==1 and b0_s[j[0]]==1 and b0_s[j[1]]==1)
        {column_sum= a0[i[2]]^a0[i[3]]^b0[j[0]]^b0[j[1]];state=2;marked_bit=i[4];return;}
    
    if( a0_s[i[0]]==1 and a0_s[i[1]]==1 and b0_s[j[2]]==1 and b0_s[j[4]]==1)
        {column_sum= a0[i[0]]^a0[i[1]]^b0[j[2]]^b0[j[4]];state=2;marked_bit=i[3];return;}
    if( a0_s[i[0]]==1 and a0_s[i[2]]==1 and b0_s[j[1]]==1 and b0_s[j[4]]==1)
        {column_sum= a0[i[0]]^a0[i[2]]^b0[j[1]]^b0[j[4]];state=2;marked_bit=i[3];return;}
    if( a0_s[i[0]]==1 and a0_s[i[4]]==1 and b0_s[j[1]]==1 and b0_s[j[2]]==1)
        {column_sum= a0[i[0]]^a0[i[4]]^b0[j[1]]^b0[j[2]];state=2;marked_bit=i[3];return;}
    if( a0_s[i[1]]==1 and a0_s[i[2]]==1 and b0_s[j[0]]==1 and b0_s[j[4]]==1)
        {column_sum= a0[i[1]]^a0[i[2]]^b0[j[0]]^b0[j[4]];state=2;marked_bit=i[3];return;}
    if( a0_s[i[1]]==1 and a0_s[i[4]]==1 and b0_s[j[0]]==1 and b0_s[j[2]]==1)
        {column_sum= a0[i[1]]^a0[i[4]]^b0[j[0]]^b0[j[2]];state=2;marked_bit=i[3];return;}
    if( a0_s[i[2]]==1 and a0_s[i[4]]==1 and b0_s[j[0]]==1 and b0_s[j[1]]==1)
        {column_sum= a0[i[2]]^a0[i[4]]^b0[j[0]]^b0[j[1]];state=2;marked_bit=i[3];return;}
    
    if( a0_s[i[0]]==1 and a0_s[i[1]]==1 and b0_s[j[4]]==1 and b0_s[j[3]]==1)
        {column_sum= a0[i[0]]^a0[i[1]]^b0[j[4]]^b0[j[3]];state=2;marked_bit=i[2];return;}
    if( a0_s[i[0]]==1 and a0_s[i[4]]==1 and b0_s[j[1]]==1 and b0_s[j[3]]==1)
        {column_sum= a0[i[0]]^a0[i[4]]^b0[j[1]]^b0[j[3]];state=2;marked_bit=i[2];return;}
    if( a0_s[i[0]]==1 and a0_s[i[3]]==1 and b0_s[j[1]]==1 and b0_s[j[4]]==1)
        {column_sum= a0[i[0]]^a0[i[3]]^b0[j[1]]^b0[j[4]];state=2;marked_bit=i[2];return;}
    if( a0_s[i[1]]==1 and a0_s[i[4]]==1 and b0_s[j[0]]==1 and b0_s[j[3]]==1)
        {column_sum= a0[i[1]]^a0[i[4]]^b0[j[0]]^b0[j[3]];state=2;marked_bit=i[2];return;}
    if( a0_s[i[1]]==1 and a0_s[i[3]]==1 and b0_s[j[0]]==1 and b0_s[j[4]]==1)
        {column_sum= a0[i[1]]^a0[i[3]]^b0[j[0]]^b0[j[4]];state=2;marked_bit=i[2];return;}
    if( a0_s[i[4]]==1 and a0_s[i[3]]==1 and b0_s[j[0]]==1 and b0_s[j[1]]==1)
        {column_sum= a0[i[4]]^a0[i[3]]^b0[j[0]]^b0[j[1]];state=2;marked_bit=i[2];return;}
    
    if( a0_s[i[0]]==1 and a0_s[i[4]]==1 and b0_s[j[2]]==1 and b0_s[j[3]]==1)
        {column_sum= a0[i[0]]^a0[i[4]]^b0[j[2]]^b0[j[3]];state=2;marked_bit=i[1];return;}
    if( a0_s[i[0]]==1 and a0_s[i[2]]==1 and b0_s[j[4]]==1 and b0_s[j[3]]==1)
        {column_sum= a0[i[0]]^a0[i[2]]^b0[j[4]]^b0[j[3]];state=2;marked_bit=i[1];return;}
    if( a0_s[i[0]]==1 and a0_s[i[3]]==1 and b0_s[j[4]]==1 and b0_s[j[2]]==1)
        {column_sum= a0[i[0]]^a0[i[3]]^b0[j[4]]^b0[j[2]];state=2;marked_bit=i[1];return;}
    if( a0_s[i[4]]==1 and a0_s[i[2]]==1 and b0_s[j[0]]==1 and b0_s[j[3]]==1)
        {column_sum= a0[i[4]]^a0[i[2]]^b0[j[0]]^b0[j[3]];state=2;marked_bit=i[1];return;}
    if( a0_s[i[4]]==1 and a0_s[i[3]]==1 and b0_s[j[0]]==1 and b0_s[j[2]]==1)
        {column_sum= a0[i[4]]^a0[i[3]]^b0[j[0]]^b0[j[2]];state=2;marked_bit=i[1];return;}
    if( a0_s[i[2]]==1 and a0_s[i[3]]==1 and b0_s[j[0]]==1 and b0_s[j[4]]==1)
        {column_sum= a0[i[2]]^a0[i[3]]^b0[j[0]]^b0[j[4]];state=2;marked_bit=i[1];return;}
    
    if( a0_s[i[4]]==1 and a0_s[i[1]]==1 and b0_s[j[2]]==1 and b0_s[j[3]]==1)
        {column_sum= a0[i[4]]^a0[i[1]]^b0[j[2]]^b0[j[3]];state=2;marked_bit=i[0];return;}
    if( a0_s[i[4]]==1 and a0_s[i[2]]==1 and b0_s[j[1]]==1 and b0_s[j[3]]==1)
        {column_sum= a0[i[4]]^a0[i[2]]^b0[j[1]]^b0[j[3]];state=2;marked_bit=i[0];return;}
    if( a0_s[i[4]]==1 and a0_s[i[3]]==1 and b0_s[j[1]]==1 and b0_s[j[2]]==1)
        {column_sum= a0[i[4]]^a0[i[3]]^b0[j[1]]^b0[j[2]];state=2;marked_bit=i[0];return;}
    if( a0_s[i[1]]==1 and a0_s[i[2]]==1 and b0_s[j[4]]==1 and b0_s[j[3]]==1)
        {column_sum= a0[i[1]]^a0[i[2]]^b0[j[4]]^b0[j[3]];state=2;marked_bit=i[0];return;}
    if( a0_s[i[1]]==1 and a0_s[i[3]]==1 and b0_s[j[4]]==1 and b0_s[j[2]]==1)
        {column_sum= a0[i[1]]^a0[i[3]]^b0[j[4]]^b0[j[2]];state=2;marked_bit=i[0];return;}
    if( a0_s[i[2]]==1 and a0_s[i[3]]==1 and b0_s[j[4]]==1 and b0_s[j[1]]==1)
        {column_sum= a0[i[2]]^a0[i[3]]^b0[j[4]]^b0[j[1]];state=2;marked_bit=i[0];return;}
    //state0
    column_sum=2;state=0;marked_bit=-1;
    return;
}


void COLUMNSUM(uint8_t* a0,uint8_t* b0,uint8_t* a0_s,uint8_t* b0_s,int column_index,uint8_t* Sigma,uint8_t* Sigma_s,int* MarkedBit)
{
    int x=column_index/64;
    int z=mod(column_index,64);
    int i[5]={0};
    i[0]=get_locate(x,0,z);i[1]=get_locate(x,1,z);i[2]=get_locate(x,2,z);i[3]=get_locate(x,3,z);i[4]=get_locate(x,4,z);
    uint8_t state_sum=0;
    int state_case=0;
    int marked_bit=0;
    judge_state(a0,b0,a0_s,b0_s,column_index,state_sum,state_case,marked_bit);
    if(state_case==1)
    {
        Sigma[column_index]=state_sum;
        Sigma_s[column_index]=1;
        return;
    }
    if(state_case==2)
    {
        Sigma[column_index]=state_sum;
        Sigma_s[column_index]=2;
        MarkedBit[column_index]=marked_bit;
        return;
    }
    if(state_case==0)
    {
        Sigma_s[column_index]=0;
        return;
    }
}

int LINEARTRANS(uint8_t* a0,uint8_t* b0,uint8_t* a0_s,uint8_t* b0_s,uint8_t* Sigma,uint8_t* Sigma_s,int* MarkedBit)
{
    int bit_cnt=0;
    for(int i0=0;i0<1600;i0++)
    {
        int i=sigma(i0);
        int i1=phi1(i0);
        int i2=phi2(i0);
        if(b0_s[i]==0 and a0_s[i0]==1 and Sigma_s[i1]==1 and Sigma_s[i2]==1)
        {
            b0[i]=a0[i0]^Sigma[i1]^Sigma[i2];b0_s[i]=1;bit_cnt+=1;
            bit_cnt+=CPKERNEL(a0,b0,a0_s,b0_s,phi0(i0));
        }
        if(b0_s[i]==1 and a0_s[i0]==0 and Sigma_s[i1]==1 and Sigma_s[i2]==1)
        {   
            a0[i0]=b0[i]^Sigma[i1]^Sigma[i2];a0_s[i0]=1;bit_cnt+=1;
            bit_cnt+=CPKERNEL(a0,b0,a0_s,b0_s,phi0(i0));
        }
        if(b0_s[i]==1 and a0_s[i0]==1 and Sigma_s[i1]==0 and Sigma_s[i2]==1)
        {    
            Sigma[i1]=b0[i]^a0[i0]^Sigma[i2];Sigma_s[i1]=1;
        }
        if(b0_s[i]==1 and a0_s[i0]==1 and Sigma_s[i1]==2 and Sigma_s[i2]==1)
        {
            a0[MarkedBit[i1]]=b0[i]^a0[i0]^Sigma[i1]^Sigma[i2];
            a0_s[MarkedBit[i1]]=1;bit_cnt+=1;
            Sigma[i1]=a0[MarkedBit[i1]]^Sigma[i1];Sigma_s[i1]=1;
            bit_cnt+=CPKERNEL(a0,b0,a0_s,b0_s,phi0(MarkedBit[i1]));
        }
        if(b0_s[i]==1 and a0_s[i0]==1 and Sigma_s[i1]==1 and Sigma_s[i2]==0)
        {
            Sigma[i2]=b0[i]^a0[i0]^Sigma[i1];Sigma_s[i2]=1;
        }
        if(b0_s[i]==1 and a0_s[i0]==1 and Sigma_s[i1]==1 and Sigma_s[i2]==2)
        {
            a0[MarkedBit[i2]]=b0[i]^a0[i0]^Sigma[i1]^Sigma[i2];
            a0_s[MarkedBit[i2]]=1;bit_cnt+=1;
            Sigma[i2]=a0[MarkedBit[i2]]^Sigma[i2];Sigma_s[i2]=1;
            bit_cnt+=CPKERNEL(a0,b0,a0_s,b0_s,phi0(MarkedBit[i2]));
        }
    }
    return bit_cnt;
}

int SIEVE(uint8_t* a0,uint8_t* b0,uint8_t* a1,uint8_t* a0_s,uint8_t* b0_s)
{
    int bit_cnt=0;
    for(int y=0;y<5;y++)
    {
        for(int z=0;z<64;z++)
        {
            int deta_out=0;
            for(int x_=0;x_<5;x_++)
            {
                deta_out+=((a1[get_locate(x_,y,z)])<<(4-x_));
            }
            int deta_in=0;
            int Tranc=0;
            for(int x_=0;x_<5;x_++)
            {
                deta_in+=((b0[get_locate(x_,y,z)])<<(4-x_));
                Tranc+=((b0_s[get_locate(x_,y,z)])<<(4-x_));
            }
            int T_deta_in=((deta_in<<5)+Tranc);
            uint8_t T_deta_in_vector[10]={0};
            bin_to_vector(T_deta_in,10,T_deta_in_vector);
            uint8_t* T=TDDT[T_deta_in][deta_out];
            uint8_t discard_vector[10]={2,2,2,2,2,2,2,2,2,2};
            //====================================================
            // if(y==0 and z==25)
            // {
            //     printf("%d,%d,[",y,z);
            //     for(int i=0;i<10;i++)
            //     {
            //         printf("%d,",T_deta_in_vector[i]);
            //     }
            //     printf("],%d\n",deta_out);
            // }
            //====================================================
            if(comp_vector(T,discard_vector,10)==1)
            {
                return -1;
            }
            if(comp_vector(T,T_deta_in_vector,10)==1)
            {
                continue;
            }
            else
            {
                int i[5]={0};
                i[0]=get_locate(0,y,z);i[1]=get_locate(1,y,z);i[2]=get_locate(2,y,z);i[3]=get_locate(3,y,z);i[4]=get_locate(4,y,z);
                for(int j=0;j<5;j++)
                {
                    if(T_deta_in_vector[j+5]==0 && T[j+5]==1)
                    {
                        b0[i[j]]=T[j];
                        b0_s[i[j]]=1;
                        bit_cnt++;
                        bit_cnt+=CPKERNEL(a0,b0,a0_s,b0_s,phi0(sigma_inverse(i[j])));
                    }
                }
            }
        }
    }
    return bit_cnt;
}

int SIEVE_back(uint8_t* a0,uint8_t* b0,uint8_t* a1,uint8_t* a0_s,uint8_t* b0_s)
{
    int bit_cnt=0;
    for(int y=0;y<5;y++)
    {
        for(int z=0;z<64;z++)
        {
            int deta_out=0;
            for(int x_=0;x_<5;x_++)
            {
                deta_out+=((a1[get_locate(x_,y,z)])<<(4-x_));
            }
            int deta_in=0;
            int Tranc=0;
            for(int x_=0;x_<5;x_++)
            {
                deta_in+=((b0[get_locate(x_,y,z)])<<(4-x_));
                Tranc+=((b0_s[get_locate(x_,y,z)])<<(4-x_));
            }
            int T_deta_in=((deta_in<<5)+Tranc);
            uint8_t T_deta_in_vector[10]={0};
            bin_to_vector(T_deta_in,10,T_deta_in_vector);
            uint8_t* T=TDDT[T_deta_in][deta_out];
            uint8_t discard_vector[10]={2,2,2,2,2,2,2,2,2,2};
            // if (y==0 and z==40)
            // {
            //     printf("%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n",y,z,T_deta_in_vector[0],
            //     T_deta_in_vector[1],T_deta_in_vector[2],T_deta_in_vector[3],T_deta_in_vector[4],
            //     T_deta_in_vector[5],T_deta_in_vector[6],T_deta_in_vector[7],T_deta_in_vector[8],
            //     T_deta_in_vector[9],deta_out);
            // }
            if(comp_vector(T,discard_vector,10)==1)
            {
                return -1;
            }
            if(comp_vector(T,T_deta_in_vector,10)==1)
            {
                continue;
            }
            else
            {
                int i[5]={0};
                i[0]=get_locate(0,y,z);i[1]=get_locate(1,y,z);i[2]=get_locate(2,y,z);i[3]=get_locate(3,y,z);i[4]=get_locate(4,y,z);
                //check back
                int update_bit[5]={-1,-1,-1,-1,-1};
                int update_bit_cnt=0;
                for(int j=0;j<5;j++)
                {
                    if(T_deta_in_vector[j+5]==0 && T[j+5]==1)
                    {
                        b0[i[j]]=T[j];
                        b0_s[i[j]]=1;
                        bit_cnt++;
                        update_bit[update_bit_cnt]=i[j];
                        update_bit_cnt++;
                        // bit_cnt+=CPKERNEL(a0,b0,a0_s,b0_s,phi0(sigma_inverse(i[j])));
                    }
                }
                if(update_bit_cnt>=2)
                {
                    int a0_bit[5]={-1,-1,-1,-1,-1};
                    int a0_bit_cnt=0;
                    for(int k=0;k<update_bit_cnt;k++)
                    {
                        if(a0_s[sigma_inverse(update_bit[k])]==1)
                        {
                            a0_bit[a0_bit_cnt]=sigma_inverse(update_bit[k]);
                            a0_bit_cnt++;
                            // printf("11111111\n");
                        }
                    }
                    if(a0_bit_cnt>=2)
                    {
                        if((a0[a0_bit[0]]^a0[a0_bit[1]])!=(b0[sigma(a0_bit[0])]^b0[sigma(a0_bit[1])]))
                        {
                            // for(int k=0;k<update_bit_cnt;k++)
                            // {
                            //     printf("%d,",update_bit[k]);
                            // }
                            // for(int k=0;k<a0_bit_cnt;k++)
                            // {
                            //     printf("%d,",a0_bit[k]);
                            // }
                            return -1;
                        }
                    }
                }
                for(int k=0;k<update_bit_cnt;k++)
                {
                    bit_cnt+=CPKERNEL(a0,b0,a0_s,b0_s,phi0(sigma_inverse(update_bit[k])));
                }
            }
        }
    }
    return bit_cnt;
}

int INITIALISEDP(uint8_t* a0,uint8_t* b0,uint8_t* a1,uint8_t* a0_s,uint8_t* b0_s)
{
    int bit_cnt=0;
    // printf("IN INITIUALDP,a0[829]=%d,a0_s[829]=%d\n",a0[829],a0_s[829]);
    for(int i=0;i<320;i++)
    {
        int x=i/64;
        int z=mod(i,64);
        bit_cnt+=CPKERNEL(a0,b0,a0_s,b0_s,i);
    }
    return bit_cnt;
}

int DP_back(uint8_t* a0,uint8_t* b0,uint8_t* a1,uint8_t* a0_s,uint8_t* b0_s,uint8_t* Sigma,uint8_t* Sigma_s,int* MarkedBit)
{
    // printf("INTO DP\n");
    int a=INITIALISEDP(a0,b0,a1,a0_s,b0_s);
    // printf("AFTER INITIAL,b0_s[168]=%d,b0[168]=%d\n",b0_s[168],b0[168]);
    while(a!=0)
    {
        for(int i=0;i<320;i++)
        {
            if(Sigma_s[i]==0)
            {
                COLUMNSUM(a0,b0,a0_s,b0_s,i,Sigma,Sigma_s,MarkedBit);
            }
            else if(Sigma_s[i]==2 && a0_s[MarkedBit[i]]==1)
            {
                Sigma_s[i]=1;
                Sigma[i]=Sigma[i]^a0[MarkedBit[i]];
            }
        }
        // printf("AFTER COLUMNSUM,b0_s[168]=%d,b0[168]=%d\n",b0_s[168],b0[168]);
        int linear_cnt=LINEARTRANS(a0,b0,a0_s,b0_s,Sigma,Sigma_s,MarkedBit);
        // printf("AFTER LINEAR,b0_s[168]=%d,b0[168]=%d\n",b0_s[168],b0[168]);
        int flag=SIEVE_back(a0,b0,a1,a0_s,b0_s);
        if(flag==-1)
        {
            return -1;
        }
        else
        {
            a=flag;
        }
    }
    return 1;
}

int DP(uint8_t* a0,uint8_t* b0,uint8_t* a1,uint8_t* a0_s,uint8_t* b0_s,uint8_t* Sigma,uint8_t* Sigma_s,int* MarkedBit)
{
    // printf("INTO DP\n");
    int a=INITIALISEDP(a0,b0,a1,a0_s,b0_s);
    // printf("AFTER INITIAL,b0_s[168]=%d,b0[168]=%d\n",b0_s[168],b0[168]);
    while(a!=0)
    {
        for(int i=0;i<320;i++)
        {
            if(Sigma_s[i]==0)
            {
                COLUMNSUM(a0,b0,a0_s,b0_s,i,Sigma,Sigma_s,MarkedBit);
            }
            else if(Sigma_s[i]==2 && a0_s[MarkedBit[i]]==1)
            {
                Sigma_s[i]=1;
                Sigma[i]=Sigma[i]^a0[MarkedBit[i]];
            }
        }
        // printf("AFTER COLUMNSUM,b0_s[345]=%d,b0[345]=%d\n",b0_s[345],b0[345]);
        int linear_cnt=LINEARTRANS(a0,b0,a0_s,b0_s,Sigma,Sigma_s,MarkedBit);
        // printf("AFTER LINEAR,b0_s[345]=%d,b0[345]=%d\n",b0_s[345],b0[345]);
        int flag=SIEVE(a0,b0,a1,a0_s,b0_s);
        if(flag==-1)
        {
            return -1;
        }
        else
        {
            a=flag;
        }
    }
    return 1;
}

void INITIALVP(uint8_t* a0,uint8_t* b0,uint8_t* a1,uint8_t* a0_s,uint8_t* b0_s,uint8_t* B,uint8_t* B_p,uint8_t* Bs,uint8_t* Bs_p)
{
    for(int y=0;y<5;y++)
    {
        for(int z=0;z<64;z++)
        {
            int deta_out=0;
            for(int x_=0;x_<5;x_++)
            {
                deta_out+=((a1[get_locate(x_,y,z)])<<(4-x_));
            }
            int deta_in=0;
            int Tranc=0;
            for(int x_=0;x_<5;x_++)
            {
                deta_in+=((b0[get_locate(x_,y,z)])<<(4-x_));
                Tranc+=((b0_s[get_locate(x_,y,z)])<<(4-x_));
            }
            int T_deta_in=((deta_in<<5)+Tranc);
            uint8_t T_deta_in_vector[10]={0};
            bin_to_vector(T_deta_in,10,T_deta_in_vector);
            uint8_t* v=FVDT[T_deta_in][deta_out];
            uint8_t discard_vector[10]={2,2,2,2,2,2,2,2,2,2};
            if(comp_vector(v,discard_vector,10)==1)
            {
                continue;
            }
            else
            {
                int i[5]={0};
                i[0]=get_locate(0,y,z);i[1]=get_locate(1,y,z);i[2]=get_locate(2,y,z);i[3]=get_locate(3,y,z);i[4]=get_locate(4,y,z);
                for(int j=0;j<5;j++)
                {
                    if(v[j+5]==1)
                    {
                        B[i[j]]=v[j];B_p[i[j]]=v[j];Bs[i[j]]=1;Bs_p[i[j]]=1;
                    }
                }
            }
        }
    }
}

int UPDATE(uint8_t* a0,uint8_t* b0,uint8_t* a1,uint8_t* a0_s,uint8_t* b0_s,uint8_t* B,uint8_t* B_p,uint8_t* Bs,uint8_t* Bs_p)
{
    int flag=0;
    for(int i=0;i<1600;i++)
    {
        if(b0_s[i]==0 and Bs[i]==1 and Bs_p[i]==1)
        {
            b0_s[i]=1;b0[i]=B[i]^B_p[i];flag=1;
        }
    }
    // printf("After 1600,b0_s=%d\n",b0_s[1305]);
    for(int x=0;x<5;x++)
    {
        for(int y=0;y<5;y++)
        {
            for(int z=0;z<64;z++)
            {
                int i0=get_locate(x,y,z);int i1=get_locate(mod(x+1,5),y,z);int i3=get_locate(mod(x+3,5),y,z);int i4=get_locate(mod(x+4,5),y,z);
                // if(x==2 and y==4 and z==25)
                // {
                //     printf("Before floop%d,%d,%d\n",i1,i3,i4);
                //     printf("%d,%d,%d\n",x,y,z);
                //     printf("%d\n---1112222211----------------\n",b0_s[1305]);
                // }
                if(Bs[i0]==1 and Bs_p[i0]==1)
                {
                    if(B[i0]==0 and B_p[i0]==0)
                    {
                        if(b0_s[i3]==0)
                        {
                            b0_s[i3]=1;b0[i3]=a1[i3];flag=1;
                        }
                        if(b0_s[i1]==0 and b0_s[i4]==1)
                        {
                            b0_s[i1]=1;b0[i1]=a1[i4]^b0[i4];flag=1;
                        }
                        else if(b0_s[i1]==1 and b0_s[i4]==0)
                        {
                            b0_s[i4]=1;b0[i4]=a1[i4]^b0[i1];flag=1;
                        }
                        // if(x==2 and y==4 and z==25)
                        // {
                        //     printf("In floop%d,%d,%d\n",i1,i3,i4);
                        //     printf("%d,%d,%d\n",x,y,z);
                        //     printf("%d\n---11111111111----------------\n",b0_s[1305]);
                        // }
                        continue;
                    }
                    if(B[i0]==1 and B_p[i0]==1)
                    {
                        if(b0_s[i4]==0)
                        {
                            b0_s[i4]=1;b0[i4]=a1[i4];flag=1;
                        }
                        if(b0_s[i3]==0 and b0_s[i4]==1)
                        {
                            b0_s[i3]=1;b0[i3]=a1[i3]^b0[i4];flag=1;
                        }
                        else if(b0_s[i3]==1 and b0_s[i4]==0)
                        {
                            b0_s[i4]=1;b0[i4]=a1[i3]^b0[i3];flag=1;
                        }
                        // if(x==2 and y==4 and z==25)
                        // {
                        //     printf("In floop%d,%d,%d\n",i1,i3,i4);
                        //     printf("%d,%d,%d\n",x,y,z);
                        //     printf("%d\n---222222222----------------\n",b0_s[1305]);
                        // }
                        continue;
                    }
                    if(B[i0]==1 and B_p[i0]==0)
                    {
                        if(b0_s[i4]==0 and Bs_p[i1]==1)
                        {
                            b0_s[i4]=1;b0[i4]=a1[i4]^B_p[i1];flag=1;
                        }
                        if(b0_s[i3]==0 and Bs[i4]==1)
                        {
                            b0_s[i3]=1;b0[i3]=a1[i3]^B[i4]^1;flag=1;
                        }
                        // if(x==2 and y==4 and z==25)
                        // {
                        //     printf("In floop,%d,%d,%d,%d\n",i0,i1,i3,i4);
                        //     printf("%d,%d,%d\n",x,y,z);
                        //     printf("%d\n---333333333----------------\n",b0_s[1305]);
                        // }
                        continue;
                    }
                    if(B[i0]==0 and B_p[i0]==1)
                    {
                        if(b0_s[i4]==0 and Bs[i1]==1)
                        {
                            b0_s[i4]=1;b0[i4]=a1[i4]^B[i1];flag=1;
                        }
                        if(b0_s[i3]==0 and Bs_p[i4]==1)
                        {
                            b0_s[i3]=1;b0[i3]=a1[i3]^B_p[i4]^1;flag=1;
                        }
                        // if(x==2 and y==4 and z==25)
                        // {
                        //     printf("In floop,%d,%d,%d,%d\n",i0,i1,i3,i4);
                        //     printf("%d,%d,%d\n",x,y,z);
                        //     printf("%d\n---444444444----------------\n",b0_s[1305]);
                        //     printf("%d\n",Bs_p[1369]);
                        // }
                        continue;
                    }
                }
            }
        }
    }
    return flag;
}

int VP(uint8_t* a0,uint8_t* b0,uint8_t* a1,uint8_t* a0_s,uint8_t* b0_s,uint8_t* A,uint8_t* A_p,uint8_t* As,uint8_t* As_p,uint8_t* B,uint8_t* B_p,uint8_t* Bs,uint8_t* Bs_p)
{
    INITIALVP(a0,b0,a1,a0_s,b0_s,B,B_p,Bs,Bs_p);
    // printf("After Initial,b0_s[1305]=%d\n",b0_s[1305]);
    // printf("After Initial,Bs_p[1369]=%d\n",Bs_p[1369]);
    // printf("This is A_p\n");
    // for(int i=0;i<1600;i++)
    // {
    //     printf("%d,",A_p[i]);
    // }
    // printf("\nThis is As_p\n");
    // for(int i=0;i<1600;i++)
    // {
    //     printf("%d,",As_p[i]);
    // }
    // printf("\nThis is A\n");
    // for(int i=0;i<1600;i++)
    // {
    //     printf("%d,",A[i]);
    // }
    // printf("\nThis is As\n");
    // for(int i=0;i<1600;i++)
    // {
    //     printf("%d,",As[i]);
    // }
    // printf("\n");
    for(int i=0;i<320;i++)
    {
        CPKERNEL(A,B,As,Bs,i);
        CPKERNEL(A_p,B_p,As_p,Bs_p,i);
    }
    // printf("After Kernerl,b0_s[1305]=%d\n",b0_s[1305]);
    // printf("Before update,Bs_p[1369]=%d\n",Bs_p[1369]);
    int a=UPDATE(a0,b0,a1,a0_s,b0_s,B,B_p,Bs,Bs_p);
    // printf("After update,b0_s[1305]=%d\n",b0_s[1305]);
    if(a==0)
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

int DS_41bit_back(uint8_t* M1,uint8_t* M1_p,uint8_t* _A,uint8_t* _A_p)
{
    uint8_t A[1600]={0};
    uint8_t A_p[1600]={0};
    uint8_t As[1600]={0};
    uint8_t As_p[1600]={0};
    uint8_t B[1600]={0};
    uint8_t B_p[1600]={0};
    uint8_t Bs[1600]={0};
    uint8_t Bs_p[1600]={0};
    uint8_t a0[1600]={0};
    uint8_t a0_s[1600]={0};
    uint8_t b0[1600]={0};
    uint8_t b0_s[1600]={0};
    uint8_t Sigma[320]={0};
    uint8_t Sigma_s[320]={0};
    int MarkedBit[320]={0};
    set<int> S0=get_S0();
    set<int> S1=get_S1();
    INITIAL_41bit(M1,M1_p,S0,S1,A,A_p,As,As_p,B,B_p,Bs,Bs_p,a0,b0,a0_s,b0_s,Sigma,Sigma_s,MarkedBit);
    int flag=1;
    while(flag>0)
    {
        flag=DP_back(a0,b0,a1,a0_s,b0_s,Sigma,Sigma_s,MarkedBit);
        if(flag>-1)
        {
            flag=VP(a0,b0,a1,a0_s,b0_s,A,A_p,As,As_p,B,B_p,Bs,Bs_p);
            if(flag==0)
            {
                memcpy(_A,A,1600);
                memcpy(_A_p,A_p,1600);
                return 1;
            }
        }
        else
        {
            return 0;
        }
    }
    return flag;
}

int DS_39bit(uint8_t* M1,uint8_t* M1_p,uint8_t* _A,uint8_t* _A_p)
{
    uint8_t A[1600]={0};
    uint8_t A_p[1600]={0};
    uint8_t As[1600]={0};
    uint8_t As_p[1600]={0};
    uint8_t B[1600]={0};
    uint8_t B_p[1600]={0};
    uint8_t Bs[1600]={0};
    uint8_t Bs_p[1600]={0};
    uint8_t a0[1600]={0};
    uint8_t a0_s[1600]={0};
    uint8_t b0[1600]={0};
    uint8_t b0_s[1600]={0};
    uint8_t Sigma[320]={0};
    uint8_t Sigma_s[320]={0};
    int MarkedBit[320]={0};
    set<int> S0=get_S0();
    set<int> S1=get_S1();
    uint8_t initial_flag=INITIAL_discard(M1,M1_p,S0,S1,A,A_p,As,As_p,B,B_p,Bs,Bs_p,a0,b0,a0_s,b0_s,Sigma,Sigma_s,MarkedBit);
    if(initial_flag==1)//discard
    {
        return 0;
    }
    int flag=1;
    while(flag>0)
    {
        // printf("1111111111111111111111\n");
        flag=DP(a0,b0,a1,a0_s,b0_s,Sigma,Sigma_s,MarkedBit);
        if(flag>-1)
        {
            // printf("update_bit:%d\n",flag);
            // printf("%d,%d\n",b0_s[1305],b0[1305]);
            flag=VP(a0,b0,a1,a0_s,b0_s,A,A_p,As,As_p,B,B_p,Bs,Bs_p);
            // printf("update_bit:%d\n",flag);
            // printf("%d,%d\n",b0_s[1305],b0[1305]);
            // for(int i=0;i<1600;i++)
            // {
            //     printf("%d,",b0[i]);
            // }
            if(flag==0)
            {
                memcpy(_A,A,1600);
                memcpy(_A_p,A_p,1600);
                return 1;
            }
        }
        else
        {
            return 0;
        }
    }
    return flag;
}

int DS_inversestate(uint8_t* M1,uint8_t* M1_p,uint8_t* _A,uint8_t* _A_p)
{
    uint8_t A[1600]={0};
    uint8_t A_p[1600]={0};
    uint8_t As[1600]={0};
    uint8_t As_p[1600]={0};
    uint8_t B[1600]={0};
    uint8_t B_p[1600]={0};
    uint8_t Bs[1600]={0};
    uint8_t Bs_p[1600]={0};
    uint8_t a0[1600]={0};
    uint8_t a0_s[1600]={0};
    uint8_t b0[1600]={0};
    uint8_t b0_s[1600]={0};
    uint8_t Sigma[320]={0};
    uint8_t Sigma_s[320]={0};
    int MarkedBit[320]={0};
    set<int> S0=get_S0();
    set<int> S1=get_S1();
    uint8_t discard_flag=INITIAL_inversestate(M1,M1_p,S0,S1,A,A_p,As,As_p,B,B_p,Bs,Bs_p,a0,b0,a0_s,b0_s,Sigma,Sigma_s,MarkedBit);
    printf("Discard_condition:%x\n",discard_flag);
    if(discard_flag==1)
    {
        return 0;
    }
    int flag=1;
    while(flag>0)
    {
        flag=DP(a0,b0,a1,a0_s,b0_s,Sigma,Sigma_s,MarkedBit);
        if(flag>-1)
        {
            flag=VP(a0,b0,a1,a0_s,b0_s,A,A_p,As,As_p,B,B_p,Bs,Bs_p);
            if(flag==0)
            {
                memcpy(_A,A,1600);
                memcpy(_A_p,A_p,1600);
                return 1;
            }
        }
        else
        {
            return 0;
        }
    }
    return flag;
}

int DS_only_DP(uint8_t* M1,uint8_t* M1_p)
{
    uint8_t A[1600]={0};
    uint8_t A_p[1600]={0};
    uint8_t As[1600]={0};
    uint8_t As_p[1600]={0};
    uint8_t B[1600]={0};
    uint8_t B_p[1600]={0};
    uint8_t Bs[1600]={0};
    uint8_t Bs_p[1600]={0};
    uint8_t a0[1600]={0};
    uint8_t a0_s[1600]={0};
    uint8_t b0[1600]={0};
    uint8_t b0_s[1600]={0};
    uint8_t Sigma[320]={0};
    uint8_t Sigma_s[320]={0};
    int MarkedBit[320]={0};
    set<int> S0=get_S0();
    set<int> S1=get_S1();
    INITIAL_39bit(M1,M1_p,S0,S1,A,A_p,As,As_p,B,B_p,Bs,Bs_p,a0,b0,a0_s,b0_s,Sigma,Sigma_s,MarkedBit);
    int flag=1;
    return DP(a0,b0,a1,a0_s,b0_s,Sigma,Sigma_s,MarkedBit);
}