#include "keccak_permutation.h"
#include <stdio.h>
#include <unistd.h>
uint64_t hash_table[16777216][13]={0};
uint8_t hash_flag[16777216]={0};
uint64_t get_index(uint8_t* M)
{
    uint8_t condition[24]={0};
    condition[0]=M[863]^M[1183];
    condition[1]=M[864]^M[1184];
    condition[2]=M[877]^M[1197];
    condition[3]=M[878]^M[1198];
    condition[4]=M[879]^M[1199];
    condition[5]=M[880]^M[1200];
    condition[6]=M[881]^M[1201];
    condition[7]=M[882]^M[1202];
    condition[8]=M[883]^M[1203];
    condition[9]=M[884]^M[1204];
    condition[10]=M[885]^M[1205];
    condition[11]=M[934]^M[1574];
    condition[12]=M[935]^M[1575];
    condition[13]=M[995]^M[1315];
    condition[14]=M[996]^M[1316];
    condition[15]=M[997]^M[1317];
    condition[16]=M[1032]^M[1352];
    condition[17]=M[1033]^M[1353];
    condition[18]=M[1034]^M[1354];
    condition[19]=M[1035]^M[1355];
    condition[20]=M[1036]^M[1356];
    condition[21]=M[1144]^M[1464];
    condition[22]=M[1145]^M[1465];
    condition[23]=M[1146]^M[1466];
    uint64_t index=0;
    for(int i=0;i<24;i++)
    {
        index+=(uint64_t(condition[i]))*(uint64_t(pow(2,i)));
    }
    return index;
}

int main(int argc,char* argv[])
{
    // int var_list[3]={1,2,3};
    // int point[3]={1,-1,0};
    // string ip=add_impossible_point(var_list,point,3);
    // printf("%s\n",ip.data());

    // int var_list[3]={1,0,3};
    // int result_bit=1;
    // string Clause[2];
    // Xor_sum(var_list,result_bit,3,Clause);
    // for(int i=0;i<2;i++)
    // {
    //     printf("%s\n",Clause[i].data());
    // }
    

    srand(time(NULL));
    for(uint64_t i=0;i<pow(2,1);i++)
    {
        FILE *f=fopen(("message"+to_string(i)).data(),"w");
        for(uint64_t j=0;j<pow(2,25);j++)
        {
            if(j%(uint64_t)(pow(2,15))==0)
            {
                printf("%ld\n",j);
            }
            uint64_t rand_num_list[13]={0};
            uint8_t M1[1600]={0};
            uint8_t hash_value[1600]={0};
            for(uint64_t k=0;k<13;k++)
            {
                uint64_t rand_num0=(uint64_t)(rand());
                uint64_t rand_num1=(uint64_t)(rand());
                uint64_t rand_num2=(uint64_t)(rand());
                uint64_t rand_num=(rand_num0<<33)+(rand_num1<<2)+(rand_num2%4);
                for(int l=0;l<64;l++)
                {
                    bin64_to_vector(rand_num,M1+k*64);
                }
                rand_num_list[k]=rand_num;
            }
            // printf("M1=[");
            // for(int l=0;l<1600;l++)
            // {
            //     printf("%d,",M1[l]);
            // }
            // printf("]\n");
            f_4(M1,hash_value);
            // printf("hash_value=[");
            // for(int l=0;l<1600;l++)
            // {
            //     printf("%d,",hash_value[l]);
            // }
            // printf("]\n");
            uint64_t index=get_index(hash_value);
            if(hash_flag[index]==1)
            {
                for(int l=0;l<13;l++)
                {
                    fprintf(f,"%lu,",hash_table[index][l]);
                }
                fprintf(f,"\n");
                for(int l=0;l<13;l++)
                {
                    fprintf(f,"%lu,",rand_num_list[l]);
                }
                fprintf(f,"\n");
            }
            memcpy(hash_table[index],rand_num_list,13*8);
            hash_flag[index]=1;
            // printf("index=%ld\n",index);
        }
        fclose(f);
    }

    // string a="";
    // string b[3];
    // b[0]=a+"123";
    // b[1]=a+"456";
    // b[2]=a+"789";
    // for(int i=0;i<3;i++)
    // {
    //     printf("%s\n",b[i].data());
    // }

        // uint64_t a[13]={1,2,3,4,5,6,7,8,9,1,2,3,4};
        // memcpy(hash_table[3],a,13*8);
        // for(int i=0;i<13;i++)
        // {
        //     printf("%ld\n",hash_table[3][i]);
        // }


    // uint64_t a=pow(2,63)+pow(2,62);
    // uint8_t result[64]={0};
    // bin64_to_vector(a,result);
    // for(int i=0;i<64;i++)
    // {
    //     printf("%d,",result[i]);
    // }
    // printf("\n");

    // uint64_t x=10000000000000;
    // cout<<(uint64_t(pow(2,63))<<1)<<endl;
    
    // cout<<time(NULL)<<endl;
    // srand(time(NULL));
    // for(int i=0;i<10;i++)
    // {
    //     cout<<rand()<<endl;
    // }
    
    // string file="message";
    // for (int i=0;i<3;i++)
    // {
    //     FILE *f=fopen((file+to_string(i)).data(),"w");
    //     fprintf(f,"123");
    //     fclose(f);
    //     sleep(10);
    // }
}