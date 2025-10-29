import numpy as np
l=64
p=2
n=5
Field=GF(p)
def get_locate(x,y,z):
    return 64*(5*y+x)+z

def theta(a):
    result=[0]*25*l
    for x in range(5):
        for y in range(5):
            for z in range(l):
                _sum=0
                for y_pie in range(5):
                    _sum=_sum+a[get_locate((x-1)%5,y_pie,z)]+a[get_locate((x+1)%5,y_pie,(z-1)%l)]
                result[get_locate(x,y,z)]=a[get_locate(x,y,z)]+_sum
    return result

t_step=[[0%l,36%l,3%l,105%l,210%l],
        [1%l,300%l,10%l,45%l,66%l],
        [190%l,6%l,171%l,15%l,253%l],
        [28%l,55%l,153%l,21%l,120%l],
        [91%l,276%l,231%l,136%l,78%l]]

def rho(a):
    result=[0]*25*l
    for x in range(5):
        for y in range(5):
            for z in range(l):
                result[get_locate(x,y,z)]=a[get_locate(x,y,(z-t_step[x][y])%l)]
    return result

pi_step=[[(0,0),(3,0),(1,0),(4,0),(2,0)],
         [(1,1),(4,1),(2,1),(0,1),(3,1)],
         [(2,2),(0,2),(3,2),(1,2),(4,2)],
         [(3,3),(1,3),(4,3),(2,3),(0,3)],
         [(4,4),(2,4),(0,4),(3,4),(1,4)]]
pi_inverse=[[(0,0) for _ in range(5)] for _ in range(5)]
for i in range(5):
    for j in range(5):
        new_i=pi_step[i][j][0]
        new_j=pi_step[i][j][1]
        pi_inverse[new_i][new_j]=(i,j)

def pi(a):
    result=[0]*25*l
    for x in range(5):
        for y in range(5):
            for z in range(l):
                result[get_locate(x,y,z)]=a[get_locate(pi_step[x][y][0],pi_step[x][y][1],z)]
    return result
    

def chi(a):
    result=[0]*25*l
    for x in range(5):
        for y in range(5):
            for z in range(l):
                result[get_locate(x,y,z)]=a[get_locate(x,y,z)]+a[get_locate((x+1)%5,y,z)]*a[get_locate((x+2)%5,y,z)]+a[get_locate((x+2)%5,y,z)]
    return result

RC=[0x0000000000000001,0x0000000000008082,0x800000000000808A,
    0x8000000080008000,0x000000000000808B,0x0000000080000001,
    0x8000000080008081,0x8000000000008009,0x000000000000008A,
    0x0000000000000088,0x0000000080008009,0x000000008000000A,
    0x000000008000808B,0x800000000000008B,0x8000000000008089,
    0x8000000000008003,0x8000000000008002,0x8000000000000080,
    0x000000000000800A,0x800000008000000A,0x8000000080008081,
    0x8000000000008080,0x0000000080000001,0x8000000080008008]

def bin_to_vector(x,length):
    result=[]
    for i in range(length):
        b=x//2**(length-1)
        result.append(Field(b))
        x=(x<<1)%(2**length)
    return result

def vector_to_bin(x,start,end):
    result=0
    length=end-start+1
    for i in range(length):
        result=result+Integer(x[start+i])*2**(length-1-i)
    return result

def iota(a,ir):
    round_key=bin_to_vector(RC[ir],l)
    for z in range(l):
        a[get_locate(0,0,z)]=a[get_locate(0,0,z)]+round_key[z]
    return a

def keccak_f(s,ir):
    return iota(chi(pi(rho(theta(s)))),ir)

def chi_5(x):
    result=[Field(0) for _ in range(5)]
    for i in range(5):
        result[i]=x[i]+x[(i+1)%5]*x[(i+2)%5]+x[(i+2)%5]
    return vector(result)
DDT=[[0 for _ in range(32)] for _ in range(32)]
DDT_value=[[[] for _ in range(32)] for _ in range(32)]
DDT_out=[[] for _ in range(32)]
DDT_in=[[] for _ in range(32)]
p=2
n=5
r=4
Field=GF(p)
def p_decimal_to_Decimal(x):
    SUM=0
    for i in range(n):
        SUM+=Integer(x[i])*(p)**(n-i-1)
    return SUM

def Decimal_to_p_decimal(x):
    result=[Field(0) for _ in range(n)]
    for i in range(n):
        remainder=x%p
        result[n-i-1]=Field(remainder)
        x=x//p
    return result

for i in range(32):
    vector_i=vector(Decimal_to_p_decimal(i))
    for deta_in in range(32):
        vector_deta_in=vector(Decimal_to_p_decimal(deta_in))
        vector_j=vector_i+vector_deta_in
        vector_i_out=chi_5(vector_i)
        vector_j_out=chi_5(vector_j)
        vector_deta_out=vector_i_out+vector_j_out
        deta_out=p_decimal_to_Decimal(vector_deta_out)
        DDT[deta_in][deta_out]+=1
        DDT_value[deta_in][deta_out].append(list(vector_i))
        DDT_out[deta_out].append(list(vector_deta_in))
        DDT_in[deta_in].append(list(vector_deta_out))

def print_DDT(DDT):
    for i in range(32):
        for j in range(32):
            print(DDT[i][j],end=",")
        print()
        
def pusai0(i):
    return (i//64)%5

def pusai1(i):
    return i//320

def pusai2(i):
    return i%64

def phi0(i):
    return 64*pusai0(i)+pusai2(i)

def phi1(i):
    return 64*((pusai0(i)-1)%5)+pusai2(i)

def phi2(i):
    return 64*((pusai0(i)+1)%5)+(pusai2(i)-1)%l

def sigma(i):
    x=pusai0(i)
    y=pusai1(i)
    z=pusai2(i)
    z_1=(z+t_step[x][y])%l
    x_1=pi_inverse[x][y][0]
    y_1=pi_inverse[x][y][1]
    return get_locate(x_1,y_1,z_1)

def sigma_inverse(i):
    x=pusai0(i)
    y=pusai1(i)
    z=pusai2(i)
    x_1=pi_step[x][y][0]
    y_1=pi_step[x][y][1]
    z_1=(z-t_step[x_1][y_1])%l
    return get_locate(x_1,y_1,z_1)

alpha1=[0x7c0bc4f5b4398002,0x2407de4bc9668001,0xac02095d32eb8000,0xd402e98975068000,0x3c05706a07f58000
,0x7c0bccf5b4398002,0x240fde4bc9e68001,0xac02095d32ef8000,0xc40ae98975068000,0x3414706a05f58000
,0x7c0bc4f5b4398000,0x240fda4bc9e68001,0xac02095d32eb8000,0xc40ae9897d068000,0x3c15706a25f58000
,0x7c0bc4f5bc398002,0x240fde4fc9668001,0xac02095d32eb8000,0xc40ae98975068000,0x3c15706a05f48000
,0x7c0bc4f1b4398002,0x240fde4bc9e68001,0xac02095d3aeb8000,0xd40ae98975868000,0x3c15706a05f58000]

a1=[0 for _ in range(1600)]
for i in range(1600):
    a1[i]=Field(alpha_i(i))

def alpha_i(i):
    x=pusai0(i)
    y=pusai1(i)
    z=pusai2(i)
    index=y*5+x
    return (alpha1[index]>>z)&0x0000000000000001

def regular(a_vector):
    for i in range(n):
        if a_vector[i]==1 and a_vector[n+i]==0:
            return 0
    return 1

TDTT=[[-1 for _ in range(2**n)] for _ in range(2**(2*n))]
def construct_TDTT(DDT):
    for a in range(2**n):
        a_vector=bin_to_vector(a,n)
        for tranc in range(2**n):
            tranc_vector=bin_to_vector(tranc,n)
            if regular(a_vector+tranc_vector)==0:
                continue
            for b in range(2**n):
                b_vector=bin_to_vector(b,n)
                solve_value=[0 for _ in range(n)]
                first_time=1
                have_solve=0
                solve_tranc=[1 for _ in range(n)]
                for j in range(len(DDT_out[b])):
                    flag=1
                    for k in range(5):
                        if tranc_vector[k]==1 and DDT_out[b][j][k]!=a_vector[k]:
                            flag=0;break
                    if flag==0:
                        continue
                    have_solve=1
                    if first_time==1:
                        solve_value=DDT_out[b][j].copy()
                        first_time=0
                    else:
                        for k in range(5):
                            if tranc_vector[k]==0 and solve_tranc[k]==1 and DDT_out[b][j][k]!=solve_value[k]:
                                solve_tranc[k]=0;solve_value[k]=0
                if have_solve==1:
                    TDTT[(a<<5)+tranc][b]=solve_value+solve_tranc
                    
FVDT=[[-1 for _ in range(2**n)] for _ in range(2**(2*n))]
def construct_FVDT():
    for a in range(2**n):
        a_vector=bin_to_vector(a,n)
        for tranc in range(2**n):
            tranc_vector=bin_to_vector(tranc,n)
            if regular(a_vector+tranc_vector)==0:
                continue
            for b in range(2**n):
                b_vector=bin_to_vector(b,n)
                solve_value=[0 for _ in range(n)]
                first_time=1
                have_solve=0
                solve_tranc=[1 for _ in range(n)]
                for j in range(len(DDT_out[b])):
                    flag=1
                    for k in range(5):
                        if tranc_vector[k]==1 and DDT_out[b][j][k]!=a_vector[k]:
                            flag=0;break
                    if flag==0:
                        continue
                    have_solve=1
                    for ele in DDT_value[vector_to_bin(DDT_out[b][j],0,4)][b]:
                        if first_time==1:
                            solve_value=ele.copy()
                            first_time=0
                        else:
                            for k in range(5):
                                if solve_tranc[k]==1 and ele[k]!=solve_value[k]:
                                    solve_tranc[k]=0;solve_value[k]=0
                if have_solve==1:
                    FVDT[(a<<5)+tranc][b]=solve_value+solve_tranc 
construct_TDTT(DDT)
construct_FVDT()
def Type1_condition():
    S0=[]
    S1=[]
    SA=[]
    matrix0=[]
    matrix1=[]
    for i in range(1600):
        x=pusai0(i)
        y=pusai1(i)
        z=pusai2(i)
        deta_out=0
        for x_ in range(5):
            deta_out+=(alpha_i(get_locate(x_,y,z))<<(4-x_))
        if deta_out==0:
            S0.append(i)
        else:
            if deta_out in [0x1,0x2,0x4,0x8,0x10] and alpha_i(i)==1:
                S1.append(i)
    for i in range(828,1600):
        for j in range(i+1,1600):
            if ((sigma(i) in S0 and sigma(j) in S0) or (sigma(i) in S1 and sigma(j) in S1)) and phi0(i)==phi0(j):
                vec=[0 for _ in range(1600-828)]
                vec[i-828]=1;vec[j-828]=1
                matrix0.append(vec+[0])
            else:
                if ((sigma(i) in S0 and sigma(j) in S1) or (sigma(i) in S1 and sigma(j) in S0)) and phi0(i)==phi0(j):
                    vec=[0 for _ in range(1600-828)]
                    vec[i-828]=1;vec[j-828]=1
                    matrix0.append(vec+[1])
    for i in range(828):
        for j in range(828,1600):
            if ((sigma(i) in S0 and sigma(j) in S0) or (sigma(i) in S1 and sigma(j) in S1)) and phi0(i)==phi0(j):
                vec=[0 for _ in range(1600)]
                vec[i]=1;vec[j]=1
                matrix1.append(vec+[0])
            else:
                if ((sigma(i) in S0 and sigma(j) in S1) or (sigma(i) in S1 and sigma(j) in S0)) and phi0(i)==phi0(j):
                    vec=[0 for _ in range(1600)]
                    vec[i]=1;vec[j]=1
                    matrix1.append(vec+[1])
    return S0,S1,matrix0,matrix1


def print_state_with_sbox(state):
    for y in range(5):
        for z in range(64):
            print(y,z,":",state[get_locate(0,y,z)],state[get_locate(1,y,z)],state[get_locate(2,y,z)],state[get_locate(3,y,z)],state[get_locate(4,y,z)])

def f_4(state):
    for i in range(4):
        state=theta(state)
        state=rho(state)
        state=pi(state)
        state=chi(state)
        state=iota(state,i)
    return state

def Update_Column(a0_s,b0_s,column_index):
    i=[0 for _ in range(5)]
    x=column_index//64
    z=column_index%64
    i[0]=get_locate(x,0,z);i[1]=get_locate(x,1,z);i[2]=get_locate(x,2,z);i[3]=get_locate(x,3,z);i[4]=get_locate(x,4,z)
    flag=0
    for j in range(5):
        if a0_s[i[j]]==1 and b0_s[sigma(i[j])]==1:
            flag=1
            break
    if flag==1:
        for j in range(5):
            if a0_s[i[j]]==1 and b0_s[sigma(i[j])]==0:
                b0_s[sigma(i[j])]=1
            else:
                if a0_s[i[j]]==0 and b0_s[sigma(i[j])]==1:
                    a0_s[i[j]]=1

def Deduce_known_bit_position():
    a0_s=[0 for _ in range(1600)];b0_s=[0 for _ in range(1600)]
    for i in range(1600):
        if i>=828:
            a0_s[i]=1
        else:
            a0_s[i]=0
        if i in S0:
            b0_s[i]=1
        else:
            if i in S1:
                b0_s[i]=1
            else:
                b0_s[i]=0
    for i in range(320):
        Update_Column(a0_s,b0_s,i)
    return a0_s,b0_s

def Vector_equal(v1,v2,T,i):
    for j in range(len(T)):
        if T[j]==1 and j!=i:
            if v1[j]!=v2[j]:
                return 0
    return 1

def print_condition(discard_condition):
    for i in range(len(discard_condition)-1):
        print(discard_condition[i],end=' or ')
    print(discard_condition[len(discard_condition)-1])

def Derive_new_conditions():
    a0_s,b0_s=Deduce_known_bit_position()
    discard_position=[]
    for y in range(4):
        for z in range(64):
            discard_vector=set()
            #Get S_discard
            bit0=get_locate(0,y,z);bit1=get_locate(1,y,z);bit2=get_locate(2,y,z);bit3=get_locate(3,y,z);bit4=get_locate(4,y,z)
            output_diff=vector_to_bin([a1[bit0],a1[bit1],a1[bit2],a1[bit3],a1[bit4]],0,4)
            S_discard=[]
            flag=-1
            flag_value=-1
            T=[b0_s[bit0],b0_s[bit1],b0_s[bit2],b0_s[bit3],b0_s[bit4]]
            T_num=vector_to_bin(T,0,4)
            for input_diff_num in range(32):
                if regular(bin_to_vector(input_diff_num,5)+T)==0:
                    continue
                input_diff_T=input_diff_num*32+T_num
                input_diff_vec=bin_to_vector(input_diff_num,5)
                dis_flag=1
                for x in range(5):
                    if get_locate(x,y,z) in S0 and input_diff_vec[x]!=0:
                        dis_flag=0
                    if get_locate(x,y,z) in S1 and input_diff_vec[x]!=1:
                        dis_flag=0
                if TDTT[input_diff_T][output_diff]==-1 and dis_flag==1:
                    S_discard.append(input_diff_vec)
            #Extract checkbits
            check_bits=[]
            hm=Integer(T[0])+Integer(T[1])+Integer(T[2])+Integer(T[3])+Integer(T[4])
            for x in range(5):
                if T[x]==0:
                    continue
                if get_locate(x,y,z) in S0 or get_locate(x,y,z) in S1:
                    continue
                pair_cnt=0
                for i in range(len(S_discard)):
                    for j in range(i+1,len(S_discard)):
                        if Vector_equal(S_discard[i],S_discard[j],T,x)==1:
                            pair_cnt+=1
                            continue
                if pair_cnt!=len(S_discard)/2:
                    check_bits.append(get_locate(x,y,z))
            for i in range(len(S_discard)):
                discard_v=[0 for _ in range(5)]
                for j in range(len(check_bits)):
                    discard_v[pusai0(check_bits[j])]=S_discard[i][pusai0(check_bits[j])]
                discard_vector.add(tuple(discard_v))
            #check properties
            discard_vector=list(discard_vector)
            for j in range(len(discard_vector)):
                discard_condition=[]
                position_flag=1
                for i in check_bits:
                    dv_i=discard_vector[j][pusai0(i)]
                    if sigma_inverse(i)>=828:
                        x_=pusai0(sigma_inverse(i))
                        z_=pusai2(sigma_inverse(i))
                        for y_ in range(5):
                            if get_locate(x_,y_,z_)>=828 and sigma(get_locate(x_,y_,z_)) in S0:
                                flag=get_locate(x_,y_,z_)
                                flag_value=0
                            if get_locate(x_,y_,z_)>=828 and sigma(get_locate(x_,y_,z_)) in S1:
                                flag=get_locate(x_,y_,z_)
                                flag_value=1
                        if flag<0:
                            position_flag=0
                        else:
                            if GF(2)(flag_value)+GF(2)(dv_i)==0:
                                discard_condition.append("a0[{}]!=a0[{}]".format(sigma_inverse(i),flag))
                            else:
                                discard_condition.append("a0[{}]=a0[{}]".format(sigma_inverse(i),flag))
                    else:
                        position_flag=0
                if position_flag==1 and len(check_bits)!=0:
                    discard_position.append([y,z])
                    print(y,z)
                    print_condition(discard_condition)

S0,S1,matrix0,matrix1=Type1_condition()
a0_s,b0_s=Deduce_known_bit_position()
Derive_new_conditions()