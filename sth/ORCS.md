### Mitigating the impact of selfish routing: An optimal-ratio control  scheme (ORCS) inspired by autonomous driving(2018)

#### 1 Network Data

|         Attribute         |     Value      |
| :-----------------------: | :------------: |
|          Network          |  Sioux Falls   |
| Links, nodes and OD pairs | 76, 24 and 528 |
|       TSTT under UE       |    7480225     |
|       TSTT under SO       |    7194922     |
|     Total flow demand     |     360600     |

#### 2 Model

Control a proportion of travelers for each OD pair to gain TSTT savings.

<img src="C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20250408151845811.png" alt="image-20250408151845811" style="zoom: 80%;" />

where $\gamma$ is control intensity. The objective includes the total controlled demand(abbreviated as "First term") and TSTT(Abbreviated as "Second term").

#### 3 Questions

**Q1: Solutions in section 5.1.3(ORCS with full control potential) may not be optimal.**

+ Control potential $C_{max}=1$.
+ Control coefficient $\gamma=0.1$.

<img src="C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20250408152107157.png" alt="image-20250408152107157" style="zoom: 67%;" />(First paragraph in page 9)

<img src="C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20250408191609736.png" alt="image-20250408191609736" style="zoom: 80%;" />(Fig 3 in page 10)



**The results in the paper are:**

+ First term: $\gamma||\tilde{q}||_1=0.1*360600*0.07 =2524$.

+ Second term: $TSTT=7480225-(7480225-7194922)*0.85=7237718$.
+ ==Objective values==: $7237718+2524=7240242$

**However, if we assume that all travelers are controlled(which equals to the SO state), the results are:**

+ First term: $\gamma||\tilde{q_{so}}||_1=0.1*360600 =36060$.
+ Second term: $TSTT_{SO}=7194922$.
+ ==Objective values==: $7194922+36060=7230982<7240242$

This means that the objective values can be decreased.

**Replicated results:**

+ First term: $\gamma||\tilde{q_{*}}||_1=0.1*154461 =15446$.
+ Second term: $TSTT_{*}=7195742$.
+ ==Objective values==: $7195742+15446=7211188<7230982$.

<img src="C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20250408171757995.png" alt="image-20250408171757995" style="zoom:67%;" />

<img src="C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20250408173616842.png" alt="image-20250408173616842" style="zoom:67%;" />

**Q2: Inconsistency in context**:

In Fig 3, when $\gamma=0.1$, the percent of SO users is more than ==7%==. However, in Fig 4, it is lower than ==7%==.

![image-20250408192250500](C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20250408192250500.png)(Fig 3, where percent of SO users is around 7.2%)

![image-20250408192433729](C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20250408192433729.png)(Fig 4, where percent of SO users is around 6.7%)