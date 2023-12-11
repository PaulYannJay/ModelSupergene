Script for the paper: "The interplay between local adaptation and gene flow contributes to the formation of supergenes"

Usage:

julia MainSupergene.jl NameOutput(string) NumberOfRepetition(integer) Identifier(integer)
e.g. julia MainSupergene.jl TestScript 10 1

The parameter space explored is defined by the arrays in MainSupergene.jl. Change the values to explore other parameters. 
To compute a single simulation (e.g. figure 2B-G), uncomment the section corresponding in MainSupergene.jl and comment the "Sensitivity analysis" section.


The sensitivity analyse function produce two tables:

1: SensitivityAnalysis\_Name.txt, which looks like that:
*param1;param2;param3;param4;param5;rep;state;genotype;frequency 
0.05;0.0;0.0;0.2;0.5;7;0;1;0.0 
0.05;0.0;0.0;0.2;0.5;7;0;2;0.0 
0.05;0.0;0.0;0.2;0.5;7;0;3;0.0 
0.05;0.0;0.0;0.2;0.5;7;0;4;0.0 
0.05;0.0;0.0;0.2;0.5;7;0;5;0.0 
0.05;0.0;0.0;0.2;0.5;7;0;6;0.0 
0.05;0.0;0.0;0.2;0.5;7;0;7;0.0 
0.05;0.0;0.0;0.2;0.5;7;0;8;0.0 
0.05;0.0;0.0;0.2;0.5;7;0;9;0.0* 
This file describe the state of the populations at the end of the migratory phase (state=0) and at the end of the non-migratory phase (state=1), depending on the parameter defined on the script (param1-5), and for each repetition. The state of the populations is described by the frequency of each genotype. Genotypes are numbered from 1 to 128. Please refer to the list at the end of this document to identify each genotype

and 

2: SensitivityAnalysis2\_Name.txt, which looks like that:
*param1;param2;param3;param4;param5;rep;state;deme;chromosome;frequency 
0.05;0.0;0.0;0.2;0.5;7;0;1;1;0.0 
0.05;0.0;0.0;0.2;0.5;7;0;1;2;0.0 
0.05;0.0;0.0;0.2;0.5;7;0;1;3;0.0 
0.05;0.0;0.0;0.2;0.5;7;0;1;4;0.0 
0.05;0.0;0.0;0.2;0.5;7;0;1;5;0.7261452290458091 
0.05;0.0;0.0;0.2;0.5;7;0;1;6;0.0 
0.05;0.0;0.0;0.2;0.5;7;0;1;7;0.0 
0.05;0.0;0.0;0.2;0.5;7;0;1;8;0.2738547709541908 
0.05;0.0;0.0;0.2;0.5;7;0;2;1;0.0*
This file describe the state of the populations at the end of the migratory phase (state=0) and at the end of the non-migratory phase (state=1), depending on the parameter defined the the script (param1-5), and for each repetition. The state of the population is described by the frequency of each haplotype ("chromosome") in each population (deme 1 and deme 2). Haplotypes are numbered from 1 to 8:  
 1: B-A (inversion)  2: b-A (inversion)  3: B-a (inversion)  4: b-a (inversion)  5: A-B  6: A-b 7: a-B 8: a-b


The single simulation function produce two tables:
-	1: SingleTimeSeries\_Name.txt, which looks like that:
*time;genotype;frequency 
0;1;0.0 
0;2;0.0 
0;3;0.0 
0;4;0.0 
0;5;0.0 
0;6;0.0 
0;7;0.0 
0;8;0.0 
0;9;0.0*

and 
- 2: SingleTimeSeries2\_Name.txt, which looks like that:
*time;deme;chromosome;frequency 
0;1;1;0.0 
0;1;2;0.0 
0;1;3;0.0 
0;1;4;0.0 
0;1;5;0.5 
0;1;6;0.0 
0;1;7;0.0 
0;1;8;0.5 
0;2;1;0.0*

These tables describe the state of the populations every generation. The output are either the genotype frequency (in SingleTimeSeries_\*; Please refer to the list at the end of this document to identify each genotype) or in haplotype frequency. 

*Genotype list:
1--> Population=1, Genotype=BA/BA
2--> Population=1, Genotype=BA/bA
3--> Population=1, Genotype=BA/Ba
4--> Population=1, Genotype=BA/ba
5--> Population=1, Genotype=BA/AB
6--> Population=1, Genotype=BA/Ab
7--> Population=1, Genotype=BA/aB
8--> Population=1, Genotype=BA/ab
9--> Population=1, Genotype=bA/BA
10--> Population=1, Genotype=bA/bA
11--> Population=1, Genotype=bA/Ba
12--> Population=1, Genotype=bA/ba
13--> Population=1, Genotype=bA/AB
14--> Population=1, Genotype=bA/Ab
15--> Population=1, Genotype=bA/aB
16--> Population=1, Genotype=bA/ab
17--> Population=1, Genotype=Ba/BA
18--> Population=1, Genotype=Ba/bA
19--> Population=1, Genotype=Ba/Ba
20--> Population=1, Genotype=Ba/ba
21--> Population=1, Genotype=Ba/AB
22--> Population=1, Genotype=Ba/Ab
23--> Population=1, Genotype=Ba/aB
24--> Population=1, Genotype=Ba/ab
25--> Population=1, Genotype=ba/BA
26--> Population=1, Genotype=ba/bA
27--> Population=1, Genotype=ba/Ba
28--> Population=1, Genotype=ba/ba
29--> Population=1, Genotype=ba/AB
30--> Population=1, Genotype=ba/Ab
31--> Population=1, Genotype=ba/aB
32--> Population=1, Genotype=ba/ab
33--> Population=1, Genotype=AB/BA
34--> Population=1, Genotype=AB/bA
35--> Population=1, Genotype=AB/Ba
36--> Population=1, Genotype=AB/ba
37--> Population=1, Genotype=AB/AB
38--> Population=1, Genotype=AB/Ab
39--> Population=1, Genotype=AB/aB
40--> Population=1, Genotype=AB/ab
41--> Population=1, Genotype=Ab/BA
42--> Population=1, Genotype=Ab/bA
43--> Population=1, Genotype=Ab/Ba
44--> Population=1, Genotype=Ab/ba
45--> Population=1, Genotype=Ab/AB
46--> Population=1, Genotype=Ab/Ab
47--> Population=1, Genotype=Ab/aB
48--> Population=1, Genotype=Ab/ab
49--> Population=1, Genotype=aB/BA
50--> Population=1, Genotype=aB/bA
51--> Population=1, Genotype=aB/Ba
52--> Population=1, Genotype=aB/ba
53--> Population=1, Genotype=aB/AB
54--> Population=1, Genotype=aB/Ab
55--> Population=1, Genotype=aB/aB
56--> Population=1, Genotype=aB/ab
57--> Population=1, Genotype=ab/BA
58--> Population=1, Genotype=ab/bA
59--> Population=1, Genotype=ab/Ba
60--> Population=1, Genotype=ab/ba
61--> Population=1, Genotype=ab/AB
62--> Population=1, Genotype=ab/Ab
63--> Population=1, Genotype=ab/aB
64--> Population=1, Genotype=ab/ab
65--> Population=2, Genotype=BA/BA
66--> Population=2, Genotype=BA/bA
67--> Population=2, Genotype=BA/Ba
68--> Population=2, Genotype=BA/ba
69--> Population=2, Genotype=BA/AB
70--> Population=2, Genotype=BA/Ab
71--> Population=2, Genotype=BA/aB
72--> Population=2, Genotype=BA/ab
73--> Population=2, Genotype=bA/BA
74--> Population=2, Genotype=bA/bA
75--> Population=2, Genotype=bA/Ba
76--> Population=2, Genotype=bA/ba
77--> Population=2, Genotype=bA/AB
78--> Population=2, Genotype=bA/Ab
79--> Population=2, Genotype=bA/aB
80--> Population=2, Genotype=bA/ab
81--> Population=2, Genotype=Ba/BA
82--> Population=2, Genotype=Ba/bA
83--> Population=2, Genotype=Ba/Ba
84--> Population=2, Genotype=Ba/ba
85--> Population=2, Genotype=Ba/AB
86--> Population=2, Genotype=Ba/Ab
87--> Population=2, Genotype=Ba/aB
88--> Population=2, Genotype=Ba/ab
89--> Population=2, Genotype=ba/BA
90--> Population=2, Genotype=ba/bA
91--> Population=2, Genotype=ba/Ba
92--> Population=2, Genotype=ba/ba
93--> Population=2, Genotype=ba/AB
94--> Population=2, Genotype=ba/Ab
95--> Population=2, Genotype=ba/aB
96--> Population=2, Genotype=ba/ab
97--> Population=2, Genotype=AB/BA
98--> Population=2, Genotype=AB/bA
99--> Population=2, Genotype=AB/Ba
100--> Population=2, Genotype=AB/ba
101--> Population=2, Genotype=AB/AB
102--> Population=2, Genotype=AB/Ab
103--> Population=2, Genotype=AB/aB
104--> Population=2, Genotype=AB/ab
105--> Population=2, Genotype=Ab/BA
106--> Population=2, Genotype=Ab/bA
107--> Population=2, Genotype=Ab/Ba
108--> Population=2, Genotype=Ab/ba
109--> Population=2, Genotype=Ab/AB
110--> Population=2, Genotype=Ab/Ab
111--> Population=2, Genotype=Ab/aB
112--> Population=2, Genotype=Ab/ab
113--> Population=2, Genotype=aB/BA
114--> Population=2, Genotype=aB/bA
115--> Population=2, Genotype=aB/Ba
116--> Population=2, Genotype=aB/ba
117--> Population=2, Genotype=aB/AB
118--> Population=2, Genotype=aB/Ab
119--> Population=2, Genotype=aB/aB
120--> Population=2, Genotype=aB/ab
121--> Population=2, Genotype=ab/BA
122--> Population=2, Genotype=ab/bA
123--> Population=2, Genotype=ab/Ba
124--> Population=2, Genotype=ab/ba
125--> Population=2, Genotype=ab/AB
126--> Population=2, Genotype=ab/Ab
127--> Population=2, Genotype=ab/aB
128--> Population=2, Genotype=ab/ab
